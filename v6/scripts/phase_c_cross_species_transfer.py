#!/usr/bin/env python3
"""
Phase C: Cross-Species Transfer Learning

Train PCA-LR on mouse tissue gene expression → predict human pre/post flight
via orthologs. Tests whether mouse spaceflight classifiers transfer to human.
"""

import sys
import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score
from datetime import datetime

sys.path.insert(0, os.path.dirname(__file__))
from v6_utils import (
    load_ortholog_map, load_ensembl_to_symbol, load_cfrna_feature_matrix,
    ALL_TISSUES, V6_EVAL, GENELAB_BASE, save_json
)

# Add v4 scripts to path for data loading
sys.path.insert(0, str(GENELAB_BASE / "v4" / "scripts"))


def load_mouse_expression(tissue):
    """Load mouse tissue expression + labels using v4 infrastructure."""
    from v4_utils import load_metadata, load_gene_features, align_features_with_meta, encode_labels

    meta = load_metadata(tissue)
    genes = load_gene_features(tissue)
    genes, meta = align_features_with_meta(genes, meta)
    y, valid = encode_labels(meta)

    return genes[valid], y[valid], meta[valid]


def find_common_genes(mouse_genes_df, human_genes, orth_map, ens2sym):
    """Find genes present in both mouse expression and human feature matrix.

    Returns:
        common_mouse_cols: list of mouse column names (Ensembl or symbol)
        common_human_cols: list of corresponding human gene names
    """
    common_mouse = []
    common_human = []
    human_gene_set = set(human_genes)

    for col in mouse_genes_df.columns:
        # Map to symbol first
        symbol = ens2sym.get(col, col) if col.startswith("ENSMUSG") else col

        # Map to human
        human = orth_map.get(symbol)
        if human is None:
            human = orth_map.get(symbol.capitalize())
        if human is None:
            human = orth_map.get(symbol.upper())

        if human and human in human_gene_set:
            common_mouse.append(col)
            common_human.append(human)

    return common_mouse, common_human


def main():
    print("=" * 60)
    print("Phase C: Cross-Species Transfer Learning")
    print("=" * 60)

    # Load human cfRNA feature matrix
    print("\n1. Loading human cfRNA data...")
    X_human, meta_human = load_cfrna_feature_matrix()
    human_genes = list(X_human.columns)
    print(f"  cfRNA: {X_human.shape[0]} samples × {X_human.shape[1]} genes")
    print(f"  Phases: {meta_human['phase'].value_counts().to_dict()}")

    # Binary labels: Pre=0, Post=1 (exclude Recovery)
    # Also try Pre=0, Post+Recovery=1
    analyses = {
        "pre_vs_post": {
            "mask": meta_human["phase"].isin(["Pre", "Post"]),
            "label_fn": lambda phase: 1 if phase == "Post" else 0,
        },
        "pre_vs_post_recovery": {
            "mask": meta_human["phase"].isin(["Pre", "Post", "Recovery"]),
            "label_fn": lambda phase: 0 if phase == "Pre" else 1,
        },
    }

    orth_map = load_ortholog_map()
    ens2sym = load_ensembl_to_symbol()

    results = {}

    for analysis_name, analysis_cfg in analyses.items():
        print(f"\n{'='*40}")
        print(f"Analysis: {analysis_name}")
        print(f"{'='*40}")

        mask = analysis_cfg["mask"]
        X_h = X_human[mask].values
        y_h = meta_human[mask]["phase"].map(analysis_cfg["label_fn"]).values
        crews_h = meta_human[mask]["crew_id"].values
        n_pos = int(y_h.sum())
        n_neg = int(len(y_h) - n_pos)
        print(f"  Human samples: {len(y_h)} (pos={n_pos}, neg={n_neg})")

        tissue_results = {}

        for tissue in ALL_TISSUES:
            print(f"\n  Tissue: {tissue}")
            try:
                mouse_genes, mouse_y, mouse_meta = load_mouse_expression(tissue)
            except Exception as e:
                print(f"    ERROR loading: {e}")
                continue

            # Find common genes
            common_mouse, common_human = find_common_genes(
                mouse_genes, human_genes, orth_map, ens2sym)

            n_common = len(common_mouse)
            if n_common < 10:
                print(f"    Only {n_common} common genes, SKIPPING")
                continue

            print(f"    Common genes: {n_common}")

            # Mouse feature matrix (common genes only)
            X_mouse = mouse_genes[common_mouse].values.astype(np.float32)
            y_mouse = mouse_y.values.astype(int)

            # Human feature matrix (common genes, same order)
            human_col_idx = [human_genes.index(h) for h in common_human]
            X_h_common = X_h[:, human_col_idx].astype(np.float32)

            # Build model: StandardScaler → PCA → LR
            n_components = min(50, n_common - 1, X_mouse.shape[0] - 1)

            # Fit StandardScaler on MOUSE data only
            scaler = StandardScaler()
            X_mouse_scaled = scaler.fit_transform(X_mouse)
            X_h_scaled = scaler.transform(X_h_common)

            # Fit PCA on MOUSE data only
            pca = PCA(n_components=n_components, random_state=42)
            X_mouse_pca = pca.fit_transform(X_mouse_scaled)
            X_h_pca = pca.transform(X_h_scaled)

            # Fit LR on MOUSE data
            lr = LogisticRegression(max_iter=1000, random_state=42)
            lr.fit(X_mouse_pca, y_mouse)

            # Predict on HUMAN data
            y_h_score = lr.predict_proba(X_h_pca)[:, 1]

            # Overall AUROC
            try:
                auroc = roc_auc_score(y_h, y_h_score)
            except ValueError:
                auroc = np.nan

            # LOCO-CV on human (leave-one-crew-out)
            unique_crews = np.unique(crews_h)
            loco_aurocs = []
            for test_crew in unique_crews:
                test_mask = crews_h == test_crew
                y_test = y_h[test_mask]
                y_score_test = y_h_score[test_mask]
                if len(np.unique(y_test)) >= 2:
                    try:
                        auc = roc_auc_score(y_test, y_score_test)
                        loco_aurocs.append({
                            "crew": str(test_crew), "auroc": round(auc, 4),
                            "n_test": int(len(y_test))
                        })
                    except ValueError:
                        pass

            # Permutation test: shuffle mouse labels, retrain, predict human
            rng = np.random.default_rng(42)
            n_perm = 1000
            perm_aurocs = []
            for _ in range(n_perm):
                y_perm = rng.permutation(y_mouse)
                lr_perm = LogisticRegression(max_iter=1000, random_state=42)
                lr_perm.fit(X_mouse_pca, y_perm)
                y_perm_score = lr_perm.predict_proba(X_h_pca)[:, 1]
                try:
                    perm_auc = roc_auc_score(y_h, y_perm_score)
                    perm_aurocs.append(perm_auc)
                except ValueError:
                    pass

            if perm_aurocs and not np.isnan(auroc):
                perm_p = (np.sum(np.array(perm_aurocs) >= auroc) + 1) / (len(perm_aurocs) + 1)
            else:
                perm_p = 1.0

            # Mouse internal AUROC (sanity check via LOMO-CV)
            # Quick 3-fold CV for speed
            from sklearn.model_selection import StratifiedKFold
            mouse_aurocs = []
            skf = StratifiedKFold(n_splits=min(3, len(np.unique(y_mouse))),
                                  shuffle=True, random_state=42)
            for train_idx, test_idx in skf.split(X_mouse_pca, y_mouse):
                lr_cv = LogisticRegression(max_iter=1000, random_state=42)
                lr_cv.fit(X_mouse_pca[train_idx], y_mouse[train_idx])
                try:
                    auc_cv = roc_auc_score(
                        y_mouse[test_idx],
                        lr_cv.predict_proba(X_mouse_pca[test_idx])[:, 1])
                    mouse_aurocs.append(auc_cv)
                except ValueError:
                    pass

            tissue_results[tissue] = {
                "n_common_genes": n_common,
                "n_pca_components": n_components,
                "n_mouse_train": int(X_mouse.shape[0]),
                "n_human_test": int(len(y_h)),
                "transfer_auroc": round(float(auroc), 4) if not np.isnan(auroc) else None,
                "loco_aurocs": loco_aurocs,
                "mean_loco_auroc": round(float(np.mean([f["auroc"] for f in loco_aurocs])), 4)
                    if loco_aurocs else None,
                "perm_p": round(float(perm_p), 4),
                "significant": perm_p < 0.05,
                "mouse_internal_auroc": round(float(np.mean(mouse_aurocs)), 4)
                    if mouse_aurocs else None,
                "perm_auroc_mean": round(float(np.mean(perm_aurocs)), 4)
                    if perm_aurocs else None,
            }

            sig = "*" if perm_p < 0.05 else ""
            print(f"    Transfer AUROC: {auroc:.4f}, perm_p={perm_p:.4f}{sig}")
            print(f"    Mouse internal: {np.mean(mouse_aurocs):.4f}" if mouse_aurocs else "")

        results[analysis_name] = tissue_results

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    for analysis_name, tissue_results in results.items():
        print(f"\n  {analysis_name}:")
        for tissue, r in sorted(tissue_results.items(),
                                key=lambda x: x[1].get("transfer_auroc", 0) or 0,
                                reverse=True):
            auc = r.get("transfer_auroc", "N/A")
            sig = "*" if r.get("significant") else ""
            print(f"    {tissue:>15s}: AUROC={auc}{sig} (n_genes={r['n_common_genes']})")

    # Assemble output
    output = {
        "analyses": results,
        "ortholog_map_size": len(orth_map),
        "cfRNA_genes": len(human_genes),
        "timestamp": datetime.now().isoformat(),
    }

    save_json(output, V6_EVAL / "V6_C_cross_species_transfer.json")
    print("\nDone!")


if __name__ == "__main__":
    main()

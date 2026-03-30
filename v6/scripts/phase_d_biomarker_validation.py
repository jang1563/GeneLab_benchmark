#!/usr/bin/env python3
"""
Phase D: Biomarker Panel Validation — Check v5 20-gene panel in human cfRNA

For each panel gene:
- Check if detected in cfRNA (26,845 genes)
- Check if DE (FDR<0.05)
- Check if DRR (in 466 DRR list)
- Check direction concordance with mouse
"""

import sys
import os
import numpy as np
from scipy import stats
from datetime import datetime

sys.path.insert(0, os.path.dirname(__file__))
from v6_utils import (
    load_ortholog_map, load_ensembl_to_symbol, map_mouse_to_human,
    load_cfrna_de, load_cfrna_drr, load_cfrna_feature_matrix,
    load_biomarker_panel, V6_EVAL, save_json
)


def main():
    print("=" * 60)
    print("Phase D: Biomarker Panel Validation")
    print("=" * 60)

    # Load data
    print("\n1. Loading data...")
    orth_map = load_ortholog_map()
    ens2sym = load_ensembl_to_symbol()
    cfrna_de = load_cfrna_de()
    drr_genes = load_cfrna_drr()
    biomarker = load_biomarker_panel()

    # Try loading feature matrix for classification
    try:
        X_cfrna, meta_cfrna = load_cfrna_feature_matrix()
        has_feature_matrix = True
        print(f"  cfRNA feature matrix: {X_cfrna.shape[0]} samples × {X_cfrna.shape[1]} genes")
    except Exception as e:
        has_feature_matrix = False
        print(f"  cfRNA feature matrix not available: {e}")

    panel = biomarker.get("panel", [])
    print(f"  Biomarker panel: {len(panel)} genes")

    # Map each panel gene to human
    print("\n2. Checking panel genes in human cfRNA...")
    gene_results = []

    for g in panel:
        mouse_gene = g.get("gene", "")
        human_symbol = g.get("human_symbol", "")
        score = g.get("score", 0)
        tissues = g.get("tissues", [])
        drugs = g.get("drugs", [])

        # Try to find valid human ortholog
        human_gene = None

        # First try the stored human_symbol
        if human_symbol and human_symbol in cfrna_de.index:
            human_gene = human_symbol
        else:
            # Try ortholog mapping
            mapped, _ = map_mouse_to_human([mouse_gene], orth_map, ens2sym)
            if mapped:
                h = list(mapped.values())[0]
                if h in cfrna_de.index:
                    human_gene = h

        # If still not found, try case variations
        if human_gene is None:
            for candidate in [human_symbol, mouse_gene.upper(), mouse_gene.capitalize()]:
                if candidate and candidate in cfrna_de.index:
                    human_gene = candidate
                    break

        result = {
            "rank": g.get("rank", 0),
            "mouse_gene": mouse_gene,
            "human_symbol_stored": human_symbol,
            "human_gene_matched": human_gene,
            "panel_score": score,
            "tissues": tissues,
            "n_drugs": len(drugs),
            "detected_in_cfrna": human_gene is not None,
            "is_de": False,
            "de_fdr": None,
            "de_diff": None,
            "is_drr": False,
            "direction": None,
            "concordant_with_mouse": None,
        }

        if human_gene is not None:
            row = cfrna_de.loc[human_gene]
            fdr = row.get("edge_pre_vs_flight_fdr", 1.0)
            diff = row.get("edge_pre_vs_flight_diff", 0.0)

            result["is_de"] = float(fdr) < 0.05
            result["de_fdr"] = round(float(fdr), 6)
            result["de_diff"] = round(float(diff), 4)
            result["is_drr"] = human_gene in drr_genes
            result["direction"] = "up_in_flight" if diff > 0 else "down_in_flight"

        gene_results.append(result)

        status = "✓" if result["detected_in_cfrna"] else "✗"
        de_str = f"DE(FDR={result['de_fdr']:.4f})" if result["is_de"] else "not_DE"
        drr_str = "DRR" if result["is_drr"] else ""
        matched = result["human_gene_matched"] or "NO_ORTHOLOG"
        print(f"  {status} {mouse_gene:>12s} → {matched:>12s} | {de_str} {drr_str}")

    # Summary statistics
    n_detected = sum(1 for r in gene_results if r["detected_in_cfrna"])
    n_de = sum(1 for r in gene_results if r["is_de"])
    n_drr = sum(1 for r in gene_results if r["is_drr"])
    n_total = len(gene_results)

    print(f"\n3. Summary:")
    print(f"  Detected in cfRNA: {n_detected}/{n_total} ({n_detected/n_total:.0%})")
    print(f"  DE (FDR<0.05):     {n_de}/{n_total}")
    print(f"  DRR (466 list):    {n_drr}/{n_total}")

    # Enrichment vs random: are panel genes enriched among DRR?
    print("\n4. DRR enrichment vs random 20-gene draws...")
    # Build universe of human genes present in cfRNA with mouse orthologs
    human_genes_with_orth = set(orth_map.values()) & set(cfrna_de.index)
    universe_size = len(human_genes_with_orth)
    drr_in_universe = drr_genes & human_genes_with_orth
    panel_human_genes = {r["human_gene_matched"] for r in gene_results
                        if r["detected_in_cfrna"]}
    panel_drr = panel_human_genes & drr_in_universe
    n_panel_drr = len(panel_drr)

    # Hypergeometric test
    p_hyper = stats.hypergeom.sf(
        n_panel_drr - 1, universe_size, len(drr_in_universe), n_detected)

    expected = n_detected * len(drr_in_universe) / universe_size if universe_size > 0 else 0
    fold_enrich = n_panel_drr / expected if expected > 0 else 0

    # Permutation test (10K random draws of same size)
    rng = np.random.default_rng(42)
    universe_list = list(human_genes_with_orth)
    n_perm = 10000
    count_ge = 0
    for _ in range(n_perm):
        random_panel = set(rng.choice(universe_list, size=n_detected, replace=False))
        if len(random_panel & drr_in_universe) >= n_panel_drr:
            count_ge += 1
    p_perm = (count_ge + 1) / (n_perm + 1)

    print(f"  Panel DRR overlap: {n_panel_drr}/{n_detected} detected genes")
    print(f"  Expected: {expected:.1f}, Fold enrichment: {fold_enrich:.2f}")
    print(f"  Hypergeometric p: {p_hyper:.4f}")
    print(f"  Permutation p: {p_perm:.4f}")

    # Panel classification on cfRNA (if feature matrix available)
    panel_classification = None
    if has_feature_matrix and n_detected >= 5:
        print("\n5. Panel classification on cfRNA...")
        panel_genes_in_matrix = [r["human_gene_matched"] for r in gene_results
                                 if r["detected_in_cfrna"]
                                 and r["human_gene_matched"] in X_cfrna.columns]

        if len(panel_genes_in_matrix) >= 3:
            from sklearn.linear_model import LogisticRegression
            from sklearn.preprocessing import StandardScaler
            from sklearn.pipeline import Pipeline
            from sklearn.metrics import roc_auc_score

            # Binary: Pre=0, Post=1 (exclude Recovery for cleaner comparison)
            mask = meta_cfrna["phase"].isin(["Pre", "Post"])
            X_panel = X_cfrna.loc[mask, panel_genes_in_matrix].values
            y = (meta_cfrna.loc[mask, "phase"] == "Post").astype(int).values

            if len(np.unique(y)) == 2 and X_panel.shape[0] >= 6:
                # LOCO-CV (leave-one-crew-out)
                crews = meta_cfrna.loc[mask, "crew_id"].values
                unique_crews = np.unique(crews)
                fold_aurocs = []

                for test_crew in unique_crews:
                    test_mask = crews == test_crew
                    train_mask = ~test_mask

                    if len(np.unique(y[train_mask])) < 2:
                        continue
                    if len(np.unique(y[test_mask])) < 2:
                        continue

                    model = Pipeline([
                        ("scaler", StandardScaler()),
                        ("lr", LogisticRegression(max_iter=1000, random_state=42)),
                    ])
                    model.fit(X_panel[train_mask], y[train_mask])
                    y_score = model.predict_proba(X_panel[test_mask])[:, 1]
                    try:
                        auc = roc_auc_score(y[test_mask], y_score)
                        fold_aurocs.append({"crew": str(test_crew), "auroc": round(auc, 4)})
                    except ValueError:
                        pass

                if fold_aurocs:
                    mean_auc = np.mean([f["auroc"] for f in fold_aurocs])
                    panel_classification = {
                        "n_panel_genes_in_matrix": len(panel_genes_in_matrix),
                        "panel_genes_used": panel_genes_in_matrix,
                        "n_samples": int(X_panel.shape[0]),
                        "n_pre": int(np.sum(y == 0)),
                        "n_post": int(np.sum(y == 1)),
                        "cv_strategy": "leave-one-crew-out",
                        "fold_aurocs": fold_aurocs,
                        "mean_auroc": round(float(mean_auc), 4),
                    }
                    print(f"  Panel genes in matrix: {len(panel_genes_in_matrix)}")
                    print(f"  LOCO-CV mean AUROC: {mean_auc:.4f}")
                    for f in fold_aurocs:
                        print(f"    Crew {f['crew']}: AUROC={f['auroc']}")
            else:
                print(f"  Insufficient samples for classification (n={X_panel.shape[0]})")
        else:
            print(f"  Only {len(panel_genes_in_matrix)} panel genes in 500-gene matrix")

    # Assemble output
    output = {
        "panel_size": n_total,
        "n_detected_in_cfrna": n_detected,
        "detection_rate": round(n_detected / n_total, 4) if n_total > 0 else 0,
        "n_de_fdr05": n_de,
        "n_drr": n_drr,
        "gene_results": gene_results,
        "enrichment": {
            "universe_size": universe_size,
            "drr_in_universe": len(drr_in_universe),
            "panel_drr_overlap": n_panel_drr,
            "expected_overlap": round(expected, 2),
            "fold_enrichment": round(fold_enrich, 3),
            "hypergeometric_p": float(p_hyper),
            "permutation_p": float(p_perm),
            "significant": p_hyper < 0.05,
        },
        "panel_classification": panel_classification,
        "timestamp": datetime.now().isoformat(),
    }

    save_json(output, V6_EVAL / "V6_D_biomarker_validation.json")
    print("\nDone!")


if __name__ == "__main__":
    main()

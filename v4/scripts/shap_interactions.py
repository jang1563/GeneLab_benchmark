#!/usr/bin/env python3
"""
shap_interactions.py — GeneLabBench v4 Phase 4: SHAP interaction terms

Computes SHAP interaction values for tree-based methods (RF, XGB) only.
TreeExplainer.shap_interaction_values() gives pairwise gene-gene interactions.

Pre-filters to top-2000 genes by variance to avoid OOM on 17K × 17K matrix.

Usage:
  python shap_interactions.py --tissue kidney --method xgb
  python shap_interactions.py  # all tree-method tissue combinations
"""

import json
import sys
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

warnings.filterwarnings("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parent))

from v4_utils import V4_EVAL_DIR, get_folds
from classifier_registry import (
    get_classifier, adapt_pca_components, SHAP_STRATEGY, CLASSIFIERS,
)

# Ensembl → Symbol map (same as shap_multi_method.py)
BASE_DIR = Path(__file__).resolve().parent.parent.parent
SYMBOL_MAP_PATH = BASE_DIR / "processed" / "ensembl_symbol_map.csv"


def load_symbol_map():
    """Load ENSEMBL → gene symbol map."""
    if not SYMBOL_MAP_PATH.exists():
        return {}
    df = pd.read_csv(SYMBOL_MAP_PATH)
    return dict(zip(df["ENSEMBL"], df["SYMBOL"]))

# Top-2 methods per tissue — filter to tree methods only
TOP2_METHODS = {
    "liver":         ["svm_rbf", "mlp"],
    "gastrocnemius": ["elasticnet_lr", "pca_lr"],
    "kidney":        ["xgb", "rf"],
    "thymus":        ["knn", "pca_lr"],
    "eye":           ["pca_lr", "xgb"],
    "skin":          ["elasticnet_lr", "pca_lr"],
    "lung":          ["elasticnet_lr", "pca_lr"],
    "colon":         ["pca_lr", "xgb"],
}

TREE_COMBOS = [(t, m) for t, methods in TOP2_METHODS.items()
               for m in methods if SHAP_STRATEGY.get(m) == "tree"]

# Top-K genes for interaction analysis (avoid n_features² memory)
INTERACTION_TOP_K = 2000


def run_interactions(tissue, method_key, top_k=INTERACTION_TOP_K, seed=42, verbose=True):
    """Compute SHAP interaction values for a tree-based method.

    Returns top gene-gene interaction pairs aggregated across folds.
    """
    import shap

    folds = get_folds(tissue)

    if verbose:
        print(f"\nInteractions: {tissue} × {method_key}")
        print(f"  Folds: {len(folds)}, top_k: {top_k}")

    fold_interactions = []
    gene_names_sub = None

    for fi, fold in enumerate(folds):
        X_train = fold["train_X"]
        y_train = fold["train_y"]
        X_test = fold["test_X"]
        gene_names = fold["gene_names"]

        if len(np.unique(fold["test_y"])) < 2:
            continue

        # Variance-based gene subset (computed on training data)
        n_features = X_train.shape[1]
        if n_features > top_k:
            variances = np.var(X_train, axis=0)
            top_idx = np.argsort(variances)[-top_k:]
            top_idx = np.sort(top_idx)
            X_train_sub = X_train[:, top_idx]
            X_test_sub = X_test[:, top_idx]
            gene_names_sub = [gene_names[i] for i in top_idx]
        else:
            X_train_sub = X_train
            X_test_sub = X_test
            gene_names_sub = list(gene_names)

        # Train model on subset
        _, model = get_classifier(method_key, seed=seed)
        model.fit(X_train_sub, y_train)

        # SHAP interaction values
        explainer = shap.TreeExplainer(model)

        # Subsample test set (interaction matrix is large: n_test × K × K)
        max_test = min(20, X_test_sub.shape[0])
        if X_test_sub.shape[0] > max_test:
            rng = np.random.default_rng(seed + fi)
            idx = rng.choice(X_test_sub.shape[0], max_test, replace=False)
            X_test_sub = X_test_sub[idx]

        if verbose:
            print(f"  Fold {fold['fold_name']}: {X_test_sub.shape[0]} test × {X_test_sub.shape[1]} genes", end="")

        try:
            interaction_values = explainer.shap_interaction_values(X_test_sub)

            # Handle different shap return formats:
            # - list [class0_3d, class1_3d] for RF (older shap)
            # - 4D array (n_samples, n_features, n_features, n_classes) for RF (newer shap)
            # - 3D array (n_samples, n_features, n_features) for XGB
            if isinstance(interaction_values, list):
                iv = interaction_values[1]  # class 1 = FLT
            elif interaction_values.ndim == 4:
                iv = interaction_values[:, :, :, 1]  # class 1
            else:
                iv = interaction_values

            assert iv.ndim == 3, f"Expected 3D interactions, got {iv.ndim}D"

            # iv shape: (n_test, n_genes, n_genes)
            # Mean absolute interaction across test samples
            mean_abs_interaction = np.mean(np.abs(iv), axis=0)  # (K, K)
            fold_interactions.append(mean_abs_interaction)

            if verbose:
                print(f" — OK (interaction matrix {mean_abs_interaction.shape})")
        except Exception as e:
            if verbose:
                print(f" — FAILED: {e}")
            continue

    if not fold_interactions:
        return {"error": f"No valid folds for interactions {tissue}/{method_key}"}

    # Resolve Ensembl → gene symbols
    sym_map = load_symbol_map()
    gene_symbols = [sym_map.get(g, g) for g in gene_names_sub]

    # Average interaction matrix across folds
    mean_interaction = np.mean(fold_interactions, axis=0)

    # Extract top interaction pairs (off-diagonal only)
    n_genes = mean_interaction.shape[0]
    pairs = []
    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            pairs.append({
                "gene_a": gene_symbols[i],
                "gene_b": gene_symbols[j],
                "interaction_strength": round(float(mean_interaction[i, j]), 6),
            })

    # Sort by interaction strength (descending)
    pairs.sort(key=lambda x: -x["interaction_strength"])

    # Also get diagonal (main effects in interaction space)
    main_effects = {gene_symbols[i]: round(float(mean_interaction[i, i]), 6)
                    for i in range(n_genes)}
    main_effects = dict(sorted(main_effects.items(), key=lambda x: -x[1]))

    result = {
        "tissue": tissue,
        "method": method_key,
        "n_genes": n_genes,
        "gene_filter": f"top_{top_k}_variance" if n_genes <= top_k else "all",
        "n_folds_used": len(fold_interactions),
        "top_50_interactions": pairs[:50],
        "top_20_main_effects": dict(list(main_effects.items())[:20]),
        "n_total_pairs": len(pairs),
        "timestamp": datetime.now().isoformat(),
    }

    return result


def main():
    parser = argparse.ArgumentParser(description="SHAP interaction terms (tree methods)")
    parser.add_argument("--tissue", type=str, default=None)
    parser.add_argument("--method", type=str, default=None)
    parser.add_argument("--top-k", type=int, default=INTERACTION_TOP_K)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    V4_EVAL_DIR.mkdir(parents=True, exist_ok=True)

    if args.tissue and args.method:
        combos = [(args.tissue, args.method)]
    else:
        combos = TREE_COMBOS

    print(f"SHAP Interactions: {len(combos)} tree-method combinations")
    for t, m in combos:
        print(f"  {t} × {m}")

    for tissue, method in combos:
        output_path = V4_EVAL_DIR / f"SHAP_interactions_{tissue}_{method}.json"

        if output_path.exists():
            print(f"\n  SKIP {tissue}/{method}: exists")
            continue

        result = run_interactions(tissue, method, top_k=args.top_k, seed=args.seed)

        with open(output_path, 'w') as f:
            json.dump(result, f, indent=2)

        if "error" not in result:
            top_pair = result["top_50_interactions"][0] if result["top_50_interactions"] else {}
            print(f"  → {output_path.name}: {result['n_genes']} genes, "
                  f"top pair: {top_pair.get('gene_a', '?')}↔{top_pair.get('gene_b', '?')}")

    print(f"\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

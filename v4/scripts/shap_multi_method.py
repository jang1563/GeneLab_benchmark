#!/usr/bin/env python3
"""
shap_multi_method.py — GeneLabBench v4 Phase 4: SHAP interpretability

Computes feature importance for top-2 methods per tissue (16 combinations).
SHAP strategy by classifier family (DD-23):
  - Tree (RF, XGB): shap.TreeExplainer
  - Linear (ElasticNet-LR): |coef| in standardized space × scaler.scale_
  - PCA-Linear (PCA-LR): |lr_coef @ PCA_components| × scaler.scale_
  - Kernel (SVM-RBF, kNN, MLP, TabNet): shap.KernelExplainer(kmeans=100)
    - Kernel methods pre-filter to top-2000 genes by training-set variance

Output per tissue × method:
  v4/evaluation/SHAP_{tissue}_{method}.json

Usage:
  python shap_multi_method.py --tissue liver --method svm_rbf
  python shap_multi_method.py  # all 16 top-2 combinations
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

from v4_utils import TISSUE_MISSIONS, V4_EVAL_DIR, get_folds
from classifier_registry import (
    get_classifier, adapt_pca_components, fit_tabnet_with_eval,
    SHAP_STRATEGY, CLASSIFIERS,
)

# Ensembl → Symbol map
BASE_DIR = Path(__file__).resolve().parent.parent.parent
SYMBOL_MAP_PATH = BASE_DIR / "processed" / "ensembl_symbol_map.csv"


def load_symbol_map():
    """Load ENSEMBL → gene symbol map."""
    if not SYMBOL_MAP_PATH.exists():
        return {}
    df = pd.read_csv(SYMBOL_MAP_PATH)
    return dict(zip(df["ENSEMBL"], df["SYMBOL"]))


def resolve_gene_names(gene_ids, sym_map):
    """Convert Ensembl IDs to gene symbols where possible."""
    return [sym_map.get(g, g) for g in gene_ids]

# Top-2 methods per tissue (from Phase 1 M1 results, gene features)
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

# Max genes for KernelExplainer (DD-23: avoid OOM on 17K features)
KERNEL_MAX_GENES = 2000
KERNEL_NSAMPLES = 500  # KernelExplainer nsamples (default=2*K+2048 is too many)
KERNEL_KMEANS_K = 100


def compute_shap_tree(model, X_train, X_test, gene_names):
    """TreeExplainer for RF / XGBoost (not in Pipeline — raw features)."""
    import shap

    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(X_test)

    # Handle different shap return formats:
    # - shap <0.40: list [class0_2d, class1_2d] for RF
    # - shap >=0.40: may return 3D array (n_samples, n_features, n_classes) for RF
    # - XGB: 2D array (n_samples, n_features), positive = class 1
    if isinstance(shap_values, list):
        sv = shap_values[1]  # class 1 = FLT
    elif shap_values.ndim == 3:
        sv = shap_values[:, :, 1]  # class 1 = FLT
    else:
        sv = shap_values  # XGB: 2D, positive = FLT

    # sv shape: (n_test, n_features)
    assert sv.ndim == 2, f"Expected 2D SHAP values, got {sv.ndim}D shape {sv.shape}"
    importance = np.mean(np.abs(sv), axis=0)  # mean |SHAP|
    direction = np.mean(sv, axis=0)  # mean signed SHAP

    return importance, direction, gene_names, None


def compute_shap_linear(model, X_train, X_test, gene_names):
    """Linear importance for ElasticNet-LR (Pipeline: scaler + clf)."""
    scaler = model.named_steps["scaler"]
    clf = model.named_steps["clf"]

    coef = clf.coef_.ravel()  # shape: (n_features,)
    scale = scaler.scale_  # std of each feature in original space

    # Importance in original feature space: |coef| × std
    # This gives expected |change in log-odds| per 1-SD change in original feature
    importance = np.abs(coef) * scale
    # Direction: positive coef × scale → gene higher in FLT increases P(FLT)
    direction = coef * scale

    return importance, direction, gene_names, None


def compute_shap_pca_linear(model, X_train, X_test, gene_names):
    """PCA-LR: back-project PCA-space coefficients to gene space."""
    scaler = model.named_steps["scaler"]
    pca = model.named_steps["pca"]
    clf = model.named_steps["clf"]

    lr_coef = clf.coef_.ravel()  # shape: (n_pca_components,)
    pca_components = pca.components_  # shape: (n_pca_components, n_features_scaled)
    scale = scaler.scale_  # std of each original feature

    # Back-project: gene importance in scaled space
    gene_coef_scaled = lr_coef @ pca_components  # shape: (n_features,)

    # Convert to original space
    importance = np.abs(gene_coef_scaled) * scale
    direction = gene_coef_scaled * scale

    # Also store explained variance ratio
    evr = float(np.sum(pca.explained_variance_ratio_))

    return importance, direction, gene_names, {"explained_variance_ratio": round(evr, 4)}


def compute_shap_kernel_full_pipeline(model, X_train, X_test, y_train,
                                       gene_names, method_key, seed=42,
                                       max_genes=KERNEL_MAX_GENES,
                                       n_samples=KERNEL_NSAMPLES):
    """KernelExplainer using re-trained model on gene subset.

    Retrains a fresh model on top-{max_genes} genes to ensure consistency.
    """
    import shap

    # Variance-based feature selection on training data
    n_features = X_train.shape[1]
    if n_features > max_genes:
        variances = np.var(X_train, axis=0)
        top_idx = np.argsort(variances)[-max_genes:]
        top_idx = np.sort(top_idx)
        X_train_sub = X_train[:, top_idx]
        X_test_sub = X_test[:, top_idx]
        gene_names_sub = [gene_names[i] for i in top_idx]
        gene_filter = f"top_{max_genes}_variance"
    else:
        X_train_sub = X_train
        X_test_sub = X_test
        gene_names_sub = list(gene_names)
        gene_filter = "all"

    # Build and train fresh model on subset
    _, fresh_model = get_classifier(method_key, seed=seed)

    if method_key == "tabnet":
        n_val = max(2, int(0.1 * len(y_train)))
        fit_tabnet_with_eval(
            fresh_model,
            X_train_sub[:-n_val].astype(np.float32),
            y_train[:-n_val].astype(np.int64),
            X_train_sub[-n_val:].astype(np.float32),
            y_train[-n_val:].astype(np.int64),
            max_epochs=200, patience=20,
        )
        predict_fn = lambda x: fresh_model.predict_proba(x.astype(np.float32))[:, 1]
        background_data = X_train_sub.astype(np.float32)
    else:
        fresh_model.fit(X_train_sub, y_train)
        predict_fn = lambda x: fresh_model.predict_proba(x)[:, 1]
        background_data = X_train_sub

    # Background dataset
    background = shap.kmeans(background_data, min(KERNEL_KMEANS_K, len(background_data)))

    explainer = shap.KernelExplainer(predict_fn, background)

    # Subsample test set if too large (>30 samples)
    if X_test_sub.shape[0] > 30:
        rng = np.random.default_rng(seed)
        test_idx = rng.choice(X_test_sub.shape[0], 30, replace=False)
        X_test_shap = X_test_sub[test_idx]
    else:
        X_test_shap = X_test_sub

    if method_key == "tabnet":
        X_test_shap = X_test_shap.astype(np.float32)

    sv = explainer.shap_values(X_test_shap, nsamples=n_samples, silent=True)

    importance = np.mean(np.abs(sv), axis=0)
    direction = np.mean(sv, axis=0)

    return importance, direction, gene_names_sub, {"gene_filter": gene_filter}


def run_shap_single(tissue, method_key, seed=42, verbose=True):
    """Compute SHAP importance for a single tissue × method combination.

    Returns result dict with gene importance rankings.
    """
    strategy = SHAP_STRATEGY[method_key]
    folds = get_folds(tissue)

    if verbose:
        print(f"\n{'='*60}")
        print(f"SHAP: {tissue} × {method_key} (strategy={strategy})")
        print(f"  Folds: {len(folds)}")

    fold_importances = []
    fold_directions = []
    fold_gene_names_list = []  # per-fold gene names (may differ for kernel methods)
    fold_extras = []

    for fi, fold in enumerate(folds):
        X_train = fold["train_X"]
        y_train = fold["train_y"]
        X_test = fold["test_X"]
        y_test = fold["test_y"]
        gene_names = fold["gene_names"]

        if len(np.unique(y_test)) < 2:
            if verbose:
                print(f"  Fold {fold['fold_name']}: single-class test, SKIP")
            continue

        if verbose:
            print(f"  Fold {fold['fold_name']}: n_train={len(y_train)}, n_test={len(y_test)}", end="")

        # Build and train model (same as multi_method_eval.py)
        _, model = get_classifier(method_key, seed=seed)
        adapt_pca_components(model, X_train.shape[0], X_train.shape[1])

        try:
            if method_key == "tabnet":
                n_val = max(2, int(0.1 * len(y_train)))
                fit_tabnet_with_eval(
                    model,
                    X_train[:-n_val].astype(np.float32),
                    y_train[:-n_val].astype(np.int64),
                    X_train[-n_val:].astype(np.float32),
                    y_train[-n_val:].astype(np.int64),
                    max_epochs=200, patience=20,
                )
            else:
                model.fit(X_train, y_train)
        except Exception as e:
            if verbose:
                print(f" — training FAILED: {e}")
            continue

        # Compute SHAP based on strategy
        try:
            if strategy == "tree":
                imp, dirn, gn, extra = compute_shap_tree(model, X_train, X_test, gene_names)
            elif strategy == "linear":
                imp, dirn, gn, extra = compute_shap_linear(model, X_train, X_test, gene_names)
            elif strategy == "pca_linear":
                imp, dirn, gn, extra = compute_shap_pca_linear(model, X_train, X_test, gene_names)
            elif strategy == "kernel":
                imp, dirn, gn, extra = compute_shap_kernel_full_pipeline(
                    model, X_train, X_test, y_train, gene_names,
                    method_key, seed=seed,
                )
            else:
                raise ValueError(f"Unknown SHAP strategy: {strategy}")
        except Exception as e:
            if verbose:
                print(f" — SHAP FAILED: {e}")
            continue

        fold_importances.append(imp)
        fold_directions.append(dirn)
        fold_gene_names_list.append(gn)
        fold_extras.append(extra)

        if verbose:
            top3 = np.argsort(imp)[-3:][::-1]
            top3_str = ", ".join(f"{gn[i]}={imp[i]:.4f}" for i in top3)
            print(f" — top3: {top3_str}")

    if not fold_importances:
        return {"error": f"No valid folds for {tissue}/{method_key}"}

    # Load symbol map for Ensembl → gene symbol conversion
    sym_map = load_symbol_map()

    # Aggregate across folds
    # For tree/linear/pca_linear: all folds have identical gene sets
    # For kernel: folds may select different top-K genes by variance
    # → align genes across folds using union + zero-fill
    all_same = all(list(gn) == list(fold_gene_names_list[0])
                   for gn in fold_gene_names_list)

    if all_same:
        # Fast path: same genes in all folds (tree/linear/pca_linear)
        mean_importance = np.mean(fold_importances, axis=0)
        mean_direction = np.mean(fold_directions, axis=0)
        final_gene_names = list(fold_gene_names_list[0])
    else:
        # Kernel methods: align by gene name across folds
        # Build union of all gene names preserving order
        seen = set()
        final_gene_names = []
        for gn_list in fold_gene_names_list:
            for g in gn_list:
                if g not in seen:
                    seen.add(g)
                    final_gene_names.append(g)

        gene_to_idx = {g: i for i, g in enumerate(final_gene_names)}
        n_genes_union = len(final_gene_names)

        # Accumulate importance/direction with per-gene fold counts
        sum_importance = np.zeros(n_genes_union)
        sum_direction = np.zeros(n_genes_union)
        count = np.zeros(n_genes_union)

        for imp, dirn, gn_list in zip(fold_importances, fold_directions,
                                       fold_gene_names_list):
            for j, g in enumerate(gn_list):
                idx = gene_to_idx[g]
                sum_importance[idx] += imp[j]
                sum_direction[idx] += dirn[j]
                count[idx] += 1

        count = np.maximum(count, 1)  # avoid division by zero
        mean_importance = sum_importance / count
        mean_direction = sum_direction / count

        if verbose:
            print(f"  Gene alignment: {n_genes_union} union genes "
                  f"from {len(fold_importances)} folds")

    # Resolve gene names (Ensembl → Symbol)
    gene_symbols = resolve_gene_names(final_gene_names, sym_map)

    # Rank genes
    rank_order = np.argsort(mean_importance)[::-1]  # descending importance
    top_100 = [gene_symbols[i] for i in rank_order[:100]]

    # Full gene importance dict (sorted by importance, using symbols)
    gene_importance = {}
    gene_direction = {}
    for i in rank_order:
        gene_importance[gene_symbols[i]] = round(float(mean_importance[i]), 6)
        gene_direction[gene_symbols[i]] = round(float(mean_direction[i]), 6)

    n_mapped = sum(1 for s, e in zip(gene_symbols, final_gene_names) if s != e)

    result = {
        "tissue": tissue,
        "method": method_key,
        "method_label": CLASSIFIERS[method_key][0],
        "shap_strategy": strategy,
        "n_genes": len(final_gene_names),
        "n_symbols_mapped": n_mapped,
        "n_folds_used": len(fold_importances),
        "n_folds_total": len(folds),
        "top_100_genes": top_100,
        "gene_importance": gene_importance,
        "gene_direction": gene_direction,
        "seed": seed,
        "timestamp": datetime.now().isoformat(),
    }

    # Add extras from first fold (gene_filter, explained_variance_ratio, etc.)
    for extra in fold_extras:
        if extra:
            result.update(extra)

    return result


def main():
    parser = argparse.ArgumentParser(description="SHAP multi-method interpretability")
    parser.add_argument("--tissue", type=str, default=None,
                        help="Single tissue (default: all top-2 combinations)")
    parser.add_argument("--method", type=str, default=None,
                        help="Single method (default: top-2 for tissue)")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    V4_EVAL_DIR.mkdir(parents=True, exist_ok=True)

    # Determine which combinations to run
    if args.tissue and args.method:
        combos = [(args.tissue, args.method)]
    elif args.tissue:
        combos = [(args.tissue, m) for m in TOP2_METHODS.get(args.tissue, [])]
    else:
        combos = []
        for tissue, methods in TOP2_METHODS.items():
            for method in methods:
                combos.append((tissue, method))

    print(f"SHAP Phase 4: {len(combos)} tissue × method combinations")

    for tissue, method in combos:
        output_path = V4_EVAL_DIR / f"SHAP_{tissue}_{method}.json"

        if output_path.exists():
            print(f"\n  SKIP {tissue}/{method}: {output_path} exists")
            continue

        result = run_shap_single(tissue, method, seed=args.seed)

        with open(output_path, 'w') as f:
            json.dump(result, f, indent=2)

        if "error" not in result:
            print(f"  → {output_path.name}: {result['n_genes']} genes, "
                  f"{result['n_folds_used']}/{result['n_folds_total']} folds")
        else:
            print(f"  → {output_path.name}: ERROR — {result['error']}")

    print(f"\nDone. Outputs in {V4_EVAL_DIR}/SHAP_*.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())

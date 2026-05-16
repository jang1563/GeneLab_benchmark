#!/usr/bin/env python3
"""
ablation_pca_components.py — GeneLabBench v4: PCA Component Ablation

Optimal PCA dimensionality for spaceflight classification.

Design:
- n_components: [5, 10, 20, 50, 100, 200]
- PCA-LR only (other methods don't use PCA)
- 8 tissues
- NO variance filter (matches Phase 1)
- n_components=50 is Phase 1 default → must match M1 exactly (regression test)
- Records explained_variance_ratio_sum per fold
- Scope: 6 levels × 8 tissues = 48 evaluations

Usage:
  python ablation_pca_components.py --tissue liver --n-components 20
"""

import json
import time
import argparse
import warnings
import traceback
import numpy as np
from pathlib import Path
from datetime import datetime
from sklearn.metrics import roc_auc_score
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

warnings.filterwarnings("ignore")

import sys
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from v4_utils import (
    TISSUE_MISSIONS, LOMO_TISSUES,
    V4_EVAL_DIR,
    bootstrap_auroc, permutation_pvalue,
    get_folds,
)

NC_VALUES = [5, 10, 20, 50, 100, 200]


def evaluate_pca_ablation(tissue, n_components,
                          n_bootstrap=2000, n_perm=1000,
                          seed=42, verbose=True):
    """Run PCA-LR with specific n_components."""
    t_start = time.time()

    if verbose:
        print(f"\n{'='*60}")
        print(f"PCA Ablation | Tissue: {tissue} | n_components={n_components}")
        print(f"{'='*60}")

    folds = get_folds(tissue)
    fold_results = []

    for fold in folds:
        X_train = fold["train_X"]
        y_train = fold["train_y"]
        X_test = fold["test_X"]
        y_test = fold["test_y"]
        test_mission = fold["test_mission"]

        if len(np.unique(y_test)) < 2:
            continue

        # Scale
        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)

        # PCA with clamping
        nc_actual = min(n_components, X_train_s.shape[0] - 1, X_train_s.shape[1])
        pca = PCA(n_components=nc_actual, random_state=seed)
        X_train_p = pca.fit_transform(X_train_s)
        X_test_p = pca.transform(X_test_s)

        explained_var_sum = float(np.sum(pca.explained_variance_ratio_))

        # Logistic Regression (matches Phase 1 PCA-LR pipeline)
        clf = LogisticRegression(
            C=1.0, class_weight="balanced",
            max_iter=5000, random_state=seed,
        )
        t0 = time.time()
        clf.fit(X_train_p, y_train)
        train_time = time.time() - t0

        y_score = clf.predict_proba(X_test_p)[:, 1]
        auroc = float(roc_auc_score(y_test, y_score))
        mean_bs, lower_bs, upper_bs = bootstrap_auroc(y_test, y_score, n_bootstrap)
        pval = permutation_pvalue(y_test, y_score, n_perm)

        result = {
            "test_mission": test_mission,
            "fold_name": fold["fold_name"],
            "auroc": round(auroc, 4),
            "bootstrap_mean": round(mean_bs, 4),
            "ci_lower": round(lower_bs, 4),
            "ci_upper": round(upper_bs, 4),
            "perm_pvalue": round(pval, 4),
            "n_train": fold["n_train"],
            "n_test": fold["n_test"],
            "n_components_requested": n_components,
            "n_components_actual": nc_actual,
            "explained_variance_ratio_sum": round(explained_var_sum, 4),
            "train_time_s": round(train_time, 2),
        }
        fold_results.append(result)

        if verbose:
            print(f"  {test_mission}: AUROC={auroc:.3f} "
                  f"(nc={nc_actual}, var={explained_var_sum:.3f})")

    if not fold_results:
        return None

    aurocs = [r["auroc"] for r in fold_results]
    total_time = time.time() - t_start

    summary = {
        "experiment": "ablation_pca",
        "model": "PCA-LR",
        "model_key": "pca_lr",
        "tissue": tissue,
        "feature_type": "gene",
        "ablation_param": "n_components",
        "ablation_value": n_components,
        "cv_strategy": "LOMO" if tissue in LOMO_TISSUES else "5-fold_stratified",
        "n_folds": len(fold_results),
        "mean_auroc": round(float(np.mean(aurocs)), 4),
        "std_auroc": round(float(np.std(aurocs)), 4),
        "min_auroc": round(float(np.min(aurocs)), 4),
        "max_auroc": round(float(np.max(aurocs)), 4),
        "mean_explained_var": round(float(np.mean(
            [r["explained_variance_ratio_sum"] for r in fold_results])), 4),
        "mean_perm_pvalue": round(float(np.mean(
            [r["perm_pvalue"] for r in fold_results])), 4),
        "n_bootstrap": n_bootstrap,
        "n_perm": n_perm,
        "seed": seed,
        "wall_time_sec": round(total_time, 1),
        "timestamp": datetime.now().isoformat(),
        "folds": fold_results,
    }

    if verbose:
        print(f"\n  Mean AUROC = {summary['mean_auroc']:.4f} ± {summary['std_auroc']:.4f}")
        print(f"  Mean explained variance = {summary['mean_explained_var']:.4f}")

    return summary


def main():
    parser = argparse.ArgumentParser(description="PCA Component Ablation")
    parser.add_argument("--tissue", required=True, choices=list(TISSUE_MISSIONS.keys()))
    parser.add_argument("--n-components", required=True, type=int,
                        choices=NC_VALUES)
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=str, default=None)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    verbose = not args.quiet
    output_dir = Path(args.output_dir) if args.output_dir else V4_EVAL_DIR

    try:
        result = evaluate_pca_ablation(
            tissue=args.tissue,
            n_components=args.n_components,
            n_bootstrap=args.n_bootstrap,
            n_perm=args.n_perm,
            seed=args.seed,
            verbose=verbose,
        )
        if result:
            fname = f"ABL_pca_{args.tissue}_nc{args.n_components}.json"
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / fname).write_text(json.dumps(result, indent=2))
            if verbose:
                print(f"  Saved: {output_dir / fname}")
    except Exception as e:
        print(f"[ERROR] {args.tissue}/nc{args.n_components}: {e}")
        traceback.print_exc()


if __name__ == "__main__":
    main()

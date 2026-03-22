#!/usr/bin/env python3
"""
ablation_feature_count.py — GeneLabBench v4: Feature Count Ablation

How many genes are needed? AUROC vs top-K genes by variance.

Design:
- K values: [100, 500, 1000, 2000, 5000, 10000, all]
- Top-2 methods: PCA-LR, ElasticNet-LR
- 8 tissues
- Variance ranking on TRAINING data only (per fold, no leakage)
- NO baseline variance filter (matches Phase 1)
- K=all uses ALL genes (identical to Phase 1)
- Scope: 7 levels × 2 methods × 8 tissues = 112 evaluations

Usage:
  python ablation_feature_count.py --tissue liver --method pca_lr --top-k 1000
  python ablation_feature_count.py --tissue liver --method pca_lr --top-k all
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
from classifier_registry import (
    CLASSIFIERS,
    get_classifier, adapt_pca_components, fit_tabnet_with_eval,
)

TOP_K_VALUES = [100, 500, 1000, 2000, 5000, 10000, "all"]


def evaluate_feature_ablation(tissue, method_key, top_k,
                              n_bootstrap=2000, n_perm=1000,
                              seed=42, verbose=True):
    """Run evaluation with top-K genes by variance."""
    t_start = time.time()
    label, _ = get_classifier(method_key, seed)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Feature Ablation | Tissue: {tissue} | Method: {label} | K={top_k}")
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

        # Feature selection: top-K by TRAINING variance (no leakage)
        if top_k == "all":
            X_train_f, X_test_f = X_train, X_test
            top_k_actual = X_train.shape[1]
        else:
            var = np.var(X_train, axis=0)
            top_k_actual = min(int(top_k), X_train.shape[1])
            top_indices = np.argsort(var)[-top_k_actual:]
            X_train_f = X_train[:, top_indices]
            X_test_f = X_test[:, top_indices]

        _, model = get_classifier(method_key, seed)
        adapt_pca_components(model, X_train_f.shape[0], X_train_f.shape[1])

        t0 = time.time()
        try:
            if method_key == "tabnet":
                n_val = max(2, int(0.1 * len(y_train)))
                fit_tabnet_with_eval(
                    model, X_train_f[:-n_val], y_train[:-n_val].astype(np.int64),
                    X_train_f[-n_val:], y_train[-n_val:].astype(np.int64),
                    max_epochs=200, patience=20
                )
            else:
                model.fit(X_train_f, y_train)
        except Exception as e:
            if verbose:
                print(f"  [ERROR] {test_mission}: {e}")
            continue
        train_time = time.time() - t0

        if hasattr(model, "predict_proba"):
            y_score = model.predict_proba(X_test_f)[:, 1]
        else:
            y_score = model.decision_function(X_test_f)

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
            "n_features_used": top_k_actual,
            "train_time_s": round(train_time, 2),
        }
        fold_results.append(result)

        if verbose:
            print(f"  {test_mission}: AUROC={auroc:.3f} (K={top_k_actual})")

    if not fold_results:
        return None

    aurocs = [r["auroc"] for r in fold_results]
    total_time = time.time() - t_start

    summary = {
        "experiment": "ablation_feature_count",
        "model": label,
        "model_key": method_key,
        "tissue": tissue,
        "feature_type": "gene",
        "ablation_param": "top_k",
        "ablation_value": top_k if top_k == "all" else int(top_k),
        "top_k_actual": int(fold_results[0]["n_features_used"]),
        "cv_strategy": "LOMO" if tissue in LOMO_TISSUES else "5-fold_stratified",
        "n_folds": len(fold_results),
        "mean_auroc": round(float(np.mean(aurocs)), 4),
        "std_auroc": round(float(np.std(aurocs)), 4),
        "min_auroc": round(float(np.min(aurocs)), 4),
        "max_auroc": round(float(np.max(aurocs)), 4),
        "mean_perm_pvalue": round(float(np.mean([r["perm_pvalue"] for r in fold_results])), 4),
        "n_bootstrap": n_bootstrap,
        "n_perm": n_perm,
        "seed": seed,
        "wall_time_sec": round(total_time, 1),
        "timestamp": datetime.now().isoformat(),
        "folds": fold_results,
    }

    if verbose:
        print(f"\n  Mean AUROC = {summary['mean_auroc']:.4f} ± {summary['std_auroc']:.4f}")

    return summary


def main():
    parser = argparse.ArgumentParser(description="Feature Count Ablation")
    parser.add_argument("--tissue", required=True, choices=list(TISSUE_MISSIONS.keys()))
    parser.add_argument("--method", required=True, choices=["pca_lr", "elasticnet_lr"])
    parser.add_argument("--top-k", required=True,
                        help="Number of top-variance genes (or 'all')")
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=str, default=None)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    verbose = not args.quiet
    output_dir = Path(args.output_dir) if args.output_dir else V4_EVAL_DIR

    top_k = args.top_k if args.top_k == "all" else int(args.top_k)

    try:
        result = evaluate_feature_ablation(
            tissue=args.tissue,
            method_key=args.method,
            top_k=top_k,
            n_bootstrap=args.n_bootstrap,
            n_perm=args.n_perm,
            seed=args.seed,
            verbose=verbose,
        )
        if result:
            k_str = top_k if top_k == "all" else f"{int(top_k)}"
            fname = f"ABL_feat_{args.tissue}_{args.method}_k{k_str}.json"
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / fname).write_text(json.dumps(result, indent=2))
            if verbose:
                print(f"  Saved: {output_dir / fname}")
    except Exception as e:
        print(f"[ERROR] {args.tissue}/{args.method}/k{top_k}: {e}")
        traceback.print_exc()


if __name__ == "__main__":
    main()

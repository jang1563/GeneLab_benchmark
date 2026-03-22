#!/usr/bin/env python3
"""
nc_random_features.py — GeneLabBench v4: NC2 Random Features Negative Control

Confirms classifiers don't find spurious patterns in noise.
Real labels but random Gaussian features matching original dimensionality.
Expected AUROC ≈ 0.5.

Design:
- Real labels (not shuffled)
- X = random Gaussian noise, shape matched to original gene features
- Same fold structure as Phase 1 (via get_folds() with real labels)
- Gene-dimensionality matching only
- Scope: 8 tissues × 8 methods = 64 evaluations

Usage:
  python nc_random_features.py --tissue liver --method pca_lr
  python nc_random_features.py --tissue liver --method all
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
    TISSUE_MISSIONS, LOMO_TISSUES, KFOLD_TISSUES,
    BASE_DIR, V4_EVAL_DIR,
    load_metadata, load_gene_features, align_features_with_meta,
    encode_labels, bootstrap_auroc, permutation_pvalue,
    get_folds,
)
from classifier_registry import (
    CLASSIFIERS, CPU_METHODS, GPU_METHODS,
    get_classifier, adapt_pca_components, fit_tabnet_with_eval,
)


def evaluate_nc2(tissue, method_key, n_bootstrap=2000, n_perm=1000,
                 seed=42, verbose=True):
    """Run NC2 evaluation: classifier on random Gaussian features with real labels."""
    t_start = time.time()
    label, _ = get_classifier(method_key, seed)

    if verbose:
        print(f"\n{'='*60}")
        print(f"NC2 Random Features | Tissue: {tissue} | Method: {label}")
        print(f"{'='*60}")

    # Get real folds (real labels, real fold structure)
    folds = get_folds(tissue)

    if not folds:
        print(f"  [ERROR] No folds for {tissue}")
        return None

    n_features = folds[0]["train_X"].shape[1]
    if verbose:
        print(f"  {len(folds)} folds, replacing {n_features} gene features with random noise")

    rng = np.random.default_rng(seed)
    fold_results = []

    for fold in folds:
        test_mission = fold["test_mission"]
        y_train = fold["train_y"]
        y_test = fold["test_y"]

        if len(np.unique(y_test)) < 2:
            continue

        # Replace features with random Gaussian noise
        X_train = rng.standard_normal((fold["n_train"], n_features)).astype(np.float32)
        X_test = rng.standard_normal((fold["n_test"], n_features)).astype(np.float32)

        _, model = get_classifier(method_key, seed)
        adapt_pca_components(model, X_train.shape[0], X_train.shape[1])

        t0 = time.time()
        try:
            if method_key == "tabnet":
                n_val = max(2, int(0.1 * len(y_train)))
                fit_tabnet_with_eval(
                    model, X_train[:-n_val], y_train[:-n_val].astype(np.int64),
                    X_train[-n_val:], y_train[-n_val:].astype(np.int64),
                    max_epochs=200, patience=20
                )
            else:
                model.fit(X_train, y_train)
        except Exception as e:
            if verbose:
                print(f"  [ERROR] {test_mission}: {e}")
            continue
        train_time = time.time() - t0

        if hasattr(model, "predict_proba"):
            y_score = model.predict_proba(X_test)[:, 1]
        else:
            y_score = model.decision_function(X_test)

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
            "n_features": n_features,
            "train_time_s": round(train_time, 2),
        }
        fold_results.append(result)

        if verbose:
            print(f"  {test_mission}: AUROC={auroc:.3f} p={pval:.3f}")

    if not fold_results:
        return None

    aurocs = [r["auroc"] for r in fold_results]
    total_time = time.time() - t_start

    summary = {
        "experiment": "NC2_random",
        "model": label,
        "model_key": method_key,
        "tissue": tissue,
        "feature_type": "random_gaussian",
        "n_original_features": n_features,
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
        verdict = "PASS" if 0.35 <= summary["mean_auroc"] <= 0.65 else "WARN"
        print(f"  Verdict: {verdict}")

    return summary


def main():
    parser = argparse.ArgumentParser(description="NC2: Random Features Negative Control")
    parser.add_argument("--tissue", required=True, choices=list(TISSUE_MISSIONS.keys()))
    parser.add_argument("--method", required=True,
                        help="Classifier key (or 'all' / 'cpu' / 'gpu')")
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=str, default=None)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    verbose = not args.quiet
    output_dir = Path(args.output_dir) if args.output_dir else V4_EVAL_DIR

    if args.method == "all":
        methods = list(CLASSIFIERS.keys())
    elif args.method == "cpu":
        methods = CPU_METHODS
    elif args.method == "gpu":
        methods = GPU_METHODS
    else:
        methods = [args.method]

    for method_key in methods:
        if method_key not in CLASSIFIERS:
            print(f"[WARN] Unknown method: {method_key}")
            continue

        try:
            result = evaluate_nc2(
                tissue=args.tissue,
                method_key=method_key,
                n_bootstrap=args.n_bootstrap,
                n_perm=args.n_perm,
                seed=args.seed,
                verbose=verbose,
            )
            if result:
                fname = f"NC2_random_{args.tissue}_{method_key}.json"
                output_dir.mkdir(parents=True, exist_ok=True)
                (output_dir / fname).write_text(json.dumps(result, indent=2))
                if verbose:
                    print(f"  Saved: {output_dir / fname}")
        except Exception as e:
            print(f"[ERROR] {args.tissue}/{method_key}: {e}")
            traceback.print_exc()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
ablation_sample_size.py — GeneLabBench v4: Sample Size Ablation

Learning curve via subsampled training data.

Design:
- Fractions: [0.2, 0.4, 0.6, 0.8, 1.0]
- Top-2 methods: PCA-LR, ElasticNet-LR
- 4 largest LOMO tissues: liver (~261), thymus (~80), skin (~65), gastro (~50)
- 10 independent repeats per fraction (stratified subsampling)
- LOMO fold structure preserved: subsample WITHIN each fold's training set
- fraction=1.0: bypass subsampling entirely (regression test vs M1)
- NO variance filter (matches Phase 1)
- Scope: 5 fractions × 2 methods × 4 tissues × 10 repeats = 400 evaluations

Usage:
  python ablation_sample_size.py --tissue liver --method pca_lr --fraction 0.6 --repeat 0
  python ablation_sample_size.py --tissue liver --method pca_lr --fraction 1.0 --repeat 0
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
from sklearn.model_selection import StratifiedShuffleSplit

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

SAMPLE_TISSUES = ["liver", "thymus", "skin", "gastrocnemius"]
FRACTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]


def evaluate_sample_ablation(tissue, method_key, fraction, repeat,
                             n_bootstrap=2000, n_perm=1000,
                             seed=42, verbose=True):
    """Run evaluation with subsampled training data."""
    t_start = time.time()
    label, _ = get_classifier(method_key, seed)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Sample Ablation | Tissue: {tissue} | Method: {label}")
        print(f"Fraction: {fraction} | Repeat: {repeat}")
        print(f"{'='*60}")

    folds = get_folds(tissue)
    fold_results = []

    for fold_idx, fold in enumerate(folds):
        X_train = fold["train_X"]
        y_train = fold["train_y"]
        X_test = fold["test_X"]
        y_test = fold["test_y"]
        test_mission = fold["test_mission"]

        if len(np.unique(y_test)) < 2:
            continue

        # Subsample training data
        if fraction < 1.0:
            n_sub = max(int(len(y_train) * fraction), 2)
            subsample_seed = seed + repeat * 1000 + fold_idx

            # Check if stratified split is possible
            from collections import Counter
            class_counts = Counter(y_train)
            min_class = min(class_counts.values())

            if min_class < 2:
                # Can't stratify, use random subsample
                rng = np.random.default_rng(subsample_seed)
                sub_idx = rng.choice(len(y_train), size=n_sub, replace=False)
            else:
                try:
                    sss = StratifiedShuffleSplit(
                        n_splits=1, train_size=n_sub,
                        random_state=subsample_seed
                    )
                    sub_idx, _ = next(sss.split(X_train, y_train))
                except ValueError:
                    # Fallback if stratification fails
                    rng = np.random.default_rng(subsample_seed)
                    sub_idx = rng.choice(len(y_train), size=n_sub, replace=False)

            X_train_sub = X_train[sub_idx]
            y_train_sub = y_train[sub_idx]

            # Skip if only one class after subsampling
            if len(np.unique(y_train_sub)) < 2:
                if verbose:
                    print(f"  [SKIP] {test_mission}: single class after subsampling")
                continue
        else:
            # fraction=1.0: bypass subsampling entirely (exact Phase 1 match)
            X_train_sub = X_train
            y_train_sub = y_train
            n_sub = len(y_train)

        _, model = get_classifier(method_key, seed)
        adapt_pca_components(model, X_train_sub.shape[0], X_train_sub.shape[1])

        t0 = time.time()
        try:
            if method_key == "tabnet":
                n_val = max(2, int(0.1 * len(y_train_sub)))
                fit_tabnet_with_eval(
                    model, X_train_sub[:-n_val],
                    y_train_sub[:-n_val].astype(np.int64),
                    X_train_sub[-n_val:],
                    y_train_sub[-n_val:].astype(np.int64),
                    max_epochs=200, patience=20
                )
            else:
                model.fit(X_train_sub, y_train_sub)
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
            "n_train_full": fold["n_train"],
            "n_train_used": n_sub,
            "n_test": fold["n_test"],
            "train_time_s": round(train_time, 2),
        }
        fold_results.append(result)

        if verbose:
            print(f"  {test_mission}: AUROC={auroc:.3f} "
                  f"(train: {n_sub}/{fold['n_train']})")

    if not fold_results:
        return None

    aurocs = [r["auroc"] for r in fold_results]
    total_time = time.time() - t_start

    summary = {
        "experiment": "ablation_sample_size",
        "model": label,
        "model_key": method_key,
        "tissue": tissue,
        "feature_type": "gene",
        "ablation_param": "fraction",
        "ablation_value": fraction,
        "repeat": repeat,
        "cv_strategy": "LOMO" if tissue in LOMO_TISSUES else "5-fold_stratified",
        "n_folds": len(fold_results),
        "mean_auroc": round(float(np.mean(aurocs)), 4),
        "std_auroc": round(float(np.std(aurocs)), 4),
        "min_auroc": round(float(np.min(aurocs)), 4),
        "max_auroc": round(float(np.max(aurocs)), 4),
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

    return summary


def main():
    parser = argparse.ArgumentParser(description="Sample Size Ablation")
    parser.add_argument("--tissue", required=True, choices=SAMPLE_TISSUES)
    parser.add_argument("--method", required=True, choices=["pca_lr", "elasticnet_lr"])
    parser.add_argument("--fraction", required=True, type=float,
                        choices=FRACTIONS)
    parser.add_argument("--repeat", required=True, type=int,
                        help="Repeat index (0-9)")
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=str, default=None)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    verbose = not args.quiet
    output_dir = Path(args.output_dir) if args.output_dir else V4_EVAL_DIR

    try:
        result = evaluate_sample_ablation(
            tissue=args.tissue,
            method_key=args.method,
            fraction=args.fraction,
            repeat=args.repeat,
            n_bootstrap=args.n_bootstrap,
            n_perm=args.n_perm,
            seed=args.seed,
            verbose=verbose,
        )
        if result:
            frac_str = f"{args.fraction:.1f}".replace(".", "")
            fname = f"ABL_sample_{args.tissue}_{args.method}_f{frac_str}_r{args.repeat}.json"
            output_dir.mkdir(parents=True, exist_ok=True)
            (output_dir / fname).write_text(json.dumps(result, indent=2))
            if verbose:
                print(f"  Saved: {output_dir / fname}")
    except Exception as e:
        print(f"[ERROR] {args.tissue}/{args.method}/f{args.fraction}/r{args.repeat}: {e}")
        traceback.print_exc()


if __name__ == "__main__":
    main()

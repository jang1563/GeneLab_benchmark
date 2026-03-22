#!/usr/bin/env python3
"""
nc_shuffled_labels.py — GeneLabBench v4: NC1 Shuffled Labels Negative Control

Confirms classifiers don't exploit data structure artifacts.
Shuffles flight/ground labels → expected AUROC ≈ 0.5.

Design:
- Identical preprocessing to Phase 1 (no variance filter, same fold structure)
- Shuffles y at 0/1 level AFTER encode_labels() filtering
- Manually builds LOMO/kfold folds using shuffled y (cannot reuse get_folds()
  since it hardcodes real labels internally)
- Gene features only
- Scope: 8 tissues × 8 methods = 64 evaluations

Usage:
  python nc_shuffled_labels.py --tissue liver --method pca_lr
  python nc_shuffled_labels.py --tissue liver --method all
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
)
from classifier_registry import (
    CLASSIFIERS, CPU_METHODS, GPU_METHODS,
    get_classifier, adapt_pca_components, fit_tabnet_with_eval,
)


def build_shuffled_lomo_folds(tissue, seed=42):
    """Build LOMO folds with shuffled labels.

    Cannot use get_folds() because it calls encode_labels() internally
    with real labels. Instead: load → encode → filter → shuffle → split.
    """
    meta = load_metadata(tissue)
    genes = load_gene_features(tissue)
    genes, meta = align_features_with_meta(genes, meta)

    y, valid_mask = encode_labels(meta)
    genes = genes[valid_mask]
    meta = meta[valid_mask]
    y = y[valid_mask].values.astype(int)
    X = genes.values.astype(np.float32)

    # Shuffle labels (preserves global class balance)
    rng = np.random.default_rng(seed)
    y_shuffled = rng.permutation(y)

    folds = []
    for test_mission in TISSUE_MISSIONS[tissue]:
        test_mask = (meta["mission"] == test_mission).values
        train_mask = ~test_mask

        if test_mask.sum() < 2 or train_mask.sum() < 5:
            continue

        y_test = y_shuffled[test_mask]
        # Skip if test set has only one class after shuffle
        if len(np.unique(y_test)) < 2:
            continue

        folds.append({
            "train_X": X[train_mask],
            "train_y": y_shuffled[train_mask],
            "test_X": X[test_mask],
            "test_y": y_test,
            "test_mission": test_mission,
            "fold_name": f"fold_{test_mission}_test",
            "n_train": int(train_mask.sum()),
            "n_test": int(test_mask.sum()),
        })

    return folds


def build_shuffled_kfold_folds(tissue, n_splits=5, seed=42):
    """Build stratified k-fold folds with shuffled labels."""
    from sklearn.model_selection import StratifiedKFold

    meta = load_metadata(tissue)
    genes = load_gene_features(tissue)
    genes, meta = align_features_with_meta(genes, meta)

    y, valid_mask = encode_labels(meta)
    genes = genes[valid_mask]
    y = y[valid_mask].values.astype(int)
    X = genes.values.astype(np.float32)

    # Shuffle labels
    rng = np.random.default_rng(seed)
    y_shuffled = rng.permutation(y)

    # Use REAL y for stratification (to match Phase 1 fold structure),
    # but train/test with SHUFFLED y
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    folds = []

    for i, (train_idx, test_idx) in enumerate(skf.split(X, y)):
        y_test = y_shuffled[test_idx]
        if len(np.unique(y_test)) < 2:
            continue

        folds.append({
            "train_X": X[train_idx],
            "train_y": y_shuffled[train_idx],
            "test_X": X[test_idx],
            "test_y": y_test,
            "test_mission": f"fold_{i+1}",
            "fold_name": f"kfold_{i+1}",
            "n_train": len(train_idx),
            "n_test": len(test_idx),
        })

    return folds


def evaluate_nc1(tissue, method_key, n_bootstrap=2000, n_perm=1000,
                 seed=42, verbose=True):
    """Run NC1 evaluation: one method on one tissue with shuffled labels."""
    t_start = time.time()
    label, _ = get_classifier(method_key, seed)

    if verbose:
        print(f"\n{'='*60}")
        print(f"NC1 Shuffled Labels | Tissue: {tissue} | Method: {label}")
        print(f"{'='*60}")

    # Build folds with shuffled labels
    if tissue in LOMO_TISSUES:
        folds = build_shuffled_lomo_folds(tissue, seed)
    else:
        folds = build_shuffled_kfold_folds(tissue, seed=seed)

    if verbose:
        print(f"  {len(folds)} folds (some may be skipped due to single-class after shuffle)")

    fold_results = []

    for fold in folds:
        X_train = fold["train_X"]
        y_train = fold["train_y"]
        X_test = fold["test_X"]
        y_test = fold["test_y"]
        test_mission = fold["test_mission"]

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
            "train_time_s": round(train_time, 2),
        }
        fold_results.append(result)

        if verbose:
            print(f"  {test_mission}: AUROC={auroc:.3f} p={pval:.3f}")

    if not fold_results:
        if verbose:
            print(f"  [ERROR] No valid folds")
        return None

    aurocs = [r["auroc"] for r in fold_results]
    total_time = time.time() - t_start

    summary = {
        "experiment": "NC1_shuffled",
        "model": label,
        "model_key": method_key,
        "tissue": tissue,
        "feature_type": "gene",
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
    parser = argparse.ArgumentParser(description="NC1: Shuffled Labels Negative Control")
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
            result = evaluate_nc1(
                tissue=args.tissue,
                method_key=method_key,
                n_bootstrap=args.n_bootstrap,
                n_perm=args.n_perm,
                seed=args.seed,
                verbose=verbose,
            )
            if result:
                fname = f"NC1_shuffled_{args.tissue}_{method_key}.json"
                output_dir.mkdir(parents=True, exist_ok=True)
                (output_dir / fname).write_text(json.dumps(result, indent=2))
                if verbose:
                    print(f"  Saved: {output_dir / fname}")
        except Exception as e:
            print(f"[ERROR] {args.tissue}/{method_key}: {e}")
            traceback.print_exc()


if __name__ == "__main__":
    main()

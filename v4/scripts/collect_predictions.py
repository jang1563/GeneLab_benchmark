#!/usr/bin/env python3
"""
collect_predictions.py — GeneLabBench v4 Phase 3: Collect raw predictions for DeLong test

Re-runs all classifiers on gene features to save (y_true, y_score) per fold.
Exact same preprocessing chain as multi_method_eval.py (lines 150-218).
Verifies AUROC matches M1 JSONs (±0.001).

Usage:
  python collect_predictions.py                           # All 64 tissue×method combos
  python collect_predictions.py --tissue liver            # All methods for liver
  python collect_predictions.py --tissue liver --method pca_lr  # Single combo (debug)
"""

import json
import sys
import time
import argparse
import traceback
import numpy as np
from pathlib import Path
from sklearn.metrics import roc_auc_score

# Add parent dirs to path
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "scripts"))

from v4_utils import (
    TISSUE_MISSIONS, V4_EVAL_DIR, get_folds
)

ALL_TISSUES = list(TISSUE_MISSIONS.keys())
from classifier_registry import (
    CLASSIFIERS, get_classifier, adapt_pca_components, fit_tabnet_with_eval
)


def collect_single(tissue, method_key, seed=42, verbose=True):
    """Collect (y_true, y_score) for one tissue × method.

    Mirrors multi_method_eval.py evaluate_single() preprocessing chain exactly.
    Returns list of fold dicts with y_true, y_score, auroc.
    """
    folds = get_folds(tissue)

    fold_predictions = []
    for fold in folds:
        X_train = fold["train_X"]
        y_train = fold["train_y"]
        X_test = fold["test_X"]
        y_test = fold["test_y"]
        fold_name = fold["fold_name"]

        # Skip single-class test folds (same as Phase 1)
        if len(np.unique(y_test)) < 2:
            if verbose:
                print(f"  [SKIP] {fold_name}: only one class in test set")
            continue

        # Build fresh Pipeline (includes StandardScaler + optional PCA)
        _, model = get_classifier(method_key, seed)

        # Adapt PCA components for small training sets
        adapt_pca_components(model, X_train.shape[0], X_train.shape[1])

        # Train — TabNet requires special handling
        t0 = time.time()
        try:
            if method_key == "tabnet":
                n_val = max(2, int(0.1 * len(y_train)))
                fit_tabnet_with_eval(
                    model,
                    X_train[:-n_val], y_train[:-n_val].astype(np.int64),
                    X_train[-n_val:], y_train[-n_val:].astype(np.int64),
                    max_epochs=200, patience=20
                )
            else:
                model.fit(X_train, y_train)
        except Exception as e:
            if verbose:
                print(f"  [ERROR] {fold_name}: training failed — {e}")
                traceback.print_exc()
            continue
        train_time = time.time() - t0

        # Predict
        if hasattr(model, "predict_proba"):
            y_score = model.predict_proba(X_test)[:, 1]
        else:
            y_score = model.decision_function(X_test)

        auroc = float(roc_auc_score(y_test, y_score))

        fold_predictions.append({
            "fold_name": fold_name,
            "test_mission": fold["test_mission"],
            "y_true": y_test.tolist(),
            "y_score": y_score.tolist(),
            "auroc": round(auroc, 4),
            "n_test": int(len(y_test)),
            "n_flight_test": int(y_test.sum()),
            "n_ground_test": int((y_test == 0).sum()),
            "train_time_s": round(train_time, 2),
        })

        if verbose:
            print(f"  {fold_name}: AUROC={auroc:.3f} (n={len(y_test)}, t={train_time:.1f}s)")

    return fold_predictions


def verify_against_m1(tissue, method_key, fold_predictions, tolerance=0.001):
    """Verify re-computed AUROCs match M1 JSONs."""
    m1_path = V4_EVAL_DIR / f"M1_{tissue}_gene_{method_key}.json"
    if not m1_path.exists():
        return None, "M1 JSON not found"

    with open(m1_path) as f:
        m1 = json.load(f)

    m1_folds = {fd["fold_name"]: fd["auroc"] for fd in m1["folds"]}
    mismatches = []

    for fp in fold_predictions:
        fn = fp["fold_name"]
        if fn in m1_folds:
            diff = abs(fp["auroc"] - m1_folds[fn])
            if diff > tolerance:
                mismatches.append({
                    "fold": fn,
                    "collected": fp["auroc"],
                    "m1": m1_folds[fn],
                    "diff": round(diff, 4)
                })

    if mismatches:
        return False, mismatches
    return True, "All folds match"


def main():
    parser = argparse.ArgumentParser(description="Collect raw predictions for DeLong test")
    parser.add_argument("--tissue", type=str, default=None,
                        help="Tissue to process (default: all)")
    parser.add_argument("--method", type=str, default=None,
                        help="Method key (default: all)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=str, default=None,
                        help="Output JSON path (default: META_predictions.json)")
    parser.add_argument("--no-verify", action="store_true",
                        help="Skip M1 verification")
    args = parser.parse_args()

    tissues = [args.tissue] if args.tissue else ALL_TISSUES
    methods = [args.method] if args.method else list(CLASSIFIERS.keys())

    output_path = Path(args.output) if args.output else V4_EVAL_DIR / "META_predictions.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Load existing results if incremental
    if output_path.exists() and (args.tissue or args.method):
        with open(output_path) as f:
            results = json.load(f)
    else:
        results = {}

    total_combos = len(tissues) * len(methods)
    completed = 0
    all_verified = True
    verification_report = []

    t_start = time.time()

    for tissue in tissues:
        if tissue not in results:
            results[tissue] = {}

        for method_key in methods:
            completed += 1
            label = CLASSIFIERS[method_key][0]
            print(f"\n[{completed}/{total_combos}] {tissue} / {label}")

            fold_preds = collect_single(tissue, method_key, seed=args.seed)
            results[tissue][method_key] = fold_preds

            # Verify against M1
            if not args.no_verify:
                ok, detail = verify_against_m1(tissue, method_key, fold_preds)
                status = "PASS" if ok is True else ("SKIP" if ok is None else "FAIL")
                verification_report.append({
                    "tissue": tissue, "method": method_key,
                    "status": status, "detail": detail
                })
                if ok is False:
                    all_verified = False
                    print(f"  *** VERIFICATION FAILED: {detail}")
                elif ok is True:
                    print(f"  ✓ Verified against M1")

    elapsed = time.time() - t_start

    # Save results
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n{'='*60}")
    print(f"Saved: {output_path}")
    print(f"Total time: {elapsed:.1f}s ({elapsed/60:.1f}m)")
    print(f"Combos: {completed}, Tissues: {len(tissues)}, Methods: {len(methods)}")

    # Verification summary
    if not args.no_verify:
        n_pass = sum(1 for v in verification_report if v["status"] == "PASS")
        n_fail = sum(1 for v in verification_report if v["status"] == "FAIL")
        n_skip = sum(1 for v in verification_report if v["status"] == "SKIP")
        print(f"\nVerification: {n_pass} PASS, {n_fail} FAIL, {n_skip} SKIP")

        if n_fail > 0:
            print("\n*** FAILURES (AUROC mismatch > 0.001) ***")
            for v in verification_report:
                if v["status"] == "FAIL":
                    print(f"  {v['tissue']}/{v['method']}: {v['detail']}")

        # Save verification report
        verify_path = output_path.with_name("META_predictions_verification.json")
        with open(verify_path, 'w') as f:
            json.dump({
                "summary": {"pass": n_pass, "fail": n_fail, "skip": n_skip},
                "all_verified": all_verified,
                "details": verification_report,
                "tolerance": 0.001,
                "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
            }, f, indent=2)
        print(f"Verification report: {verify_path}")

    return 0 if all_verified else 1


if __name__ == "__main__":
    sys.exit(main())

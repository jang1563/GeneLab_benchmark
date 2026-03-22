#!/usr/bin/env python3
"""
ablation_bootstrap.py — GeneLabBench v4: Bootstrap Stability Analysis

When does bootstrap CI converge? Validates n_bootstrap=2000 is sufficient.

Design:
- n_bootstrap values: [100, 500, 1000, 2000, 5000]
- Gene features, PCA-LR, all 8 tissues
- One training pass per tissue → store (y_test, y_score) → re-bootstrap at all levels
- Outputs ALL 5 n_bootstrap levels in a single JSON per tissue
- No SLURM needed (~10s per tissue)
- AUROC at each level must match M1 gene PCA-LR (regression test)

Usage:
  python ablation_bootstrap.py --tissue liver
  python ablation_bootstrap.py --tissue all
"""

import json
import time
import argparse
import warnings
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
    bootstrap_auroc,
    get_folds,
)
from classifier_registry import get_classifier, adapt_pca_components

BOOTSTRAP_LEVELS = [100, 500, 1000, 2000, 5000]


def evaluate_bootstrap_ablation(tissue, seed=42, verbose=True):
    """Train PCA-LR once, then compute bootstrap CI at multiple levels."""
    t_start = time.time()

    if verbose:
        print(f"\n{'='*60}")
        print(f"Bootstrap Ablation | Tissue: {tissue} | PCA-LR")
        print(f"{'='*60}")

    folds = get_folds(tissue)

    # Train once, store predictions
    fold_predictions = []

    for fold in folds:
        X_train = fold["train_X"]
        y_train = fold["train_y"]
        X_test = fold["test_X"]
        y_test = fold["test_y"]
        test_mission = fold["test_mission"]

        if len(np.unique(y_test)) < 2:
            continue

        _, model = get_classifier("pca_lr", seed)
        adapt_pca_components(model, X_train.shape[0], X_train.shape[1])
        model.fit(X_train, y_train)
        y_score = model.predict_proba(X_test)[:, 1]
        auroc = float(roc_auc_score(y_test, y_score))

        fold_predictions.append({
            "test_mission": test_mission,
            "fold_name": fold["fold_name"],
            "y_test": y_test,
            "y_score": y_score,
            "auroc": auroc,
            "n_test": fold["n_test"],
        })

        if verbose:
            print(f"  {test_mission}: AUROC={auroc:.3f}")

    if not fold_predictions:
        return None

    # Compute bootstrap CI at all levels
    levels_results = {}

    for n_boot in BOOTSTRAP_LEVELS:
        fold_results = []
        for fp in fold_predictions:
            mean_bs, lower_bs, upper_bs = bootstrap_auroc(
                fp["y_test"], fp["y_score"], n_bootstrap=n_boot)
            ci_width = upper_bs - lower_bs if not (np.isnan(upper_bs) or np.isnan(lower_bs)) else float("nan")

            fold_results.append({
                "test_mission": fp["test_mission"],
                "auroc": round(fp["auroc"], 4),
                "bootstrap_mean": round(mean_bs, 4),
                "ci_lower": round(lower_bs, 4),
                "ci_upper": round(upper_bs, 4),
                "ci_width": round(ci_width, 4) if not np.isnan(ci_width) else None,
                "n_test": fp["n_test"],
            })

        aurocs = [r["auroc"] for r in fold_results]
        ci_widths = [r["ci_width"] for r in fold_results if r["ci_width"] is not None]

        levels_results[str(n_boot)] = {
            "n_bootstrap": n_boot,
            "mean_auroc": round(float(np.mean(aurocs)), 4),
            "std_auroc": round(float(np.std(aurocs)), 4),
            "mean_ci_width": round(float(np.mean(ci_widths)), 4) if ci_widths else None,
            "std_ci_width": round(float(np.std(ci_widths)), 4) if ci_widths else None,
            "mean_ci_lower": round(float(np.mean([r["ci_lower"] for r in fold_results])), 4),
            "mean_ci_upper": round(float(np.mean([r["ci_upper"] for r in fold_results])), 4),
            "folds": fold_results,
        }

        if verbose:
            mean_width = float(np.mean(ci_widths)) if ci_widths else float("nan")
            print(f"  n_boot={n_boot:5d}: mean_CI_width={mean_width:.4f}")

    total_time = time.time() - t_start

    summary = {
        "experiment": "ablation_bootstrap",
        "model": "PCA-LR",
        "model_key": "pca_lr",
        "tissue": tissue,
        "feature_type": "gene",
        "cv_strategy": "LOMO" if tissue in LOMO_TISSUES else "5-fold_stratified",
        "n_folds": len(fold_predictions),
        "bootstrap_levels": BOOTSTRAP_LEVELS,
        "levels": levels_results,
        "seed": seed,
        "wall_time_sec": round(total_time, 1),
        "timestamp": datetime.now().isoformat(),
    }

    if verbose:
        print(f"\n  Wall time: {total_time:.1f}s")

    return summary


def main():
    parser = argparse.ArgumentParser(description="Bootstrap Stability Analysis")
    parser.add_argument("--tissue", required=True,
                        help="Tissue name (or 'all' for all 8 tissues)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=str, default=None)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    verbose = not args.quiet
    output_dir = Path(args.output_dir) if args.output_dir else V4_EVAL_DIR

    if args.tissue == "all":
        tissues = list(TISSUE_MISSIONS.keys())
    else:
        if args.tissue not in TISSUE_MISSIONS:
            print(f"[ERROR] Unknown tissue: {args.tissue}")
            return
        tissues = [args.tissue]

    for tissue in tissues:
        try:
            result = evaluate_bootstrap_ablation(
                tissue=tissue, seed=args.seed, verbose=verbose)
            if result:
                fname = f"ABL_boot_{tissue}.json"
                output_dir.mkdir(parents=True, exist_ok=True)
                (output_dir / fname).write_text(json.dumps(result, indent=2))
                if verbose:
                    print(f"  Saved: {output_dir / fname}")
        except Exception as e:
            print(f"[ERROR] {tissue}: {e}")
            import traceback
            traceback.print_exc()


if __name__ == "__main__":
    main()

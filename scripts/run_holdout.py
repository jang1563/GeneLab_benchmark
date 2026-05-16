#!/usr/bin/env python3
"""
run_holdout.py — GeneLab_benchmark: Held-Out Evaluation

Runs baseline classifiers on held-out test sets and reports AUROC with
Bootstrap 95% CI and permutation p-values.

Reuses model definitions and evaluation logic from run_baselines.py.

Usage:
    python scripts/run_holdout.py --tissue skin --holdout-mission RR-7
    python scripts/run_holdout.py --tissue thymus --holdout-mission RR-23
    python scripts/run_holdout.py --all
"""

import json
import sys
import argparse
import warnings
from pathlib import Path
from datetime import datetime

warnings.filterwarnings("ignore")

BASE_DIR = Path(__file__).resolve().parent.parent
TASKS_DIR = BASE_DIR / "tasks"
EVAL_DIR = BASE_DIR / "evaluation"

# Add scripts to path for imports
sys.path.insert(0, str(BASE_DIR / "scripts"))
from run_baselines import (
    evaluate_fold,
    build_lr, build_rf, build_pca_lr,
)

# ── Holdout configurations ────────────────────────────────────────────────────
HOLDOUTS = {
    "thymus_RR-23": {
        "tissue": "thymus",
        "mission": "RR-23",
        "task": "A4",
        "fold_dir": TASKS_DIR / "A4_thymus_lomo" / "fold_RR-23_holdout",
        "output": "A4_holdout_results.json",
    },
    "skin_RR-7": {
        "tissue": "skin",
        "mission": "RR-7",
        "task": "A5",
        "fold_dir": TASKS_DIR / "A5_skin_lomo" / "fold_RR-7_holdout",
        "output": "A5_holdout_results.json",
    },
}

MODELS = {
    "lr": ("Logistic Regression (ElasticNet)", build_lr),
    "rf": ("Random Forest", build_rf),
    "pca_lr": ("PCA-50 + LogReg", build_pca_lr),
}


class _NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        import numpy as np
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj) if np.isfinite(obj) else None
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def run_holdout(config, n_bootstrap=2000, n_perm=10000):
    """Run all models on a held-out fold."""
    fold_dir = config["fold_dir"]
    tissue = config["tissue"]
    mission = config["mission"]

    if not fold_dir.exists():
        print(f"ERROR: {fold_dir} not found")
        return None

    print(f"\n{'='*50}")
    print(f"Held-Out: {tissue} / {mission}")
    print(f"Fold: {fold_dir}")
    print(f"{'='*50}")

    results = {}
    for model_name, (desc, builder) in MODELS.items():
        print(f"\n  Model: {desc}")
        model = builder()
        result = evaluate_fold(
            fold_dir, model_name, model,
            n_bootstrap=n_bootstrap, n_perm=n_perm,
            verbose=True,
        )
        if result is not None:
            results[model_name] = result

    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tissue", default=None,
                        help="Tissue (skin or thymus)")
    parser.add_argument("--holdout-mission", default=None,
                        help="Held-out mission (e.g., RR-7, RR-23)")
    parser.add_argument("--all", action="store_true",
                        help="Run all configured holdouts")
    parser.add_argument("--n-bootstrap", type=int, default=2000,
                        help="Number of bootstrap iterations (default: 2000)")
    parser.add_argument("--n-perm", type=int, default=10000,
                        help="Number of permutation iterations (default: 10000)")
    args = parser.parse_args()

    configs_to_run = []

    if args.all:
        configs_to_run = list(HOLDOUTS.values())
    elif args.tissue and args.holdout_mission:
        key = f"{args.tissue}_{args.holdout_mission}"
        if key not in HOLDOUTS:
            print(f"ERROR: Unknown holdout '{key}'. Available: {list(HOLDOUTS.keys())}")
            sys.exit(1)
        configs_to_run = [HOLDOUTS[key]]
    else:
        print("Specify --tissue and --holdout-mission, or --all")
        sys.exit(1)

    all_results = {}
    for config in configs_to_run:
        results = run_holdout(config, args.n_bootstrap, args.n_perm)
        if results:
            # Write individual output
            out_path = EVAL_DIR / config["output"]
            with open(out_path, "w") as f:
                json.dump(results, f, indent=2, cls=_NumpyEncoder)
            print(f"\nWrote {out_path}")
            all_results[f"{config['tissue']}_{config['mission']}"] = results

    # Print comparison table if multiple holdouts
    if len(all_results) > 1:
        print(f"\n{'='*60}")
        print("Cross-Tissue Held-Out Comparison")
        print(f"{'='*60}")
        print(f"{'Holdout':<20} {'Model':<10} {'AUROC':<8} {'CI':<20} {'n_test':<8}")
        for key, results in all_results.items():
            for model_name, r in results.items():
                print(f"{key:<20} {model_name:<10} {r['auroc']:.3f}  "
                      f"[{r['ci_lower']:.3f}, {r['ci_upper']:.3f}]  {r['n_test']}")


if __name__ == "__main__":
    main()

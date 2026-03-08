#!/usr/bin/env python3
"""
aggregate_geneformer_results.py — Collect Geneformer LOMO results across tissues.

Collects Geneformer LOMO results across tissues and compares with Tier 1 baselines.

Supports two input modes:
  1. Pre-aggregated: evaluation/geneformer_{model}_{task}_lomo_results.json
  2. Per-fold (HPC array jobs): tasks/A*_lomo/fold_*/geneformer_tokens/{model}/finetune_result.json
     → auto-aggregates and saves evaluation JSON for future use

Input:
  evaluation/geneformer_{model}_{task}_lomo_results.json  (per-tissue, from finetune.py --fold lomo)
  OR tasks/A*_lomo/fold_*/geneformer_tokens/{model}/finetune_result.json  (per-fold, from HPC)
  evaluation/A{1..6}_baseline_results.json                (Tier 1 baselines)

Output:
  evaluation/geneformer_{model}_all_tissues_summary.json
  stdout: comparison table (markdown)

Usage:
  python scripts/aggregate_geneformer_results.py
  python scripts/aggregate_geneformer_results.py --model-version mouse_gf
  python scripts/aggregate_geneformer_results.py --model-version v1
"""

import json
import argparse
from pathlib import Path
from datetime import datetime

BASE_DIR = Path(__file__).resolve().parent.parent
EVAL_DIR = BASE_DIR / "evaluation"

# Task → tissue mapping
TASK_TISSUE = {
    "A1": "liver",
    "A2": "gastrocnemius",
    "A3": "kidney",
    "A4": "thymus",
    "A5": "skin",
    "A6": "eye",
}

# Best baseline model per task (from Phase 1 results)
BEST_BASELINE = {
    "A1": "lr",           # LR ElasticNet
    "A2": "lr",           # LR ElasticNet (converged)
    "A3": "lr",           # LR
    "A4": "pca_lr",       # PCA-LR
    "A5": "lr",           # LR ElasticNet
    "A6": "pca_lr",       # PCA-LR (pathway)
}


TASKS_DIR = BASE_DIR / "tasks"

# Task → task_dir mapping
TASK_DIRS = {
    "A1": "A1_liver_lomo",
    "A2": "A2_gastrocnemius_lomo",
    "A3": "A3_kidney_lomo",
    "A4": "A4_thymus_lomo",
    "A5": "A5_skin_lomo",
    "A6": "A6_eye_lomo",
}


def load_geneformer_results(model_version: str) -> dict:
    """Load Geneformer LOMO results.

    First checks for pre-aggregated evaluation JSONs.
    Falls back to scanning per-fold finetune_result.json files
    in the task directories (produced by HPC array jobs).
    """
    import numpy as np

    results = {}
    for task in TASK_TISSUE:
        # 1. Try pre-aggregated LOMO JSON
        path = EVAL_DIR / f"geneformer_{model_version}_{task}_lomo_results.json"
        if path.exists():
            with open(path) as f:
                results[task] = json.load(f)
            continue

        # 2. Fallback: collect per-fold results from task directory
        task_dir_name = TASK_DIRS.get(task)
        if not task_dir_name:
            continue
        task_dir = TASKS_DIR / task_dir_name
        if not task_dir.exists():
            continue

        fold_results = []
        for fold_dir in sorted(task_dir.glob("fold_*_test")):
            result_file = fold_dir / "geneformer_tokens" / model_version / "finetune_result.json"
            if result_file.exists():
                with open(result_file) as f:
                    fold_results.append(json.load(f))

        if not fold_results:
            continue

        aurocs = [r["best_test_auroc"] for r in fold_results if r.get("status") == "ok"]
        if not aurocs:
            continue

        lomo_result = {
            "task": task,
            "model": f"Geneformer-{model_version}",
            "tissue": TASK_TISSUE[task],
            "mean_auroc": float(np.mean(aurocs)),
            "std_auroc": float(np.std(aurocs)),
            "n_folds": len(aurocs),
            "fold_results": fold_results,
        }

        # Save aggregated JSON for future use
        EVAL_DIR.mkdir(exist_ok=True)
        out_path = EVAL_DIR / f"geneformer_{model_version}_{task}_lomo_results.json"
        with open(out_path, "w") as f:
            json.dump(lomo_result, f, indent=2)

        results[task] = lomo_result
    return results


def load_baseline_results() -> dict:
    """Load Tier 1 baseline results for comparison."""
    baselines = {}
    for task, tissue in TASK_TISSUE.items():
        path = EVAL_DIR / f"{task}_baseline_results.json"
        if not path.exists():
            continue
        with open(path) as f:
            data = json.load(f)
        # Get the best model for this task
        best_key = BEST_BASELINE.get(task, "pca_lr")
        if best_key in data:
            baselines[task] = {
                "model": data[best_key].get("model", best_key),
                "mean_auroc": data[best_key].get("mean_auroc", None),
                "n_folds": data[best_key].get("n_folds", None),
            }
        else:
            # Fallback: pick first available model
            for key, val in data.items():
                if isinstance(val, dict) and "mean_auroc" in val:
                    baselines[task] = {
                        "model": val.get("model", key),
                        "mean_auroc": val["mean_auroc"],
                        "n_folds": val.get("n_folds", None),
                    }
                    break
    return baselines


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model-version", default="mouse_gf",
                        choices=["v1", "v2", "mouse_gf"],
                        help="Geneformer model version (default: mouse_gf)")
    args = parser.parse_args()

    model = args.model_version
    gf_results = load_geneformer_results(model)
    baselines = load_baseline_results()

    if not gf_results:
        print(f"No Geneformer results found for model '{model}'.")
        print(f"Expected files: evaluation/geneformer_{model}_A*_lomo_results.json")
        print(f"\nRun fine-tuning first:")
        print(f"  bash scripts/hpc_submit_all_tissues.sh")
        return

    # Build summary
    summary = {
        "model": f"Geneformer-{model}",
        "timestamp": datetime.now().isoformat(),
        "n_tissues": len(gf_results),
        "tissues": {},
    }

    # Print comparison table
    print()
    print(f"## Geneformer ({model}) vs Tier 1 Baseline — LOMO AUROC")
    print()
    print("| Task | Tissue | Geneformer AUROC | Baseline AUROC | Baseline Model | Delta | Winner |")
    print("|------|--------|-----------------|---------------|----------------|-------|--------|")

    gf_aurocs = []
    bl_aurocs = []

    for task in sorted(TASK_TISSUE.keys()):
        tissue = TASK_TISSUE[task]

        if task in gf_results:
            gf = gf_results[task]
            gf_auroc = gf.get("mean_auroc", None)
            gf_std = gf.get("std_auroc", 0)
            gf_n = gf.get("n_folds", 0)
        else:
            gf_auroc = None
            gf_std = 0
            gf_n = 0

        if task in baselines:
            bl = baselines[task]
            bl_auroc = bl.get("mean_auroc", None)
            bl_model = bl.get("model", "?")
        else:
            bl_auroc = None
            bl_model = "N/A"

        # Format
        gf_str = f"{gf_auroc:.3f}" if gf_auroc is not None else "—"
        bl_str = f"{bl_auroc:.3f}" if bl_auroc is not None else "—"

        if gf_auroc is not None and bl_auroc is not None:
            delta = gf_auroc - bl_auroc
            delta_str = f"{delta:+.3f}"
            winner = "Geneformer" if delta > 0.01 else ("Baseline" if delta < -0.01 else "Tie")
            gf_aurocs.append(gf_auroc)
            bl_aurocs.append(bl_auroc)
        else:
            delta_str = "—"
            winner = "—"

        print(f"| {task} | {tissue.capitalize()} | {gf_str} | {bl_str} | {bl_model} | {delta_str} | {winner} |")

        # Add to summary
        summary["tissues"][task] = {
            "tissue": tissue,
            "geneformer_mean_auroc": gf_auroc,
            "geneformer_std_auroc": gf_std,
            "geneformer_n_folds": gf_n,
            "baseline_mean_auroc": bl_auroc,
            "baseline_model": bl_model,
            "delta": (gf_auroc - bl_auroc) if (gf_auroc and bl_auroc) else None,
        }

    # Overall
    if gf_aurocs:
        import numpy as np
        gf_mean = float(np.mean(gf_aurocs))
        bl_mean = float(np.mean(bl_aurocs)) if bl_aurocs else None
        print(f"|---|---|---|---|---|---|---|")
        print(f"| **Mean** | **{len(gf_aurocs)} tissues** | **{gf_mean:.3f}** | "
              f"**{bl_mean:.3f}** | — | **{gf_mean - bl_mean:+.3f}** | "
              f"**{'Geneformer' if gf_mean > bl_mean else 'Baseline'}** |")

        summary["overall"] = {
            "geneformer_mean": gf_mean,
            "baseline_mean": bl_mean,
            "n_compared": len(gf_aurocs),
        }

    # Save
    out_path = EVAL_DIR / f"geneformer_{model}_all_tissues_summary.json"
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved to {out_path}")

    # Per-fold details
    print(f"\n### Per-Fold Details")
    for task in sorted(gf_results.keys()):
        gf = gf_results[task]
        tissue = TASK_TISSUE[task]
        print(f"\n**{task} ({tissue})**: mean={gf['mean_auroc']:.3f} ± {gf.get('std_auroc', 0):.3f}")
        for fold in gf.get("fold_results", []):
            auroc = fold.get("best_test_auroc", fold.get("auroc", None))
            ci_low = fold.get("ci_low", None)
            ci_high = fold.get("ci_high", None)
            fname = fold.get("fold", "?")
            auroc_str = f"{auroc:.3f}" if isinstance(auroc, (int, float)) else "?"
            if isinstance(ci_low, (int, float)) and isinstance(ci_high, (int, float)):
                ci_str = f"[{ci_low:.3f}, {ci_high:.3f}]"
            else:
                ci_str = "[?, ?]"
            print(f"  {fname}: AUROC={auroc_str} CI={ci_str}")


if __name__ == "__main__":
    main()

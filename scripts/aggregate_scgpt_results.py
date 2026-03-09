#!/usr/bin/env python3
"""
aggregate_scgpt_results.py
Collect per-fold scGPT result JSONs → combined summary JSON.

Usage:
    python scripts/aggregate_scgpt_results.py \
        --eval-base evaluation/ \
        --model-version whole_human \
        --output evaluation/scgpt_whole_human_all_tissues_summary.json
"""
import argparse
import json
import glob
import numpy as np
from pathlib import Path

TISSUE_MAP = {
    "A1": "liver",
    "A2": "gastrocnemius",
    "A3": "kidney",
    "A4": "thymus",
    "A5": "skin",
    "A6": "eye",
}

# Baseline (PCA-LR LOMO) mean AUROCs from RESULTS_SUMMARY.md
BASELINE_AUROC = {
    "A1": 0.588,   # Liver LR ElasticNet best
    "A2": 0.801,   # Gastrocnemius
    "A3": 0.538,   # Kidney
    "A4": 0.923,   # Thymus
    "A5": 0.821,   # Skin
    "A6": 0.789,   # Eye
}

# Geneformer mean AUROCs from existing results
GENEFORMER_AUROC = {
    "A1": 0.486,
    "A2": 0.432,
    "A3": 0.432,
    "A4": 0.476,
    "A5": 0.532,
    "A6": 0.478,
}


def load_fold_results(eval_base: Path, model_version: str, task: str) -> list[dict]:
    pattern = str(eval_base / f"scgpt_{model_version}_{task}_*_result.json")
    files = sorted(glob.glob(pattern))
    results = []
    for f in files:
        with open(f) as fh:
            results.append(json.load(fh))
    return results


def summarize_task(task: str, folds: list[dict]) -> dict:
    aurocs = [f["auroc"] for f in folds]
    mean_auroc = float(np.mean(aurocs))
    std_auroc = float(np.std(aurocs))
    n_folds = len(folds)
    baseline = BASELINE_AUROC.get(task, None)
    gf = GENEFORMER_AUROC.get(task, None)
    return {
        "task": task,
        "tissue": TISSUE_MAP.get(task, task),
        "model": f"scgpt_{args.model_version}",
        "mean_auroc": round(mean_auroc, 4),
        "std_auroc": round(std_auroc, 4),
        "n_folds": n_folds,
        "baseline_auroc": baseline,
        "geneformer_auroc": gf,
        "delta_vs_baseline": round(mean_auroc - baseline, 4) if baseline else None,
        "delta_vs_geneformer": round(mean_auroc - gf, 4) if gf else None,
        "fold_results": [
            {
                "fold": f["test_mission"],
                "auroc": f["auroc"],
                "best_epoch": f["best_epoch"],
                "n_train": f["n_train"],
                "n_test": f["n_test"],
                "epoch_aurocs": f.get("epoch_aurocs", []),
            }
            for f in sorted(folds, key=lambda x: x["test_mission"])
        ],
    }


def main():
    global args
    parser = argparse.ArgumentParser()
    parser.add_argument("--eval-base", default="evaluation")
    parser.add_argument("--model-version", default="whole_human")
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    eval_base = Path(args.eval_base)
    output = Path(args.output) if args.output else \
        eval_base / f"scgpt_{args.model_version}_all_tissues_summary.json"

    all_tasks = {}
    all_aurocs = []
    print(f"\nscGPT {args.model_version} — Results Aggregation")
    print("=" * 60)

    for task in sorted(TISSUE_MAP.keys()):
        folds = load_fold_results(eval_base, args.model_version, task)
        if not folds:
            print(f"  {task} ({TISSUE_MAP[task]}): no results found — skipping")
            continue
        summary = summarize_task(task, folds)
        all_tasks[task] = summary
        all_aurocs.append(summary["mean_auroc"])

        gf = summary["geneformer_auroc"]
        bl = summary["baseline_auroc"]
        print(f"  {task} {TISSUE_MAP[task]:14s}: scGPT={summary['mean_auroc']:.3f}"
              f"  GF={gf:.3f}  Baseline={bl:.3f}"
              f"  Δ_GF={summary['delta_vs_geneformer']:+.3f}"
              f"  Δ_BL={summary['delta_vs_baseline']:+.3f}"
              f"  ({summary['n_folds']} folds)")

    overall_mean = float(np.mean(all_aurocs)) if all_aurocs else None
    gf_overall = 0.476
    bl_overall = 0.758

    print("=" * 60)
    print(f"  Overall mean AUROC: {overall_mean:.3f}  (GF={gf_overall:.3f}, Baseline={bl_overall:.3f})")
    print(f"  Δ vs Geneformer:   {overall_mean - gf_overall:+.3f}")
    print(f"  Δ vs Baseline:     {overall_mean - bl_overall:+.3f}")

    output_dict = {
        "model": f"scgpt_{args.model_version}",
        "description": "scGPT whole_human fine-tuned on mouse spaceflight LOMO folds (10 epochs, batch=8, lr=1e-4, freeze=10/12 layers)",
        "overall_mean_auroc": round(overall_mean, 4),
        "geneformer_mean_auroc": gf_overall,
        "baseline_mean_auroc": bl_overall,
        "delta_vs_geneformer": round(overall_mean - gf_overall, 4),
        "delta_vs_baseline": round(overall_mean - bl_overall, 4),
        "n_tasks": len(all_tasks),
        "n_folds_total": sum(t["n_folds"] for t in all_tasks.values()),
        "tasks": all_tasks,
    }

    with open(output, "w") as f:
        json.dump(output_dict, f, indent=2)
    print(f"\nSaved: {output}")


if __name__ == "__main__":
    main()

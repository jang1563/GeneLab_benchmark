#!/usr/bin/env python3
"""
aggregate_scgpt_results.py
Collect per-fold scGPT result JSONs → combined summary JSON.

Usage:
    python scripts/aggregate_scgpt_results.py \
        --eval-base evaluation/ \
        --model-version whole_human

This discovers fold results in evaluation/scgpt/ when invoked from the repo root
and writes the combined summary alongside those per-fold result JSONs by default.
"""
import argparse
import json
import glob
import numpy as np
from pathlib import Path
from typing import Optional

REPO_ROOT = Path(__file__).resolve().parents[1]
TISSUE_MAP = {
    "A1": "liver",
    "A2": "gastrocnemius",
    "A3": "kidney",
    "A4": "thymus",
    "A5": "skin",
    "A6": "eye",
}

REFERENCE_SUMMARY_NAME = "geneformer_mouse_gf_all_tissues_summary.json"


def resolve_results_dir(eval_base: Path, model_version: str) -> Path:
    candidates = [eval_base, eval_base / "scgpt"]
    for candidate in candidates:
        pattern = str(candidate / f"scgpt_{model_version}_*_result.json")
        if glob.glob(pattern):
            return candidate
    return eval_base / "scgpt" if (eval_base / "scgpt").exists() else eval_base


def resolve_reference_summary(results_dir: Path, eval_base: Path, reference_summary: Optional[str]) -> Path:
    if reference_summary:
        return Path(reference_summary)

    candidates = [
        eval_base / REFERENCE_SUMMARY_NAME,
        results_dir / REFERENCE_SUMMARY_NAME,
        results_dir.parent / REFERENCE_SUMMARY_NAME,
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate

    searched = ", ".join(str(candidate) for candidate in candidates)
    raise FileNotFoundError(
        "Could not find Geneformer/baseline reference summary. "
        f"Searched: {searched}"
    )


def load_reference_metrics(summary_path: Path) -> dict:
    with open(summary_path) as fh:
        payload = json.load(fh)

    task_metrics = {}
    for task, tissue_summary in payload.get("tissues", {}).items():
        task_metrics[task] = {
            "baseline": float(tissue_summary["baseline_mean_auroc"]),
            "geneformer": float(tissue_summary["geneformer_mean_auroc"]),
        }

    overall = payload.get("overall", {})
    return {
        "tasks": task_metrics,
        "overall": {
            "baseline": float(overall["baseline_mean"]),
            "geneformer": float(overall["geneformer_mean"]),
        },
        "source": str(summary_path),
    }


def load_fold_results(results_dir: Path, model_version: str, task: str) -> list[dict]:
    pattern = str(results_dir / f"scgpt_{model_version}_{task}_*_result.json")
    files = sorted(glob.glob(pattern))
    results = []
    for f in files:
        with open(f) as fh:
            results.append(json.load(fh))
    return results


def summarize_task(task: str, folds: list[dict], model_version: str, reference_metrics: dict) -> dict:
    aurocs = [f["auroc"] for f in folds]
    mean_auroc = float(np.mean(aurocs))
    std_auroc = float(np.std(aurocs))
    n_folds = len(folds)
    task_reference = reference_metrics.get(task, {})
    baseline = task_reference.get("baseline", None)
    gf = task_reference.get("geneformer", None)
    return {
        "task": task,
        "tissue": TISSUE_MAP.get(task, task),
        "model": f"scgpt_{model_version}",
        "mean_auroc": round(mean_auroc, 4),
        "std_auroc": round(std_auroc, 4),
        "n_folds": n_folds,
        "baseline_auroc": round(baseline, 4) if baseline is not None else None,
        "geneformer_auroc": round(gf, 4) if gf is not None else None,
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


def format_metric(value: Optional[float]) -> str:
    return f"{value:.3f}" if value is not None else "n/a"


def repo_relative_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--eval-base", default="evaluation")
    parser.add_argument("--model-version", default="whole_human")
    parser.add_argument("--output", default=None)
    parser.add_argument("--reference-summary", default=None)
    args = parser.parse_args()

    eval_base = Path(args.eval_base)
    results_dir = resolve_results_dir(eval_base, args.model_version)
    output = (
        Path(args.output)
        if args.output
        else results_dir / f"scgpt_{args.model_version}_all_tissues_summary.json"
    )
    reference_summary = resolve_reference_summary(results_dir, eval_base, args.reference_summary)
    reference_metrics = load_reference_metrics(reference_summary)

    all_tasks = {}
    all_aurocs = []
    print(f"\nscGPT {args.model_version} — Results Aggregation")
    print("=" * 60)
    print(f"  Result files: {repo_relative_path(results_dir)}")
    print(f"  Reference summary: {repo_relative_path(reference_summary)}")

    for task in sorted(TISSUE_MAP.keys()):
        folds = load_fold_results(results_dir, args.model_version, task)
        if not folds:
            print(f"  {task} ({TISSUE_MAP[task]}): no results found — skipping")
            continue
        summary = summarize_task(task, folds, args.model_version, reference_metrics["tasks"])
        all_tasks[task] = summary
        all_aurocs.append(summary["mean_auroc"])

        gf = summary["geneformer_auroc"]
        bl = summary["baseline_auroc"]
        print(
            f"  {task} {TISSUE_MAP[task]:14s}: scGPT={summary['mean_auroc']:.3f}"
            f"  GF={format_metric(gf)}  Baseline={format_metric(bl)}"
            f"  Δ_GF={summary['delta_vs_geneformer']:+.3f}"
            f"  Δ_BL={summary['delta_vs_baseline']:+.3f}"
            f"  ({summary['n_folds']} folds)"
        )

    overall_mean = float(np.mean(all_aurocs)) if all_aurocs else None
    gf_overall = reference_metrics["overall"]["geneformer"]
    bl_overall = reference_metrics["overall"]["baseline"]

    print("=" * 60)
    print(f"  Overall mean AUROC: {overall_mean:.3f}  (GF={gf_overall:.3f}, Baseline={bl_overall:.3f})")
    print(f"  Δ vs Geneformer:   {overall_mean - gf_overall:+.3f}")
    print(f"  Δ vs Baseline:     {overall_mean - bl_overall:+.3f}")

    output_dict = {
        "model": f"scgpt_{args.model_version}",
        "description": "scGPT whole_human fine-tuned on mouse spaceflight LOMO folds (10 epochs, batch=8, lr=1e-4, freeze=10/12 layers)",
        "overall_mean_auroc": round(overall_mean, 4),
        "geneformer_mean_auroc": round(gf_overall, 4),
        "baseline_mean_auroc": round(bl_overall, 4),
        "delta_vs_geneformer": round(overall_mean - gf_overall, 4),
        "delta_vs_baseline": round(overall_mean - bl_overall, 4),
        "n_tasks": len(all_tasks),
        "n_folds_total": sum(t["n_folds"] for t in all_tasks.values()),
        "reference_summary": repo_relative_path(reference_summary),
        "tasks": all_tasks,
    }

    with open(output, "w") as f:
        json.dump(output_dict, f, indent=2)
    print(f"\nSaved: {output}")


if __name__ == "__main__":
    main()

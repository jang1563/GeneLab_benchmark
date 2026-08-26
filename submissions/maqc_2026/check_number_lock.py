#!/usr/bin/env python3
"""Recompute every numerical result reported in the MAQC 2026 abstract."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def load_json(relative_path: str) -> dict:
    with (ROOT / relative_path).open(encoding="utf-8") as handle:
        return json.load(handle)


scgpt = load_json("evaluation/scgpt/scgpt_whole_human_all_tissues_summary.json")
geneformer_summary = load_json("evaluation/geneformer_mouse_gf_all_tissues_summary.json")

geneformer_tissues = [
    load_json(f"evaluation/geneformer_mouse_gf_A{task}_lomo_results.json")
    for task in range(1, 7)
]
baseline_tissues = [
    load_json(f"evaluation/A{task}_baseline_results.json")
    for task in range(1, 7)
]
thymus = load_json("evaluation/submission_PCA-LR_baseline_v1_A4_eval.json")

scgpt_tasks = list(scgpt["tasks"].values())
scgpt_fixed_tissue_means = [
    sum(fold["epoch_aurocs"][-1] for fold in task["fold_results"])
    / len(task["fold_results"])
    for task in scgpt_tasks
]

geneformer_fixed_tissue_means = []
geneformer_folds = []
for task in geneformer_tissues:
    valid_folds = [fold for fold in task["fold_results"] if fold.get("status") == "ok"]
    geneformer_folds.extend(valid_folds)
    geneformer_fixed_tissue_means.append(
        sum(fold["history"][-1]["test_auroc"] for fold in valid_folds)
        / len(valid_folds)
    )

values = {
    "profiles_task_inventory": sum(fold["n_test"] for fold in geneformer_folds),
    "task_folds": len(geneformer_folds),
    "scgpt_result_folds": sum(len(task["fold_results"]) for task in scgpt_tasks),
    "scgpt_result_profiles": sum(
        fold["n_test"] for task in scgpt_tasks for fold in task["fold_results"]
    ),
    "geneformer_result_folds": len(geneformer_folds),
    "scgpt_initial_macro_auroc": scgpt["overall_mean_auroc"],
    "scgpt_epoch10_macro_auroc": sum(scgpt_fixed_tissue_means)
    / len(scgpt_fixed_tissue_means),
    "geneformer_initial_macro_auroc": geneformer_summary["overall"][
        "geneformer_mean"
    ],
    "geneformer_epoch10_macro_auroc": sum(geneformer_fixed_tissue_means)
    / len(geneformer_fixed_tissue_means),
    "pca_lr_macro_auroc": sum(task["pca_lr"]["mean_auroc"] for task in baseline_tissues)
    / len(baseline_tissues),
    "thymus_pca_lr_mean_mission_auroc": thymus["summary"]["mean_auroc"],
    "thymus_pca_lr_pooled_oof_auroc": thymus["summary"]["overall_auroc_pooled"],
}

print(json.dumps(values, indent=2, sort_keys=True))


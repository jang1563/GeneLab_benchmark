#!/usr/bin/env python3
"""Recompute every quantitative value used in the MAQC 2026 abstract.

The script reads only the frozen local evaluation artifacts and prints a
machine-readable JSON record.  The tissue, rather than the fold or sample, is
the averaging unit for the three six-tissue macro AUROCs.
"""

from __future__ import annotations

import hashlib
import json
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TASKS = tuple(f"A{i}" for i in range(1, 7))


def read_json(relative_path: str) -> dict:
    return json.loads((ROOT / relative_path).read_text())


def sha256(relative_path: str) -> str:
    return hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()


scgpt_path = "evaluation/scgpt/scgpt_whole_human_all_tissues_summary.json"
scgpt = read_json(scgpt_path)

scgpt_initial_tissue_aurocs = [
    scgpt["tasks"][task]["mean_auroc"] for task in TASKS
]
scgpt_epoch10_tissue_aurocs = [
    statistics.mean(
        fold["epoch_aurocs"][9]
        for fold in scgpt["tasks"][task]["fold_results"]
    )
    for task in TASKS
]
scgpt_folds = [
    fold
    for task in TASKS
    for fold in scgpt["tasks"][task]["fold_results"]
]

geneformer_paths = [
    f"evaluation/geneformer_mouse_gf_{task}_lomo_results.json"
    for task in TASKS
]
geneformer_results = [read_json(path) for path in geneformer_paths]
geneformer_folds_by_task = [
    [fold for fold in result["fold_results"] if fold.get("status") == "ok"]
    for result in geneformer_results
]
geneformer_folds = [fold for folds in geneformer_folds_by_task for fold in folds]
geneformer_initial_tissue_aurocs = [
    statistics.mean(fold["best_test_auroc"] for fold in folds)
    for folds in geneformer_folds_by_task
]
geneformer_epoch10_tissue_aurocs = [
    statistics.mean(fold["history"][-1]["test_auroc"] for fold in folds)
    for folds in geneformer_folds_by_task
]

baseline_paths = [
    f"evaluation/{task}_baseline_results.json" for task in TASKS
]
pca_lr_tissue_aurocs = [
    read_json(path)["pca_lr"]["mean_auroc"] for path in baseline_paths
]

thymus_path = "evaluation/submission_PCA-LR_baseline_v1_A4_eval.json"
thymus = read_json(thymus_path)

source_paths = [scgpt_path, *geneformer_paths, *baseline_paths, thymus_path]
record = {
    "averaging_unit_for_macro_aurocs": "tissue",
    "n_tissues": len(TASKS),
    "scheduled_task_folds": len(geneformer_folds),
    "scheduled_profiles": sum(fold["n_test"] for fold in geneformer_folds),
    "scgpt": {
        "covered_folds": len(scgpt_folds),
        "covered_profiles": sum(fold["n_test"] for fold in scgpt_folds),
        "initial_tissue_aurocs": scgpt_initial_tissue_aurocs,
        "initial_macro_auroc": statistics.mean(scgpt_initial_tissue_aurocs),
        "epoch10_tissue_aurocs": scgpt_epoch10_tissue_aurocs,
        "epoch10_macro_auroc": statistics.mean(scgpt_epoch10_tissue_aurocs),
    },
    "mouse_geneformer": {
        "covered_folds": len(geneformer_folds),
        "covered_profiles": sum(fold["n_test"] for fold in geneformer_folds),
        "initial_tissue_aurocs": geneformer_initial_tissue_aurocs,
        "initial_macro_auroc": statistics.mean(geneformer_initial_tissue_aurocs),
        "epoch10_tissue_aurocs": geneformer_epoch10_tissue_aurocs,
        "epoch10_macro_auroc": statistics.mean(geneformer_epoch10_tissue_aurocs),
    },
    "pca_logistic_regression": {
        "tissue_aurocs": pca_lr_tissue_aurocs,
        "macro_auroc": statistics.mean(pca_lr_tissue_aurocs),
        "thymus_mean_mission_auroc": thymus["summary"]["mean_auroc"],
        "thymus_pooled_oof_auroc": thymus["summary"]["overall_auroc_pooled"],
    },
    "source_sha256": {path: sha256(path) for path in source_paths},
}

print(json.dumps(record, indent=2, sort_keys=True))

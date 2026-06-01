#!/usr/bin/env python3
"""Build premium static visuals for the v1-v9 deck/manuscript.

The first pass intentionally rebuilds a small number of main-story figures from
canonical values instead of screenshotting legacy HTML figures.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("output/.matplotlib").resolve()))
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "output" / "premium_figures"
MANIFEST_DIR = OUT_DIR / "manifests"
MANUSCRIPT_DIR = OUT_DIR / "manuscript_variants"
TABLE_DIR = ROOT / "output" / "premium_tables"


PALETTE = {
    "paper": "#FBFAF7",
    "ink": "#1F2933",
    "muted": "#6B7280",
    "grid": "#E5E7EB",
    "chance": "#9CA3AF",
    "flight": "#E85D5A",
    "classical": "#2F6C9F",
    "evidence": "#2F6C9F",
    "pathway": "#1F9D8A",
    "biology": "#1F9D8A",
    "gene": "#3B82F6",
    "artifact": "#3B82F6",
    "highlight": "#B45309",
    "tissue_highlight": "#B45309",
    "caution": "#D99A2B",
    "low": "#A7B0BA",
    "blocker": "#B91C1C",
    "ready": "#1B8A5A",
    "draft": "#D99A2B",
    "diagnostic": "#2F6C9F",
}


TISSUE_LABELS = {
    "thymus": "Thymus",
    "gastrocnemius": "Gastrocnemius",
    "skin": "Skin",
    "eye": "Eye",
    "liver": "Liver",
    "kidney": "Kidney",
}


MANUSCRIPT_TISSUE_LABELS = {
    **TISSUE_LABELS,
    "gastrocnemius": "Gastroc.",
}


CATEGORY_B_TRANSFER = [
    {
        "tissue": "thymus",
        "auroc": 0.860,
        "ci_low": 0.763,
        "ci_high": 0.953,
        "tier": 1,
        "source": "evaluation/RESULTS_SUMMARY.md:Category B",
    },
    {
        "tissue": "gastrocnemius",
        "auroc": 0.801,
        "ci_low": 0.653,
        "ci_high": 0.944,
        "tier": 1,
        "source": "evaluation/RESULTS_SUMMARY.md:Category B",
    },
    {
        "tissue": "skin",
        "auroc": 0.772,
        "ci_low": 0.691,
        "ci_high": 0.834,
        "tier": 2,
        "source": "evaluation/RESULTS_SUMMARY.md:Category B",
    },
    {
        "tissue": "eye",
        "auroc": 0.754,
        "ci_low": 0.688,
        "ci_high": 0.838,
        "tier": 2,
        "source": "evaluation/RESULTS_SUMMARY.md:Category B",
    },
    {
        "tissue": "liver",
        "auroc": 0.577,
        "ci_low": 0.492,
        "ci_high": 0.666,
        "tier": 3,
        "source": "evaluation/RESULTS_SUMMARY.md:Category B",
    },
    {
        "tissue": "kidney",
        "auroc": 0.555,
        "ci_low": 0.397,
        "ci_high": 0.681,
        "tier": 3,
        "source": "evaluation/RESULTS_SUMMARY.md:Category B",
    },
]


CATEGORY_A_DETECTION = [
    {
        "tissue": "skin",
        "auroc": 0.821,
        "p": 0.002,
        "q": 0.012,
        "significant": True,
        "source": "evaluation/RESULTS_SUMMARY.md:Category A Detection Significance",
    },
    {
        "tissue": "gastrocnemius",
        "auroc": 0.824,
        "p": 0.026,
        "q": 0.074,
        "significant": False,
        "source": "evaluation/RESULTS_SUMMARY.md:Category A Detection Significance",
    },
    {
        "tissue": "thymus",
        "auroc": 0.923,
        "p": 0.037,
        "q": 0.074,
        "significant": False,
        "source": "evaluation/RESULTS_SUMMARY.md:Category A Detection Significance",
    },
    {
        "tissue": "eye",
        "auroc": 0.789,
        "p": 0.063,
        "q": 0.095,
        "significant": False,
        "source": "evaluation/RESULTS_SUMMARY.md:Category A Detection Significance",
    },
    {
        "tissue": "liver",
        "auroc": 0.670,
        "p": 0.091,
        "q": 0.109,
        "significant": False,
        "source": "evaluation/RESULTS_SUMMARY.md:Category A Detection Significance",
    },
    {
        "tissue": "kidney",
        "auroc": 0.432,
        "p": 0.281,
        "q": 0.281,
        "significant": False,
        "source": "evaluation/RESULTS_SUMMARY.md:Category A Detection Significance",
    },
]


GENE_PATHWAY_DETECTION = [
    {
        "tissue": "liver",
        "gene": 0.670,
        "pathway": 0.574,
        "source": "evaluation/RESULTS_SUMMARY.md:J5 Category A",
    },
    {
        "tissue": "gastrocnemius",
        "gene": 0.824,
        "pathway": 0.688,
        "source": "evaluation/RESULTS_SUMMARY.md:J5 Category A",
    },
    {
        "tissue": "kidney",
        "gene": 0.432,
        "pathway": 0.743,
        "source": "evaluation/RESULTS_SUMMARY.md:J5 Category A",
    },
    {
        "tissue": "thymus",
        "gene": 0.923,
        "pathway": 0.879,
        "source": "evaluation/RESULTS_SUMMARY.md:J5 Category A",
    },
    {
        "tissue": "eye",
        "gene": 0.789,
        "pathway": 0.915,
        "source": "evaluation/RESULTS_SUMMARY.md:J5 Category A",
    },
]


CONFOUNDER_FEATURE_COMPARISON = [
    {
        "task": "Mission ID\nliver",
        "gene": 1.000,
        "pathway": 0.056,
        "note": "D3",
        "source": "evaluation/RESULTS_SUMMARY.md:Category D",
    },
    {
        "task": "Hardware\nliver",
        "gene": 1.000,
        "pathway": 0.386,
        "note": "D5",
        "source": "evaluation/RESULTS_SUMMARY.md:Category D",
    },
    {
        "task": "Hardware\nthymus",
        "gene": 1.000,
        "pathway": 0.352,
        "note": "D5",
        "source": "evaluation/RESULTS_SUMMARY.md:Category D",
    },
    {
        "task": "Gravity\nliver",
        "gene": 0.886,
        "pathway": 0.413,
        "note": "D6",
        "source": "evaluation/RESULTS_SUMMARY.md:Category D",
    },
    {
        "task": "Gravity\nthymus",
        "gene": 0.657,
        "pathway": 0.641,
        "note": "D6",
        "source": "evaluation/RESULTS_SUMMARY.md:Category D",
    },
]


NES_TRANSFER = [
    {
        "tissue": "thymus",
        "nes_r": 0.619,
        "transfer_auroc": 0.860,
        "source": "evaluation/RESULTS_SUMMARY.md:NES Conservation vs Cross-Mission Transfer",
    },
    {
        "tissue": "eye",
        "nes_r": 0.335,
        "transfer_auroc": 0.754,
        "source": "evaluation/RESULTS_SUMMARY.md:NES Conservation vs Cross-Mission Transfer",
    },
    {
        "tissue": "skin",
        "nes_r": 0.147,
        "transfer_auroc": 0.772,
        "source": "evaluation/RESULTS_SUMMARY.md:NES Conservation vs Cross-Mission Transfer",
    },
    {
        "tissue": "liver",
        "nes_r": 0.059,
        "transfer_auroc": 0.577,
        "source": "evaluation/RESULTS_SUMMARY.md:NES Conservation vs Cross-Mission Transfer",
    },
    {
        "tissue": "gastrocnemius",
        "nes_r": 0.057,
        "transfer_auroc": 0.801,
        "source": "evaluation/RESULTS_SUMMARY.md:NES Conservation vs Cross-Mission Transfer",
    },
    {
        "tissue": "kidney",
        "nes_r": 0.048,
        "transfer_auroc": 0.555,
        "source": "evaluation/RESULTS_SUMMARY.md:NES Conservation vs Cross-Mission Transfer",
    },
]


MODEL_TIER_MEANS = [
    {
        "model": "PCA-LR",
        "score": 0.758,
        "family": "Classical ML",
        "surface": "6-tissue v1 mean",
        "source": "docs/CANONICAL_RESULTS_V7_1.md:Canonical FM / LLM Snapshot",
    },
    {
        "model": "scGPT",
        "score": 0.666,
        "family": "Single-cell foundation model",
        "surface": "6-tissue v1 mean",
        "source": "docs/CANONICAL_RESULTS_V7_1.md:Canonical FM / LLM Snapshot",
    },
    {
        "model": "Mouse-Geneformer",
        "score": 0.476,
        "family": "Single-cell foundation model",
        "surface": "6-tissue v1 mean",
        "source": "docs/CANONICAL_RESULTS_V7_1.md:Canonical FM / LLM Snapshot",
    },
    {
        "model": "DeepSeek-V3",
        "score": 0.471,
        "family": "Text LLM",
        "surface": "zero-shot, 6 tasks",
        "source": "evaluation/RESULTS_SUMMARY.md:Tier 3 LLM Zero-Shot Classification",
    },
    {
        "model": "Llama-3.3-70B",
        "score": 0.484,
        "family": "Text LLM",
        "surface": "zero-shot, 6 tasks",
        "source": "evaluation/RESULTS_SUMMARY.md:Tier 3 LLM Zero-Shot Classification",
    },
    {
        "model": "Gemini-2.5-Flash",
        "score": 0.505,
        "family": "Text LLM",
        "surface": "zero-shot, 6 tasks",
        "source": "evaluation/RESULTS_SUMMARY.md:Tier 3 LLM Zero-Shot Classification",
    },
]


MODEL_TIER_DELTAS = [
    {"tissue": "liver", "model": "scGPT", "delta": 0.040},
    {"tissue": "gastrocnemius", "model": "scGPT", "delta": -0.222},
    {"tissue": "kidney", "model": "scGPT", "delta": 0.035},
    {"tissue": "thymus", "model": "scGPT", "delta": -0.141},
    {"tissue": "skin", "model": "scGPT", "delta": -0.130},
    {"tissue": "eye", "model": "scGPT", "delta": -0.139},
    {"tissue": "liver", "model": "Mouse-Geneformer", "delta": -0.102},
    {"tissue": "gastrocnemius", "model": "Mouse-Geneformer", "delta": -0.525},
    {"tissue": "kidney", "model": "Mouse-Geneformer", "delta": -0.069},
    {"tissue": "thymus", "model": "Mouse-Geneformer", "delta": -0.428},
    {"tissue": "skin", "model": "Mouse-Geneformer", "delta": -0.264},
    {"tissue": "eye", "model": "Mouse-Geneformer", "delta": -0.305},
]


MODEL_EXTENSION_ROWS = [
    {
        "model": "scFoundation",
        "score": 0.635,
        "surface": "best single-tissue extension row",
        "source": "docs/CANONICAL_RESULTS_V7_1.md:Canonical FM / LLM Snapshot",
    },
    {
        "model": "UCE",
        "score": 0.632,
        "surface": "best single-tissue extension row",
        "source": "docs/CANONICAL_RESULTS_V7_1.md:Canonical FM / LLM Snapshot",
    },
]


V9_STATUS_ROWS = [
    {
        "lane": "Public mouse bulk",
        "provenance": "22 sources\nchecksums parsed",
        "evaluation": "8 tasks / 33 folds\n24 baseline rows",
        "release": "metadata preview\nallowed",
        "constraint": "local data copy/hash\nblocked",
        "release_kind": "ready",
        "constraint_kind": "blocked",
        "source": "v9/reports/public_bulk_alpha_snapshot_decision/snapshot_decision_summary.csv",
    },
    {
        "lane": "Human organoid",
        "provenance": "2 sources\n42 samples",
        "evaluation": "pilot baseline\ngene checks",
        "release": "draft with\nbiology checks",
        "constraint": "donor/source\ncoupling",
        "release_kind": "draft",
        "constraint_kind": "draft",
        "source": "v9/human_organoid/reports/ORGANOID_DIAGNOSTIC_CONSOLIDATION_AND_RELEASE_BOUNDARY.md",
    },
    {
        "lane": "Multispecies",
        "provenance": "3 sources\n124 samples",
        "evaluation": "plant/fly plus\nOSD-120 pilots",
        "release": "feasibility only\nno leaderboard",
        "constraint": "species namespace\nand strata",
        "release_kind": "draft",
        "constraint_kind": "draft",
        "source": "v9/multispecies/reports/MULTISPECIES_BASELINE_SENSITIVITY_REVIEW.md",
    },
    {
        "lane": "Single-cell",
        "provenance": "54 legacy assets\nindexed",
        "evaluation": "metric plan\npresent",
        "release": "no runnable\ndata yet",
        "constraint": "missing data files\n17 blockers",
        "release_kind": "blocked",
        "constraint_kind": "blocked",
        "source": "v9/sc_spaceflight/obs_var_audit_summary.csv",
    },
]


V9_FOOTPRINT_ROWS = [
    {
        "group": "Bulk",
        "metric": "public source rows",
        "value": 22,
        "display": "22",
        "kind": "ready",
        "source": "v9/source_inventory.csv",
    },
    {
        "group": "Bulk",
        "metric": "held-out mission folds",
        "value": 33,
        "display": "33",
        "kind": "ready",
        "source": "v9/reports/public_bulk_alpha_gap_matrix/public_bulk_alpha_gap_summary.csv",
    },
    {
        "group": "Bulk",
        "metric": "baseline rows",
        "value": 24,
        "display": "24",
        "kind": "ready",
        "source": "v9/reports/public_bulk_alpha_gap_matrix/public_bulk_alpha_gap_summary.csv",
    },
    {
        "group": "Organoid",
        "metric": "public samples",
        "value": 42,
        "display": "42",
        "kind": "draft",
        "source": "v9/human_organoid/sample_factors.draft.csv",
    },
    {
        "group": "Organoid",
        "metric": "gene contrasts",
        "value": 8,
        "display": "8",
        "kind": "draft",
        "source": "v9/human_organoid/de_references/human_organoid_de_reference_manifest.draft.json",
    },
    {
        "group": "Organoid",
        "metric": "gene reference rows",
        "value": 242708,
        "display": "242.7k",
        "kind": "draft",
        "source": "v9/human_organoid/de_references/human_organoid_de_reference_manifest.draft.json",
    },
    {
        "group": "Multispecies",
        "metric": "public samples",
        "value": 124,
        "display": "124",
        "kind": "draft",
        "source": "v9/multispecies/sample_factors.draft.csv",
    },
    {
        "group": "Single-cell",
        "metric": "legacy assets indexed",
        "value": 54,
        "display": "54",
        "kind": "diagnostic",
        "source": "v9/sc_spaceflight/asset_inventory_summary.csv",
    },
    {
        "group": "Single-cell",
        "metric": "data blockers",
        "value": 17,
        "display": "17",
        "kind": "blocked",
        "source": "v9/sc_spaceflight/obs_var_audit_summary.csv",
    },
]


BULK_ALPHA_GATES = [
    {
        "gate": "Task table registry",
        "status": "pass",
        "detail": "8 task tables / 33 folds indexed",
        "source": "v9/task_manifest_index.csv",
    },
    {
        "gate": "Source inventory scope",
        "status": "pass",
        "detail": "22 public mouse-bulk source rows",
        "source": "v9/source_inventory.csv",
    },
    {
        "gate": "Checksum evidence",
        "status": "pass",
        "detail": "22/22 sources have parsed OSDR checksum-manifest evidence",
        "source": "v9/source_checksum_audit.csv",
    },
    {
        "gate": "Baseline output evidence",
        "status": "pass",
        "detail": "24/24 simple baseline rows evaluated",
        "source": "v9/reports/bulk_lomo_baseline_summary.csv",
    },
    {
        "gate": "Data Package boundary",
        "status": "needs update",
        "detail": "metadata and data-file evidence need linked descriptors",
        "source": "v9/datapackage.draft.json",
    },
    {
        "gate": "Dataset description",
        "status": "needs update",
        "detail": "description should separate indexed metadata from data files",
        "source": "docs/v9_hf_dataset_card.md",
    },
    {
        "gate": "Data-file hash verification",
        "status": "pending",
        "detail": "0/22 sources have local data-file hash evidence",
        "source": "v9/source_checksum_audit.csv",
    },
]


BULK_ALPHA_DECISIONS = [
    {
        "path": "Metadata records",
        "status": "available",
        "detail": "source, task, fold, checksum, and baseline metadata indexed",
        "kind": "ready",
        "source": "v9/reports/public_bulk_alpha_snapshot_decision/snapshot_option_matrix.csv",
    },
    {
        "path": "Data-file mirror",
        "status": "pending",
        "detail": "local file copy and hash verification still separate",
        "kind": "diagnostic",
        "source": "v9/reports/public_bulk_alpha_snapshot_decision/snapshot_option_matrix.csv",
    },
    {
        "path": "Integrated package",
        "status": "future",
        "detail": "metadata and data-file evidence can be joined after verification",
        "kind": "diagnostic",
        "source": "v9/reports/public_bulk_alpha_snapshot_decision/snapshot_option_matrix.csv",
    },
]


BULK_ALPHA_LANGUAGE = [
    {
        "side": "Metadata evidence",
        "text": "v9 public bulk metadata records",
        "kind": "ready",
        "source": "v9/reports/public_bulk_alpha_snapshot_decision/snapshot_claim_boundary.csv",
    },
    {
        "side": "Metadata evidence",
        "text": "Task/source/checksum/baseline metadata",
        "kind": "ready",
        "source": "v9/reports/public_bulk_alpha_snapshot_decision/snapshot_claim_boundary.csv",
    },
    {
        "side": "Metadata evidence",
        "text": "Baseline outputs as reproducibility checks",
        "kind": "ready",
        "source": "v9/reports/public_bulk_alpha_snapshot_decision/snapshot_claim_boundary.csv",
    },
    {
        "side": "Data-file evidence",
        "text": "Local data-file mirror",
        "kind": "diagnostic",
        "source": "v9/reports/public_bulk_alpha_snapshot_decision/snapshot_claim_boundary.csv",
    },
    {
        "side": "Data-file evidence",
        "text": "Local hash-verified data bundle",
        "kind": "diagnostic",
        "source": "v9/reports/public_bulk_alpha_snapshot_decision/snapshot_claim_boundary.csv",
    },
    {
        "side": "Data-file evidence",
        "text": "Complete public benchmark package",
        "kind": "diagnostic",
        "source": "v9/reports/public_bulk_alpha_snapshot_decision/snapshot_claim_boundary.csv",
    },
]


ORGANOID_FOOTPRINT_ROWS = [
    {
        "metric": "public sources",
        "value": 2,
        "display": "2",
        "detail": "OSD-863 / OSD-871",
        "kind": "draft",
        "source": "v9/human_organoid/source_inventory.draft.csv",
    },
    {
        "metric": "samples",
        "value": 42,
        "display": "42",
        "detail": "public neural organoid RNA-seq",
        "kind": "draft",
        "source": "v9/human_organoid/sample_factors.draft.csv",
    },
    {
        "metric": "flight-ground gene contrasts",
        "value": 8,
        "display": "8",
        "detail": "Ground vs LEO/ISS",
        "kind": "diagnostic",
        "source": "v9/human_organoid/de_references/human_organoid_de_reference_manifest.draft.json",
    },
    {
        "metric": "gene reference rows",
        "value": 242708,
        "display": "242.7k",
        "detail": "gene/contrast rows",
        "kind": "diagnostic",
        "source": "v9/human_organoid/de_references/human_organoid_de_reference_manifest.draft.json",
    },
    {
        "metric": "significant gene rows",
        "value": 2368,
        "display": "2,368",
        "detail": "significant reference rows",
        "kind": "diagnostic",
        "source": "v9/human_organoid/de_references/human_organoid_de_reference_manifest.draft.json",
    },
    {
        "metric": "pilot status",
        "value": 0,
        "display": "draft",
        "detail": "not frozen; diagnostics only",
        "kind": "blocked",
        "source": "v9/human_organoid/reports/ORGANOID_DIAGNOSTIC_CONSOLIDATION_AND_RELEASE_BOUNDARY.md",
    },
]


ORGANOID_METRIC_ROWS = [
    {
        "family": "Primary task",
        "metric": "AUROC",
        "value": 0.6147727273,
        "role": "prediction",
        "kind": "ready",
        "source": "v9/human_organoid/reports/nearest_centroid/human_organoid_baseline_summary.csv",
    },
    {
        "family": "Primary task",
        "metric": "Macro-F1",
        "value": 0.5194508009,
        "role": "prediction",
        "kind": "ready",
        "source": "v9/human_organoid/reports/nearest_centroid/human_organoid_baseline_summary.csv",
    },
    {
        "family": "Flight-response pattern",
        "metric": "Direction match",
        "value": 0.7706734868,
        "role": "biology check",
        "kind": "diagnostic",
        "source": "v9/human_organoid/reports/source_transfer_signature/metrics.json",
    },
    {
        "family": "Flight-response pattern",
        "metric": "Rank correlation",
        "value": 0.1760078660,
        "role": "biology check",
        "kind": "diagnostic",
        "source": "v9/human_organoid/reports/source_transfer_signature/metrics.json",
    },
    {
        "family": "Model gene effects",
        "metric": "Direction match",
        "value": 0.6078431373,
        "role": "secondary check",
        "kind": "draft",
        "source": "v9/human_organoid/reports/logistic_feature_effect/metrics.json",
    },
    {
        "family": "Model gene effects",
        "metric": "Rank correlation",
        "value": 0.0867280024,
        "role": "secondary check",
        "kind": "draft",
        "source": "v9/human_organoid/reports/logistic_feature_effect/metrics.json",
    },
]


ORGANOID_TOPK_ROWS = [
    {
        "k": 50,
        "observed": 1,
        "expected": 0.6375,
        "enrichment": 1.5686,
        "p_value": 0.474071,
        "source": "v9/human_organoid/reports/ORGANOID_FEATURE_EFFECT_NULL_CALIBRATION_REVIEW.md",
    },
    {
        "k": 100,
        "observed": 5,
        "expected": 1.2750,
        "enrichment": 3.9216,
        "p_value": 0.009090,
        "source": "v9/human_organoid/reports/ORGANOID_FEATURE_EFFECT_NULL_CALIBRATION_REVIEW.md",
    },
    {
        "k": 250,
        "observed": 10,
        "expected": 3.1875,
        "enrichment": 3.1373,
        "p_value": 0.001414,
        "source": "v9/human_organoid/reports/ORGANOID_FEATURE_EFFECT_NULL_CALIBRATION_REVIEW.md",
    },
    {
        "k": 500,
        "observed": 14,
        "expected": 6.3750,
        "enrichment": 2.1961,
        "p_value": 0.004966,
        "source": "v9/human_organoid/reports/ORGANOID_FEATURE_EFFECT_NULL_CALIBRATION_REVIEW.md",
    },
]


ORGANOID_DECISION_ROWS = [
    {
        "family": "Task definition",
        "artifact": "sample table and source metadata",
        "decision": "pilot task",
        "kind": "draft",
        "source": "v9/human_organoid/reports/ORGANOID_DIAGNOSTIC_CONSOLIDATION_AND_RELEASE_BOUNDARY.md",
    },
    {
        "family": "Flight-response pattern",
        "artifact": "cross-source response genes",
        "decision": "show by default",
        "kind": "ready",
        "source": "v9/human_organoid/reports/source_transfer_signature/metrics.json",
    },
    {
        "family": "Model gene effects",
        "artifact": "model gene-weight check",
        "decision": "show as secondary",
        "kind": "draft",
        "source": "v9/human_organoid/reports/logistic_feature_effect/metrics.json",
    },
    {
        "family": "Microglia-matched pattern",
        "artifact": "microglia-aware source check",
        "decision": "keep secondary",
        "kind": "diagnostic",
        "source": "v9/human_organoid/reports/microglia_source_transfer_signature/metrics.json",
    },
    {
        "family": "Shared-control subset",
        "artifact": "partial disease and microglia match",
        "decision": "partial/negative",
        "kind": "blocked",
        "source": "v9/human_organoid/reports/shared_control_source_transfer_signature/metrics.json",
    },
    {
        "family": "PCA-LR gene effects",
        "artifact": "reconstructed gene weights",
        "decision": "negative comparison",
        "kind": "blocked",
        "source": "v9/human_organoid/reports/pca_lr_feature_effect/metrics.json",
    },
]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": PALETTE["paper"],
            "axes.facecolor": "white",
            "axes.edgecolor": "#D1D5DB",
            "axes.labelcolor": PALETTE["ink"],
            "xtick.color": PALETTE["ink"],
            "ytick.color": PALETTE["ink"],
            "text.color": PALETTE["ink"],
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titleweight": "bold",
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def write_manifest(rows: list[dict[str, object]], stem: str) -> None:
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
    json_path = MANIFEST_DIR / f"{stem}.json"
    csv_path = MANIFEST_DIR / f"{stem}.csv"
    with json_path.open("w", encoding="utf-8") as fh:
        json.dump(rows, fh, indent=2)
        fh.write("\n")
    keys = sorted({key for row in rows for key in row})
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def write_table(rows: list[dict[str, object]], stem: str) -> list[Path]:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    json_path = TABLE_DIR / f"{stem}.json"
    csv_path = TABLE_DIR / f"{stem}.csv"
    with json_path.open("w", encoding="utf-8") as fh:
        json.dump(rows, fh, indent=2)
        fh.write("\n")
    keys = sorted({key for row in rows for key in row})
    with csv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)
    return [csv_path, json_path]


def strip_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.8)
    ax.set_axisbelow(True)


def draw_tissue_transfer_panel(ax: plt.Axes) -> None:
    rows = sorted(CATEGORY_B_TRANSFER, key=lambda row: row["auroc"], reverse=True)
    y = np.arange(len(rows))
    values = np.array([float(row["auroc"]) for row in rows])
    lows = np.array([float(row["ci_low"]) for row in rows])
    highs = np.array([float(row["ci_high"]) for row in rows])
    tier_style = {
        1: {"marker": "o", "label": "High transfer", "color": PALETTE["tissue_highlight"]},
        2: {"marker": "D", "label": "Mid-range", "color": PALETTE["evidence"]},
        3: {"marker": "s", "label": "Near chance", "color": PALETTE["low"]},
    }
    ax.hlines(y, lows, highs, color="#CBD5E1", linewidth=2.0, zorder=1)
    seen_tiers: set[int] = set()
    for yi, value, row in zip(y, values, rows):
        tier = int(row["tier"])
        style = tier_style[tier]
        label = str(style["label"]) if tier not in seen_tiers else None
        ax.scatter(
            value,
            yi,
            s=88 if tier == 1 else 78,
            color=str(style["color"]),
            marker=str(style["marker"]),
            edgecolor="white",
            linewidth=0.9,
            label=label,
            zorder=3,
        )
        seen_tiers.add(tier)
    for yi, row in zip(y, rows):
        ax.text(
            float(row["auroc"]) + 0.012,
            yi,
            f"{float(row['auroc']):.3f}",
            va="center",
            ha="left",
            fontsize=9,
            fontweight="bold" if row["tissue"] in {"thymus", "liver"} else "normal",
        )
    ax.axvline(0.5, color=PALETTE["chance"], linestyle="--", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels([TISSUE_LABELS[str(row["tissue"])] for row in rows])
    ax.set_xlim(0.35, 1.02)
    ax.set_ylim(len(rows) - 0.45, -0.55)
    ax.set_xlabel("Held-out mission AUROC")
    ax.set_title("Held-out mission prediction by tissue", loc="left")
    ax.hlines([1.5, 3.5], 0.35, 1.02, color="#E5E7EB", linewidth=0.9, zorder=0)
    ax.text(0.37, 0.5, "High", fontsize=8, color=PALETTE["muted"], va="center")
    ax.text(0.37, 2.5, "Mid", fontsize=8, color=PALETTE["muted"], va="center")
    ax.text(0.37, 4.5, "Near chance", fontsize=8, color=PALETTE["muted"], va="center")
    ax.text(
        0.54,
        0.08,
        "Thymus - liver: +0.283 AUROC (p=0.001)",
        fontsize=9,
        color=PALETTE["tissue_highlight"],
        fontweight="bold",
        transform=ax.transAxes,
    )
    ax.legend(loc="lower right", fontsize=8, ncol=3, handletextpad=0.35, columnspacing=0.9)
    strip_axes(ax)


def draw_tissue_claim_panel(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(0.0, 1.02, "B. Main readout", fontsize=12, fontweight="bold", transform=ax.transAxes)
    rows = [
        ("Highest transfer", "Thymus", "0.860 AUROC"),
        ("Strong transfer", "Gastrocnemius", "0.801 AUROC"),
        ("Mid-range", "Skin / eye", "0.772 / 0.754"),
        ("Near chance", "Liver / kidney", "0.577 / 0.555"),
    ]
    header_y = 0.87
    ax.text(0.00, header_y, "Observation", fontsize=8.4, fontweight="bold", transform=ax.transAxes)
    ax.text(0.45, header_y, "Tissue", fontsize=8.4, fontweight="bold", transform=ax.transAxes)
    ax.text(0.76, header_y, "Value", fontsize=8.4, fontweight="bold", transform=ax.transAxes)
    ax.hlines(header_y - 0.055, 0.0, 0.98, color="#D1D5DB", linewidth=1.0, transform=ax.transAxes)
    for i, (observation, tissue, value) in enumerate(rows):
        y = 0.73 - i * 0.13
        color = PALETTE["tissue_highlight"] if i == 0 else PALETTE["low"] if i == 3 else PALETTE["evidence"]
        ax.hlines(y - 0.055, 0.0, 0.98, color="#E5E7EB", linewidth=0.8, transform=ax.transAxes)
        ax.scatter(0.018, y, s=48, color=color, transform=ax.transAxes)
        ax.text(0.05, y, observation, fontsize=8.6, va="center", transform=ax.transAxes)
        ax.text(0.45, y, tissue, fontsize=8.6, va="center", fontweight="bold", transform=ax.transAxes)
        ax.text(0.76, y, value, fontsize=8.6, va="center", color=color, fontweight="bold", transform=ax.transAxes)
    ax.text(
        0.0,
        0.18,
        "Interpretation: transfer is tissue-specific.\nScope: mouse bulk RNA-seq; not a universal model or tissue claim.",
        fontsize=8.8,
        color=PALETTE["ink"],
        transform=ax.transAxes,
        linespacing=1.28,
    )


def tissue_transfer_summary_table() -> list[dict[str, object]]:
    return [
        {
            "observation": "highest_transfer",
            "tissue": "thymus",
            "auroc": 0.860,
            "ci_low": 0.763,
            "ci_high": 0.953,
            "display": "Thymus 0.860 AUROC",
            "source": "evaluation/RESULTS_SUMMARY.md:Category B",
        },
        {
            "observation": "strong_transfer",
            "tissue": "gastrocnemius",
            "auroc": 0.801,
            "ci_low": 0.653,
            "ci_high": 0.944,
            "display": "Gastrocnemius 0.801 AUROC",
            "source": "evaluation/RESULTS_SUMMARY.md:Category B",
        },
        {
            "observation": "mid_range_transfer",
            "tissue": "skin_eye",
            "auroc": "0.772 / 0.754",
            "ci_low": "0.691 / 0.688",
            "ci_high": "0.834 / 0.838",
            "display": "Skin / eye 0.772 / 0.754",
            "source": "evaluation/RESULTS_SUMMARY.md:Category B",
        },
        {
            "observation": "near_chance_transfer",
            "tissue": "liver_kidney",
            "auroc": "0.577 / 0.555",
            "ci_low": "0.492 / 0.397",
            "ci_high": "0.666 / 0.681",
            "display": "Liver / kidney 0.577 / 0.555",
            "source": "evaluation/RESULTS_SUMMARY.md:Category B",
        },
        {
            "observation": "lead_difference",
            "tissue": "thymus_minus_liver",
            "auroc": 0.283,
            "ci_low": "",
            "ci_high": "",
            "display": "Thymus - liver +0.283 AUROC, p=0.001",
            "source": "evaluation/RESULTS_SUMMARY.md:Category B",
        },
    ]


def draw_detection_panel(ax: plt.Axes) -> None:
    rows = sorted(CATEGORY_A_DETECTION, key=lambda row: float(row["auroc"]), reverse=True)
    x = np.arange(len(rows))
    colors = [
        PALETTE["highlight"] if row["tissue"] == "thymus" else
        PALETTE["flight"] if bool(row["significant"]) else
        PALETTE["classical"]
        for row in rows
    ]
    values = [float(row["auroc"]) for row in rows]
    ax.bar(x, values, color=colors, width=0.68)
    ax.axhline(0.5, color=PALETTE["chance"], linestyle="--", linewidth=1.0)
    ax.axhline(0.7, color=PALETTE["grid"], linestyle="-", linewidth=1.0)
    for xi, row in zip(x, rows):
        ax.text(
            xi,
            float(row["auroc"]) + 0.025,
            f"{float(row['auroc']):.3f}",
            ha="center",
            va="bottom",
            fontsize=8,
            fontweight="bold" if row["tissue"] in {"thymus", "skin"} else "normal",
        )
        if bool(row["significant"]):
            ax.text(xi, 0.97, "FDR<0.05", ha="center", va="top", fontsize=7, color=PALETTE["flight"])
    ax.set_xticks(x)
    ax.set_xticklabels([TISSUE_LABELS[str(row["tissue"])] for row in rows], rotation=35, ha="right")
    ax.set_ylim(0.35, 1.02)
    ax.set_ylabel("Held-out mission AUROC")
    ax.set_title("B. Detection is strong but significance is fold-limited", loc="left")
    ax.grid(axis="y", color=PALETTE["grid"], linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def draw_pathway_rescue_panel(ax: plt.Axes) -> None:
    rows = sorted(
        GENE_PATHWAY_DETECTION,
        key=lambda row: float(row["pathway"]) - float(row["gene"]),
    )
    y = np.arange(len(rows))
    gene = np.array([float(row["gene"]) for row in rows])
    pathway = np.array([float(row["pathway"]) for row in rows])
    ax.hlines(y, gene, pathway, color=PALETTE["muted"], linewidth=1.6, alpha=0.75)
    ax.scatter(gene, y, s=56, color=PALETTE["gene"], label="Gene")
    ax.scatter(pathway, y, s=66, color=PALETTE["pathway"], label="Pathway")
    for yi, row in zip(y, rows):
        delta = float(row["pathway"]) - float(row["gene"])
        color = PALETTE["pathway"] if delta > 0 else PALETTE["muted"]
        label = f"{delta:+.3f}"
        ax.text(max(float(row["gene"]), float(row["pathway"])) + 0.018, yi, label, va="center", fontsize=9, color=color)
    ax.axvline(0.5, color=PALETTE["chance"], linestyle="--", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels([TISSUE_LABELS[str(row["tissue"])] for row in rows])
    ax.set_xlim(0.35, 1.02)
    ax.set_xlabel("Detection AUROC")
    ax.set_title("C. Pathway summaries help selected weak gene tasks", loc="left")
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.18), ncol=2, fontsize=8)
    strip_axes(ax)


def build_core_figure() -> list[Path]:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(14.8, 8.3), constrained_layout=False)
    ax_a = fig.add_axes([0.12, 0.16, 0.83, 0.64])
    draw_tissue_transfer_panel(ax_a)
    fig.suptitle(
        "Some mouse tissues generalize across missions better than others",
        x=0.12,
        y=0.955,
        ha="left",
        fontsize=20,
        fontweight="bold",
    )
    fig.text(
        0.12,
        0.905,
        "Six mouse bulk RNA-seq tasks tested by holding out one mission at a time. Chance AUROC = 0.5; intervals show bootstrap uncertainty.",
        ha="left",
        fontsize=11,
        color=PALETTE["muted"],
    )
    outputs = [
        OUT_DIR / "premium_fig1_tissue_transfer_hierarchy.png",
        OUT_DIR / "premium_fig1_tissue_transfer_hierarchy.pdf",
    ]
    fig.savefig(outputs[0], dpi=220)
    fig.savefig(outputs[1])
    legacy_outputs = [
        OUT_DIR / "premium_fig1_core_tissue_pathway.png",
        OUT_DIR / "premium_fig1_core_tissue_pathway.pdf",
    ]
    fig.savefig(legacy_outputs[0], dpi=220)
    fig.savefig(legacy_outputs[1])
    plt.close(fig)
    manifest_rows: list[dict[str, object]] = []
    for row in CATEGORY_B_TRANSFER:
        item = dict(row)
        item["panel"] = "held_out_mission_tissue_transfer"
        manifest_rows.append(item)
    manifest_rows.extend(
        [
            {
                "panel": "claim_boundary",
                "claim": "transfer_is_tissue_specific",
                "scope": "mouse_bulk_rna_seq_six_tissue_tasks",
                "source": "evaluation/RESULTS_SUMMARY.md:Category B",
            }
        ]
    )
    write_manifest(manifest_rows, "premium_fig1_tissue_transfer_hierarchy_manifest")
    write_manifest(manifest_rows, "premium_fig1_core_tissue_pathway_manifest")
    table_outputs = write_table(tissue_transfer_summary_table(), "table1_tissue_transfer_summary")
    return outputs + legacy_outputs + table_outputs


def draw_confounder_panel(ax: plt.Axes) -> None:
    rows = list(reversed(CONFOUNDER_FEATURE_COMPARISON))
    y = np.arange(len(rows))
    height = 0.34
    gene = np.array([float(row["gene"]) for row in rows])
    pathway = np.array([float(row["pathway"]) for row in rows])
    ax.barh(y + height / 2, gene, height=height, color=PALETTE["gene"], label="Gene")
    ax.barh(y - height / 2, pathway, height=height, color=PALETTE["pathway"], label="Pathway")
    for yi, row in zip(y, rows):
        ax.text(1.02, yi, str(row["note"]), va="center", ha="left", fontsize=8, color=PALETTE["muted"])
    ax.axvline(0.5, color=PALETTE["chance"], linestyle="--", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels([str(row["task"]) for row in rows])
    ax.set_xlim(0, 1.12)
    ax.set_xlabel("Prediction score (macro-F1)")
    ax.set_title("A. Gene inputs expose mission labels; pathway summaries reduce them", loc="left")
    ax.legend(loc="lower right", ncol=2, fontsize=8)
    strip_axes(ax)


def draw_rescue_delta_panel(ax: plt.Axes) -> None:
    rows = sorted(
        GENE_PATHWAY_DETECTION,
        key=lambda row: float(row["pathway"]) - float(row["gene"]),
        reverse=True,
    )
    x = np.arange(len(rows))
    deltas = np.array([float(row["pathway"]) - float(row["gene"]) for row in rows])
    colors = [PALETTE["pathway"] if value > 0 else PALETTE["low"] for value in deltas]
    ax.bar(x, deltas, color=colors, width=0.68)
    ax.axhline(0, color=PALETTE["ink"], linewidth=1.0)
    for xi, value in zip(x, deltas):
        va = "bottom" if value >= 0 else "top"
        offset = 0.018 if value >= 0 else -0.018
        ax.text(xi, value + offset, f"{value:+.3f}", ha="center", va=va, fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([TISSUE_LABELS[str(row["tissue"])] for row in rows], rotation=35, ha="right")
    ax.set_ylabel("Pathway - gene AUROC")
    ax.set_ylim(-0.18, 0.36)
    ax.set_title("B. Selected weak tasks improve", loc="left")
    ax.grid(axis="y", color=PALETTE["grid"], linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def draw_nes_transfer_panel(ax: plt.Axes) -> None:
    rows = NES_TRANSFER
    x = np.array([float(row["nes_r"]) for row in rows])
    y = np.array([float(row["transfer_auroc"]) for row in rows])
    colors = [
        PALETTE["tissue_highlight"] if row["tissue"] == "thymus" else
        PALETTE["low"] if row["tissue"] == "gastrocnemius" else
        PALETTE["evidence"]
        for row in rows
    ]
    ax.scatter(x, y, s=90, color=colors, edgecolor="white", linewidth=1.0, zorder=3)
    for row in rows:
        tx = float(row["nes_r"])
        ty = float(row["transfer_auroc"])
        label = TISSUE_LABELS[str(row["tissue"])]
        dx = 0.012
        dy = 0.012
        if row["tissue"] == "gastrocnemius":
            dy = -0.035
        ax.text(tx + dx, ty + dy, label, fontsize=8, color=PALETTE["ink"])
    ax.axhline(0.5, color=PALETTE["chance"], linestyle="--", linewidth=1.0)
    ax.set_xlim(0.0, 0.68)
    ax.set_ylim(0.48, 0.91)
    ax.set_xlabel("Pathway activity agreement")
    ax.set_ylabel("Held-out mission AUROC")
    ax.set_title("C. Activity agreement across tissues", loc="left")
    ax.text(
        0.04,
        0.875,
        "Six tissues; gastrocnemius remains discordant",
        fontsize=9,
        color=PALETTE["muted"],
    )
    ax.grid(color=PALETTE["grid"], linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def build_pathway_figure() -> list[Path]:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(16, 9), constrained_layout=False)
    grid = fig.add_gridspec(
        nrows=2,
        ncols=2,
        height_ratios=[1.0, 0.95],
        width_ratios=[1.05, 0.95],
        left=0.075,
        right=0.965,
        top=0.82,
        bottom=0.13,
        wspace=0.28,
        hspace=0.38,
    )
    ax_a = fig.add_subplot(grid[:, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 1])
    draw_confounder_panel(ax_a)
    draw_rescue_delta_panel(ax_b)
    draw_nes_transfer_panel(ax_c)
    fig.suptitle(
        "Pathway summaries reduce unwanted mission signals and help selected weak tasks",
        x=0.075,
        y=0.955,
        ha="left",
        fontsize=20,
        fontweight="bold",
    )
    fig.text(
        0.075,
        0.905,
        "Gene-level inputs can identify mission and hardware labels; pathway summaries often dampen those signals, but benefits are tissue-specific.",
        ha="left",
        fontsize=11,
        color=PALETTE["muted"],
    )
    outputs = [
        OUT_DIR / "premium_fig2_pathway_artifact_rescue.png",
        OUT_DIR / "premium_fig2_pathway_artifact_rescue.pdf",
    ]
    fig.savefig(outputs[0], dpi=220)
    fig.savefig(outputs[1])
    plt.close(fig)
    manifest_rows: list[dict[str, object]] = []
    for label, rows in [
        ("confounder_feature_comparison", CONFOUNDER_FEATURE_COMPARISON),
        ("gene_pathway_detection_delta", GENE_PATHWAY_DETECTION),
        ("nes_transfer", NES_TRANSFER),
    ]:
        for row in rows:
            item = dict(row)
            item["panel"] = label
            manifest_rows.append(item)
    write_manifest(manifest_rows, "premium_fig2_pathway_artifact_rescue_manifest")
    return outputs


def draw_confounder_lollipop_panel(ax: plt.Axes) -> None:
    rows = CONFOUNDER_FEATURE_COMPARISON
    y = np.arange(len(rows))
    gene = np.array([float(row["gene"]) for row in rows])
    pathway = np.array([float(row["pathway"]) for row in rows])
    ax.hlines(y, pathway, gene, color="#CBD5E1", linewidth=2.2, zorder=1)
    ax.scatter(gene, y, s=54, color=PALETTE["gene"], edgecolor="white", linewidth=0.8, label="Gene", zorder=3)
    ax.scatter(pathway, y, s=62, color=PALETTE["pathway"], edgecolor="white", linewidth=0.8, label="Pathway", zorder=3)
    for yi, row in zip(y, rows):
        ax.text(1.025, yi, str(row["note"]), va="center", ha="left", fontsize=7.3, color=PALETTE["muted"])
    ax.axvline(0.5, color=PALETTE["chance"], linestyle="--", linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([str(row["task"]) for row in rows])
    ax.set_xlim(0, 1.08)
    ax.set_ylim(len(rows) - 0.45, -0.55)
    ax.set_xlabel("Macro-F1 for coupled-label prediction")
    ax.set_title("A. Unwanted label signals fall in pathway summaries", loc="left")
    strip_axes(ax)


def draw_rescue_delta_lollipop_panel(ax: plt.Axes) -> None:
    rows = sorted(
        GENE_PATHWAY_DETECTION,
        key=lambda row: float(row["pathway"]) - float(row["gene"]),
    )
    y = np.arange(len(rows))
    deltas = np.array([float(row["pathway"]) - float(row["gene"]) for row in rows])
    colors = [PALETTE["pathway"] if value > 0 else PALETTE["low"] for value in deltas]
    ax.axvline(0, color=PALETTE["ink"], linewidth=0.9)
    for yi, value, color in zip(y, deltas, colors):
        ax.hlines(yi, min(0, value), max(0, value), color="#CBD5E1", linewidth=2.0, zorder=1)
        ax.scatter(value, yi, s=56, color=color, edgecolor="white", linewidth=0.8, zorder=3)
        ha = "left"
        dx = 0.025
        ax.text(value + dx, yi, f"{value:+.3f}", va="center", ha=ha, fontsize=7.4, color=PALETTE["ink"])
    ax.set_yticks(y)
    ax.set_yticklabels([MANUSCRIPT_TISSUE_LABELS[str(row["tissue"])] for row in rows])
    ax.set_xlim(-0.18, 0.35)
    ax.set_xlabel("Pathway - gene AUROC")
    ax.set_title("B. Detection-task change", loc="left")
    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def draw_nes_transfer_manuscript_panel(ax: plt.Axes) -> None:
    rows = NES_TRANSFER
    x = np.array([float(row["nes_r"]) for row in rows])
    y = np.array([float(row["transfer_auroc"]) for row in rows])
    colors = [
        PALETTE["tissue_highlight"] if row["tissue"] == "thymus" else
        PALETTE["low"] if row["tissue"] == "gastrocnemius" else
        PALETTE["evidence"]
        for row in rows
    ]
    ax.scatter(x, y, s=58, color=colors, edgecolor="white", linewidth=0.8, zorder=3)
    label_offsets = {
        "thymus": (-0.095, -0.018),
        "eye": (0.014, -0.005),
        "skin": (0.020, 0.022),
        "liver": (0.014, -0.005),
        "kidney": (0.014, -0.010),
        "gastrocnemius": (0.014, -0.045),
    }
    for row in rows:
        tx = float(row["nes_r"])
        ty = float(row["transfer_auroc"])
        dx, dy = label_offsets[str(row["tissue"])]
        ax.text(tx + dx, ty + dy, MANUSCRIPT_TISSUE_LABELS[str(row["tissue"])], fontsize=7.2, color=PALETTE["ink"])
    ax.axhline(0.5, color=PALETTE["chance"], linestyle="--", linewidth=0.9)
    ax.set_xlim(0.0, 0.68)
    ax.set_ylim(0.48, 0.90)
    ax.set_xlabel("Pathway activity agreement")
    ax.set_ylabel("Held-out mission AUROC")
    ax.set_title("C. Activity agreement check", loc="left")
    ax.text(0.04, 0.855, "No fitted trend; n=6 tissues", fontsize=7.4, color=PALETTE["muted"])
    ax.grid(color=PALETTE["grid"], linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def build_pathway_manuscript_figure() -> list[Path]:
    setup_style()
    MANUSCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.2, 5.65), constrained_layout=False)
    grid = fig.add_gridspec(
        nrows=2,
        ncols=2,
        height_ratios=[1.18, 1.0],
        width_ratios=[0.96, 1.04],
        left=0.16,
        right=0.975,
        top=0.79,
        bottom=0.15,
        wspace=0.46,
        hspace=0.64,
    )
    ax_a = fig.add_subplot(grid[0, :])
    ax_b = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[1, 1])
    draw_confounder_lollipop_panel(ax_a)
    draw_rescue_delta_lollipop_panel(ax_b)
    draw_nes_transfer_manuscript_panel(ax_c)
    fig.suptitle(
        "Pathway summaries reduce unwanted label signals",
        x=0.16,
        y=0.962,
        ha="left",
        fontsize=13.2,
        fontweight="bold",
    )
    fig.text(
        0.16,
        0.915,
        "Diagnostic label prediction falls for pathway inputs; selected prediction gains remain tissue-specific.",
        ha="left",
        fontsize=8.6,
        color=PALETTE["muted"],
    )
    fig.text(
        0.70,
        0.845,
        "Gene",
        ha="left",
        fontsize=7.4,
        color=PALETTE["gene"],
        fontweight="bold",
    )
    fig.text(
        0.77,
        0.845,
        "Pathway",
        ha="left",
        fontsize=7.4,
        color=PALETTE["pathway"],
        fontweight="bold",
    )
    outputs = [
        MANUSCRIPT_DIR / "premium_fig2_pathway_rescue_manuscript.png",
        MANUSCRIPT_DIR / "premium_fig2_pathway_rescue_manuscript.pdf",
    ]
    fig.savefig(outputs[0], dpi=320)
    fig.savefig(outputs[1])
    plt.close(fig)
    manifest_rows: list[dict[str, object]] = []
    for label, rows in [
        ("confounder_feature_comparison", CONFOUNDER_FEATURE_COMPARISON),
        ("gene_pathway_detection_delta", GENE_PATHWAY_DETECTION),
        ("pathway_activity_agreement", NES_TRANSFER),
    ]:
        for row in rows:
            item = dict(row)
            item["panel"] = label
            item["variant"] = "manuscript"
            manifest_rows.append(item)
    write_manifest(manifest_rows, "premium_fig2_pathway_rescue_manuscript_manifest")
    return outputs


def family_color(family: str) -> str:
    if family == "Classical ML":
        return PALETTE["classical"]
    if family == "Single-cell foundation model":
        return "#7C3AED"
    if family == "Text LLM":
        return "#C026D3"
    return PALETTE["muted"]


def draw_model_mean_panel(ax: plt.Axes) -> None:
    rows = sorted(MODEL_TIER_MEANS, key=lambda row: float(row["score"]))
    y = np.arange(len(rows))
    values = [float(row["score"]) for row in rows]
    colors = [family_color(str(row["family"])) for row in rows]
    ax.barh(y, values, color=colors, height=0.62)
    ax.axvline(0.5, color=PALETTE["chance"], linestyle="--", linewidth=1.0)
    for yi, row in zip(y, rows):
        ax.text(float(row["score"]) + 0.012, yi, f"{float(row['score']):.3f}", va="center", fontsize=9)
    ax.set_yticks(y)
    ax.set_yticklabels([str(row["model"]) for row in rows])
    ax.set_xlim(0.4, 0.82)
    ax.set_xlabel("Mean AUROC / score")
    ax.set_title("A. Classical baseline is strongest across the shared 6-task set", loc="left")
    strip_axes(ax)
    legend_items = [
        ("Classical ML", family_color("Classical ML")),
        ("Single-cell model", family_color("Single-cell foundation model")),
        ("Text LLM", family_color("Text LLM")),
    ]
    for i, (name, color) in enumerate(legend_items):
        ax.scatter(0.42 + i * 0.13, len(rows) - 0.15, s=60, color=color)
        ax.text(0.435 + i * 0.13, len(rows) - 0.15, name, va="center", fontsize=8)


def draw_model_delta_panel(ax: plt.Axes) -> None:
    tissues = ["liver", "gastrocnemius", "kidney", "thymus", "skin", "eye"]
    y = np.arange(len(tissues))
    scgpt = {row["tissue"]: float(row["delta"]) for row in MODEL_TIER_DELTAS if row["model"] == "scGPT"}
    geneformer = {
        row["tissue"]: float(row["delta"])
        for row in MODEL_TIER_DELTAS
        if row["model"] == "Mouse-Geneformer"
    }
    ax.axvline(0, color=PALETTE["ink"], linewidth=1.0)
    ax.scatter([scgpt[t] for t in tissues], y + 0.13, s=64, color="#7C3AED", label="scGPT")
    ax.scatter([geneformer[t] for t in tissues], y - 0.13, s=64, color="#A78BFA", label="Mouse-Geneformer")
    for yi, tissue in zip(y, tissues):
        ax.hlines(yi, min(scgpt[tissue], geneformer[tissue], 0), max(scgpt[tissue], geneformer[tissue], 0), color=PALETTE["grid"], linewidth=2)
    ax.set_yticks(y)
    ax.set_yticklabels([TISSUE_LABELS[t] for t in tissues])
    ax.set_xlim(-0.58, 0.08)
    ax.set_xlabel("AUROC change vs matched classical baseline")
    ax.set_title("B. Foundation-model gains are local; most tissues are lower", loc="left")
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def draw_model_boundary_panel(ax: plt.Axes) -> None:
    ax.axis("off")
    rows = [
        ("Shared 6-task mean", "PCA-LR, scGPT, Mouse-Geneformer", "same task set"),
        ("Text-only tasks", "DeepSeek, Llama, Gemini", "zero-shot text input"),
        ("Single-tissue extensions", "UCE 0.632; scFoundation 0.635", "single-tissue extension rows"),
    ]
    y0 = 0.86
    ax.text(0.0, 0.98, "C. Evaluation scope", fontsize=12, fontweight="bold", transform=ax.transAxes)
    for i, (surface, models, status) in enumerate(rows):
        y = y0 - (i + 1) * 0.21
        color = PALETTE["classical"] if i == 0 else "#7C3AED" if i == 2 else "#C026D3"
        ax.add_patch(plt.Rectangle((0.0, y - 0.055), 0.025, 0.11, color=color, transform=ax.transAxes))
        ax.text(0.045, y + 0.032, surface, fontsize=10, fontweight="bold", transform=ax.transAxes)
        ax.text(0.045, y - 0.005, models, fontsize=9, color=PALETTE["ink"], transform=ax.transAxes)
        ax.text(0.045, y - 0.044, status, fontsize=8, color=PALETTE["muted"], transform=ax.transAxes)

def build_model_figure() -> list[Path]:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(16, 9), constrained_layout=False)
    grid = fig.add_gridspec(
        nrows=2,
        ncols=2,
        height_ratios=[1.0, 0.9],
        width_ratios=[1.05, 0.95],
        left=0.085,
        right=0.965,
        top=0.82,
        bottom=0.13,
        wspace=0.32,
        hspace=0.36,
    )
    ax_a = fig.add_subplot(grid[:, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 1])
    draw_model_mean_panel(ax_a)
    draw_model_delta_panel(ax_b)
    draw_model_boundary_panel(ax_c)
    fig.suptitle(
        "Larger models do not solve the small spaceflight training-set problem",
        x=0.085,
        y=0.955,
        ha="left",
        fontsize=20,
        fontweight="bold",
    )
    fig.text(
        0.085,
        0.905,
        "Current single-cell foundation models and text LLMs do not reliably beat a tuned classical baseline on these held-out mission tasks.",
        ha="left",
        fontsize=11,
        color=PALETTE["muted"],
    )
    outputs = [
        OUT_DIR / "premium_fig3_model_tier_comparison.png",
        OUT_DIR / "premium_fig3_model_tier_comparison.pdf",
    ]
    fig.savefig(outputs[0], dpi=220)
    fig.savefig(outputs[1])
    plt.close(fig)
    manifest_rows: list[dict[str, object]] = []
    for label, rows in [
        ("model_tier_means", MODEL_TIER_MEANS),
        ("model_tier_deltas", MODEL_TIER_DELTAS),
        ("model_extension_rows", MODEL_EXTENSION_ROWS),
    ]:
        for row in rows:
            item = dict(row)
            item["panel"] = label
            manifest_rows.append(item)
    write_manifest(manifest_rows, "premium_fig3_model_tier_comparison_manifest")
    return outputs


def model_surface_coverage_table() -> list[dict[str, object]]:
    return [
        {
            "surface": "shared_6_task_mean",
            "models": "PCA-LR; scGPT; Mouse-Geneformer",
            "interpretation": "main like-for-like comparison",
            "source": "docs/CANONICAL_RESULTS_V7_1.md:Canonical FM / LLM Snapshot",
        },
        {
            "surface": "text_only_zero_shot_tasks",
            "models": "DeepSeek-V3; Llama-3.3-70B; Gemini-2.5-Flash",
            "interpretation": "separate zero-shot text task surface",
            "source": "evaluation/RESULTS_SUMMARY.md:Tier 3 LLM Zero-Shot Classification",
        },
        {
            "surface": "best_single_tissue_extension_rows",
            "models": "scFoundation; UCE",
            "interpretation": "reported separately; not a shared leaderboard row",
            "source": "docs/CANONICAL_RESULTS_V7_1.md:Canonical FM / LLM Snapshot",
        },
    ]


def draw_model_score_lollipop_panel(ax: plt.Axes) -> None:
    rows = sorted(MODEL_TIER_MEANS, key=lambda row: float(row["score"]), reverse=True)
    y = np.arange(len(rows))
    values = np.array([float(row["score"]) for row in rows])
    colors = [family_color(str(row["family"])) for row in rows]
    ax.axvline(0.5, color=PALETTE["chance"], linestyle="--", linewidth=0.9, zorder=1)
    for yi, value, color in zip(y, values, colors):
        ax.hlines(yi, min(0.5, value), max(0.5, value), color="#CBD5E1", linewidth=2.1, zorder=1)
        ax.scatter(value, yi, s=60, color=color, edgecolor="white", linewidth=0.8, zorder=3)
        if value < 0.5:
            ax.text(value - 0.014, yi, f"{value:.3f}", va="center", ha="right", fontsize=7.7)
        else:
            ax.text(value + 0.014, yi, f"{value:.3f}", va="center", ha="left", fontsize=7.7)
    ax.set_yticks(y)
    ax.set_yticklabels([str(row["model"]) for row in rows])
    ax.set_xlim(0.0, 0.83)
    ax.set_ylim(len(rows) - 0.45, -0.55)
    ax.set_xlabel("Mean AUROC / task score")
    ax.set_title("A. Classical baseline remains highest on the shared task set", loc="left")
    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def draw_model_delta_manuscript_panel(ax: plt.Axes) -> None:
    tissues = ["liver", "gastrocnemius", "kidney", "thymus", "skin", "eye"]
    y = np.arange(len(tissues))
    scgpt = {row["tissue"]: float(row["delta"]) for row in MODEL_TIER_DELTAS if row["model"] == "scGPT"}
    geneformer = {
        row["tissue"]: float(row["delta"])
        for row in MODEL_TIER_DELTAS
        if row["model"] == "Mouse-Geneformer"
    }
    ax.axvline(0, color=PALETTE["ink"], linewidth=0.9)
    for yi, tissue in zip(y, tissues):
        ax.hlines(yi, min(scgpt[tissue], geneformer[tissue], 0), max(scgpt[tissue], geneformer[tissue], 0), color="#CBD5E1", linewidth=2.1, zorder=1)
    ax.scatter([scgpt[t] for t in tissues], y + 0.12, s=54, color="#7C3AED", edgecolor="white", linewidth=0.8, label="scGPT", zorder=3)
    ax.scatter([geneformer[t] for t in tissues], y - 0.12, s=54, color="#A78BFA", edgecolor="white", linewidth=0.8, label="Mouse-Geneformer", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([MANUSCRIPT_TISSUE_LABELS[t] for t in tissues])
    ax.set_xlim(-0.58, 0.08)
    ax.set_xlabel("AUROC change vs matched classical baseline")
    ax.set_title("B. Foundation-model gains are local and usually negative", loc="left")
    ax.legend(loc="upper left", fontsize=7.2, handletextpad=0.3)
    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def build_model_manuscript_figure() -> list[Path]:
    setup_style()
    MANUSCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.2, 5.75), constrained_layout=False)
    grid = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=[1.08, 0.92],
        left=0.20,
        right=0.975,
        top=0.78,
        bottom=0.17,
        hspace=0.72,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[1, 0])
    draw_model_score_lollipop_panel(ax_a)
    draw_model_delta_manuscript_panel(ax_b)
    legend_items = [
        ("Classical ML", family_color("Classical ML")),
        ("Single-cell model", family_color("Single-cell foundation model")),
        ("Text LLM", family_color("Text LLM")),
    ]
    for idx, (label, color) in enumerate(legend_items):
        x = 0.20 + idx * 0.22
        fig.text(x, 0.855, label, fontsize=7.5, color=color, ha="left", va="center", fontweight="bold")
    fig.suptitle(
        "Model scale does not solve held-out mission prediction",
        x=0.20,
        y=0.962,
        ha="left",
        fontsize=13.2,
        fontweight="bold",
    )
    fig.text(
        0.20,
        0.915,
        "Only shared 6-task rows are direct model comparisons; other rows use separate task sets.",
        ha="left",
        fontsize=8.6,
        color=PALETTE["muted"],
    )
    outputs = [
        MANUSCRIPT_DIR / "premium_fig3_model_tier_comparison_manuscript.png",
        MANUSCRIPT_DIR / "premium_fig3_model_tier_comparison_manuscript.pdf",
    ]
    fig.savefig(outputs[0], dpi=320)
    fig.savefig(outputs[1])
    plt.close(fig)
    manifest_rows: list[dict[str, object]] = []
    for label, rows in [
        ("model_tier_means", MODEL_TIER_MEANS),
        ("model_tier_deltas", MODEL_TIER_DELTAS),
    ]:
        for row in rows:
            item = dict(row)
            item["panel"] = label
            item["variant"] = "manuscript"
            manifest_rows.append(item)
    write_manifest(manifest_rows, "premium_fig3_model_tier_comparison_manuscript_manifest")
    table_outputs = write_table(model_surface_coverage_table(), "table3_model_surface_coverage")
    return outputs + table_outputs


def status_color(kind: str) -> str:
    if kind == "ready":
        return PALETTE["ready"]
    if kind == "blocked":
        return PALETTE["blocker"]
    if kind == "diagnostic":
        return PALETTE["diagnostic"]
    if kind == "draft":
        return PALETTE["draft"]
    return PALETTE["muted"]


def draw_v9_status_matrix(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(0.0, 1.02, "A. v9 release lanes are at different readiness levels", fontsize=12, fontweight="bold", transform=ax.transAxes)
    columns = [
        ("lane", "Lane"),
        ("provenance", "Evidence trail"),
        ("evaluation", "Evaluation"),
        ("release", "Release status"),
        ("constraint", "Main constraint"),
    ]
    widths = [0.18, 0.2, 0.22, 0.19, 0.21]
    x_positions = np.cumsum([0] + widths[:-1])
    header_y = 0.9
    row_h = 0.19
    for x, (_, label) in zip(x_positions, columns):
        ax.text(x + 0.012, header_y + 0.04, label, va="center", fontsize=8, fontweight="bold", transform=ax.transAxes)
    ax.hlines(header_y, 0.0, 0.99, color="#D1D5DB", linewidth=1.0, transform=ax.transAxes)
    ax.hlines(header_y + 0.08, 0.0, 0.99, color="#D1D5DB", linewidth=1.0, transform=ax.transAxes)
    for row_i, row in enumerate(V9_STATUS_ROWS):
        y = header_y - (row_i + 1) * row_h
        ax.hlines(y - 0.02, 0.0, 0.99, color="#E5E7EB", linewidth=0.8, transform=ax.transAxes)
        for col_i, (key, _) in enumerate(columns):
            x = x_positions[col_i]
            text = str(row[key])
            color = PALETTE["ink"]
            weight = "bold" if key == "lane" else "normal"
            if key == "release":
                color = status_color(str(row["release_kind"]))
                weight = "bold"
            if key == "constraint":
                color = status_color(str(row["constraint_kind"]))
                weight = "bold" if row["constraint_kind"] == "blocked" else "normal"
            ax.text(
                x + 0.012,
                y + row_h * 0.5 - 0.006,
                text,
                va="center",
                fontsize=8.2,
                color=color,
                fontweight=weight,
                transform=ax.transAxes,
                linespacing=1.15,
            )


def draw_v9_footprint_panel(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(0.0, 1.02, "B. Evidence footprint now spans four lanes", fontsize=12, fontweight="bold", transform=ax.transAxes)
    for idx, row in enumerate(V9_FOOTPRINT_ROWS):
        y = 0.88 - idx * 0.09
        color = status_color(str(row["kind"]))
        ax.hlines(y - 0.04, 0.0, 0.98, color="#E5E7EB", linewidth=0.8, transform=ax.transAxes)
        ax.scatter(0.02, y, s=48, color=color, transform=ax.transAxes)
        ax.text(0.06, y, str(row["group"]), va="center", fontsize=7.6, color=PALETTE["muted"], transform=ax.transAxes)
        ax.text(0.26, y, str(row["display"]), va="center", fontsize=12, fontweight="bold", color=color, transform=ax.transAxes)
        ax.text(0.43, y, str(row["metric"]), va="center", fontsize=7.8, color=PALETTE["ink"], transform=ax.transAxes)


def draw_v9_architecture_panel(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(0.0, 1.02, "A. Data-to-evaluation architecture", fontsize=12, fontweight="bold", transform=ax.transAxes)
    steps = [
        ("Public sources", "OSDR / GEO\nsample tables"),
        ("Task table", "tasks, folds,\nsource inventory"),
        ("Score/check layer", "classification,\nresponse checks"),
        ("Reproducibility layer", "checksums,\ndata hashes,\nsource trace"),
    ]
    xs = np.linspace(0.08, 0.92, len(steps))
    y = 0.6
    ax.plot([xs[0], xs[-1]], [y, y], color="#9CA3AF", linewidth=1.2, transform=ax.transAxes)
    for i, (title, body) in enumerate(steps):
        x = xs[i]
        ax.scatter(x, y, s=115, color="white", edgecolor=PALETTE["classical"], linewidth=1.4, transform=ax.transAxes, zorder=3)
        ax.text(x, y, str(i + 1), ha="center", va="center", fontsize=8, fontweight="bold", color=PALETTE["classical"], transform=ax.transAxes)
        ax.text(x, y + 0.16, title, ha="center", fontsize=8.2, fontweight="bold", transform=ax.transAxes)
        ax.text(x, y - 0.18, body, ha="center", fontsize=7.3, color=PALETTE["muted"], transform=ax.transAxes, linespacing=1.15)
        if i < len(steps) - 1:
            ax.annotate(
                "",
                xy=(xs[i + 1] - 0.04, y),
                xytext=(x + 0.04, y),
                xycoords=ax.transAxes,
                textcoords=ax.transAxes,
                arrowprops={"arrowstyle": "->", "color": PALETTE["muted"], "lw": 1.2},
            )
    lane_y = 0.12
    lane_items = [
        ("indexed metadata", PALETTE["ready"]),
        ("pilot evidence", PALETTE["draft"]),
        ("data pending", PALETTE["blocker"]),
    ]
    ax.text(0.0, lane_y + 0.15, "Evidence depth differs by data layer.", fontsize=8.5, color=PALETTE["ink"], transform=ax.transAxes)
    for i, (label, color) in enumerate(lane_items):
        x = 0.02 + i * 0.28
        ax.scatter(x, lane_y + 0.055, s=80, color=color, transform=ax.transAxes)
        ax.text(x + 0.035, lane_y + 0.055, label, va="center", fontsize=8, transform=ax.transAxes)


def draw_v9_lane_schematic(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(0.0, 1.02, "B. Extension datasets differ in evidence depth", fontsize=12, fontweight="bold", transform=ax.transAxes)
    lanes = [
        {
            "label": "Public mouse bulk",
            "evidence": "22 sources; 33 held-out folds; 24 baselines",
            "status": "metadata indexed",
            "constraint": "local data copy pending",
            "status_kind": "ready",
            "constraint_kind": "blocked",
        },
        {
            "label": "Human organoid",
            "evidence": "2 sources; 42 samples; gene checks",
            "status": "pilot biology checks",
            "constraint": "donor/source coupling",
            "status_kind": "draft",
            "constraint_kind": "draft",
        },
        {
            "label": "Multispecies",
            "evidence": "3 sources; 124 samples; feasibility pilots",
            "status": "feasibility only",
            "constraint": "species namespace and strata",
            "status_kind": "draft",
            "constraint_kind": "draft",
        },
        {
            "label": "Single-cell",
            "evidence": "54 legacy assets indexed",
            "status": "metric plan present",
            "constraint": "missing data files",
            "status_kind": "diagnostic",
            "constraint_kind": "blocked",
        },
    ]
    y_values = np.linspace(0.78, 0.22, len(lanes))
    for y, lane in zip(y_values, lanes):
        status_color_value = status_color(str(lane["status_kind"]))
        constraint_color_value = status_color(str(lane["constraint_kind"]))
        ax.hlines(y, 0.16, 0.88, color="#D1D5DB", linewidth=1.2, transform=ax.transAxes, zorder=1)
        ax.text(0.0, y + 0.035, str(lane["label"]), fontsize=9.2, fontweight="bold", transform=ax.transAxes)
        ax.text(0.0, y - 0.02, str(lane["evidence"]), fontsize=7.7, color=PALETTE["muted"], transform=ax.transAxes)
        ax.scatter(0.32, y, s=92, color=status_color_value, edgecolor="white", linewidth=0.8, transform=ax.transAxes, zorder=3)
        ax.text(0.35, y + 0.018, "evidence", fontsize=6.8, color=PALETTE["muted"], transform=ax.transAxes)
        ax.text(0.35, y - 0.025, str(lane["status"]), fontsize=8.0, color=status_color_value, fontweight="bold", transform=ax.transAxes)
        ax.scatter(0.72, y, s=92, color=constraint_color_value, edgecolor="white", linewidth=0.8, transform=ax.transAxes, zorder=3)
        ax.text(0.75, y + 0.018, "limitation", fontsize=6.8, color=PALETTE["muted"], transform=ax.transAxes)
        ax.text(0.75, y - 0.025, str(lane["constraint"]), fontsize=8.0, color=constraint_color_value, fontweight="bold", transform=ax.transAxes)


def build_v9_status_figure() -> list[Path]:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(16, 9), constrained_layout=False)
    ax_a = fig.add_axes([0.08, 0.45, 0.86, 0.36])
    ax_b = fig.add_axes([0.08, 0.13, 0.86, 0.23])
    draw_v9_architecture_panel(ax_a)
    draw_v9_lane_schematic(ax_b)
    fig.suptitle(
        "SpaceBio-Bench organizes public space-biology evidence layers",
        x=0.08,
        y=0.955,
        ha="left",
        fontsize=20,
        fontweight="bold",
    )
    fig.text(
        0.08,
        0.905,
        "Public-bulk evidence is most complete; organoid, multispecies, and single-cell tracks are smaller extension datasets.",
        ha="left",
        fontsize=11,
        color=PALETTE["muted"],
    )
    outputs = [
        OUT_DIR / "premium_fig4_v9_platform_architecture.png",
        OUT_DIR / "premium_fig4_v9_platform_architecture.pdf",
    ]
    legacy_outputs = [
        OUT_DIR / "premium_fig4_v9_platform_status.png",
        OUT_DIR / "premium_fig4_v9_platform_status.pdf",
    ]
    fig.savefig(outputs[0], dpi=220)
    fig.savefig(outputs[1])
    fig.savefig(legacy_outputs[0], dpi=220)
    fig.savefig(legacy_outputs[1])
    plt.close(fig)
    manifest_rows: list[dict[str, object]] = []
    for row in V9_STATUS_ROWS:
        item = dict(row)
        item["panel"] = "v9_status_matrix"
        manifest_rows.append(item)
    for row in V9_FOOTPRINT_ROWS:
        item = dict(row)
        item["panel"] = "v9_evidence_footprint"
        manifest_rows.append(item)
    manifest_rows.extend(
        [
            {
                "panel": "v9_architecture",
                "step": "public_sources_to_release_gates",
                "source": "docs/V9_LONG_HORIZON_EXECUTION_PLAN.md;v9/README.md",
                "claim_boundary": "platform_status_not_final_benchmark_results",
            }
        ]
    )
    write_manifest(manifest_rows, "premium_fig4_v9_platform_status_manifest")
    table_outputs = []
    table_outputs.extend(write_table(V9_STATUS_ROWS, "table4_v9_lane_readiness"))
    table_outputs.extend(write_table(V9_FOOTPRINT_ROWS, "table4_v9_evidence_footprint"))
    return outputs + legacy_outputs + table_outputs


def bulk_gate_color(status: str) -> str:
    if status == "pass":
        return PALETTE["ready"]
    if status == "blocked":
        return PALETTE["blocker"]
    return PALETTE["draft"]


def draw_bulk_gate_panel(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(0.0, 1.02, "A. Metadata gates are mostly ready; frozen data files are not", fontsize=12, fontweight="bold", transform=ax.transAxes)
    ax.text(0.055, 0.91, "Gate and evidence", fontsize=8, fontweight="bold", transform=ax.transAxes)
    ax.text(0.86, 0.91, "Status", fontsize=8, fontweight="bold", transform=ax.transAxes)
    ax.hlines(0.87, 0.0, 0.98, color="#D1D5DB", linewidth=1.0, transform=ax.transAxes)
    y_start = 0.81
    row_h = 0.105
    for i, row in enumerate(BULK_ALPHA_GATES):
        y = y_start - i * row_h
        color = bulk_gate_color(str(row["status"]))
        ax.hlines(y - 0.05, 0.0, 0.98, color="#E5E7EB", linewidth=0.8, transform=ax.transAxes)
        ax.scatter(0.025, y, s=55, color=color, transform=ax.transAxes)
        ax.text(0.055, y + 0.018, str(row["gate"]), fontsize=8.2, fontweight="bold", va="center", transform=ax.transAxes)
        ax.text(0.055, y - 0.018, str(row["detail"]), fontsize=7.0, color=PALETTE["muted"], va="center", transform=ax.transAxes)
        ax.text(
            0.86,
            y,
            str(row["status"]).upper(),
            fontsize=7.4,
            color=color,
            fontweight="bold",
            va="center",
            transform=ax.transAxes,
        )


def draw_bulk_decision_panel(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(0.0, 1.02, "B. Data-package components", fontsize=12, fontweight="bold", transform=ax.transAxes)
    ax.plot([0.06, 0.06], [0.2, 0.86], color="#D1D5DB", linewidth=1.0, transform=ax.transAxes)
    for i, row in enumerate(BULK_ALPHA_DECISIONS):
        y = 0.78 - i * 0.25
        color = status_color(str(row["kind"]))
        ax.scatter(0.06, y, s=80, color=color, edgecolor="white", linewidth=0.8, transform=ax.transAxes, zorder=3)
        ax.text(0.12, y + 0.035, str(row["path"]), fontsize=9.0, fontweight="bold", transform=ax.transAxes)
        ax.text(0.12, y - 0.015, str(row["detail"]), fontsize=7.8, color=PALETTE["muted"], transform=ax.transAxes)
        ax.text(0.78, y + 0.02, str(row["status"]).upper(), fontsize=7.2, color=color, fontweight="bold", transform=ax.transAxes)


def draw_bulk_language_panel(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(0.0, 1.02, "C. Evidence wording for metadata and data files", fontsize=12, fontweight="bold", transform=ax.transAxes)
    columns = [("Metadata evidence", PALETTE["ready"], 0.0), ("Data-file evidence", PALETTE["diagnostic"], 0.51)]
    for title, color, x in columns:
        ax.text(x + 0.0, 0.89, title, color=color, fontsize=9, fontweight="bold", va="center", transform=ax.transAxes)
        ax.hlines(0.84, x, x + 0.46, color=color, linewidth=1.5, transform=ax.transAxes)
        items = [row for row in BULK_ALPHA_LANGUAGE if row["side"] == title]
        for i, row in enumerate(items):
            y = 0.72 - i * 0.18
            ax.hlines(y - 0.06, x, x + 0.46, color="#E5E7EB", linewidth=0.8, transform=ax.transAxes)
            marker = "+" if title == "Metadata evidence" else "o"
            ax.text(x + 0.025, y, marker, fontsize=12, color=color, fontweight="bold", va="center", transform=ax.transAxes)
            ax.text(x + 0.065, y, str(row["text"]), fontsize=8.1, va="center", transform=ax.transAxes, linespacing=1.15)


def draw_bulk_boundary_schematic(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(0.0, 1.02, "A. Public bulk data-package state", fontsize=12, fontweight="bold", transform=ax.transAxes)
    nodes = [
        {
            "title": "Metadata records",
            "body": "Sources, tasks, folds, and baselines are indexed",
            "color": PALETTE["ready"],
            "x": 0.12,
        },
        {
            "title": "Data-file mirror",
            "body": "Local copy and hash checks remain pending",
            "color": PALETTE["blocker"],
            "x": 0.50,
        },
        {
            "title": "Benchmark package",
            "body": "Metadata and data-file evidence stay separated",
            "color": PALETTE["diagnostic"],
            "x": 0.86,
        },
    ]
    y = 0.58
    ax.hlines(y, nodes[0]["x"], nodes[-1]["x"], color="#D1D5DB", linewidth=1.2, transform=ax.transAxes)
    for i, node in enumerate(nodes):
        x = float(node["x"])
        color = str(node["color"])
        ax.scatter(x, y, s=150, color="white", edgecolor=color, linewidth=1.8, transform=ax.transAxes, zorder=3)
        ax.text(x, y, str(i + 1), ha="center", va="center", fontsize=9, fontweight="bold", color=color, transform=ax.transAxes)
        ax.text(x, y + 0.16, str(node["title"]), ha="center", fontsize=10, fontweight="bold", color=color, transform=ax.transAxes)
        ax.text(x, y - 0.17, str(node["body"]), ha="center", fontsize=8.4, color=PALETTE["muted"], transform=ax.transAxes, linespacing=1.16)
        if i < len(nodes) - 1:
            ax.annotate(
                "",
                xy=(float(nodes[i + 1]["x"]) - 0.07, y),
                xytext=(x + 0.07, y),
                xycoords=ax.transAxes,
                textcoords=ax.transAxes,
                arrowprops={"arrowstyle": "->", "color": PALETTE["muted"], "lw": 1.2},
            )


def build_public_bulk_alpha_figure() -> list[Path]:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(16, 9), constrained_layout=False)
    ax_a = fig.add_axes([0.08, 0.19, 0.86, 0.58])
    draw_bulk_boundary_schematic(ax_a)
    fig.suptitle(
        "Public bulk v9 has indexed metadata but pending data-file mirroring",
        x=0.08,
        y=0.955,
        ha="left",
        fontsize=20,
        fontweight="bold",
    )
    fig.text(
        0.08,
        0.905,
        "Source records, task tables, checksum evidence, and baseline outputs are present; local file-copy verification is still separate.",
        ha="left",
        fontsize=11,
        color=PALETTE["muted"],
    )
    outputs = [
        OUT_DIR / "premium_fig5_public_bulk_release_boundary_schematic.png",
        OUT_DIR / "premium_fig5_public_bulk_release_boundary_schematic.pdf",
    ]
    legacy_outputs = [
        OUT_DIR / "premium_fig5_public_bulk_alpha_boundary.png",
        OUT_DIR / "premium_fig5_public_bulk_alpha_boundary.pdf",
    ]
    fig.savefig(outputs[0], dpi=220)
    fig.savefig(outputs[1])
    fig.savefig(legacy_outputs[0], dpi=220)
    fig.savefig(legacy_outputs[1])
    plt.close(fig)
    manifest_rows: list[dict[str, object]] = []
    for label, rows in [
        ("bulk_alpha_gates", BULK_ALPHA_GATES),
        ("bulk_alpha_decisions", BULK_ALPHA_DECISIONS),
        ("bulk_alpha_language", BULK_ALPHA_LANGUAGE),
    ]:
        for row in rows:
            item = dict(row)
            item["panel"] = label
            manifest_rows.append(item)
    write_manifest(manifest_rows, "premium_fig5_public_bulk_alpha_boundary_manifest")
    table_outputs = []
    table_outputs.extend(write_table(BULK_ALPHA_GATES, "table5_public_bulk_alpha_gates"))
    table_outputs.extend(write_table(BULK_ALPHA_DECISIONS, "table5_public_bulk_release_options"))
    table_outputs.extend(write_table(BULK_ALPHA_LANGUAGE, "table5_public_bulk_claim_language"))
    return outputs + legacy_outputs + table_outputs


def draw_organoid_footprint_panel(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(0.0, 1.02, "A. Organoid dataset is public but small", fontsize=12, fontweight="bold", transform=ax.transAxes)
    visible_rows = [row for row in ORGANOID_FOOTPRINT_ROWS if row["metric"] != "pilot status"]
    xs = np.linspace(0.06, 0.94, len(visible_rows))
    y_mid = 0.58
    ax.plot([xs[0], xs[-1]], [y_mid, y_mid], color="#D1D5DB", linewidth=1.0, transform=ax.transAxes, zorder=1)
    for idx, (x, row) in enumerate(zip(xs, visible_rows)):
        color = status_color(str(row["kind"]))
        ax.scatter(x, y_mid, s=85, color=color, edgecolor="white", linewidth=0.9, transform=ax.transAxes, zorder=3)
        y_num = 0.78 if idx % 2 == 0 else 0.32
        y_lab = 0.68 if idx % 2 == 0 else 0.2
        va = "bottom" if idx % 2 == 0 else "top"
        ax.plot([x, x], [min(y_mid, y_num - 0.03), max(y_mid, y_num + 0.02)], color="#D1D5DB", linewidth=0.8, transform=ax.transAxes)
        ax.text(x, y_num, str(row["display"]), ha="center", va=va, fontsize=16, fontweight="bold", color=color, transform=ax.transAxes)
        ax.text(x, y_lab, str(row["metric"]), ha="center", va="center", fontsize=7.6, fontweight="bold", transform=ax.transAxes)
        ax.text(x, y_lab - 0.07 if idx % 2 == 0 else y_lab - 0.06, str(row["detail"]), ha="center", va="center", fontsize=6.8, color=PALETTE["muted"], transform=ax.transAxes)


def draw_organoid_metric_panel(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(0.0, 1.02, "B. Keep prediction and biology checks separate", fontsize=12, fontweight="bold", transform=ax.transAxes)
    groups = [
        ("Primary task", "prediction", 0.02, PALETTE["ready"]),
        ("Flight-response pattern", "biology check", 0.36, PALETTE["diagnostic"]),
        ("Model gene effects", "secondary check", 0.70, PALETTE["draft"]),
    ]
    width = 0.25
    for family, role, x0, color in groups:
        group_rows = [row for row in ORGANOID_METRIC_ROWS if row["family"] == family]
        ax.text(x0, 0.88, family, fontsize=8.8, fontweight="bold", color=color, transform=ax.transAxes)
        ax.text(x0, 0.82, role, fontsize=7.4, color=PALETTE["muted"], transform=ax.transAxes)
        for j, row in enumerate(group_rows):
            y = 0.62 - j * 0.22
            value = float(row["value"])
            ax.hlines(y, x0, x0 + width, color="#D1D5DB", linewidth=1.0, transform=ax.transAxes)
            ax.scatter(x0 + value * width, y, s=70, color=color, edgecolor="white", linewidth=0.8, transform=ax.transAxes, zorder=3)
            ax.text(x0, y + 0.055, str(row["metric"]), fontsize=7.4, transform=ax.transAxes)
            ax.text(x0 + width + 0.012, y, f"{value:.3f}", va="center", fontsize=7.6, color=color, fontweight="bold", transform=ax.transAxes)
        ax.text(x0, 0.18, "0", fontsize=6.8, color=PALETTE["muted"], transform=ax.transAxes)
        ax.text(x0 + width, 0.18, "1", fontsize=6.8, color=PALETTE["muted"], ha="right", transform=ax.transAxes)


def draw_organoid_topk_panel(ax: plt.Axes) -> None:
    rows = ORGANOID_TOPK_ROWS
    x = np.arange(len(rows))
    enrich = np.array([float(row["enrichment"]) for row in rows])
    colors = [PALETTE["ready"] if float(row["p_value"]) < 0.05 else PALETTE["low"] for row in rows]
    ax.bar(x, enrich, color=colors, width=0.64)
    ax.axhline(1.0, color=PALETTE["chance"], linestyle="--", linewidth=1.0)
    for xi, row in zip(x, rows):
        ax.text(
            xi,
            float(row["enrichment"]) + 0.12,
            f"{row['observed']} vs {row['expected']:.2f}",
            ha="center",
            va="bottom",
            fontsize=7.5,
        )
        p_label = "p<0.05" if float(row["p_value"]) < 0.05 else "not sig."
        ax.text(xi, 0.08, p_label, ha="center", va="bottom", fontsize=7, color=colors[xi], fontweight="bold")
    ax.set_xticks(x)
    ax.set_xticklabels([f"Top {row['k']}" for row in rows])
    ax.set_ylabel("Significant-gene enrichment")
    ax.set_ylim(0, 4.4)
    ax.set_title("C. Top-ranked model genes are enriched, but unevenly", loc="left")
    ax.grid(axis="y", color=PALETTE["grid"], linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def draw_organoid_decision_panel(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.text(0.0, 1.02, "D. Biology-check evidence roles", fontsize=12, fontweight="bold", transform=ax.transAxes)
    header_y = 0.89
    ax.text(0.03, header_y, "Check", fontsize=8, fontweight="bold", transform=ax.transAxes)
    ax.text(0.36, header_y, "Evidence", fontsize=8, fontweight="bold", transform=ax.transAxes)
    ax.text(0.76, header_y, "Role", fontsize=8, fontweight="bold", transform=ax.transAxes)
    ax.hlines(header_y - 0.045, 0.0, 0.98, color="#D1D5DB", linewidth=1.0, transform=ax.transAxes)
    y0 = 0.78
    row_h = 0.115
    for i, row in enumerate(ORGANOID_DECISION_ROWS):
        y = y0 - i * row_h
        color = status_color(str(row["kind"]))
        ax.hlines(y - 0.048, 0.0, 0.98, color="#E5E7EB", linewidth=0.8, transform=ax.transAxes)
        ax.scatter(0.012, y, s=42, color=color, transform=ax.transAxes)
        ax.text(0.03, y, str(row["family"]), fontsize=8.0, fontweight="bold", va="center", transform=ax.transAxes)
        ax.text(0.36, y, str(row["artifact"]), fontsize=7.1, color=PALETTE["muted"], va="center", transform=ax.transAxes)
        ax.text(0.76, y, str(row["decision"]), fontsize=7.4, color=color, fontweight="bold", va="center", transform=ax.transAxes)


def build_organoid_figure() -> list[Path]:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(16, 9), constrained_layout=False)
    grid = fig.add_gridspec(
        nrows=2,
        ncols=2,
        height_ratios=[0.95, 1.05],
        width_ratios=[0.96, 1.04],
        left=0.07,
        right=0.965,
        top=0.82,
        bottom=0.13,
        wspace=0.34,
        hspace=0.38,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, :])
    draw_organoid_footprint_panel(ax_a)
    draw_organoid_metric_panel(ax_b)
    draw_organoid_topk_panel(ax_c)
    fig.suptitle(
        "Human organoids provide a small biology-check dataset",
        x=0.07,
        y=0.955,
        ha="left",
        fontsize=20,
        fontweight="bold",
    )
    fig.text(
        0.07,
        0.905,
        "Two public neural-organoid RNA-seq studies support a 42-sample pilot; study, disease, donor, and microglia factors remain coupled.",
        ha="left",
        fontsize=11,
        color=PALETTE["muted"],
    )
    outputs = [
        OUT_DIR / "premium_fig6_organoid_diagnostic_surface.png",
        OUT_DIR / "premium_fig6_organoid_diagnostic_surface.pdf",
    ]
    fig.savefig(outputs[0], dpi=220)
    fig.savefig(outputs[1])
    plt.close(fig)
    manifest_rows: list[dict[str, object]] = []
    for label, rows in [
        ("organoid_footprint", ORGANOID_FOOTPRINT_ROWS),
        ("organoid_metrics", ORGANOID_METRIC_ROWS),
        ("organoid_topk", ORGANOID_TOPK_ROWS),
    ]:
        for row in rows:
            item = dict(row)
            item["panel"] = label
            manifest_rows.append(item)
    write_manifest(manifest_rows, "premium_fig6_organoid_diagnostic_surface_manifest")
    table_outputs = write_table(ORGANOID_DECISION_ROWS, "table6_organoid_biology_check_decisions")
    return outputs + table_outputs


def draw_organoid_footprint_manuscript_panel(ax: plt.Axes) -> None:
    ax.axis("off")
    rows = [
        ("2", "public\nsources", PALETTE["draft"]),
        ("42", "samples", PALETTE["draft"]),
        ("8", "flight-ground\ngene contrasts", PALETTE["diagnostic"]),
        ("2,368", "significant\ngene rows", PALETTE["diagnostic"]),
    ]
    xs = np.linspace(0.08, 0.92, len(rows))
    ax.hlines(0.50, xs[0], xs[-1], color="#D1D5DB", linewidth=1.0, transform=ax.transAxes)
    for x, (value, label, color) in zip(xs, rows):
        ax.scatter(x, 0.50, s=58, color=color, edgecolor="white", linewidth=0.8, transform=ax.transAxes, zorder=3)
        ax.text(x, 0.72, value, ha="center", va="center", fontsize=13, fontweight="bold", color=color, transform=ax.transAxes)
        ax.text(x, 0.30, label, ha="center", va="center", fontsize=7.0, color=PALETTE["ink"], transform=ax.transAxes, linespacing=1.08)
    ax.text(0.0, 1.02, "A. Small public organoid dataset", fontsize=10.5, fontweight="bold", transform=ax.transAxes)


def draw_organoid_metric_manuscript_panel(ax: plt.Axes) -> None:
    rows = [
        ("Primary AUROC", 0.6147727273, PALETTE["ready"], "o"),
        ("Primary macro-F1", 0.5194508009, PALETTE["ready"], "o"),
        ("Flight direction", 0.7706734868, PALETTE["diagnostic"], "D"),
        ("Flight rank corr.", 0.1760078660, PALETTE["diagnostic"], "D"),
        ("Model gene direction", 0.6078431373, PALETTE["draft"], "s"),
        ("Model gene rank corr.", 0.0867280024, PALETTE["draft"], "s"),
    ]
    y = np.arange(len(rows))
    ax.hlines(y, 0, 1, color="#D1D5DB", linewidth=0.9, zorder=1)
    for yi, (label, value, color, marker) in zip(y, rows):
        ax.scatter(value, yi, s=52, color=color, marker=marker, edgecolor="white", linewidth=0.7, zorder=3)
        ax.text(value + 0.035, yi, f"{value:.3f}", va="center", fontsize=7.3, color=color, fontweight="bold")
    ax.set_yticks(y)
    ax.set_yticklabels([str(row[0]) for row in rows])
    ax.set_xlim(-0.02, 1.05)
    ax.set_ylim(len(rows) - 0.45, -0.55)
    ax.set_xlabel("Value")
    ax.set_title("B. Prediction and biology checks remain separate", loc="left")
    ax.grid(axis="x", color=PALETTE["grid"], linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def draw_organoid_topk_manuscript_panel(ax: plt.Axes) -> None:
    rows = ORGANOID_TOPK_ROWS
    x = np.arange(len(rows))
    enrich = np.array([float(row["enrichment"]) for row in rows])
    colors = [PALETTE["ready"] if float(row["p_value"]) < 0.05 else PALETTE["low"] for row in rows]
    ax.bar(x, enrich, color=colors, width=0.58)
    ax.axhline(1.0, color=PALETTE["chance"], linestyle="--", linewidth=0.9)
    for xi, row in zip(x, rows):
        ax.text(
            xi,
            float(row["enrichment"]) + 0.12,
            f"{row['observed']} / {row['expected']:.2f}",
            ha="center",
            va="bottom",
            fontsize=7.0,
        )
    ax.set_xticks(x)
    ax.set_xticklabels([f"Top {row['k']}" for row in rows])
    ax.set_ylabel("Enrichment")
    ax.set_ylim(0, 4.35)
    ax.set_title("C. Top-ranked model genes overlap flight-response genes", loc="left")
    ax.grid(axis="y", color=PALETTE["grid"], linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def build_organoid_manuscript_figure() -> list[Path]:
    setup_style()
    MANUSCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.2, 5.7), constrained_layout=False)
    grid = fig.add_gridspec(
        nrows=3,
        ncols=1,
        height_ratios=[0.72, 1.18, 1.02],
        left=0.245,
        right=0.975,
        top=0.80,
        bottom=0.12,
        hspace=0.70,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[1, 0])
    ax_c = fig.add_subplot(grid[2, 0])
    draw_organoid_footprint_manuscript_panel(ax_a)
    draw_organoid_metric_manuscript_panel(ax_b)
    draw_organoid_topk_manuscript_panel(ax_c)
    fig.suptitle(
        "Human organoid biology checks",
        x=0.245,
        y=0.965,
        ha="left",
        fontsize=13.2,
        fontweight="bold",
    )
    fig.text(
        0.245,
        0.917,
        "Public neural-organoid RNA-seq adds flight-response gene checks to prediction.",
        ha="left",
        fontsize=8.6,
        color=PALETTE["muted"],
    )
    outputs = [
        MANUSCRIPT_DIR / "premium_fig6_organoid_biology_check_manuscript.png",
        MANUSCRIPT_DIR / "premium_fig6_organoid_biology_check_manuscript.pdf",
    ]
    fig.savefig(outputs[0], dpi=320)
    fig.savefig(outputs[1])
    plt.close(fig)
    manifest_rows: list[dict[str, object]] = []
    for label, rows in [
        ("organoid_footprint_manuscript", ORGANOID_FOOTPRINT_ROWS),
        ("organoid_metrics_manuscript", ORGANOID_METRIC_ROWS),
        ("organoid_topk_manuscript", ORGANOID_TOPK_ROWS),
    ]:
        for row in rows:
            item = dict(row)
            item["panel"] = label
            item["variant"] = "manuscript"
            manifest_rows.append(item)
    write_manifest(manifest_rows, "premium_fig6_organoid_biology_check_manuscript_manifest")
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--figure",
        choices=[
            "core",
            "pathway",
            "pathway_manuscript",
            "model",
            "model_manuscript",
            "manuscript",
            "v9",
            "bulk_alpha",
            "organoid",
            "organoid_manuscript",
            "all",
        ],
        default="all",
        help="Figure family to build.",
    )
    args = parser.parse_args()
    generated: list[Path] = []
    if args.figure in {"core", "all"}:
        generated.extend(build_core_figure())
    if args.figure in {"pathway", "all"}:
        generated.extend(build_pathway_figure())
    if args.figure in {"pathway_manuscript", "manuscript", "all"}:
        generated.extend(build_pathway_manuscript_figure())
    if args.figure in {"model", "all"}:
        generated.extend(build_model_figure())
    if args.figure in {"model_manuscript", "manuscript", "all"}:
        generated.extend(build_model_manuscript_figure())
    if args.figure in {"v9", "all"}:
        generated.extend(build_v9_status_figure())
    if args.figure in {"bulk_alpha", "all"}:
        generated.extend(build_public_bulk_alpha_figure())
    if args.figure in {"organoid", "all"}:
        generated.extend(build_organoid_figure())
    if args.figure in {"organoid_manuscript", "manuscript", "all"}:
        generated.extend(build_organoid_manuscript_figure())
    print(json.dumps({"generated": [str(path.relative_to(ROOT)) for path in generated]}, indent=2))


if __name__ == "__main__":
    main()

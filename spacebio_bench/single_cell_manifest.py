"""Single-cell AnnData task-manifest draft helpers for SpaceBio-Bench v9."""

from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any

from .manifests import validate_task_manifest


DEFAULT_SC_ANNDATA_MANIFEST_DRAFT_ID = "v9_sc_002_rrrm1_blood_anndata_manifest_draft"
DEFAULT_SC_ANNDATA_TASK_ID = "draft_rrrm1_blood_single_cell_spaceflight"

SC_ANNDATA_MANIFEST_DRAFT_SUMMARY_FIELDS = [
    "manifest_draft_id",
    "task_id",
    "decision_status",
    "selected_source_id",
    "selected_tissue",
    "sample_count",
    "flight_sample_count",
    "ground_sample_count",
    "qc_cell_count",
    "qc_gene_count",
    "local_h5ad_count",
    "local_payload_status",
    "v9_manifest_validation_status",
    "metric_profile",
    "next_required_block",
    "claim_boundary",
]

SC_ANNDATA_MANIFEST_BLOCKER_FIELDS = [
    "manifest_draft_id",
    "blocker_id",
    "blocker_status",
    "required_before_runnable",
    "evidence_path",
    "owner_action",
]


def _resolve_path(path: str | Path, repo_root: Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return repo_root / candidate


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _animal_id(source_name: str) -> str:
    # Keep the condition prefix because RRRM-1 names can otherwise collide
    # between matched Ground and Flight samples such as HC_LAR_09/FL_LAR_09.
    return source_name


def _condition_label(condition: str) -> str:
    if condition == "FLT":
        return "Flight"
    if condition == "GC":
        return "Ground"
    return condition


def _candidate_folds(sample_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    folds: list[dict[str, Any]] = []
    for row in sample_rows:
        heldout_animal = _animal_id(row["source_name"])
        heldout_label = _condition_label(row["condition"])
        train_rows = [other for other in sample_rows if other is not row]
        label_counts = Counter(_condition_label(other["condition"]) for other in train_rows)
        folds.append(
            {
                "fold_id": f"holdout_{heldout_animal}",
                "status": "draft_blocked_pending_local_anndata_payload",
                "heldout_animal_id": heldout_animal,
                "heldout_sample_id": row["gsm"],
                "heldout_srx": row["srx"],
                "heldout_label": heldout_label,
                "n_train_samples": len(train_rows),
                "n_test_samples": 1,
                "train_label_distribution": dict(sorted(label_counts.items())),
                "test_label_distribution": {heldout_label: 1},
            }
        )
    return folds


def build_sc_anndata_manifest_draft(
    *,
    repo_root: str | Path = ".",
    manifest_draft_id: str = DEFAULT_SC_ANNDATA_MANIFEST_DRAFT_ID,
    task_id: str = DEFAULT_SC_ANNDATA_TASK_ID,
    selected_source_id: str = "OSD-918",
) -> dict[str, Any]:
    """Build a non-runnable RRRM-1 blood AnnData task-manifest draft."""

    root = Path(repo_root).resolve()
    inventory_summary_path = root / "v9" / "sc_spaceflight" / "asset_inventory_summary.csv"
    inventory_path = root / "v9" / "sc_spaceflight" / "asset_inventory.csv"
    payload_scan_path = root / "v9" / "sc_spaceflight" / "local_payload_scan.csv"
    benchmark_ready_path = root / "v2" / "docs" / "RRRM1_BENCHMARK_READY_MANIFEST_2026-03-12.csv"
    sample_qc_path = root / "v2" / "docs" / "RRRM1_SAMPLE_QC_SUMMARY_2026-03-12.csv"
    condition_map_path = root / "v2" / "docs" / "RRRM1_SRX_CONDITION_MAP.csv"
    annotation_summary_path = root / "v2" / "docs" / "RRRM1_BROAD_ANNOTATION_SUMMARY_2026-03-12.md"

    benchmark_rows = _read_csv_rows(benchmark_ready_path)
    qc_rows = _read_csv_rows(sample_qc_path)
    condition_rows = _read_csv_rows(condition_map_path)
    payload_rows = _read_csv_rows(payload_scan_path)
    selected_benchmark = next(
        row for row in benchmark_rows if row["osd"] == selected_source_id
    )
    selected_qc = next(row for row in qc_rows if row["osd"] == selected_source_id)
    selected_samples = [
        row for row in condition_rows if row["osd"] == selected_source_id
    ]
    label_counts = Counter(_condition_label(row["condition"]) for row in selected_samples)
    local_h5ad_count = sum(1 for row in payload_rows if row["payload_type"] == "h5ad")
    local_payload_status = (
        "local_anndata_payload_present_requires_review"
        if local_h5ad_count
        else "blocked_no_local_anndata_payload"
    )
    claim_boundary = "anndata_manifest_contract_only_no_local_payload_or_score_claim"

    manifest = {
        "schema_version": "0.1.0",
        "task_id": task_id,
        "task_family": "sc_spaceflight",
        "title": "Draft RRRM-1 blood single-cell Flight vs Ground AnnData task",
        "release_status": "draft_non_runnable_contract",
        "runnable_status": local_payload_status,
        "claim_boundary": claim_boundary,
        "organism": "Mus musculus",
        "species_common_name": "mouse",
        "material_type": "single_cells",
        "biospecimen_type": "blood",
        "model_system": "RRRM-1_mouse_spaceflight",
        "assay_modality": "single_cell_rna_seq",
        "platform": "10x_3prime_v3_STARsolo_legacy",
        "feature_namespace": "mouse_gene",
        "ground_control_type": "matched_ground_control",
        "spaceflight_environment": "LEO_or_ISS_spaceflight",
        "source_records": [
            {
                "source_id": selected_source_id,
                "url_or_accession": selected_source_id,
                "source_url": f"https://osdr.nasa.gov/bio/repo/data/studies/{selected_source_id}",
                "access_status": "public_processed_reference_pending_review",
                "privacy_class": "public",
                "checksum_status": "legacy_rrrm_inventory_only_no_local_payload_hash",
                "organism": "Mus musculus",
                "species_common_name": "mouse",
                "taxon_id": "10090",
                "assay_modality": "single_cell_rna_seq",
                "platform": selected_benchmark["assay"],
                "biospecimen_type": selected_benchmark["tissue"],
                "material_type": "single_cells",
                "mission": "RRRM-1",
                "legacy_h5ad_path": selected_benchmark["annotated_h5ad_path"],
                "legacy_h5ad_status": selected_benchmark["status"],
                "benchmark_readiness": selected_benchmark["benchmark_readiness"],
                "local_payload_status": local_payload_status,
                "notes": "Selected as first v9 sc_spaceflight manifest contract because RRRM-1 blood has compact sample/QC/annotation evidence.",
            }
        ],
        "anndata_contract": {
            "input_format": "h5ad",
            "status": "contract_only_payload_not_local",
            "required_matrix": "X or raw counts; exact normalized/count layer must be frozen in V9-SC-003",
            "required_obs_fields": [
                "cell_id",
                "sample_id",
                "source_id",
                "tissue",
                "condition",
                "flight_ground_label",
                "animal_id",
                "age_months",
                "broad_celltype",
            ],
            "recommended_obs_fields": [
                "srx",
                "gsm",
                "cluster_id",
                "n_genes_by_counts",
                "total_counts",
                "pct_counts_mt",
            ],
            "required_var_fields": [
                "gene_symbol",
                "feature_id",
            ],
            "recommended_var_fields": [
                "ensembl_id",
                "highly_variable",
                "mitochondrial",
            ],
            "required_uns_fields": [
                "spacebio_bench_source_id",
                "spacebio_bench_task_id",
                "preprocessing_summary",
                "annotation_source",
            ],
        },
        "split": {
            "name": "rrrm1_blood_leave_animal_out_draft",
            "unit": "animal_id",
            "strategy": "leave_one_animal_out_within_rrrm1_blood",
            "status": "draft_blocked_pending_local_anndata_payload",
            "candidate_folds": _candidate_folds(selected_samples),
            "sample_count": len(selected_samples),
            "label_distribution": dict(sorted(label_counts.items())),
            "blocking_factors": [
                "animal_id",
                "age_months",
                "broad_celltype",
                "library_or_sample_id",
            ],
        },
        "metrics": [
            {
                "metric_id": "state_label_auroc",
                "profile": "genelab_sc",
                "interpretation": "Flight versus Ground ranking within the held-out animal; formula must be finalized in V9-SC-003.",
            },
            {
                "metric_id": "state_label_auprc",
                "profile": "genelab_sc",
                "interpretation": "Class-imbalance-aware Flight versus Ground ranking; non-primary until labels and folds are frozen.",
            },
            {
                "metric_id": "mission_discrimination",
                "profile": "genelab_sc",
                "interpretation": "Representation diagnostic for mission/state separation, not a biological mechanism claim.",
            },
            {
                "metric_id": "de_overlap_at_n",
                "profile": "genelab_sc",
                "interpretation": "Cell-type-specific DE recovery metric; blocked until DE reference and cell-type strata are frozen.",
            },
            {
                "metric_id": "de_direction_match",
                "profile": "genelab_sc",
                "interpretation": "Direction agreement for response genes; diagnostic-only until V9-SC-003 metric specification.",
            },
        ],
        "output": {
            "prediction_format": "csv_or_h5ad_contract",
            "primary_artifacts": [
                "predictions.csv",
                "metrics.json",
                "run_manifest.json",
                "response_signature.csv",
            ],
            "required_columns": [
                "cell_id",
                "sample_id",
                "true_label",
                "predicted_label",
            ],
            "optional_columns": [
                "flight_probability",
                "ground_probability",
                "predicted_broad_celltype",
                "embedding_*",
            ],
            "label_domain": [
                "Ground",
                "Flight",
            ],
            "positive_label": "Flight",
        },
        "provenance": {
            "asset_inventory_summary": inventory_summary_path.relative_to(root).as_posix(),
            "asset_inventory": inventory_path.relative_to(root).as_posix(),
            "local_payload_scan": payload_scan_path.relative_to(root).as_posix(),
            "benchmark_ready_manifest": benchmark_ready_path.relative_to(root).as_posix(),
            "sample_qc_summary": sample_qc_path.relative_to(root).as_posix(),
            "condition_map": condition_map_path.relative_to(root).as_posix(),
            "annotation_summary": annotation_summary_path.relative_to(root).as_posix(),
            "local_payload_status": local_payload_status,
            "legacy_score_status": "reference_only_not_v9_benchmark_result",
            "manifest_status": "draft_contract_only",
        },
        "sample_qc": {
            "n_cells": int(selected_qc["n_cells"]),
            "n_genes": int(selected_qc["n_genes"]),
            "median_total_counts": float(selected_qc["median_total_counts"]),
            "median_n_genes_by_counts": float(selected_qc["median_n_genes_by_counts"]),
            "median_pct_counts_mt": float(selected_qc["median_pct_counts_mt"]),
        },
        "limitations": [
            "Draft manifest contract only; no v9 single-cell benchmark result is claimed.",
            "No local .h5ad, .loom, or .mtx payload is present in this repository scan.",
            "Legacy h5ad paths point to scratch locations and are not distributable package resources.",
            "AnnData obs/var/uns fields have not been verified locally.",
            "Metric formulas and DE-reference policy must be finalized in V9-SC-003.",
            "Legacy RRRM figures, scripts, and evaluation outputs remain reference-only.",
        ],
    }
    validate_task_manifest(manifest)

    blockers = [
        {
            "manifest_draft_id": manifest_draft_id,
            "blocker_id": "local_anndata_payload_absent",
            "blocker_status": "blocking_runnable_task",
            "required_before_runnable": "stage or regenerate the selected RRRM-1 blood AnnData object locally",
            "evidence_path": payload_scan_path.relative_to(root).as_posix(),
            "owner_action": "do not claim runnable sc_spaceflight task until payload exists",
        },
        {
            "manifest_draft_id": manifest_draft_id,
            "blocker_id": "obs_var_contract_unverified",
            "blocker_status": "blocking_evaluator",
            "required_before_runnable": "inspect AnnData obs/var/uns fields against the manifest contract",
            "evidence_path": benchmark_ready_path.relative_to(root).as_posix(),
            "owner_action": "verify required obs/var fields after local payload staging",
        },
        {
            "manifest_draft_id": manifest_draft_id,
            "blocker_id": "metric_formula_spec_pending",
            "blocker_status": "blocking_leaderboard",
            "required_before_runnable": "write V9-SC-003 metric formulas and skip policy",
            "evidence_path": "spacebio_bench/profiles.py",
            "owner_action": "define genelab_sc metrics before evaluator implementation",
        },
        {
            "manifest_draft_id": manifest_draft_id,
            "blocker_id": "legacy_scores_reference_only",
            "blocker_status": "claim_guard",
            "required_before_runnable": "regenerate any scores from a v9 task manifest and run manifest",
            "evidence_path": "v3/evaluation/",
            "owner_action": "do not promote legacy RRRM score surfaces as v9 benchmark results",
        },
    ]
    summary = {
        "manifest_draft_id": manifest_draft_id,
        "task_id": task_id,
        "decision_status": "draft_anndata_manifest_contract_ready_payload_blocked",
        "selected_source_id": selected_source_id,
        "selected_tissue": selected_benchmark["tissue"],
        "sample_count": str(len(selected_samples)),
        "flight_sample_count": str(label_counts["Flight"]),
        "ground_sample_count": str(label_counts["Ground"]),
        "qc_cell_count": selected_qc["n_cells"],
        "qc_gene_count": selected_qc["n_genes"],
        "local_h5ad_count": str(local_h5ad_count),
        "local_payload_status": local_payload_status,
        "v9_manifest_validation_status": "passes_minimal_manifest_validator",
        "metric_profile": "genelab_sc",
        "next_required_block": "V9-SC-003: genelab_sc metric specification",
        "claim_boundary": claim_boundary,
    }
    review_md = f"""# V9-SC-002 AnnData Task Manifest Draft

Status: `{summary["decision_status"]}`

Task id: `{task_id}`

Claim boundary: `{claim_boundary}`

## Draft Decision

Draft one non-runnable `sc_spaceflight` manifest for RRRM-1 blood
(`{selected_source_id}`). This is the most compact first single-cell flagship
contract because the repo already has RRRM-1 sample, QC, annotation, and
benchmark-ready evidence for blood, while local AnnData payload files are still
absent.

## Evidence

- Samples in condition map: {len(selected_samples)}
- Flight samples: {label_counts["Flight"]}
- Ground samples: {label_counts["Ground"]}
- QC cell count: {selected_qc["n_cells"]}
- QC gene count: {selected_qc["n_genes"]}
- Local `.h5ad` payload count: {local_h5ad_count}
- Minimal manifest validator: pass

## Not Claimed

- This is not a runnable v9 single-cell benchmark task.
- It does not claim local AnnData payload availability.
- It does not promote legacy RRRM scripts, figures, or scores as v9 results.
- It does not define final `genelab_sc` formulas.

## Next Block

Run `V9-SC-003: genelab_sc metric specification`. Define metric formulas,
required inputs, skip policy, and which metrics remain diagnostic before any
single-cell evaluator implementation.
"""
    return {
        "summary": [summary],
        "manifest": manifest,
        "blockers": blockers,
        "review_md": review_md,
    }


def write_sc_anndata_manifest_draft(
    *,
    output_dir: str | Path = "v9/sc_spaceflight",
    repo_root: str | Path = ".",
    manifest_draft_id: str = DEFAULT_SC_ANNDATA_MANIFEST_DRAFT_ID,
    task_id: str = DEFAULT_SC_ANNDATA_TASK_ID,
    selected_source_id: str = "OSD-918",
) -> dict[str, Path]:
    """Write the non-runnable RRRM-1 blood AnnData manifest draft."""

    root = Path(repo_root).resolve()
    package = build_sc_anndata_manifest_draft(
        repo_root=root,
        manifest_draft_id=manifest_draft_id,
        task_id=task_id,
        selected_source_id=selected_source_id,
    )
    outdir = _resolve_path(output_dir, root)
    manifest_dir = outdir / "task_manifests"
    outputs = {
        "summary_csv": outdir / "anndata_manifest_draft_summary.csv",
        "summary_json": outdir / "anndata_manifest_draft_summary.json",
        "manifest_json": manifest_dir / f"{task_id}.json",
        "blockers_csv": outdir / "anndata_manifest_blockers.csv",
        "blockers_json": outdir / "anndata_manifest_blockers.json",
        "review_md": outdir / "ANNDATA_MANIFEST_DRAFT_REVIEW.md",
    }
    _write_csv(
        outputs["summary_csv"],
        package["summary"],
        SC_ANNDATA_MANIFEST_DRAFT_SUMMARY_FIELDS,
    )
    _write_json(outputs["summary_json"], package["summary"])
    _write_json(outputs["manifest_json"], package["manifest"])
    _write_csv(
        outputs["blockers_csv"],
        package["blockers"],
        SC_ANNDATA_MANIFEST_BLOCKER_FIELDS,
    )
    _write_json(outputs["blockers_json"], package["blockers"])
    outputs["review_md"].parent.mkdir(parents=True, exist_ok=True)
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs

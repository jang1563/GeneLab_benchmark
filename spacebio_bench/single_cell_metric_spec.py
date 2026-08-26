"""Single-cell metric-specification helpers for SpaceBio-Bench v9."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable

from .profiles import get_metric_profile


DEFAULT_SC_METRIC_SPEC_ID = "v9_sc_003_genelab_sc_metric_spec"

SC_METRIC_SPEC_SUMMARY_FIELDS = [
    "metric_spec_id",
    "decision_status",
    "profile",
    "profile_metric_count",
    "primary_metric_count",
    "diagnostic_metric_count",
    "optional_metric_count",
    "manifest_task_id",
    "manifest_runnable_status",
    "local_payload_status",
    "next_required_block",
    "claim_boundary",
]

SC_METRIC_SPEC_METRIC_FIELDS = [
    "metric_spec_id",
    "metric_id",
    "metric_role",
    "formula",
    "required_inputs",
    "minimum_requirements",
    "aggregation",
    "skip_policy_id",
    "claim_status",
]

SC_METRIC_SPEC_INPUT_FIELDS = [
    "metric_spec_id",
    "input_id",
    "input_status",
    "artifact",
    "required_fields",
    "used_by_metrics",
    "skip_if_missing",
    "notes",
]

SC_METRIC_SPEC_SKIP_FIELDS = [
    "metric_spec_id",
    "metric_id",
    "skip_policy_id",
    "skip_conditions",
    "reported_value",
    "reported_reason",
]


def _resolve_path(path: str | Path, repo_root: Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return repo_root / candidate


def _write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _pipe_join(values: Iterable[str]) -> str:
    return "|".join(str(value) for value in values if str(value))


def _metric_rows(metric_spec_id: str) -> list[dict[str, str]]:
    return [
        {
            "metric_spec_id": metric_spec_id,
            "metric_id": "state_label_auroc",
            "metric_role": "primary_after_payload_freeze",
            "formula": (
                "AUROC(true_flight_label_binary, flight_probability) computed "
                "within each eligible fold or stratum."
            ),
            "required_inputs": "predictions_csv_cell_labels",
            "minimum_requirements": "at least one Flight and one Ground observation with finite flight_probability",
            "aggregation": "macro mean across eligible folds; also report per-fold values",
            "skip_policy_id": "skip_state_label_auc_if_single_class_or_no_probability",
            "claim_status": "formula_spec_only_no_score_claim",
        },
        {
            "metric_spec_id": metric_spec_id,
            "metric_id": "state_label_auprc",
            "metric_role": "primary_after_payload_freeze",
            "formula": (
                "Average precision with Flight as the positive class using "
                "flight_probability."
            ),
            "required_inputs": "predictions_csv_cell_labels",
            "minimum_requirements": "at least one positive Flight observation and finite flight_probability",
            "aggregation": "macro mean across eligible folds; also report per-fold values",
            "skip_policy_id": "skip_state_label_auprc_if_no_positive_or_probability",
            "claim_status": "formula_spec_only_no_score_claim",
        },
        {
            "metric_spec_id": metric_spec_id,
            "metric_id": "mission_discrimination",
            "metric_role": "diagnostic_representation",
            "formula": (
                "Nearest-centroid retrieval accuracy from embedding_* columns: "
                "each test cell is assigned to the closest train centroid by "
                "flight_ground_label."
            ),
            "required_inputs": "embedding_columns",
            "minimum_requirements": "embedding_* columns and at least two labels with train/test cells",
            "aggregation": "macro mean across eligible folds; report label-level confusion where possible",
            "skip_policy_id": "skip_mission_discrimination_if_no_embeddings_or_two_labels",
            "claim_status": "diagnostic_only_no_biological_mechanism_claim",
        },
        {
            "metric_spec_id": metric_spec_id,
            "metric_id": "de_overlap_at_n",
            "metric_role": "diagnostic_de_recovery",
            "formula": (
                "For N in {50,100,200}, compute |topN_predicted_genes "
                "intersect topN_reference_genes| / N within each frozen "
                "cell-type/contrast stratum."
            ),
            "required_inputs": "response_signature_csv|de_reference_table",
            "minimum_requirements": "same contrast and at least N overlapping gene identifiers",
            "aggregation": "report each N and macro mean across eligible strata",
            "skip_policy_id": "skip_de_overlap_if_no_reference_or_insufficient_gene_overlap",
            "claim_status": "diagnostic_only_until_de_reference_frozen",
        },
        {
            "metric_spec_id": metric_spec_id,
            "metric_id": "de_direction_match",
            "metric_role": "diagnostic_de_recovery",
            "formula": (
                "Mean indicator sign(predicted_log2fc_flight_minus_ground) "
                "equals sign(reference_log2fc_flight_minus_ground) over "
                "nonzero overlapping genes."
            ),
            "required_inputs": "response_signature_csv|de_reference_table",
            "minimum_requirements": "at least one overlapping gene with nonzero predicted and reference signs",
            "aggregation": "macro mean across eligible cell-type/contrast strata",
            "skip_policy_id": "skip_de_direction_if_no_signed_overlap",
            "claim_status": "diagnostic_only_until_de_reference_frozen",
        },
        {
            "metric_spec_id": metric_spec_id,
            "metric_id": "expression_mae_when_applicable",
            "metric_role": "optional_reconstruction",
            "formula": (
                "Mean absolute error between predicted and reference expression "
                "values after aligning cells or pseudobulk groups and genes in "
                "the declared normalized layer."
            ),
            "required_inputs": "expression_prediction_layer|anndata_obs_contract|anndata_var_contract",
            "minimum_requirements": "matching cells or groups, matching genes, and declared normalization layer",
            "aggregation": "mean across aligned values; report per-cell-type or per-group when available",
            "skip_policy_id": "skip_expression_mae_unless_reconstruction_artifact_declared",
            "claim_status": "optional_only_not_required_for_classifier_tasks",
        },
    ]


def _input_rows(metric_spec_id: str, manifest: dict[str, Any]) -> list[dict[str, str]]:
    anndata_contract = manifest["anndata_contract"]
    return [
        {
            "metric_spec_id": metric_spec_id,
            "input_id": "anndata_obs_contract",
            "input_status": "required_before_runnable",
            "artifact": "task h5ad obs",
            "required_fields": _pipe_join(anndata_contract["required_obs_fields"]),
            "used_by_metrics": "all cell-level metrics",
            "skip_if_missing": "all metrics skip because cells cannot be aligned to labels and folds",
            "notes": "Must be verified after local AnnData payload staging.",
        },
        {
            "metric_spec_id": metric_spec_id,
            "input_id": "anndata_var_contract",
            "input_status": "required_for_gene_metrics",
            "artifact": "task h5ad var",
            "required_fields": _pipe_join(anndata_contract["required_var_fields"]),
            "used_by_metrics": "de_overlap_at_n|de_direction_match|expression_mae_when_applicable",
            "skip_if_missing": "gene-level metrics skip because gene ids cannot be aligned",
            "notes": "Gene id namespace must be frozen before DE metrics run.",
        },
        {
            "metric_spec_id": metric_spec_id,
            "input_id": "predictions_csv_cell_labels",
            "input_status": "required_for_state_label_metrics",
            "artifact": "predictions.csv",
            "required_fields": "cell_id|sample_id|true_label|predicted_label|flight_probability",
            "used_by_metrics": "state_label_auroc|state_label_auprc",
            "skip_if_missing": "state label ranking metrics skip",
            "notes": "true_label must use the manifest label domain Ground/Flight.",
        },
        {
            "metric_spec_id": metric_spec_id,
            "input_id": "embedding_columns",
            "input_status": "optional_diagnostic",
            "artifact": "predictions.csv",
            "required_fields": "embedding_*",
            "used_by_metrics": "mission_discrimination",
            "skip_if_missing": "mission_discrimination skip with reason no_embedding_columns",
            "notes": "Embeddings are diagnostic and must not be required for classifiers.",
        },
        {
            "metric_spec_id": metric_spec_id,
            "input_id": "response_signature_csv",
            "input_status": "optional_diagnostic_pending_reference",
            "artifact": "response_signature.csv",
            "required_fields": (
                "task_id|source_id|fold_id|cell_type|gene_symbol|"
                "predicted_log2fc_flight_minus_ground|predicted_rank"
            ),
            "used_by_metrics": "de_overlap_at_n|de_direction_match",
            "skip_if_missing": "DE recovery metrics skip",
            "notes": "Only score after the reference contrast and cell-type strata are frozen.",
        },
        {
            "metric_spec_id": metric_spec_id,
            "input_id": "de_reference_table",
            "input_status": "blocked_pending_v9_sc_reference_freeze",
            "artifact": "future v9/sc_spaceflight/de_reference/*.csv",
            "required_fields": (
                "source_id|contrast_id|cell_type|gene_symbol|"
                "reference_log2fc_flight_minus_ground|reference_rank"
            ),
            "used_by_metrics": "de_overlap_at_n|de_direction_match",
            "skip_if_missing": "DE recovery metrics skip with reason de_reference_missing",
            "notes": "Do not infer reference from legacy figures; freeze a table first.",
        },
        {
            "metric_spec_id": metric_spec_id,
            "input_id": "expression_prediction_layer",
            "input_status": "optional_only",
            "artifact": "predicted h5ad layer or aligned expression table",
            "required_fields": "cell_or_group_id|gene_symbol|predicted_expression|reference_expression",
            "used_by_metrics": "expression_mae_when_applicable",
            "skip_if_missing": "expression_mae_when_applicable skip",
            "notes": "Not required for the RRRM-1 blood classifier contract.",
        },
    ]


def _skip_rows(metric_spec_id: str) -> list[dict[str, str]]:
    return [
        {
            "metric_spec_id": metric_spec_id,
            "metric_id": "state_label_auroc",
            "skip_policy_id": "skip_state_label_auc_if_single_class_or_no_probability",
            "skip_conditions": "missing flight_probability|nonfinite scores|only one true label class in evaluated fold",
            "reported_value": "NA",
            "reported_reason": "state_label_auroc_skipped",
        },
        {
            "metric_spec_id": metric_spec_id,
            "metric_id": "state_label_auprc",
            "skip_policy_id": "skip_state_label_auprc_if_no_positive_or_probability",
            "skip_conditions": "missing flight_probability|nonfinite scores|no Flight positives in evaluated fold",
            "reported_value": "NA",
            "reported_reason": "state_label_auprc_skipped",
        },
        {
            "metric_spec_id": metric_spec_id,
            "metric_id": "mission_discrimination",
            "skip_policy_id": "skip_mission_discrimination_if_no_embeddings_or_two_labels",
            "skip_conditions": "no embedding_* columns|fewer than two labels|no train centroid for a label",
            "reported_value": "NA",
            "reported_reason": "mission_discrimination_skipped",
        },
        {
            "metric_spec_id": metric_spec_id,
            "metric_id": "de_overlap_at_n",
            "skip_policy_id": "skip_de_overlap_if_no_reference_or_insufficient_gene_overlap",
            "skip_conditions": "missing response_signature.csv|missing DE reference|fewer than N overlapping genes",
            "reported_value": "NA",
            "reported_reason": "de_overlap_at_n_skipped",
        },
        {
            "metric_spec_id": metric_spec_id,
            "metric_id": "de_direction_match",
            "skip_policy_id": "skip_de_direction_if_no_signed_overlap",
            "skip_conditions": "missing response_signature.csv|missing DE reference|no nonzero signed overlapping genes",
            "reported_value": "NA",
            "reported_reason": "de_direction_match_skipped",
        },
        {
            "metric_spec_id": metric_spec_id,
            "metric_id": "expression_mae_when_applicable",
            "skip_policy_id": "skip_expression_mae_unless_reconstruction_artifact_declared",
            "skip_conditions": "no expression prediction artifact|no declared normalization layer|no aligned cells or genes",
            "reported_value": "NA",
            "reported_reason": "expression_mae_when_applicable_skipped",
        },
    ]


def build_sc_metric_spec(
    *,
    repo_root: str | Path = ".",
    metric_spec_id: str = DEFAULT_SC_METRIC_SPEC_ID,
    manifest_path: str | Path = (
        "v9/sc_spaceflight/task_manifests/"
        "draft_rrrm1_blood_single_cell_spaceflight.json"
    ),
) -> dict[str, Any]:
    """Build a machine-readable `genelab_sc` metric specification."""

    root = Path(repo_root).resolve()
    resolved_manifest_path = _resolve_path(manifest_path, root)
    manifest = json.loads(resolved_manifest_path.read_text())
    profile = get_metric_profile("genelab_sc")
    metrics = _metric_rows(metric_spec_id)
    metric_ids = [row["metric_id"] for row in metrics]
    profile_metric_ids = profile["metrics"]
    if set(metric_ids) != set(profile_metric_ids):
        raise ValueError("metric specification does not match genelab_sc profile")

    inputs = _input_rows(metric_spec_id, manifest)
    skip_rows = _skip_rows(metric_spec_id)
    primary_count = sum(
        1 for row in metrics if row["metric_role"] == "primary_after_payload_freeze"
    )
    optional_count = sum(1 for row in metrics if row["metric_role"] == "optional_reconstruction")
    diagnostic_count = len(metrics) - primary_count - optional_count
    claim_boundary = "genelab_sc_metric_spec_only_no_evaluator_or_score_claim"
    summary = {
        "metric_spec_id": metric_spec_id,
        "decision_status": "genelab_sc_metric_spec_ready_no_evaluator",
        "profile": "genelab_sc",
        "profile_metric_count": str(len(profile_metric_ids)),
        "primary_metric_count": str(primary_count),
        "diagnostic_metric_count": str(diagnostic_count),
        "optional_metric_count": str(optional_count),
        "manifest_task_id": manifest["task_id"],
        "manifest_runnable_status": manifest["runnable_status"],
        "local_payload_status": manifest["provenance"]["local_payload_status"],
        "next_required_block": "V9-SC-004: AnnData payload staging and obs/var audit plan",
        "claim_boundary": claim_boundary,
    }
    review_md = f"""# V9-SC-003 `genelab_sc` Metric Specification

Status: `{summary["decision_status"]}`

Profile: `genelab_sc`

Claim boundary: `{claim_boundary}`

## Scope

This is a metric specification only. It defines formulas, required inputs, and
skip policy for the first non-runnable RRRM-1 blood `sc_spaceflight` manifest:
`{manifest["task_id"]}`.

No evaluator, local AnnData payload, benchmark score, or legacy RRRM score
promotion is claimed by this block.

## Metric Roles

- Primary after payload freeze: `state_label_auroc`, `state_label_auprc`
- Diagnostic representation: `mission_discrimination`
- Diagnostic DE recovery: `de_overlap_at_n`, `de_direction_match`
- Optional reconstruction: `expression_mae_when_applicable`

## Required Boundaries

- State-label metrics require `predictions.csv` with cell labels and
  `flight_probability`.
- Embedding diagnostics require `embedding_*` columns and are skippable.
- DE recovery metrics require both `response_signature.csv` and a frozen
  v9 single-cell DE reference table.
- Expression MAE is optional and must skip unless a reconstruction artifact and
  normalization layer are declared.
- Every skipped metric must report `NA` plus a machine-readable skip reason.

## Next Block

Run `V9-SC-004: AnnData payload staging and obs/var audit plan`. The metric
spec is now fixed enough to decide exactly which local AnnData object and
metadata fields must be staged or regenerated before evaluator implementation.
"""
    return {
        "summary": [summary],
        "metrics": metrics,
        "inputs": inputs,
        "skip_policy": skip_rows,
        "review_md": review_md,
    }


def write_sc_metric_spec(
    *,
    output_dir: str | Path = "v9/sc_spaceflight",
    repo_root: str | Path = ".",
    metric_spec_id: str = DEFAULT_SC_METRIC_SPEC_ID,
    manifest_path: str | Path = (
        "v9/sc_spaceflight/task_manifests/"
        "draft_rrrm1_blood_single_cell_spaceflight.json"
    ),
) -> dict[str, Path]:
    """Write `genelab_sc` metric-specification artifacts."""

    root = Path(repo_root).resolve()
    package = build_sc_metric_spec(
        repo_root=root,
        metric_spec_id=metric_spec_id,
        manifest_path=manifest_path,
    )
    outdir = _resolve_path(output_dir, root)
    outputs = {
        "summary_csv": outdir / "metric_spec_summary.csv",
        "summary_json": outdir / "metric_spec_summary.json",
        "metrics_csv": outdir / "metric_spec_metrics.csv",
        "metrics_json": outdir / "metric_spec_metrics.json",
        "inputs_csv": outdir / "metric_spec_required_inputs.csv",
        "inputs_json": outdir / "metric_spec_required_inputs.json",
        "skip_csv": outdir / "metric_spec_skip_policy.csv",
        "skip_json": outdir / "metric_spec_skip_policy.json",
        "review_md": root / "docs" / "V9_SC_METRIC_SPEC.md",
    }
    _write_csv(
        outputs["summary_csv"],
        package["summary"],
        SC_METRIC_SPEC_SUMMARY_FIELDS,
    )
    _write_json(outputs["summary_json"], package["summary"])
    _write_csv(outputs["metrics_csv"], package["metrics"], SC_METRIC_SPEC_METRIC_FIELDS)
    _write_json(outputs["metrics_json"], package["metrics"])
    _write_csv(outputs["inputs_csv"], package["inputs"], SC_METRIC_SPEC_INPUT_FIELDS)
    _write_json(outputs["inputs_json"], package["inputs"])
    _write_csv(outputs["skip_csv"], package["skip_policy"], SC_METRIC_SPEC_SKIP_FIELDS)
    _write_json(outputs["skip_json"], package["skip_policy"])
    outputs["review_md"].parent.mkdir(parents=True, exist_ok=True)
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs

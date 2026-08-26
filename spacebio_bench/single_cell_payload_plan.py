"""Single-cell AnnData payload-staging plan helpers for SpaceBio-Bench v9."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


DEFAULT_SC_PAYLOAD_STAGING_PLAN_ID = "v9_sc_004_rrrm1_blood_payload_staging_plan"

SC_PAYLOAD_STAGING_PLAN_SUMMARY_FIELDS = [
    "payload_plan_id",
    "decision_status",
    "task_id",
    "selected_source_id",
    "selected_tissue",
    "canonical_payload_path",
    "local_payload_present",
    "required_obs_field_count",
    "required_var_field_count",
    "required_uns_field_count",
    "audit_requirement_count",
    "staging_action_count",
    "next_required_block",
    "claim_boundary",
]

SC_PAYLOAD_STAGING_CANDIDATE_FIELDS = [
    "payload_plan_id",
    "candidate_id",
    "candidate_role",
    "candidate_path",
    "path_status",
    "promotion_status",
    "required_evidence",
    "notes",
]

SC_OBS_VAR_AUDIT_REQUIREMENT_FIELDS = [
    "payload_plan_id",
    "field_id",
    "axis",
    "requirement_status",
    "accepted_aliases",
    "validation_rule",
    "required_for",
    "blocker_if_missing",
]

SC_PAYLOAD_STAGING_ACTION_FIELDS = [
    "payload_plan_id",
    "action_id",
    "action_order",
    "action_type",
    "inputs",
    "outputs",
    "gate_status",
    "success_criteria",
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


def _path_status(path_text: str, repo_root: Path) -> str:
    if "$" in path_text:
        return "unresolved_external_path"
    path = _resolve_path(path_text, repo_root)
    if path.exists():
        return "local_exists"
    if path_text.startswith("v9/"):
        return "planned_repo_path_absent"
    return "local_absent"


def _pipe_join(values: list[str]) -> str:
    return "|".join(values)


def _candidate_rows(
    *,
    payload_plan_id: str,
    manifest: dict[str, Any],
    canonical_payload_path: str,
    repo_root: Path,
) -> list[dict[str, str]]:
    source = manifest["source_records"][0]
    source_id = source["source_id"]
    tissue = source["biospecimen_type"]
    legacy_annotated = source["legacy_h5ad_path"]
    legacy_labeled = (
        f"$SCRATCH_DIR/rrrm1_scrna/{source_id}/"
        f"{source_id}_{tissue}_labeled.h5ad"
    )
    return [
        {
            "payload_plan_id": payload_plan_id,
            "candidate_id": "canonical_v9_payload_target",
            "candidate_role": "target_payload_path",
            "candidate_path": canonical_payload_path,
            "path_status": _path_status(canonical_payload_path, repo_root),
            "promotion_status": "target_blocked_until_payload_staged_and_audited",
            "required_evidence": "local h5ad exists|sha256 recorded|obs_var_audit_passed",
            "notes": "Canonical location for the first runnable RRRM-1 blood v9 payload.",
        },
        {
            "payload_plan_id": payload_plan_id,
            "candidate_id": "legacy_annotated_h5ad",
            "candidate_role": "preferred_source_if_available",
            "candidate_path": legacy_annotated,
            "path_status": _path_status(legacy_annotated, repo_root),
            "promotion_status": "usable_only_after_local_copy_hash_and_contract_audit",
            "required_evidence": "copy provenance|sha256|broad_celltype present|cell barcode stability",
            "notes": "Preferred because the legacy annotation summary names broad cell-type labels.",
        },
        {
            "payload_plan_id": payload_plan_id,
            "candidate_id": "legacy_labeled_h5ad",
            "candidate_role": "metadata_source_if_annotation_missing",
            "candidate_path": legacy_labeled,
            "path_status": _path_status(legacy_labeled, repo_root),
            "promotion_status": "insufficient_without_broad_celltype_annotation",
            "required_evidence": "condition labels|animal_id|age_months|annotation join plan",
            "notes": "Useful fallback for sample metadata but not enough for cell-type diagnostics.",
        },
        {
            "payload_plan_id": payload_plan_id,
            "candidate_id": "starsolo_per_srx_regeneration",
            "candidate_role": "regeneration_source",
            "candidate_path": "$SCRATCH_DIR/rrrm1_scrna/OSD-918/starsolo_per_srx/",
            "path_status": "unresolved_external_path",
            "promotion_status": "valid_regeneration_route_if_raw_or_starsolo_outputs_available",
            "required_evidence": "per-SRX matrix files|condition map|script version|output hash",
            "notes": "Regenerate with the legacy per-SRX merge script, then normalize to the v9 contract.",
        },
    ]


def _obs_var_requirement_rows(
    *,
    payload_plan_id: str,
    manifest: dict[str, Any],
) -> list[dict[str, str]]:
    contract = manifest["anndata_contract"]
    required_obs = contract["required_obs_fields"]
    recommended_obs = contract["recommended_obs_fields"]
    required_var = contract["required_var_fields"]
    recommended_var = contract["recommended_var_fields"]
    required_uns = contract["required_uns_fields"]
    rows: list[dict[str, str]] = []

    obs_aliases = {
        "cell_id": ["obs_names", "barcode", "cell_barcode"],
        "sample_id": ["srx", "gsm", "source_name"],
        "source_id": ["osd", "study_accession"],
        "condition": ["FLT/GC", "flight_ground_label"],
        "flight_ground_label": ["condition", "exp"],
        "broad_celltype": ["predicted.id", "celltype", "cell_type"],
    }
    obs_rules = {
        "cell_id": "unique non-empty cell identifiers; if absent derive from obs_names",
        "sample_id": "maps each cell to exactly one SRX/GSM sample",
        "source_id": "constant OSD-918 for this first payload",
        "tissue": "constant blood and matches manifest biospecimen_type",
        "condition": "domain must be FLT or GC before label normalization",
        "flight_ground_label": "domain must be Flight or Ground after normalization",
        "animal_id": "exactly eight animal/source_name groups for RRRM-1 blood",
        "age_months": "numeric values expected from the condition map",
        "broad_celltype": "non-empty broad label per cell for cell-type diagnostics",
    }
    for field in required_obs:
        rows.append(
            {
                "payload_plan_id": payload_plan_id,
                "field_id": field,
                "axis": "obs",
                "requirement_status": "required_before_runnable",
                "accepted_aliases": _pipe_join(obs_aliases.get(field, [])),
                "validation_rule": obs_rules.get(field, "non-empty field"),
                "required_for": "folds|state_label_metrics|cell_alignment",
                "blocker_if_missing": "payload_not_runnable",
            }
        )
    for field in recommended_obs:
        rows.append(
            {
                "payload_plan_id": payload_plan_id,
                "field_id": field,
                "axis": "obs",
                "requirement_status": "recommended_audit_field",
                "accepted_aliases": "",
                "validation_rule": "record presence, completeness, and value domain",
                "required_for": "diagnostics|provenance",
                "blocker_if_missing": "warn_not_block_unless_needed_by_metric",
            }
        )
    for field in required_var:
        rows.append(
            {
                "payload_plan_id": payload_plan_id,
                "field_id": field,
                "axis": "var",
                "requirement_status": "required_before_gene_metrics",
                "accepted_aliases": "var_names|gene_ids|gene_symbols",
                "validation_rule": "non-empty gene identifiers with stable namespace and duplicate report",
                "required_for": "DE recovery|expression alignment",
                "blocker_if_missing": "gene_metrics_not_runnable",
            }
        )
    for field in recommended_var:
        rows.append(
            {
                "payload_plan_id": payload_plan_id,
                "field_id": field,
                "axis": "var",
                "requirement_status": "recommended_audit_field",
                "accepted_aliases": "",
                "validation_rule": "record presence and completeness",
                "required_for": "diagnostics|feature namespace review",
                "blocker_if_missing": "warn_not_block",
            }
        )
    for field in required_uns:
        rows.append(
            {
                "payload_plan_id": payload_plan_id,
                "field_id": field,
                "axis": "uns",
                "requirement_status": "required_before_release",
                "accepted_aliases": "",
                "validation_rule": "must record v9 task/source/preprocessing/annotation provenance",
                "required_for": "provenance|release_boundary",
                "blocker_if_missing": "payload_not_releasable",
            }
        )
    rows.extend(
        [
            {
                "payload_plan_id": payload_plan_id,
                "field_id": "X",
                "axis": "matrix",
                "requirement_status": "required_before_runnable",
                "accepted_aliases": "raw.X|layers['counts']",
                "validation_rule": "matrix shape matches n_obs x n_vars and finite numeric values",
                "required_for": "model input|feature alignment",
                "blocker_if_missing": "payload_not_runnable",
            },
            {
                "payload_plan_id": payload_plan_id,
                "field_id": "layers.counts",
                "axis": "layers",
                "requirement_status": "preferred_before_release",
                "accepted_aliases": "raw counts|GeneFull filtered counts",
                "validation_rule": "record whether integer count layer exists and which layer evaluators use",
                "required_for": "preprocessing reproducibility",
                "blocker_if_missing": "release_requires_explicit_matrix_policy",
            },
            {
                "payload_plan_id": payload_plan_id,
                "field_id": "raw",
                "axis": "raw",
                "requirement_status": "optional_audit_field",
                "accepted_aliases": "",
                "validation_rule": "record presence, shape, and feature namespace if available",
                "required_for": "diagnostics",
                "blocker_if_missing": "no_blocker",
            },
        ]
    )
    return rows


def _action_rows(
    *,
    payload_plan_id: str,
    canonical_payload_path: str,
) -> list[dict[str, str]]:
    return [
        {
            "payload_plan_id": payload_plan_id,
            "action_id": "locate_or_regenerate_rrrm1_blood_h5ad",
            "action_order": "1",
            "action_type": "payload_staging",
            "inputs": "legacy annotated h5ad or STARsolo per-SRX matrices",
            "outputs": canonical_payload_path,
            "gate_status": "blocked_until_payload_available",
            "success_criteria": "local h5ad exists at canonical path",
        },
        {
            "payload_plan_id": payload_plan_id,
            "action_id": "normalize_obs_contract",
            "action_order": "2",
            "action_type": "metadata_normalization",
            "inputs": "AnnData.obs|RRRM1_SRX_CONDITION_MAP.csv",
            "outputs": "required obs fields with Flight/Ground labels",
            "gate_status": "pending_payload",
            "success_criteria": "all required obs fields present and complete",
        },
        {
            "payload_plan_id": payload_plan_id,
            "action_id": "normalize_var_contract",
            "action_order": "3",
            "action_type": "feature_namespace_normalization",
            "inputs": "AnnData.var|var_names",
            "outputs": "gene_symbol and feature_id columns",
            "gate_status": "pending_payload",
            "success_criteria": "gene namespace and duplicate policy are recorded",
        },
        {
            "payload_plan_id": payload_plan_id,
            "action_id": "write_uns_provenance",
            "action_order": "4",
            "action_type": "provenance_embedding",
            "inputs": "task manifest|metric spec|preprocessing notes",
            "outputs": "spacebio_bench_* uns keys",
            "gate_status": "pending_payload",
            "success_criteria": "uns provenance keys match manifest contract",
        },
        {
            "payload_plan_id": payload_plan_id,
            "action_id": "compute_payload_hash_and_shape_manifest",
            "action_order": "5",
            "action_type": "integrity_audit",
            "inputs": canonical_payload_path,
            "outputs": "future payload_manifest.csv",
            "gate_status": "pending_payload",
            "success_criteria": "sha256, n_obs, n_vars, and file size recorded",
        },
        {
            "payload_plan_id": payload_plan_id,
            "action_id": "run_obs_var_audit",
            "action_order": "6",
            "action_type": "contract_audit",
            "inputs": canonical_payload_path,
            "outputs": "future obs_var_audit.csv",
            "gate_status": "pending_payload",
            "success_criteria": "required fields pass or emit machine-readable blockers",
        },
        {
            "payload_plan_id": payload_plan_id,
            "action_id": "promote_manifest_to_runnable_candidate",
            "action_order": "7",
            "action_type": "manifest_gate",
            "inputs": "payload manifest|obs_var_audit|metric spec",
            "outputs": "updated runnable-status decision",
            "gate_status": "blocked_until_actions_1_to_6_pass",
            "success_criteria": "no score claim; only runnable candidate status changes",
        },
    ]


def build_sc_payload_staging_plan(
    *,
    repo_root: str | Path = ".",
    payload_plan_id: str = DEFAULT_SC_PAYLOAD_STAGING_PLAN_ID,
    manifest_path: str | Path = (
        "v9/sc_spaceflight/task_manifests/"
        "draft_rrrm1_blood_single_cell_spaceflight.json"
    ),
    metric_spec_summary_path: str | Path = "v9/sc_spaceflight/metric_spec_summary.csv",
    canonical_payload_path: str = (
        "v9/sc_spaceflight/payloads/rrrm1_blood/"
        "OSD-918_blood_rrrm1_bench.h5ad"
    ),
) -> dict[str, Any]:
    """Build a machine-readable staging and obs/var audit plan for RRRM-1 blood."""

    root = Path(repo_root).resolve()
    resolved_manifest_path = _resolve_path(manifest_path, root)
    resolved_metric_summary = _resolve_path(metric_spec_summary_path, root)
    manifest = json.loads(resolved_manifest_path.read_text())
    candidate_rows = _candidate_rows(
        payload_plan_id=payload_plan_id,
        manifest=manifest,
        canonical_payload_path=canonical_payload_path,
        repo_root=root,
    )
    audit_rows = _obs_var_requirement_rows(
        payload_plan_id=payload_plan_id,
        manifest=manifest,
    )
    action_rows = _action_rows(
        payload_plan_id=payload_plan_id,
        canonical_payload_path=canonical_payload_path,
    )
    contract = manifest["anndata_contract"]
    source = manifest["source_records"][0]
    local_payload_present = (
        candidate_rows[0]["path_status"] == "local_exists"
        and manifest["runnable_status"] != "blocked_no_local_anndata_payload"
    )
    claim_boundary = "payload_staging_plan_only_no_local_payload_or_score_claim"
    summary = {
        "payload_plan_id": payload_plan_id,
        "decision_status": "payload_staging_plan_ready_no_local_payload",
        "task_id": manifest["task_id"],
        "selected_source_id": source["source_id"],
        "selected_tissue": source["biospecimen_type"],
        "canonical_payload_path": canonical_payload_path,
        "local_payload_present": str(local_payload_present).lower(),
        "required_obs_field_count": str(len(contract["required_obs_fields"])),
        "required_var_field_count": str(len(contract["required_var_fields"])),
        "required_uns_field_count": str(len(contract["required_uns_fields"])),
        "audit_requirement_count": str(len(audit_rows)),
        "staging_action_count": str(len(action_rows)),
        "next_required_block": "V9-SC-005: AnnData obs/var audit implementation",
        "claim_boundary": claim_boundary,
    }
    review_md = f"""# V9-SC-004 AnnData Payload Staging And Obs/Var Audit Plan

Status: `{summary["decision_status"]}`

Task id: `{manifest["task_id"]}`

Claim boundary: `{claim_boundary}`

## Decision

The first `sc_spaceflight` payload target is `{canonical_payload_path}` for
RRRM-1 blood (`{source["source_id"]}`). This block does not stage, copy,
download, or regenerate the payload. It fixes the gate that must be satisfied
before the manifest can become runnable.

## Required Gate

- Required obs fields: {summary["required_obs_field_count"]}
- Required var fields: {summary["required_var_field_count"]}
- Required uns fields: {summary["required_uns_field_count"]}
- Audit requirement rows: {summary["audit_requirement_count"]}
- Staging actions: {summary["staging_action_count"]}
- Local payload present now: `{summary["local_payload_present"]}`

## Source Preference

Use the annotated legacy RRRM-1 blood h5ad first if it can be staged locally,
because it is the candidate with broad cell-type labels. If it is unavailable,
regenerate from per-SRX STARsolo matrices and then reapply the v9 obs/var/uns
contract. A labeled but unannotated h5ad can help recover sample metadata but
is not enough for the `genelab_sc` diagnostic surface.

## Not Claimed

- No local h5ad payload is claimed.
- No payload checksum is claimed.
- No obs/var audit pass is claimed.
- No evaluator, leaderboard result, or legacy RRRM score promotion is claimed.

## Inputs Used

- `{resolved_manifest_path.relative_to(root).as_posix()}`
- `{resolved_metric_summary.relative_to(root).as_posix()}`
- `v2/docs/RRRM1_BENCHMARK_READY_MANIFEST_2026-03-12.csv`
- `v2/docs/RRRM1_BROAD_ANNOTATION_SUMMARY_2026-03-12.md`

## Next Block

Run `V9-SC-005: AnnData obs/var audit implementation`. That block should add
the actual h5ad reader/auditor, but still skip cleanly until a local payload is
available.
"""
    return {
        "summary": [summary],
        "candidates": candidate_rows,
        "audit_requirements": audit_rows,
        "staging_actions": action_rows,
        "review_md": review_md,
    }


def write_sc_payload_staging_plan(
    *,
    output_dir: str | Path = "v9/sc_spaceflight",
    repo_root: str | Path = ".",
    payload_plan_id: str = DEFAULT_SC_PAYLOAD_STAGING_PLAN_ID,
    manifest_path: str | Path = (
        "v9/sc_spaceflight/task_manifests/"
        "draft_rrrm1_blood_single_cell_spaceflight.json"
    ),
    metric_spec_summary_path: str | Path = "v9/sc_spaceflight/metric_spec_summary.csv",
    canonical_payload_path: str = (
        "v9/sc_spaceflight/payloads/rrrm1_blood/"
        "OSD-918_blood_rrrm1_bench.h5ad"
    ),
) -> dict[str, Path]:
    """Write the RRRM-1 blood payload staging and obs/var audit plan."""

    root = Path(repo_root).resolve()
    package = build_sc_payload_staging_plan(
        repo_root=root,
        payload_plan_id=payload_plan_id,
        manifest_path=manifest_path,
        metric_spec_summary_path=metric_spec_summary_path,
        canonical_payload_path=canonical_payload_path,
    )
    outdir = _resolve_path(output_dir, root)
    outputs = {
        "summary_csv": outdir / "payload_staging_plan_summary.csv",
        "summary_json": outdir / "payload_staging_plan_summary.json",
        "candidates_csv": outdir / "payload_staging_candidates.csv",
        "candidates_json": outdir / "payload_staging_candidates.json",
        "audit_csv": outdir / "obs_var_audit_requirements.csv",
        "audit_json": outdir / "obs_var_audit_requirements.json",
        "actions_csv": outdir / "payload_staging_actions.csv",
        "actions_json": outdir / "payload_staging_actions.json",
        "review_md": root / "docs" / "V9_SC_PAYLOAD_STAGING_PLAN.md",
    }
    _write_csv(
        outputs["summary_csv"],
        package["summary"],
        SC_PAYLOAD_STAGING_PLAN_SUMMARY_FIELDS,
    )
    _write_json(outputs["summary_json"], package["summary"])
    _write_csv(
        outputs["candidates_csv"],
        package["candidates"],
        SC_PAYLOAD_STAGING_CANDIDATE_FIELDS,
    )
    _write_json(outputs["candidates_json"], package["candidates"])
    _write_csv(
        outputs["audit_csv"],
        package["audit_requirements"],
        SC_OBS_VAR_AUDIT_REQUIREMENT_FIELDS,
    )
    _write_json(outputs["audit_json"], package["audit_requirements"])
    _write_csv(
        outputs["actions_csv"],
        package["staging_actions"],
        SC_PAYLOAD_STAGING_ACTION_FIELDS,
    )
    _write_json(outputs["actions_json"], package["staging_actions"])
    outputs["review_md"].parent.mkdir(parents=True, exist_ok=True)
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs

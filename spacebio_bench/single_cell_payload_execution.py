"""Single-cell payload staging execution helpers for SpaceBio-Bench v9."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


DEFAULT_SC_PAYLOAD_STAGING_EXECUTION_ID = (
    "v9_sc_006_rrrm1_blood_payload_staging_execution"
)

SC_PAYLOAD_STAGING_EXECUTION_SUMMARY_FIELDS = [
    "execution_id",
    "decision_status",
    "task_id",
    "canonical_payload_path",
    "canonical_payload_status",
    "selected_route",
    "payload_manifest_status",
    "obs_var_audit_status",
    "local_payload_staged",
    "candidate_count",
    "regeneration_step_count",
    "next_required_block",
    "claim_boundary",
]

SC_PAYLOAD_STAGING_EXECUTION_CANDIDATE_FIELDS = [
    "execution_id",
    "candidate_id",
    "candidate_path",
    "check_scope",
    "path_status",
    "evidence_status",
    "action_decision",
    "notes",
]

SC_PAYLOAD_REGENERATION_STEP_FIELDS = [
    "execution_id",
    "step_id",
    "step_order",
    "required_input",
    "input_status",
    "output",
    "gate_status",
    "notes",
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


def _path_status(path_text: str, repo_root: Path) -> str:
    if "$" in path_text:
        return "external_path_unresolved"
    path = _resolve_path(path_text, repo_root)
    if path.exists():
        return "local_exists"
    if path_text.startswith("v9/"):
        return "planned_repo_path_absent"
    return "local_absent"


def _candidate_rows(
    *,
    execution_id: str,
    repo_root: Path,
    candidate_table_path: Path,
    external_probe_note: str,
) -> list[dict[str, str]]:
    source_rows = _read_csv_rows(candidate_table_path)
    rows: list[dict[str, str]] = []
    for row in source_rows:
        candidate_id = row["candidate_id"]
        path_status = _path_status(row["candidate_path"], repo_root)
        if candidate_id == "canonical_v9_payload_target":
            action_decision = (
                "audit_existing_payload"
                if path_status == "local_exists"
                else "blocked_until_source_payload_available"
            )
            evidence_status = "local_file_hashable" if path_status == "local_exists" else "no_local_payload"
        elif candidate_id == "legacy_annotated_h5ad":
            action_decision = (
                "copy_to_canonical_and_audit"
                if path_status == "local_exists"
                else "preferred_route_blocked"
            )
            evidence_status = (
                "preferred_source_available"
                if path_status == "local_exists"
                else "preferred_source_not_available_in_repo_context"
            )
        elif candidate_id == "legacy_labeled_h5ad":
            action_decision = (
                "use_only_for_metadata_recovery"
                if path_status == "local_exists"
                else "fallback_route_blocked"
            )
            evidence_status = (
                "metadata_source_available"
                if path_status == "local_exists"
                else "metadata_source_not_available_in_repo_context"
            )
        else:
            action_decision = "regeneration_route_requires_per_srx_matrices"
            evidence_status = "not_executed_in_this_block"
        notes = row["notes"]
        if external_probe_note and candidate_id in {
            "legacy_annotated_h5ad",
            "legacy_labeled_h5ad",
            "starsolo_per_srx_regeneration",
        }:
            notes = f"{notes} External probe note: {external_probe_note}"
        rows.append(
            {
                "execution_id": execution_id,
                "candidate_id": candidate_id,
                "candidate_path": row["candidate_path"],
                "check_scope": "repo_resolved_path_or_unresolved_external_path",
                "path_status": path_status,
                "evidence_status": evidence_status,
                "action_decision": action_decision,
                "notes": notes,
            }
        )
    return rows


def _regeneration_rows(
    *,
    execution_id: str,
    canonical_payload_path: str,
) -> list[dict[str, str]]:
    return [
        {
            "execution_id": execution_id,
            "step_id": "confirm_per_srx_starsolo_matrices",
            "step_order": "1",
            "required_input": "$SCRATCH_DIR/rrrm1_scrna/OSD-918/starsolo_per_srx/*/Solo.out/GeneFull/filtered",
            "input_status": "unverified_external_payload",
            "output": "per-SRX matrix availability table",
            "gate_status": "blocked_until_external_payload_available",
            "notes": "Required before regenerating the labeled blood AnnData object.",
        },
        {
            "execution_id": execution_id,
            "step_id": "merge_per_srx_h5ad",
            "step_order": "2",
            "required_input": "v2/scripts/rrrm1_merge_per_srx.py|v2/docs/RRRM1_SRX_CONDITION_MAP.csv",
            "input_status": "repo_script_and_condition_map_available",
            "output": "$SCRATCH_DIR/rrrm1_scrna/OSD-918/OSD-918_blood_labeled.h5ad",
            "gate_status": "blocked_until_step_1_passes",
            "notes": "Legacy script can recover condition, animal, age, SRX, OSD, and tissue metadata.",
        },
        {
            "execution_id": execution_id,
            "step_id": "apply_broad_celltype_annotation",
            "step_order": "3",
            "required_input": "v2/scripts/rrrm1_broad_annotate.py|marker sets|processed h5ad",
            "input_status": "repo_script_available_payload_missing",
            "output": "$SCRATCH_DIR/rrrm1_scrna/downstream_initial/annotations/objects/RRRM1_blood_annotated.h5ad",
            "gate_status": "blocked_until_labeled_or_processed_h5ad_available",
            "notes": "Preferred output because the v9 contract requires broad_celltype before runnable status.",
        },
        {
            "execution_id": execution_id,
            "step_id": "normalize_v9_contract_and_copy",
            "step_order": "4",
            "required_input": "annotated h5ad|v9 obs/var/uns contract",
            "input_status": "blocked_no_annotated_h5ad",
            "output": canonical_payload_path,
            "gate_status": "blocked_until_step_3_passes",
            "notes": "Add v9 obs aliases, gene_symbol/feature_id var columns, and spacebio_bench_* uns provenance keys.",
        },
        {
            "execution_id": execution_id,
            "step_id": "rerun_obs_var_audit",
            "step_order": "5",
            "required_input": canonical_payload_path,
            "input_status": "blocked_no_canonical_payload",
            "output": "v9/sc_spaceflight/obs_var_audit_summary.csv",
            "gate_status": "blocked_until_canonical_payload_exists",
            "notes": "Only this step can move the task from skip-only audit to pass/blocker audit.",
        },
    ]


def build_sc_payload_staging_execution(
    *,
    repo_root: str | Path = ".",
    execution_id: str = DEFAULT_SC_PAYLOAD_STAGING_EXECUTION_ID,
    manifest_path: str | Path = (
        "v9/sc_spaceflight/task_manifests/"
        "draft_rrrm1_blood_single_cell_spaceflight.json"
    ),
    candidate_table_path: str | Path = "v9/sc_spaceflight/payload_staging_candidates.csv",
    payload_manifest_path: str | Path = "v9/sc_spaceflight/payload_manifest.csv",
    obs_var_summary_path: str | Path = "v9/sc_spaceflight/obs_var_audit_summary.csv",
    canonical_payload_path: str = (
        "v9/sc_spaceflight/payloads/rrrm1_blood/"
        "OSD-918_blood_rrrm1_bench.h5ad"
    ),
    external_probe_note: str = "",
) -> dict[str, Any]:
    """Build the V9-SC-006 staging execution decision package."""

    root = Path(repo_root).resolve()
    manifest = json.loads(_resolve_path(manifest_path, root).read_text())
    payload_rows = _read_csv_rows(_resolve_path(payload_manifest_path, root))
    obs_var_rows = _read_csv_rows(_resolve_path(obs_var_summary_path, root))
    candidates = _candidate_rows(
        execution_id=execution_id,
        repo_root=root,
        candidate_table_path=_resolve_path(candidate_table_path, root),
        external_probe_note=external_probe_note,
    )
    regeneration_rows = _regeneration_rows(
        execution_id=execution_id,
        canonical_payload_path=canonical_payload_path,
    )
    canonical_status = _path_status(canonical_payload_path, root)
    local_payload_staged = canonical_status == "local_exists"
    preferred_available = any(
        row["candidate_id"] == "legacy_annotated_h5ad"
        and row["path_status"] == "local_exists"
        for row in candidates
    )
    if local_payload_staged:
        decision_status = "canonical_payload_present_ready_for_obs_var_audit"
        selected_route = "audit_existing_canonical_payload"
        next_required_block = "V9-SC-007: single-cell evaluator smoke implementation"
    elif preferred_available:
        decision_status = "legacy_annotated_payload_available_not_copied"
        selected_route = "copy_legacy_annotated_h5ad_then_audit"
        next_required_block = "V9-SC-006b: copy preferred h5ad and rerun audit"
    else:
        decision_status = "payload_staging_execution_blocked_no_candidate_payload"
        selected_route = "prepare_regeneration_from_starsolo_per_srx"
        next_required_block = "V9-SC-006b: locate STARsolo per-SRX matrices or restage annotated h5ad"
    claim_boundary = "payload_staging_execution_no_payload_or_score_claim"
    summary = {
        "execution_id": execution_id,
        "decision_status": decision_status,
        "task_id": manifest["task_id"],
        "canonical_payload_path": canonical_payload_path,
        "canonical_payload_status": canonical_status,
        "selected_route": selected_route,
        "payload_manifest_status": payload_rows[0]["path_status"],
        "obs_var_audit_status": obs_var_rows[0]["decision_status"],
        "local_payload_staged": str(local_payload_staged).lower(),
        "candidate_count": str(len(candidates)),
        "regeneration_step_count": str(len(regeneration_rows)),
        "next_required_block": next_required_block,
        "claim_boundary": claim_boundary,
    }
    review_md = f"""# V9-SC-006 Payload Staging Execution

Status: `{decision_status}`

Task id: `{manifest["task_id"]}`

Claim boundary: `{claim_boundary}`

## Decision

The canonical RRRM-1 blood AnnData payload is still not staged at
`{canonical_payload_path}`. The preferred annotated legacy h5ad is the correct
first route if it can be found, but this execution package does not claim that
any h5ad was copied, downloaded, regenerated, or audited.

Selected route: `{selected_route}`

## Current Evidence

- Canonical payload status: `{canonical_status}`
- Existing payload manifest status: `{payload_rows[0]["path_status"]}`
- Existing obs/var audit status: `{obs_var_rows[0]["decision_status"]}`
- Candidate rows checked: {len(candidates)}
- Regeneration steps fixed: {len(regeneration_rows)}

## Regeneration Gate

If the annotated h5ad cannot be restaged, regeneration must start by confirming
per-SRX STARsolo matrices for `OSD-918` and then rerun the legacy RRRM-1 merge
and broad annotation scripts before applying the v9 obs/var/uns contract.

## Not Claimed

- No local h5ad payload is claimed.
- No payload checksum or shape is claimed beyond the existing missing-payload
  manifest.
- No obs/var pass is claimed while the canonical payload remains absent.
- No evaluator, leaderboard result, or legacy RRRM score promotion is claimed.

## Next Block

Run `{next_required_block}`.
"""
    return {
        "summary": [summary],
        "candidates": candidates,
        "regeneration_steps": regeneration_rows,
        "review_md": review_md,
    }


def write_sc_payload_staging_execution(
    *,
    output_dir: str | Path = "v9/sc_spaceflight",
    repo_root: str | Path = ".",
    execution_id: str = DEFAULT_SC_PAYLOAD_STAGING_EXECUTION_ID,
    manifest_path: str | Path = (
        "v9/sc_spaceflight/task_manifests/"
        "draft_rrrm1_blood_single_cell_spaceflight.json"
    ),
    candidate_table_path: str | Path = "v9/sc_spaceflight/payload_staging_candidates.csv",
    payload_manifest_path: str | Path = "v9/sc_spaceflight/payload_manifest.csv",
    obs_var_summary_path: str | Path = "v9/sc_spaceflight/obs_var_audit_summary.csv",
    canonical_payload_path: str = (
        "v9/sc_spaceflight/payloads/rrrm1_blood/"
        "OSD-918_blood_rrrm1_bench.h5ad"
    ),
    external_probe_note: str = "",
) -> dict[str, Path]:
    """Write the V9-SC-006 RRRM-1 blood payload staging execution package."""

    root = Path(repo_root).resolve()
    package = build_sc_payload_staging_execution(
        repo_root=root,
        execution_id=execution_id,
        manifest_path=manifest_path,
        candidate_table_path=candidate_table_path,
        payload_manifest_path=payload_manifest_path,
        obs_var_summary_path=obs_var_summary_path,
        canonical_payload_path=canonical_payload_path,
        external_probe_note=external_probe_note,
    )
    outdir = _resolve_path(output_dir, root)
    outputs = {
        "summary_csv": outdir / "payload_staging_execution_summary.csv",
        "summary_json": outdir / "payload_staging_execution_summary.json",
        "candidates_csv": outdir / "payload_staging_execution_candidates.csv",
        "candidates_json": outdir / "payload_staging_execution_candidates.json",
        "regeneration_csv": outdir / "payload_regeneration_steps.csv",
        "regeneration_json": outdir / "payload_regeneration_steps.json",
        "review_md": root / "docs" / "V9_SC_PAYLOAD_STAGING_EXECUTION.md",
    }
    _write_csv(
        outputs["summary_csv"],
        package["summary"],
        SC_PAYLOAD_STAGING_EXECUTION_SUMMARY_FIELDS,
    )
    _write_json(outputs["summary_json"], package["summary"])
    _write_csv(
        outputs["candidates_csv"],
        package["candidates"],
        SC_PAYLOAD_STAGING_EXECUTION_CANDIDATE_FIELDS,
    )
    _write_json(outputs["candidates_json"], package["candidates"])
    _write_csv(
        outputs["regeneration_csv"],
        package["regeneration_steps"],
        SC_PAYLOAD_REGENERATION_STEP_FIELDS,
    )
    _write_json(outputs["regeneration_json"], package["regeneration_steps"])
    outputs["review_md"].parent.mkdir(parents=True, exist_ok=True)
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs

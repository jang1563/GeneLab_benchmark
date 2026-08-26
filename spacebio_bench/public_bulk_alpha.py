"""Public bulk alpha readiness gap matrix helpers for SpaceBio-Bench v9."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable


DEFAULT_PUBLIC_BULK_ALPHA_GAP_MATRIX_ID = "v9_public_bulk_alpha_freeze_gap_matrix"
DEFAULT_PUBLIC_BULK_ALPHA_SNAPSHOT_DECISION_ID = (
    "v9_public_bulk_alpha_metadata_snapshot_decision"
)

PUBLIC_BULK_ALPHA_GAP_SUMMARY_FIELDS = [
    "gap_matrix_id",
    "decision_status",
    "bulk_task_count",
    "bulk_fold_count",
    "bulk_source_count",
    "checksum_parsed_source_count",
    "freeze_ready_source_count",
    "baseline_row_count",
    "pass_count",
    "blocker_count",
    "needs_update_count",
    "allowed_claim_count",
    "prohibited_claim_count",
    "next_required_block",
    "claim_boundary",
]

PUBLIC_BULK_ALPHA_GAP_ROW_FIELDS = [
    "gap_matrix_id",
    "gap_id",
    "gap_category",
    "readiness_status",
    "severity",
    "evidence_paths",
    "evidence_summary",
    "owner_action",
    "allowed_current_claim",
    "prohibited_release_claim",
]

PUBLIC_BULK_ALPHA_PAYLOAD_BOUNDARY_FIELDS = [
    "gap_matrix_id",
    "source_id",
    "api_status",
    "audit_status",
    "checksum_manifest_count",
    "parsed_checksum_entries",
    "checksum_payload_matches",
    "freeze_ready",
    "current_boundary",
    "release_blocker",
]

PUBLIC_BULK_ALPHA_CLAIM_BOUNDARY_FIELDS = [
    "gap_matrix_id",
    "claim_id",
    "claim_status",
    "statement",
    "supporting_evidence",
    "disallowed_language",
    "next_allowed_action",
]

PUBLIC_BULK_ALPHA_PACKAGE_UPDATE_FIELDS = [
    "gap_matrix_id",
    "update_id",
    "target_artifact",
    "update_status",
    "required_change",
    "evidence_paths",
    "blocked_by",
]

PUBLIC_BULK_ALPHA_SNAPSHOT_DECISION_SUMMARY_FIELDS = [
    "snapshot_decision_id",
    "decision_status",
    "selected_path",
    "deferred_path",
    "bulk_task_count",
    "bulk_fold_count",
    "bulk_source_count",
    "checksum_parsed_source_count",
    "freeze_ready_source_count",
    "baseline_row_count",
    "gap_blocker_count",
    "metadata_only_allowed",
    "payload_release_allowed",
    "dataset_card_update_required",
    "datapackage_update_required",
    "next_required_block",
    "claim_boundary",
]

PUBLIC_BULK_ALPHA_SNAPSHOT_OPTION_FIELDS = [
    "snapshot_decision_id",
    "path_id",
    "decision",
    "status",
    "rationale",
    "required_wording",
    "remaining_blockers",
    "next_action",
]

PUBLIC_BULK_ALPHA_SNAPSHOT_CLAIM_FIELDS = [
    "snapshot_decision_id",
    "claim_id",
    "claim_status",
    "allowed_language",
    "disallowed_language",
    "supporting_evidence",
    "next_action",
]

PUBLIC_BULK_ALPHA_SNAPSHOT_LANGUAGE_FIELDS = [
    "snapshot_decision_id",
    "snippet_id",
    "target_artifact",
    "section_hint",
    "language_status",
    "snippet_text",
    "must_not_imply",
]

PUBLIC_BULK_ALPHA_SNAPSHOT_ACTION_FIELDS = [
    "snapshot_decision_id",
    "action_id",
    "action_status",
    "target_artifact",
    "required_change",
    "evidence_paths",
    "blocked_by",
    "unblocks",
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


def _pipe_join(values: Iterable[str]) -> str:
    return "|".join(str(value) for value in values if str(value))


def _truthy(value: str) -> bool:
    return value.strip().lower() in {"true", "1", "yes"}


def _gap_row(
    *,
    gap_matrix_id: str,
    gap_id: str,
    gap_category: str,
    readiness_status: str,
    severity: str,
    evidence_paths: Iterable[str],
    evidence_summary: str,
    owner_action: str,
    allowed_current_claim: str,
    prohibited_release_claim: str,
) -> dict[str, str]:
    return {
        "gap_matrix_id": gap_matrix_id,
        "gap_id": gap_id,
        "gap_category": gap_category,
        "readiness_status": readiness_status,
        "severity": severity,
        "evidence_paths": _pipe_join(evidence_paths),
        "evidence_summary": evidence_summary,
        "owner_action": owner_action,
        "allowed_current_claim": allowed_current_claim,
        "prohibited_release_claim": prohibited_release_claim,
    }


def build_public_bulk_alpha_gap_matrix(
    *,
    repo_root: str | Path = ".",
    gap_matrix_id: str = DEFAULT_PUBLIC_BULK_ALPHA_GAP_MATRIX_ID,
) -> dict[str, Any]:
    """Build a public bulk alpha readiness pass/blocker matrix."""

    root = Path(repo_root)
    task_index_path = _resolve_path("v9/task_manifest_index.csv", root)
    source_inventory_path = _resolve_path("v9/source_inventory.csv", root)
    checksum_audit_path = _resolve_path("v9/source_checksum_audit.csv", root)
    baseline_summary_path = _resolve_path(
        "v9/reports/bulk_lomo_baseline_summary.csv",
        root,
    )
    datapackage_path = _resolve_path("v9/datapackage.draft.json", root)
    dataset_card_path = _resolve_path("docs/v9_hf_dataset_card.md", root)
    recenter_summary_path = _resolve_path(
        "v9/reports/recenter_decision/recenter_decision_summary.csv",
        root,
    )

    task_rows = _read_csv_rows(task_index_path)
    source_rows = _read_csv_rows(source_inventory_path)
    checksum_rows = _read_csv_rows(checksum_audit_path)
    baseline_rows = _read_csv_rows(baseline_summary_path)
    recenter_rows = _read_csv_rows(recenter_summary_path)
    datapackage = json.loads(datapackage_path.read_text())
    dataset_card = dataset_card_path.read_text()

    task_count = len(task_rows)
    fold_count = sum(int(row.get("n_folds", "0") or "0") for row in task_rows)
    source_count = len(source_rows)
    checksum_parsed_count = sum(
        1 for row in checksum_rows if row.get("audit_status") == "checksum_manifest_parsed"
    )
    freeze_ready_count = sum(
        1 for row in checksum_rows if _truthy(row.get("freeze_ready", ""))
    )
    baseline_row_count = len(baseline_rows)
    evaluated_baseline_count = sum(
        1 for row in baseline_rows if row.get("status") == "evaluated"
    )
    datapackage_resource_count = len(datapackage.get("resources", []))
    dataset_card_declares_draft = "draft" in dataset_card.lower()
    recenter_selected_bulk = (
        bool(recenter_rows)
        and recenter_rows[0].get("selected_next_lane") == "public_bulk_alpha"
    )

    common_paths = [
        task_index_path.relative_to(root).as_posix(),
        source_inventory_path.relative_to(root).as_posix(),
        checksum_audit_path.relative_to(root).as_posix(),
        baseline_summary_path.relative_to(root).as_posix(),
        datapackage_path.relative_to(root).as_posix(),
        dataset_card_path.relative_to(root).as_posix(),
    ]
    gap_rows = [
        _gap_row(
            gap_matrix_id=gap_matrix_id,
            gap_id="task_manifest_registry",
            gap_category="task_registry",
            readiness_status="pass",
            severity="P0",
            evidence_paths=[task_index_path.relative_to(root).as_posix()],
            evidence_summary=f"{task_count} public bulk task manifests and {fold_count} folds indexed.",
            owner_action="keep generated task index in alpha package",
            allowed_current_claim="public bulk task manifests are indexed",
            prohibited_release_claim="frozen benchmark release",
        ),
        _gap_row(
            gap_matrix_id=gap_matrix_id,
            gap_id="source_inventory_scope",
            gap_category="source_inventory",
            readiness_status="pass",
            severity="P0",
            evidence_paths=[source_inventory_path.relative_to(root).as_posix()],
            evidence_summary=f"{source_count} public bulk source rows scoped to mouse bulk LOMO.",
            owner_action="keep organoid and multispecies draft tracks out of public bulk core",
            allowed_current_claim="public mouse bulk source inventory exists",
            prohibited_release_claim="all v9 draft tracks are in alpha core",
        ),
        _gap_row(
            gap_matrix_id=gap_matrix_id,
            gap_id="checksum_manifest_evidence",
            gap_category="source_integrity",
            readiness_status="pass",
            severity="P0",
            evidence_paths=[checksum_audit_path.relative_to(root).as_posix()],
            evidence_summary=f"{checksum_parsed_count}/{source_count} sources have parsed checksum-manifest evidence.",
            owner_action="retain API/checksum-manifest evidence in alpha package",
            allowed_current_claim="checksum manifests were parsed for public bulk sources",
            prohibited_release_claim="payload files are locally hash-verified",
        ),
        _gap_row(
            gap_matrix_id=gap_matrix_id,
            gap_id="payload_hash_verification",
            gap_category="source_integrity",
            readiness_status="blocker",
            severity="P0",
            evidence_paths=[checksum_audit_path.relative_to(root).as_posix()],
            evidence_summary=f"{freeze_ready_count}/{source_count} sources are freeze_ready=true.",
            owner_action="choose payload-mirror hash plan or metadata-only alpha wording",
            allowed_current_claim="payload hash verification remains pending",
            prohibited_release_claim="frozen payload mirror",
        ),
        _gap_row(
            gap_matrix_id=gap_matrix_id,
            gap_id="baseline_output_evidence",
            gap_category="baseline_outputs",
            readiness_status="pass",
            severity="P1",
            evidence_paths=[baseline_summary_path.relative_to(root).as_posix()],
            evidence_summary=f"{evaluated_baseline_count}/{baseline_row_count} public bulk baseline rows are evaluated.",
            owner_action="keep baseline outputs linked through run manifests",
            allowed_current_claim="simple public bulk baselines are reproduced",
            prohibited_release_claim="state-of-the-art leaderboard",
        ),
        _gap_row(
            gap_matrix_id=gap_matrix_id,
            gap_id="datapackage_descriptor",
            gap_category="package_metadata",
            readiness_status="needs_update",
            severity="P1",
            evidence_paths=[datapackage_path.relative_to(root).as_posix()],
            evidence_summary=f"draft Data Package has {datapackage_resource_count} resources.",
            owner_action="add alpha gap-matrix resources and payload-hash boundary after this block",
            allowed_current_claim="draft Data Package descriptor exists",
            prohibited_release_claim="complete release metadata descriptor",
        ),
        _gap_row(
            gap_matrix_id=gap_matrix_id,
            gap_id="dataset_card_language",
            gap_category="release_language",
            readiness_status="needs_update",
            severity="P1",
            evidence_paths=[dataset_card_path.relative_to(root).as_posix()],
            evidence_summary=(
                "dataset card declares draft status"
                if dataset_card_declares_draft
                else "dataset card draft wording not detected"
            ),
            owner_action="add explicit alpha gap and payload-boundary wording",
            allowed_current_claim="dataset-card draft exists",
            prohibited_release_claim="frozen public release card",
        ),
        _gap_row(
            gap_matrix_id=gap_matrix_id,
            gap_id="track_scope_separation",
            gap_category="scope_boundary",
            readiness_status="pass",
            severity="P1",
            evidence_paths=[
                "v9/human_organoid/",
                "v9/multispecies/",
                source_inventory_path.relative_to(root).as_posix(),
            ],
            evidence_summary="public bulk core remains separate from draft organoid and multispecies tracks.",
            owner_action="do not pull draft extension artifacts into the public bulk alpha core",
            allowed_current_claim="public bulk alpha is mouse bulk LOMO scoped",
            prohibited_release_claim="organoid or multispecies draft tasks are alpha leaderboard tasks",
        ),
        _gap_row(
            gap_matrix_id=gap_matrix_id,
            gap_id="recenter_alignment",
            gap_category="program_sequence",
            readiness_status="pass" if recenter_selected_bulk else "needs_update",
            severity="P1",
            evidence_paths=[recenter_summary_path.relative_to(root).as_posix()],
            evidence_summary="V9-REFOCUS-001 selected public bulk alpha as the next active lane.",
            owner_action="continue V9-BULK-ALPHA blocks before V9-SC-001",
            allowed_current_claim="program has recentered on public bulk alpha readiness",
            prohibited_release_claim="single-cell flagship has started",
        ),
        _gap_row(
            gap_matrix_id=gap_matrix_id,
            gap_id="minimal_alpha_snapshot_decision",
            gap_category="release_decision",
            readiness_status="blocker",
            severity="P0",
            evidence_paths=common_paths,
            evidence_summary="metadata-only alpha snapshot decision is not yet recorded.",
            owner_action="decide metadata-only alpha versus payload-mirror prerequisite",
            allowed_current_claim="alpha decision is pending",
            prohibited_release_claim="alpha release approved",
        ),
    ]

    payload_rows = [
        {
            "gap_matrix_id": gap_matrix_id,
            "source_id": row["source_id"],
            "api_status": row["api_status"],
            "audit_status": row["audit_status"],
            "checksum_manifest_count": row["checksum_manifest_count"],
            "parsed_checksum_entries": row["parsed_checksum_entries"],
            "checksum_payload_matches": row["checksum_payload_matches"],
            "freeze_ready": row["freeze_ready"],
            "current_boundary": (
                "checksum_manifest_evidence_only"
                if not _truthy(row["freeze_ready"])
                else "payload_hash_verified"
            ),
            "release_blocker": (
                "local_payload_hash_verification_pending"
                if not _truthy(row["freeze_ready"])
                else ""
            ),
        }
        for row in checksum_rows
    ]
    claim_rows = [
        {
            "gap_matrix_id": gap_matrix_id,
            "claim_id": "public_bulk_metadata_scaffold_exists",
            "claim_status": "allowed_current_alpha_gap",
            "statement": "The public bulk scaffold has task, source, package, and baseline metadata.",
            "supporting_evidence": _pipe_join(common_paths),
            "disallowed_language": "frozen public benchmark release",
            "next_allowed_action": "describe as alpha scaffold with explicit blockers",
        },
        {
            "gap_matrix_id": gap_matrix_id,
            "claim_id": "checksum_manifests_parsed",
            "claim_status": "allowed_current_alpha_gap",
            "statement": "OSDR checksum manifests have been parsed for all public bulk sources.",
            "supporting_evidence": checksum_audit_path.relative_to(root).as_posix(),
            "disallowed_language": "payload hash verified",
            "next_allowed_action": "separate manifest parsing from local payload hashing",
        },
        {
            "gap_matrix_id": gap_matrix_id,
            "claim_id": "baseline_outputs_reproduced",
            "claim_status": "allowed_current_alpha_gap",
            "statement": "Simple baseline outputs are present for the public bulk task set.",
            "supporting_evidence": baseline_summary_path.relative_to(root).as_posix(),
            "disallowed_language": "leaderboard or model superiority claim",
            "next_allowed_action": "use as reproducibility evidence only",
        },
        {
            "gap_matrix_id": gap_matrix_id,
            "claim_id": "no_payload_freeze_claim",
            "claim_status": "prohibited_current_alpha_gap",
            "statement": "The current package must not claim local payload hash freeze.",
            "supporting_evidence": checksum_audit_path.relative_to(root).as_posix(),
            "disallowed_language": "frozen payload mirror",
            "next_allowed_action": "add payload hash plan or metadata-only alpha decision",
        },
        {
            "gap_matrix_id": gap_matrix_id,
            "claim_id": "no_extension_track_promotion",
            "claim_status": "prohibited_current_alpha_gap",
            "statement": "Draft organoid and multispecies outputs are not public bulk alpha core tasks.",
            "supporting_evidence": "v9/human_organoid/|v9/multispecies/",
            "disallowed_language": "all v9 draft tracks are public alpha leaderboard tasks",
            "next_allowed_action": "keep extension tracks in separate draft lanes",
        },
        {
            "gap_matrix_id": gap_matrix_id,
            "claim_id": "no_alpha_release_approval",
            "claim_status": "prohibited_current_alpha_gap",
            "statement": "A public alpha release approval is not recorded by this gap matrix.",
            "supporting_evidence": "v9/reports/public_bulk_alpha_gap_matrix/",
            "disallowed_language": "release approved",
            "next_allowed_action": "run V9-BULK-ALPHA-002 decision block",
        },
    ]
    update_rows = [
        {
            "gap_matrix_id": gap_matrix_id,
            "update_id": "add_gap_matrix_to_datapackage",
            "target_artifact": datapackage_path.relative_to(root).as_posix(),
            "update_status": "pending_after_gap_matrix",
            "required_change": "add public bulk alpha gap-matrix resources after review",
            "evidence_paths": "v9/reports/public_bulk_alpha_gap_matrix/",
            "blocked_by": "V9-BULK-ALPHA-001 review completion",
        },
        {
            "gap_matrix_id": gap_matrix_id,
            "update_id": "add_payload_boundary_to_dataset_card",
            "target_artifact": dataset_card_path.relative_to(root).as_posix(),
            "update_status": "pending_after_gap_matrix",
            "required_change": "state checksum-manifest evidence versus local payload hash verification explicitly",
            "evidence_paths": checksum_audit_path.relative_to(root).as_posix(),
            "blocked_by": "V9-BULK-ALPHA-002 metadata-only alpha decision",
        },
        {
            "gap_matrix_id": gap_matrix_id,
            "update_id": "keep_extension_tracks_out_of_bulk_alpha",
            "target_artifact": source_inventory_path.relative_to(root).as_posix(),
            "update_status": "guarded_no_change",
            "required_change": "do not merge organoid or multispecies draft sources into public bulk inventory",
            "evidence_paths": "v9/human_organoid/|v9/multispecies/",
            "blocked_by": "",
        },
    ]

    pass_count = sum(1 for row in gap_rows if row["readiness_status"] == "pass")
    blocker_count = sum(1 for row in gap_rows if row["readiness_status"] == "blocker")
    needs_update_count = sum(
        1 for row in gap_rows if row["readiness_status"] == "needs_update"
    )
    allowed_claim_count = sum(
        1 for row in claim_rows if row["claim_status"] == "allowed_current_alpha_gap"
    )
    prohibited_claim_count = sum(
        1 for row in claim_rows if row["claim_status"] == "prohibited_current_alpha_gap"
    )
    claim_boundary = "public_bulk_alpha_gap_matrix_no_release_approval"
    summary = {
        "gap_matrix_id": gap_matrix_id,
        "decision_status": "metadata_alpha_gap_matrix_ready_payload_hash_blocked",
        "bulk_task_count": str(task_count),
        "bulk_fold_count": str(fold_count),
        "bulk_source_count": str(source_count),
        "checksum_parsed_source_count": str(checksum_parsed_count),
        "freeze_ready_source_count": str(freeze_ready_count),
        "baseline_row_count": str(baseline_row_count),
        "pass_count": str(pass_count),
        "blocker_count": str(blocker_count),
        "needs_update_count": str(needs_update_count),
        "allowed_claim_count": str(allowed_claim_count),
        "prohibited_claim_count": str(prohibited_claim_count),
        "next_required_block": (
            "V9-BULK-ALPHA-002: metadata-only alpha snapshot decision"
        ),
        "claim_boundary": claim_boundary,
    }
    review_md = f"""# V9-BULK-ALPHA-001 Public Bulk Alpha Freeze-Gap Matrix Review

Status: `{summary["decision_status"]}`

Claim boundary: `{claim_boundary}`

## Decision

The public bulk alpha scaffold is close enough to stay as the active lane, but
not ready for frozen release language. The strongest blocker is payload hash
verification: `{freeze_ready_count}/{source_count}` public bulk sources are
`freeze_ready=true`, even though `{checksum_parsed_count}/{source_count}` have
parsed OSDR checksum-manifest evidence.

## Matrix Summary

- Pass rows: {pass_count}
- Blocker rows: {blocker_count}
- Needs-update rows: {needs_update_count}
- Allowed current claims: {allowed_claim_count}
- Prohibited release claims: {prohibited_claim_count}

## Current Allowed Language

It is safe to say that the public bulk scaffold has indexed task manifests,
source inventory, OSDR API/checksum-manifest evidence, draft Data Package
metadata, dataset-card draft language, and reproduced simple baselines.

## Current Prohibited Language

Do not claim a frozen public benchmark release, locally verified payload mirror,
complete release metadata descriptor, extension-track leaderboard, or approved
alpha release.

## Next Block

Run `V9-BULK-ALPHA-002: metadata-only alpha snapshot decision`. That block
should decide whether the project can publish a metadata-only alpha snapshot
with explicit payload-hash blockers, or whether payload mirroring and local hash
verification must precede any alpha wording.

## External Guidance Anchors

- OSDR API evidence can establish source/file-list traceability, but local
  payload verification is a separate benchmark-release claim.
- Frictionless Data Package resources can describe current tables and reports,
  but package metadata must not imply release completeness while payload and
  license/citation boundaries remain unresolved.
"""
    return {
        "summary": [summary],
        "gap_rows": gap_rows,
        "payload_rows": payload_rows,
        "claim_rows": claim_rows,
        "update_rows": update_rows,
        "review_md": review_md,
    }


def write_public_bulk_alpha_gap_matrix(
    *,
    output_dir: str | Path = "v9/reports/public_bulk_alpha_gap_matrix",
    repo_root: str | Path = ".",
    gap_matrix_id: str = DEFAULT_PUBLIC_BULK_ALPHA_GAP_MATRIX_ID,
) -> dict[str, Path]:
    """Write public bulk alpha gap matrix artifacts."""

    root = Path(repo_root)
    package = build_public_bulk_alpha_gap_matrix(
        repo_root=root,
        gap_matrix_id=gap_matrix_id,
    )
    outdir = _resolve_path(output_dir, root)
    outputs = {
        "summary_csv": outdir / "public_bulk_alpha_gap_summary.csv",
        "summary_json": outdir / "public_bulk_alpha_gap_summary.json",
        "gap_csv": outdir / "public_bulk_alpha_gap_matrix.csv",
        "gap_json": outdir / "public_bulk_alpha_gap_matrix.json",
        "payload_csv": outdir / "payload_hash_boundary.csv",
        "payload_json": outdir / "payload_hash_boundary.json",
        "claim_csv": outdir / "public_bulk_alpha_claim_boundary.csv",
        "claim_json": outdir / "public_bulk_alpha_claim_boundary.json",
        "update_csv": outdir / "package_update_plan.csv",
        "update_json": outdir / "package_update_plan.json",
        "review_md": root / "docs" / "V9_PUBLIC_BULK_ALPHA_FREEZE_GAP_MATRIX_REVIEW.md",
    }
    _write_csv(
        outputs["summary_csv"],
        package["summary"],
        PUBLIC_BULK_ALPHA_GAP_SUMMARY_FIELDS,
    )
    _write_json(outputs["summary_json"], package["summary"])
    _write_csv(
        outputs["gap_csv"],
        package["gap_rows"],
        PUBLIC_BULK_ALPHA_GAP_ROW_FIELDS,
    )
    _write_json(outputs["gap_json"], package["gap_rows"])
    _write_csv(
        outputs["payload_csv"],
        package["payload_rows"],
        PUBLIC_BULK_ALPHA_PAYLOAD_BOUNDARY_FIELDS,
    )
    _write_json(outputs["payload_json"], package["payload_rows"])
    _write_csv(
        outputs["claim_csv"],
        package["claim_rows"],
        PUBLIC_BULK_ALPHA_CLAIM_BOUNDARY_FIELDS,
    )
    _write_json(outputs["claim_json"], package["claim_rows"])
    _write_csv(
        outputs["update_csv"],
        package["update_rows"],
        PUBLIC_BULK_ALPHA_PACKAGE_UPDATE_FIELDS,
    )
    _write_json(outputs["update_json"], package["update_rows"])
    outputs["review_md"].parent.mkdir(parents=True, exist_ok=True)
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs


def build_public_bulk_alpha_snapshot_decision(
    *,
    repo_root: str | Path = ".",
    gap_report_dir: str | Path = "v9/reports/public_bulk_alpha_gap_matrix",
    snapshot_decision_id: str = DEFAULT_PUBLIC_BULK_ALPHA_SNAPSHOT_DECISION_ID,
) -> dict[str, Any]:
    """Build the metadata-only alpha versus payload-mirror-first decision."""

    root = Path(repo_root)
    report_dir = _resolve_path(gap_report_dir, root)
    gap_summary_path = report_dir / "public_bulk_alpha_gap_summary.csv"
    gap_matrix_path = report_dir / "public_bulk_alpha_gap_matrix.csv"
    payload_boundary_path = report_dir / "payload_hash_boundary.csv"
    claim_boundary_path = report_dir / "public_bulk_alpha_claim_boundary.csv"
    package_update_path = report_dir / "package_update_plan.csv"
    gap_review_path = root / "docs" / "V9_PUBLIC_BULK_ALPHA_FREEZE_GAP_MATRIX_REVIEW.md"
    dataset_card_path = root / "docs" / "v9_hf_dataset_card.md"
    datapackage_path = root / "v9" / "datapackage.draft.json"

    gap_summary_rows = _read_csv_rows(gap_summary_path)
    if not gap_summary_rows:
        raise ValueError(f"no gap summary rows found in {gap_summary_path}")
    gap_summary = gap_summary_rows[0]
    gap_rows = _read_csv_rows(gap_matrix_path)
    payload_rows = _read_csv_rows(payload_boundary_path)
    prior_claim_rows = _read_csv_rows(claim_boundary_path)
    package_update_rows = _read_csv_rows(package_update_path)

    source_count = int(gap_summary["bulk_source_count"])
    freeze_ready_count = int(gap_summary["freeze_ready_source_count"])
    checksum_parsed_count = int(gap_summary["checksum_parsed_source_count"])
    blocker_count = int(gap_summary["blocker_count"])
    all_payload_hashes_pending = bool(payload_rows) and all(
        row["current_boundary"] == "checksum_manifest_evidence_only"
        and row["release_blocker"] == "local_payload_hash_verification_pending"
        for row in payload_rows
    )
    has_payload_blocker = any(
        row["gap_id"] == "payload_hash_verification"
        and row["readiness_status"] == "blocker"
        for row in gap_rows
    )
    metadata_only_allowed = (
        checksum_parsed_count == source_count
        and freeze_ready_count == 0
        and all_payload_hashes_pending
        and has_payload_blocker
    )
    payload_release_allowed = freeze_ready_count == source_count and source_count > 0

    selected_path = (
        "metadata_only_alpha_snapshot"
        if metadata_only_allowed and not payload_release_allowed
        else "payload_mirror_first"
    )
    deferred_path = (
        "payload_mirror_first"
        if selected_path == "metadata_only_alpha_snapshot"
        else "metadata_only_alpha_snapshot"
    )
    claim_boundary = "metadata_only_public_bulk_alpha_no_payload_release"
    next_required_block = (
        "V9-BULK-ALPHA-003: dataset card and Data Package alpha boundary update"
    )

    option_rows = [
        {
            "snapshot_decision_id": snapshot_decision_id,
            "path_id": "metadata_only_alpha_snapshot",
            "decision": "selected",
            "status": "allowed_with_explicit_blockers",
            "rationale": (
                "The public bulk task, source, checksum-manifest, baseline, and "
                "package metadata spine is indexed, while payload hashes remain "
                "clearly blocked rather than silently implied."
            ),
            "required_wording": (
                "metadata-only alpha snapshot; checksum-manifest evidence parsed; "
                "no local payload mirror; no frozen payload hash verification"
            ),
            "remaining_blockers": (
                "dataset-card alpha-boundary wording|Data Package decision "
                "resources|payload release still requires local hash verification"
            ),
            "next_action": "update dataset-card and Data Package boundary language",
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "path_id": "payload_mirror_first",
            "decision": "deferred",
            "status": "valid_but_not_required_before_metadata_alpha",
            "rationale": (
                "Payload mirroring is the right prerequisite for a frozen payload "
                "release, but it is not required to publish a transparent metadata "
                "alpha that says the payload hashes are not locally verified."
            ),
            "required_wording": (
                "use only for a future frozen payload release after local mirroring "
                "and hash comparison are complete"
            ),
            "remaining_blockers": (
                "mirror upstream payloads|compare local hashes for 22 public bulk "
                "sources|record citation/license/version boundary"
            ),
            "next_action": "keep as a later payload-freeze block after metadata alpha",
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "path_id": "no_alpha_until_payload_frozen",
            "decision": "rejected_for_current_sequence",
            "status": "too_conservative_for_metadata_scaffold",
            "rationale": (
                "This would stall the v9 alpha spine even though the current "
                "artifact can be bounded as metadata only and can help reviewers "
                "inspect task/source/provenance decisions."
            ),
            "required_wording": "not selected",
            "remaining_blockers": "would defer all alpha wording until payload freeze",
            "next_action": "do not use this as the active v9 sequence",
        },
    ]

    evidence_paths = [
        gap_summary_path.relative_to(root).as_posix(),
        gap_matrix_path.relative_to(root).as_posix(),
        payload_boundary_path.relative_to(root).as_posix(),
        claim_boundary_path.relative_to(root).as_posix(),
        package_update_path.relative_to(root).as_posix(),
        gap_review_path.relative_to(root).as_posix(),
    ]
    claim_rows = [
        {
            "snapshot_decision_id": snapshot_decision_id,
            "claim_id": "metadata_alpha_label",
            "claim_status": "allowed_metadata_alpha",
            "allowed_language": "SpaceBio-Bench v9 public bulk metadata alpha",
            "disallowed_language": "SpaceBio-Bench v9 frozen public benchmark release",
            "supporting_evidence": _pipe_join(evidence_paths),
            "next_action": "use as the short public status label",
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "claim_id": "metadata_scaffold_contents",
            "claim_status": "allowed_metadata_alpha",
            "allowed_language": (
                "The snapshot includes task manifests, source inventory, "
                "checksum-manifest evidence, baseline summaries, and provenance "
                "reports for the public mouse bulk LOMO scaffold."
            ),
            "disallowed_language": "The snapshot distributes all upstream payload data.",
            "supporting_evidence": _pipe_join(evidence_paths),
            "next_action": "mirror this wording in the dataset card overview",
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "claim_id": "checksum_manifest_evidence",
            "claim_status": "allowed_metadata_alpha",
            "allowed_language": (
                "OSDR file-list and checksum-manifest evidence has been parsed "
                "for all 22 public bulk source rows."
            ),
            "disallowed_language": "All payload files are locally hash verified.",
            "supporting_evidence": payload_boundary_path.relative_to(root).as_posix(),
            "next_action": "keep checksum-manifest evidence separate from payload freeze",
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "claim_id": "upstream_source_credit",
            "claim_status": "allowed_metadata_alpha",
            "allowed_language": (
                "The snapshot references upstream NASA OSDR datasets and should "
                "credit OSDR/source studies rather than redistributing payloads."
            ),
            "disallowed_language": "SpaceBio-Bench is the authoritative source archive.",
            "supporting_evidence": "https://science.nasa.gov/reference/osdr-faq/",
            "next_action": "keep OSDR citation/source-credit wording visible",
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "claim_id": "payload_release_claim",
            "claim_status": "prohibited_payload_release",
            "allowed_language": "Payload release remains pending local mirroring and hash verification.",
            "disallowed_language": "frozen payload mirror",
            "supporting_evidence": payload_boundary_path.relative_to(root).as_posix(),
            "next_action": "block release language until freeze_ready=true sources exist",
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "claim_id": "local_hash_verified_claim",
            "claim_status": "prohibited_payload_release",
            "allowed_language": "Local payload-hash verification is not yet complete.",
            "disallowed_language": "locally hash-verified data bundle",
            "supporting_evidence": payload_boundary_path.relative_to(root).as_posix(),
            "next_action": "defer to a future payload-mirror freeze block",
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "claim_id": "leaderboard_or_model_superiority",
            "claim_status": "prohibited_payload_release",
            "allowed_language": "Baseline outputs are reproducibility checks for the scaffold.",
            "disallowed_language": "state-of-the-art leaderboard",
            "supporting_evidence": "v9/reports/bulk_lomo_baseline_summary.csv",
            "next_action": "keep baseline language descriptive",
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "claim_id": "extension_track_promotion",
            "claim_status": "prohibited_payload_release",
            "allowed_language": "Organoid and multispecies artifacts remain draft extension lanes.",
            "disallowed_language": "all v9 draft tracks are public bulk alpha tasks",
            "supporting_evidence": "v9/human_organoid/|v9/multispecies/",
            "next_action": "keep extension tracks out of the public bulk alpha core",
        },
    ]

    language_rows = [
        {
            "snapshot_decision_id": snapshot_decision_id,
            "snippet_id": "short_label",
            "target_artifact": dataset_card_path.relative_to(root).as_posix(),
            "section_hint": "title/status",
            "language_status": "allowed",
            "snippet_text": "SpaceBio-Bench v9 public bulk metadata alpha",
            "must_not_imply": "frozen payload release",
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "snippet_id": "status_paragraph",
            "target_artifact": dataset_card_path.relative_to(root).as_posix(),
            "section_hint": "Dataset Status",
            "language_status": "allowed",
            "snippet_text": (
                "This metadata-only alpha snapshot documents task manifests, "
                "source inventory, checksum-manifest evidence, baseline summaries, "
                "and provenance reports for the public mouse bulk LOMO scaffold. "
                "It is not a frozen payload release and does not include a local "
                "payload mirror or local payload-hash verification."
            ),
            "must_not_imply": "payload files are bundled or hash verified",
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "snippet_id": "limitations_payload_boundary",
            "target_artifact": dataset_card_path.relative_to(root).as_posix(),
            "section_hint": "Limitations",
            "language_status": "required",
            "snippet_text": (
                "The current evidence boundary is OSDR API file-list and "
                "checksum-manifest parsing. Local payload mirroring and payload "
                "hash verification remain future work before any frozen data "
                "release wording."
            ),
            "must_not_imply": "release-ready benchmark payload",
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "snippet_id": "datapackage_resource_boundary",
            "target_artifact": datapackage_path.relative_to(root).as_posix(),
            "section_hint": "resources/descriptions",
            "language_status": "required",
            "snippet_text": (
                "Add gap-matrix and snapshot-decision resources as metadata "
                "resources; do not add upstream payload files as locally bundled "
                "resources until a payload mirror and hash manifest exist."
            ),
            "must_not_imply": "complete Data Package release descriptor",
        },
    ]

    action_rows = [
        {
            "snapshot_decision_id": snapshot_decision_id,
            "action_id": "update_dataset_card_alpha_boundary",
            "action_status": "pending_next_block",
            "target_artifact": dataset_card_path.relative_to(root).as_posix(),
            "required_change": (
                "add metadata-only alpha status paragraph and explicit payload "
                "hash limitation wording"
            ),
            "evidence_paths": _pipe_join(evidence_paths),
            "blocked_by": "",
            "unblocks": next_required_block,
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "action_id": "update_datapackage_snapshot_resources",
            "action_status": "pending_next_block",
            "target_artifact": datapackage_path.relative_to(root).as_posix(),
            "required_change": (
                "add public bulk gap-matrix and snapshot-decision tables as "
                "metadata resources"
            ),
            "evidence_paths": "v9/reports/public_bulk_alpha_gap_matrix/|v9/reports/public_bulk_alpha_snapshot_decision/",
            "blocked_by": "",
            "unblocks": next_required_block,
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "action_id": "defer_payload_mirror_freeze",
            "action_status": "deferred_after_metadata_alpha",
            "target_artifact": "future payload mirror manifest",
            "required_change": (
                "mirror upstream payloads and compare local file hashes before "
                "any payload release claim"
            ),
            "evidence_paths": payload_boundary_path.relative_to(root).as_posix(),
            "blocked_by": "metadata alpha boundary should be recorded first",
            "unblocks": "future payload-freeze release candidate",
        },
        {
            "snapshot_decision_id": snapshot_decision_id,
            "action_id": "keep_extension_tracks_separate",
            "action_status": "guarded_no_change",
            "target_artifact": "v9/source_inventory.csv",
            "required_change": (
                "do not merge organoid or multispecies draft sources into the "
                "public bulk alpha core"
            ),
            "evidence_paths": "v9/human_organoid/|v9/multispecies/",
            "blocked_by": "",
            "unblocks": "purpose-drift guard",
        },
    ]

    allowed_claim_count = sum(
        1 for row in claim_rows if row["claim_status"] == "allowed_metadata_alpha"
    )
    prohibited_claim_count = sum(
        1 for row in claim_rows if row["claim_status"] == "prohibited_payload_release"
    )
    decision_status = (
        "metadata_only_alpha_snapshot_allowed_with_payload_blockers"
        if selected_path == "metadata_only_alpha_snapshot"
        else "payload_mirror_first_required_before_alpha"
    )
    summary = {
        "snapshot_decision_id": snapshot_decision_id,
        "decision_status": decision_status,
        "selected_path": selected_path,
        "deferred_path": deferred_path,
        "bulk_task_count": gap_summary["bulk_task_count"],
        "bulk_fold_count": gap_summary["bulk_fold_count"],
        "bulk_source_count": gap_summary["bulk_source_count"],
        "checksum_parsed_source_count": gap_summary["checksum_parsed_source_count"],
        "freeze_ready_source_count": gap_summary["freeze_ready_source_count"],
        "baseline_row_count": gap_summary["baseline_row_count"],
        "gap_blocker_count": str(blocker_count),
        "metadata_only_allowed": str(metadata_only_allowed).lower(),
        "payload_release_allowed": str(payload_release_allowed).lower(),
        "dataset_card_update_required": "true",
        "datapackage_update_required": "true",
        "next_required_block": next_required_block,
        "claim_boundary": claim_boundary,
    }

    review_md = f"""# V9-BULK-ALPHA-002 Metadata-Only Alpha Snapshot Decision

Status: `{summary["decision_status"]}`

Selected path: `{selected_path}`

Deferred path: `{deferred_path}`

Claim boundary: `{claim_boundary}`

## Decision

Proceed with a metadata-only alpha snapshot for the public bulk lane, with
explicit payload-hash blockers. This is not a frozen payload release. It is a
bounded review surface for task manifests, source inventory, OSDR API and
checksum-manifest evidence, baseline summaries, and provenance reports.

The payload-release path remains blocked because `{freeze_ready_count}/{source_count}`
public bulk sources are locally payload-hash verified, while
`{checksum_parsed_count}/{source_count}` sources have parsed checksum-manifest
evidence.

## Option Comparison

| Path | Decision | Status |
| --- | --- | --- |
| metadata-only alpha snapshot | selected | allowed with explicit blockers |
| payload mirror first | deferred | valid for future payload release, not required before metadata alpha |
| no alpha until payload frozen | rejected | too conservative for metadata scaffold |

## Allowed Language

- `SpaceBio-Bench v9 public bulk metadata alpha`
- The snapshot documents public mouse bulk LOMO task/source/provenance metadata.
- OSDR file-list and checksum-manifest evidence has been parsed for all 22
  public bulk source rows.
- Payload mirroring and local payload-hash verification remain pending.

## Prohibited Language

- Frozen public benchmark release.
- Frozen payload mirror.
- Locally hash-verified data bundle.
- DOI/archive release, complete release Data Package, or leaderboard claim.
- Organoid or multispecies draft tracks as public bulk alpha core tasks.

## External Guidance Anchors

- Hugging Face dataset cards are README/metadata surfaces meant to help users
  understand dataset contents, context, and responsible use:
  https://huggingface.co/docs/hub/datasets-cards
- Frictionless Data Package descriptors separate package metadata from resource
  entries and can describe metadata resources without implying a local payload
  mirror:
  https://specs.frictionlessdata.io/data-package/
- OSDR API file-list and metadata endpoints support source/file traceability,
  while local benchmark payload hashing remains a separate project claim:
  https://visualization.osdr.nasa.gov/biodata/api/
- NASA OSDR should remain the credited upstream source for space biology data:
  https://science.nasa.gov/reference/osdr-faq/

## Next Block

Run `{next_required_block}`. That block should update
`docs/v9_hf_dataset_card.md` and `v9/datapackage.draft.json` using the claim
boundary in this decision package.
"""
    return {
        "summary": [summary],
        "option_rows": option_rows,
        "claim_rows": claim_rows,
        "language_rows": language_rows,
        "action_rows": action_rows,
        "prior_claim_rows": prior_claim_rows,
        "package_update_rows": package_update_rows,
        "review_md": review_md,
    }


def write_public_bulk_alpha_snapshot_decision(
    *,
    output_dir: str | Path = "v9/reports/public_bulk_alpha_snapshot_decision",
    repo_root: str | Path = ".",
    gap_report_dir: str | Path = "v9/reports/public_bulk_alpha_gap_matrix",
    snapshot_decision_id: str = DEFAULT_PUBLIC_BULK_ALPHA_SNAPSHOT_DECISION_ID,
) -> dict[str, Path]:
    """Write metadata-only alpha snapshot decision artifacts."""

    root = Path(repo_root)
    package = build_public_bulk_alpha_snapshot_decision(
        repo_root=root,
        gap_report_dir=gap_report_dir,
        snapshot_decision_id=snapshot_decision_id,
    )
    outdir = _resolve_path(output_dir, root)
    outputs = {
        "summary_csv": outdir / "snapshot_decision_summary.csv",
        "summary_json": outdir / "snapshot_decision_summary.json",
        "option_csv": outdir / "snapshot_option_matrix.csv",
        "option_json": outdir / "snapshot_option_matrix.json",
        "claim_csv": outdir / "snapshot_claim_boundary.csv",
        "claim_json": outdir / "snapshot_claim_boundary.json",
        "language_csv": outdir / "snapshot_language_snippets.csv",
        "language_json": outdir / "snapshot_language_snippets.json",
        "action_csv": outdir / "snapshot_next_actions.csv",
        "action_json": outdir / "snapshot_next_actions.json",
        "review_md": root
        / "docs"
        / "V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md",
    }
    _write_csv(
        outputs["summary_csv"],
        package["summary"],
        PUBLIC_BULK_ALPHA_SNAPSHOT_DECISION_SUMMARY_FIELDS,
    )
    _write_json(outputs["summary_json"], package["summary"])
    _write_csv(
        outputs["option_csv"],
        package["option_rows"],
        PUBLIC_BULK_ALPHA_SNAPSHOT_OPTION_FIELDS,
    )
    _write_json(outputs["option_json"], package["option_rows"])
    _write_csv(
        outputs["claim_csv"],
        package["claim_rows"],
        PUBLIC_BULK_ALPHA_SNAPSHOT_CLAIM_FIELDS,
    )
    _write_json(outputs["claim_json"], package["claim_rows"])
    _write_csv(
        outputs["language_csv"],
        package["language_rows"],
        PUBLIC_BULK_ALPHA_SNAPSHOT_LANGUAGE_FIELDS,
    )
    _write_json(outputs["language_json"], package["language_rows"])
    _write_csv(
        outputs["action_csv"],
        package["action_rows"],
        PUBLIC_BULK_ALPHA_SNAPSHOT_ACTION_FIELDS,
    )
    _write_json(outputs["action_json"], package["action_rows"])
    outputs["review_md"].parent.mkdir(parents=True, exist_ok=True)
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs

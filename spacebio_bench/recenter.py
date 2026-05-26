"""Recenter decision helpers for SpaceBio-Bench v9."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable


DEFAULT_V9_RECENTER_DECISION_ID = "v9_refocus_001_post_osd120_recenter"

RECENTER_DECISION_SUMMARY_FIELDS = [
    "recenter_id",
    "decision_status",
    "selected_next_lane",
    "selected_next_block",
    "selected_next_block_goal",
    "public_bulk_readiness_score",
    "single_cell_readiness_score",
    "bulk_task_count",
    "bulk_fold_count",
    "bulk_source_count",
    "bulk_api_ok_source_count",
    "bulk_checksum_parsed_source_count",
    "bulk_freeze_ready_source_count",
    "bulk_baseline_row_count",
    "single_cell_asset_count",
    "single_cell_h5ad_count",
    "single_cell_v9_manifest_count",
    "osd120_branch_status",
    "claim_boundary",
]

RECENTER_CANDIDATE_FIELDS = [
    "recenter_id",
    "candidate_id",
    "lane_label",
    "decision",
    "readiness_status",
    "readiness_score",
    "evidence_items",
    "blocking_gaps",
    "external_method_anchor",
    "next_block_if_selected",
    "rationale",
]

RECENTER_NEXT_BLOCK_ACTION_FIELDS = [
    "recenter_id",
    "block_id",
    "action_id",
    "lane_id",
    "objective",
    "expected_files",
    "expected_tests",
    "done_when",
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


def _single_cell_asset_paths(repo_root: Path) -> list[str]:
    keywords = (
        "singlecell",
        "single-cell",
        "single_cell",
        "scrna",
        "snrnaseq",
        "rrrm",
        "anndata",
        "h5ad",
        "scanpy",
    )
    search_roots = [
        repo_root / "v2",
        repo_root / "v3",
        repo_root / "docs",
        repo_root / "scripts",
    ]
    paths: list[str] = []
    for search_root in search_roots:
        if not search_root.exists():
            continue
        for path in search_root.rglob("*"):
            if not path.is_file():
                continue
            rel = path.relative_to(repo_root).as_posix()
            lowered = rel.lower()
            if any(keyword in lowered for keyword in keywords):
                paths.append(rel)
    return sorted(set(paths))


def build_v9_recenter_decision(
    *,
    repo_root: str | Path = ".",
    recenter_id: str = DEFAULT_V9_RECENTER_DECISION_ID,
) -> dict[str, Any]:
    """Build the V9-REFOCUS-001 public-bulk versus single-cell decision."""

    root = Path(repo_root)
    task_index_path = _resolve_path("v9/task_manifest_index.csv", root)
    source_inventory_path = _resolve_path("v9/source_inventory.csv", root)
    checksum_audit_path = _resolve_path("v9/source_checksum_audit.csv", root)
    baseline_summary_path = _resolve_path(
        "v9/reports/bulk_lomo_baseline_summary.csv",
        root,
    )
    osd120_closeout_path = _resolve_path(
        "v9/multispecies/reports/interaction_diagnostic_metadata_release_note/"
        "diagnostic_metadata_release_summary.csv",
        root,
    )
    datapackage_path = _resolve_path("v9/datapackage.draft.json", root)
    dataset_card_path = _resolve_path("docs/v9_hf_dataset_card.md", root)
    sc_manifest_dir = _resolve_path("v9/sc_spaceflight/task_manifests", root)
    profile_path = _resolve_path("spacebio_bench/profiles.py", root)

    task_rows = _read_csv_rows(task_index_path)
    source_rows = _read_csv_rows(source_inventory_path)
    checksum_rows = _read_csv_rows(checksum_audit_path)
    baseline_rows = _read_csv_rows(baseline_summary_path)
    osd120_rows = _read_csv_rows(osd120_closeout_path)
    osd120_status = (
        osd120_rows[0]["branch_closeout_status"]
        if osd120_rows
        else "missing_osd120_closeout_summary"
    )

    bulk_task_count = len(task_rows)
    bulk_fold_count = sum(int(row.get("n_folds", "0") or "0") for row in task_rows)
    bulk_source_count = len(source_rows)
    bulk_api_ok_count = sum(1 for row in checksum_rows if row.get("api_status") == "ok")
    bulk_checksum_parsed_count = sum(
        1 for row in checksum_rows if row.get("audit_status") == "checksum_manifest_parsed"
    )
    bulk_freeze_ready_count = sum(
        1 for row in checksum_rows if _truthy(row.get("freeze_ready", ""))
    )
    bulk_baseline_row_count = len(baseline_rows)
    bulk_evaluated_baseline_count = sum(
        1 for row in baseline_rows if row.get("status") == "evaluated"
    )
    single_cell_assets = _single_cell_asset_paths(root)
    single_cell_h5ad_count = sum(
        1
        for pattern in ("*.h5ad", "*.loom", "*.mtx")
        for _path in root.glob(f"**/{pattern}")
    )
    single_cell_v9_manifest_count = (
        len(list(sc_manifest_dir.glob("*.json"))) if sc_manifest_dir.exists() else 0
    )

    bulk_score = 0
    bulk_score += 20 if bulk_task_count else 0
    bulk_score += 10 if bulk_fold_count else 0
    bulk_score += 15 if bulk_source_count else 0
    bulk_score += 10 if bulk_api_ok_count == bulk_source_count else 0
    bulk_score += 10 if bulk_checksum_parsed_count == bulk_source_count else 0
    bulk_score += 15 if bulk_evaluated_baseline_count == bulk_baseline_row_count else 0
    bulk_score += 10 if datapackage_path.exists() else 0
    bulk_score += 10 if dataset_card_path.exists() else 0
    bulk_score -= 10 if bulk_freeze_ready_count == 0 else 0

    single_cell_score = 0
    single_cell_score += 20 if single_cell_assets else 0
    single_cell_score += 10 if "genelab_sc" in profile_path.read_text() else 0
    single_cell_score += 10 if single_cell_h5ad_count else 0
    single_cell_score += 10 if single_cell_v9_manifest_count else 0
    single_cell_score += 5 if _resolve_path("docs/CELL_EVAL_INTEGRATION_NOTE.md", root).exists() else 0

    selected_next_block = "V9-BULK-ALPHA-001: public bulk alpha freeze-gap matrix"
    selected_next_goal = (
        "Turn the existing public bulk scaffold into a release-readiness gap "
        "matrix before any single-cell flagship implementation."
    )
    claim_boundary = "recenter_decision_only_no_new_benchmark_or_release_claim"
    summary = {
        "recenter_id": recenter_id,
        "decision_status": "selected_public_bulk_alpha_recenter",
        "selected_next_lane": "public_bulk_alpha",
        "selected_next_block": selected_next_block,
        "selected_next_block_goal": selected_next_goal,
        "public_bulk_readiness_score": str(bulk_score),
        "single_cell_readiness_score": str(single_cell_score),
        "bulk_task_count": str(bulk_task_count),
        "bulk_fold_count": str(bulk_fold_count),
        "bulk_source_count": str(bulk_source_count),
        "bulk_api_ok_source_count": str(bulk_api_ok_count),
        "bulk_checksum_parsed_source_count": str(bulk_checksum_parsed_count),
        "bulk_freeze_ready_source_count": str(bulk_freeze_ready_count),
        "bulk_baseline_row_count": str(bulk_baseline_row_count),
        "single_cell_asset_count": str(len(single_cell_assets)),
        "single_cell_h5ad_count": str(single_cell_h5ad_count),
        "single_cell_v9_manifest_count": str(single_cell_v9_manifest_count),
        "osd120_branch_status": osd120_status,
        "claim_boundary": claim_boundary,
    }

    bulk_evidence = [
        f"{bulk_task_count} public bulk task manifests indexed",
        f"{bulk_fold_count} fold definitions across public bulk tasks",
        f"{bulk_source_count} public OSDR source rows",
        f"{bulk_api_ok_count}/{bulk_source_count} source rows API-ok",
        f"{bulk_checksum_parsed_count}/{bulk_source_count} checksum-manifest parsed",
        f"{bulk_evaluated_baseline_count} evaluated baseline rows",
        "draft Data Package present" if datapackage_path.exists() else "draft Data Package missing",
        "dataset card present" if dataset_card_path.exists() else "dataset card missing",
    ]
    bulk_gaps = [
        f"{bulk_freeze_ready_count}/{bulk_source_count} sources freeze_ready=true",
        "payload files are indexed but not locally hash-verified for release",
        "public bulk release claim language still needs a final alpha boundary",
    ]
    sc_evidence = [
        f"{len(single_cell_assets)} legacy RRRM/single-cell-related files found",
        "genelab_sc metric profile exists",
        "cell-eval/OpenProblems alignment notes exist",
    ]
    sc_gaps = [
        f"{single_cell_h5ad_count} local h5ad/loom/mtx files found in repo scan",
        f"{single_cell_v9_manifest_count} v9 sc_spaceflight task manifests found",
        "RRRM assets need inventory before AnnData task-card promotion",
        "single-cell evaluator and loader are not yet implemented",
    ]
    candidates = [
        {
            "recenter_id": recenter_id,
            "candidate_id": "public_bulk_alpha",
            "lane_label": "Public bulk alpha readiness",
            "decision": "selected_next",
            "readiness_status": "platform_spine_ready_freeze_gap_incomplete",
            "readiness_score": str(bulk_score),
            "evidence_items": _pipe_join(bulk_evidence),
            "blocking_gaps": _pipe_join(bulk_gaps),
            "external_method_anchor": (
                "OSDR/API checksum evidence and Data Package-style metadata are "
                "already represented locally"
            ),
            "next_block_if_selected": selected_next_block,
            "rationale": (
                "This lane is closer to a coherent v9-alpha platform checkpoint "
                "because task manifests, source rows, baselines, package draft, "
                "and dataset-card draft already exist."
            ),
        },
        {
            "recenter_id": recenter_id,
            "candidate_id": "single_cell_flagship",
            "lane_label": "First single-cell flagship inventory",
            "decision": "defer_until_after_bulk_alpha_gap_matrix",
            "readiness_status": "legacy_assets_present_no_v9_manifest_or_h5ad_inventory",
            "readiness_score": str(single_cell_score),
            "evidence_items": _pipe_join(sc_evidence),
            "blocking_gaps": _pipe_join(sc_gaps),
            "external_method_anchor": (
                "AnnData/cell-eval/OpenProblems imply that a v9 single-cell lane "
                "needs explicit AnnData, obs/var metadata, and task API contracts"
            ),
            "next_block_if_selected": "V9-SC-001: RRRM asset inventory",
            "rationale": (
                "The lane remains strategically important, but current local "
                "evidence is inventory-level rather than v9 task-manifest-ready."
            ),
        },
    ]
    actions = [
        {
            "recenter_id": recenter_id,
            "block_id": "V9-BULK-ALPHA-001",
            "action_id": "bulk_alpha_gap_matrix",
            "lane_id": "public_bulk_alpha",
            "objective": "Create a machine-readable pass/blocker matrix for the public bulk alpha scaffold.",
            "expected_files": "v9/reports/public_bulk_alpha_gap_matrix/",
            "expected_tests": "tests for row counts, blocker classes, and claim boundary",
            "done_when": "Every current public bulk release blocker has an owner, evidence file, and next action.",
        },
        {
            "recenter_id": recenter_id,
            "block_id": "V9-BULK-ALPHA-001",
            "action_id": "payload_hash_boundary",
            "lane_id": "public_bulk_alpha",
            "objective": "Separate checksum-manifest evidence from locally verified payload hash evidence.",
            "expected_files": "v9/reports/public_bulk_alpha_gap_matrix/payload_hash_boundary.csv",
            "expected_tests": "tests confirm freeze_ready=false is not promoted to release_ready",
            "done_when": "Dataset-card language can explain what is verified now and what remains pending.",
        },
        {
            "recenter_id": recenter_id,
            "block_id": "V9-BULK-ALPHA-001",
            "action_id": "minimal_alpha_snapshot_decision",
            "lane_id": "public_bulk_alpha",
            "objective": "Decide whether a metadata-only alpha snapshot is acceptable before full payload mirroring.",
            "expected_files": "docs/V9_PUBLIC_BULK_ALPHA_SNAPSHOT_DECISION.md",
            "expected_tests": "tests cover allowed and disallowed alpha claims",
            "done_when": "The next user-facing phrase is bounded and not a frozen data release claim.",
        },
        {
            "recenter_id": recenter_id,
            "block_id": "V9-BULK-ALPHA-001",
            "action_id": "bulk_package_update_plan",
            "lane_id": "public_bulk_alpha",
            "objective": "Name the package/card fields that need updates after the gap matrix.",
            "expected_files": "v9/reports/public_bulk_alpha_gap_matrix/package_update_plan.csv",
            "expected_tests": "tests assert no OSD-120 or organoid draft surfaces are pulled into public bulk core",
            "done_when": "Public bulk core, organoid draft, and multispecies draft tracks remain separated.",
        },
    ]
    review_md = f"""# V9-REFOCUS-001 Post-OSD-120 Recenter Decision

Status: `{summary["decision_status"]}`

Selected next lane: `{summary["selected_next_lane"]}`

Selected next block: `{selected_next_block}`

## Decision

Choose public bulk alpha readiness before the first single-cell flagship
implementation.

This is not a retreat from single-cell work. It is a sequencing correction:
the public bulk scaffold already has task manifests, source inventory,
checksum-manifest evidence, baseline outputs, a draft Data Package, and a
dataset-card draft. It still lacks a clean alpha gap matrix and payload hash
boundary. The single-cell lane has valuable RRRM legacy assets and a
`genelab_sc` metric profile, but no v9 `sc_spaceflight` task manifests and no
local h5ad/loom/mtx files found by the current repo scan.

## Readiness Snapshot

| Lane | Score | Decision | Main gap |
|---|---:|---|---|
| Public bulk alpha | {bulk_score} | selected next | payload hash/release claim boundary |
| Single-cell flagship | {single_cell_score} | deferred | RRRM asset inventory before AnnData task cards |

## Selected Next Block

`{selected_next_block}` should produce a machine-readable gap matrix for the
public bulk alpha scaffold, with explicit pass/blocker rows for source
inventory, checksum evidence, payload hash verification, package metadata,
dataset-card language, and baseline output evidence.

## External Method Anchors

- AnnData stores expression matrix `X` with observation metadata in `obs`,
  variable metadata in `var`, and h5ad-backed on-disk storage. A single-cell
  flagship should therefore start from explicit AnnData/obs/var inventory.
- cell-eval expects predicted and real AnnData inputs for perturbation-response
  evaluation, so v9 should not add model adapters before the real/predicted
  AnnData contract is defined.
- OpenProblems emphasizes task APIs and benchmark contracts; v9 should follow
  that discipline for `sc_spaceflight` after the asset inventory.

## Claim Boundary

This is a recenter decision only. It does not create a new benchmark release,
does not promote OSD-120 to the v9-alpha core, and does not claim a
single-cell flagship result.
"""
    return {
        "summary": [summary],
        "candidates": candidates,
        "actions": actions,
        "review_md": review_md,
        "single_cell_asset_paths": single_cell_assets,
    }


def write_v9_recenter_decision(
    *,
    output_dir: str | Path = "v9/reports/recenter_decision",
    repo_root: str | Path = ".",
    recenter_id: str = DEFAULT_V9_RECENTER_DECISION_ID,
) -> dict[str, Path]:
    """Write the V9-REFOCUS-001 decision artifacts."""

    root = Path(repo_root)
    package = build_v9_recenter_decision(
        repo_root=root,
        recenter_id=recenter_id,
    )
    outdir = _resolve_path(output_dir, root)
    outputs = {
        "summary_csv": outdir / "recenter_decision_summary.csv",
        "summary_json": outdir / "recenter_decision_summary.json",
        "candidate_csv": outdir / "recenter_candidate_matrix.csv",
        "candidate_json": outdir / "recenter_candidate_matrix.json",
        "action_csv": outdir / "recenter_next_block_actions.csv",
        "action_json": outdir / "recenter_next_block_actions.json",
        "asset_json": outdir / "single_cell_asset_paths.json",
        "review_md": root / "docs" / "V9_REFOCUS_001_POST_OSD120_RECENTER_DECISION.md",
    }
    _write_csv(
        outputs["summary_csv"],
        package["summary"],
        RECENTER_DECISION_SUMMARY_FIELDS,
    )
    _write_json(outputs["summary_json"], package["summary"])
    _write_csv(
        outputs["candidate_csv"],
        package["candidates"],
        RECENTER_CANDIDATE_FIELDS,
    )
    _write_json(outputs["candidate_json"], package["candidates"])
    _write_csv(
        outputs["action_csv"],
        package["actions"],
        RECENTER_NEXT_BLOCK_ACTION_FIELDS,
    )
    _write_json(outputs["action_json"], package["actions"])
    _write_json(outputs["asset_json"], package["single_cell_asset_paths"])
    outputs["review_md"].parent.mkdir(parents=True, exist_ok=True)
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs

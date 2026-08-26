"""External single-cell payload availability helpers for SpaceBio-Bench v9."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


DEFAULT_SC_EXTERNAL_PAYLOAD_AVAILABILITY_ID = (
    "v9_sc_006b_rrrm1_blood_external_payload_availability"
)

DEFAULT_SC_EXTERNAL_CHECKED_BASES = [
    "/home/fs01/jak4013/rrrm1_scrna",
    "/home/fs01/jak4013/scratch/rrrm1_scrna",
]

SC_EXTERNAL_PAYLOAD_AVAILABILITY_SUMMARY_FIELDS = [
    "availability_id",
    "decision_status",
    "task_id",
    "source_id",
    "tissue",
    "checked_base_count",
    "expected_blood_srx_count",
    "checked_matrix_row_count",
    "annotated_h5ad_found",
    "labeled_h5ad_found",
    "complete_starsolo_srx_count",
    "canonical_copy_allowed",
    "regeneration_allowed",
    "next_required_block",
    "claim_boundary",
]

SC_EXTERNAL_PAYLOAD_CANDIDATE_FIELDS = [
    "availability_id",
    "candidate_id",
    "candidate_type",
    "expected_path",
    "checked_scope",
    "path_status",
    "required_evidence_status",
    "action_decision",
    "notes",
]

SC_EXTERNAL_STARSOLO_MATRIX_FIELDS = [
    "availability_id",
    "source_id",
    "tissue",
    "srx",
    "gsm",
    "animal_id",
    "condition",
    "age_months",
    "checked_base",
    "expected_filtered_dir",
    "filtered_dir_status",
    "matrix_status",
    "features_status",
    "barcodes_status",
    "complete_matrix_bundle",
    "notes",
]

SC_CANONICAL_PAYLOAD_COPY_DECISION_FIELDS = [
    "availability_id",
    "decision_id",
    "source_payload_status",
    "canonical_payload_path",
    "copy_allowed",
    "reason",
    "required_before_copy",
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


def _path_status(path: Path) -> str:
    return "local_exists" if path.exists() else "not_found_in_checked_scope"


def _file_status(paths: list[Path]) -> str:
    return "present" if any(path.exists() for path in paths) else "missing"


def _blood_condition_rows(condition_map_path: Path) -> list[dict[str, str]]:
    rows = [
        row
        for row in _read_csv_rows(condition_map_path)
        if row["osd"] == "OSD-918" and row["tissue"] == "blood"
    ]
    if not rows:
        raise ValueError(f"no OSD-918 blood rows in {condition_map_path}")
    return rows


def _candidate_rows(
    *,
    availability_id: str,
    checked_bases: list[str],
    canonical_payload_path: str,
    repo_root: Path,
    external_probe_note: str,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = [
        {
            "availability_id": availability_id,
            "candidate_id": "canonical_v9_payload_target",
            "candidate_type": "canonical_target",
            "expected_path": canonical_payload_path,
            "checked_scope": "repo",
            "path_status": _path_status(_resolve_path(canonical_payload_path, repo_root)),
            "required_evidence_status": "missing_local_payload",
            "action_decision": "do_not_audit_until_source_payload_exists",
            "notes": "Canonical v9 target remains the only allowed local payload path.",
        }
    ]
    for index, base in enumerate(checked_bases, start=1):
        base_path = Path(base)
        annotated = (
            base_path
            / "downstream_initial"
            / "annotations"
            / "objects"
            / "RRRM1_blood_annotated.h5ad"
        )
        labeled = base_path / "OSD-918" / "OSD-918_blood_labeled.h5ad"
        for candidate_id, candidate_type, expected_path, action in [
            (
                f"legacy_annotated_h5ad_checked_base_{index}",
                "preferred_annotated_h5ad",
                annotated,
                "copy_to_canonical_if_found_and_hashable",
            ),
            (
                f"legacy_labeled_h5ad_checked_base_{index}",
                "metadata_labeled_h5ad",
                labeled,
                "use_only_if_annotation_regeneration_is_possible",
            ),
        ]:
            status = _path_status(expected_path)
            rows.append(
                {
                    "availability_id": availability_id,
                    "candidate_id": candidate_id,
                    "candidate_type": candidate_type,
                    "expected_path": expected_path.as_posix(),
                    "checked_scope": base,
                    "path_status": status,
                    "required_evidence_status": (
                        "source_payload_available" if status == "local_exists" else "source_payload_not_found"
                    ),
                    "action_decision": action if status == "local_exists" else "blocked_not_found",
                    "notes": external_probe_note,
                }
            )
    return rows


def _starsolo_rows(
    *,
    availability_id: str,
    checked_bases: list[str],
    condition_rows: list[dict[str, str]],
    external_probe_note: str,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for condition in condition_rows:
        for base in checked_bases:
            filtered_dir = (
                Path(base)
                / condition["osd"]
                / "starsolo_per_srx"
                / condition["srx"]
                / "Solo.out"
                / "GeneFull"
                / "filtered"
            )
            matrix_status = _file_status(
                [filtered_dir / "matrix.mtx", filtered_dir / "matrix.mtx.gz"]
            )
            features_status = _file_status(
                [filtered_dir / "features.tsv", filtered_dir / "features.tsv.gz"]
            )
            barcodes_status = _file_status(
                [filtered_dir / "barcodes.tsv", filtered_dir / "barcodes.tsv.gz"]
            )
            complete = (
                _path_status(filtered_dir) == "local_exists"
                and matrix_status == "present"
                and features_status == "present"
                and barcodes_status == "present"
            )
            rows.append(
                {
                    "availability_id": availability_id,
                    "source_id": condition["osd"],
                    "tissue": condition["tissue"],
                    "srx": condition["srx"],
                    "gsm": condition["gsm"],
                    "animal_id": condition["source_name"],
                    "condition": condition["condition"],
                    "age_months": condition["age_months"],
                    "checked_base": base,
                    "expected_filtered_dir": filtered_dir.as_posix(),
                    "filtered_dir_status": _path_status(filtered_dir),
                    "matrix_status": matrix_status,
                    "features_status": features_status,
                    "barcodes_status": barcodes_status,
                    "complete_matrix_bundle": str(complete).lower(),
                    "notes": external_probe_note,
                }
            )
    return rows


def build_sc_external_payload_availability(
    *,
    repo_root: str | Path = ".",
    availability_id: str = DEFAULT_SC_EXTERNAL_PAYLOAD_AVAILABILITY_ID,
    manifest_path: str | Path = (
        "v9/sc_spaceflight/task_manifests/"
        "draft_rrrm1_blood_single_cell_spaceflight.json"
    ),
    condition_map_path: str | Path = "v2/docs/RRRM1_SRX_CONDITION_MAP.csv",
    canonical_payload_path: str = (
        "v9/sc_spaceflight/payloads/rrrm1_blood/"
        "OSD-918_blood_rrrm1_bench.h5ad"
    ),
    checked_bases: list[str] | None = None,
    external_probe_note: str = (
        "Exact checked paths are absent in the current repo context; no payload "
        "copy or audit is claimed."
    ),
) -> dict[str, Any]:
    """Build the V9-SC-006b external payload availability package."""

    root = Path(repo_root).resolve()
    bases = list(checked_bases or DEFAULT_SC_EXTERNAL_CHECKED_BASES)
    manifest = json.loads(_resolve_path(manifest_path, root).read_text())
    condition_rows = _blood_condition_rows(_resolve_path(condition_map_path, root))
    candidates = _candidate_rows(
        availability_id=availability_id,
        checked_bases=bases,
        canonical_payload_path=canonical_payload_path,
        repo_root=root,
        external_probe_note=external_probe_note,
    )
    matrix_rows = _starsolo_rows(
        availability_id=availability_id,
        checked_bases=bases,
        condition_rows=condition_rows,
        external_probe_note=external_probe_note,
    )
    annotated_found = any(
        row["candidate_type"] == "preferred_annotated_h5ad"
        and row["path_status"] == "local_exists"
        for row in candidates
    )
    labeled_found = any(
        row["candidate_type"] == "metadata_labeled_h5ad"
        and row["path_status"] == "local_exists"
        for row in candidates
    )
    complete_srxs = {
        row["srx"]
        for row in matrix_rows
        if row["complete_matrix_bundle"] == "true"
    }
    complete_starsolo_srx_count = len(complete_srxs)
    all_srx_complete = complete_starsolo_srx_count == len(condition_rows)
    canonical_copy_allowed = annotated_found
    regeneration_allowed = all_srx_complete
    if canonical_copy_allowed:
        decision_status = "external_payload_found_ready_for_canonical_copy_decision"
        next_required_block = "V9-SC-006c: canonical h5ad copy and obs/var audit rerun"
    elif regeneration_allowed:
        decision_status = "starsolo_matrices_found_ready_for_regeneration"
        next_required_block = "V9-SC-006c: regenerate RRRM-1 blood h5ad"
    else:
        decision_status = "external_payload_availability_blocked_no_h5ad_or_starsolo_matrices"
        next_required_block = "V9-SC-006c: OSDR processed payload discovery or owner scratch path request"
    claim_boundary = "external_payload_availability_no_payload_or_score_claim"
    source = manifest["source_records"][0]
    summary = {
        "availability_id": availability_id,
        "decision_status": decision_status,
        "task_id": manifest["task_id"],
        "source_id": source["source_id"],
        "tissue": source["biospecimen_type"],
        "checked_base_count": str(len(bases)),
        "expected_blood_srx_count": str(len(condition_rows)),
        "checked_matrix_row_count": str(len(matrix_rows)),
        "annotated_h5ad_found": str(annotated_found).lower(),
        "labeled_h5ad_found": str(labeled_found).lower(),
        "complete_starsolo_srx_count": str(complete_starsolo_srx_count),
        "canonical_copy_allowed": str(canonical_copy_allowed).lower(),
        "regeneration_allowed": str(regeneration_allowed).lower(),
        "next_required_block": next_required_block,
        "claim_boundary": claim_boundary,
    }
    source_payload_status = (
        "annotated_h5ad_found"
        if annotated_found
        else "complete_starsolo_matrices_found"
        if all_srx_complete
        else "no_source_payload_found"
    )
    copy_decision = [
        {
            "availability_id": availability_id,
            "decision_id": "canonical_copy_gate",
            "source_payload_status": source_payload_status,
            "canonical_payload_path": canonical_payload_path,
            "copy_allowed": str(canonical_copy_allowed).lower(),
            "reason": (
                "Only a hashable annotated h5ad can be copied directly to the canonical path."
            ),
            "required_before_copy": (
                "source h5ad exists|sha256 recorded|n_obs/n_vars recorded|obs_var_audit rerun planned"
            ),
        }
    ]
    review_md = f"""# V9-SC-006b External Payload Availability

Status: `{decision_status}`

Task id: `{manifest["task_id"]}`

Claim boundary: `{claim_boundary}`

## Decision

No canonical or external RRRM-1 blood source payload is available in the checked
scope. Direct canonical copy is therefore blocked, and STARsolo regeneration is
also blocked until all eight OSD-918 blood SRX matrix bundles are found.

## Checked Scope

- Checked bases: {len(bases)}
- Expected OSD-918 blood SRX rows: {len(condition_rows)}
- Matrix availability rows: {len(matrix_rows)}
- Annotated h5ad found: `{summary["annotated_h5ad_found"]}`
- Labeled h5ad found: `{summary["labeled_h5ad_found"]}`
- Complete STARsolo SRX bundles: {complete_starsolo_srx_count}/{len(condition_rows)}

## Not Claimed

- No h5ad was copied into the canonical v9 payload path.
- No source payload hash, shape, or obs/var pass is claimed.
- No RRRM-1 h5ad regeneration was run.
- No evaluator, leaderboard result, or legacy RRRM score promotion is claimed.

## Next Block

Run `{next_required_block}`.
"""
    return {
        "summary": [summary],
        "candidates": candidates,
        "matrix_availability": matrix_rows,
        "copy_decision": copy_decision,
        "review_md": review_md,
    }


def write_sc_external_payload_availability(
    *,
    output_dir: str | Path = "v9/sc_spaceflight",
    repo_root: str | Path = ".",
    availability_id: str = DEFAULT_SC_EXTERNAL_PAYLOAD_AVAILABILITY_ID,
    manifest_path: str | Path = (
        "v9/sc_spaceflight/task_manifests/"
        "draft_rrrm1_blood_single_cell_spaceflight.json"
    ),
    condition_map_path: str | Path = "v2/docs/RRRM1_SRX_CONDITION_MAP.csv",
    canonical_payload_path: str = (
        "v9/sc_spaceflight/payloads/rrrm1_blood/"
        "OSD-918_blood_rrrm1_bench.h5ad"
    ),
    checked_bases: list[str] | None = None,
    external_probe_note: str = (
        "Exact checked paths are absent in the current repo context; no payload "
        "copy or audit is claimed."
    ),
) -> dict[str, Path]:
    """Write the V9-SC-006b external payload availability package."""

    root = Path(repo_root).resolve()
    package = build_sc_external_payload_availability(
        repo_root=root,
        availability_id=availability_id,
        manifest_path=manifest_path,
        condition_map_path=condition_map_path,
        canonical_payload_path=canonical_payload_path,
        checked_bases=checked_bases,
        external_probe_note=external_probe_note,
    )
    outdir = _resolve_path(output_dir, root)
    outputs = {
        "summary_csv": outdir / "external_payload_availability_summary.csv",
        "summary_json": outdir / "external_payload_availability_summary.json",
        "candidates_csv": outdir / "external_payload_candidates.csv",
        "candidates_json": outdir / "external_payload_candidates.json",
        "matrix_csv": outdir / "external_starsolo_matrix_availability.csv",
        "matrix_json": outdir / "external_starsolo_matrix_availability.json",
        "copy_decision_csv": outdir / "canonical_payload_copy_decision.csv",
        "copy_decision_json": outdir / "canonical_payload_copy_decision.json",
        "review_md": root / "docs" / "V9_SC_EXTERNAL_PAYLOAD_AVAILABILITY.md",
    }
    _write_csv(
        outputs["summary_csv"],
        package["summary"],
        SC_EXTERNAL_PAYLOAD_AVAILABILITY_SUMMARY_FIELDS,
    )
    _write_json(outputs["summary_json"], package["summary"])
    _write_csv(
        outputs["candidates_csv"],
        package["candidates"],
        SC_EXTERNAL_PAYLOAD_CANDIDATE_FIELDS,
    )
    _write_json(outputs["candidates_json"], package["candidates"])
    _write_csv(
        outputs["matrix_csv"],
        package["matrix_availability"],
        SC_EXTERNAL_STARSOLO_MATRIX_FIELDS,
    )
    _write_json(outputs["matrix_json"], package["matrix_availability"])
    _write_csv(
        outputs["copy_decision_csv"],
        package["copy_decision"],
        SC_CANONICAL_PAYLOAD_COPY_DECISION_FIELDS,
    )
    _write_json(outputs["copy_decision_json"], package["copy_decision"])
    outputs["review_md"].parent.mkdir(parents=True, exist_ok=True)
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs

"""OSDR processed-payload discovery helpers for SpaceBio-Bench v9 single-cell."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Callable

from .source_audit import (
    BIODATA_API_BASE,
    FetchResult,
    extract_files_from_listing,
    fetch_url,
    source_file_api_url,
)


DEFAULT_SC_PROCESSED_PAYLOAD_DISCOVERY_ID = (
    "v9_sc_006c_osd918_processed_payload_discovery"
)

SC_PROCESSED_PAYLOAD_DISCOVERY_SUMMARY_FIELDS = [
    "discovery_id",
    "decision_status",
    "task_id",
    "source_id",
    "glds_prefix",
    "tissue",
    "api_url",
    "api_status",
    "api_response_sha256",
    "osdr_file_count",
    "metadata_file_count",
    "raw_fastq_count",
    "raw_checksum_count",
    "raw_multiqc_count",
    "processed_h5ad_count",
    "processed_starsolo_count",
    "processed_checksum_manifest_count",
    "expected_blood_srx_count",
    "complete_expected_fastq_pair_count",
    "missing_expected_fastq_pair_count",
    "canonical_copy_allowed",
    "regeneration_allowed",
    "owner_scratch_request_required",
    "next_required_block",
    "claim_boundary",
]

SC_OSDR_FILE_DISCOVERY_FIELDS = [
    "discovery_id",
    "source_id",
    "filename",
    "rest_url",
    "download_url",
    "file_role",
    "processed_payload_candidate",
    "expected_srx",
    "claim_status",
    "notes",
]

SC_OSDR_EXPECTED_SRX_COVERAGE_FIELDS = [
    "discovery_id",
    "source_id",
    "tissue",
    "srx",
    "gsm",
    "animal_id",
    "condition",
    "expected_r1_file",
    "expected_r1_status",
    "expected_r2_file",
    "expected_r2_status",
    "fastq_pair_complete",
    "processed_matrix_status",
    "action_decision",
]

SC_OWNER_SCRATCH_REQUEST_FIELDS = [
    "discovery_id",
    "request_id",
    "recipient_role",
    "request_status",
    "required_path_or_artifact",
    "required_evidence",
    "reason",
    "unblock_step",
]

SC_PAYLOAD_DISCOVERY_DEFERRAL_FIELDS = [
    "discovery_id",
    "decision_id",
    "osdr_processed_payload_status",
    "canonical_payload_path",
    "copy_allowed",
    "regeneration_allowed",
    "deferral_reason",
    "next_action",
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


def _blood_condition_rows(condition_map_path: Path) -> list[dict[str, str]]:
    rows = [
        row
        for row in _read_csv_rows(condition_map_path)
        if row["osd"] == "OSD-918" and row["tissue"] == "blood"
    ]
    if not rows:
        raise ValueError(f"no OSD-918 blood rows in {condition_map_path}")
    return rows


def _load_listing_result(
    *,
    source_id: str,
    api_base: str,
    api_listing_json: str | Path | None,
    fetcher: Callable[[str], FetchResult] | None,
) -> tuple[str, FetchResult]:
    api_url = source_file_api_url(source_id, api_base=api_base)
    if api_listing_json is not None:
        body = Path(api_listing_json).read_bytes()
        return api_url, FetchResult(ok=True, url=str(api_listing_json), body=body)
    if fetcher is not None:
        return api_url, fetcher(api_url)
    return api_url, fetch_url(api_url, timeout=30, max_bytes=None)


def _classify_file(filename: str) -> tuple[str, str, str]:
    lower = filename.lower()
    if lower.endswith(".h5ad") or "anndata" in lower:
        return "processed_h5ad", "true", "remote_processed_payload_candidate"
    if any(token in lower for token in ["starsolo", "solo.out", "matrix.mtx"]):
        return "processed_starsolo_matrix", "true", "remote_processed_payload_candidate"
    if "processed" in lower and any(
        token in lower for token in ["md5", "checksum", "sha256", "sha1"]
    ):
        return "processed_checksum_manifest", "false", "checksum_evidence_only"
    if lower.endswith("_raw.fastq.gz") or "_raw.fastq.gz" in lower:
        return "raw_fastq", "false", "raw_input_only_not_processed_payload"
    if "raw_md5sum" in lower or ("raw" in lower and "md5" in lower):
        return "raw_checksum_manifest", "false", "checksum_evidence_only"
    if "multiqc" in lower:
        return "raw_multiqc_report", "false", "qc_report_not_payload"
    if "metadata" in lower or lower.endswith("-isa.zip"):
        return "metadata", "false", "metadata_only"
    return "other", "false", "not_a_processed_payload"


def _expected_srx_from_filename(filename: str, expected_srxs: set[str]) -> str:
    matches = sorted(srx for srx in expected_srxs if srx in filename)
    return "|".join(matches)


def _file_rows(
    *,
    discovery_id: str,
    source_id: str,
    files: dict[str, dict[str, Any]],
    expected_srxs: set[str],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for filename, metadata in sorted(files.items()):
        role, candidate, claim_status = _classify_file(filename)
        rows.append(
            {
                "discovery_id": discovery_id,
                "source_id": source_id,
                "filename": filename,
                "rest_url": str(metadata.get("REST_URL", "")),
                "download_url": str(metadata.get("URL", "")),
                "file_role": role,
                "processed_payload_candidate": candidate,
                "expected_srx": _expected_srx_from_filename(filename, expected_srxs),
                "claim_status": claim_status,
                "notes": "OSDR file-list discovery row; no payload was downloaded.",
            }
        )
    return rows


def _srx_coverage_rows(
    *,
    discovery_id: str,
    source_id: str,
    condition_rows: list[dict[str, str]],
    file_names: set[str],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for condition in condition_rows:
        srx = condition["srx"]
        r1 = next(
            (name for name in sorted(file_names) if srx in name and "_R1_raw.fastq.gz" in name),
            "",
        )
        r2 = next(
            (name for name in sorted(file_names) if srx in name and "_R2_raw.fastq.gz" in name),
            "",
        )
        processed_matrix_found = any(
            srx in name
            and any(token in name.lower() for token in ["starsolo", "matrix.mtx", "h5ad"])
            for name in file_names
        )
        complete = bool(r1 and r2)
        rows.append(
            {
                "discovery_id": discovery_id,
                "source_id": source_id,
                "tissue": condition["tissue"],
                "srx": srx,
                "gsm": condition["gsm"],
                "animal_id": condition["source_name"],
                "condition": condition["condition"],
                "expected_r1_file": r1 or f"GLDS-746_scRNA-Seq_{srx}_R1_raw.fastq.gz",
                "expected_r1_status": "listed" if r1 else "missing",
                "expected_r2_file": r2 or f"GLDS-746_scRNA-Seq_{srx}_R2_raw.fastq.gz",
                "expected_r2_status": "listed" if r2 else "missing",
                "fastq_pair_complete": str(complete).lower(),
                "processed_matrix_status": (
                    "listed" if processed_matrix_found else "not_listed_in_osdr_file_api"
                ),
                "action_decision": (
                    "raw_fastq_pair_available_but_regeneration_not_started"
                    if complete
                    else "raw_fastq_pair_missing"
                ),
            }
        )
    return rows


def _owner_request_rows(discovery_id: str) -> list[dict[str, str]]:
    return [
        {
            "discovery_id": discovery_id,
            "request_id": "preferred_annotated_h5ad_owner_path",
            "recipient_role": "dataset_owner_or_local_pipeline_operator",
            "request_status": "required",
            "required_path_or_artifact": (
                "RRRM1_blood_annotated.h5ad or exact scratch path to that object"
            ),
            "required_evidence": "absolute_path|sha256|n_obs|n_vars|annotation_source",
            "reason": "Preferred route to canonical v9 payload without re-running raw FASTQ processing.",
            "unblock_step": "canonical h5ad copy and obs/var audit rerun",
        },
        {
            "discovery_id": discovery_id,
            "request_id": "all_eight_starsolo_matrix_bundles",
            "recipient_role": "dataset_owner_or_local_pipeline_operator",
            "request_status": "required_if_annotated_h5ad_absent",
            "required_path_or_artifact": (
                "OSD-918/starsolo_per_srx/<SRX>/Solo.out/GeneFull/filtered for all 8 blood SRX"
            ),
            "required_evidence": "matrix.mtx(.gz)|features.tsv(.gz)|barcodes.tsv(.gz)|per_srx_paths",
            "reason": "Required to regenerate the labeled h5ad from existing STARsolo output.",
            "unblock_step": "RRRM-1 blood h5ad regeneration",
        },
        {
            "discovery_id": discovery_id,
            "request_id": "raw_fastq_regeneration_feasibility",
            "recipient_role": "hpc_pipeline_operator",
            "request_status": "deferred_planning_only",
            "required_path_or_artifact": "raw FASTQ pairs for 8 OSD-918 blood SRX",
            "required_evidence": "OSDR raw FASTQ listing|raw md5 manifest|HPC storage and compute plan",
            "reason": "OSDR lists raw inputs, but regeneration from raw FASTQ is a separate HPC block.",
            "unblock_step": "raw FASTQ to STARsolo regeneration feasibility decision",
        },
    ]


def build_sc_processed_payload_discovery(
    *,
    repo_root: str | Path = ".",
    discovery_id: str = DEFAULT_SC_PROCESSED_PAYLOAD_DISCOVERY_ID,
    manifest_path: str | Path = (
        "v9/sc_spaceflight/task_manifests/"
        "draft_rrrm1_blood_single_cell_spaceflight.json"
    ),
    condition_map_path: str | Path = "v2/docs/RRRM1_SRX_CONDITION_MAP.csv",
    canonical_payload_path: str = (
        "v9/sc_spaceflight/payloads/rrrm1_blood/"
        "OSD-918_blood_rrrm1_bench.h5ad"
    ),
    api_base: str = BIODATA_API_BASE,
    api_listing_json: str | Path | None = None,
    fetcher: Callable[[str], FetchResult] | None = None,
) -> dict[str, Any]:
    """Build the V9-SC-006c OSDR processed-payload discovery package."""

    root = Path(repo_root).resolve()
    manifest = json.loads(_resolve_path(manifest_path, root).read_text())
    source = manifest["source_records"][0]
    source_id = source["source_id"]
    condition_rows = _blood_condition_rows(_resolve_path(condition_map_path, root))
    api_url, listing_result = _load_listing_result(
        source_id=source_id,
        api_base=api_base,
        api_listing_json=api_listing_json,
        fetcher=fetcher,
    )
    expected_srxs = {row["srx"] for row in condition_rows}
    files: dict[str, dict[str, Any]] = {}
    api_status = "ok" if listing_result.ok else "error"
    parse_error = ""
    if listing_result.ok:
        try:
            listing = json.loads(listing_result.body.decode("utf-8"))
            files = {
                filename: dict(metadata)
                for filename, metadata in extract_files_from_listing(source_id, listing).items()
            }
        except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
            api_status = "invalid_json"
            parse_error = str(exc)
    else:
        parse_error = listing_result.error

    file_rows = _file_rows(
        discovery_id=discovery_id,
        source_id=source_id,
        files=files,
        expected_srxs=expected_srxs,
    )
    coverage_rows = _srx_coverage_rows(
        discovery_id=discovery_id,
        source_id=source_id,
        condition_rows=condition_rows,
        file_names=set(files),
    )
    owner_requests = _owner_request_rows(discovery_id)
    role_counts = {
        role: sum(1 for row in file_rows if row["file_role"] == role)
        for role in {
            "metadata",
            "raw_fastq",
            "raw_checksum_manifest",
            "raw_multiqc_report",
            "processed_h5ad",
            "processed_starsolo_matrix",
            "processed_checksum_manifest",
        }
    }
    complete_fastq_pairs = sum(
        1 for row in coverage_rows if row["fastq_pair_complete"] == "true"
    )
    missing_fastq_pairs = len(coverage_rows) - complete_fastq_pairs
    processed_h5ad_count = role_counts["processed_h5ad"]
    processed_starsolo_count = role_counts["processed_starsolo_matrix"]
    processed_checksum_count = role_counts["processed_checksum_manifest"]

    if api_status != "ok":
        decision_status = "osdr_processed_payload_discovery_api_error"
        owner_scratch_request_required = "true"
        next_required_block = "V9-SC-006d: retry OSDR discovery or owner scratch path intake"
    elif processed_h5ad_count:
        decision_status = "osdr_processed_h5ad_listed_download_hash_required"
        owner_scratch_request_required = "false"
        next_required_block = "V9-SC-006d: download h5ad hash shape and obs/var audit"
    elif processed_starsolo_count:
        decision_status = "osdr_processed_starsolo_listed_download_bundle_audit_required"
        owner_scratch_request_required = "false"
        next_required_block = "V9-SC-006d: download STARsolo bundles and regenerate h5ad"
    else:
        decision_status = "osdr_processed_payload_unavailable_owner_scratch_request_required"
        owner_scratch_request_required = "true"
        next_required_block = (
            "V9-SC-006d: owner scratch path intake or raw FASTQ regeneration feasibility decision"
        )

    canonical_copy_allowed = "false"
    regeneration_allowed = "false"
    claim_boundary = "osdr_processed_payload_discovery_only_no_payload_copy_or_score_claim"
    deferral_reason = (
        parse_error
        if parse_error
        else (
            "OSDR lists raw FASTQ inputs and raw QC/checksum files, but no processed h5ad, "
            "processed STARsolo matrix bundle, or processed checksum manifest for the v9 "
            "RRRM-1 blood payload."
        )
    )
    summary = {
        "discovery_id": discovery_id,
        "decision_status": decision_status,
        "task_id": manifest["task_id"],
        "source_id": source_id,
        "glds_prefix": "GLDS-746",
        "tissue": source["biospecimen_type"],
        "api_url": api_url,
        "api_status": api_status,
        "api_response_sha256": listing_result.sha256,
        "osdr_file_count": str(len(file_rows)),
        "metadata_file_count": str(role_counts["metadata"]),
        "raw_fastq_count": str(role_counts["raw_fastq"]),
        "raw_checksum_count": str(role_counts["raw_checksum_manifest"]),
        "raw_multiqc_count": str(role_counts["raw_multiqc_report"]),
        "processed_h5ad_count": str(processed_h5ad_count),
        "processed_starsolo_count": str(processed_starsolo_count),
        "processed_checksum_manifest_count": str(processed_checksum_count),
        "expected_blood_srx_count": str(len(condition_rows)),
        "complete_expected_fastq_pair_count": str(complete_fastq_pairs),
        "missing_expected_fastq_pair_count": str(missing_fastq_pairs),
        "canonical_copy_allowed": canonical_copy_allowed,
        "regeneration_allowed": regeneration_allowed,
        "owner_scratch_request_required": owner_scratch_request_required,
        "next_required_block": next_required_block,
        "claim_boundary": claim_boundary,
    }
    deferral = [
        {
            "discovery_id": discovery_id,
            "decision_id": "canonical_copy_or_regeneration_deferral",
            "osdr_processed_payload_status": decision_status,
            "canonical_payload_path": canonical_payload_path,
            "copy_allowed": canonical_copy_allowed,
            "regeneration_allowed": regeneration_allowed,
            "deferral_reason": deferral_reason,
            "next_action": next_required_block,
        }
    ]
    review_md = f"""# V9-SC-006c OSDR Processed Payload Discovery

Status: `{decision_status}`

Task id: `{manifest["task_id"]}`

Claim boundary: `{claim_boundary}`

## Decision

The OSDR Biodata file-list endpoint for `{source_id}` was queried at
`{api_url}`. The current file list records {len(file_rows)} files: {role_counts["metadata"]}
metadata file, {role_counts["raw_checksum_manifest"]} raw checksum manifest,
{role_counts["raw_multiqc_report"]} raw MultiQC report, and
{role_counts["raw_fastq"]} raw FASTQ files. It does not list a processed h5ad,
processed STARsolo matrix bundle, or processed checksum manifest for the
canonical RRRM-1 blood AnnData target.

## Coverage

- Expected OSD-918 blood SRX rows: {len(condition_rows)}
- Complete expected raw FASTQ pairs listed by OSDR: {complete_fastq_pairs}/{len(condition_rows)}
- Processed h5ad files listed: {processed_h5ad_count}
- Processed STARsolo matrix rows listed: {processed_starsolo_count}
- Canonical copy allowed: `{canonical_copy_allowed}`
- Regeneration allowed: `{regeneration_allowed}`
- Owner scratch request required: `{owner_scratch_request_required}`

## Not Claimed

- No OSDR payload was downloaded.
- No source payload hash, AnnData shape, or obs/var pass is claimed.
- No canonical v9 h5ad copy was made.
- No raw FASTQ-to-STARsolo regeneration was started.
- No evaluator, leaderboard result, or legacy RRRM score promotion is claimed.

## Next Block

Run `{next_required_block}`.
"""
    return {
        "summary": [summary],
        "files": file_rows,
        "expected_srx_coverage": coverage_rows,
        "owner_requests": owner_requests,
        "deferral": deferral,
        "review_md": review_md,
    }


def write_sc_processed_payload_discovery(
    *,
    output_dir: str | Path = "v9/sc_spaceflight",
    repo_root: str | Path = ".",
    discovery_id: str = DEFAULT_SC_PROCESSED_PAYLOAD_DISCOVERY_ID,
    manifest_path: str | Path = (
        "v9/sc_spaceflight/task_manifests/"
        "draft_rrrm1_blood_single_cell_spaceflight.json"
    ),
    condition_map_path: str | Path = "v2/docs/RRRM1_SRX_CONDITION_MAP.csv",
    canonical_payload_path: str = (
        "v9/sc_spaceflight/payloads/rrrm1_blood/"
        "OSD-918_blood_rrrm1_bench.h5ad"
    ),
    api_base: str = BIODATA_API_BASE,
    api_listing_json: str | Path | None = None,
    fetcher: Callable[[str], FetchResult] | None = None,
) -> dict[str, Path]:
    """Write the V9-SC-006c OSDR processed-payload discovery package."""

    root = Path(repo_root).resolve()
    package = build_sc_processed_payload_discovery(
        repo_root=root,
        discovery_id=discovery_id,
        manifest_path=manifest_path,
        condition_map_path=condition_map_path,
        canonical_payload_path=canonical_payload_path,
        api_base=api_base,
        api_listing_json=api_listing_json,
        fetcher=fetcher,
    )
    outdir = _resolve_path(output_dir, root)
    outputs = {
        "summary_csv": outdir / "osdr_processed_payload_discovery_summary.csv",
        "summary_json": outdir / "osdr_processed_payload_discovery_summary.json",
        "files_csv": outdir / "osdr_file_discovery.csv",
        "files_json": outdir / "osdr_file_discovery.json",
        "coverage_csv": outdir / "osdr_expected_srx_coverage.csv",
        "coverage_json": outdir / "osdr_expected_srx_coverage.json",
        "owner_request_csv": outdir / "owner_scratch_request.csv",
        "owner_request_json": outdir / "owner_scratch_request.json",
        "deferral_csv": outdir / "processed_payload_deferral_decision.csv",
        "deferral_json": outdir / "processed_payload_deferral_decision.json",
        "review_md": root / "docs" / "V9_SC_OSDR_PROCESSED_PAYLOAD_DISCOVERY.md",
    }
    _write_csv(
        outputs["summary_csv"],
        package["summary"],
        SC_PROCESSED_PAYLOAD_DISCOVERY_SUMMARY_FIELDS,
    )
    _write_json(outputs["summary_json"], package["summary"])
    _write_csv(outputs["files_csv"], package["files"], SC_OSDR_FILE_DISCOVERY_FIELDS)
    _write_json(outputs["files_json"], package["files"])
    _write_csv(
        outputs["coverage_csv"],
        package["expected_srx_coverage"],
        SC_OSDR_EXPECTED_SRX_COVERAGE_FIELDS,
    )
    _write_json(outputs["coverage_json"], package["expected_srx_coverage"])
    _write_csv(
        outputs["owner_request_csv"],
        package["owner_requests"],
        SC_OWNER_SCRATCH_REQUEST_FIELDS,
    )
    _write_json(outputs["owner_request_json"], package["owner_requests"])
    _write_csv(
        outputs["deferral_csv"],
        package["deferral"],
        SC_PAYLOAD_DISCOVERY_DEFERRAL_FIELDS,
    )
    _write_json(outputs["deferral_json"], package["deferral"])
    outputs["review_md"].parent.mkdir(parents=True, exist_ok=True)
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs

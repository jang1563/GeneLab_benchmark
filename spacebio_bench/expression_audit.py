"""Expression-matrix audit helpers for draft extension tasks."""

from __future__ import annotations

import csv
import hashlib
import io
import json
import re
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

from .source_audit import (
    BIODATA_API_BASE,
    FetchResult,
    extract_files_from_listing,
    fetch_url,
    source_file_api_url,
)


NORMALIZED_COUNT_RE = re.compile(r"rna_seq_normalized_counts_glbulkrnaseq\.csv$", re.I)
EXPRESSION_MATRIX_AUDIT_FIELDS = [
    "source_id",
    "glds_prefix",
    "task_ids",
    "api_url",
    "api_status",
    "n_files",
    "matrix_file_count",
    "matrix_files",
    "matrix_urls",
    "downloaded_matrix_count",
    "local_matrix_paths",
    "matrix_sha256",
    "matrix_bytes",
    "n_feature_rows",
    "n_matrix_columns",
    "sample_columns",
    "n_sample_columns",
    "expected_sample_count",
    "matching_sample_count",
    "missing_samples",
    "extra_samples",
    "audit_status",
    "pending_reason",
]


def _join(values: Sequence[str]) -> str:
    return ";".join(sorted({value for value in values if value}))


def _download_url(metadata: Mapping[str, Any]) -> str:
    return str(
        metadata.get("URL")
        or metadata.get("url")
        or metadata.get("download_url")
        or metadata.get("REST_URL")
        or ""
    )


def is_normalized_count_matrix(filename: str) -> bool:
    normalized = filename.replace("-", "_").lower()
    if "unnormalized" in normalized or "rrnarm" in normalized:
        return False
    return bool(NORMALIZED_COUNT_RE.search(normalized))


def _source_expected_samples(
    sample_factor_rows: Sequence[Mapping[str, str]],
    source_id: str,
) -> list[str]:
    return sorted(
        {
            str(row.get("sample_id", "") or "")
            for row in sample_factor_rows
            if row.get("source_id") == source_id and row.get("sample_id")
        }
    )


def inspect_expression_matrix(text: str) -> dict[str, object]:
    """Inspect a normalized-count CSV without assuming a specific gene id name."""

    reader = csv.reader(io.StringIO(text.lstrip("\ufeff")))
    try:
        header = next(reader)
    except StopIteration:
        return {
            "n_feature_rows": 0,
            "n_matrix_columns": 0,
            "sample_columns": [],
        }
    sample_columns = [column.strip() for column in header[1:] if column.strip()]
    n_feature_rows = 0
    for row in reader:
        if any(str(value or "").strip() for value in row):
            n_feature_rows += 1
    return {
        "n_feature_rows": n_feature_rows,
        "n_matrix_columns": len(header),
        "sample_columns": sample_columns,
    }


def audit_expression_matrix_row(
    row: Mapping[str, str],
    *,
    sample_factor_rows: Sequence[Mapping[str, str]],
    output_dir: str | Path,
    fetcher: Callable[[str], FetchResult] | None = None,
    api_base: str = BIODATA_API_BASE,
) -> dict[str, str]:
    source_id = str(row["source_id"])
    glds_prefix = str(row.get("glds_prefix", ""))
    api_url = source_file_api_url(source_id, api_base=api_base)
    fetch = fetcher or (lambda url: fetch_url(url, timeout=60, max_bytes=None))
    listing_result = fetch(api_url)
    output = {field: "" for field in EXPRESSION_MATRIX_AUDIT_FIELDS}
    output.update(
        {
            "source_id": source_id,
            "glds_prefix": glds_prefix,
            "task_ids": str(row.get("task_ids", "")),
            "api_url": api_url,
        }
    )
    if not listing_result.ok:
        output.update(
            {
                "api_status": "error",
                "audit_status": "api_error",
                "pending_reason": listing_result.error,
            }
        )
        return output

    try:
        listing = json.loads(listing_result.body.decode("utf-8"))
        files = extract_files_from_listing(source_id, listing)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        output.update(
            {
                "api_status": "invalid_json",
                "audit_status": "api_error",
                "pending_reason": str(exc),
            }
        )
        return output

    matrix_files = [
        (filename, metadata)
        for filename, metadata in sorted(files.items())
        if is_normalized_count_matrix(filename)
    ]
    output["api_status"] = "ok"
    output["n_files"] = str(len(files))
    output["matrix_file_count"] = str(len(matrix_files))
    output["matrix_files"] = _join([filename for filename, _ in matrix_files])
    output["matrix_urls"] = _join([_download_url(metadata) for _, metadata in matrix_files])
    if not matrix_files:
        output.update(
            {
                "audit_status": "no_normalized_count_matrix_listed",
                "pending_reason": "OSDR file listing returned no normalized-count matrix.",
            }
        )
        return output

    expected_samples = _source_expected_samples(sample_factor_rows, source_id)
    output_path_dir = Path(output_dir)
    output_path_dir.mkdir(parents=True, exist_ok=True)

    downloaded = 0
    local_paths: list[str] = []
    hashes: list[str] = []
    byte_sizes: list[str] = []
    feature_rows: list[str] = []
    matrix_columns: list[str] = []
    all_sample_columns: list[str] = []
    errors: list[str] = []

    for filename, metadata in matrix_files:
        url = _download_url(metadata)
        if not url:
            errors.append(f"{filename}: missing download URL")
            continue
        matrix_result = fetch(url)
        if not matrix_result.ok:
            errors.append(f"{filename}: {matrix_result.error}")
            continue
        downloaded += 1
        file_path = output_path_dir / filename
        file_path.write_bytes(matrix_result.body)
        local_paths.append(file_path.as_posix())
        hashes.append(hashlib.sha256(matrix_result.body).hexdigest())
        byte_sizes.append(str(len(matrix_result.body)))
        inspected = inspect_expression_matrix(matrix_result.body.decode("utf-8", errors="replace"))
        feature_rows.append(str(inspected["n_feature_rows"]))
        matrix_columns.append(str(inspected["n_matrix_columns"]))
        all_sample_columns.extend(str(value) for value in inspected["sample_columns"])

    sample_columns = sorted({value for value in all_sample_columns if value})
    expected_set = set(expected_samples)
    observed_set = set(sample_columns)
    missing = sorted(expected_set - observed_set)
    extra = sorted(observed_set - expected_set)

    output["downloaded_matrix_count"] = str(downloaded)
    output["local_matrix_paths"] = _join(local_paths)
    output["matrix_sha256"] = _join(hashes)
    output["matrix_bytes"] = _join(byte_sizes)
    output["n_feature_rows"] = _join(feature_rows)
    output["n_matrix_columns"] = _join(matrix_columns)
    output["sample_columns"] = _join(sample_columns)
    output["n_sample_columns"] = str(len(sample_columns))
    output["expected_sample_count"] = str(len(expected_samples))
    output["matching_sample_count"] = str(len(expected_set & observed_set))
    output["missing_samples"] = _join(missing)
    output["extra_samples"] = _join(extra)

    if downloaded and not missing and not extra and expected_samples:
        output.update(
            {
                "audit_status": "matrix_downloaded_sample_aligned",
                "pending_reason": "Matrix downloaded and sample columns align with parsed sample factors.",
            }
        )
    elif downloaded:
        output.update(
            {
                "audit_status": "matrix_downloaded_sample_mismatch",
                "pending_reason": "Matrix downloaded, but sample columns do not exactly match parsed sample factors.",
            }
        )
    else:
        output.update(
            {
                "audit_status": "matrix_download_failed",
                "pending_reason": " | ".join(errors) or "No matrix files were downloaded.",
            }
        )
    return output


def audit_expression_matrix_inventory(
    rows: Sequence[Mapping[str, str]],
    *,
    sample_factor_rows: Sequence[Mapping[str, str]],
    output_dir: str | Path,
    fetcher: Callable[[str], FetchResult] | None = None,
    api_base: str = BIODATA_API_BASE,
) -> list[dict[str, str]]:
    return [
        audit_expression_matrix_row(
            row,
            sample_factor_rows=sample_factor_rows,
            output_dir=output_dir,
            fetcher=fetcher,
            api_base=api_base,
        )
        for row in rows
    ]


def write_expression_matrix_audit(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    if not rows:
        raise ValueError("cannot write an empty expression matrix audit")
    output_csv = Path(csv_path)
    output_json = Path(json_path)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    normalized = [
        {field: str(row.get(field, "") or "") for field in EXPRESSION_MATRIX_AUDIT_FIELDS}
        for row in rows
    ]
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=EXPRESSION_MATRIX_AUDIT_FIELDS)
        writer.writeheader()
        writer.writerows(normalized)
    output_json.write_text(json.dumps(normalized, indent=2, sort_keys=True) + "\n")
    return output_csv, output_json

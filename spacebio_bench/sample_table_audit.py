"""Sample-table audit helpers for draft extension sources."""

from __future__ import annotations

import csv
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


SAMPLE_TABLE_NAME_RE = re.compile(r"(sample.?table|sample.?metadata|metadata)", re.I)
CONDITION_COLUMN_RE = re.compile(
    r"(factor|condition|spaceflight|microgravity|ground|microglia|organoid|"
    r"donor|subject|genotype|treatment|group|sample.?source|characteristics)",
    re.I,
)
SAMPLE_TABLE_AUDIT_FIELDS = [
    "source_id",
    "glds_prefix",
    "task_ids",
    "api_url",
    "api_status",
    "api_response_sha256",
    "n_files",
    "sample_table_file_count",
    "sample_table_files",
    "sample_table_urls",
    "fetched_sample_table_count",
    "sample_rows",
    "sample_table_columns",
    "candidate_condition_columns",
    "condition_value_summary",
    "audit_status",
    "pending_reason",
]
SAMPLE_FACTOR_FIELDS = [
    "source_id",
    "glds_prefix",
    "sample_id",
    "raw_condition",
    "disease_context",
    "spaceflight_condition",
    "true_label",
    "microglia_condition",
    "organoid_type",
    "sample_table_file",
    "parse_status",
    "donor_or_line_id",
    "ipsc_line_id",
    "geo_sample_title",
    "geo_cell_type",
    "geo_treatment",
    "geo_metadata_status",
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


def is_sample_table_file(filename: str) -> bool:
    lower = filename.lower()
    if not SAMPLE_TABLE_NAME_RE.search(lower):
        return False
    return lower.endswith((".csv", ".tsv", ".txt"))


def _delimiter_for(filename: str, text: str) -> str:
    if filename.lower().endswith(".tsv"):
        return "\t"
    first_line = text.splitlines()[0] if text.splitlines() else ""
    return "\t" if first_line.count("\t") > first_line.count(",") else ","


def parse_sample_table(text: str, *, filename: str = "") -> tuple[list[str], list[dict[str, str]]]:
    """Parse a CSV/TSV sample table into column names and rows."""

    stripped = text.lstrip("\ufeff")
    if not stripped.strip():
        return [], []
    delimiter = _delimiter_for(filename, stripped)
    reader = csv.DictReader(io.StringIO(stripped), delimiter=delimiter)
    columns = [str(column or "").strip() for column in (reader.fieldnames or [])]
    rows = [
        {str(key or "").strip(): str(value or "").strip() for key, value in row.items()}
        for row in reader
        if any(str(value or "").strip() for value in row.values())
    ]
    return columns, rows


def _normalize_token(value: str) -> str:
    normalized = re.sub(r"[^A-Za-z0-9]+", "_", value.strip().lower()).strip("_")
    return normalized or "unknown"


def parse_condition_factors(condition: str) -> dict[str, str]:
    """Parse OSD-863/OSD-871 compact condition strings into factors."""

    raw = str(condition or "").strip()
    output = {
        "raw_condition": raw,
        "disease_context": "",
        "spaceflight_condition": "",
        "true_label": "",
        "microglia_condition": "",
        "parse_status": "unparsed",
    }
    parts = raw.split("...")
    if len(parts) != 3:
        output["parse_status"] = "unexpected_condition_format"
        return output

    disease_raw, spaceflight_raw, microglia_raw = parts
    disease_lookup = {
        "no.known.diseases": "no_known_diseases",
        "primary.progressive.multiple.sclerosis": "primary_progressive_multiple_sclerosis",
        "sporadic.parkinson.disease": "sporadic_parkinson_disease",
    }
    spaceflight_lookup = {
        "ground.control": ("Ground", "Ground"),
        "space.flight": ("LEO_or_ISS", "LEO_or_ISS"),
    }
    microglia_lookup = {
        "with.microglia": "with_microglia",
        "without.microglia": "without_microglia",
    }

    disease_key = disease_raw.lower()
    spaceflight_key = spaceflight_raw.lower()
    microglia_key = microglia_raw.lower()
    output["disease_context"] = disease_lookup.get(disease_key, _normalize_token(disease_raw))
    if spaceflight_key in spaceflight_lookup:
        output["spaceflight_condition"], output["true_label"] = spaceflight_lookup[spaceflight_key]
    else:
        output["spaceflight_condition"] = _normalize_token(spaceflight_raw)
        output["true_label"] = output["spaceflight_condition"]
    output["microglia_condition"] = microglia_lookup.get(
        microglia_key,
        _normalize_token(microglia_raw),
    )
    output["parse_status"] = "parsed"
    return output


def _sample_id_from_row(row: Mapping[str, str]) -> str:
    for key in ("", "Sample Name", "sample_name", "sample_id", "Sample ID", "Run"):
        value = str(row.get(key, "") or "").strip()
        if value:
            return value
    for value in row.values():
        text = str(value or "").strip()
        if text:
            return text
    return ""


def _condition_columns(columns: Sequence[str]) -> list[str]:
    return [column for column in columns if CONDITION_COLUMN_RE.search(column)]


def _condition_summary(rows: Sequence[Mapping[str, str]], columns: Sequence[str]) -> str:
    parts: list[str] = []
    for column in columns:
        values = sorted({str(row.get(column, "") or "") for row in rows if row.get(column)})
        if not values:
            continue
        shown = values[:8]
        suffix = f"+{len(values) - len(shown)}more" if len(values) > len(shown) else ""
        parts.append(f"{column}={'/'.join(shown)}{suffix}")
    return " | ".join(parts)


def audit_sample_table_row(
    row: Mapping[str, str],
    *,
    fetcher: Callable[[str], FetchResult] | None = None,
    api_base: str = BIODATA_API_BASE,
) -> dict[str, str]:
    source_id = str(row["source_id"])
    api_url = source_file_api_url(source_id, api_base=api_base)
    fetch = fetcher or (lambda url: fetch_url(url, timeout=30, max_bytes=None))
    listing_result = fetch(api_url)
    output = {field: "" for field in SAMPLE_TABLE_AUDIT_FIELDS}
    output.update(
        {
            "source_id": source_id,
            "glds_prefix": str(row.get("glds_prefix", "")),
            "task_ids": str(row.get("task_ids", "")),
            "api_url": api_url,
            "api_response_sha256": listing_result.sha256,
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

    sample_files = [
        (filename, metadata)
        for filename, metadata in sorted(files.items())
        if is_sample_table_file(filename)
    ]
    output["api_status"] = "ok"
    output["n_files"] = str(len(files))
    output["sample_table_file_count"] = str(len(sample_files))
    output["sample_table_files"] = _join([filename for filename, _ in sample_files])
    output["sample_table_urls"] = _join([_download_url(metadata) for _, metadata in sample_files])

    if not sample_files:
        output.update(
            {
                "audit_status": "no_sample_table_listed",
                "pending_reason": "OSDR file listing returned no sample table-like files.",
            }
        )
        return output

    all_columns: list[str] = []
    all_rows: list[dict[str, str]] = []
    fetched_count = 0
    errors: list[str] = []
    for filename, metadata in sample_files:
        url = _download_url(metadata)
        if not url:
            errors.append(f"{filename}: missing download URL")
            continue
        sample_result = fetch(url)
        if not sample_result.ok:
            errors.append(f"{filename}: {sample_result.error}")
            continue
        fetched_count += 1
        columns, parsed_rows = parse_sample_table(
            sample_result.body.decode("utf-8", errors="replace"),
            filename=filename,
        )
        all_columns.extend(columns)
        all_rows.extend(parsed_rows)

    unique_columns = sorted({column for column in all_columns if column})
    condition_columns = _condition_columns(unique_columns)
    output["fetched_sample_table_count"] = str(fetched_count)
    output["sample_rows"] = str(len(all_rows))
    output["sample_table_columns"] = _join(unique_columns)
    output["candidate_condition_columns"] = _join(condition_columns)
    output["condition_value_summary"] = _condition_summary(all_rows, condition_columns)

    if all_rows:
        output.update(
            {
                "audit_status": "sample_table_parsed",
                "pending_reason": (
                    "Sample table parsed; manual factor mapping is still required before "
                    "final fold generation."
                ),
            }
        )
    elif errors:
        output.update(
            {
                "audit_status": "sample_table_fetch_or_parse_error",
                "pending_reason": " | ".join(errors),
            }
        )
    else:
        output.update(
            {
                "audit_status": "sample_table_listed_unparsed",
                "pending_reason": "Sample table files were listed but no rows were parsed.",
            }
        )
    return output


def audit_sample_table_inventory(
    rows: Sequence[Mapping[str, str]],
    *,
    fetcher: Callable[[str], FetchResult] | None = None,
    api_base: str = BIODATA_API_BASE,
) -> list[dict[str, str]]:
    return [
        audit_sample_table_row(row, fetcher=fetcher, api_base=api_base)
        for row in rows
    ]


def build_sample_factor_rows(
    rows: Sequence[Mapping[str, str]],
    *,
    fetcher: Callable[[str], FetchResult] | None = None,
    api_base: str = BIODATA_API_BASE,
) -> list[dict[str, str]]:
    """Fetch sample tables and expand compact condition strings into factors."""

    fetch = fetcher or (lambda url: fetch_url(url, timeout=30, max_bytes=None))
    factor_rows: list[dict[str, str]] = []
    for source_row in rows:
        source_id = str(source_row["source_id"])
        listing_result = fetch(source_file_api_url(source_id, api_base=api_base))
        if not listing_result.ok:
            factor_rows.append(
                {
                    "source_id": source_id,
                    "glds_prefix": str(source_row.get("glds_prefix", "")),
                    "parse_status": f"api_error: {listing_result.error}",
                }
            )
            continue
        try:
            listing = json.loads(listing_result.body.decode("utf-8"))
            files = extract_files_from_listing(source_id, listing)
        except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
            factor_rows.append(
                {
                    "source_id": source_id,
                    "glds_prefix": str(source_row.get("glds_prefix", "")),
                    "parse_status": f"api_error: {exc}",
                }
            )
            continue
        sample_files = [
            (filename, metadata)
            for filename, metadata in sorted(files.items())
            if is_sample_table_file(filename)
        ]
        for filename, metadata in sample_files:
            url = _download_url(metadata)
            if not url:
                factor_rows.append(
                    {
                        "source_id": source_id,
                        "glds_prefix": str(source_row.get("glds_prefix", "")),
                        "sample_table_file": filename,
                        "parse_status": "missing_download_url",
                    }
                )
                continue
            sample_result = fetch(url)
            if not sample_result.ok:
                factor_rows.append(
                    {
                        "source_id": source_id,
                        "glds_prefix": str(source_row.get("glds_prefix", "")),
                        "sample_table_file": filename,
                        "parse_status": f"sample_table_fetch_error: {sample_result.error}",
                    }
                )
                continue
            _, sample_rows = parse_sample_table(
                sample_result.body.decode("utf-8", errors="replace"),
                filename=filename,
            )
            for sample_row in sample_rows:
                factors = parse_condition_factors(str(sample_row.get("condition", "") or ""))
                factor_rows.append(
                    {
                        "source_id": source_id,
                        "glds_prefix": str(source_row.get("glds_prefix", "")),
                        "sample_id": _sample_id_from_row(sample_row),
                        "organoid_type": str(source_row.get("organoid_type", "")),
                        "sample_table_file": filename,
                        **factors,
                    }
                )
    return [
        {field: str(row.get(field, "") or "") for field in SAMPLE_FACTOR_FIELDS}
        for row in factor_rows
    ]


def write_sample_table_audit(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    if not rows:
        raise ValueError("cannot write an empty sample table audit")
    output_csv = Path(csv_path)
    output_json = Path(json_path)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    normalized = [
        {field: str(row.get(field, "") or "") for field in SAMPLE_TABLE_AUDIT_FIELDS}
        for row in rows
    ]
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SAMPLE_TABLE_AUDIT_FIELDS)
        writer.writeheader()
        writer.writerows(normalized)
    output_json.write_text(json.dumps(normalized, indent=2, sort_keys=True) + "\n")
    return output_csv, output_json


def write_sample_factor_table(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    if not rows:
        raise ValueError("cannot write an empty sample factor table")
    output_csv = Path(csv_path)
    output_json = Path(json_path)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    normalized = [
        {field: str(row.get(field, "") or "") for field in SAMPLE_FACTOR_FIELDS}
        for row in rows
    ]
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SAMPLE_FACTOR_FIELDS)
        writer.writeheader()
        writer.writerows(normalized)
    output_json.write_text(json.dumps(normalized, indent=2, sort_keys=True) + "\n")
    return output_csv, output_json

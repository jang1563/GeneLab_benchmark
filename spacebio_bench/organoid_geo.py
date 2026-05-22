"""GEO metadata helpers for the draft human organoid extension task."""

from __future__ import annotations

import csv
import gzip
import json
import re
from pathlib import Path
from typing import Mapping, Sequence

from .sample_table_audit import write_sample_factor_table


GSE259421_SERIES_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE259nnn/GSE259421/matrix/"
    "GSE259421_series_matrix.txt.gz"
)
GEO_SAMPLE_METADATA_FIELDS = [
    "sample_id",
    "geo_accession",
    "geo_series_accession",
    "geo_pubmed_id",
    "geo_sample_title",
    "geo_subject_id",
    "donor_or_line_id",
    "ipsc_line_id",
    "geo_source_name_ch1",
    "geo_cell_line",
    "geo_cell_type",
    "geo_treatment",
    "geo_title_organoid_type",
    "geo_title_microglia_condition",
    "geo_title_spaceflight_condition",
    "geo_title_replicate",
    "geo_biosample_accession",
    "geo_sra_accession",
    "geo_metadata_status",
    "geo_metadata_source",
]
GEO_SAMPLE_FACTOR_FIELDS = [
    "donor_or_line_id",
    "ipsc_line_id",
    "geo_sample_title",
    "geo_cell_type",
    "geo_treatment",
    "geo_metadata_status",
]
_TITLE_RE = re.compile(
    r"^(?P<subject>Subject\d+),\s*"
    r"(?P<organoid>.+?\s+organoid),\s*"
    r"(?P<microglia>with|without)\s+microglia\s*,\s*"
    r"(?P<condition>Ground|LEO)\s*(?P<replicate>\d+)\s*$",
    re.IGNORECASE,
)


def _read_text(path: str | Path) -> str:
    input_path = Path(path)
    if input_path.suffix == ".gz":
        with gzip.open(input_path, "rt", encoding="utf-8", errors="replace") as handle:
            return handle.read()
    return input_path.read_text(encoding="utf-8", errors="replace")


def _parse_matrix_line(line: str) -> list[str]:
    return next(csv.reader([line], delimiter="\t", quotechar='"'))


def _series_scalar(sample_lines: Mapping[str, list[list[str]]], key: str) -> str:
    rows = sample_lines.get(key, [])
    if not rows:
        return ""
    values = rows[0]
    return values[0] if values else ""


def _normalize_organoid_type(value: str) -> str:
    lower = value.strip().lower()
    if "dopaminergic" in lower:
        return "dopaminergic_neural_organoid"
    if "cortical" in lower:
        return "cortical_neural_organoid"
    return re.sub(r"[^a-z0-9]+", "_", lower).strip("_")


def _normalize_microglia(value: str) -> str:
    lower = " ".join(value.strip().lower().split())
    if lower == "with":
        return "with_microglia"
    if lower == "without":
        return "without_microglia"
    if lower.startswith("with "):
        return "with_microglia"
    if lower.startswith("without "):
        return "without_microglia"
    return re.sub(r"[^a-z0-9]+", "_", lower).strip("_")


def _normalize_condition(value: str) -> str:
    lower = value.strip().lower()
    if lower == "leo":
        return "LEO_or_ISS"
    if lower == "ground":
        return "Ground"
    return value.strip()


def _parse_title(title: str) -> dict[str, str]:
    match = _TITLE_RE.match(" ".join(title.strip().split()))
    if not match:
        return {
            "geo_subject_id": "",
            "donor_or_line_id": "",
            "geo_title_organoid_type": "",
            "geo_title_microglia_condition": "",
            "geo_title_spaceflight_condition": "",
            "geo_title_replicate": "",
        }
    subject = match.group("subject")
    return {
        "geo_subject_id": subject,
        "donor_or_line_id": subject,
        "geo_title_organoid_type": _normalize_organoid_type(match.group("organoid")),
        "geo_title_microglia_condition": _normalize_microglia(match.group("microglia")),
        "geo_title_spaceflight_condition": _normalize_condition(match.group("condition")),
        "geo_title_replicate": match.group("replicate"),
    }


def _parse_key_value(value: str) -> tuple[str, str]:
    if ":" not in value:
        return "", value.strip()
    key, parsed_value = value.split(":", 1)
    return key.strip().lower(), parsed_value.strip()


def _parse_relation(value: str) -> tuple[str, str]:
    if ":" not in value:
        return "", value.strip()
    key, url = value.split(":", 1)
    text = url.strip()
    accession = text.rstrip("/").split("/")[-1]
    if "term=" in text:
        accession = text.split("term=", 1)[1].split("&", 1)[0]
    return key.strip().lower(), accession


def parse_geo_series_matrix(text: str) -> list[dict[str, str]]:
    """Parse sample-level metadata from a GEO series matrix text payload."""

    matrix_lines: dict[str, list[list[str]]] = {}
    series_fields: dict[str, str] = {}
    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        if line == "!series_matrix_table_begin":
            break
        if not line:
            continue
        cells = _parse_matrix_line(line)
        if not cells:
            continue
        key = cells[0]
        values = cells[1:]
        if key.startswith("!Series_"):
            if values and key not in series_fields:
                series_fields[key] = values[0]
        if key.startswith("!Sample_"):
            matrix_lines.setdefault(key, []).append(values)

    accessions = matrix_lines.get("!Sample_geo_accession", [[]])[0]
    if not accessions:
        raise ValueError("GEO series matrix is missing !Sample_geo_accession")
    n_samples = len(accessions)
    rows = [
        {
            "sample_id": accession,
            "geo_accession": accession,
            "geo_series_accession": series_fields.get("!Series_geo_accession", ""),
            "geo_pubmed_id": series_fields.get("!Series_pubmed_id", ""),
            "geo_metadata_source": GSE259421_SERIES_MATRIX_URL,
        }
        for accession in accessions
    ]

    def assign_once(key: str, output_field: str) -> None:
        values = matrix_lines.get(key, [[]])[0]
        if values and len(values) != n_samples:
            raise ValueError(f"{key} has {len(values)} values for {n_samples} samples")
        for index, value in enumerate(values):
            rows[index][output_field] = value

    assign_once("!Sample_title", "geo_sample_title")
    assign_once("!Sample_source_name_ch1", "geo_source_name_ch1")

    for values in matrix_lines.get("!Sample_characteristics_ch1", []):
        if len(values) != n_samples:
            raise ValueError(
                f"!Sample_characteristics_ch1 has {len(values)} values for {n_samples} samples"
            )
        for index, raw_value in enumerate(values):
            key, value = _parse_key_value(raw_value)
            if key == "cell line":
                rows[index]["geo_cell_line"] = value
            elif key == "cell type":
                rows[index]["geo_cell_type"] = value
            elif key == "treatment":
                rows[index]["geo_treatment"] = value

    for values in matrix_lines.get("!Sample_relation", []):
        if len(values) != n_samples:
            raise ValueError(f"!Sample_relation has {len(values)} values for {n_samples} samples")
        for index, raw_value in enumerate(values):
            key, accession = _parse_relation(raw_value)
            if key == "biosample":
                rows[index]["geo_biosample_accession"] = accession
            elif key == "sra":
                rows[index]["geo_sra_accession"] = accession

    normalized_rows: list[dict[str, str]] = []
    for row in rows:
        title_fields = _parse_title(row.get("geo_sample_title", ""))
        cell_line = row.get("geo_cell_line", "") or row.get("geo_source_name_ch1", "")
        row.update(title_fields)
        row["ipsc_line_id"] = cell_line
        row["geo_metadata_status"] = (
            "parsed_from_geo_series_matrix"
            if row.get("donor_or_line_id") and row.get("ipsc_line_id")
            else "geo_series_matrix_missing_donor_or_line"
        )
        normalized_rows.append(
            {field: str(row.get(field, "") or "") for field in GEO_SAMPLE_METADATA_FIELDS}
        )
    return normalized_rows


def read_geo_series_matrix(path: str | Path) -> list[dict[str, str]]:
    """Read a plain-text or gzipped GEO series matrix and parse sample metadata."""

    return parse_geo_series_matrix(_read_text(path))


def merge_geo_metadata_with_sample_factors(
    sample_factor_rows: Sequence[Mapping[str, str]],
    geo_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, str]]:
    """Attach GEO donor/iPSC-line metadata to existing sample-factor rows."""

    by_sample = {str(row.get("sample_id", "") or ""): row for row in geo_rows}
    merged_rows: list[dict[str, str]] = []
    for sample_factor in sample_factor_rows:
        sample_id = str(sample_factor.get("sample_id", "") or "")
        geo = by_sample.get(sample_id, {})
        merged = dict(sample_factor)
        for field in GEO_SAMPLE_FACTOR_FIELDS:
            merged[field] = str(geo.get(field, "") or "")
        if not geo:
            merged["geo_metadata_status"] = "geo_metadata_not_found_for_sample"
        merged_rows.append(merged)
    return merged_rows


def write_geo_sample_metadata(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    """Write parsed GEO sample metadata as CSV and JSON."""

    if not rows:
        raise ValueError("cannot write an empty GEO sample metadata table")
    output_csv = Path(csv_path)
    output_json = Path(json_path)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    normalized = [
        {field: str(row.get(field, "") or "") for field in GEO_SAMPLE_METADATA_FIELDS}
        for row in rows
    ]
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GEO_SAMPLE_METADATA_FIELDS)
        writer.writeheader()
        writer.writerows(normalized)
    output_json.write_text(json.dumps(normalized, indent=2, sort_keys=True) + "\n")
    return output_csv, output_json


def read_csv_rows(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_merged_sample_factors(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    return write_sample_factor_table(rows, csv_path=csv_path, json_path=json_path)

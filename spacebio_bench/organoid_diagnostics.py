"""Diagnostics for draft human organoid extension tasks."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from statistics import median
from typing import Mapping, Sequence

import numpy as np

from .data import HumanOrganoidTaskData


SAMPLE_DIAGNOSTIC_FIELDS = [
    "source_id",
    "glds_prefix",
    "sample_id",
    "true_label",
    "organoid_type",
    "microglia_condition",
    "disease_context",
    "donor_or_line_id",
    "ipsc_line_id",
    "matrix_path",
    "n_features",
    "sample_sum",
    "sample_mean",
    "sample_median",
    "nonzero_features",
    "zero_fraction",
]
GROUP_DIAGNOSTIC_FIELDS = [
    "group_field",
    "group_value",
    "n_samples",
    "sample_sum_min",
    "sample_sum_median",
    "sample_sum_max",
    "zero_fraction_median",
    "label_distribution",
]
DONOR_AUDIT_FIELDS = [
    "source_id",
    "glds_prefix",
    "sample_table_file",
    "sample_table_columns",
    "sample_rows",
    "donor_field_status",
    "donor_or_line_count",
    "donor_or_line_ids",
    "ipsc_line_ids",
    "blocking_status",
    "pending_reason",
]


def _format_float(value: float) -> str:
    return f"{value:.10g}"


def _factor_by_sample(task: HumanOrganoidTaskData) -> dict[str, dict[str, str]]:
    return {
        str(row["sample_id"]): row
        for row in task.sample_factors
        if row.get("parse_status") == "parsed" and row.get("sample_id")
    }


def build_organoid_sample_diagnostics(task: HumanOrganoidTaskData) -> list[dict[str, str]]:
    """Summarize per-sample expression scale diagnostics."""

    by_sample = _factor_by_sample(task)
    rows: list[dict[str, str]] = []
    source_matrix_paths = {source_id: path.as_posix() for source_id, path in task.matrix_paths.items()}
    for sample_id, values in task.features.iterrows():
        factor = by_sample[str(sample_id)]
        array = values.to_numpy(dtype=np.float64)
        nonzero = int(np.count_nonzero(array))
        rows.append(
            {
                "source_id": str(factor.get("source_id", "")),
                "glds_prefix": str(factor.get("glds_prefix", "")),
                "sample_id": str(sample_id),
                "true_label": str(factor.get("true_label", "")),
                "organoid_type": str(factor.get("organoid_type", "")),
                "microglia_condition": str(factor.get("microglia_condition", "")),
                "disease_context": str(factor.get("disease_context", "")),
                "donor_or_line_id": str(factor.get("donor_or_line_id", "")),
                "ipsc_line_id": str(factor.get("ipsc_line_id", "")),
                "matrix_path": source_matrix_paths.get(str(factor.get("source_id", "")), ""),
                "n_features": str(len(array)),
                "sample_sum": _format_float(float(array.sum())),
                "sample_mean": _format_float(float(array.mean())),
                "sample_median": _format_float(float(np.median(array))),
                "nonzero_features": str(nonzero),
                "zero_fraction": _format_float(1.0 - (nonzero / len(array))),
            }
        )
    return rows


def _label_distribution(rows: Sequence[Mapping[str, str]]) -> str:
    counts: dict[str, int] = {}
    for row in rows:
        label = str(row.get("true_label", "") or "")
        if label:
            counts[label] = counts.get(label, 0) + 1
    return ";".join(f"{label}:{count}" for label, count in sorted(counts.items()))


def build_organoid_group_diagnostics(
    sample_rows: Sequence[Mapping[str, str]],
    *,
    group_fields: Sequence[str] = (
        "source_id",
        "true_label",
        "organoid_type",
        "microglia_condition",
        "disease_context",
        "donor_or_line_id",
        "ipsc_line_id",
    ),
) -> list[dict[str, str]]:
    """Aggregate sample-total diagnostics by key sample factors."""

    output: list[dict[str, str]] = []
    for field in group_fields:
        values = sorted({str(row.get(field, "") or "") for row in sample_rows if row.get(field)})
        for value in values:
            group_rows = [row for row in sample_rows if row.get(field) == value]
            sample_sums = [float(row["sample_sum"]) for row in group_rows]
            zero_fractions = [float(row["zero_fraction"]) for row in group_rows]
            output.append(
                {
                    "group_field": field,
                    "group_value": value,
                    "n_samples": str(len(group_rows)),
                    "sample_sum_min": _format_float(min(sample_sums)),
                    "sample_sum_median": _format_float(median(sample_sums)),
                    "sample_sum_max": _format_float(max(sample_sums)),
                    "zero_fraction_median": _format_float(median(zero_fractions)),
                    "label_distribution": _label_distribution(group_rows),
                }
            )
    return output


def build_organoid_donor_metadata_audit(
    sample_table_audit_rows: Sequence[Mapping[str, str]],
    sample_factor_rows: Sequence[Mapping[str, str]] | None = None,
) -> list[dict[str, str]]:
    """Record whether donor/iPSC-line blocking metadata are available locally."""

    output: list[dict[str, str]] = []
    donor_keywords = ("donor", "subject", "individual", "patient", "line", "ipsc")
    factors_by_source: dict[str, list[Mapping[str, str]]] = {}
    for factor in sample_factor_rows or []:
        source_id = str(factor.get("source_id", "") or "")
        if source_id:
            factors_by_source.setdefault(source_id, []).append(factor)
    for row in sample_table_audit_rows:
        source_id = str(row.get("source_id", ""))
        source_factors = factors_by_source.get(source_id, [])
        donor_ids = sorted(
            {
                str(factor.get("donor_or_line_id", "") or "")
                for factor in source_factors
                if factor.get("donor_or_line_id")
            }
        )
        ipsc_line_ids = sorted(
            {
                str(factor.get("ipsc_line_id", "") or "")
                for factor in source_factors
                if factor.get("ipsc_line_id")
            }
        )
        columns = str(row.get("sample_table_columns", "") or "")
        lower_columns = columns.lower()
        has_donor_column = any(keyword in lower_columns for keyword in donor_keywords)
        if donor_ids:
            donor_status = "recovered_from_geo_series_matrix"
            blocking_status = "available_but_confounded_pilot_only"
            reason = (
                "GEO series matrix provides Subject-level and iPSC-line identifiers; "
                "donor holdout remains a pilot-only diagnostic because donor, source, "
                "organoid fate, and disease context are not independently crossed."
            )
        elif has_donor_column:
            donor_status = "candidate_donor_field_present"
            blocking_status = "pending_donor_field_parser"
            reason = "SampleTable has donor-like columns that need parser review."
        else:
            donor_status = "not_available_in_osdr_sample_table"
            blocking_status = "donor_block_named_but_unresolved"
            reason = "OSDR SampleTable exposes condition fields but no donor/iPSC-line field."
        output.append(
            {
                "source_id": source_id,
                "glds_prefix": str(row.get("glds_prefix", "")),
                "sample_table_file": str(row.get("sample_table_files", "")),
                "sample_table_columns": columns,
                "sample_rows": str(row.get("sample_rows", "")),
                "donor_field_status": donor_status,
                "donor_or_line_count": str(len(donor_ids)) if donor_ids else "",
                "donor_or_line_ids": ";".join(donor_ids),
                "ipsc_line_ids": ";".join(ipsc_line_ids),
                "blocking_status": blocking_status,
                "pending_reason": reason,
            }
        )
    return output


def _write_rows(
    rows: Sequence[Mapping[str, str]],
    *,
    fields: Sequence[str],
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    if not rows:
        raise ValueError("cannot write an empty diagnostics table")
    output_csv = Path(csv_path)
    output_json = Path(json_path)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    normalized = [{field: str(row.get(field, "") or "") for field in fields} for row in rows]
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(normalized)
    output_json.write_text(json.dumps(normalized, indent=2, sort_keys=True) + "\n")
    return output_csv, output_json


def write_organoid_sample_diagnostics(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    return _write_rows(
        rows,
        fields=SAMPLE_DIAGNOSTIC_FIELDS,
        csv_path=csv_path,
        json_path=json_path,
    )


def write_organoid_group_diagnostics(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    return _write_rows(
        rows,
        fields=GROUP_DIAGNOSTIC_FIELDS,
        csv_path=csv_path,
        json_path=json_path,
    )


def write_organoid_donor_metadata_audit(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    return _write_rows(
        rows,
        fields=DONOR_AUDIT_FIELDS,
        csv_path=csv_path,
        json_path=json_path,
    )

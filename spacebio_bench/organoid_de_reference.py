"""Derived DE-reference extraction for the draft human organoid task."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import re
from pathlib import Path
from typing import Callable, Mapping, Sequence

from .source_audit import FetchResult, fetch_url


DE_REFERENCE_FIELDS = [
    "source_id",
    "glds_prefix",
    "task_id",
    "organoid_type",
    "disease_context",
    "microglia_condition",
    "contrast_id",
    "ensembl_id",
    "gene_symbol",
    "entrez_id",
    "log2fc_leo_or_iss_minus_ground",
    "source_log2fc_value",
    "stat",
    "p_value",
    "adj_p_value",
    "significant_fdr_0_05",
    "source_de_file",
]

BENCHMARK_LOG2FC_ORIENTATION = "LEO_or_ISS_minus_Ground"
SOURCE_LOG2FC_ORIENTATION = "group_a_minus_group_b_assumed_from_Log2fc_A_v_B_header"
ORIENTATION_TRANSFORM = "negated_ground_control_v_space_flight_to_leo_or_iss_minus_ground"
RESPONSE_SIGNATURE_REQUIRED_COLUMNS = [
    "task_id",
    "source_id",
    "contrast_id",
    "gene_symbol",
    "ensembl_id",
    "predicted_log2fc_leo_or_iss_minus_ground",
]


def _split_semicolon(value: str) -> list[str]:
    return [part.strip() for part in str(value or "").split(";") if part.strip()]


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", value.strip().lower()).strip("_")
    return slug or "unknown"


def _header_tail(value: str) -> str:
    return str(value or "").strip().lstrip("\ufeff").rsplit("/", 1)[-1]


def _as_float(value: str) -> float | None:
    text = str(value or "").strip()
    if text.lower() in {"", "na", "nan", "inf", "-inf", "null", "none"}:
        return None
    try:
        parsed = float(text)
    except ValueError:
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def _format_float(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value:.12g}"


def _condition_from_group(group: str) -> dict[str, str]:
    parts = [part.strip() for part in group.strip().strip("()").split(" & ")]
    if len(parts) != 3:
        return {
            "disease_context": "",
            "flight_condition": "",
            "microglia_condition": "",
        }
    disease_context, flight_condition, microglia_condition = parts
    return {
        "disease_context": disease_context,
        "flight_condition": flight_condition,
        "microglia_condition": microglia_condition,
    }


def parse_contrast_label(label: str) -> dict[str, str]:
    """Parse an OSDR contrast label such as ``(Ground)v(Space Flight)``."""

    text = str(label or "").strip()
    if text.startswith("Log2fc_"):
        text = text.removeprefix("Log2fc_")
    match = re.match(r"^\((?P<a>.+)\)v\((?P<b>.+)\)$", text)
    if not match:
        return {
            "contrast_label": text,
            "source_group_a": "",
            "source_group_b": "",
            "disease_context": "",
            "microglia_condition": "",
            "is_direct_spaceflight": "false",
        }
    source_group_a = match.group("a").strip()
    source_group_b = match.group("b").strip()
    parsed_a = _condition_from_group(source_group_a)
    parsed_b = _condition_from_group(source_group_b)
    same_disease = parsed_a["disease_context"] == parsed_b["disease_context"]
    same_microglia = parsed_a["microglia_condition"] == parsed_b["microglia_condition"]
    is_direct = (
        same_disease
        and same_microglia
        and parsed_a["flight_condition"] == "Ground Control"
        and parsed_b["flight_condition"] == "Space Flight"
    )
    return {
        "contrast_label": text,
        "source_group_a": source_group_a,
        "source_group_b": source_group_b,
        "disease_context": parsed_a["disease_context"] if same_disease else "",
        "microglia_condition": parsed_a["microglia_condition"] if same_microglia else "",
        "is_direct_spaceflight": "true" if is_direct else "false",
    }


def contrast_id(source_id: str, label: str) -> str:
    parsed = parse_contrast_label(label)
    parts = [
        source_id,
        parsed.get("disease_context", ""),
        parsed.get("microglia_condition", ""),
        "leo_or_iss_vs_ground",
    ]
    return _slug("_".join(part for part in parts if part))


def select_canonical_de_reference(audit_row: Mapping[str, str]) -> dict[str, str]:
    """Select the non-rRNArm OSDR bulk-RNA differential-expression table."""

    files = _split_semicolon(str(audit_row.get("de_reference_files", "")))
    urls = _split_semicolon(str(audit_row.get("de_reference_urls", "")))
    sizes = _split_semicolon(str(audit_row.get("de_reference_file_sizes", "")))
    data_types = _split_semicolon(str(audit_row.get("de_reference_data_types", "")))
    candidates: list[dict[str, str]] = []
    for index, filename in enumerate(files):
        candidates.append(
            {
                "source_de_file": filename,
                "source_de_url": urls[index] if index < len(urls) else "",
                "source_de_file_size": sizes[index] if index < len(sizes) else "",
                "source_de_data_type": data_types[index] if index < len(data_types) else "",
            }
        )

    for candidate in candidates:
        filename = candidate["source_de_file"]
        if "differential_expression_GLbulkRNAseq.csv" in filename and "rRNArm" not in filename:
            return candidate
    for candidate in candidates:
        filename = candidate["source_de_file"]
        if "differential_expression" in filename and "rRNArm" not in filename:
            return candidate
    raise ValueError(
        f"no canonical non-rRNArm DE reference listed for {audit_row.get('source_id', '')}"
    )


def _contrast_columns(headers: Sequence[str], direct_contrasts: Sequence[str]) -> dict[str, dict[str, int]]:
    direct_set = set(direct_contrasts)
    columns: dict[str, dict[str, int]] = {
        contrast: {} for contrast in direct_contrasts
    }
    prefixes = {
        "Log2fc_": "log2fc",
        "Stat_": "stat",
        "P.value_": "p_value",
        "Adj.p.value_": "adj_p_value",
    }
    for index, header in enumerate(headers):
        tail = _header_tail(header)
        for prefix, metric_key in prefixes.items():
            if not tail.startswith(prefix):
                continue
            label = tail.removeprefix(prefix)
            if label in direct_set:
                columns[label][metric_key] = index
    missing = {
        label: sorted({"log2fc", "stat", "p_value", "adj_p_value"} - set(indices))
        for label, indices in columns.items()
        if {"log2fc", "stat", "p_value", "adj_p_value"} - set(indices)
    }
    if missing:
        formatted = "; ".join(
            f"{label}: {','.join(keys)}" for label, keys in sorted(missing.items())
        )
        raise ValueError(f"DE table is missing selected direct contrast columns: {formatted}")
    return columns


def _gene_column_indices(headers: Sequence[str]) -> dict[str, int]:
    tails = {_header_tail(header): index for index, header in enumerate(headers)}
    required = {
        "ensembl_id": "ENSEMBL",
        "gene_symbol": "SYMBOL",
        "entrez_id": "ENTREZID",
    }
    return {
        output_key: tails[source_key]
        for output_key, source_key in required.items()
        if source_key in tails
    }


def _row_value(row: Sequence[str], index: int | None) -> str:
    if index is None or index >= len(row):
        return ""
    return str(row[index] or "").strip()


def _is_gene_data_row(row: Sequence[str], gene_columns: Mapping[str, int]) -> bool:
    values = [
        _row_value(row, gene_columns.get("ensembl_id")),
        _row_value(row, gene_columns.get("gene_symbol")),
        _row_value(row, gene_columns.get("entrez_id")),
    ]
    if not any(values):
        return False
    type_tokens = {"string", "number", "numeric", "integer", "double"}
    return not all(value.lower() in type_tokens for value in values if value)


def extract_organoid_de_reference_rows(
    *,
    source_row: Mapping[str, str],
    audit_row: Mapping[str, str],
    de_csv_bytes: bytes,
    source_de_file: str,
    source_de_url: str,
    source_de_sha256: str | None = None,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    """Extract selected organoid DE contrasts from one OSDR DE table."""

    direct_contrasts = _split_semicolon(str(audit_row.get("direct_spaceflight_contrasts", "")))
    if not direct_contrasts:
        raise ValueError(f"no direct spaceflight contrasts listed for {audit_row.get('source_id', '')}")
    source_id = str(audit_row.get("source_id") or source_row.get("source_id") or "")
    glds_prefix = str(audit_row.get("glds_prefix") or source_row.get("glds_prefix") or "")
    task_id = _split_semicolon(str(audit_row.get("task_ids", "")))
    sha256 = source_de_sha256 or hashlib.sha256(de_csv_bytes).hexdigest()

    text = io.TextIOWrapper(
        io.BytesIO(de_csv_bytes),
        encoding="utf-8-sig",
        errors="replace",
        newline="",
    )
    reader = csv.reader(text)
    try:
        headers = next(reader)
    except StopIteration as exc:
        raise ValueError("DE CSV has no header") from exc
    gene_columns = _gene_column_indices(headers)
    selected_columns = _contrast_columns(headers, direct_contrasts)

    output_rows: list[dict[str, str]] = []
    counts_by_contrast = {label: 0 for label in direct_contrasts}
    significant_by_contrast = {label: 0 for label in direct_contrasts}
    for row in reader:
        if not _is_gene_data_row(row, gene_columns):
            continue
        gene_values = {
            "ensembl_id": _row_value(row, gene_columns.get("ensembl_id")),
            "gene_symbol": _row_value(row, gene_columns.get("gene_symbol")),
            "entrez_id": _row_value(row, gene_columns.get("entrez_id")),
        }
        for label, columns in selected_columns.items():
            source_log2fc = _as_float(_row_value(row, columns.get("log2fc")))
            if source_log2fc is None:
                continue
            normalized_log2fc = -source_log2fc
            adj_p_value = _as_float(_row_value(row, columns.get("adj_p_value")))
            significant = adj_p_value is not None and adj_p_value <= 0.05
            parsed = parse_contrast_label(label)
            counts_by_contrast[label] += 1
            if significant:
                significant_by_contrast[label] += 1
            output_rows.append(
                {
                    "source_id": source_id,
                    "glds_prefix": glds_prefix,
                    "task_id": task_id[0] if task_id else "draft_human_organoid_spaceflight",
                    "organoid_type": str(source_row.get("organoid_type", "") or ""),
                    "disease_context": parsed["disease_context"],
                    "microglia_condition": parsed["microglia_condition"],
                    "contrast_id": contrast_id(source_id, label),
                    **gene_values,
                    "log2fc_leo_or_iss_minus_ground": _format_float(normalized_log2fc),
                    "source_log2fc_value": _format_float(source_log2fc),
                    "stat": _row_value(row, columns.get("stat")),
                    "p_value": _row_value(row, columns.get("p_value")),
                    "adj_p_value": _row_value(row, columns.get("adj_p_value")),
                    "significant_fdr_0_05": "true" if significant else "false",
                    "source_de_file": source_de_file,
                }
            )

    metadata = {
        "source_id": source_id,
        "glds_prefix": glds_prefix,
        "source_de_file": source_de_file,
        "source_de_url": source_de_url,
        "source_de_sha256": sha256,
        "direct_spaceflight_contrasts": direct_contrasts,
        "contrasts": [
            {
                "contrast_id": contrast_id(source_id, label),
                **{
                    key: value
                    for key, value in parse_contrast_label(label).items()
                    if key != "is_direct_spaceflight"
                },
                "source_log2fc_orientation": SOURCE_LOG2FC_ORIENTATION,
                "orientation_transform": ORIENTATION_TRANSFORM,
                "benchmark_log2fc_orientation": BENCHMARK_LOG2FC_ORIENTATION,
            }
            for label in direct_contrasts
        ],
        "n_contrasts": len(direct_contrasts),
        "n_rows": len(output_rows),
        "n_significant_fdr_0_05": sum(significant_by_contrast.values()),
        "rows_by_contrast": {
            contrast_id(source_id, label): counts_by_contrast[label]
            for label in direct_contrasts
        },
        "significant_fdr_0_05_by_contrast": {
            contrast_id(source_id, label): significant_by_contrast[label]
            for label in direct_contrasts
        },
    }
    return output_rows, metadata


def _default_fetcher(url: str) -> FetchResult:
    return fetch_url(url, timeout=300)


def build_organoid_de_reference(
    *,
    source_rows: Sequence[Mapping[str, str]],
    audit_rows: Sequence[Mapping[str, str]],
    fetcher: Callable[[str], FetchResult] = _default_fetcher,
) -> tuple[list[dict[str, str]], dict[str, object]]:
    """Build compact derived rows and manifest metadata for organoid DE references."""

    sources_by_id = {str(row.get("source_id", "")): row for row in source_rows}
    all_rows: list[dict[str, str]] = []
    source_metadata: list[dict[str, object]] = []
    for audit_row in sorted(audit_rows, key=lambda row: str(row.get("source_id", ""))):
        source_id = str(audit_row.get("source_id", ""))
        if not source_id:
            continue
        source_row = sources_by_id.get(source_id, {})
        canonical = select_canonical_de_reference(audit_row)
        if not canonical["source_de_url"]:
            raise ValueError(f"canonical DE reference for {source_id} has no download URL")
        fetched = fetcher(canonical["source_de_url"])
        if not fetched.ok:
            raise RuntimeError(
                f"failed to fetch canonical DE reference for {source_id}: {fetched.error}"
            )
        rows, metadata = extract_organoid_de_reference_rows(
            source_row=source_row,
            audit_row=audit_row,
            de_csv_bytes=fetched.body,
            source_de_file=canonical["source_de_file"],
            source_de_url=canonical["source_de_url"],
            source_de_sha256=fetched.sha256,
        )
        metadata.update(
            {
                "source_de_file_size": canonical.get("source_de_file_size", ""),
                "source_de_data_type": canonical.get("source_de_data_type", ""),
            }
        )
        source_metadata.append(metadata)
        all_rows.extend(rows)

    manifest = {
        "schema_version": "0.1.0",
        "artifact_type": "human_organoid_de_reference",
        "release_status": "draft_not_frozen",
        "build_policy": "v9-org-014-freeze-direct-non-rRNArm-de-contrasts",
        "canonical_de_file_rule": (
            "Select *_rna_seq_differential_expression_GLbulkRNAseq.csv and exclude rRNArm."
        ),
        "contrast_rule": (
            "Use direct matched Ground Control versus Space Flight contrasts with the same "
            "disease context and microglia condition."
        ),
        "log2fc_orientation_contract": {
            "source_orientation": SOURCE_LOG2FC_ORIENTATION,
            "benchmark_orientation": BENCHMARK_LOG2FC_ORIENTATION,
            "transform": ORIENTATION_TRANSFORM,
            "note": (
                "OSDR direct contrast labels are Ground Control v Space Flight. "
                "The derived table negates source Log2fc values so positive values "
                "mean higher in LEO/ISS than matched ground control."
            ),
        },
        "feature_key_contract": {
            "primary_key": "gene_symbol",
            "fallback_keys": ["ensembl_id", "entrez_id"],
            "feature_namespace": "human_gene",
        },
        "response_signature_contract": {
            "artifact": "response_signature.csv",
            "required_columns": RESPONSE_SIGNATURE_REQUIRED_COLUMNS,
            "predicted_log2fc_orientation": BENCHMARK_LOG2FC_ORIENTATION,
            "join_keys": ["source_id", "contrast_id", "gene_symbol"],
            "fallback_join_keys": ["source_id", "contrast_id", "ensembl_id"],
        },
        "totals": {
            "n_sources": len(source_metadata),
            "n_contrasts": sum(int(source["n_contrasts"]) for source in source_metadata),
            "n_rows": len(all_rows),
            "n_significant_fdr_0_05": sum(
                int(source["n_significant_fdr_0_05"]) for source in source_metadata
            ),
        },
        "sources": source_metadata,
    }
    return all_rows, manifest


def write_organoid_de_reference(
    *,
    source_rows: Sequence[Mapping[str, str]],
    audit_rows: Sequence[Mapping[str, str]],
    csv_path: str | Path,
    json_path: str | Path,
    fetcher: Callable[[str], FetchResult] = _default_fetcher,
) -> tuple[Path, Path]:
    """Write the derived human organoid DE reference CSV and manifest JSON."""

    rows, manifest = build_organoid_de_reference(
        source_rows=source_rows,
        audit_rows=audit_rows,
        fetcher=fetcher,
    )
    output_csv = Path(csv_path)
    output_json = Path(json_path)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    manifest["reference_table"] = output_csv.as_posix()
    manifest["reference_manifest"] = output_json.as_posix()
    if output_csv.name.endswith(".gz"):
        handle_context = gzip.open(output_csv, "wt", newline="")
    else:
        handle_context = output_csv.open("w", newline="")
    with handle_context as handle:
        writer = csv.DictWriter(handle, fieldnames=DE_REFERENCE_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in DE_REFERENCE_FIELDS})
    output_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return output_csv, output_json

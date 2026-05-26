"""Response-signature metrics for SpaceBio-Bench tasks."""

from __future__ import annotations

import csv
import gzip
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from .organoid_de_reference import RESPONSE_SIGNATURE_REQUIRED_COLUMNS


SIGNATURE_METRIC_IDS = {"de_direction_match", "signature_rank_correlation"}
REFERENCE_LOG2FC_COLUMN = "log2fc_leo_or_iss_minus_ground"
PREDICTED_LOG2FC_COLUMN = "predicted_log2fc_leo_or_iss_minus_ground"


def _open_text(path: str | Path):
    parsed = Path(path)
    if parsed.name.endswith(".gz"):
        return gzip.open(parsed, "rt", newline="")
    return parsed.open(newline="")


def _read_rows(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    with _open_text(path) as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        return fieldnames, [dict(row) for row in reader]


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


def _metric_value(value: float, details: Mapping[str, Any]) -> dict[str, Any]:
    return {"status": "computed", "value": value, "details": dict(details)}


def _metric_skip(reason: str, details: Mapping[str, Any] | None = None) -> dict[str, Any]:
    payload: dict[str, Any] = {"status": "skipped", "reason": reason}
    if details is not None:
        payload["details"] = dict(details)
    return payload


def _average_ranks(values: Sequence[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    cursor = 0
    while cursor < len(indexed):
        end = cursor + 1
        while end < len(indexed) and indexed[end][1] == indexed[cursor][1]:
            end += 1
        average_rank = (cursor + 1 + end) / 2.0
        for index, _ in indexed[cursor:end]:
            ranks[index] = average_rank
        cursor = end
    return ranks


def _spearman(xs: Sequence[float], ys: Sequence[float]) -> float | None:
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    x_ranks = _average_ranks(xs)
    y_ranks = _average_ranks(ys)
    mean_x = sum(x_ranks) / len(x_ranks)
    mean_y = sum(y_ranks) / len(y_ranks)
    numerator = sum((x - mean_x) * (y - mean_y) for x, y in zip(x_ranks, y_ranks))
    denom_x = sum((x - mean_x) ** 2 for x in x_ranks)
    denom_y = sum((y - mean_y) ** 2 for y in y_ranks)
    denominator = math.sqrt(denom_x * denom_y)
    if denominator == 0:
        return None
    return numerator / denominator


def _sign(value: float) -> int:
    if value > 0:
        return 1
    if value < 0:
        return -1
    return 0


def _key(row: Mapping[str, str], feature_key: str) -> tuple[str, str, str] | None:
    source_id = str(row.get("source_id", "") or "").strip()
    contrast_id = str(row.get("contrast_id", "") or "").strip()
    feature_value = str(row.get(feature_key, "") or "").strip()
    if not source_id or not contrast_id or not feature_value:
        return None
    return (source_id, contrast_id, feature_value)


def _validate_response_signature(
    *,
    manifest: Mapping[str, Any],
    response_signature_path: str | Path,
) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    path = Path(response_signature_path)
    try:
        columns, rows = _read_rows(path)
    except FileNotFoundError:
        return {
            "ok": False,
            "path": path.as_posix(),
            "n_rows": 0,
            "columns": [],
            "errors": [f"response_signature.csv artifact not found: {path}"],
            "warnings": [],
        }
    except OSError as exc:
        return {
            "ok": False,
            "path": path.as_posix(),
            "n_rows": 0,
            "columns": [],
            "errors": [f"could not read response_signature.csv {path}: {exc}"],
            "warnings": [],
        }

    missing = [column for column in RESPONSE_SIGNATURE_REQUIRED_COLUMNS if column not in columns]
    if missing:
        errors.append(f"missing required response_signature columns: {', '.join(missing)}")
    if columns and not rows:
        errors.append("response_signature.csv has no rows")

    task_id = str(manifest.get("task_id", "") or "")
    if "task_id" in columns:
        mismatches = sorted(
            {
                str(row.get("task_id", "") or "")
                for row in rows
                if row.get("task_id") and row.get("task_id") != task_id
            }
        )
        if mismatches:
            errors.append(
                "response_signature task_id column does not match manifest "
                f"task_id {task_id}: {', '.join(mismatches)}"
            )

    if PREDICTED_LOG2FC_COLUMN in columns:
        for row_index, row in enumerate(rows, start=2):
            if _as_float(str(row.get(PREDICTED_LOG2FC_COLUMN, ""))) is None:
                errors.append(
                    f"{PREDICTED_LOG2FC_COLUMN} is not numeric at CSV row {row_index}"
                )
                if len(errors) >= 20:
                    errors.append("additional numeric validation errors suppressed")
                    break

    if {"gene_symbol", "ensembl_id"} <= set(columns):
        no_feature_rows = [
            index
            for index, row in enumerate(rows, start=2)
            if not str(row.get("gene_symbol", "") or "").strip()
            and not str(row.get("ensembl_id", "") or "").strip()
        ]
        if no_feature_rows:
            preview = ", ".join(str(index) for index in no_feature_rows[:10])
            errors.append(
                "response_signature rows must include gene_symbol or ensembl_id; "
                f"missing at CSV rows: {preview}"
            )
            if len(no_feature_rows) > 10:
                warnings.append("additional feature-key validation errors suppressed")

    return {
        "ok": not errors,
        "path": path.as_posix(),
        "n_rows": len(rows),
        "columns": columns,
        "errors": errors,
        "warnings": warnings,
    }


def _resolve_reference_path(
    *,
    manifest: Mapping[str, Any],
    reference_signature_path: str | Path | None,
) -> Path | None:
    if reference_signature_path is not None:
        return Path(reference_signature_path)
    reference_signatures = manifest.get("reference_signatures", {})
    if isinstance(reference_signatures, Mapping) and reference_signatures.get("de_reference_table"):
        return Path(str(reference_signatures["de_reference_table"]))
    return None


def _load_reference_indices(path: Path) -> tuple[dict[tuple[str, str, str], Mapping[str, str]], dict[tuple[str, str, str], Mapping[str, str]], dict[str, Any]]:
    columns, rows = _read_rows(path)
    primary: dict[tuple[str, str, str], Mapping[str, str]] = {}
    fallback: dict[tuple[str, str, str], Mapping[str, str]] = {}
    duplicate_primary = 0
    duplicate_fallback = 0
    for row in rows:
        primary_key = _key(row, "gene_symbol")
        if primary_key is not None:
            if primary_key in primary:
                duplicate_primary += 1
            else:
                primary[primary_key] = row
        fallback_key = _key(row, "ensembl_id")
        if fallback_key is not None:
            if fallback_key in fallback:
                duplicate_fallback += 1
            else:
                fallback[fallback_key] = row
    summary = {
        "reference_path": path.as_posix(),
        "reference_columns": columns,
        "n_reference_rows": len(rows),
        "n_primary_keys": len(primary),
        "n_fallback_keys": len(fallback),
        "duplicate_primary_keys": duplicate_primary,
        "duplicate_fallback_keys": duplicate_fallback,
    }
    return primary, fallback, summary


def _joined_signature_rows(
    *,
    response_rows: Sequence[Mapping[str, str]],
    reference_primary: Mapping[tuple[str, str, str], Mapping[str, str]],
    reference_fallback: Mapping[tuple[str, str, str], Mapping[str, str]],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    joined: list[dict[str, Any]] = []
    counts = {
        "n_response_rows": len(response_rows),
        "n_joined_rows": 0,
        "n_primary_join_rows": 0,
        "n_fallback_join_rows": 0,
        "n_unmatched_rows": 0,
        "n_duplicate_response_keys": 0,
        "n_non_numeric_reference_rows": 0,
    }
    seen_response_keys: set[tuple[str, str, str, str]] = set()
    for response_row in response_rows:
        predicted = _as_float(str(response_row.get(PREDICTED_LOG2FC_COLUMN, "")))
        if predicted is None:
            continue
        join_kind = ""
        reference_row: Mapping[str, str] | None = None
        response_key = _key(response_row, "gene_symbol")
        if response_key is not None and response_key in reference_primary:
            reference_row = reference_primary[response_key]
            join_kind = "gene_symbol"
        else:
            response_key = _key(response_row, "ensembl_id")
            if response_key is not None and response_key in reference_fallback:
                reference_row = reference_fallback[response_key]
                join_kind = "ensembl_id"
        if response_key is None or reference_row is None:
            counts["n_unmatched_rows"] += 1
            continue
        duplicate_key = (
            join_kind,
            response_key[0],
            response_key[1],
            response_key[2],
        )
        if duplicate_key in seen_response_keys:
            counts["n_duplicate_response_keys"] += 1
            continue
        seen_response_keys.add(duplicate_key)
        reference = _as_float(str(reference_row.get(REFERENCE_LOG2FC_COLUMN, "")))
        if reference is None:
            counts["n_non_numeric_reference_rows"] += 1
            continue
        counts["n_joined_rows"] += 1
        if join_kind == "gene_symbol":
            counts["n_primary_join_rows"] += 1
        elif join_kind == "ensembl_id":
            counts["n_fallback_join_rows"] += 1
        joined.append(
            {
                "source_id": response_key[0],
                "contrast_id": response_key[1],
                "feature_id": response_key[2],
                "join_key": join_kind,
                "predicted_log2fc": predicted,
                "reference_log2fc": reference,
                "significant_fdr_0_05": (
                    str(reference_row.get("significant_fdr_0_05", "")).lower() == "true"
                ),
            }
        )
    return joined, counts


def _direction_match(joined: Sequence[Mapping[str, Any]]) -> tuple[float | None, dict[str, Any]]:
    significant = [row for row in joined if row.get("significant_fdr_0_05")]
    scored = [
        row
        for row in significant
        if _sign(float(row["predicted_log2fc"])) != 0
        and _sign(float(row["reference_log2fc"])) != 0
    ]
    matches = sum(
        _sign(float(row["predicted_log2fc"])) == _sign(float(row["reference_log2fc"]))
        for row in scored
    )
    details = {
        "n_joined_rows": len(joined),
        "n_significant_reference_rows": len(significant),
        "n_direction_scored": len(scored),
        "n_direction_matches": matches,
        "reference_filter": "significant_fdr_0_05",
    }
    if not scored:
        return None, details
    return matches / len(scored), details


def _rank_correlation(joined: Sequence[Mapping[str, Any]]) -> tuple[float | None, dict[str, Any]]:
    predicted = [float(row["predicted_log2fc"]) for row in joined]
    reference = [float(row["reference_log2fc"]) for row in joined]
    value = _spearman(predicted, reference)
    details = {
        "n_joined_rows": len(joined),
        "n_rank_scored": len(predicted),
        "correlation": "spearman",
        "reference_filter": "all_joined_rows",
    }
    return value, details


def _per_contrast_details(joined: Sequence[Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    by_contrast: dict[str, list[Mapping[str, Any]]] = {}
    for row in joined:
        by_contrast.setdefault(str(row["contrast_id"]), []).append(row)
    details: dict[str, dict[str, Any]] = {}
    for contrast_id, rows in sorted(by_contrast.items()):
        direction_value, direction_details = _direction_match(rows)
        rank_value, rank_details = _rank_correlation(rows)
        details[contrast_id] = {
            "n_joined_rows": len(rows),
            "de_direction_match": direction_value,
            "de_direction_match_details": direction_details,
            "signature_rank_correlation": rank_value,
            "signature_rank_correlation_details": rank_details,
        }
    return details


def compute_response_signature_metrics(
    *,
    manifest: Mapping[str, Any],
    response_signature_path: str | Path,
    reference_signature_path: str | Path | None = None,
) -> dict[str, Any]:
    """Validate and score a response-signature artifact against a DE reference."""

    validation = _validate_response_signature(
        manifest=manifest,
        response_signature_path=response_signature_path,
    )
    metrics: dict[str, Any] = {}
    if not validation["ok"]:
        reason = "response_signature.csv validation failed: " + "; ".join(validation["errors"])
        for metric_id in sorted(SIGNATURE_METRIC_IDS):
            metrics[metric_id] = _metric_skip(reason, {"validation": validation})
        return {"validation": validation, "metrics": metrics}

    reference_path = _resolve_reference_path(
        manifest=manifest,
        reference_signature_path=reference_signature_path,
    )
    if reference_path is None:
        reason = "DE reference table path missing from manifest reference_signatures"
        for metric_id in sorted(SIGNATURE_METRIC_IDS):
            metrics[metric_id] = _metric_skip(reason, {"validation": validation})
        return {"validation": validation, "metrics": metrics}
    if not reference_path.exists():
        reason = f"DE reference table not found: {reference_path}"
        for metric_id in sorted(SIGNATURE_METRIC_IDS):
            metrics[metric_id] = _metric_skip(reason, {"validation": validation})
        return {"validation": validation, "metrics": metrics}

    _, response_rows = _read_rows(response_signature_path)
    reference_primary, reference_fallback, reference_summary = _load_reference_indices(reference_path)
    joined, join_counts = _joined_signature_rows(
        response_rows=response_rows,
        reference_primary=reference_primary,
        reference_fallback=reference_fallback,
    )
    shared_details = {
        "validation": validation,
        "reference": reference_summary,
        "join": join_counts,
        "per_contrast": _per_contrast_details(joined),
    }
    if not joined:
        reason = "response_signature.csv had no rows that joined to the DE reference"
        for metric_id in sorted(SIGNATURE_METRIC_IDS):
            metrics[metric_id] = _metric_skip(reason, shared_details)
        return {"validation": validation, "metrics": metrics}

    direction_value, direction_details = _direction_match(joined)
    rank_value, rank_details = _rank_correlation(joined)
    if direction_value is None:
        metrics["de_direction_match"] = _metric_skip(
            "de_direction_match requires at least one joined significant reference gene "
            "with non-zero predicted and reference log2FC",
            {**shared_details, "aggregate": direction_details},
        )
    else:
        metrics["de_direction_match"] = _metric_value(
            direction_value,
            {**shared_details, "aggregate": direction_details},
        )

    if rank_value is None:
        metrics["signature_rank_correlation"] = _metric_skip(
            "signature_rank_correlation requires at least two joined genes with "
            "variation in predicted and reference log2FC",
            {**shared_details, "aggregate": rank_details},
        )
    else:
        metrics["signature_rank_correlation"] = _metric_value(
            rank_value,
            {**shared_details, "aggregate": rank_details},
        )
    return {"validation": validation, "metrics": metrics}

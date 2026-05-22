"""Submission validation for SpaceBio-Bench prediction files."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping


@dataclass
class SubmissionValidationReport:
    """Structured validation report for one prediction file."""

    ok: bool
    submission_path: str
    task_id: str
    n_rows: int = 0
    columns: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    extra_columns: list[str] = field(default_factory=list)
    embedding_columns: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "ok": self.ok,
            "submission_path": self.submission_path,
            "task_id": self.task_id,
            "n_rows": self.n_rows,
            "columns": self.columns,
            "errors": self.errors,
            "warnings": self.warnings,
            "extra_columns": self.extra_columns,
            "embedding_columns": self.embedding_columns,
        }


def _read_submission(path: str | Path) -> tuple[list[str], list[dict[str, str]]]:
    submission_path = Path(path)
    with submission_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]
    return fieldnames, rows


def _expected_labels(manifest: Mapping[str, Any]) -> set[str]:
    output = manifest.get("output", {})
    if isinstance(output, Mapping) and output.get("label_domain"):
        return {str(label) for label in output["label_domain"]}

    labels: set[str] = set()
    split = manifest.get("split", {})
    folds = split.get("folds", []) if isinstance(split, Mapping) else []
    for fold in folds:
        if not isinstance(fold, Mapping):
            continue
        for key in ("train_label_distribution", "test_label_distribution"):
            distribution = fold.get(key, {})
            if isinstance(distribution, Mapping):
                labels.update(str(label) for label in distribution)
    return labels


def _required_columns(manifest: Mapping[str, Any]) -> list[str]:
    output = manifest.get("output", {})
    if not isinstance(output, Mapping):
        return []
    return [str(column) for column in output.get("required_columns", [])]


def validate_submission(
    manifest: Mapping[str, Any],
    submission_path: str | Path,
) -> SubmissionValidationReport:
    """Validate a prediction CSV against a task manifest."""

    path = Path(submission_path)
    errors: list[str] = []
    warnings: list[str] = []
    columns: list[str] = []
    rows: list[dict[str, str]] = []

    try:
        columns, rows = _read_submission(path)
    except FileNotFoundError:
        errors.append(f"submission file not found: {path}")
    except OSError as exc:
        errors.append(f"could not read submission file {path}: {exc}")

    task_id = str(manifest.get("task_id", ""))
    required = _required_columns(manifest)
    if not columns and not errors:
        errors.append("submission file has no header")

    missing = [column for column in required if column not in columns]
    if missing:
        errors.append(f"missing required columns: {', '.join(missing)}")

    if columns and not rows:
        errors.append("submission file has no prediction rows")

    if "task_id" in columns:
        mismatches = sorted(
            {
                row.get("task_id", "")
                for row in rows
                if row.get("task_id", "") and row.get("task_id", "") != task_id
            }
        )
        if mismatches:
            errors.append(
                "task_id column does not match manifest task_id "
                f"{task_id}: {', '.join(mismatches)}"
            )

    label_domain = _expected_labels(manifest)
    for label_column in ("true_label", "predicted_label"):
        if label_domain and label_column in columns:
            unknown = sorted(
                {
                    row.get(label_column, "")
                    for row in rows
                    if row.get(label_column, "") not in label_domain
                }
            )
            if unknown:
                errors.append(
                    f"{label_column} contains labels outside expected domain "
                    f"{sorted(label_domain)}: {', '.join(unknown)}"
                )

    embedding_columns = [column for column in columns if column.startswith("embedding_")]
    allowed_optional = {"task_id", "fold_id", "mission", "flight_probability"}
    extra_columns = [
        column
        for column in columns
        if column not in required
        and column not in allowed_optional
        and not column.startswith("embedding_")
    ]
    if extra_columns:
        warnings.append(f"unrecognized extra columns preserved: {', '.join(extra_columns)}")

    if "flight_probability" in columns:
        for row_index, row in enumerate(rows, start=2):
            value = row.get("flight_probability", "")
            if value == "":
                errors.append(f"flight_probability missing at CSV row {row_index}")
                continue
            try:
                probability = float(value)
            except ValueError:
                errors.append(f"flight_probability is not numeric at CSV row {row_index}")
                continue
            if probability < 0.0 or probability > 1.0:
                errors.append(
                    f"flight_probability outside [0, 1] at CSV row {row_index}"
                )

    for column in embedding_columns:
        for row_index, row in enumerate(rows, start=2):
            value = row.get(column, "")
            if value == "":
                errors.append(f"{column} missing at CSV row {row_index}")
                continue
            try:
                float(value)
            except ValueError:
                errors.append(f"{column} is not numeric at CSV row {row_index}")

    return SubmissionValidationReport(
        ok=not errors,
        submission_path=path.as_posix(),
        task_id=task_id,
        n_rows=len(rows),
        columns=columns,
        errors=errors,
        warnings=warnings,
        extra_columns=extra_columns,
        embedding_columns=embedding_columns,
    )

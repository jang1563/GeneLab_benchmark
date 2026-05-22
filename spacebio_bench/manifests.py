"""Manifest helpers for SpaceBio-Bench task contracts."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping


REQUIRED_TASK_FIELDS = (
    "schema_version",
    "task_id",
    "task_family",
    "title",
    "source_records",
    "split",
    "metrics",
    "output",
)

REQUIRED_SOURCE_FIELDS = (
    "source_id",
    "url_or_accession",
    "access_status",
    "checksum_status",
)

REQUIRED_SPLIT_FIELDS = ("name", "unit", "strategy")
REQUIRED_METRIC_FIELDS = ("metric_id", "profile", "interpretation")
REQUIRED_OUTPUT_FIELDS = ("prediction_format", "primary_artifacts")


def _require_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a mapping")
    return value


def _missing_fields(record: Mapping[str, Any], required: tuple[str, ...]) -> list[str]:
    return [field for field in required if field not in record]


def validate_task_manifest(manifest: Mapping[str, Any]) -> bool:
    """Validate the minimal v9 task-manifest contract.

    The validator is intentionally lightweight. It enforces fields that are
    critical for v9 design decisions without requiring a full JSON Schema
    dependency in the repository's test path.
    """

    _require_mapping(manifest, "manifest")
    missing = _missing_fields(manifest, REQUIRED_TASK_FIELDS)
    if missing:
        raise ValueError(f"task manifest missing required fields: {', '.join(missing)}")

    sources = manifest["source_records"]
    if not isinstance(sources, list) or not sources:
        raise ValueError("source_records must be a non-empty list")
    for index, source in enumerate(sources):
        source_map = _require_mapping(source, f"source_records[{index}]")
        missing = _missing_fields(source_map, REQUIRED_SOURCE_FIELDS)
        if missing:
            raise ValueError(
                f"source_records[{index}] missing required fields: {', '.join(missing)}"
            )

    split = _require_mapping(manifest["split"], "split")
    missing = _missing_fields(split, REQUIRED_SPLIT_FIELDS)
    if missing:
        raise ValueError(f"split missing required fields: {', '.join(missing)}")

    metrics = manifest["metrics"]
    if not isinstance(metrics, list) or not metrics:
        raise ValueError("metrics must be a non-empty list")
    for index, metric in enumerate(metrics):
        metric_map = _require_mapping(metric, f"metrics[{index}]")
        missing = _missing_fields(metric_map, REQUIRED_METRIC_FIELDS)
        if missing:
            raise ValueError(
                f"metrics[{index}] missing required fields: {', '.join(missing)}"
            )

    output = _require_mapping(manifest["output"], "output")
    missing = _missing_fields(output, REQUIRED_OUTPUT_FIELDS)
    if missing:
        raise ValueError(f"output missing required fields: {', '.join(missing)}")

    return True


def load_task_manifest(path: str | Path) -> Mapping[str, Any]:
    """Load and validate a task manifest from JSON."""

    manifest_path = Path(path)
    manifest = json.loads(manifest_path.read_text())
    validate_task_manifest(manifest)
    return manifest

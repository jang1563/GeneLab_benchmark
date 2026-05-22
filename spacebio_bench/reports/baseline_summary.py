"""Cross-baseline report aggregation helpers."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Mapping, Sequence


BASELINE_SUMMARY_FIELDS = [
    "baseline_id",
    "task_id",
    "status",
    "validation_ok",
    "n_predictions",
    "macro_f1",
    "balanced_accuracy",
    "auroc",
    "calibration_error",
    "mission_discrimination",
    "output_dir",
    "predictions",
    "metrics",
    "run_manifest",
    "error",
]


def read_baseline_summary(
    path: str | Path,
    *,
    baseline_id: str | None = None,
) -> list[dict[str, str]]:
    """Read one baseline summary CSV and normalize its columns."""

    with Path(path).open(newline="") as handle:
        rows = [dict(row) for row in csv.DictReader(handle)]

    normalized: list[dict[str, str]] = []
    for row in rows:
        output_row = {field: str(row.get(field, "") or "") for field in BASELINE_SUMMARY_FIELDS}
        if baseline_id is not None:
            output_row["baseline_id"] = baseline_id
        if not output_row["baseline_id"]:
            raise ValueError(f"{path} has no baseline_id column; pass baseline_id explicitly")
        normalized.append(output_row)
    return normalized


def write_baseline_summary(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> dict[str, Path]:
    """Write normalized cross-baseline CSV and JSON summaries."""

    out_csv = Path(csv_path)
    out_json = Path(json_path)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    sorted_rows = sorted(
        ({field: str(row.get(field, "") or "") for field in BASELINE_SUMMARY_FIELDS} for row in rows),
        key=lambda row: (row["baseline_id"], row["task_id"]),
    )
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=BASELINE_SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(sorted_rows)
    out_json.write_text(json.dumps(sorted_rows, indent=2, sort_keys=True) + "\n")
    return {"csv": out_csv, "json": out_json}

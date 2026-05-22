"""Task registry helpers for SpaceBio-Bench manifests."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Mapping

from .manifests import load_task_manifest


@dataclass(frozen=True)
class TaskRegistry:
    """In-memory registry of validated v9 task manifests."""

    manifests: Mapping[str, Mapping[str, Any]]

    @classmethod
    def from_dir(cls, manifest_dir: str | Path) -> "TaskRegistry":
        path = Path(manifest_dir)
        manifests: dict[str, Mapping[str, Any]] = {}
        for manifest_path in sorted(path.glob("*.json")):
            manifest = load_task_manifest(manifest_path)
            task_id = str(manifest["task_id"])
            if task_id in manifests:
                raise ValueError(f"duplicate task_id in registry: {task_id}")
            manifests[task_id] = manifest
        if not manifests:
            raise ValueError(f"no task manifests found in {path}")
        return cls(manifests=manifests)

    def __len__(self) -> int:
        return len(self.manifests)

    def __iter__(self) -> Iterator[Mapping[str, Any]]:
        for task_id in self.task_ids():
            yield self.manifests[task_id]

    def task_ids(self) -> list[str]:
        return sorted(self.manifests)

    def get(self, task_id: str) -> Mapping[str, Any]:
        try:
            return self.manifests[task_id]
        except KeyError as exc:
            raise KeyError(f"unknown task_id {task_id!r}") from exc

    def summary_rows(self) -> list[dict[str, str]]:
        return [summarize_task(manifest) for manifest in self]


def _join(values: list[Any]) -> str:
    return ";".join(str(value) for value in values)


def summarize_task(manifest: Mapping[str, Any]) -> dict[str, str]:
    """Return a compact index row for one task manifest."""

    split = manifest["split"]
    sources = manifest["source_records"]
    metrics = manifest["metrics"]
    return {
        "task_id": str(manifest["task_id"]),
        "legacy_task_id": str(manifest.get("legacy_task_id", "")),
        "task_family": str(manifest["task_family"]),
        "tissue": str(manifest.get("tissue", "")),
        "organism": str(manifest.get("organism", "")),
        "material_type": str(manifest.get("material_type", "")),
        "model_system": str(manifest.get("model_system", "")),
        "assay_modality": str(manifest.get("assay_modality", "")),
        "feature_namespace": str(manifest.get("feature_namespace", "")),
        "variant": str(manifest.get("variant", "")),
        "n_missions": str(split.get("n_missions", "")),
        "n_folds": str(split.get("n_folds", "")),
        "missions": _join(list(split.get("missions", []))),
        "n_sources": str(len(sources)),
        "source_ids": _join([source["source_id"] for source in sources]),
        "metric_ids": _join([metric["metric_id"] for metric in metrics]),
        "manifest_schema_version": str(manifest["schema_version"]),
    }


def write_task_index(
    registry: TaskRegistry,
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    """Write CSV and JSON task-index summaries."""

    rows = registry.summary_rows()
    output_csv = Path(csv_path)
    output_json = Path(json_path)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = list(rows[0])
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    output_json.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
    return output_csv, output_json

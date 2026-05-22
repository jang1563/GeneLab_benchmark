"""Draft Frictionless Data Package helpers for SpaceBio-Bench v9 artifacts."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Sequence


PACKAGE_CREATED_DATE = "2026-05-21"
CSV_RESOURCE_SPECS = [
    {
        "name": "task_manifest_index",
        "path": "v9/task_manifest_index.csv",
        "title": "Task Manifest Index",
        "description": "Compact index of generated v9 bulk LOMO task manifests.",
        "bundle_part": "metadata_spine",
        "primary_key": ["task_id"],
    },
    {
        "name": "task_data_index",
        "path": "v9/task_data_index.csv",
        "title": "Task Data Index",
        "description": (
            "Fold-level index of legacy processed matrices used by the current "
            "public bulk LOMO scaffold."
        ),
        "bundle_part": "data_index",
        "primary_key": ["task_id", "fold_id"],
    },
    {
        "name": "source_inventory",
        "path": "v9/source_inventory.csv",
        "title": "Source Inventory",
        "description": "Deduplicated OSDR source inventory for generated public bulk tasks.",
        "bundle_part": "metadata_spine",
        "primary_key": ["source_id"],
    },
    {
        "name": "source_checksum_audit",
        "path": "v9/source_checksum_audit.csv",
        "title": "Source Checksum Audit",
        "description": "OSDR API and checksum-manifest evidence for source inventory rows.",
        "bundle_part": "provenance_evidence",
        "primary_key": ["source_id"],
    },
    {
        "name": "bulk_lomo_baseline_summary",
        "path": "v9/reports/bulk_lomo_baseline_summary.csv",
        "title": "Bulk LOMO Baseline Summary",
        "description": "Cross-baseline summary for all current public bulk LOMO tasks.",
        "bundle_part": "benchmark_outputs",
        "primary_key": ["baseline_id", "task_id"],
    },
    {
        "name": "nearest_centroid_summary",
        "path": "v9/reports/nearest_centroid/bulk_lomo_summary.csv",
        "title": "Nearest-Centroid Summary",
        "description": "Nearest-centroid per-task summary table.",
        "bundle_part": "benchmark_outputs",
        "primary_key": ["task_id"],
    },
    {
        "name": "sklearn_baseline_summary",
        "path": "v9/reports/sklearn_baselines/bulk_lomo_summary.csv",
        "title": "Sklearn Baseline Summary",
        "description": "PCA-LR and L2 logistic-regression per-task summary table.",
        "bundle_part": "benchmark_outputs",
        "primary_key": ["baseline_id", "task_id"],
    },
]


def _repo_path(repo_root: str | Path, relative_path: str | Path) -> Path:
    return Path(repo_root) / Path(relative_path)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _infer_field_type(values: Iterable[str]) -> str:
    observed = [value for value in values if value not in {"", None}]
    if not observed:
        return "string"
    lowered = {str(value).lower() for value in observed}
    if lowered <= {"true", "false"}:
        return "boolean"
    try:
        for value in observed:
            int(str(value))
        return "integer"
    except ValueError:
        pass
    try:
        for value in observed:
            float(str(value))
        return "number"
    except ValueError:
        return "string"


def _csv_schema(path: Path) -> dict[str, Any]:
    rows = _read_csv_rows(path)
    if not rows:
        with path.open(newline="") as handle:
            fieldnames = csv.DictReader(handle).fieldnames or []
    else:
        fieldnames = list(rows[0])
    fields = [
        {
            "name": field,
            "type": _infer_field_type(row.get(field, "") for row in rows),
        }
        for field in fieldnames
    ]
    return {"fields": fields, "missingValues": [""]}


def _file_entry(repo_root: Path, path: Path) -> dict[str, Any]:
    relative_path = path.relative_to(repo_root).as_posix()
    return {
        "path": relative_path,
        "bytes": path.stat().st_size,
        "hash": _sha256_file(path),
        "hash_algorithm": "sha256",
    }


def _tabular_resource(
    repo_root: Path,
    *,
    name: str,
    relative_path: str,
    title: str,
    description: str,
    bundle_part: str,
    primary_key: Sequence[str],
) -> dict[str, Any]:
    path = _repo_path(repo_root, relative_path)
    resource = {
        "profile": "tabular-data-resource",
        "name": name,
        "title": title,
        "description": description,
        "path": relative_path,
        "format": "csv",
        "mediatype": "text/csv",
        "bytes": path.stat().st_size,
        "hash": _sha256_file(path),
        "hash_algorithm": "sha256",
        "schema": _csv_schema(path),
        "spacebio_bench:bundle_part": bundle_part,
    }
    if primary_key:
        resource["schema"]["primaryKey"] = list(primary_key)
    return resource


def _collection_resource(
    repo_root: Path,
    *,
    name: str,
    title: str,
    description: str,
    paths: Sequence[Path],
    mediatype: str,
    bundle_part: str,
    schema_from: Path | None = None,
) -> dict[str, Any]:
    entries = [_file_entry(repo_root, path) for path in paths]
    resource: dict[str, Any] = {
        "name": name,
        "title": title,
        "description": description,
        "path": [entry["path"] for entry in entries],
        "mediatype": mediatype,
        "bytes": sum(entry["bytes"] for entry in entries),
        "spacebio_bench:file_count": len(entries),
        "spacebio_bench:files": entries,
        "spacebio_bench:bundle_part": bundle_part,
    }
    if schema_from is not None:
        resource.update(
            {
                "profile": "tabular-data-resource",
                "format": "csv",
                "schema": _csv_schema(schema_from),
            }
        )
    return resource


def _glob_existing(repo_root: Path, pattern: str) -> list[Path]:
    return sorted(path for path in repo_root.glob(pattern) if path.is_file())


def build_v9_public_bulk_datapackage(
    *,
    repo_root: str | Path = ".",
    version: str = "0.1.0-draft",
) -> dict[str, Any]:
    """Build a draft Data Package descriptor for current v9 public bulk artifacts."""

    root = Path(repo_root).resolve()
    resources = [
        _tabular_resource(
            root,
            name=str(spec["name"]),
            relative_path=str(spec["path"]),
            title=str(spec["title"]),
            description=str(spec["description"]),
            bundle_part=str(spec["bundle_part"]),
            primary_key=spec["primary_key"],
        )
        for spec in CSV_RESOURCE_SPECS
    ]

    task_manifest_paths = _glob_existing(root, "v9/task_manifests/*.json")
    prediction_paths = _glob_existing(root, "v9/reports/**/predictions.csv")
    metrics_paths = _glob_existing(root, "v9/reports/**/metrics.json")
    run_manifest_paths = _glob_existing(root, "v9/reports/**/run_manifest.json")

    if task_manifest_paths:
        resources.append(
            _collection_resource(
                root,
                name="task_manifests",
                title="Task Manifests",
                description="Generated v9 bulk LOMO task manifest JSON files.",
                paths=task_manifest_paths,
                mediatype="application/json",
                bundle_part="metadata_spine",
            )
        )
    if prediction_paths:
        resources.append(
            _collection_resource(
                root,
                name="baseline_predictions",
                title="Baseline Predictions",
                description="Per-task baseline prediction CSV outputs.",
                paths=prediction_paths,
                mediatype="text/csv",
                bundle_part="benchmark_outputs",
                schema_from=prediction_paths[0],
            )
        )
    if metrics_paths:
        resources.append(
            _collection_resource(
                root,
                name="baseline_metrics",
                title="Baseline Metrics",
                description="Per-task baseline metrics JSON outputs.",
                paths=metrics_paths,
                mediatype="application/json",
                bundle_part="benchmark_outputs",
            )
        )
    if run_manifest_paths:
        resources.append(
            _collection_resource(
                root,
                name="baseline_run_manifests",
                title="Baseline Run Manifests",
                description="Per-task baseline run provenance JSON outputs.",
                paths=run_manifest_paths,
                mediatype="application/json",
                bundle_part="benchmark_outputs",
            )
        )

    return {
        "profile": "data-package",
        "name": "spacebio-bench-v9-public-bulk-draft",
        "title": "SpaceBio-Bench v9 Public Bulk Draft Package",
        "description": (
            "Draft machine-readable package descriptor for SpaceBio-Bench v9 "
            "public bulk LOMO metadata, provenance evidence, and baseline outputs."
        ),
        "version": version,
        "created": PACKAGE_CREATED_DATE,
        "licenses": [
            {
                "name": "not-yet-finalized",
                "title": "Draft: OSDR source reuse and benchmark artifact license pending release review",
            }
        ],
        "sources": [
            {
                "title": "NASA Open Science Data Repository Biological Data API",
                "path": "https://visualization.osdr.nasa.gov/biodata/api/",
            },
            {
                "title": "SpaceBio-Bench v9 source inventory",
                "path": "v9/source_inventory.csv",
            },
            {
                "title": "SpaceBio-Bench v9 source checksum audit",
                "path": "v9/source_checksum_audit.csv",
            },
        ],
        "keywords": [
            "space-biology",
            "benchmark",
            "gene-expression",
            "mission-held-out",
            "OSDR",
        ],
        "resources": resources,
        "spacebio_bench:release_status": "draft_not_frozen",
        "spacebio_bench:payload_verification_status": "checksum_manifests_parsed_payloads_not_hashed",
        "spacebio_bench:public_core_only": True,
    }


def write_datapackage(package: dict[str, Any], path: str | Path) -> Path:
    """Write a Data Package JSON descriptor."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(package, indent=2, sort_keys=True) + "\n")
    return output_path

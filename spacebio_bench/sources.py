"""Source inventory builders for SpaceBio-Bench manifests."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from .registry import TaskRegistry


OSDR_STUDY_URL = "https://osdr.nasa.gov/bio/repo/data/studies/"
SOURCE_INVENTORY_FIELDS = [
    "source_id",
    "glds_prefix",
    "osd_url",
    "url_or_accession",
    "mission",
    "tissue",
    "organism",
    "taxon_id",
    "species_common_name",
    "material_type",
    "model_system",
    "biospecimen_type",
    "assay_modality",
    "platform",
    "spaceflight_environment",
    "ground_control_type",
    "donor_or_strain_block",
    "orthology_strategy",
    "feature_namespace",
    "task_ids",
    "task_families",
    "variants",
    "access_status",
    "privacy_class",
    "checksum_status",
    "release_target",
    "notes",
]


def _add(target: dict[str, set[str]], key: str, value: Any) -> None:
    text = str(value or "").strip()
    if text:
        target[key].add(text)


def _join(values: set[str]) -> str:
    return ";".join(sorted(values))


def _source_url(source_id: str, source: Mapping[str, Any]) -> str:
    url = str(source.get("source_url") or "")
    if url:
        return url
    if source_id.startswith("OSD-"):
        return f"{OSDR_STUDY_URL}{source_id}"
    return ""


def _metadata_value(
    source: Mapping[str, Any],
    manifest: Mapping[str, Any],
    key: str,
) -> Any:
    return source.get(key) or manifest.get(key)


def build_source_inventory(
    registry: TaskRegistry,
    *,
    release_target: str = "v9_alpha_public_bulk_candidate",
) -> list[dict[str, str]]:
    """Build one source-level inventory row per public source accession."""

    grouped: dict[str, dict[str, set[str]]] = {}
    for manifest in registry:
        task_id = str(manifest["task_id"])
        for source in manifest.get("source_records", []):
            source_id = str(source.get("source_id") or source.get("url_or_accession") or "")
            if not source_id:
                raise ValueError(f"{task_id} has a source record without source_id")
            row = grouped.setdefault(source_id, {field: set() for field in SOURCE_INVENTORY_FIELDS})
            _add(row, "source_id", source_id)
            _add(row, "glds_prefix", source.get("glds_prefix"))
            _add(row, "osd_url", _source_url(source_id, source))
            _add(row, "url_or_accession", source.get("url_or_accession") or source_id)
            _add(row, "mission", source.get("mission"))
            _add(row, "tissue", manifest.get("tissue"))
            _add(row, "organism", _metadata_value(source, manifest, "organism"))
            _add(row, "taxon_id", _metadata_value(source, manifest, "taxon_id"))
            _add(
                row,
                "species_common_name",
                _metadata_value(source, manifest, "species_common_name"),
            )
            _add(row, "material_type", _metadata_value(source, manifest, "material_type"))
            _add(row, "model_system", _metadata_value(source, manifest, "model_system"))
            _add(
                row,
                "biospecimen_type",
                _metadata_value(source, manifest, "biospecimen_type"),
            )
            _add(
                row,
                "assay_modality",
                _metadata_value(source, manifest, "assay_modality"),
            )
            _add(row, "platform", _metadata_value(source, manifest, "platform"))
            _add(
                row,
                "spaceflight_environment",
                _metadata_value(source, manifest, "spaceflight_environment"),
            )
            _add(
                row,
                "ground_control_type",
                _metadata_value(source, manifest, "ground_control_type"),
            )
            _add(
                row,
                "donor_or_strain_block",
                _metadata_value(source, manifest, "donor_or_strain_block"),
            )
            _add(
                row,
                "orthology_strategy",
                _metadata_value(source, manifest, "orthology_strategy"),
            )
            _add(
                row,
                "feature_namespace",
                _metadata_value(source, manifest, "feature_namespace"),
            )
            _add(row, "task_ids", task_id)
            _add(row, "task_families", manifest.get("task_family"))
            _add(row, "variants", manifest.get("variant"))
            _add(row, "access_status", source.get("access_status"))
            _add(row, "privacy_class", source.get("privacy_class"))
            _add(row, "checksum_status", source.get("checksum_status"))
            _add(row, "release_target", release_target)
            _add(row, "notes", source.get("notes"))

    rows = [
        {field: _join(values[field]) for field in SOURCE_INVENTORY_FIELDS}
        for _, values in sorted(grouped.items())
    ]
    return rows


def write_source_inventory(
    rows: Sequence[Mapping[str, str]],
    *,
    csv_path: str | Path,
    json_path: str | Path,
) -> tuple[Path, Path]:
    """Write source inventory rows as CSV and JSON."""

    if not rows:
        raise ValueError("cannot write an empty source inventory")
    output_csv = Path(csv_path)
    output_json = Path(json_path)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_json.parent.mkdir(parents=True, exist_ok=True)

    normalized = [
        {field: str(row.get(field, "") or "") for field in SOURCE_INVENTORY_FIELDS}
        for row in rows
    ]
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SOURCE_INVENTORY_FIELDS)
        writer.writeheader()
        writer.writerows(normalized)
    output_json.write_text(json.dumps(normalized, indent=2, sort_keys=True) + "\n")
    return output_csv, output_json

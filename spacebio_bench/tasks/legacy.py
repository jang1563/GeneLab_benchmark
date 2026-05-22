"""Converters from existing GeneLab task metadata to v9 task manifests."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Mapping

from spacebio_bench.manifests import validate_task_manifest


OSDR_STUDY_URL = "https://osdr.nasa.gov/bio/repo/data/studies/"
MOUSE_BULK_CONTEXT = {
    "organism": "Mus musculus",
    "taxon_id": "10090",
    "species_common_name": "mouse",
    "material_type": "animal_tissue",
    "model_system": "rodent_spaceflight",
    "assay_modality": "bulk_rna_seq",
    "platform": "OSDR_GLbulkRNAseq_processed_counts",
    "spaceflight_environment": "LEO_or_ISS_spaceflight",
    "ground_control_type": "matched_ground_control",
    "orthology_strategy": "not_applicable_within_mouse_task",
    "feature_namespace": "mouse_gene",
}


def _load_download_manifest() -> Mapping[str, Mapping[str, Any]]:
    """Load the existing repository OSDR manifest without adding a new registry."""

    try:
        from scripts.fetch_osdr import DOWNLOAD_MANIFEST
    except ImportError as exc:
        raise RuntimeError("could not import scripts.fetch_osdr.DOWNLOAD_MANIFEST") from exc
    return DOWNLOAD_MANIFEST


def _task_variant(task_dir: Path, task_id: str) -> str:
    prefix = f"{task_id}_"
    name = task_dir.name
    if not name.startswith(prefix):
        return ""
    suffix = name[len(prefix) :]
    if suffix.endswith("_lomo"):
        return ""
    marker = "_lomo_"
    if marker in suffix:
        return suffix.split(marker, 1)[1]
    return ""


def _source_records_for_task(
    task_id: str,
    tissue: str,
    download_manifest: Mapping[str, Mapping[str, Any]],
    allowed_missions: set[str] | None = None,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for osd_id, info in sorted(download_manifest.items()):
        if info.get("task") != task_id or info.get("tissue") != tissue:
            continue
        source_mission = str(info.get("mission", ""))
        if allowed_missions is not None and not _mission_matches(source_mission, allowed_missions):
            continue
        glds_prefix = str(info.get("glds_prefix", ""))
        record = {
            "source_id": osd_id,
            "url_or_accession": osd_id,
            "access_status": "public",
            "checksum_status": "legacy_task_source_unfrozen",
            "privacy_class": "public",
            "source_url": f"{OSDR_STUDY_URL}{osd_id}",
            "glds_prefix": glds_prefix,
            "mission": source_mission,
            "notes": info.get("notes", ""),
            **MOUSE_BULK_CONTEXT,
            "biospecimen_type": tissue,
        }
        records.append(record)
    if not records:
        raise ValueError(f"no OSDR source records found for {task_id} / {tissue}")
    return records


def _mission_matches(source_mission: str, allowed_missions: set[str]) -> bool:
    if source_mission in allowed_missions:
        return True
    return any(
        source_mission.startswith(f"{mission} ")
        or source_mission.startswith(f"{mission}_")
        for mission in allowed_missions
    )


def _all_missions(folds: Iterable[Mapping[str, Any]]) -> list[str]:
    missions: set[str] = set()
    for fold in folds:
        missions.add(str(fold["test_mission"]))
        missions.update(str(mission) for mission in fold.get("train_missions", []))
    return sorted(missions)


def legacy_task_info_to_manifest(
    task_info: Mapping[str, Any],
    *,
    task_dir: str | Path,
    download_manifest: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Convert an existing `tasks/*/task_info.json` record to a v9 manifest."""

    task_path = Path(task_dir)
    legacy_task_id = str(task_info.get("task_id") or task_path.name.split("_", 1)[0])
    task_id = legacy_task_id.split("_", 1)[0]
    tissue = str(task_info["tissue"])
    folds = task_info.get("folds", [])
    if not isinstance(folds, list) or not folds:
        raise ValueError("legacy task_info must include a non-empty folds list")
    missions = _all_missions(folds)

    source_manifest = download_manifest or _load_download_manifest()
    variant = _task_variant(task_path, task_id)
    if not variant and legacy_task_id != task_id and legacy_task_id.startswith(f"{task_id}_"):
        variant = legacy_task_id[len(task_id) + 1 :]
    manifest_task_id = f"{task_id}_{tissue}_bulk_lomo"
    if variant:
        manifest_task_id = f"{manifest_task_id}_{variant}"

    manifest = {
        "schema_version": "0.1.0",
        "task_id": manifest_task_id,
        "legacy_task_id": legacy_task_id,
        "legacy_task_dir": task_path.as_posix(),
        "task_family": "bulk_lomo",
        "title": f"GeneLab {task_id} {tissue} bulk leave-one-mission-out",
        "tissue": tissue,
        "variant": variant or "canonical",
        **MOUSE_BULK_CONTEXT,
        "biospecimen_type": tissue,
        "source_records": _source_records_for_task(
            task_id,
            tissue,
            source_manifest,
            allowed_missions=set(missions),
        ),
        "split": {
            "name": "leave_one_mission_out",
            "unit": "mission",
            "strategy": "train on all missions except the held-out mission",
            "n_missions": task_info.get("n_missions"),
            "n_folds": task_info.get("n_folds_generated", len(folds)),
            "missions": missions,
            "folds": folds,
        },
        "metrics": [
            {
                "metric_id": "macro_f1",
                "profile": "genelab_minimal",
                "interpretation": "Primary class-balance-aware Flight vs Ground classification score.",
            },
            {
                "metric_id": "balanced_accuracy",
                "profile": "genelab_minimal",
                "interpretation": "Mean recall across Flight and Ground labels.",
            },
            {
                "metric_id": "auroc",
                "profile": "genelab_minimal",
                "interpretation": "Ranking quality for binary Flight vs Ground probabilities.",
            },
            {
                "metric_id": "mission_discrimination",
                "profile": "genelab_minimal",
                "interpretation": "Ranks whether representations are nearest to their own mission centroid.",
            },
        ],
        "output": {
            "prediction_format": "csv",
            "primary_artifacts": ["predictions.csv", "metrics.json"],
            "required_columns": ["sample_id", "true_label", "predicted_label"],
            "optional_columns": ["flight_probability", "embedding_*"],
        },
        "provenance": {
            "derived_from": task_path.joinpath("task_info.json").as_posix(),
            "legacy_generated_at": task_info.get("generated_at"),
            "legacy_split": task_info.get("split"),
            "dry_run": task_info.get("dry_run"),
        },
    }
    validate_task_manifest(manifest)
    return manifest


def load_legacy_task_info(path: str | Path) -> Mapping[str, Any]:
    """Load a legacy `task_info.json` file."""

    return json.loads(Path(path).read_text())


def write_manifest(manifest: Mapping[str, Any], path: str | Path) -> Path:
    """Write a v9 manifest with stable key order and ASCII-safe JSON."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return output_path

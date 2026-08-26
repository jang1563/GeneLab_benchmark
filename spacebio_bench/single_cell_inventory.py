"""RRRM single-cell asset inventory helpers for SpaceBio-Bench v9."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

from .profiles import METRIC_PROFILES


DEFAULT_SINGLE_CELL_ASSET_INVENTORY_ID = "v9_sc_001_rrrm_asset_inventory"

SINGLE_CELL_ASSET_INVENTORY_SUMMARY_FIELDS = [
    "inventory_id",
    "decision_status",
    "total_asset_count",
    "rrrm1_asset_count",
    "rrrm2_asset_count",
    "documentation_count",
    "script_count",
    "evaluation_count",
    "figure_count",
    "processed_summary_count",
    "generated_cache_count",
    "local_h5ad_count",
    "local_loom_count",
    "local_mtx_count",
    "local_anndata_payload_count",
    "v9_sc_manifest_count",
    "metric_profile_status",
    "next_required_block",
    "claim_boundary",
]

SINGLE_CELL_ASSET_INVENTORY_ROW_FIELDS = [
    "inventory_id",
    "asset_path",
    "asset_family",
    "asset_class",
    "file_type",
    "exists",
    "bytes",
    "sha256",
    "current_use",
    "promotion_status",
    "notes",
]

SINGLE_CELL_LOCAL_PAYLOAD_SCAN_FIELDS = [
    "inventory_id",
    "payload_path",
    "payload_type",
    "bytes",
    "sha256",
    "claim_boundary",
]


def _resolve_path(path: str | Path, repo_root: Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return repo_root / candidate


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _asset_family(relative_path: str) -> str:
    if relative_path.startswith("v2/"):
        return "RRRM-1"
    if relative_path.startswith("v3/"):
        return "RRRM-2"
    return "single_cell_related"


def _asset_class(relative_path: str) -> str:
    path = Path(relative_path)
    parts = set(path.parts)
    if "__pycache__" in parts:
        return "generated_cache"
    if "docs" in parts:
        return "documentation"
    if "scripts" in parts:
        return "script"
    if "evaluation" in parts:
        return "evaluation_config"
    if "figures" in parts:
        return "figure"
    if "processed" in parts:
        return "processed_summary"
    if path.suffix.lower() in {".h5ad", ".loom", ".mtx"}:
        return "local_anndata_payload"
    return "other"


def _asset_use_and_status(asset_class: str) -> tuple[str, str, str]:
    if asset_class == "documentation":
        return (
            "legacy evidence for RRRM scope, samples, QC, and task design",
            "review_for_v9_manifest",
            "safe to cite as inventory evidence, not as a v9 task manifest",
        )
    if asset_class == "script":
        return (
            "legacy conversion, annotation, benchmark, or HPC pipeline code",
            "inspect_before_reuse",
            "requires v9 API/contract review before promotion",
        )
    if asset_class == "evaluation_config":
        return (
            "legacy RRRM-2 benchmark result/config surface",
            "reference_only",
            "do not promote as v9 score without regenerated task contract",
        )
    if asset_class == "figure":
        return (
            "legacy result visualization",
            "reference_only",
            "do not promote as v9 result",
        )
    if asset_class == "processed_summary":
        return (
            "legacy processed summary table",
            "evidence_only_no_raw_payload",
            "not a local AnnData payload",
        )
    if asset_class == "generated_cache":
        return (
            "interpreter cache artifact",
            "exclude_from_v9",
            "cache files must not enter v9 artifacts",
        )
    if asset_class == "local_anndata_payload":
        return (
            "local single-cell matrix payload",
            "payload_requires_manifest_review",
            "must be checked for obs/var and provenance before v9 use",
        )
    return ("single-cell related artifact", "review_needed", "unclassified")


def _payload_scan_paths(repo_root: Path) -> list[Path]:
    excluded_dirs = {
        ".git",
        ".venv",
        "__pycache__",
        "submissions",
        "node_modules",
    }
    payloads: list[Path] = []
    for path in repo_root.rglob("*"):
        if any(part in excluded_dirs for part in path.parts):
            continue
        if path.is_file() and path.suffix.lower() in {".h5ad", ".loom", ".mtx"}:
            payloads.append(path)
    return sorted(payloads)


def _read_asset_paths(path: Path) -> list[str]:
    return [str(item) for item in json.loads(path.read_text())]


def _count(rows: Iterable[dict[str, str]], key: str, value: str) -> int:
    return sum(1 for row in rows if row[key] == value)


def build_single_cell_asset_inventory(
    *,
    repo_root: str | Path = ".",
    inventory_id: str = DEFAULT_SINGLE_CELL_ASSET_INVENTORY_ID,
    seed_asset_paths: str | Path = "v9/reports/recenter_decision/single_cell_asset_paths.json",
) -> dict[str, Any]:
    """Build a conservative RRRM single-cell asset inventory."""

    root = Path(repo_root).resolve()
    seed_path = _resolve_path(seed_asset_paths, root)
    relative_paths = _read_asset_paths(seed_path)
    asset_rows: list[dict[str, str]] = []
    for relative_path in relative_paths:
        path = _resolve_path(relative_path, root)
        exists = path.exists()
        asset_class = _asset_class(relative_path)
        current_use, promotion_status, notes = _asset_use_and_status(asset_class)
        asset_rows.append(
            {
                "inventory_id": inventory_id,
                "asset_path": relative_path,
                "asset_family": _asset_family(relative_path),
                "asset_class": asset_class,
                "file_type": path.suffix.lower().lstrip(".") or "directory",
                "exists": str(exists).lower(),
                "bytes": str(path.stat().st_size) if exists and path.is_file() else "",
                "sha256": _sha256_file(path) if exists and path.is_file() else "",
                "current_use": current_use,
                "promotion_status": promotion_status,
                "notes": notes,
            }
        )

    payload_paths = _payload_scan_paths(root)
    payload_rows = [
        {
            "inventory_id": inventory_id,
            "payload_path": path.relative_to(root).as_posix(),
            "payload_type": path.suffix.lower().lstrip("."),
            "bytes": str(path.stat().st_size),
            "sha256": _sha256_file(path),
            "claim_boundary": "local_payload_detected_requires_obs_var_review",
        }
        for path in payload_paths
    ]
    v9_sc_manifest_paths = sorted(
        (root / "v9" / "sc_spaceflight" / "task_manifests").glob("*.json")
        if (root / "v9" / "sc_spaceflight" / "task_manifests").exists()
        else []
    )
    h5ad_count = sum(1 for row in payload_rows if row["payload_type"] == "h5ad")
    loom_count = sum(1 for row in payload_rows if row["payload_type"] == "loom")
    mtx_count = sum(1 for row in payload_rows if row["payload_type"] == "mtx")
    local_payload_count = h5ad_count + loom_count + mtx_count
    metric_profile_status = (
        "genelab_sc_profile_present"
        if "genelab_sc" in METRIC_PROFILES
        else "genelab_sc_profile_missing"
    )
    claim_boundary = "single_cell_asset_inventory_only_no_v9_sc_task_or_payload_claim"
    summary = {
        "inventory_id": inventory_id,
        "decision_status": "legacy_rrrm_assets_indexed_no_local_anndata_payload",
        "total_asset_count": str(len(asset_rows)),
        "rrrm1_asset_count": str(_count(asset_rows, "asset_family", "RRRM-1")),
        "rrrm2_asset_count": str(_count(asset_rows, "asset_family", "RRRM-2")),
        "documentation_count": str(_count(asset_rows, "asset_class", "documentation")),
        "script_count": str(_count(asset_rows, "asset_class", "script")),
        "evaluation_count": str(_count(asset_rows, "asset_class", "evaluation_config")),
        "figure_count": str(_count(asset_rows, "asset_class", "figure")),
        "processed_summary_count": str(_count(asset_rows, "asset_class", "processed_summary")),
        "generated_cache_count": str(_count(asset_rows, "asset_class", "generated_cache")),
        "local_h5ad_count": str(h5ad_count),
        "local_loom_count": str(loom_count),
        "local_mtx_count": str(mtx_count),
        "local_anndata_payload_count": str(local_payload_count),
        "v9_sc_manifest_count": str(len(v9_sc_manifest_paths)),
        "metric_profile_status": metric_profile_status,
        "next_required_block": "V9-SC-002: AnnData task manifest draft",
        "claim_boundary": claim_boundary,
    }
    review_md = f"""# V9-SC-001 RRRM Asset Inventory

Status: `{summary["decision_status"]}`

Claim boundary: `{claim_boundary}`

## Inventory Summary

- Total legacy RRRM/single-cell asset paths: {len(asset_rows)}
- RRRM-1 asset paths: {summary["rrrm1_asset_count"]}
- RRRM-2 asset paths: {summary["rrrm2_asset_count"]}
- Documentation rows: {summary["documentation_count"]}
- Script rows: {summary["script_count"]}
- Evaluation/config rows: {summary["evaluation_count"]}
- Figure rows: {summary["figure_count"]}
- Processed-summary rows: {summary["processed_summary_count"]}
- Generated-cache rows excluded from promotion: {summary["generated_cache_count"]}
- Local `.h5ad` files found: {h5ad_count}
- Local `.loom` files found: {loom_count}
- Local `.mtx` files found: {mtx_count}
- v9 `sc_spaceflight` task manifests found: {len(v9_sc_manifest_paths)}
- Metric profile status: `{metric_profile_status}`

## Decision

The single-cell lane has enough legacy RRRM documentation, scripts, figures,
and evaluation surfaces to justify a v9 flagship manifest draft, but not enough
local matrix payload evidence to claim a runnable v9 single-cell task. This
inventory is therefore evidence only: it does not create a v9 `sc_spaceflight`
task manifest, does not claim local AnnData payload availability, and does not
promote legacy RRRM scores as v9 benchmark results.

## Strongest Candidate Surface

- RRRM-1 has compact sample/QC/benchmark-ready manifests in `v2/docs/`.
- RRRM-2 has legacy benchmark scripts and evaluation JSON surfaces under
  `v3/scripts/` and `v3/evaluation/`.
- The next manifest draft should choose one candidate task and explicitly list
  required AnnData `obs`, `var`, label, split, and provenance fields before any
  evaluator implementation.

## Next Block

Run `V9-SC-002: AnnData task manifest draft`. Start from this inventory and
draft one non-runnable `sc_spaceflight` manifest with explicit payload blockers
instead of trying to regenerate RRRM matrices in the same block.
"""
    return {
        "summary": [summary],
        "asset_rows": asset_rows,
        "payload_rows": payload_rows,
        "review_md": review_md,
    }


def write_single_cell_asset_inventory(
    *,
    output_dir: str | Path = "v9/sc_spaceflight",
    repo_root: str | Path = ".",
    inventory_id: str = DEFAULT_SINGLE_CELL_ASSET_INVENTORY_ID,
    seed_asset_paths: str | Path = "v9/reports/recenter_decision/single_cell_asset_paths.json",
) -> dict[str, Path]:
    """Write RRRM single-cell asset inventory artifacts."""

    root = Path(repo_root).resolve()
    package = build_single_cell_asset_inventory(
        repo_root=root,
        inventory_id=inventory_id,
        seed_asset_paths=seed_asset_paths,
    )
    outdir = _resolve_path(output_dir, root)
    outputs = {
        "summary_csv": outdir / "asset_inventory_summary.csv",
        "summary_json": outdir / "asset_inventory_summary.json",
        "asset_csv": outdir / "asset_inventory.csv",
        "asset_json": outdir / "asset_inventory.json",
        "payload_csv": outdir / "local_payload_scan.csv",
        "payload_json": outdir / "local_payload_scan.json",
        "review_md": outdir / "asset_inventory.md",
    }
    _write_csv(
        outputs["summary_csv"],
        package["summary"],
        SINGLE_CELL_ASSET_INVENTORY_SUMMARY_FIELDS,
    )
    _write_json(outputs["summary_json"], package["summary"])
    _write_csv(
        outputs["asset_csv"],
        package["asset_rows"],
        SINGLE_CELL_ASSET_INVENTORY_ROW_FIELDS,
    )
    _write_json(outputs["asset_json"], package["asset_rows"])
    _write_csv(
        outputs["payload_csv"],
        package["payload_rows"],
        SINGLE_CELL_LOCAL_PAYLOAD_SCAN_FIELDS,
    )
    _write_json(outputs["payload_json"], package["payload_rows"])
    outputs["review_md"].parent.mkdir(parents=True, exist_ok=True)
    outputs["review_md"].write_text(str(package["review_md"]))
    return outputs

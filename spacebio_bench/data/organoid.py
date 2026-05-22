"""Data adapters for draft human organoid extension tasks."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from spacebio_bench.manifests import load_task_manifest


DEFAULT_ORGANOID_MANIFEST = (
    "v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json"
)


@dataclass(frozen=True)
class OrganoidFoldData:
    """Sample-level split for one draft organoid blocking factor."""

    task_id: str
    fold_id: str
    heldout_factor: str
    heldout_value: str
    train_sample_ids: list[str]
    test_sample_ids: list[str]
    train_label_distribution: Mapping[str, int]
    test_label_distribution: Mapping[str, int]
    status: str


@dataclass(frozen=True)
class HumanOrganoidTaskData:
    """Loaded draft human organoid task backed by normalized count matrices."""

    manifest: Mapping[str, Any]
    features: pd.DataFrame
    sample_factors: list[dict[str, str]]
    folds: list[OrganoidFoldData]
    diagnostic_folds: list[OrganoidFoldData]
    matrix_paths: Mapping[str, Path]
    feature_namespace: str

    @property
    def task_id(self) -> str:
        return str(self.manifest["task_id"])

    @property
    def n_samples(self) -> int:
        return int(self.features.shape[0])

    @property
    def n_features(self) -> int:
        return int(self.features.shape[1])


def _relative_or_absolute(path: str | Path, repo_root: Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return repo_root / candidate


def _read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _label_distribution(rows: list[dict[str, str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        label = str(row.get("true_label", "") or "")
        if not label:
            continue
        counts[label] = counts.get(label, 0) + 1
    return dict(sorted(counts.items()))


def _source_missions(manifest: Mapping[str, Any]) -> dict[str, str]:
    missions: dict[str, str] = {}
    for source in manifest.get("source_records", []):
        if not isinstance(source, Mapping):
            continue
        source_id = str(source.get("source_id", "") or "")
        if source_id:
            missions[source_id] = str(source.get("mission", "") or "")
    return missions


def _read_source_matrix(path: Path, sample_ids: list[str]) -> pd.DataFrame:
    matrix = pd.read_csv(path, index_col=0)
    if matrix.empty:
        raise ValueError(f"expression matrix is empty: {path}")
    matrix.index = matrix.index.map(str)
    sample_frame = matrix.T
    sample_frame.index = sample_frame.index.map(str)
    missing = sorted(set(sample_ids) - set(sample_frame.index))
    if missing:
        preview = ", ".join(missing[:5])
        raise ValueError(f"{path} is missing expected sample columns: {preview}")
    return sample_frame.loc[sample_ids].astype(float)


def _load_feature_matrix(
    *,
    manifest: Mapping[str, Any],
    sample_factors: list[dict[str, str]],
    repo_root: Path,
) -> tuple[pd.DataFrame, dict[str, Path]]:
    split = manifest.get("split", {})
    if not isinstance(split, Mapping):
        raise ValueError("human organoid manifest is missing split metadata")
    matrix_sources = split.get("expression_matrix_sources", {})
    if not isinstance(matrix_sources, Mapping) or not matrix_sources:
        raise ValueError("human organoid manifest is missing expression matrix sources")

    by_source: dict[str, list[str]] = {}
    for row in sample_factors:
        if row.get("parse_status") != "parsed":
            continue
        source_id = str(row.get("source_id", "") or "")
        sample_id = str(row.get("sample_id", "") or "")
        if source_id and sample_id:
            by_source.setdefault(source_id, []).append(sample_id)

    frames: list[pd.DataFrame] = []
    matrix_paths: dict[str, Path] = {}
    for source_id, source in sorted(matrix_sources.items()):
        if not isinstance(source, Mapping):
            continue
        sample_ids = by_source.get(str(source_id), [])
        if not sample_ids:
            raise ValueError(f"no parsed sample factors for expression source {source_id}")
        matrix_path = _relative_or_absolute(str(source["local_matrix_path"]), repo_root)
        if not matrix_path.exists():
            raise FileNotFoundError(f"expression matrix not found: {matrix_path}")
        matrix_paths[str(source_id)] = matrix_path
        frames.append(_read_source_matrix(matrix_path, sample_ids))

    if not frames:
        raise ValueError("no human organoid expression matrices were loaded")
    features = pd.concat(frames, axis=0, join="inner")
    if features.shape[1] == 0:
        raise ValueError("human organoid matrices have no common feature columns")

    sample_order = [
        str(row["sample_id"])
        for row in sample_factors
        if row.get("parse_status") == "parsed" and row.get("sample_id")
    ]
    missing_ordered = sorted(set(sample_order) - set(features.index))
    if missing_ordered:
        preview = ", ".join(missing_ordered[:5])
        raise ValueError(f"combined expression matrix is missing samples: {preview}")
    return features.loc[sample_order].astype(float), matrix_paths


def _build_folds(
    *,
    manifest: Mapping[str, Any],
    sample_factors: list[dict[str, str]],
    candidate_key: str = "candidate_folds",
    allowed_statuses: set[str] | None = None,
) -> list[OrganoidFoldData]:
    parsed_rows = [row for row in sample_factors if row.get("parse_status") == "parsed"]
    folds: list[OrganoidFoldData] = []
    split = manifest.get("split", {})
    candidates = split.get(candidate_key, []) if isinstance(split, Mapping) else []
    statuses = allowed_statuses or {"sample_count_backed_draft"}
    for candidate in candidates:
        if not isinstance(candidate, Mapping):
            continue
        if str(candidate.get("status", "") or "") not in statuses:
            continue
        factor = str(candidate.get("heldout_factor", "") or "")
        value = str(candidate.get("heldout_value", "") or "")
        if not factor or not value:
            continue
        test_rows = [row for row in parsed_rows if row.get(factor) == value]
        train_rows = [row for row in parsed_rows if row.get(factor) != value]
        if len(train_rows) != int(candidate.get("n_train", len(train_rows))):
            raise ValueError(f"{candidate.get('fold_id')} train sample count mismatch")
        if len(test_rows) != int(candidate.get("n_test", len(test_rows))):
            raise ValueError(f"{candidate.get('fold_id')} test sample count mismatch")
        folds.append(
            OrganoidFoldData(
                task_id=str(manifest["task_id"]),
                fold_id=str(candidate["fold_id"]),
                heldout_factor=factor,
                heldout_value=value,
                train_sample_ids=[str(row["sample_id"]) for row in train_rows],
                test_sample_ids=[str(row["sample_id"]) for row in test_rows],
                train_label_distribution=_label_distribution(train_rows),
                test_label_distribution=_label_distribution(test_rows),
                status=str(candidate["status"]),
            )
        )
    if not folds:
        raise ValueError(f"human organoid manifest has no supported folds in {candidate_key}")
    return folds


def load_human_organoid_task(
    *,
    manifest_path: str | Path = DEFAULT_ORGANOID_MANIFEST,
    repo_root: str | Path = ".",
) -> HumanOrganoidTaskData:
    """Load the draft human organoid task and aligned normalized-count features."""

    root = Path(repo_root)
    resolved_manifest = _relative_or_absolute(manifest_path, root)
    manifest = load_task_manifest(resolved_manifest)
    if manifest.get("task_family") != "human_organoid_spaceflight":
        raise ValueError(f"{resolved_manifest} is not a human_organoid_spaceflight task")

    split = manifest.get("split", {})
    if not isinstance(split, Mapping):
        raise ValueError("human organoid manifest is missing split metadata")
    sample_factor_path = _relative_or_absolute(
        str(split.get("sample_factor_table", "v9/human_organoid/sample_factors.draft.csv")),
        root,
    )
    sample_factors = _read_csv_rows(sample_factor_path)
    mission_by_source = _source_missions(manifest)
    for row in sample_factors:
        row["mission"] = mission_by_source.get(row.get("source_id", ""), "")

    features, matrix_paths = _load_feature_matrix(
        manifest=manifest,
        sample_factors=sample_factors,
        repo_root=root,
    )
    folds = _build_folds(manifest=manifest, sample_factors=sample_factors)
    diagnostic_folds = _build_folds(
        manifest=manifest,
        sample_factors=sample_factors,
        candidate_key="donor_diagnostic_folds",
        allowed_statuses={"metadata_backed_pilot_only_not_default"},
    )
    return HumanOrganoidTaskData(
        manifest=manifest,
        features=features,
        sample_factors=sample_factors,
        folds=folds,
        diagnostic_folds=diagnostic_folds,
        matrix_paths=matrix_paths,
        feature_namespace=str(manifest.get("feature_namespace", "")),
    )

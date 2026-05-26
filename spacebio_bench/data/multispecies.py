"""Data adapters for draft multispecies extension tasks."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import pandas as pd

from spacebio_bench.manifests import load_task_manifest


DEFAULT_MULTISPECIES_MANIFEST_DIR = "v9/multispecies/task_manifests"
DEFAULT_MULTISPECIES_INTERACTION_MANIFEST_DIR = (
    "v9/multispecies/interaction_task_manifests"
)


@dataclass(frozen=True)
class MultispeciesFoldData:
    """Sample-level split for one multispecies condition-stratum fold."""

    task_id: str
    fold_id: str
    fold_family: str
    heldout_factor: str
    heldout_value: str
    train_sample_ids: list[str]
    test_sample_ids: list[str]
    train_label_distribution: Mapping[str, int]
    test_label_distribution: Mapping[str, int]
    status: str

    @property
    def n_train(self) -> int:
        return len(self.train_sample_ids)

    @property
    def n_test(self) -> int:
        return len(self.test_sample_ids)


@dataclass(frozen=True)
class MultispeciesTaskData:
    """Loaded draft multispecies task backed by a local normalized-count matrix."""

    manifest: Mapping[str, Any]
    features: pd.DataFrame
    sample_factors: list[dict[str, str]]
    folds: list[MultispeciesFoldData]
    fold_families: Mapping[str, list[MultispeciesFoldData]]
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
        if label:
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


def _manifest_source_ids(manifest: Mapping[str, Any]) -> set[str]:
    source_ids: set[str] = set()
    for source in manifest.get("source_records", []):
        if isinstance(source, Mapping) and source.get("source_id"):
            source_ids.add(str(source["source_id"]))
    if not source_ids:
        raise ValueError("multispecies manifest has no source ids")
    return source_ids


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
        raise ValueError("multispecies manifest is missing split metadata")
    matrix_sources = split.get("expression_matrix_sources", {})
    if not isinstance(matrix_sources, Mapping) or not matrix_sources:
        raise ValueError("multispecies manifest is missing expression matrix sources")

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
        raise ValueError("no multispecies expression matrices were loaded")
    features = pd.concat(frames, axis=0, join="inner")
    if features.shape[1] == 0:
        raise ValueError("multispecies matrices have no common feature columns")

    source_ids = set(matrix_sources)
    sample_order = [
        str(row["sample_id"])
        for row in sample_factors
        if row.get("parse_status") == "parsed"
        and row.get("sample_id")
        and row.get("source_id") in source_ids
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
) -> list[MultispeciesFoldData]:
    split = manifest.get("split", {})
    candidates = split.get("candidate_folds", []) if isinstance(split, Mapping) else []
    return _build_folds_from_candidates(
        manifest=manifest,
        sample_factors=sample_factors,
        candidates=candidates,
        supported_statuses={"sample_count_backed_candidate_not_default"},
        fold_family="condition_stratum_candidate_folds",
    )


def _build_folds_from_candidates(
    *,
    manifest: Mapping[str, Any],
    sample_factors: list[dict[str, str]],
    candidates: Any,
    supported_statuses: set[str],
    fold_family: str,
) -> list[MultispeciesFoldData]:
    parsed_rows = [row for row in sample_factors if row.get("parse_status") == "parsed"]
    folds: list[MultispeciesFoldData] = []
    for candidate in candidates:
        if not isinstance(candidate, Mapping):
            continue
        if str(candidate.get("status", "") or "") not in supported_statuses:
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
            MultispeciesFoldData(
                task_id=str(manifest["task_id"]),
                fold_id=str(candidate["fold_id"]),
                fold_family=fold_family,
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
        raise ValueError(f"multispecies manifest has no supported {fold_family} folds")
    return folds


def load_multispecies_task(
    manifest_path: str | Path,
    *,
    repo_root: str | Path = ".",
) -> MultispeciesTaskData:
    """Load one draft multispecies species-native task."""

    root = Path(repo_root)
    resolved_manifest = _relative_or_absolute(manifest_path, root)
    manifest = load_task_manifest(resolved_manifest)
    if manifest.get("task_family") != "multispecies_species_native_spaceflight":
        raise ValueError(f"{resolved_manifest} is not a multispecies species-native task")

    split = manifest.get("split", {})
    if not isinstance(split, Mapping):
        raise ValueError("multispecies manifest is missing split metadata")
    sample_factor_path = _relative_or_absolute(
        str(split.get("sample_factor_table", "v9/multispecies/sample_factors.draft.csv")),
        root,
    )
    source_ids = _manifest_source_ids(manifest)
    all_sample_factors = _read_csv_rows(sample_factor_path)
    sample_factors = [
        row
        for row in all_sample_factors
        if str(row.get("source_id", "") or "") in source_ids
    ]
    mission_by_source = _source_missions(manifest)
    for row in sample_factors:
        row["mission"] = mission_by_source.get(row.get("source_id", ""), "")

    features, matrix_paths = _load_feature_matrix(
        manifest=manifest,
        sample_factors=sample_factors,
        repo_root=root,
    )
    folds = _build_folds(manifest=manifest, sample_factors=sample_factors)
    return MultispeciesTaskData(
        manifest=manifest,
        features=features,
        sample_factors=sample_factors,
        folds=folds,
        fold_families={"condition_stratum_candidate_folds": folds},
        matrix_paths=matrix_paths,
        feature_namespace=str(manifest.get("feature_namespace", "")),
    )


def _build_interaction_fold_families(
    *,
    manifest: Mapping[str, Any],
    sample_factors: list[dict[str, str]],
) -> dict[str, list[MultispeciesFoldData]]:
    split = manifest.get("split", {})
    if not isinstance(split, Mapping):
        raise ValueError("multispecies interaction manifest is missing split metadata")
    return {
        "primary_genotype_or_ecotype_holdout": _build_folds_from_candidates(
            manifest=manifest,
            sample_factors=sample_factors,
            candidates=split.get("primary_candidate_folds", []),
            supported_statuses={"sample_count_backed_primary_candidate"},
            fold_family="primary_genotype_or_ecotype_holdout",
        ),
        "secondary_light_treatment_holdout": _build_folds_from_candidates(
            manifest=manifest,
            sample_factors=sample_factors,
            candidates=split.get("secondary_light_treatment_folds", []),
            supported_statuses={"sample_count_backed_secondary_candidate"},
            fold_family="secondary_light_treatment_holdout",
        ),
        "diagnostic_condition_stratum_holdout": _build_folds_from_candidates(
            manifest=manifest,
            sample_factors=sample_factors,
            candidates=split.get("condition_stratum_diagnostic_folds", []),
            supported_statuses={"sample_count_backed_diagnostic_candidate"},
            fold_family="diagnostic_condition_stratum_holdout",
        ),
    }


def load_multispecies_interaction_task(
    manifest_path: str | Path,
    *,
    repo_root: str | Path = ".",
) -> MultispeciesTaskData:
    """Load one draft OSD-120 multispecies interaction task."""

    root = Path(repo_root)
    resolved_manifest = _relative_or_absolute(manifest_path, root)
    manifest = load_task_manifest(resolved_manifest)
    if manifest.get("task_family") != "multispecies_light_interaction_spaceflight":
        raise ValueError(f"{resolved_manifest} is not a multispecies interaction task")

    split = manifest.get("split", {})
    if not isinstance(split, Mapping):
        raise ValueError("multispecies interaction manifest is missing split metadata")
    sample_factor_path = _relative_or_absolute(
        str(split.get("sample_factor_table", "v9/multispecies/sample_factors.draft.csv")),
        root,
    )
    source_ids = _manifest_source_ids(manifest)
    all_sample_factors = _read_csv_rows(sample_factor_path)
    sample_factors = [
        row
        for row in all_sample_factors
        if str(row.get("source_id", "") or "") in source_ids
    ]
    mission_by_source = _source_missions(manifest)
    for row in sample_factors:
        row["mission"] = mission_by_source.get(row.get("source_id", ""), "")

    features, matrix_paths = _load_feature_matrix(
        manifest=manifest,
        sample_factors=sample_factors,
        repo_root=root,
    )
    fold_families = _build_interaction_fold_families(
        manifest=manifest,
        sample_factors=sample_factors,
    )
    primary_folds = fold_families["primary_genotype_or_ecotype_holdout"]
    return MultispeciesTaskData(
        manifest=manifest,
        features=features,
        sample_factors=sample_factors,
        folds=primary_folds,
        fold_families=fold_families,
        matrix_paths=matrix_paths,
        feature_namespace=str(manifest.get("feature_namespace", "")),
    )


def load_all_multispecies_tasks(
    *,
    manifest_dir: str | Path = DEFAULT_MULTISPECIES_MANIFEST_DIR,
    repo_root: str | Path = ".",
) -> list[MultispeciesTaskData]:
    """Load all draft multispecies species-native manifests from a directory."""

    root = Path(repo_root)
    resolved_dir = _relative_or_absolute(manifest_dir, root)
    manifests = sorted(resolved_dir.glob("*.json"))
    if not manifests:
        raise ValueError(f"no multispecies task manifests found in {resolved_dir}")
    return [
        load_multispecies_task(manifest_path=manifest_path, repo_root=root)
        for manifest_path in manifests
    ]


def load_all_multispecies_interaction_tasks(
    *,
    manifest_dir: str | Path = DEFAULT_MULTISPECIES_INTERACTION_MANIFEST_DIR,
    repo_root: str | Path = ".",
) -> list[MultispeciesTaskData]:
    """Load all draft multispecies interaction manifests from a directory."""

    root = Path(repo_root)
    resolved_dir = _relative_or_absolute(manifest_dir, root)
    manifests = sorted(resolved_dir.glob("*.json"))
    if not manifests:
        raise ValueError(f"no multispecies interaction task manifests found in {resolved_dir}")
    return [
        load_multispecies_interaction_task(manifest_path=manifest_path, repo_root=root)
        for manifest_path in manifests
    ]

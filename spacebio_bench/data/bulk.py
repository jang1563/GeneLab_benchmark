"""Read-only adapters for v9 bulk LOMO task data."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from spacebio_bench.registry import TaskRegistry


REQUIRED_FOLD_FILES = (
    "train_X.csv",
    "train_y.csv",
    "train_meta.csv",
    "test_X.csv",
    "test_y.csv",
    "test_meta.csv",
    "selected_genes.txt",
    "fold_info.json",
)


@dataclass(frozen=True)
class BulkFoldData:
    """File-level adapter for one legacy bulk LOMO fold."""

    task_id: str
    fold_id: str
    test_mission: str
    train_missions: list[str]
    fold_dir: Path
    paths: Mapping[str, Path]
    fold_info: Mapping[str, Any]
    selected_gene_count: int
    train_row_count: int
    test_row_count: int

    def to_summary_row(self) -> dict[str, str]:
        return {
            "task_id": self.task_id,
            "fold_id": self.fold_id,
            "test_mission": self.test_mission,
            "train_missions": ";".join(self.train_missions),
            "n_train": str(self.fold_info.get("n_train", "")),
            "n_test": str(self.fold_info.get("n_test", "")),
            "train_y_rows": str(self.train_row_count),
            "test_y_rows": str(self.test_row_count),
            "selected_gene_count": str(self.selected_gene_count),
            "fold_dir": self.fold_dir.as_posix(),
            "train_X": self.paths["train_X.csv"].as_posix(),
            "test_X": self.paths["test_X.csv"].as_posix(),
        }


@dataclass(frozen=True)
class BulkTaskData:
    """Read-only view of one v9 bulk task backed by legacy fold files."""

    manifest: Mapping[str, Any]
    task_dir: Path
    folds: list[BulkFoldData]

    @property
    def task_id(self) -> str:
        return str(self.manifest["task_id"])

    def to_summary_rows(self) -> list[dict[str, str]]:
        return [fold.to_summary_row() for fold in self.folds]


def _count_csv_rows(path: Path) -> int:
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        try:
            next(reader)
        except StopIteration:
            return 0
        return sum(1 for _ in reader)


def _count_lines(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for line in handle if line.strip())


def _fold_dir_name(test_mission: str) -> str:
    return f"fold_{test_mission}_test"


def _relative_or_absolute(path: str | Path, repo_root: Path) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    return repo_root / candidate


def _load_fold(
    *,
    manifest: Mapping[str, Any],
    task_dir: Path,
    fold_info_from_manifest: Mapping[str, Any],
) -> BulkFoldData:
    test_mission = str(fold_info_from_manifest["test_mission"])
    fold_dir = task_dir / _fold_dir_name(test_mission)
    paths = {filename: fold_dir / filename for filename in REQUIRED_FOLD_FILES}
    missing = [path.as_posix() for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            f"missing required files for {manifest['task_id']} / {test_mission}: "
            + ", ".join(missing)
        )

    fold_info = json.loads(paths["fold_info.json"].read_text())
    selected_gene_count = _count_lines(paths["selected_genes.txt"])
    train_row_count = _count_csv_rows(paths["train_y.csv"])
    test_row_count = _count_csv_rows(paths["test_y.csv"])

    expected_train = int(fold_info.get("n_train", train_row_count))
    expected_test = int(fold_info.get("n_test", test_row_count))
    if train_row_count != expected_train:
        raise ValueError(
            f"{paths['train_y.csv']} has {train_row_count} rows, expected {expected_train}"
        )
    if test_row_count != expected_test:
        raise ValueError(
            f"{paths['test_y.csv']} has {test_row_count} rows, expected {expected_test}"
        )

    expected_genes = fold_info.get("n_genes_after_var_filter")
    if expected_genes is not None and selected_gene_count != int(expected_genes):
        raise ValueError(
            f"{paths['selected_genes.txt']} has {selected_gene_count} genes, "
            f"expected {expected_genes}"
        )

    return BulkFoldData(
        task_id=str(manifest["task_id"]),
        fold_id=_fold_dir_name(test_mission),
        test_mission=test_mission,
        train_missions=[str(mission) for mission in fold_info.get("train_missions", [])],
        fold_dir=fold_dir,
        paths=paths,
        fold_info=fold_info,
        selected_gene_count=selected_gene_count,
        train_row_count=train_row_count,
        test_row_count=test_row_count,
    )


def load_bulk_task(
    task_id: str,
    *,
    manifest_dir: str | Path = "v9/task_manifests",
    repo_root: str | Path = ".",
) -> BulkTaskData:
    """Load a v9 bulk task as paths and validated fold metadata."""

    root = Path(repo_root)
    registry = TaskRegistry.from_dir(_relative_or_absolute(manifest_dir, root))
    manifest = registry.get(task_id)
    if manifest.get("task_family") != "bulk_lomo":
        raise ValueError(f"{task_id} is not a bulk_lomo task")

    task_dir = _relative_or_absolute(str(manifest["legacy_task_dir"]), root)
    if not task_dir.exists():
        raise FileNotFoundError(f"legacy task directory not found: {task_dir}")

    folds = [
        _load_fold(
            manifest=manifest,
            task_dir=task_dir,
            fold_info_from_manifest=fold_info,
        )
        for fold_info in manifest["split"]["folds"]
    ]
    return BulkTaskData(manifest=manifest, task_dir=task_dir, folds=folds)


def load_all_bulk_tasks(
    *,
    manifest_dir: str | Path = "v9/task_manifests",
    repo_root: str | Path = ".",
) -> list[BulkTaskData]:
    """Load all bulk LOMO task manifests in a directory."""

    root = Path(repo_root)
    registry = TaskRegistry.from_dir(_relative_or_absolute(manifest_dir, root))
    tasks: list[BulkTaskData] = []
    for task_id in registry.task_ids():
        manifest = registry.get(task_id)
        if manifest.get("task_family") == "bulk_lomo":
            tasks.append(
                load_bulk_task(
                    task_id,
                    manifest_dir=manifest_dir,
                    repo_root=root,
                )
            )
    return tasks


def bulk_task_summary_rows(tasks: list[BulkTaskData]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for task in tasks:
        rows.extend(task.to_summary_rows())
    return rows

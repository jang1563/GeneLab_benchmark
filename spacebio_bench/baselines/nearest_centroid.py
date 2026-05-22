"""Nearest-centroid baseline for v9 bulk leave-one-mission-out tasks."""

from __future__ import annotations

import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from spacebio_bench.data import BulkFoldData, BulkTaskData, load_all_bulk_tasks, load_bulk_task
from spacebio_bench.evaluate import evaluate_submission
from spacebio_bench.reports import write_evaluation_report


BASELINE_ID = "nearest_centroid"
NEGATIVE_LABEL = "Ground"
POSITIVE_LABEL = "Flight"
PREDICTION_COLUMNS = [
    "task_id",
    "fold_id",
    "sample_id",
    "mission",
    "true_label",
    "predicted_label",
    "flight_probability",
    "embedding_0",
    "embedding_1",
]


@dataclass(frozen=True)
class NearestCentroidConfig:
    """Configuration for the nearest-centroid reference baseline."""

    scaling: str = "none"

    def __post_init__(self) -> None:
        valid = {"none", "zscore"}
        if self.scaling not in valid:
            raise ValueError(f"scaling must be one of {sorted(valid)}, got {self.scaling!r}")


@dataclass(frozen=True)
class FoldPredictionResult:
    """Predictions and fold-level metadata for one held-out mission."""

    task_id: str
    fold_id: str
    test_mission: str
    n_train: int
    n_test: int
    predictions: list[dict[str, str]]


@dataclass(frozen=True)
class NearestCentroidTaskResult:
    """Filesystem outputs and metrics for one baseline task run."""

    task_id: str
    output_dir: Path
    predictions_path: Path
    metrics_path: Path
    run_manifest_path: Path
    evaluation: Mapping[str, Any]
    n_predictions: int


def _label_from_binary(value: object) -> str:
    if isinstance(value, str):
        normalized = value.strip().lower()
        if normalized in {"1", "1.0", "flight"}:
            return POSITIVE_LABEL
        if normalized in {"0", "0.0", "ground"}:
            return NEGATIVE_LABEL
    numeric = float(value)
    if numeric == 1.0:
        return POSITIVE_LABEL
    if numeric == 0.0:
        return NEGATIVE_LABEL
    raise ValueError(f"expected a binary label encoded as 0/1, got {value!r}")


def _read_features(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    if frame.empty:
        raise ValueError(f"feature matrix is empty: {path}")
    return frame.astype(float)


def _read_labels(path: Path, sample_index: pd.Index) -> list[str]:
    frame = pd.read_csv(path, index_col=0)
    if frame.shape[1] < 1:
        raise ValueError(f"label file has no label column: {path}")
    labels = frame.iloc[:, 0].map(_label_from_binary)
    missing = sample_index.difference(labels.index)
    if len(missing) > 0:
        preview = ", ".join(str(value) for value in missing[:5])
        raise ValueError(f"{path} is missing labels for samples: {preview}")
    return [str(value) for value in labels.reindex(sample_index)]


def _read_missions(path: Path, sample_index: pd.Index, fallback: str) -> list[str]:
    try:
        frame = pd.read_csv(path, index_col=0)
    except FileNotFoundError:
        return [fallback] * len(sample_index)
    if "mission" not in frame.columns:
        return [fallback] * len(sample_index)
    missions = frame["mission"].reindex(sample_index).fillna(fallback)
    return [str(value) for value in missions]


def _align_test_features(train: pd.DataFrame, test: pd.DataFrame, test_path: Path) -> pd.DataFrame:
    missing = train.columns.difference(test.columns)
    if len(missing) > 0:
        preview = ", ".join(str(value) for value in missing[:5])
        raise ValueError(f"{test_path} is missing train feature columns: {preview}")
    return test.loc[:, train.columns].astype(float)


def _scale_train_test(
    train: np.ndarray,
    test: np.ndarray,
    *,
    scaling: str,
) -> tuple[np.ndarray, np.ndarray]:
    if scaling == "none":
        return train, test
    if scaling != "zscore":
        raise ValueError(f"unsupported scaling mode: {scaling}")
    mean = train.mean(axis=0)
    std = train.std(axis=0)
    std[std == 0.0] = 1.0
    return (train - mean) / std, (test - mean) / std


def _class_centroids(features: np.ndarray, labels: Sequence[str]) -> np.ndarray:
    centroids = []
    for label in (NEGATIVE_LABEL, POSITIVE_LABEL):
        mask = np.asarray([value == label for value in labels], dtype=bool)
        if not mask.any():
            raise ValueError(f"training fold has no {label} samples")
        centroids.append(features[mask].mean(axis=0))
    return np.vstack(centroids)


def _distances_to_centroids(features: np.ndarray, centroids: np.ndarray) -> np.ndarray:
    squared = np.square(features[:, np.newaxis, :] - centroids[np.newaxis, :, :]).sum(axis=2)
    return np.sqrt(squared)


def _flight_probability(ground_distance: np.ndarray, flight_distance: np.ndarray) -> np.ndarray:
    total = ground_distance + flight_distance
    return np.divide(
        ground_distance,
        total,
        out=np.full_like(total, 0.5, dtype=float),
        where=total > 0.0,
    )


def _format_float(value: float) -> str:
    return f"{value:.10g}"


def predict_fold(
    fold: BulkFoldData,
    *,
    config: NearestCentroidConfig | None = None,
) -> FoldPredictionResult:
    """Fit centroids on one training fold and predict the held-out mission."""

    cfg = config or NearestCentroidConfig()
    train_x = _read_features(fold.paths["train_X.csv"])
    test_x = _align_test_features(
        train_x,
        _read_features(fold.paths["test_X.csv"]),
        fold.paths["test_X.csv"],
    )
    train_labels = _read_labels(fold.paths["train_y.csv"], train_x.index)
    test_labels = _read_labels(fold.paths["test_y.csv"], test_x.index)
    missions = _read_missions(fold.paths["test_meta.csv"], test_x.index, fold.test_mission)

    train_values, test_values = _scale_train_test(
        train_x.to_numpy(dtype=float),
        test_x.to_numpy(dtype=float),
        scaling=cfg.scaling,
    )
    centroids = _class_centroids(train_values, train_labels)
    distances = _distances_to_centroids(test_values, centroids)
    ground_distance = distances[:, 0]
    flight_distance = distances[:, 1]
    probabilities = _flight_probability(ground_distance, flight_distance)
    predicted = np.where(flight_distance <= ground_distance, POSITIVE_LABEL, NEGATIVE_LABEL)

    rows: list[dict[str, str]] = []
    for index, sample_id in enumerate(test_x.index):
        rows.append(
            {
                "task_id": fold.task_id,
                "fold_id": fold.fold_id,
                "sample_id": str(sample_id),
                "mission": missions[index],
                "true_label": test_labels[index],
                "predicted_label": str(predicted[index]),
                "flight_probability": _format_float(float(probabilities[index])),
                "embedding_0": _format_float(float(ground_distance[index])),
                "embedding_1": _format_float(float(flight_distance[index])),
            }
        )

    return FoldPredictionResult(
        task_id=fold.task_id,
        fold_id=fold.fold_id,
        test_mission=fold.test_mission,
        n_train=len(train_x),
        n_test=len(test_x),
        predictions=rows,
    )


def _write_predictions(path: Path, rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PREDICTION_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def run_nearest_centroid_task(
    task: BulkTaskData,
    *,
    output_dir: str | Path,
    config: NearestCentroidConfig | None = None,
    task_manifest_path: str | Path | None = None,
    command: Sequence[str] | None = None,
) -> NearestCentroidTaskResult:
    """Run nearest-centroid predictions and evaluation for one bulk task."""

    cfg = config or NearestCentroidConfig()
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    prediction_rows: list[dict[str, str]] = []
    fold_summaries: list[dict[str, Any]] = []
    for fold in task.folds:
        fold_result = predict_fold(fold, config=cfg)
        prediction_rows.extend(fold_result.predictions)
        fold_summaries.append(
            {
                "fold_id": fold_result.fold_id,
                "test_mission": fold_result.test_mission,
                "n_train": fold_result.n_train,
                "n_test": fold_result.n_test,
            }
        )

    predictions_path = outdir / "predictions.csv"
    _write_predictions(predictions_path, prediction_rows)
    evaluation = evaluate_submission(task.manifest, predictions_path)
    evaluation["baseline"] = {
        "baseline_id": BASELINE_ID,
        "scaling": cfg.scaling,
        "score_column": "flight_probability",
        "embedding_columns": ["embedding_0", "embedding_1"],
        "embedding_interpretation": "Euclidean distances to Ground and Flight centroids.",
    }
    evaluation["folds"] = fold_summaries
    outputs = write_evaluation_report(
        evaluation_result=evaluation,
        task_manifest=task.manifest,
        task_manifest_path=task_manifest_path,
        submission_path=predictions_path,
        output_dir=outdir,
        command=command if command is not None else sys.argv,
    )
    return NearestCentroidTaskResult(
        task_id=task.task_id,
        output_dir=outdir,
        predictions_path=predictions_path,
        metrics_path=outputs["metrics"],
        run_manifest_path=outputs["run_manifest"],
        evaluation=evaluation,
        n_predictions=len(prediction_rows),
    )


def _metric_value(evaluation: Mapping[str, Any], metric_id: str) -> str:
    metric = evaluation.get("metrics", {}).get(metric_id, {})
    if metric.get("status") != "computed":
        return ""
    return _format_float(float(metric["value"]))


def task_result_summary_row(result: NearestCentroidTaskResult) -> dict[str, str]:
    evaluation = result.evaluation
    return {
        "task_id": result.task_id,
        "status": str(evaluation.get("status", "")),
        "validation_ok": str(evaluation.get("validation", {}).get("ok", "")),
        "n_predictions": str(result.n_predictions),
        "macro_f1": _metric_value(evaluation, "macro_f1"),
        "balanced_accuracy": _metric_value(evaluation, "balanced_accuracy"),
        "auroc": _metric_value(evaluation, "auroc"),
        "calibration_error": _metric_value(evaluation, "calibration_error"),
        "mission_discrimination": _metric_value(evaluation, "mission_discrimination"),
        "output_dir": result.output_dir.as_posix(),
        "predictions": result.predictions_path.as_posix(),
        "metrics": result.metrics_path.as_posix(),
        "run_manifest": result.run_manifest_path.as_posix(),
        "error": "",
    }


def failure_summary_row(task_id: str, output_dir: Path, error: Exception) -> dict[str, str]:
    return {
        "task_id": task_id,
        "status": "failed",
        "validation_ok": "",
        "n_predictions": "",
        "macro_f1": "",
        "balanced_accuracy": "",
        "auroc": "",
        "calibration_error": "",
        "mission_discrimination": "",
        "output_dir": output_dir.as_posix(),
        "predictions": "",
        "metrics": "",
        "run_manifest": "",
        "error": f"{type(error).__name__}: {error}",
    }


def write_summary(output_root: str | Path, rows: Sequence[Mapping[str, str]]) -> dict[str, Path]:
    outdir = Path(output_root)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "bulk_lomo_summary.csv"
    json_path = outdir / "bulk_lomo_summary.json"
    fieldnames = [
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
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(list(rows), indent=2, sort_keys=True) + "\n")
    return {"csv": csv_path, "json": json_path}


def run_nearest_centroid_all(
    *,
    task_ids: Sequence[str] | None = None,
    manifest_dir: str | Path = "v9/task_manifests",
    repo_root: str | Path = ".",
    output_root: str | Path = "v9/reports/nearest_centroid",
    config: NearestCentroidConfig | None = None,
    command: Sequence[str] | None = None,
    fail_fast: bool = False,
) -> tuple[list[dict[str, str]], dict[str, Path]]:
    """Run nearest-centroid reports for one or more bulk tasks."""

    cfg = config or NearestCentroidConfig()
    if task_ids:
        tasks = [
            load_bulk_task(task_id, manifest_dir=manifest_dir, repo_root=repo_root)
            for task_id in task_ids
        ]
    else:
        tasks = load_all_bulk_tasks(manifest_dir=manifest_dir, repo_root=repo_root)

    manifest_base = Path(manifest_dir)
    rows: list[dict[str, str]] = []
    for task in tasks:
        task_output = Path(output_root) / task.task_id
        try:
            result = run_nearest_centroid_task(
                task,
                output_dir=task_output,
                config=cfg,
                task_manifest_path=manifest_base / f"{task.task_id}.json",
                command=command,
            )
            rows.append(task_result_summary_row(result))
        except Exception as exc:
            task_output.mkdir(parents=True, exist_ok=True)
            failure_path = task_output / "failure.json"
            failure_path.write_text(
                json.dumps(
                    {
                        "task_id": task.task_id,
                        "baseline_id": BASELINE_ID,
                        "error_type": type(exc).__name__,
                        "error": str(exc),
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n"
            )
            rows.append(failure_summary_row(task.task_id, task_output, exc))
            if fail_fast:
                raise

    return rows, write_summary(output_root, rows)

"""sklearn classifier baselines for v9 bulk leave-one-mission-out tasks."""

from __future__ import annotations

import csv
import json
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from spacebio_bench.data import BulkFoldData, BulkTaskData, load_all_bulk_tasks, load_bulk_task
from spacebio_bench.evaluate import evaluate_submission
from spacebio_bench.reports import write_evaluation_report

from .nearest_centroid import (
    NEGATIVE_LABEL,
    POSITIVE_LABEL,
    _align_test_features,
    _format_float,
    _read_features,
    _read_labels,
    _read_missions,
)


BASELINE_ALIASES = {
    "logistic": "logistic_regression_l2",
    "logistic_l2": "logistic_regression_l2",
    "logistic_regression_l2": "logistic_regression_l2",
    "pca_lr": "pca_logistic_regression",
    "pca_logistic": "pca_logistic_regression",
    "pca_logistic_regression": "pca_logistic_regression",
}
DEFAULT_BASELINES = ("pca_logistic_regression", "logistic_regression_l2")
BASE_PREDICTION_COLUMNS = [
    "task_id",
    "fold_id",
    "sample_id",
    "mission",
    "true_label",
    "predicted_label",
    "flight_probability",
]


@dataclass(frozen=True)
class SklearnBaselineConfig:
    """Configuration for sklearn bulk classifier baselines."""

    baseline_id: str = "pca_logistic_regression"
    pca_components: int = 50
    max_iter: int = 5000
    random_state: int = 42
    c: float = 1.0

    def __post_init__(self) -> None:
        normalized = normalize_baseline_id(self.baseline_id)
        object.__setattr__(self, "baseline_id", normalized)
        if self.pca_components < 1:
            raise ValueError("pca_components must be at least 1")
        if self.max_iter < 1:
            raise ValueError("max_iter must be at least 1")
        if self.c <= 0:
            raise ValueError("c must be positive")


@dataclass(frozen=True)
class SklearnFoldPredictionResult:
    """Predictions and fold-level metadata for one sklearn baseline fold."""

    task_id: str
    fold_id: str
    test_mission: str
    n_train: int
    n_test: int
    n_features: int
    train_time_s: float
    fit_details: Mapping[str, Any]
    predictions: list[dict[str, str]]


@dataclass(frozen=True)
class SklearnBaselineTaskResult:
    """Filesystem outputs and metrics for one sklearn baseline task run."""

    baseline_id: str
    task_id: str
    output_dir: Path
    predictions_path: Path
    metrics_path: Path
    run_manifest_path: Path
    evaluation: Mapping[str, Any]
    n_predictions: int


def normalize_baseline_id(value: str) -> str:
    try:
        return BASELINE_ALIASES[value]
    except KeyError as exc:
        valid = ", ".join(sorted(BASELINE_ALIASES))
        raise ValueError(f"unknown sklearn baseline {value!r}; valid values: {valid}") from exc


def _numeric_labels(labels: Sequence[str]) -> np.ndarray:
    mapping = {NEGATIVE_LABEL: 0, POSITIVE_LABEL: 1}
    try:
        return np.asarray([mapping[label] for label in labels], dtype=int)
    except KeyError as exc:
        raise ValueError(f"unexpected label for binary classifier: {exc.args[0]!r}") from exc


def _prediction_label(value: int) -> str:
    return POSITIVE_LABEL if int(value) == 1 else NEGATIVE_LABEL


def _adaptive_pca_components(config: SklearnBaselineConfig, n_train: int, n_features: int) -> int:
    return max(1, min(config.pca_components, n_train - 1, n_features))


def _build_model(config: SklearnBaselineConfig, n_train: int, n_features: int):
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    if config.baseline_id == "pca_logistic_regression":
        n_components = _adaptive_pca_components(config, n_train, n_features)
        return Pipeline(
            [
                ("scaler", StandardScaler()),
                ("pca", PCA(n_components=n_components, random_state=config.random_state)),
                (
                    "classifier",
                    LogisticRegression(
                        C=config.c,
                        class_weight="balanced",
                        max_iter=config.max_iter,
                        random_state=config.random_state,
                    ),
                ),
            ]
        )

    if config.baseline_id == "logistic_regression_l2":
        return Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "classifier",
                    LogisticRegression(
                        C=config.c,
                        class_weight="balanced",
                        max_iter=config.max_iter,
                        random_state=config.random_state,
                        solver="liblinear",
                    ),
                ),
            ]
        )

    raise ValueError(f"unsupported sklearn baseline: {config.baseline_id}")


def _positive_class_index(model: Any) -> int:
    classes = list(model.named_steps["classifier"].classes_)
    try:
        return classes.index(1)
    except ValueError as exc:
        raise ValueError(f"classifier did not learn positive class 1: {classes}") from exc


def _pca_embeddings(model: Any, test_values: np.ndarray) -> list[tuple[float, float]]:
    if "pca" not in model.named_steps:
        return []
    scaled = model.named_steps["scaler"].transform(test_values)
    transformed = model.named_steps["pca"].transform(scaled)
    if transformed.shape[1] == 1:
        return [(float(value), 0.0) for value in transformed[:, 0]]
    return [(float(row[0]), float(row[1])) for row in transformed]


def _prediction_fieldnames(rows: Sequence[Mapping[str, str]]) -> list[str]:
    fields = list(BASE_PREDICTION_COLUMNS)
    if rows and "embedding_0" in rows[0]:
        fields.extend(["embedding_0", "embedding_1"])
    return fields


def _write_predictions(path: Path, rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_prediction_fieldnames(rows))
        writer.writeheader()
        writer.writerows(rows)


def predict_sklearn_fold(
    fold: BulkFoldData,
    *,
    config: SklearnBaselineConfig | None = None,
) -> SklearnFoldPredictionResult:
    """Fit one sklearn classifier on a training fold and predict the held-out mission."""

    cfg = config or SklearnBaselineConfig()
    train_x = _read_features(fold.paths["train_X.csv"])
    test_x = _align_test_features(
        train_x,
        _read_features(fold.paths["test_X.csv"]),
        fold.paths["test_X.csv"],
    )
    train_labels = _read_labels(fold.paths["train_y.csv"], train_x.index)
    test_labels = _read_labels(fold.paths["test_y.csv"], test_x.index)
    missions = _read_missions(fold.paths["test_meta.csv"], test_x.index, fold.test_mission)
    y_train = _numeric_labels(train_labels)

    train_values = train_x.to_numpy(dtype=np.float32)
    test_values = test_x.to_numpy(dtype=np.float32)
    model = _build_model(cfg, len(train_x), train_x.shape[1])

    start = time.perf_counter()
    model.fit(train_values, y_train)
    train_time_s = time.perf_counter() - start

    positive_index = _positive_class_index(model)
    probabilities = model.predict_proba(test_values)[:, positive_index]
    predicted_numeric = model.predict(test_values)
    embeddings = _pca_embeddings(model, test_values)

    rows: list[dict[str, str]] = []
    for index, sample_id in enumerate(test_x.index):
        row = {
            "task_id": fold.task_id,
            "fold_id": fold.fold_id,
            "sample_id": str(sample_id),
            "mission": missions[index],
            "true_label": test_labels[index],
            "predicted_label": _prediction_label(int(predicted_numeric[index])),
            "flight_probability": _format_float(float(probabilities[index])),
        }
        if embeddings:
            row["embedding_0"] = _format_float(embeddings[index][0])
            row["embedding_1"] = _format_float(embeddings[index][1])
        rows.append(row)

    fit_details: dict[str, Any] = {
        "model_type": cfg.baseline_id,
        "n_features": int(train_x.shape[1]),
        "train_time_s": round(train_time_s, 6),
    }
    if cfg.baseline_id == "pca_logistic_regression":
        fit_details["pca_components"] = int(model.named_steps["pca"].n_components)

    return SklearnFoldPredictionResult(
        task_id=fold.task_id,
        fold_id=fold.fold_id,
        test_mission=fold.test_mission,
        n_train=len(train_x),
        n_test=len(test_x),
        n_features=train_x.shape[1],
        train_time_s=train_time_s,
        fit_details=fit_details,
        predictions=rows,
    )


def _embedding_metadata(config: SklearnBaselineConfig) -> dict[str, Any]:
    if config.baseline_id == "pca_logistic_regression":
        return {
            "embedding_columns": ["embedding_0", "embedding_1"],
            "embedding_interpretation": "Train-fitted PCA coordinates before logistic regression.",
        }
    return {
        "embedding_columns": [],
        "embedding_interpretation": "No embedding columns are emitted for plain L2 logistic regression.",
    }


def run_sklearn_baseline_task(
    task: BulkTaskData,
    *,
    output_dir: str | Path,
    config: SklearnBaselineConfig | None = None,
    task_manifest_path: str | Path | None = None,
    command: Sequence[str] | None = None,
) -> SklearnBaselineTaskResult:
    """Run one sklearn baseline and evaluation report for one bulk task."""

    cfg = config or SklearnBaselineConfig()
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    prediction_rows: list[dict[str, str]] = []
    fold_summaries: list[dict[str, Any]] = []
    for fold in task.folds:
        fold_result = predict_sklearn_fold(fold, config=cfg)
        prediction_rows.extend(fold_result.predictions)
        fold_summaries.append(
            {
                "fold_id": fold_result.fold_id,
                "test_mission": fold_result.test_mission,
                "n_train": fold_result.n_train,
                "n_test": fold_result.n_test,
                "n_features": fold_result.n_features,
                "train_time_s": round(fold_result.train_time_s, 6),
                "fit_details": dict(fold_result.fit_details),
            }
        )

    predictions_path = outdir / "predictions.csv"
    _write_predictions(predictions_path, prediction_rows)
    evaluation = evaluate_submission(task.manifest, predictions_path)
    evaluation["baseline"] = {
        "baseline_id": cfg.baseline_id,
        "framework": "sklearn",
        "score_column": "flight_probability",
        "random_state": cfg.random_state,
        "max_iter": cfg.max_iter,
        "c": cfg.c,
        **_embedding_metadata(cfg),
    }
    if cfg.baseline_id == "pca_logistic_regression":
        evaluation["baseline"]["requested_pca_components"] = cfg.pca_components
    evaluation["folds"] = fold_summaries
    outputs = write_evaluation_report(
        evaluation_result=evaluation,
        task_manifest=task.manifest,
        task_manifest_path=task_manifest_path,
        submission_path=predictions_path,
        output_dir=outdir,
        command=command if command is not None else sys.argv,
    )
    return SklearnBaselineTaskResult(
        baseline_id=cfg.baseline_id,
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


def task_result_summary_row(result: SklearnBaselineTaskResult) -> dict[str, str]:
    evaluation = result.evaluation
    return {
        "baseline_id": result.baseline_id,
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


def failure_summary_row(
    *,
    baseline_id: str,
    task_id: str,
    output_dir: Path,
    error: Exception,
) -> dict[str, str]:
    return {
        "baseline_id": baseline_id,
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
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(list(rows), indent=2, sort_keys=True) + "\n")
    return {"csv": csv_path, "json": json_path}


def run_sklearn_baseline_all(
    *,
    baseline_ids: Sequence[str] | None = None,
    task_ids: Sequence[str] | None = None,
    manifest_dir: str | Path = "v9/task_manifests",
    repo_root: str | Path = ".",
    output_root: str | Path = "v9/reports/sklearn_baselines",
    pca_components: int = 50,
    max_iter: int = 5000,
    random_state: int = 42,
    c: float = 1.0,
    command: Sequence[str] | None = None,
    fail_fast: bool = False,
) -> tuple[list[dict[str, str]], dict[str, Path]]:
    """Run sklearn baseline reports for one or more bulk tasks."""

    normalized_baselines = [
        normalize_baseline_id(baseline)
        for baseline in (baseline_ids if baseline_ids else DEFAULT_BASELINES)
    ]
    if task_ids:
        tasks = [
            load_bulk_task(task_id, manifest_dir=manifest_dir, repo_root=repo_root)
            for task_id in task_ids
        ]
    else:
        tasks = load_all_bulk_tasks(manifest_dir=manifest_dir, repo_root=repo_root)

    manifest_base = Path(manifest_dir)
    rows: list[dict[str, str]] = []
    for baseline_id in normalized_baselines:
        config = SklearnBaselineConfig(
            baseline_id=baseline_id,
            pca_components=pca_components,
            max_iter=max_iter,
            random_state=random_state,
            c=c,
        )
        for task in tasks:
            task_output = Path(output_root) / config.baseline_id / task.task_id
            try:
                result = run_sklearn_baseline_task(
                    task,
                    output_dir=task_output,
                    config=config,
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
                            "baseline_id": config.baseline_id,
                            "task_id": task.task_id,
                            "error_type": type(exc).__name__,
                            "error": str(exc),
                        },
                        indent=2,
                        sort_keys=True,
                    )
                    + "\n"
                )
                rows.append(
                    failure_summary_row(
                        baseline_id=config.baseline_id,
                        task_id=task.task_id,
                        output_dir=task_output,
                        error=exc,
                    )
                )
                if fail_fast:
                    raise

    return rows, write_summary(output_root, rows)

"""Draft human organoid baseline runners."""

from __future__ import annotations

import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from spacebio_bench.data import HumanOrganoidTaskData, OrganoidFoldData
from spacebio_bench.data import load_human_organoid_task
from spacebio_bench.evaluate import evaluate_submission
from spacebio_bench.reports import write_evaluation_report


BASELINE_ID = "organoid_nearest_centroid"
NEGATIVE_LABEL = "Ground"
POSITIVE_LABEL = "LEO_or_ISS"
PREDICTION_COLUMNS = [
    "task_id",
    "fold_id",
    "sample_id",
    "mission",
    "true_label",
    "predicted_label",
    "flight_probability",
]


@dataclass(frozen=True)
class OrganoidBaselineConfig:
    """Configuration for the draft human organoid centroid baseline."""

    transform: str = "log1p"
    scaling: str = "zscore"
    top_variable_genes: int = 2000

    def __post_init__(self) -> None:
        if self.transform not in {"none", "log1p"}:
            raise ValueError("transform must be 'none' or 'log1p'")
        if self.scaling not in {"none", "zscore"}:
            raise ValueError("scaling must be 'none' or 'zscore'")
        if self.top_variable_genes < 1:
            raise ValueError("top_variable_genes must be at least 1")


@dataclass(frozen=True)
class OrganoidFoldPredictionResult:
    """Predictions and fold metadata for one organoid draft fold."""

    task_id: str
    fold_id: str
    heldout_factor: str
    heldout_value: str
    n_train: int
    n_test: int
    n_features: int
    predictions: list[dict[str, str]]


@dataclass(frozen=True)
class OrganoidBaselineTaskResult:
    """Filesystem outputs and metrics for one organoid baseline run."""

    baseline_id: str
    task_id: str
    output_dir: Path
    predictions_path: Path
    metrics_path: Path
    run_manifest_path: Path
    evaluation: Mapping[str, Any]
    n_predictions: int


def _format_float(value: float) -> str:
    return f"{value:.10g}"


def _sample_factor_by_id(task: HumanOrganoidTaskData) -> dict[str, dict[str, str]]:
    return {
        str(row["sample_id"]): row
        for row in task.sample_factors
        if row.get("parse_status") == "parsed" and row.get("sample_id")
    }


def _transform(values: np.ndarray, mode: str) -> np.ndarray:
    if mode == "none":
        return values.astype(np.float64, copy=True)
    if mode == "log1p":
        clipped = np.clip(values.astype(np.float64, copy=True), a_min=0.0, a_max=None)
        return np.log1p(clipped)
    raise ValueError(f"unsupported transform: {mode}")


def _select_train_variable_features(
    train: np.ndarray,
    test: np.ndarray,
    *,
    top_variable_genes: int,
) -> tuple[np.ndarray, np.ndarray]:
    n_features = train.shape[1]
    if top_variable_genes >= n_features:
        return train, test
    variances = np.var(train, axis=0)
    selected = np.argsort(-variances, kind="mergesort")[:top_variable_genes]
    selected.sort()
    return train[:, selected], test[:, selected]


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


def _positive_probability(ground_distance: np.ndarray, leo_distance: np.ndarray) -> np.ndarray:
    total = ground_distance + leo_distance
    return np.divide(
        ground_distance,
        total,
        out=np.full_like(total, 0.5, dtype=float),
        where=total > 0.0,
    )


def predict_organoid_fold(
    task: HumanOrganoidTaskData,
    fold: OrganoidFoldData,
    *,
    config: OrganoidBaselineConfig | None = None,
) -> OrganoidFoldPredictionResult:
    """Fit a transparent centroid classifier on one draft organoid fold."""

    cfg = config or OrganoidBaselineConfig()
    sample_factors = _sample_factor_by_id(task)
    train_x = task.features.loc[fold.train_sample_ids]
    test_x = task.features.loc[fold.test_sample_ids]
    train_labels = [sample_factors[sample_id]["true_label"] for sample_id in train_x.index]
    test_labels = [sample_factors[sample_id]["true_label"] for sample_id in test_x.index]

    train_values = _transform(train_x.to_numpy(dtype=np.float64), cfg.transform)
    test_values = _transform(test_x.to_numpy(dtype=np.float64), cfg.transform)
    train_values, test_values = _select_train_variable_features(
        train_values,
        test_values,
        top_variable_genes=cfg.top_variable_genes,
    )
    train_values, test_values = _scale_train_test(
        train_values,
        test_values,
        scaling=cfg.scaling,
    )
    centroids = _class_centroids(train_values, train_labels)
    distances = _distances_to_centroids(test_values, centroids)
    ground_distance = distances[:, 0]
    leo_distance = distances[:, 1]
    probabilities = _positive_probability(ground_distance, leo_distance)
    predicted = np.where(leo_distance <= ground_distance, POSITIVE_LABEL, NEGATIVE_LABEL)

    rows: list[dict[str, str]] = []
    for index, sample_id in enumerate(test_x.index):
        row = sample_factors[str(sample_id)]
        rows.append(
            {
                "task_id": task.task_id,
                "fold_id": fold.fold_id,
                "sample_id": str(sample_id),
                "mission": str(row.get("mission", "")),
                "true_label": test_labels[index],
                "predicted_label": str(predicted[index]),
                "flight_probability": _format_float(float(probabilities[index])),
            }
        )

    return OrganoidFoldPredictionResult(
        task_id=task.task_id,
        fold_id=fold.fold_id,
        heldout_factor=fold.heldout_factor,
        heldout_value=fold.heldout_value,
        n_train=len(train_x),
        n_test=len(test_x),
        n_features=int(train_values.shape[1]),
        predictions=rows,
    )


def _write_predictions(path: Path, rows: Sequence[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PREDICTION_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def run_organoid_nearest_centroid(
    task: HumanOrganoidTaskData,
    *,
    output_dir: str | Path,
    config: OrganoidBaselineConfig | None = None,
    folds: Sequence[OrganoidFoldData] | None = None,
    fold_family: str = "default_sample_count_backed_folds",
    claim_boundary: str = "pilot_baseline_only_not_leaderboard",
    task_manifest_path: str | Path | None = None,
    command: Sequence[str] | None = None,
) -> OrganoidBaselineTaskResult:
    """Run the draft organoid nearest-centroid baseline and evaluation report."""

    cfg = config or OrganoidBaselineConfig()
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    prediction_rows: list[dict[str, str]] = []
    fold_summaries: list[dict[str, Any]] = []
    selected_folds = list(folds or task.folds)
    for fold in selected_folds:
        fold_result = predict_organoid_fold(task, fold, config=cfg)
        prediction_rows.extend(fold_result.predictions)
        fold_summaries.append(
            {
                "fold_id": fold_result.fold_id,
                "heldout_factor": fold_result.heldout_factor,
                "heldout_value": fold_result.heldout_value,
                "n_train": fold_result.n_train,
                "n_test": fold_result.n_test,
                "n_features": fold_result.n_features,
            }
        )

    predictions_path = outdir / "predictions.csv"
    _write_predictions(predictions_path, prediction_rows)
    evaluation = evaluate_submission(task.manifest, predictions_path)
    evaluation["baseline"] = {
        "baseline_id": BASELINE_ID,
        "release_status": "draft_not_frozen",
        "transform": cfg.transform,
        "scaling": cfg.scaling,
        "top_variable_genes": cfg.top_variable_genes,
        "score_column": "flight_probability",
        "positive_label": POSITIVE_LABEL,
        "fold_family": fold_family,
        "claim_boundary": claim_boundary,
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
    return OrganoidBaselineTaskResult(
        baseline_id=BASELINE_ID,
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


def organoid_result_summary_row(result: OrganoidBaselineTaskResult) -> dict[str, str]:
    evaluation = result.evaluation
    baseline = evaluation.get("baseline", {})
    claim_boundary = str(
        baseline.get("claim_boundary", "pilot_baseline_only_not_leaderboard")
    )
    return {
        "baseline_id": result.baseline_id,
        "variant_id": organoid_config_id(
            OrganoidBaselineConfig(
                transform=str(baseline.get("transform", "log1p")),
                scaling=str(baseline.get("scaling", "zscore")),
                top_variable_genes=int(baseline.get("top_variable_genes", 2000)),
            )
        ),
        "transform": str(baseline.get("transform", "")),
        "scaling": str(baseline.get("scaling", "")),
        "top_variable_genes": str(baseline.get("top_variable_genes", "")),
        "task_id": result.task_id,
        "status": str(evaluation.get("status", "")),
        "validation_ok": str(evaluation.get("validation", {}).get("ok", "")),
        "release_status": "draft_not_frozen",
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
        "fold_family": str(baseline.get("fold_family", "")),
        "claim_boundary": claim_boundary,
    }


def write_organoid_baseline_summary(
    output_root: str | Path,
    rows: Sequence[Mapping[str, str]],
) -> dict[str, Path]:
    outdir = Path(output_root)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "human_organoid_baseline_summary.csv"
    json_path = outdir / "human_organoid_baseline_summary.json"
    fieldnames = [
        "baseline_id",
        "variant_id",
        "transform",
        "scaling",
        "top_variable_genes",
        "task_id",
        "status",
        "validation_ok",
        "release_status",
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
        "fold_family",
        "claim_boundary",
    ]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(list(rows), indent=2, sort_keys=True) + "\n")
    return {"csv": csv_path, "json": json_path}


def organoid_config_id(config: OrganoidBaselineConfig) -> str:
    return f"tvg{config.top_variable_genes}_{config.transform}_{config.scaling}"


def run_human_organoid_baseline(
    *,
    manifest_path: str | Path = "v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json",
    repo_root: str | Path = ".",
    output_root: str | Path = "v9/human_organoid/reports/nearest_centroid",
    config: OrganoidBaselineConfig | None = None,
    command: Sequence[str] | None = None,
) -> tuple[OrganoidBaselineTaskResult, dict[str, Path]]:
    """Load the draft organoid task, run baseline, and write a summary."""

    task = load_human_organoid_task(manifest_path=manifest_path, repo_root=repo_root)
    result = run_organoid_nearest_centroid(
        task,
        output_dir=Path(output_root) / task.task_id,
        config=config,
        task_manifest_path=Path(repo_root) / manifest_path,
        command=command,
    )
    summary = write_organoid_baseline_summary(
        output_root,
        [organoid_result_summary_row(result)],
    )
    return result, summary


def run_human_organoid_donor_diagnostics(
    *,
    manifest_path: str | Path = "v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json",
    repo_root: str | Path = ".",
    output_root: str | Path = "v9/human_organoid/reports/donor_diagnostics",
    config: OrganoidBaselineConfig | None = None,
    command: Sequence[str] | None = None,
) -> tuple[OrganoidBaselineTaskResult, dict[str, Path]]:
    """Evaluate GEO-derived donor holdouts as a diagnostic-only fold family."""

    task = load_human_organoid_task(manifest_path=manifest_path, repo_root=repo_root)
    result = run_organoid_nearest_centroid(
        task,
        output_dir=Path(output_root) / task.task_id,
        config=config,
        folds=task.diagnostic_folds,
        fold_family="geo_donor_or_line_diagnostic_folds",
        claim_boundary="donor_diagnostic_only_not_leaderboard",
        task_manifest_path=Path(repo_root) / manifest_path,
        command=command,
    )
    summary = write_organoid_baseline_summary(
        output_root,
        [organoid_result_summary_row(result)],
    )
    return result, summary


def run_human_organoid_sensitivity_grid(
    *,
    manifest_path: str | Path = "v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json",
    repo_root: str | Path = ".",
    output_root: str | Path = "v9/human_organoid/reports/sensitivity",
    configs: Sequence[OrganoidBaselineConfig] | None = None,
    command: Sequence[str] | None = None,
) -> tuple[list[OrganoidBaselineTaskResult], dict[str, Path]]:
    """Run a small sensitivity grid for the draft human organoid baseline."""

    task = load_human_organoid_task(manifest_path=manifest_path, repo_root=repo_root)
    default_configs = [
        OrganoidBaselineConfig(transform=transform, scaling=scaling, top_variable_genes=top)
        for transform in ("log1p", "none")
        for scaling in ("zscore", "none")
        for top in (100, 500, 2000, 5000, task.n_features)
    ]
    selected_configs = list(configs or default_configs)
    results: list[OrganoidBaselineTaskResult] = []
    rows: list[dict[str, str]] = []
    root = Path(output_root)
    for config in selected_configs:
        variant_id = organoid_config_id(config)
        result = run_organoid_nearest_centroid(
            task,
            output_dir=root / variant_id / task.task_id,
            config=config,
            task_manifest_path=Path(repo_root) / manifest_path,
            command=command,
        )
        results.append(result)
        rows.append(organoid_result_summary_row(result))
    summary = write_organoid_baseline_summary(output_root, rows)
    return results, summary

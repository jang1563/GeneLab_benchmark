"""Evaluation helpers for SpaceBio-Bench prediction submissions."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Mapping, Sequence

from .metrics import mission_discrimination_score
from .signature_metrics import (
    RESPONSE_SIGNATURE_REQUIRED_COLUMNS,
    SIGNATURE_METRIC_IDS,
    compute_response_signature_metrics,
)
from .submissions import validate_submission


def _read_rows(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _label_domain(manifest: Mapping[str, Any], rows: Sequence[Mapping[str, str]]) -> list[str]:
    output = manifest.get("output", {})
    if isinstance(output, Mapping) and output.get("label_domain"):
        return [str(label) for label in output["label_domain"]]
    labels = sorted(
        {
            str(row[label_column])
            for row in rows
            for label_column in ("true_label", "predicted_label")
            if row.get(label_column)
        }
    )
    preferred = ["Flight", "Ground"]
    if set(labels) == set(preferred):
        return preferred
    return labels


def _positive_label(manifest: Mapping[str, Any], labels: Sequence[str]) -> str:
    output = manifest.get("output", {})
    if isinstance(output, Mapping) and output.get("positive_label"):
        return str(output["positive_label"])
    if "Flight" in labels:
        return "Flight"
    if not labels:
        raise ValueError("cannot infer positive label without labels")
    return labels[0]


def _confusion_counts(rows: Sequence[Mapping[str, str]], labels: Sequence[str]) -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = {
        label: {"tp": 0, "fp": 0, "fn": 0, "tn": 0} for label in labels
    }
    for label in labels:
        for row in rows:
            truth = row["true_label"]
            pred = row["predicted_label"]
            if truth == label and pred == label:
                counts[label]["tp"] += 1
            elif truth != label and pred == label:
                counts[label]["fp"] += 1
            elif truth == label and pred != label:
                counts[label]["fn"] += 1
            else:
                counts[label]["tn"] += 1
    return counts


def _macro_f1(rows: Sequence[Mapping[str, str]], labels: Sequence[str]) -> float:
    counts = _confusion_counts(rows, labels)
    scores: list[float] = []
    for label in labels:
        tp = counts[label]["tp"]
        fp = counts[label]["fp"]
        fn = counts[label]["fn"]
        denom = 2 * tp + fp + fn
        scores.append(0.0 if denom == 0 else (2 * tp) / denom)
    return sum(scores) / len(scores)


def _balanced_accuracy(rows: Sequence[Mapping[str, str]], labels: Sequence[str]) -> float:
    counts = _confusion_counts(rows, labels)
    recalls: list[float] = []
    for label in labels:
        tp = counts[label]["tp"]
        fn = counts[label]["fn"]
        denom = tp + fn
        recalls.append(0.0 if denom == 0 else tp / denom)
    return sum(recalls) / len(recalls)


def _average_ranks(values: Sequence[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    cursor = 0
    while cursor < len(indexed):
        end = cursor + 1
        while end < len(indexed) and indexed[end][1] == indexed[cursor][1]:
            end += 1
        average_rank = (cursor + 1 + end) / 2.0
        for index, _ in indexed[cursor:end]:
            ranks[index] = average_rank
        cursor = end
    return ranks


def _binary_auroc(
    rows: Sequence[Mapping[str, str]],
    *,
    positive_label: str,
    probability_column: str = "flight_probability",
) -> float:
    labels = [1 if row["true_label"] == positive_label else 0 for row in rows]
    scores = [float(row[probability_column]) for row in rows]
    n_pos = sum(labels)
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        raise ValueError("AUROC requires at least one positive and one negative sample")
    ranks = _average_ranks(scores)
    pos_rank_sum = sum(rank for rank, label in zip(ranks, labels) if label == 1)
    return (pos_rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def _calibration_error(
    rows: Sequence[Mapping[str, str]],
    *,
    positive_label: str,
    probability_column: str = "flight_probability",
    n_bins: int = 10,
) -> float:
    total = len(rows)
    if total == 0:
        raise ValueError("calibration_error requires at least one row")
    error = 0.0
    for bin_index in range(n_bins):
        low = bin_index / n_bins
        high = (bin_index + 1) / n_bins
        in_bin = []
        for row in rows:
            prob = float(row[probability_column])
            if (bin_index == n_bins - 1 and low <= prob <= high) or (low <= prob < high):
                in_bin.append(row)
        if not in_bin:
            continue
        mean_prob = sum(float(row[probability_column]) for row in in_bin) / len(in_bin)
        observed = sum(row["true_label"] == positive_label for row in in_bin) / len(in_bin)
        error += (len(in_bin) / total) * abs(mean_prob - observed)
    return error


def _embedding_columns(rows: Sequence[Mapping[str, str]]) -> list[str]:
    if not rows:
        return []
    return sorted(column for column in rows[0] if column.startswith("embedding_"))


def _metric_value(value: float) -> dict[str, Any]:
    return {"status": "computed", "value": value}


def _metric_skip(reason: str) -> dict[str, Any]:
    return {"status": "skipped", "reason": reason}


def _declared_metric_ids(manifest: Mapping[str, Any]) -> list[str]:
    metrics = manifest.get("metrics", [])
    if not isinstance(metrics, Sequence) or isinstance(metrics, (str, bytes)):
        return []
    metric_ids: list[str] = []
    for metric in metrics:
        if isinstance(metric, Mapping) and metric.get("metric_id"):
            metric_ids.append(str(metric["metric_id"]))
    return metric_ids


def _signature_metric_skip_reason(response_signature_path: str | Path | None) -> str:
    columns = ", ".join(RESPONSE_SIGNATURE_REQUIRED_COLUMNS)
    if response_signature_path is None:
        return (
            "response_signature.csv artifact missing; DE/signature metrics require "
            f"a gene-level artifact with columns: {columns}"
        )
    path = Path(response_signature_path)
    if not path.exists():
        return f"response_signature.csv artifact not found: {path}"
    return (
        "response_signature.csv artifact supplied, but DE/signature scoring is "
        "pending the frozen scorer implementation after V9-ORG-014"
    )


def _add_declared_signature_metric_statuses(
    result: dict[str, Any],
    manifest: Mapping[str, Any],
    *,
    response_signature_path: str | Path | None,
    reference_signature_path: str | Path | None,
) -> None:
    declared = set(_declared_metric_ids(manifest))
    declared_signature_metrics = SIGNATURE_METRIC_IDS & declared
    if not declared_signature_metrics:
        return
    if response_signature_path is not None and Path(response_signature_path).exists():
        signature_result = compute_response_signature_metrics(
            manifest=manifest,
            response_signature_path=response_signature_path,
            reference_signature_path=reference_signature_path,
        )
        result["response_signature_validation"] = signature_result["validation"]
        for metric_id in sorted(declared_signature_metrics):
            result["metrics"][metric_id] = signature_result["metrics"][metric_id]
        return
    for metric_id in sorted(declared_signature_metrics):
        if metric_id not in result["metrics"]:
            result["metrics"][metric_id] = _metric_skip(
                _signature_metric_skip_reason(response_signature_path)
            )


def evaluate_submission(
    manifest: Mapping[str, Any],
    submission_path: str | Path,
    *,
    response_signature_path: str | Path | None = None,
    reference_signature_path: str | Path | None = None,
) -> dict[str, Any]:
    """Validate and evaluate a prediction CSV for a task manifest."""

    validation = validate_submission(manifest, submission_path)
    result: dict[str, Any] = {
        "task_id": manifest["task_id"],
        "submission_path": Path(submission_path).as_posix(),
        "validation": validation.to_dict(),
        "metrics": {},
    }
    if response_signature_path is not None:
        result["response_signature_path"] = Path(response_signature_path).as_posix()
    if reference_signature_path is not None:
        result["reference_signature_path"] = Path(reference_signature_path).as_posix()
    if not validation.ok:
        result["status"] = "invalid"
        return result

    rows = _read_rows(submission_path)
    labels = _label_domain(manifest, rows)
    positive_label = _positive_label(manifest, labels)

    result["status"] = "evaluated"
    result["label_domain"] = labels
    result["positive_label"] = positive_label
    result["metrics"]["macro_f1"] = _metric_value(_macro_f1(rows, labels))
    result["metrics"]["balanced_accuracy"] = _metric_value(
        _balanced_accuracy(rows, labels)
    )

    if "flight_probability" in validation.columns:
        try:
            result["metrics"]["auroc"] = _metric_value(
                _binary_auroc(rows, positive_label=positive_label)
            )
        except ValueError as exc:
            result["metrics"]["auroc"] = _metric_skip(str(exc))
        result["metrics"]["calibration_error"] = _metric_value(
            _calibration_error(rows, positive_label=positive_label)
        )
    else:
        result["metrics"]["auroc"] = _metric_skip("flight_probability column missing")
        result["metrics"]["calibration_error"] = _metric_skip(
            "flight_probability column missing"
        )

    embedding_columns = _embedding_columns(rows)
    if embedding_columns and "mission" in validation.columns:
        embeddings = [
            [float(row[column]) for column in embedding_columns]
            for row in rows
        ]
        try:
            mission_result = mission_discrimination_score(
                embeddings,
                [row["mission"] for row in rows],
                sample_ids=[
                    row.get("sample_id", str(index)) for index, row in enumerate(rows)
                ],
            )
        except ValueError as exc:
            result["metrics"]["mission_discrimination"] = _metric_skip(str(exc))
        else:
            result["metrics"]["mission_discrimination"] = {
                "status": "computed",
                "value": mission_result["score"],
                "details": mission_result,
            }
    elif embedding_columns:
        result["metrics"]["mission_discrimination"] = _metric_skip(
            "mission column missing"
        )
    else:
        result["metrics"]["mission_discrimination"] = _metric_skip(
            "embedding_* columns missing"
        )

    _add_declared_signature_metric_statuses(
        result,
        manifest,
        response_signature_path=response_signature_path,
        reference_signature_path=reference_signature_path,
    )
    return result

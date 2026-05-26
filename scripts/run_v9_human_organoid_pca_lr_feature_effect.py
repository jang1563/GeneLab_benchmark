#!/usr/bin/env python3
"""Run the v9 human organoid PCA-LR reconstructed feature-effect pilot."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    build_pca_lr_reconstructed_feature_effect,
    compute_feature_effect_metrics,
    evaluate_submission,
    load_human_organoid_task,
    load_task_manifest,
    write_feature_effect_artifact,
)
from spacebio_bench.reports import write_evaluation_report  # noqa: E402


PREDICTION_ROWS = [
    {
        "sample_id": "pca_lr_feature_effect_smoke_ground",
        "true_label": "Ground",
        "predicted_label": "Ground",
        "flight_probability": "0.05",
    },
    {
        "sample_id": "pca_lr_feature_effect_smoke_leo_or_iss",
        "true_label": "LEO_or_ISS",
        "predicted_label": "LEO_or_ISS",
        "flight_probability": "0.95",
    },
]


def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _metric_value(evaluation: dict[str, object], metric_id: str) -> float | None:
    metrics = evaluation.get("metrics")
    metric = metrics.get(metric_id) if isinstance(metrics, dict) else None
    if not isinstance(metric, dict) or metric.get("status") != "computed":
        return None
    value = metric.get("value")
    return float(value) if isinstance(value, (int, float)) else None


def _format_metric(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.6f}"


def _top_k_summary(evaluation: dict[str, object]) -> list[tuple[str, str, str, str, str, str]]:
    metric = evaluation.get("metrics", {}).get("feature_effect_top_k_de_overlap")
    if not isinstance(metric, dict) or metric.get("status") != "computed":
        return []
    value = metric.get("value")
    if not isinstance(value, dict):
        return []
    rows = []
    for key in sorted(value, key=lambda item: int(item)):
        row = value[key]
        if not isinstance(row, dict):
            continue
        rows.append(
            (
                key,
                str(row.get("n_overlap", "")),
                _format_metric(
                    row.get("overlap_fraction")
                    if isinstance(row.get("overlap_fraction"), (int, float))
                    else None
                ),
                _format_metric(
                    row.get("expected_overlap")
                    if isinstance(row.get("expected_overlap"), (int, float))
                    else None
                ),
                _format_metric(
                    row.get("enrichment")
                    if isinstance(row.get("enrichment"), (int, float))
                    else None
                ),
                _format_metric(
                    row.get("hypergeometric_p_value_greater_equal")
                    if isinstance(row.get("hypergeometric_p_value_greater_equal"), (int, float))
                    else None
                ),
            )
        )
    return rows


def _load_comparison_metrics(path: str | Path | None) -> dict[str, object] | None:
    if not path:
        return None
    parsed = Path(path)
    if not parsed.exists():
        return None
    return json.loads(parsed.read_text())


def _comparison_value(metrics: dict[str, object] | None, metric_id: str) -> float | None:
    if not metrics:
        return None
    metric = metrics.get("metrics", {}).get(metric_id)
    if not isinstance(metric, dict) or metric.get("status") != "computed":
        return None
    value = metric.get("value")
    return float(value) if isinstance(value, (int, float)) else None


def _write_note(
    *,
    output_dir: Path,
    metadata: dict[str, object],
    evaluation: dict[str, object],
    comparison_metrics: dict[str, object] | None,
) -> Path:
    note_path = output_dir / "README.md"
    lines = [
        "# Human Organoid PCA-LR Reconstructed Feature-Effect Pilot",
        "",
        "Status: diagnostic feature-effect pilot; not a leaderboard result.",
        "",
        "This report generates `feature_effect.csv.gz` from source-transfer",
        "training samples only. Feature effects are PCA-LR coefficients",
        "reconstructed into train-standardized gene space using",
        "`pca.components_.T @ logistic.coef_[0]`. They are not log2FC",
        "response signatures and must not be interpreted as expression fold",
        "changes.",
        "",
        "The DE reference table is not used during feature-effect generation;",
        "it is used only afterward for rank, sign, and calibrated top-k",
        "diagnostic scoring.",
        "",
        "The accompanying `predictions.csv` is a two-row evaluator plumbing",
        "fixture. Classification metrics in this report are not interpretable",
        "as organoid model performance; only the feature-effect diagnostics are",
        "reviewed.",
        "",
        f"- classifier model: `{metadata['classifier_model_id']}`",
        f"- feature-effect rows: {metadata['n_effect_rows']}",
        f"- reference usage policy: `{metadata['reference_usage_policy']}`",
        "- claim boundary: feature-effect diagnostic only, not a biological",
        "  generalization claim and not a leaderboard claim.",
        "",
        "## Feature-Effect Metrics",
        "",
        "| Metric | PCA-LR Reconstructed | L2 Logistic Reference |",
        "|---|---:|---:|",
        (
            "| `feature_effect_direction_match` | "
            f"{_format_metric(_metric_value(evaluation, 'feature_effect_direction_match'))} | "
            f"{_format_metric(_comparison_value(comparison_metrics, 'feature_effect_direction_match'))} |"
        ),
        (
            "| `feature_effect_rank_correlation` | "
            f"{_format_metric(_metric_value(evaluation, 'feature_effect_rank_correlation'))} | "
            f"{_format_metric(_comparison_value(comparison_metrics, 'feature_effect_rank_correlation'))} |"
        ),
        "",
        "## Top-K DE Overlap",
        "",
        "| K | Overlap | Fraction | Expected | Enrichment | Hypergeometric p>=x |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for k, overlap, fraction, expected, enrichment, p_value in _top_k_summary(evaluation):
        lines.append(f"| {k} | {overlap} | {fraction} | {expected} | {enrichment} | {p_value} |")
    lines.append("")
    note_path.write_text("\n".join(lines))
    return note_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root for loading draft organoid matrices.",
    )
    parser.add_argument(
        "--task-manifest",
        default="v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json",
        help="Draft human organoid task manifest JSON.",
    )
    parser.add_argument(
        "--de-reference-manifest",
        default="v9/human_organoid/de_references/human_organoid_de_reference_manifest.draft.json",
        help="Derived DE reference manifest JSON used only for target contrast ids.",
    )
    parser.add_argument(
        "--reference-signature-table",
        default="v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz",
        help="Derived DE reference CSV/CSV.GZ used only for post hoc scoring.",
    )
    parser.add_argument(
        "--comparison-metrics",
        default="v9/human_organoid/reports/logistic_feature_effect/metrics.json",
        help="L2 logistic feature-effect metrics JSON for README comparison.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/human_organoid/reports/pca_lr_feature_effect",
        help="Output directory for PCA-LR reconstructed feature-effect report.",
    )
    parser.add_argument(
        "--top-variable-genes",
        type=int,
        default=2000,
        help="Number of train-fold variable genes to use per source-transfer model.",
    )
    parser.add_argument(
        "--pca-components",
        type=int,
        default=50,
        help="Requested PCA components. Adapted per source-transfer fold.",
    )
    parser.add_argument(
        "--transform",
        choices=["none", "log1p"],
        default="log1p",
        help="Feature transform before train-fold scaling.",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=5000,
        help="LogisticRegression max_iter.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Deterministic random_state.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    task = load_human_organoid_task(
        manifest_path=args.task_manifest,
        repo_root=args.repo_root,
    )
    task_manifest = load_task_manifest(args.task_manifest)
    result = build_pca_lr_reconstructed_feature_effect(
        task,
        de_reference_manifest_path=args.de_reference_manifest,
        transform=args.transform,
        top_variable_genes=args.top_variable_genes,
        pca_components=args.pca_components,
        max_iter=args.max_iter,
        random_state=args.random_state,
    )
    predictions_path = output_dir / "predictions.csv"
    feature_effect_path = output_dir / "feature_effect.csv.gz"
    metadata_path = output_dir / "feature_effect_metadata.json"
    _write_rows(predictions_path, PREDICTION_ROWS)
    _, written_metadata_path = write_feature_effect_artifact(
        result,
        feature_effect_path=feature_effect_path,
        metadata_path=metadata_path,
    )
    metadata = json.loads(written_metadata_path.read_text())
    evaluation = evaluate_submission(task_manifest, predictions_path)
    feature_metrics = compute_feature_effect_metrics(
        feature_effect_path=feature_effect_path,
        reference_signature_path=args.reference_signature_table,
        task_id=str(task_manifest["task_id"]),
    )
    evaluation["metrics"].update(feature_metrics["metrics"])
    evaluation["feature_effect_validation"] = feature_metrics["feature_effect_validation"]
    evaluation["feature_effect_path"] = feature_effect_path.as_posix()
    evaluation["reference_signature_path"] = Path(args.reference_signature_table).as_posix()
    evaluation["feature_effect_baseline"] = {
        "classifier_model_id": metadata["classifier_model_id"],
        "claim_boundary": metadata["claim_boundary"],
        "reference_usage_policy": metadata["reference_usage_policy"],
        "metadata_path": written_metadata_path.as_posix(),
        "effect_scale": "pca_lr_reconstructed_standardized_logistic_coefficient",
        "model_space": "reconstructed_gene_space_from_pca",
        "comparison_metrics_path": Path(args.comparison_metrics).as_posix(),
    }
    outputs = write_evaluation_report(
        evaluation_result=evaluation,
        task_manifest=task_manifest,
        task_manifest_path=args.task_manifest,
        submission_path=predictions_path,
        output_dir=output_dir,
        command=sys.argv,
    )
    note_path = _write_note(
        output_dir=output_dir,
        metadata=metadata,
        evaluation=evaluation,
        comparison_metrics=_load_comparison_metrics(args.comparison_metrics),
    )
    print(predictions_path)
    print(feature_effect_path)
    print(written_metadata_path)
    print(outputs["metrics"])
    print(outputs["run_manifest"])
    print(note_path)


if __name__ == "__main__":
    main()

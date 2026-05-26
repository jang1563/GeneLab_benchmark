#!/usr/bin/env python3
"""Run the v9 human organoid microglia-matched response-signature baseline."""

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
    build_microglia_matched_source_transfer_response_signature,
    evaluate_submission,
    load_human_organoid_task,
    load_task_manifest,
    write_microglia_matched_source_transfer_response_signature,
)
from spacebio_bench.reports import write_evaluation_report  # noqa: E402


PREDICTION_ROWS = [
    {
        "sample_id": "microglia_source_transfer_smoke_ground",
        "true_label": "Ground",
        "predicted_label": "Ground",
        "flight_probability": "0.05",
    },
    {
        "sample_id": "microglia_source_transfer_smoke_leo_or_iss",
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


def _metric_value(metrics: dict[str, object], metric_id: str) -> float | None:
    metric = metrics.get("metrics", {}).get(metric_id) if isinstance(metrics.get("metrics"), dict) else None
    if not isinstance(metric, dict) or metric.get("status") != "computed":
        return None
    value = metric.get("value")
    return float(value) if isinstance(value, (int, float)) else None


def _format_metric(value: float | None) -> str:
    if value is None:
        return "NA"
    return f"{value:.6f}"


def _comparison_rows(
    *,
    evaluation: dict[str, object],
    global_metrics_path: Path,
) -> list[tuple[str, float | None, float | None, float | None]]:
    global_metrics: dict[str, object] = {}
    if global_metrics_path.exists():
        global_metrics = json.loads(global_metrics_path.read_text())
    rows = []
    for metric_id in ["de_direction_match", "signature_rank_correlation"]:
        global_value = _metric_value(global_metrics, metric_id) if global_metrics else None
        microglia_value = _metric_value(evaluation, metric_id)
        delta = None
        if global_value is not None and microglia_value is not None:
            delta = microglia_value - global_value
        rows.append((metric_id, global_value, microglia_value, delta))
    return rows


def _write_note(
    *,
    output_dir: Path,
    metadata: dict[str, object],
    evaluation: dict[str, object],
    global_metrics_path: Path,
) -> Path:
    note_path = output_dir / "README.md"
    comparison = _comparison_rows(
        evaluation=evaluation,
        global_metrics_path=global_metrics_path,
    )
    lines = [
        "# Human Organoid Microglia-Matched Source-Transfer Response-Signature Baseline",
        "",
        "Status: diagnostic baseline; not a leaderboard result.",
        "",
        "This report generates `response_signature.csv.gz` from source-transfer",
        "training samples matched to each target contrast's microglia condition.",
        "The DE reference table is not used during signature generation; it is",
        "used only afterward for diagnostic scoring by the evaluator.",
        "",
        "The accompanying `predictions.csv` is a two-row evaluator plumbing",
        "fixture. Classification metrics in this report are not interpretable",
        "as organoid model performance; only the response-signature diagnostics",
        "are reviewed.",
        "",
        f"- signature model: `{metadata['signature_model_id']}`",
        f"- response rows: {metadata['n_response_rows']}",
        f"- conditioning strategy: `{metadata['conditioning_strategy']}`",
        f"- reference usage policy: `{metadata['reference_usage_policy']}`",
        "- claim boundary: microglia-matched source-transfer diagnostic only,",
        "  not a biological generalization claim and not a leaderboard claim.",
        "",
        "## Comparison To Global Source-Transfer",
        "",
        "| Metric | Global Source-Transfer | Microglia-Matched | Delta |",
        "|---|---:|---:|---:|",
    ]
    for metric_id, global_value, microglia_value, delta in comparison:
        lines.append(
            "| "
            + " | ".join(
                [
                    f"`{metric_id}`",
                    _format_metric(global_value),
                    _format_metric(microglia_value),
                    _format_metric(delta),
                ]
            )
            + " |"
        )
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
        help="Derived DE reference manifest JSON used only for target contrast metadata.",
    )
    parser.add_argument(
        "--reference-signature-table",
        default="v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz",
        help="Derived DE reference CSV/CSV.GZ used only for post hoc scoring.",
    )
    parser.add_argument(
        "--global-source-transfer-metrics",
        default="v9/human_organoid/reports/source_transfer_signature/metrics.json",
        help="Global source-transfer metrics JSON for comparison in README.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/human_organoid/reports/microglia_source_transfer_signature",
        help="Output directory for microglia-matched response-signature report.",
    )
    parser.add_argument(
        "--max-features",
        type=int,
        default=None,
        help="Optional feature limit for smoke tests; omit for the full report.",
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
    result = build_microglia_matched_source_transfer_response_signature(
        task,
        de_reference_manifest_path=args.de_reference_manifest,
        max_features=args.max_features,
    )
    predictions_path = output_dir / "predictions.csv"
    response_path = output_dir / "response_signature.csv.gz"
    metadata_path = output_dir / "response_signature_metadata.json"
    _write_rows(predictions_path, PREDICTION_ROWS)
    _, written_metadata_path = write_microglia_matched_source_transfer_response_signature(
        result,
        response_signature_path=response_path,
        metadata_path=metadata_path,
    )
    metadata = json.loads(written_metadata_path.read_text())
    evaluation = evaluate_submission(
        task_manifest,
        predictions_path,
        response_signature_path=response_path,
        reference_signature_path=args.reference_signature_table,
    )
    evaluation["response_signature_baseline"] = {
        "signature_model_id": metadata["signature_model_id"],
        "claim_boundary": metadata["claim_boundary"],
        "conditioning_strategy": metadata["conditioning_strategy"],
        "reference_usage_policy": metadata["reference_usage_policy"],
        "metadata_path": written_metadata_path.as_posix(),
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
        global_metrics_path=Path(args.global_source_transfer_metrics),
    )
    print(predictions_path)
    print(response_path)
    print(written_metadata_path)
    print(outputs["metrics"])
    print(outputs["run_manifest"])
    print(note_path)


if __name__ == "__main__":
    main()

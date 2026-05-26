#!/usr/bin/env python3
"""Run a diagnostic response-signature smoke test for the v9 organoid task."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import evaluate_submission, load_task_manifest  # noqa: E402
from spacebio_bench.reports import write_evaluation_report  # noqa: E402


PREDICTION_ROWS = [
    {
        "sample_id": "smoke_ground",
        "true_label": "Ground",
        "predicted_label": "Ground",
        "flight_probability": "0.05",
    },
    {
        "sample_id": "smoke_leo_or_iss",
        "true_label": "LEO_or_ISS",
        "predicted_label": "LEO_or_ISS",
        "flight_probability": "0.95",
    },
]


def _open_text(path: str | Path):
    parsed = Path(path)
    if parsed.name.endswith(".gz"):
        return gzip.open(parsed, "rt", newline="")
    return parsed.open(newline="")


def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _select_reference_rows(
    *,
    reference_path: str | Path,
    significant_per_contrast: int,
    background_per_contrast: int,
) -> list[dict[str, str]]:
    significant: dict[str, list[dict[str, str]]] = {}
    background: dict[str, list[dict[str, str]]] = {}
    with _open_text(reference_path) as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            contrast_id = str(row.get("contrast_id", "") or "")
            if not contrast_id:
                continue
            is_significant = str(row.get("significant_fdr_0_05", "")).lower() == "true"
            bucket = significant if is_significant else background
            limit = significant_per_contrast if is_significant else background_per_contrast
            rows = bucket.setdefault(contrast_id, [])
            if len(rows) < limit:
                rows.append(dict(row))

    selected: list[dict[str, str]] = []
    for contrast_id in sorted(set(significant) | set(background)):
        selected.extend(significant.get(contrast_id, []))
        selected.extend(background.get(contrast_id, []))
    if not selected:
        raise ValueError(f"no reference rows selected from {reference_path}")
    return selected


def _response_rows_from_reference(
    reference_rows: list[dict[str, str]],
    *,
    task_id: str,
) -> list[dict[str, str]]:
    response_rows: list[dict[str, str]] = []
    for row in reference_rows:
        response_rows.append(
            {
                "task_id": task_id,
                "source_id": str(row.get("source_id", "") or ""),
                "contrast_id": str(row.get("contrast_id", "") or ""),
                "gene_symbol": str(row.get("gene_symbol", "") or ""),
                "ensembl_id": str(row.get("ensembl_id", "") or ""),
                "predicted_log2fc_leo_or_iss_minus_ground": str(
                    row.get("log2fc_leo_or_iss_minus_ground", "") or ""
                ),
            }
        )
    return response_rows


def _write_note(
    *,
    output_dir: Path,
    reference_path: str | Path,
    n_response_rows: int,
) -> Path:
    note_path = output_dir / "README.md"
    note_path.write_text(
        "\n".join(
            [
                "# Human Organoid Response-Signature Smoke Test",
                "",
                "Status: scorer smoke test; not a biological model result.",
                "",
                "This directory exercises the response-signature evaluator against the",
                "real derived human organoid DE reference. The response signature is",
                "constructed directly from selected reference rows, so perfect metric",
                "values only mean that the evaluator join/scoring path is functioning.",
                "",
                f"- reference: `{Path(reference_path).as_posix()}`",
                f"- response rows: {n_response_rows}",
                "- claim boundary: diagnostic scorer plumbing only, not leaderboard",
                "  evidence and not a model performance claim.",
                "",
            ]
        )
    )
    return note_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--task-manifest",
        default="v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json",
        help="Draft human organoid task manifest JSON.",
    )
    parser.add_argument(
        "--reference-signature-table",
        default="v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz",
        help="Derived DE reference CSV or CSV.GZ path.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/human_organoid/reports/response_signature_smoke",
        help="Output directory for the smoke-test report.",
    )
    parser.add_argument(
        "--significant-per-contrast",
        type=int,
        default=3,
        help="Number of significant reference genes to mirror per contrast.",
    )
    parser.add_argument(
        "--background-per-contrast",
        type=int,
        default=2,
        help="Number of non-significant reference genes to mirror per contrast.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    task_manifest = load_task_manifest(args.task_manifest)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    predictions_path = output_dir / "predictions.csv"
    response_path = output_dir / "response_signature.csv"

    _write_rows(predictions_path, PREDICTION_ROWS)
    reference_rows = _select_reference_rows(
        reference_path=args.reference_signature_table,
        significant_per_contrast=args.significant_per_contrast,
        background_per_contrast=args.background_per_contrast,
    )
    response_rows = _response_rows_from_reference(
        reference_rows,
        task_id=str(task_manifest["task_id"]),
    )
    _write_rows(response_path, response_rows)
    _write_note(
        output_dir=output_dir,
        reference_path=args.reference_signature_table,
        n_response_rows=len(response_rows),
    )

    result = evaluate_submission(
        task_manifest,
        predictions_path,
        response_signature_path=response_path,
        reference_signature_path=args.reference_signature_table,
    )
    outputs = write_evaluation_report(
        evaluation_result=result,
        task_manifest=task_manifest,
        task_manifest_path=args.task_manifest,
        submission_path=predictions_path,
        output_dir=output_dir,
        command=sys.argv,
    )
    print(predictions_path)
    print(response_path)
    print(outputs["metrics"])
    print(outputs["run_manifest"])


if __name__ == "__main__":
    main()

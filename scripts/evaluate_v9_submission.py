#!/usr/bin/env python3
"""Validate and evaluate a v9 SpaceBio-Bench prediction CSV."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import evaluate_submission, load_task_manifest
from spacebio_bench.reports import write_evaluation_report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-manifest", required=True, help="v9 task manifest JSON.")
    parser.add_argument("--submission", required=True, help="Prediction CSV to evaluate.")
    parser.add_argument("--output-dir", required=True, help="Directory for evaluation outputs.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    task_manifest = load_task_manifest(args.task_manifest)
    result = evaluate_submission(task_manifest, args.submission)
    outputs = write_evaluation_report(
        evaluation_result=result,
        task_manifest=task_manifest,
        task_manifest_path=args.task_manifest,
        submission_path=args.submission,
        output_dir=args.output_dir,
        command=sys.argv,
    )
    print(outputs["metrics"])
    print(outputs["run_manifest"])
    if result["status"] == "invalid":
        raise SystemExit(1)


if __name__ == "__main__":
    main()

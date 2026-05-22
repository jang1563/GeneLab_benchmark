#!/usr/bin/env python3
"""Run the v9 nearest-centroid baseline on bulk LOMO tasks."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.baselines import NearestCentroidConfig, run_nearest_centroid_all


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--task-id",
        action="append",
        help="Task id to run. May be repeated. Defaults to all bulk LOMO tasks.",
    )
    parser.add_argument(
        "--manifest-dir",
        default="v9/task_manifests",
        help="Directory containing v9 task manifests.",
    )
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve legacy task data paths.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/reports/nearest_centroid",
        help="Directory for predictions, metrics, run manifests, and summary files.",
    )
    parser.add_argument(
        "--scaling",
        choices=["none", "zscore"],
        default="none",
        help="Feature scaling mode fit on each training fold only.",
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop on the first task failure instead of writing failure rows.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows, outputs = run_nearest_centroid_all(
        task_ids=args.task_id,
        manifest_dir=args.manifest_dir,
        repo_root=args.repo_root,
        output_root=args.output_dir,
        config=NearestCentroidConfig(scaling=args.scaling),
        command=sys.argv,
        fail_fast=args.fail_fast,
    )
    print(outputs["csv"])
    print(outputs["json"])
    failed = [row for row in rows if row["status"] == "failed"]
    if failed:
        for row in failed:
            print(f"FAILED {row['task_id']}: {row['error']}", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()

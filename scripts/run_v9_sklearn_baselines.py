#!/usr/bin/env python3
"""Run v9 sklearn classifier baselines on bulk LOMO tasks."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.baselines.sklearn_classifiers import run_sklearn_baseline_all


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--baseline",
        action="append",
        help=(
            "Baseline to run. May be repeated. Choices/aliases include "
            "pca_lr, pca_logistic_regression, logistic_l2, logistic_regression_l2. "
            "Defaults to PCA-LR and L2 logistic regression."
        ),
    )
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
        default="v9/reports/sklearn_baselines",
        help="Directory for predictions, metrics, run manifests, and summary files.",
    )
    parser.add_argument(
        "--pca-components",
        type=int,
        default=50,
        help="Requested PCA components for the PCA-LR baseline. Adapted per fold.",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=5000,
        help="Maximum logistic-regression iterations.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Deterministic sklearn random_state.",
    )
    parser.add_argument(
        "--c",
        type=float,
        default=1.0,
        help="Inverse regularization strength for LogisticRegression.",
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Stop on the first task failure instead of writing failure rows.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows, outputs = run_sklearn_baseline_all(
        baseline_ids=args.baseline,
        task_ids=args.task_id,
        manifest_dir=args.manifest_dir,
        repo_root=args.repo_root,
        output_root=args.output_dir,
        pca_components=args.pca_components,
        max_iter=args.max_iter,
        random_state=args.random_state,
        c=args.c,
        command=sys.argv,
        fail_fast=args.fail_fast,
    )
    print(outputs["csv"])
    print(outputs["json"])
    failed = [row for row in rows if row["status"] == "failed"]
    if failed:
        for row in failed:
            print(
                f"FAILED {row['baseline_id']} {row['task_id']}: {row['error']}",
                file=sys.stderr,
            )
        raise SystemExit(1)


if __name__ == "__main__":
    main()

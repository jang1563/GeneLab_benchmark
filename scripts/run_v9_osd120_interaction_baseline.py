#!/usr/bin/env python3
"""Run draft v9 OSD-120 interaction nearest-centroid baselines."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.baselines import (  # noqa: E402
    INTERACTION_FOLD_FAMILIES,
    MultispeciesBaselineConfig,
    run_multispecies_interaction_baselines,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest-dir",
        default="v9/multispecies/interaction_task_manifests",
        help="Directory containing draft multispecies interaction task manifests.",
    )
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve manifest-relative data paths.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/multispecies/reports/interaction_nearest_centroid",
        help="Directory for predictions, metrics, run manifests, and summary.",
    )
    parser.add_argument(
        "--fold-family",
        choices=INTERACTION_FOLD_FAMILIES,
        action="append",
        help="Fold family to evaluate. May be repeated. Defaults to all families.",
    )
    parser.add_argument(
        "--top-variable-genes",
        type=int,
        default=2000,
        help="Number of train-fold variable genes to use per draft fold.",
    )
    parser.add_argument(
        "--transform",
        choices=["none", "log1p"],
        default="log1p",
        help="Feature transform applied before train-fold scaling.",
    )
    parser.add_argument(
        "--scaling",
        choices=["none", "zscore"],
        default="zscore",
        help="Feature scaling fit on each training fold only.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    _, summary = run_multispecies_interaction_baselines(
        manifest_dir=args.manifest_dir,
        repo_root=args.repo_root,
        output_root=args.output_dir,
        fold_families=args.fold_family,
        config=MultispeciesBaselineConfig(
            transform=args.transform,
            scaling=args.scaling,
            top_variable_genes=args.top_variable_genes,
        ),
        command=sys.argv,
    )
    print(summary["csv"])
    print(summary["json"])


if __name__ == "__main__":
    main()

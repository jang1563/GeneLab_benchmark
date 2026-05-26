#!/usr/bin/env python3
"""Run sensitivity variants for draft v9 OSD-120 interaction baselines."""

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
    run_multispecies_interaction_sensitivity_grid,
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
        default="v9/multispecies/reports/interaction_sensitivity",
        help="Directory for per-variant predictions, metrics, run manifests, and summary.",
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
        action="append",
        help="Feature counts to test. May be repeated. Defaults to 100, 500, 2000, 5000, and all features.",
    )
    parser.add_argument(
        "--transform",
        choices=["none", "log1p"],
        action="append",
        help="Transforms to test. May be repeated. Defaults to log1p and none.",
    )
    parser.add_argument(
        "--scaling",
        choices=["none", "zscore"],
        action="append",
        help="Scaling modes to test. May be repeated. Defaults to zscore and none.",
    )
    return parser.parse_args()


def _configs(args: argparse.Namespace) -> list[MultispeciesBaselineConfig] | None:
    if not args.top_variable_genes and not args.transform and not args.scaling:
        return None
    tops = args.top_variable_genes or [100, 500, 2000, 5000]
    transforms = args.transform or ["log1p", "none"]
    scalings = args.scaling or ["zscore", "none"]
    return [
        MultispeciesBaselineConfig(
            transform=transform,
            scaling=scaling,
            top_variable_genes=top,
        )
        for transform in transforms
        for scaling in scalings
        for top in tops
    ]


def main() -> None:
    args = parse_args()
    _, summary = run_multispecies_interaction_sensitivity_grid(
        manifest_dir=args.manifest_dir,
        repo_root=args.repo_root,
        output_root=args.output_dir,
        configs=_configs(args),
        fold_families=args.fold_family,
        command=sys.argv,
    )
    print(summary["csv"])
    print(summary["json"])


if __name__ == "__main__":
    main()

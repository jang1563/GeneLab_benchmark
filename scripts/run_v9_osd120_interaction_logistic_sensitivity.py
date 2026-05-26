#!/usr/bin/env python3
"""Run compact v9 OSD-120 interaction L2 logistic sensitivity grid."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.baselines import (  # noqa: E402
    INTERACTION_FOLD_FAMILIES,
    MultispeciesLogisticBaselineConfig,
    run_multispecies_interaction_logistic_sensitivity_grid,
    write_multispecies_interaction_fold_detail_comparison,
    write_multispecies_interaction_fold_details,
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
        default="v9/multispecies/reports/interaction_logistic_l2_sensitivity",
        help="Directory for predictions, metrics, summaries, and comparison tables.",
    )
    parser.add_argument(
        "--nearest-centroid-fold-detail-csv",
        default="v9/multispecies/reports/interaction_sensitivity/fold_detail_summary.csv",
        help="Nearest-centroid fold-detail CSV to compare against.",
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
        help="Top train-fold variable-gene count. May be repeated.",
    )
    parser.add_argument(
        "--c",
        type=float,
        action="append",
        help="Inverse L2 regularization strength. May be repeated.",
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
    parser.add_argument(
        "--max-iter",
        type=int,
        default=5000,
        help="Maximum LogisticRegression iterations.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Deterministic LogisticRegression random_state.",
    )
    return parser.parse_args()


def configs_from_args(args: argparse.Namespace) -> list[MultispeciesLogisticBaselineConfig] | None:
    if args.top_variable_genes is None and args.c is None:
        return None
    top_values = args.top_variable_genes or [500, 2000]
    c_values = args.c or [0.1, 1.0, 10.0]
    return [
        MultispeciesLogisticBaselineConfig(
            transform=args.transform,
            scaling=args.scaling,
            top_variable_genes=top,
            c=c,
            max_iter=args.max_iter,
            random_state=args.random_state,
        )
        for top in top_values
        for c in c_values
    ]


def main() -> None:
    args = parse_args()
    _, summary = run_multispecies_interaction_logistic_sensitivity_grid(
        manifest_dir=args.manifest_dir,
        repo_root=args.repo_root,
        output_root=args.output_dir,
        configs=configs_from_args(args),
        fold_families=args.fold_family,
        command=sys.argv,
    )
    fold_outputs = write_multispecies_interaction_fold_details(
        summary_csv=summary["csv"],
        repo_root=args.repo_root,
        output_dir=args.output_dir,
    )
    comparison_outputs = write_multispecies_interaction_fold_detail_comparison(
        output_dir=args.output_dir,
        nearest_centroid_fold_detail_csv=args.nearest_centroid_fold_detail_csv,
        logistic_fold_detail_csv=fold_outputs["csv"],
        repo_root=args.repo_root,
        logistic_variant_id=None,
    )
    print(summary["csv"])
    print(summary["json"])
    print(fold_outputs["csv"])
    print(fold_outputs["json"])
    print(comparison_outputs["csv"])
    print(comparison_outputs["json"])


if __name__ == "__main__":
    main()

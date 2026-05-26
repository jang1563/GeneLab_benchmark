#!/usr/bin/env python3
"""Run draft v9 OSD-120 interaction sparse L1 logistic pilot."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.baselines import (  # noqa: E402
    DEFAULT_LOGISTIC_FEATURE_AUDIT_FOLDS,
    INTERACTION_FOLD_FAMILIES,
    MultispeciesLogisticBaselineConfig,
    run_multispecies_interaction_sparse_logistic_l1_grid,
    write_multispecies_interaction_fold_detail_comparison,
    write_multispecies_interaction_fold_details,
    write_multispecies_interaction_logistic_feature_audit,
)


DEFAULT_C_VALUES = [0.01, 0.03, 0.1, 0.3, 1.0]


def parse_focus_fold(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("focus fold must be fold_family=heldout_value")
    fold_family, heldout_value = value.split("=", 1)
    if not fold_family or not heldout_value:
        raise argparse.ArgumentTypeError("focus fold must be fold_family=heldout_value")
    return fold_family, heldout_value


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
        default="v9/multispecies/reports/interaction_logistic_sparse_l1",
        help="Directory for sparse L1 logistic reports.",
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
        "--focus-fold",
        type=parse_focus_fold,
        action="append",
        help="Fold to include in sparse coefficient audit.",
    )
    parser.add_argument(
        "--top-variable-genes",
        type=int,
        action="append",
        help="Top train-fold variable-gene count. Defaults to 2000.",
    )
    parser.add_argument(
        "--c",
        type=float,
        action="append",
        help="Inverse L1 regularization strength. May be repeated.",
    )
    parser.add_argument(
        "--top-n-coefficients",
        type=int,
        default=10,
        help="Number of top coefficient feature ids to include in audit summaries.",
    )
    return parser.parse_args()


def configs_from_args(args: argparse.Namespace) -> list[MultispeciesLogisticBaselineConfig] | None:
    if args.top_variable_genes is None and args.c is None:
        return None
    top_values = args.top_variable_genes or [2000]
    c_values = args.c or DEFAULT_C_VALUES
    return [
        MultispeciesLogisticBaselineConfig(
            top_variable_genes=top,
            c=c,
            penalty="l1",
            solver="liblinear",
        )
        for top in top_values
        for c in c_values
    ]


def main() -> None:
    args = parse_args()
    configs = configs_from_args(args)
    audit_configs = configs or [
        MultispeciesLogisticBaselineConfig(
            top_variable_genes=2000,
            c=c,
            penalty="l1",
            solver="liblinear",
        )
        for c in DEFAULT_C_VALUES
    ]
    _, summary = run_multispecies_interaction_sparse_logistic_l1_grid(
        manifest_dir=args.manifest_dir,
        repo_root=args.repo_root,
        output_root=args.output_dir,
        configs=configs,
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
    audit_outputs = write_multispecies_interaction_logistic_feature_audit(
        manifest_dir=args.manifest_dir,
        repo_root=args.repo_root,
        output_dir=args.output_dir,
        configs=audit_configs,
        focus_folds=args.focus_fold or DEFAULT_LOGISTIC_FEATURE_AUDIT_FOLDS,
        top_n_coefficients=args.top_n_coefficients,
    )
    print(summary["csv"])
    print(summary["json"])
    print(fold_outputs["csv"])
    print(fold_outputs["json"])
    print(comparison_outputs["csv"])
    print(comparison_outputs["json"])
    print(audit_outputs["summary_csv"])
    print(audit_outputs["summary_json"])
    print(audit_outputs["coefficient_csv"])
    print(audit_outputs["coefficient_json"])


if __name__ == "__main__":
    main()

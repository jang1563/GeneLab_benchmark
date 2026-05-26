#!/usr/bin/env python3
"""Audit selected feature sets for OSD-120 interaction L2 logistic variants."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.baselines import (  # noqa: E402
    DEFAULT_LOGISTIC_FEATURE_AUDIT_FOLDS,
    MultispeciesLogisticBaselineConfig,
    write_multispecies_interaction_logistic_feature_audit,
)


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
        default="v9/multispecies/reports/interaction_logistic_feature_audit",
        help="Directory for feature-set and coefficient audit outputs.",
    )
    parser.add_argument(
        "--focus-fold",
        type=parse_focus_fold,
        action="append",
        help=(
            "Fold to audit, formatted as fold_family=heldout_value. "
            "May be repeated. Defaults to the V9-MULTI-017 focus folds."
        ),
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
        default=1.0,
        help="Inverse L2 regularization strength for the audited variants.",
    )
    parser.add_argument(
        "--top-n-coefficients",
        type=int,
        default=10,
        help="Number of top coefficient feature ids to include in summary columns.",
    )
    return parser.parse_args()


def configs_from_args(args: argparse.Namespace) -> list[MultispeciesLogisticBaselineConfig] | None:
    if args.top_variable_genes is None:
        return None
    return [
        MultispeciesLogisticBaselineConfig(top_variable_genes=top, c=args.c)
        for top in args.top_variable_genes
    ]


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_logistic_feature_audit(
        manifest_dir=args.manifest_dir,
        repo_root=args.repo_root,
        output_dir=args.output_dir,
        configs=configs_from_args(args),
        focus_folds=args.focus_fold or DEFAULT_LOGISTIC_FEATURE_AUDIT_FOLDS,
        top_n_coefficients=args.top_n_coefficients,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["coefficient_csv"])
    print(outputs["coefficient_json"])


if __name__ == "__main__":
    main()

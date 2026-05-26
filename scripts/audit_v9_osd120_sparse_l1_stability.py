#!/usr/bin/env python3
"""Audit sparse L1 logistic feature stability for OSD-120 interaction folds."""

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
    write_multispecies_interaction_sparse_l1_stability_audit,
)


DEFAULT_C_VALUES = [0.3, 1.0]


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
        default="v9/multispecies/reports/interaction_logistic_sparse_l1_stability",
        help="Directory for sparse L1 stability audit outputs.",
    )
    parser.add_argument(
        "--focus-fold",
        type=parse_focus_fold,
        action="append",
        help="Fold to audit, formatted as fold_family=heldout_value.",
    )
    parser.add_argument(
        "--top-variable-genes",
        type=int,
        default=2000,
        help="Top train-fold variable-gene count for sparse L1 candidates.",
    )
    parser.add_argument(
        "--c",
        type=float,
        action="append",
        help="Inverse L1 regularization strength. May be repeated.",
    )
    parser.add_argument(
        "--n-subsamples",
        type=int,
        default=20,
        help="Number of deterministic balanced subsamples per fold/config.",
    )
    parser.add_argument(
        "--subsample-fraction",
        type=float,
        default=0.8,
        help="Within-class training fraction for each deterministic subsample.",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=1729,
        help="Base seed for deterministic subsample generation.",
    )
    return parser.parse_args()


def configs_from_args(args: argparse.Namespace) -> list[MultispeciesLogisticBaselineConfig]:
    return [
        MultispeciesLogisticBaselineConfig(
            top_variable_genes=args.top_variable_genes,
            c=c,
            penalty="l1",
            solver="liblinear",
        )
        for c in (args.c or DEFAULT_C_VALUES)
    ]


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_sparse_l1_stability_audit(
        manifest_dir=args.manifest_dir,
        repo_root=args.repo_root,
        output_dir=args.output_dir,
        configs=configs_from_args(args),
        focus_folds=args.focus_fold or DEFAULT_LOGISTIC_FEATURE_AUDIT_FOLDS,
        n_subsamples=args.n_subsamples,
        subsample_fraction=args.subsample_fraction,
        random_seed=args.random_seed,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["feature_csv"])
    print(outputs["feature_json"])


if __name__ == "__main__":
    main()

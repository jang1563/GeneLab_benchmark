#!/usr/bin/env python3
"""Build OSD-120 logistic fold details and nearest-centroid comparison tables."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    write_multispecies_interaction_fold_detail_comparison,
    write_multispecies_interaction_fold_details,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve relative metrics paths.",
    )
    parser.add_argument(
        "--nearest-centroid-fold-detail-csv",
        default="v9/multispecies/reports/interaction_sensitivity/fold_detail_summary.csv",
        help="Nearest-centroid fold-detail CSV to compare against.",
    )
    parser.add_argument(
        "--logistic-summary-csv",
        default="v9/multispecies/reports/interaction_logistic_l2/multispecies_baseline_summary.csv",
        help="Logistic interaction summary CSV with metrics paths.",
    )
    parser.add_argument(
        "--logistic-output-dir",
        default="v9/multispecies/reports/interaction_logistic_l2",
        help="Directory for logistic fold-detail and comparison outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logistic_outputs = write_multispecies_interaction_fold_details(
        summary_csv=args.logistic_summary_csv,
        repo_root=args.repo_root,
        output_dir=args.logistic_output_dir,
    )
    comparison_outputs = write_multispecies_interaction_fold_detail_comparison(
        output_dir=args.logistic_output_dir,
        nearest_centroid_fold_detail_csv=args.nearest_centroid_fold_detail_csv,
        logistic_fold_detail_csv=logistic_outputs["csv"],
        repo_root=args.repo_root,
    )
    print(logistic_outputs["csv"])
    print(logistic_outputs["json"])
    print(comparison_outputs["csv"])
    print(comparison_outputs["json"])


if __name__ == "__main__":
    main()

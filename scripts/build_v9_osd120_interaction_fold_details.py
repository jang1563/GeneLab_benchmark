#!/usr/bin/env python3
"""Build fold-detail summaries for OSD-120 interaction sensitivity reports."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import write_multispecies_interaction_fold_details  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--summary-csv",
        default="v9/multispecies/reports/interaction_sensitivity/multispecies_baseline_summary.csv",
        help="Interaction sensitivity summary CSV with metrics paths.",
    )
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve relative metrics paths.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/multispecies/reports/interaction_sensitivity",
        help="Directory for fold_detail_summary.csv and .json.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_fold_details(
        summary_csv=args.summary_csv,
        repo_root=args.repo_root,
        output_dir=args.output_dir,
    )
    print(outputs["csv"])
    print(outputs["json"])


if __name__ == "__main__":
    main()

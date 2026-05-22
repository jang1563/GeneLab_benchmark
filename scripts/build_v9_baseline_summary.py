#!/usr/bin/env python3
"""Build a cross-baseline v9 bulk LOMO summary table."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.reports import read_baseline_summary, write_baseline_summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--nearest-centroid-summary",
        default="v9/reports/nearest_centroid/bulk_lomo_summary.csv",
        help="Nearest-centroid summary CSV. This file predates baseline_id columns.",
    )
    parser.add_argument(
        "--sklearn-summary",
        default="v9/reports/sklearn_baselines/bulk_lomo_summary.csv",
        help="sklearn baseline summary CSV.",
    )
    parser.add_argument(
        "--csv",
        default="v9/reports/bulk_lomo_baseline_summary.csv",
        help="Output cross-baseline CSV path.",
    )
    parser.add_argument(
        "--json",
        default="v9/reports/bulk_lomo_baseline_summary.json",
        help="Output cross-baseline JSON path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = []
    nearest_path = Path(args.nearest_centroid_summary)
    sklearn_path = Path(args.sklearn_summary)
    if nearest_path.exists():
        rows.extend(read_baseline_summary(nearest_path, baseline_id="nearest_centroid"))
    if sklearn_path.exists():
        rows.extend(read_baseline_summary(sklearn_path))
    if not rows:
        raise SystemExit("no baseline summary rows found")

    outputs = write_baseline_summary(rows, csv_path=args.csv, json_path=args.json)
    print(outputs["csv"])
    print(outputs["json"])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build parsed sample-factor rows for the v9 human organoid draft task."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.sample_table_audit import (  # noqa: E402
    build_sample_factor_rows,
    write_sample_factor_table,
)
from spacebio_bench.source_audit import BIODATA_API_BASE, read_source_inventory  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-inventory",
        default="v9/human_organoid/source_inventory.draft.csv",
        help="Input human organoid source inventory CSV path.",
    )
    parser.add_argument(
        "--csv",
        default="v9/human_organoid/sample_factors.draft.csv",
        help="Output parsed sample-factor CSV path.",
    )
    parser.add_argument(
        "--json",
        default="v9/human_organoid/sample_factors.draft.json",
        help="Output parsed sample-factor JSON path.",
    )
    parser.add_argument(
        "--api-base",
        default=BIODATA_API_BASE,
        help="OSDR biological-data API base URL.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = read_source_inventory(args.source_inventory)
    factor_rows = build_sample_factor_rows(rows, api_base=args.api_base)
    csv_path, json_path = write_sample_factor_table(
        factor_rows,
        csv_path=args.csv,
        json_path=args.json,
    )
    print(csv_path)
    print(json_path)


if __name__ == "__main__":
    main()


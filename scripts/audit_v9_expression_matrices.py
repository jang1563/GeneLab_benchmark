#!/usr/bin/env python3
"""Audit/download normalized expression matrices for v9 draft extension sources."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.expression_audit import (  # noqa: E402
    audit_expression_matrix_inventory,
    write_expression_matrix_audit,
)
from spacebio_bench.source_audit import BIODATA_API_BASE, read_source_inventory  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-inventory",
        default="v9/human_organoid/source_inventory.draft.csv",
        help="Input source inventory CSV path.",
    )
    parser.add_argument(
        "--sample-factors",
        default="v9/human_organoid/sample_factors.draft.csv",
        help="Parsed sample-factor CSV path.",
    )
    parser.add_argument(
        "--matrix-dir",
        default="v9/human_organoid/matrices",
        help="Output directory for downloaded normalized-count matrices.",
    )
    parser.add_argument(
        "--csv",
        default="v9/human_organoid/expression_matrix_audit.draft.csv",
        help="Output expression matrix audit CSV path.",
    )
    parser.add_argument(
        "--json",
        default="v9/human_organoid/expression_matrix_audit.draft.json",
        help="Output expression matrix audit JSON path.",
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
    with Path(args.sample_factors).open(newline="") as handle:
        sample_factor_rows = [dict(row) for row in csv.DictReader(handle)]
    audit_rows = audit_expression_matrix_inventory(
        rows,
        sample_factor_rows=sample_factor_rows,
        output_dir=args.matrix_dir,
        api_base=args.api_base,
    )
    csv_path, json_path = write_expression_matrix_audit(
        audit_rows,
        csv_path=args.csv,
        json_path=args.json,
    )
    print(csv_path)
    print(json_path)


if __name__ == "__main__":
    main()


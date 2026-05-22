#!/usr/bin/env python3
"""Audit OSDR sample-table availability for v9 source inventories."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.sample_table_audit import (  # noqa: E402
    audit_sample_table_inventory,
    write_sample_table_audit,
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
        "--csv",
        default="v9/human_organoid/sample_table_audit.draft.csv",
        help="Output sample-table audit CSV path.",
    )
    parser.add_argument(
        "--json",
        default="v9/human_organoid/sample_table_audit.draft.json",
        help="Output sample-table audit JSON path.",
    )
    parser.add_argument(
        "--api-base",
        default=BIODATA_API_BASE,
        help="OSDR biological-data API base URL.",
    )
    parser.add_argument(
        "--source-id",
        action="append",
        default=[],
        help="Optional source accession to audit. May be repeated.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = read_source_inventory(args.source_inventory)
    if args.source_id:
        requested = set(args.source_id)
        rows = [row for row in rows if row["source_id"] in requested]
        missing = sorted(requested - {row["source_id"] for row in rows})
        if missing:
            raise SystemExit(f"source ids not found in inventory: {', '.join(missing)}")
    audit_rows = audit_sample_table_inventory(rows, api_base=args.api_base)
    csv_path, json_path = write_sample_table_audit(
        audit_rows,
        csv_path=args.csv,
        json_path=args.json,
    )
    print(csv_path)
    print(json_path)


if __name__ == "__main__":
    main()


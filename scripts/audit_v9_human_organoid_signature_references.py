#!/usr/bin/env python3
"""Audit public DE/signature reference evidence for the v9 human organoid pilot."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.organoid_signature_audit import (  # noqa: E402
    audit_organoid_signature_reference_inventory,
    write_organoid_signature_reference_audit,
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
        default="v9/human_organoid/signature_reference_audit.draft.csv",
        help="Output signature-reference audit CSV path.",
    )
    parser.add_argument(
        "--json",
        default="v9/human_organoid/signature_reference_audit.draft.json",
        help="Output signature-reference audit JSON path.",
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
    audit_rows = audit_organoid_signature_reference_inventory(
        rows,
        api_base=args.api_base,
    )
    csv_path, json_path = write_organoid_signature_reference_audit(
        audit_rows,
        csv_path=args.csv,
        json_path=args.json,
    )
    print(csv_path)
    print(json_path)


if __name__ == "__main__":
    main()

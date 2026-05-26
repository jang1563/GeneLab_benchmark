#!/usr/bin/env python3
"""Build the derived v9 human organoid DE reference table."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.organoid_de_reference import write_organoid_de_reference  # noqa: E402
from spacebio_bench.source_audit import read_source_inventory  # noqa: E402


def _read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-inventory",
        default="v9/human_organoid/source_inventory.draft.csv",
        help="Input human organoid source inventory CSV path.",
    )
    parser.add_argument(
        "--signature-reference-audit",
        default="v9/human_organoid/signature_reference_audit.draft.csv",
        help="Input human organoid signature-reference audit CSV path.",
    )
    parser.add_argument(
        "--csv",
        default="v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz",
        help="Output compact derived DE reference CSV or CSV.GZ path.",
    )
    parser.add_argument(
        "--json",
        default=(
            "v9/human_organoid/de_references/"
            "human_organoid_de_reference_manifest.draft.json"
        ),
        help="Output derived DE reference manifest JSON path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_rows = read_source_inventory(args.source_inventory)
    audit_rows = _read_csv(args.signature_reference_audit)
    csv_path, json_path = write_organoid_de_reference(
        source_rows=source_rows,
        audit_rows=audit_rows,
        csv_path=args.csv,
        json_path=args.json,
    )
    print(csv_path)
    print(json_path)


if __name__ == "__main__":
    main()

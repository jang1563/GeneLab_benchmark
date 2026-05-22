#!/usr/bin/env python3
"""Build draft source inventories for v9 organoid and multispecies extensions."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.extension_sources import (  # noqa: E402
    HUMAN_ORGANOID_DRAFT_SOURCES,
    MULTISPECIES_DRAFT_SOURCES,
    write_extension_source_inventory,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--human-organoid-csv",
        default="v9/human_organoid/source_inventory.draft.csv",
        help="Output CSV path for the human organoid draft source inventory.",
    )
    parser.add_argument(
        "--human-organoid-json",
        default="v9/human_organoid/source_inventory.draft.json",
        help="Output JSON path for the human organoid draft source inventory.",
    )
    parser.add_argument(
        "--multispecies-csv",
        default="v9/multispecies/source_inventory.draft.csv",
        help="Output CSV path for the multispecies draft source inventory.",
    )
    parser.add_argument(
        "--multispecies-json",
        default="v9/multispecies/source_inventory.draft.json",
        help="Output JSON path for the multispecies draft source inventory.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    human_csv, human_json = write_extension_source_inventory(
        HUMAN_ORGANOID_DRAFT_SOURCES,
        csv_path=args.human_organoid_csv,
        json_path=args.human_organoid_json,
    )
    multispecies_csv, multispecies_json = write_extension_source_inventory(
        MULTISPECIES_DRAFT_SOURCES,
        csv_path=args.multispecies_csv,
        json_path=args.multispecies_json,
    )
    print(human_csv)
    print(human_json)
    print(multispecies_csv)
    print(multispecies_json)


if __name__ == "__main__":
    main()


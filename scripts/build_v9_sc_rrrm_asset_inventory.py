#!/usr/bin/env python3
"""Build the V9-SC-001 RRRM single-cell asset inventory."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_SINGLE_CELL_ASSET_INVENTORY_ID,
    write_single_cell_asset_inventory,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve legacy RRRM assets.",
    )
    parser.add_argument(
        "--seed-asset-paths",
        default="v9/reports/recenter_decision/single_cell_asset_paths.json",
        help="JSON list of legacy single-cell asset paths to inventory.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/sc_spaceflight",
        help="Directory for v9 single-cell asset inventory artifacts.",
    )
    parser.add_argument(
        "--inventory-id",
        default=DEFAULT_SINGLE_CELL_ASSET_INVENTORY_ID,
        help="Identifier for the inventory artifacts.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_single_cell_asset_inventory(
        repo_root=args.repo_root,
        output_dir=args.output_dir,
        inventory_id=args.inventory_id,
        seed_asset_paths=args.seed_asset_paths,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["asset_csv"])
    print(outputs["asset_json"])
    print(outputs["payload_csv"])
    print(outputs["payload_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

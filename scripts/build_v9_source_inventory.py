#!/usr/bin/env python3
"""Build a source-level inventory for v9 task manifests."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import TaskRegistry, build_source_inventory, write_source_inventory


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest-dir",
        default="v9/task_manifests",
        help="Directory containing v9 task manifests.",
    )
    parser.add_argument(
        "--csv",
        default="v9/source_inventory.csv",
        help="Output source inventory CSV path.",
    )
    parser.add_argument(
        "--json",
        default="v9/source_inventory.json",
        help="Output source inventory JSON path.",
    )
    parser.add_argument(
        "--release-target",
        default="v9_alpha_public_bulk_candidate",
        help="Release target label to attach to each source row.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    registry = TaskRegistry.from_dir(args.manifest_dir)
    rows = build_source_inventory(registry, release_target=args.release_target)
    csv_path, json_path = write_source_inventory(
        rows,
        csv_path=args.csv,
        json_path=args.json,
    )
    print(csv_path)
    print(json_path)


if __name__ == "__main__":
    main()

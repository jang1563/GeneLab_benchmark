#!/usr/bin/env python3
"""Build a compact index over v9 SpaceBio-Bench task manifests."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.registry import TaskRegistry, write_task_index


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest-dir",
        default="v9/task_manifests",
        help="Directory containing v9 task manifest JSON files.",
    )
    parser.add_argument(
        "--csv",
        default="v9/task_manifest_index.csv",
        help="Output CSV index path.",
    )
    parser.add_argument(
        "--json",
        default="v9/task_manifest_index.json",
        help="Output JSON index path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    registry = TaskRegistry.from_dir(args.manifest_dir)
    csv_path, json_path = write_task_index(
        registry,
        csv_path=args.csv,
        json_path=args.json,
    )
    print(csv_path)
    print(json_path)


if __name__ == "__main__":
    main()

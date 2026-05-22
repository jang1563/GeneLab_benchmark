#!/usr/bin/env python3
"""Build a draft Frictionless Data Package descriptor for v9 public bulk artifacts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import build_v9_public_bulk_datapackage, write_datapackage


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=str(REPO_ROOT),
        help="Repository root containing v9 artifacts.",
    )
    parser.add_argument(
        "--json",
        default="v9/datapackage.draft.json",
        help="Output Data Package descriptor path.",
    )
    parser.add_argument(
        "--version",
        default="0.1.0-draft",
        help="Draft package version label.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    package = build_v9_public_bulk_datapackage(
        repo_root=args.repo_root,
        version=args.version,
    )
    output_path = write_datapackage(package, args.json)
    print(output_path)


if __name__ == "__main__":
    main()

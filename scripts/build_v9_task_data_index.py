#!/usr/bin/env python3
"""Build a fold-level data index for v9 bulk LOMO task manifests."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.data import bulk_task_summary_rows, load_all_bulk_tasks


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest-dir",
        default="v9/task_manifests",
        help="Directory containing v9 task manifests.",
    )
    parser.add_argument(
        "--csv",
        default="v9/task_data_index.csv",
        help="Output fold-level CSV index.",
    )
    parser.add_argument(
        "--json",
        default="v9/task_data_index.json",
        help="Output fold-level JSON index.",
    )
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve relative paths in the index.",
    )
    return parser.parse_args()


def write_rows(rows: list[dict[str, str]], csv_path: Path, json_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")


def main() -> None:
    args = parse_args()
    repo_root = Path(args.repo_root)
    tasks = load_all_bulk_tasks(manifest_dir=args.manifest_dir, repo_root=repo_root)
    rows = bulk_task_summary_rows(tasks)
    write_rows(rows, Path(args.csv), Path(args.json))
    print(args.csv)
    print(args.json)


if __name__ == "__main__":
    main()

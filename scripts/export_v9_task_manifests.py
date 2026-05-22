#!/usr/bin/env python3
"""Export existing GeneLab LOMO task metadata as v9 SpaceBio-Bench manifests."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.tasks.legacy import (
    legacy_task_info_to_manifest,
    load_legacy_task_info,
    write_manifest,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--tasks-dir",
        default="tasks",
        help="Directory containing legacy task subdirectories.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/task_manifests",
        help="Directory for generated v9 task manifests.",
    )
    parser.add_argument(
        "--task",
        action="append",
        default=[],
        help="Optional legacy task directory name to export. Can be repeated.",
    )
    return parser.parse_args()


def iter_task_info_files(tasks_dir: Path, task_names: list[str]) -> list[Path]:
    if task_names:
        return [tasks_dir / task_name / "task_info.json" for task_name in task_names]
    return sorted(tasks_dir.glob("*/task_info.json"))


def main() -> None:
    args = parse_args()
    tasks_dir = Path(args.tasks_dir)
    output_dir = Path(args.output_dir)

    written: list[Path] = []
    for task_info_path in iter_task_info_files(tasks_dir, args.task):
        if not task_info_path.exists():
            raise FileNotFoundError(task_info_path)
        task_info = load_legacy_task_info(task_info_path)
        manifest = legacy_task_info_to_manifest(
            task_info,
            task_dir=task_info_path.parent,
        )
        output_path = output_dir / f"{manifest['task_id']}.json"
        written.append(write_manifest(manifest, output_path))

    for path in written:
        print(path)


if __name__ == "__main__":
    main()

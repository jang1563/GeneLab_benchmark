#!/usr/bin/env python3
"""Preflight the OSD-120 diagnostic packaging rebuild gate."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_INTERACTION_REBUILD_GATE_ID,
    write_multispecies_interaction_rebuild_gate_manifest,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve OSD-120 report paths.",
    )
    parser.add_argument(
        "--reports-root",
        default="v9/multispecies/reports",
        help="Directory containing OSD-120 report subdirectories.",
    )
    parser.add_argument(
        "--task-manifest",
        default=(
            "v9/multispecies/interaction_task_manifests/"
            "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
        ),
        help="OSD-120 interaction task manifest.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/multispecies/reports/interaction_rebuild_gate",
        help="Directory for rebuild gate manifest outputs.",
    )
    parser.add_argument(
        "--gate-id",
        default=DEFAULT_INTERACTION_REBUILD_GATE_ID,
        help="Identifier for the emitted rebuild gate.",
    )
    parser.add_argument(
        "--mode",
        choices=["preflight"],
        default="preflight",
        help="Gate mode. Preflight hashes existing packaging outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mode = "preflight_existing_outputs"
    outputs = write_multispecies_interaction_rebuild_gate_manifest(
        repo_root=args.repo_root,
        reports_root=args.reports_root,
        task_manifest=args.task_manifest,
        output_dir=args.output_dir,
        gate_id=args.gate_id,
        mode=mode,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["step_csv"])
    print(outputs["step_json"])
    print(outputs["environment_csv"])
    print(outputs["environment_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

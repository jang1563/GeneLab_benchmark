#!/usr/bin/env python3
"""Build the OSD-120 public metadata package skeleton."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_INTERACTION_PUBLIC_METADATA_PACKAGE_ID,
    write_multispecies_interaction_public_metadata_package,
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
        default="v9/multispecies/reports/interaction_public_metadata_package",
        help="Directory for public metadata package outputs.",
    )
    parser.add_argument(
        "--package-id",
        default=DEFAULT_INTERACTION_PUBLIC_METADATA_PACKAGE_ID,
        help="Identifier for the emitted metadata package.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_public_metadata_package(
        repo_root=args.repo_root,
        reports_root=args.reports_root,
        task_manifest=args.task_manifest,
        output_dir=args.output_dir,
        package_id=args.package_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["target_csv"])
    print(outputs["target_json"])
    print(outputs["field_csv"])
    print(outputs["field_json"])
    print(outputs["reference_csv"])
    print(outputs["reference_json"])
    print(outputs["metadata_skeleton_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

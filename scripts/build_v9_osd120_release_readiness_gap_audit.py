#!/usr/bin/env python3
"""Build the OSD-120 release-readiness gap audit tables."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_INTERACTION_RELEASE_READINESS_AUDIT_ID,
    write_multispecies_interaction_release_readiness_gap_audit,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve OSD-120 input paths.",
    )
    parser.add_argument(
        "--reports-root",
        default="v9/multispecies/reports",
        help="Directory containing OSD-120 interaction report subdirectories.",
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
        default="v9/multispecies/reports/interaction_release_readiness_gap_audit",
        help="Directory for release-readiness audit outputs.",
    )
    parser.add_argument(
        "--audit-id",
        default=DEFAULT_INTERACTION_RELEASE_READINESS_AUDIT_ID,
        help="Identifier for the emitted release-readiness audit.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_release_readiness_gap_audit(
        repo_root=args.repo_root,
        reports_root=args.reports_root,
        task_manifest=args.task_manifest,
        output_dir=args.output_dir,
        audit_id=args.audit_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["gap_csv"])
    print(outputs["gap_json"])
    print(outputs["reference_csv"])
    print(outputs["reference_json"])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build the OSD-120 sparse diagnostic candidate package."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_INTERACTION_DIAGNOSTIC_CANDIDATE_ID,
    DEFAULT_INTERACTION_DIAGNOSTIC_COMPARATOR_ID,
    DEFAULT_INTERACTION_DIAGNOSTIC_PACKAGE_ID,
    write_multispecies_interaction_diagnostic_candidate_package,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve report and manifest paths.",
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
        default="v9/multispecies/reports/interaction_diagnostic_candidate_package",
        help="Directory for diagnostic candidate package outputs.",
    )
    parser.add_argument(
        "--package-id",
        default=DEFAULT_INTERACTION_DIAGNOSTIC_PACKAGE_ID,
        help="Identifier for the emitted candidate package.",
    )
    parser.add_argument(
        "--candidate-id",
        default=DEFAULT_INTERACTION_DIAGNOSTIC_CANDIDATE_ID,
        help="Candidate id from the baseline ladder summary.",
    )
    parser.add_argument(
        "--comparator-candidate-id",
        default=DEFAULT_INTERACTION_DIAGNOSTIC_COMPARATOR_ID,
        help="Comparator candidate id from the baseline ladder summary.",
    )
    parser.add_argument(
        "--stable-feature-frequency-threshold",
        type=float,
        default=0.5,
        help="Minimum subsample selection frequency for stable feature rows.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_diagnostic_candidate_package(
        repo_root=args.repo_root,
        reports_root=args.reports_root,
        task_manifest=args.task_manifest,
        output_dir=args.output_dir,
        package_id=args.package_id,
        candidate_id=args.candidate_id,
        comparator_candidate_id=args.comparator_candidate_id,
        stable_feature_frequency_threshold=args.stable_feature_frequency_threshold,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["focus_csv"])
    print(outputs["focus_json"])
    print(outputs["feature_csv"])
    print(outputs["feature_json"])
    print(outputs["claim_csv"])
    print(outputs["claim_json"])


if __name__ == "__main__":
    main()

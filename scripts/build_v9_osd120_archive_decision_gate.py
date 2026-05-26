#!/usr/bin/env python3
"""Build the OSD-120 archive identifier and license decision gate."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_INTERACTION_ARCHIVE_DECISION_GATE_ID,
    write_multispecies_interaction_archive_decision_gate,
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
        "--output-dir",
        default="v9/multispecies/reports/interaction_archive_decision_gate",
        help="Directory for archive decision gate outputs.",
    )
    parser.add_argument(
        "--decision-id",
        default=DEFAULT_INTERACTION_ARCHIVE_DECISION_GATE_ID,
        help="Identifier for the emitted archive decision gate.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_archive_decision_gate(
        repo_root=args.repo_root,
        reports_root=args.reports_root,
        output_dir=args.output_dir,
        decision_id=args.decision_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["archive_option_csv"])
    print(outputs["archive_option_json"])
    print(outputs["license_csv"])
    print(outputs["license_json"])
    print(outputs["creator_csv"])
    print(outputs["creator_json"])
    print(outputs["reference_csv"])
    print(outputs["reference_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

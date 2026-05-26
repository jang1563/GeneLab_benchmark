#!/usr/bin/env python3
"""Build the OSD-120 release-owner citation metadata fill scaffold."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_INTERACTION_CITATION_METADATA_FILL_ID,
    write_multispecies_interaction_citation_metadata_fill,
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
        default="v9/multispecies/reports/interaction_citation_metadata_fill",
        help="Directory for citation metadata fill outputs.",
    )
    parser.add_argument(
        "--owner-metadata",
        default=None,
        help=(
            "Optional CSV or JSON file with field_id, supplied_value, supplied_by, "
            "supplied_date, and supplied_evidence columns/keys."
        ),
    )
    parser.add_argument(
        "--fill-id",
        default=DEFAULT_INTERACTION_CITATION_METADATA_FILL_ID,
        help="Identifier for the emitted citation metadata fill scaffold.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_citation_metadata_fill(
        repo_root=args.repo_root,
        reports_root=args.reports_root,
        output_dir=args.output_dir,
        owner_metadata_path=args.owner_metadata,
        fill_id=args.fill_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["intake_template_csv"])
    print(outputs["intake_template_json"])
    print(outputs["fill_status_csv"])
    print(outputs["fill_status_json"])
    print(outputs["descriptor_preview_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

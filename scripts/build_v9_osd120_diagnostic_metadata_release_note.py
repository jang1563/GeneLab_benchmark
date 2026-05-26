#!/usr/bin/env python3
"""Build the OSD-120 diagnostic metadata release-note closeout."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_INTERACTION_DIAGNOSTIC_METADATA_RELEASE_NOTE_ID,
    write_multispecies_interaction_diagnostic_metadata_release_note,
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
        default="v9/multispecies/reports/interaction_diagnostic_metadata_release_note",
        help="Directory for diagnostic metadata release-note outputs.",
    )
    parser.add_argument(
        "--release-note-id",
        default=DEFAULT_INTERACTION_DIAGNOSTIC_METADATA_RELEASE_NOTE_ID,
        help="Identifier for the emitted release-note closeout.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_diagnostic_metadata_release_note(
        repo_root=args.repo_root,
        reports_root=args.reports_root,
        output_dir=args.output_dir,
        release_note_id=args.release_note_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["section_csv"])
    print(outputs["section_json"])
    print(outputs["claim_csv"])
    print(outputs["claim_json"])
    print(outputs["retry_csv"])
    print(outputs["retry_json"])
    print(outputs["release_note_md"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

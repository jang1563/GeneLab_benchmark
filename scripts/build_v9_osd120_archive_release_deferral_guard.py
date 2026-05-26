#!/usr/bin/env python3
"""Build the OSD-120 archive-release deferral/application guard."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_INTERACTION_ARCHIVE_RELEASE_DEFERRAL_GUARD_ID,
    write_multispecies_interaction_archive_release_deferral_guard,
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
        default="v9/multispecies/reports/interaction_archive_release_deferral_guard",
        help="Directory for archive-release deferral guard outputs.",
    )
    parser.add_argument(
        "--guard-id",
        default=DEFAULT_INTERACTION_ARCHIVE_RELEASE_DEFERRAL_GUARD_ID,
        help="Identifier for the emitted archive-release deferral guard.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_archive_release_deferral_guard(
        repo_root=args.repo_root,
        reports_root=args.reports_root,
        output_dir=args.output_dir,
        guard_id=args.guard_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["guard_csv"])
    print(outputs["guard_json"])
    print(outputs["action_csv"])
    print(outputs["action_json"])
    print(outputs["mutation_guard_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

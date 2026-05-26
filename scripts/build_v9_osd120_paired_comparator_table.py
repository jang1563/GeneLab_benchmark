#!/usr/bin/env python3
"""Build the OSD-120 sparse L1 c1-versus-c0.3 paired comparator table."""

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
    DEFAULT_INTERACTION_PAIRED_COMPARATOR_TABLE_ID,
    write_multispecies_interaction_paired_comparator_table,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve report paths.",
    )
    parser.add_argument(
        "--reports-root",
        default="v9/multispecies/reports",
        help="Directory containing OSD-120 interaction report subdirectories.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/multispecies/reports/interaction_paired_comparator_table",
        help="Directory for paired comparator outputs.",
    )
    parser.add_argument(
        "--paired-table-id",
        default=DEFAULT_INTERACTION_PAIRED_COMPARATOR_TABLE_ID,
        help="Identifier for the emitted paired table.",
    )
    parser.add_argument(
        "--primary-candidate-id",
        default=DEFAULT_INTERACTION_DIAGNOSTIC_CANDIDATE_ID,
        help="Primary candidate id from the ladder summary.",
    )
    parser.add_argument(
        "--compact-candidate-id",
        default=DEFAULT_INTERACTION_DIAGNOSTIC_COMPARATOR_ID,
        help="Compact comparator id from the ladder summary.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_paired_comparator_table(
        repo_root=args.repo_root,
        reports_root=args.reports_root,
        output_dir=args.output_dir,
        paired_table_id=args.paired_table_id,
        primary_candidate_id=args.primary_candidate_id,
        compact_candidate_id=args.compact_candidate_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["focus_csv"])
    print(outputs["focus_json"])
    print(outputs["decision_md"])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build the V9-BULK-ALPHA-002 public bulk metadata snapshot decision."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_PUBLIC_BULK_ALPHA_SNAPSHOT_DECISION_ID,
    write_public_bulk_alpha_snapshot_decision,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve v9 public bulk inputs.",
    )
    parser.add_argument(
        "--gap-report-dir",
        default="v9/reports/public_bulk_alpha_gap_matrix",
        help="Directory containing V9-BULK-ALPHA-001 gap-matrix artifacts.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/reports/public_bulk_alpha_snapshot_decision",
        help="Directory for public bulk alpha snapshot-decision tables.",
    )
    parser.add_argument(
        "--snapshot-decision-id",
        default=DEFAULT_PUBLIC_BULK_ALPHA_SNAPSHOT_DECISION_ID,
        help="Identifier for the snapshot-decision artifacts.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_public_bulk_alpha_snapshot_decision(
        repo_root=args.repo_root,
        gap_report_dir=args.gap_report_dir,
        output_dir=args.output_dir,
        snapshot_decision_id=args.snapshot_decision_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["option_csv"])
    print(outputs["option_json"])
    print(outputs["claim_csv"])
    print(outputs["claim_json"])
    print(outputs["language_csv"])
    print(outputs["language_json"])
    print(outputs["action_csv"])
    print(outputs["action_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

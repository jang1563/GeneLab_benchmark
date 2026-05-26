#!/usr/bin/env python3
"""Build the V9-REFOCUS-001 post-OSD-120 recenter decision."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_V9_RECENTER_DECISION_ID,
    write_v9_recenter_decision,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve v9 recenter inputs.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/reports/recenter_decision",
        help="Directory for V9-REFOCUS-001 output tables.",
    )
    parser.add_argument(
        "--recenter-id",
        default=DEFAULT_V9_RECENTER_DECISION_ID,
        help="Identifier for the recenter decision artifacts.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_v9_recenter_decision(
        repo_root=args.repo_root,
        output_dir=args.output_dir,
        recenter_id=args.recenter_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["candidate_csv"])
    print(outputs["candidate_json"])
    print(outputs["action_csv"])
    print(outputs["action_json"])
    print(outputs["asset_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

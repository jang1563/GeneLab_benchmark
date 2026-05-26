#!/usr/bin/env python3
"""Build the OSD-120 interaction baseline ladder consolidation tables."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import write_multispecies_interaction_baseline_ladder  # noqa: E402


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
        default="v9/multispecies/reports/interaction_baseline_ladder",
        help="Directory for baseline ladder outputs.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_baseline_ladder(
        repo_root=args.repo_root,
        reports_root=args.reports_root,
        output_dir=args.output_dir,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["focus_csv"])
    print(outputs["focus_json"])


if __name__ == "__main__":
    main()

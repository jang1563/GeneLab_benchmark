#!/usr/bin/env python3
"""Build human-facing OSD-120 sparse diagnostic figure/table artifacts."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import write_multispecies_interaction_figure_table_package  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve package paths.",
    )
    parser.add_argument(
        "--package-dir",
        default="v9/multispecies/reports/interaction_diagnostic_candidate_package",
        help="Directory containing diagnostic candidate package tables.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/multispecies/reports/interaction_figure_table_package",
        help="Directory for figure/table draft outputs.",
    )
    parser.add_argument(
        "--figure-table-id",
        default="osd120_sparse_l1_c1_focus_table",
        help="Identifier for the emitted figure/table rows.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_figure_table_package(
        repo_root=args.repo_root,
        package_dir=args.package_dir,
        output_dir=args.output_dir,
        figure_table_id=args.figure_table_id,
    )
    print(outputs["main_csv"])
    print(outputs["main_json"])
    print(outputs["feature_csv"])
    print(outputs["feature_json"])
    print(outputs["caption_md"])
    print(outputs["claim_box_md"])


if __name__ == "__main__":
    main()

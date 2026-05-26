#!/usr/bin/env python3
"""Build the OSD-120 RO-Crate and citation freeze scaffold."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_INTERACTION_RO_CRATE_EXPORT_SCAFFOLD_ID,
    write_multispecies_interaction_ro_crate_citation_scaffold,
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
        default="v9/multispecies/reports/interaction_ro_crate_citation_scaffold",
        help="Directory for RO-Crate/citation scaffold outputs.",
    )
    parser.add_argument(
        "--scaffold-id",
        default=DEFAULT_INTERACTION_RO_CRATE_EXPORT_SCAFFOLD_ID,
        help="Identifier for the emitted RO-Crate/citation scaffold.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_ro_crate_citation_scaffold(
        repo_root=args.repo_root,
        reports_root=args.reports_root,
        output_dir=args.output_dir,
        scaffold_id=args.scaffold_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["validation_csv"])
    print(outputs["validation_json"])
    print(outputs["citation_csv"])
    print(outputs["citation_json"])
    print(outputs["ro_crate_json"])
    print(outputs["data_package_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

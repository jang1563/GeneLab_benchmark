#!/usr/bin/env python3
"""Build the V9-BULK-ALPHA-001 public bulk alpha freeze-gap matrix."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_PUBLIC_BULK_ALPHA_GAP_MATRIX_ID,
    write_public_bulk_alpha_gap_matrix,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve v9 public bulk inputs.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/reports/public_bulk_alpha_gap_matrix",
        help="Directory for public bulk alpha gap-matrix output tables.",
    )
    parser.add_argument(
        "--gap-matrix-id",
        default=DEFAULT_PUBLIC_BULK_ALPHA_GAP_MATRIX_ID,
        help="Identifier for the gap-matrix artifacts.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_public_bulk_alpha_gap_matrix(
        repo_root=args.repo_root,
        output_dir=args.output_dir,
        gap_matrix_id=args.gap_matrix_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["gap_csv"])
    print(outputs["gap_json"])
    print(outputs["payload_csv"])
    print(outputs["payload_json"])
    print(outputs["claim_csv"])
    print(outputs["claim_json"])
    print(outputs["update_csv"])
    print(outputs["update_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

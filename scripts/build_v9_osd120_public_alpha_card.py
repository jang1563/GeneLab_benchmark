#!/usr/bin/env python3
"""Build the OSD-120 diagnostic public-alpha card draft."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_INTERACTION_PUBLIC_ALPHA_CARD_ID,
    write_multispecies_interaction_public_alpha_card,
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
        help="Directory containing OSD-120 report subdirectories.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/multispecies/reports/interaction_public_alpha_card",
        help="Directory for public-alpha card outputs.",
    )
    parser.add_argument(
        "--card-id",
        default=DEFAULT_INTERACTION_PUBLIC_ALPHA_CARD_ID,
        help="Identifier for the emitted card summary.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_public_alpha_card(
        repo_root=args.repo_root,
        reports_root=args.reports_root,
        output_dir=args.output_dir,
        card_id=args.card_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["card_md"])


if __name__ == "__main__":
    main()

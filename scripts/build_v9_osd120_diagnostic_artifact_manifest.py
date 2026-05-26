#!/usr/bin/env python3
"""Build the OSD-120 diagnostic artifact manifest and claim map."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_INTERACTION_DIAGNOSTIC_ARTIFACT_MANIFEST_ID,
    write_multispecies_interaction_diagnostic_artifact_manifest,
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
        default="v9/multispecies/reports/interaction_diagnostic_artifact_manifest",
        help="Directory for diagnostic artifact manifest outputs.",
    )
    parser.add_argument(
        "--manifest-id",
        default=DEFAULT_INTERACTION_DIAGNOSTIC_ARTIFACT_MANIFEST_ID,
        help="Identifier for the emitted artifact manifest.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_diagnostic_artifact_manifest(
        repo_root=args.repo_root,
        reports_root=args.reports_root,
        output_dir=args.output_dir,
        manifest_id=args.manifest_id,
    )
    print(outputs["artifact_csv"])
    print(outputs["artifact_json"])
    print(outputs["claim_csv"])
    print(outputs["claim_json"])


if __name__ == "__main__":
    main()

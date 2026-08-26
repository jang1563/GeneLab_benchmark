#!/usr/bin/env python3
"""Build the V9-SC-006c OSDR processed-payload discovery package."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_SC_PROCESSED_PAYLOAD_DISCOVERY_ID,
    write_sc_processed_payload_discovery,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve v9 single-cell artifacts.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/sc_spaceflight",
        help="Directory for v9 single-cell OSDR discovery artifacts.",
    )
    parser.add_argument(
        "--discovery-id",
        default=DEFAULT_SC_PROCESSED_PAYLOAD_DISCOVERY_ID,
        help="Identifier for the OSDR processed-payload discovery artifacts.",
    )
    parser.add_argument(
        "--manifest-path",
        default=(
            "v9/sc_spaceflight/task_manifests/"
            "draft_rrrm1_blood_single_cell_spaceflight.json"
        ),
        help="Draft sc_spaceflight manifest used to anchor the discovery package.",
    )
    parser.add_argument(
        "--condition-map-path",
        default="v2/docs/RRRM1_SRX_CONDITION_MAP.csv",
        help="RRRM-1 SRX condition map used to enumerate expected blood files.",
    )
    parser.add_argument(
        "--canonical-payload-path",
        default=(
            "v9/sc_spaceflight/payloads/rrrm1_blood/"
            "OSD-918_blood_rrrm1_bench.h5ad"
        ),
        help="Canonical h5ad path to stage and audit if a payload is later found.",
    )
    parser.add_argument(
        "--api-base",
        default="https://visualization.osdr.nasa.gov/biodata/api/v2",
        help="OSDR Biodata API base URL.",
    )
    parser.add_argument(
        "--api-listing-json",
        default=None,
        help="Optional local OSDR file-list JSON used instead of a live API fetch.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_sc_processed_payload_discovery(
        repo_root=args.repo_root,
        output_dir=args.output_dir,
        discovery_id=args.discovery_id,
        manifest_path=args.manifest_path,
        condition_map_path=args.condition_map_path,
        canonical_payload_path=args.canonical_payload_path,
        api_base=args.api_base,
        api_listing_json=args.api_listing_json,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["files_csv"])
    print(outputs["files_json"])
    print(outputs["coverage_csv"])
    print(outputs["coverage_json"])
    print(outputs["owner_request_csv"])
    print(outputs["owner_request_json"])
    print(outputs["deferral_csv"])
    print(outputs["deferral_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build the V9-SC-006b RRRM-1 blood external payload availability package."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_SC_EXTERNAL_CHECKED_BASES,
    DEFAULT_SC_EXTERNAL_PAYLOAD_AVAILABILITY_ID,
    write_sc_external_payload_availability,
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
        help="Directory for v9 single-cell external payload artifacts.",
    )
    parser.add_argument(
        "--availability-id",
        default=DEFAULT_SC_EXTERNAL_PAYLOAD_AVAILABILITY_ID,
        help="Identifier for the external payload availability artifacts.",
    )
    parser.add_argument(
        "--manifest-path",
        default=(
            "v9/sc_spaceflight/task_manifests/"
            "draft_rrrm1_blood_single_cell_spaceflight.json"
        ),
        help="Draft sc_spaceflight manifest used to anchor the availability package.",
    )
    parser.add_argument(
        "--condition-map-path",
        default="v2/docs/RRRM1_SRX_CONDITION_MAP.csv",
        help="RRRM-1 SRX condition map used to enumerate expected blood matrices.",
    )
    parser.add_argument(
        "--canonical-payload-path",
        default=(
            "v9/sc_spaceflight/payloads/rrrm1_blood/"
            "OSD-918_blood_rrrm1_bench.h5ad"
        ),
        help="Canonical h5ad path to stage and audit.",
    )
    parser.add_argument(
        "--checked-base",
        action="append",
        default=None,
        help=(
            "External base path checked for legacy h5ad and STARsolo matrices. "
            "May be repeated."
        ),
    )
    parser.add_argument(
        "--external-probe-note",
        default=(
            "Manual probe 2026-06-03: exact OSD-918 h5ad and 8 SRX STARsolo "
            "paths were absent in the checked local/remote scopes."
        ),
        help="Concise note from manual/remote payload probes.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_sc_external_payload_availability(
        repo_root=args.repo_root,
        output_dir=args.output_dir,
        availability_id=args.availability_id,
        manifest_path=args.manifest_path,
        condition_map_path=args.condition_map_path,
        canonical_payload_path=args.canonical_payload_path,
        checked_bases=args.checked_base or DEFAULT_SC_EXTERNAL_CHECKED_BASES,
        external_probe_note=args.external_probe_note,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["candidates_csv"])
    print(outputs["candidates_json"])
    print(outputs["matrix_csv"])
    print(outputs["matrix_json"])
    print(outputs["copy_decision_csv"])
    print(outputs["copy_decision_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

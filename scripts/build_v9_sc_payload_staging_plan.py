#!/usr/bin/env python3
"""Build the V9-SC-004 AnnData payload staging and obs/var audit plan."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_SC_PAYLOAD_STAGING_PLAN_ID,
    write_sc_payload_staging_plan,
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
        help="Directory for v9 single-cell payload-staging artifacts.",
    )
    parser.add_argument(
        "--payload-plan-id",
        default=DEFAULT_SC_PAYLOAD_STAGING_PLAN_ID,
        help="Identifier for the payload-staging plan artifacts.",
    )
    parser.add_argument(
        "--manifest-path",
        default=(
            "v9/sc_spaceflight/task_manifests/"
            "draft_rrrm1_blood_single_cell_spaceflight.json"
        ),
        help="Draft sc_spaceflight manifest used to anchor the staging plan.",
    )
    parser.add_argument(
        "--metric-spec-summary-path",
        default="v9/sc_spaceflight/metric_spec_summary.csv",
        help="Metric-spec summary used as evidence that V9-SC-003 is complete.",
    )
    parser.add_argument(
        "--canonical-payload-path",
        default=(
            "v9/sc_spaceflight/payloads/rrrm1_blood/"
            "OSD-918_blood_rrrm1_bench.h5ad"
        ),
        help="Planned canonical h5ad path for the runnable v9 payload.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_sc_payload_staging_plan(
        repo_root=args.repo_root,
        output_dir=args.output_dir,
        payload_plan_id=args.payload_plan_id,
        manifest_path=args.manifest_path,
        metric_spec_summary_path=args.metric_spec_summary_path,
        canonical_payload_path=args.canonical_payload_path,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["candidates_csv"])
    print(outputs["candidates_json"])
    print(outputs["audit_csv"])
    print(outputs["audit_json"])
    print(outputs["actions_csv"])
    print(outputs["actions_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

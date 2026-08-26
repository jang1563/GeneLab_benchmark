#!/usr/bin/env python3
"""Build the V9-SC-006 RRRM-1 blood payload staging execution package."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_SC_PAYLOAD_STAGING_EXECUTION_ID,
    write_sc_payload_staging_execution,
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
        help="Directory for v9 single-cell payload staging execution artifacts.",
    )
    parser.add_argument(
        "--execution-id",
        default=DEFAULT_SC_PAYLOAD_STAGING_EXECUTION_ID,
        help="Identifier for the payload staging execution artifacts.",
    )
    parser.add_argument(
        "--manifest-path",
        default=(
            "v9/sc_spaceflight/task_manifests/"
            "draft_rrrm1_blood_single_cell_spaceflight.json"
        ),
        help="Draft sc_spaceflight manifest used to anchor the execution package.",
    )
    parser.add_argument(
        "--candidate-table-path",
        default="v9/sc_spaceflight/payload_staging_candidates.csv",
        help="V9-SC-004 candidate table used as the source route list.",
    )
    parser.add_argument(
        "--payload-manifest-path",
        default="v9/sc_spaceflight/payload_manifest.csv",
        help="Current V9-SC-005 payload manifest.",
    )
    parser.add_argument(
        "--obs-var-summary-path",
        default="v9/sc_spaceflight/obs_var_audit_summary.csv",
        help="Current V9-SC-005 obs/var audit summary.",
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
        "--external-probe-note",
        default="",
        help="Optional concise note from manual/remote payload probes.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_sc_payload_staging_execution(
        repo_root=args.repo_root,
        output_dir=args.output_dir,
        execution_id=args.execution_id,
        manifest_path=args.manifest_path,
        candidate_table_path=args.candidate_table_path,
        payload_manifest_path=args.payload_manifest_path,
        obs_var_summary_path=args.obs_var_summary_path,
        canonical_payload_path=args.canonical_payload_path,
        external_probe_note=args.external_probe_note,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["candidates_csv"])
    print(outputs["candidates_json"])
    print(outputs["regeneration_csv"])
    print(outputs["regeneration_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

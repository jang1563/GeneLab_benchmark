#!/usr/bin/env python3
"""Audit the V9 single-cell AnnData obs/var contract."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import DEFAULT_SC_OBS_VAR_AUDIT_ID, write_sc_obs_var_audit  # noqa: E402


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
        help="Directory for v9 single-cell obs/var audit artifacts.",
    )
    parser.add_argument(
        "--audit-id",
        default=DEFAULT_SC_OBS_VAR_AUDIT_ID,
        help="Identifier for the obs/var audit artifacts.",
    )
    parser.add_argument(
        "--manifest-path",
        default=(
            "v9/sc_spaceflight/task_manifests/"
            "draft_rrrm1_blood_single_cell_spaceflight.json"
        ),
        help="Draft sc_spaceflight manifest used to anchor the audit.",
    )
    parser.add_argument(
        "--requirements-path",
        default="v9/sc_spaceflight/obs_var_audit_requirements.csv",
        help="V9-SC-004 obs/var audit requirements table.",
    )
    parser.add_argument(
        "--canonical-payload-path",
        default=(
            "v9/sc_spaceflight/payloads/rrrm1_blood/"
            "OSD-918_blood_rrrm1_bench.h5ad"
        ),
        help="Canonical h5ad path to audit.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_sc_obs_var_audit(
        repo_root=args.repo_root,
        output_dir=args.output_dir,
        audit_id=args.audit_id,
        manifest_path=args.manifest_path,
        requirements_path=args.requirements_path,
        canonical_payload_path=args.canonical_payload_path,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["results_csv"])
    print(outputs["results_json"])
    print(outputs["payload_manifest_csv"])
    print(outputs["payload_manifest_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

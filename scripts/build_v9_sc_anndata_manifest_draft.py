#!/usr/bin/env python3
"""Build the V9-SC-002 RRRM-1 blood AnnData task-manifest draft."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_SC_ANNDATA_MANIFEST_DRAFT_ID,
    DEFAULT_SC_ANNDATA_TASK_ID,
    write_sc_anndata_manifest_draft,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve RRRM-1 evidence files.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/sc_spaceflight",
        help="Directory for v9 single-cell manifest draft artifacts.",
    )
    parser.add_argument(
        "--manifest-draft-id",
        default=DEFAULT_SC_ANNDATA_MANIFEST_DRAFT_ID,
        help="Identifier for the manifest-draft artifacts.",
    )
    parser.add_argument(
        "--task-id",
        default=DEFAULT_SC_ANNDATA_TASK_ID,
        help="Draft sc_spaceflight task id.",
    )
    parser.add_argument(
        "--selected-source-id",
        default="OSD-918",
        help="RRRM-1 OSDR source id selected for the draft task.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_sc_anndata_manifest_draft(
        repo_root=args.repo_root,
        output_dir=args.output_dir,
        manifest_draft_id=args.manifest_draft_id,
        task_id=args.task_id,
        selected_source_id=args.selected_source_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["manifest_json"])
    print(outputs["blockers_csv"])
    print(outputs["blockers_json"])
    print(outputs["review_md"])


if __name__ == "__main__":
    main()

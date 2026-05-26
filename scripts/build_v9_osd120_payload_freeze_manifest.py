#!/usr/bin/env python3
"""Build the OSD-120 diagnostic payload freeze manifest."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    DEFAULT_INTERACTION_PAYLOAD_FREEZE_ID,
    write_multispecies_interaction_payload_freeze_manifest,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve OSD-120 payload paths.",
    )
    parser.add_argument(
        "--task-manifest",
        default=(
            "v9/multispecies/interaction_task_manifests/"
            "draft_osd120_arabidopsis_root_light_interaction_spaceflight.json"
        ),
        help="OSD-120 interaction task manifest.",
    )
    parser.add_argument(
        "--source-checksum-audit",
        default="v9/multispecies/source_checksum_audit.draft.csv",
        help="Source checksum audit CSV containing the OSD-120 checksum URL.",
    )
    parser.add_argument(
        "--expression-matrix-audit",
        default="v9/multispecies/expression_matrix_audit.draft.csv",
        help="Expression matrix audit CSV containing required local payload paths.",
    )
    parser.add_argument(
        "--checksum-manifest-text-path",
        default=None,
        help="Optional local checksum manifest text fixture; skips live OSDR fetch.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/multispecies/reports/interaction_payload_freeze_manifest",
        help="Directory for payload freeze outputs.",
    )
    parser.add_argument(
        "--freeze-id",
        default=DEFAULT_INTERACTION_PAYLOAD_FREEZE_ID,
        help="Identifier for the emitted payload freeze manifest.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outputs = write_multispecies_interaction_payload_freeze_manifest(
        repo_root=args.repo_root,
        task_manifest=args.task_manifest,
        source_checksum_audit=args.source_checksum_audit,
        expression_matrix_audit=args.expression_matrix_audit,
        checksum_manifest_text_path=args.checksum_manifest_text_path,
        output_dir=args.output_dir,
        freeze_id=args.freeze_id,
    )
    print(outputs["summary_csv"])
    print(outputs["summary_json"])
    print(outputs["manifest_csv"])
    print(outputs["manifest_json"])


if __name__ == "__main__":
    main()

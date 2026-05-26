#!/usr/bin/env python3
"""Build the draft v9 OSD-120 interaction task manifest and index."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.extension_sources import read_extension_source_inventory  # noqa: E402
from spacebio_bench.extension_tasks import (  # noqa: E402
    MULTISPECIES_OSD120_INTERACTION_TASK_ID,
    build_osd120_interaction_task_manifest,
    write_task_manifest,
)
from spacebio_bench.registry import TaskRegistry, write_task_index  # noqa: E402


def _read_csv_if_exists(path: str | Path) -> list[dict[str, str]]:
    csv_path = Path(path)
    if not csv_path.exists():
        return []
    with csv_path.open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-inventory",
        default="v9/multispecies/source_inventory.draft.json",
        help="Draft multispecies source inventory JSON or CSV.",
    )
    parser.add_argument(
        "--sample-factor-table",
        default="v9/multispecies/sample_factors.draft.csv",
        help="Parsed multispecies sample-factor CSV path.",
    )
    parser.add_argument(
        "--expression-matrix-audit",
        default="v9/multispecies/expression_matrix_audit.draft.csv",
        help="Expression-matrix audit CSV path.",
    )
    parser.add_argument(
        "--source-checksum-audit",
        default="v9/multispecies/source_checksum_audit.draft.csv",
        help="Source checksum audit CSV path.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/multispecies/interaction_task_manifests",
        help="Output directory for the draft OSD-120 interaction task manifest.",
    )
    parser.add_argument(
        "--index-csv",
        default="v9/multispecies/interaction_task_manifest_index.draft.csv",
        help="Output CSV path for the draft OSD-120 interaction task index.",
    )
    parser.add_argument(
        "--index-json",
        default="v9/multispecies/interaction_task_manifest_index.draft.json",
        help="Output JSON path for the draft OSD-120 interaction task index.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_rows = read_extension_source_inventory(args.source_inventory)
    manifest = build_osd120_interaction_task_manifest(
        source_rows,
        sample_factor_rows=_read_csv_if_exists(args.sample_factor_table),
        expression_matrix_audit_rows=_read_csv_if_exists(args.expression_matrix_audit),
        source_checksum_audit_rows=_read_csv_if_exists(args.source_checksum_audit),
    )
    output_dir = Path(args.output_dir)
    manifest_path = write_task_manifest(
        manifest,
        output_dir / f"{MULTISPECIES_OSD120_INTERACTION_TASK_ID}.json",
    )
    registry = TaskRegistry.from_dir(output_dir)
    index_csv, index_json = write_task_index(
        registry,
        csv_path=args.index_csv,
        json_path=args.index_json,
    )
    print(manifest_path)
    print(index_csv)
    print(index_json)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build the draft v9 human organoid task manifest and extension index."""

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
    HUMAN_ORGANOID_TASK_ID,
    build_human_organoid_task_manifest,
    write_task_manifest,
)
from spacebio_bench.registry import TaskRegistry, write_task_index  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-inventory",
        default="v9/human_organoid/source_inventory.draft.json",
        help="Draft human organoid source inventory JSON or CSV.",
    )
    parser.add_argument(
        "--output-dir",
        default="v9/human_organoid/task_manifests",
        help="Output directory for draft human organoid task manifests.",
    )
    parser.add_argument(
        "--sample-factor-table",
        default="v9/human_organoid/sample_factors.draft.csv",
        help="Optional parsed sample-factor CSV path to include sample-backed folds.",
    )
    parser.add_argument(
        "--expression-matrix-audit",
        default="v9/human_organoid/expression_matrix_audit.draft.csv",
        help="Optional expression-matrix audit CSV path to include matrix alignment status.",
    )
    parser.add_argument(
        "--signature-reference-audit",
        default="v9/human_organoid/signature_reference_audit.draft.csv",
        help="Optional signature-reference audit CSV path to include DE/signature policy.",
    )
    parser.add_argument(
        "--index-csv",
        default="v9/human_organoid/task_manifest_index.draft.csv",
        help="Output CSV path for the draft human organoid task index.",
    )
    parser.add_argument(
        "--index-json",
        default="v9/human_organoid/task_manifest_index.draft.json",
        help="Output JSON path for the draft human organoid task index.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    rows = read_extension_source_inventory(args.source_inventory)
    sample_factor_rows = []
    sample_factor_path = Path(args.sample_factor_table)
    if sample_factor_path.exists():
        with sample_factor_path.open(newline="") as handle:
            sample_factor_rows = [dict(row) for row in csv.DictReader(handle)]
    expression_matrix_audit_rows = []
    expression_matrix_audit_path = Path(args.expression_matrix_audit)
    if expression_matrix_audit_path.exists():
        with expression_matrix_audit_path.open(newline="") as handle:
            expression_matrix_audit_rows = [dict(row) for row in csv.DictReader(handle)]
    signature_reference_audit_rows = []
    signature_reference_audit_path = Path(args.signature_reference_audit)
    if signature_reference_audit_path.exists():
        with signature_reference_audit_path.open(newline="") as handle:
            signature_reference_audit_rows = [dict(row) for row in csv.DictReader(handle)]
    manifest = build_human_organoid_task_manifest(
        rows,
        sample_factor_rows=sample_factor_rows,
        expression_matrix_audit_rows=expression_matrix_audit_rows,
        signature_reference_audit_rows=signature_reference_audit_rows,
    )
    output_dir = Path(args.output_dir)
    manifest_path = write_task_manifest(
        manifest,
        output_dir / f"{HUMAN_ORGANOID_TASK_ID}.json",
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

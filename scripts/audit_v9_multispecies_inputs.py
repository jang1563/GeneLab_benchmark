#!/usr/bin/env python3
"""Build local sample-factor and expression-matrix audits for v9 multispecies inputs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.extension_sources import read_extension_source_inventory  # noqa: E402
from spacebio_bench.multispecies_audit import (  # noqa: E402
    audit_multispecies_expression_matrices,
    build_multispecies_sample_factor_rows,
    write_multispecies_expression_matrix_audit,
    write_multispecies_sample_factor_table,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source-inventory",
        default="v9/multispecies/source_inventory.draft.csv",
        help="Input multispecies source inventory CSV or JSON path.",
    )
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve local multispecies data files.",
    )
    parser.add_argument(
        "--sample-factors-csv",
        default="v9/multispecies/sample_factors.draft.csv",
        help="Output parsed sample-factor CSV path.",
    )
    parser.add_argument(
        "--sample-factors-json",
        default="v9/multispecies/sample_factors.draft.json",
        help="Output parsed sample-factor JSON path.",
    )
    parser.add_argument(
        "--matrix-audit-csv",
        default="v9/multispecies/expression_matrix_audit.draft.csv",
        help="Output expression-matrix audit CSV path.",
    )
    parser.add_argument(
        "--matrix-audit-json",
        default="v9/multispecies/expression_matrix_audit.draft.json",
        help="Output expression-matrix audit JSON path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source_rows = read_extension_source_inventory(args.source_inventory)
    sample_factor_rows = build_multispecies_sample_factor_rows(
        source_rows,
        repo_root=args.repo_root,
    )
    sample_csv, sample_json = write_multispecies_sample_factor_table(
        sample_factor_rows,
        csv_path=args.sample_factors_csv,
        json_path=args.sample_factors_json,
    )
    audit_rows = audit_multispecies_expression_matrices(
        source_rows,
        sample_factor_rows=sample_factor_rows,
        repo_root=args.repo_root,
    )
    audit_csv, audit_json = write_multispecies_expression_matrix_audit(
        audit_rows,
        csv_path=args.matrix_audit_csv,
        json_path=args.matrix_audit_json,
    )
    print(sample_csv)
    print(sample_json)
    print(audit_csv)
    print(audit_json)


if __name__ == "__main__":
    main()

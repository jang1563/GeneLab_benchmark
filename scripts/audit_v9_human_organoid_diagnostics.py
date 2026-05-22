#!/usr/bin/env python3
"""Build sample-scale and donor-metadata diagnostics for the v9 organoid pilot."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench import (  # noqa: E402
    build_organoid_donor_metadata_audit,
    build_organoid_group_diagnostics,
    build_organoid_sample_diagnostics,
    load_human_organoid_task,
    write_organoid_donor_metadata_audit,
    write_organoid_group_diagnostics,
    write_organoid_sample_diagnostics,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest-path",
        default="v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json",
        help="Draft human organoid task manifest path.",
    )
    parser.add_argument(
        "--repo-root",
        default=".",
        help="Repository root used to resolve manifest-relative data paths.",
    )
    parser.add_argument(
        "--sample-table-audit",
        default="v9/human_organoid/sample_table_audit.draft.csv",
        help="Input sample-table audit CSV path for donor metadata availability checks.",
    )
    parser.add_argument(
        "--sample-csv",
        default="v9/human_organoid/sample_scale_diagnostics.draft.csv",
        help="Output per-sample expression scale diagnostics CSV.",
    )
    parser.add_argument(
        "--sample-json",
        default="v9/human_organoid/sample_scale_diagnostics.draft.json",
        help="Output per-sample expression scale diagnostics JSON.",
    )
    parser.add_argument(
        "--group-csv",
        default="v9/human_organoid/group_scale_diagnostics.draft.csv",
        help="Output grouped expression scale diagnostics CSV.",
    )
    parser.add_argument(
        "--group-json",
        default="v9/human_organoid/group_scale_diagnostics.draft.json",
        help="Output grouped expression scale diagnostics JSON.",
    )
    parser.add_argument(
        "--donor-csv",
        default="v9/human_organoid/donor_metadata_audit.draft.csv",
        help="Output donor/iPSC-line metadata availability audit CSV.",
    )
    parser.add_argument(
        "--donor-json",
        default="v9/human_organoid/donor_metadata_audit.draft.json",
        help="Output donor/iPSC-line metadata availability audit JSON.",
    )
    return parser.parse_args()


def _read_csv_rows(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def main() -> None:
    args = parse_args()
    task = load_human_organoid_task(
        manifest_path=args.manifest_path,
        repo_root=args.repo_root,
    )
    sample_rows = build_organoid_sample_diagnostics(task)
    group_rows = build_organoid_group_diagnostics(sample_rows)
    donor_rows = build_organoid_donor_metadata_audit(
        _read_csv_rows(Path(args.repo_root) / args.sample_table_audit),
        sample_factor_rows=task.sample_factors,
    )
    sample_csv, sample_json = write_organoid_sample_diagnostics(
        sample_rows,
        csv_path=args.sample_csv,
        json_path=args.sample_json,
    )
    group_csv, group_json = write_organoid_group_diagnostics(
        group_rows,
        csv_path=args.group_csv,
        json_path=args.group_json,
    )
    donor_csv, donor_json = write_organoid_donor_metadata_audit(
        donor_rows,
        csv_path=args.donor_csv,
        json_path=args.donor_json,
    )
    for path in (sample_csv, sample_json, group_csv, group_json, donor_csv, donor_json):
        print(path)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build GEO-derived donor/iPSC-line metadata for the v9 organoid pilot."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from spacebio_bench.organoid_geo import (  # noqa: E402
    merge_geo_metadata_with_sample_factors,
    read_csv_rows,
    read_geo_series_matrix,
    write_geo_sample_metadata,
    write_merged_sample_factors,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--series-matrix-gz",
        required=True,
        help="Input GEO GSE259421 series matrix text or .gz path.",
    )
    parser.add_argument(
        "--sample-factors",
        default="v9/human_organoid/sample_factors.draft.csv",
        help="Input parsed OSDR sample-factor CSV path.",
    )
    parser.add_argument(
        "--geo-csv",
        default="v9/human_organoid/geo_sample_metadata.draft.csv",
        help="Output parsed GEO sample metadata CSV path.",
    )
    parser.add_argument(
        "--geo-json",
        default="v9/human_organoid/geo_sample_metadata.draft.json",
        help="Output parsed GEO sample metadata JSON path.",
    )
    parser.add_argument(
        "--merged-csv",
        default="v9/human_organoid/sample_factors.draft.csv",
        help="Output merged sample-factor CSV path.",
    )
    parser.add_argument(
        "--merged-json",
        default="v9/human_organoid/sample_factors.draft.json",
        help="Output merged sample-factor JSON path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    geo_rows = read_geo_series_matrix(args.series_matrix_gz)
    sample_factor_rows = read_csv_rows(args.sample_factors)
    merged_rows = merge_geo_metadata_with_sample_factors(sample_factor_rows, geo_rows)
    geo_csv, geo_json = write_geo_sample_metadata(
        geo_rows,
        csv_path=args.geo_csv,
        json_path=args.geo_json,
    )
    merged_csv, merged_json = write_merged_sample_factors(
        merged_rows,
        csv_path=args.merged_csv,
        json_path=args.merged_json,
    )
    for path in (geo_csv, geo_json, merged_csv, merged_json):
        print(path)


if __name__ == "__main__":
    main()

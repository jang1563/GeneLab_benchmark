#!/usr/bin/env python3
"""
rrrm1_merge_h5ad.py
Build a simple RRRM-1 multi-sample AnnData object and per-sample QC summary.

This is intentionally conservative:
- loads the four per-sample h5ad files produced from STARsolo
- writes a per-sample QC summary CSV
- optionally writes a concatenated h5ad for downstream clustering/annotation
"""

import argparse
import os
from pathlib import Path

import anndata as ad
import pandas as pd

SCRATCH = Path(os.environ.get("SCRATCH_DIR", ".")) / "rrrm1_scrna"
SAMPLES = {
    918: "blood",
    920: "eye",
    924: "muscle",
    934: "skin",
}


def h5ad_path(osd_n: int, tissue: str) -> Path:
    return SCRATCH / f"OSD-{osd_n}" / f"OSD-{osd_n}_{tissue}_rrrm1.h5ad"


def summarize_sample(osd_n: int, tissue: str) -> dict:
    path = h5ad_path(osd_n, tissue)
    if not path.exists():
        raise FileNotFoundError(f"Missing h5ad: {path}")

    adata = ad.read_h5ad(path, backed="r")
    obs = adata.obs
    summary = {
        "osd": f"OSD-{osd_n}",
        "tissue": tissue,
        "h5ad_path": str(path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
    }

    for key in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
        if key in obs:
            values = obs[key]
            summary[f"median_{key}"] = float(values.median())
            summary[f"mean_{key}"] = float(values.mean())

    adata.file.close()
    return summary


def load_sample(osd_n: int, tissue: str):
    path = h5ad_path(osd_n, tissue)
    if not path.exists():
        raise FileNotFoundError(f"Missing h5ad: {path}")

    adata = ad.read_h5ad(path)
    adata.obs["sample_id"] = f"OSD-{osd_n}_{tissue}"
    adata.obs["osd"] = f"OSD-{osd_n}"
    adata.obs["tissue"] = tissue
    return adata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--summary-csv",
        default=str(SCRATCH / "RRRM1_sample_qc_summary.csv"),
        help="Output CSV path for per-sample QC summary",
    )
    parser.add_argument(
        "--merged-h5ad",
        default=str(SCRATCH / "RRRM1_merged_samples.h5ad"),
        help="Output path for merged AnnData",
    )
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Write the QC summary CSV only",
    )
    args = parser.parse_args()

    summaries = []
    for osd_n, tissue in SAMPLES.items():
        summaries.append(summarize_sample(osd_n, tissue))

    summary_df = pd.DataFrame(summaries).sort_values(["tissue", "osd"])
    summary_path = Path(args.summary_csv)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(summary_path, index=False)
    print(f"Saved sample QC summary: {summary_path}")
    print(summary_df.to_string(index=False))

    if args.summary_only:
        return

    adatas = []
    keys = []
    for osd_n, tissue in SAMPLES.items():
        adatas.append(load_sample(osd_n, tissue))
        keys.append(f"OSD-{osd_n}")

    merged = ad.concat(
        adatas,
        join="outer",
        merge="same",
        label="batch",
        keys=keys,
        index_unique="-",
    )
    merged.obs["study"] = "RRRM-1"
    merged.obs["mission"] = "RR-8"

    merged_path = Path(args.merged_h5ad)
    merged_path.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(merged_path)
    print(f"Saved merged h5ad: {merged_path}")
    print(f"Merged shape: {merged.n_obs} cells x {merged.n_vars} genes")


if __name__ == "__main__":
    main()

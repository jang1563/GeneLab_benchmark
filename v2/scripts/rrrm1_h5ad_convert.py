#!/usr/bin/env python3
"""
rrrm1_h5ad_convert.py
Convert STARsolo filtered output to h5ad for RRRM-1 scRNA-seq datasets
Run after rrrm1_starsolo_job.sh completes.

Usage: python3 rrrm1_h5ad_convert.py [--osd OSD_NUMBER]
       (if --osd not specified, processes all completed OSD dirs)
"""

import argparse
import sys
from pathlib import Path

import scanpy as sc
import pandas as pd
import numpy as np

SCRATCH = Path("/athena/masonlab/scratch/users/jak4013/rrrm1_scrna")
OSD_MAP = {
    918: "blood",
    920: "eye",
    924: "muscle",
    934: "skin",
}

def load_starsolo(osd_n: int, tissue: str) -> sc.AnnData:
    """Load STARsolo filtered output into AnnData."""
    solo_dir = SCRATCH / f"OSD-{osd_n}" / "starsolo" / "Solo.out"

    # Prefer GeneFull (counts reads spanning introns+exons) over Gene
    for gene_mode in ["GeneFull", "Gene"]:
        filt_dir = solo_dir / gene_mode / "filtered"
        if filt_dir.exists():
            print(f"  Loading {gene_mode}/filtered from {filt_dir}")
            adata = sc.read_10x_mtx(
                str(filt_dir),
                var_names="gene_symbols",
                cache=False,
            )
            break
    else:
        raise FileNotFoundError(f"No filtered STARsolo output found in {solo_dir}")

    # Add metadata
    adata.obs["tissue"]    = tissue
    adata.obs["osd"]       = f"OSD-{osd_n}"
    adata.obs["mission"]   = "RR-8"
    adata.obs["study"]     = "RRRM-1"

    # Basic QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    print(f"  {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"  Median UMI: {adata.obs['total_counts'].median():.0f}")
    print(f"  Median genes: {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"  MT%: {adata.obs['pct_counts_mt'].median():.1f}%")

    return adata


def basic_filter(adata: sc.AnnData) -> sc.AnnData:
    """Basic QC filtering: remove dead cells and empty droplets."""
    n0 = adata.n_obs
    # Remove cells with >25% mitochondrial reads (likely dead cells)
    adata = adata[adata.obs["pct_counts_mt"] < 25].copy()
    # Remove cells with too few genes (likely empty drops)
    sc.pp.filter_cells(adata, min_genes=200)
    # Remove genes expressed in <3 cells
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"  After QC: {adata.n_obs}/{n0} cells ({adata.n_vars} genes)")
    return adata


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--osd", type=int, default=None,
                        help="Single OSD number to process (default: all)")
    args = parser.parse_args()

    osd_list = [args.osd] if args.osd else list(OSD_MAP.keys())

    for osd_n in osd_list:
        tissue = OSD_MAP.get(osd_n)
        if tissue is None:
            print(f"WARNING: OSD-{osd_n} not in expected list, skipping")
            continue

        out_path = SCRATCH / f"OSD-{osd_n}" / f"OSD-{osd_n}_{tissue}_rrrm1.h5ad"

        if out_path.exists():
            print(f"SKIP (exists): {out_path}")
            continue

        print(f"\n=== OSD-{osd_n} ({tissue}) ===")
        try:
            adata = load_starsolo(osd_n, tissue)
            adata = basic_filter(adata)
            adata.write_h5ad(out_path)
            print(f"  Saved: {out_path}")
        except Exception as e:
            print(f"  ERROR: {e}")
            continue

    print("\nDone.")


if __name__ == "__main__":
    main()

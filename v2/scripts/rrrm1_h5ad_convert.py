#!/usr/bin/env python3
"""
rrrm1_h5ad_convert.py
Convert STARsolo filtered output to h5ad for RRRM-1 scRNA-seq datasets
Run after rrrm1_starsolo_job.sh completes.

Usage: python3 rrrm1_h5ad_convert.py [--osd OSD_NUMBER] [--cleanup]
       --cleanup: delete raw FASTQs and STARsolo BAM after successful h5ad write
       (if --osd not specified, processes all completed OSD dirs)
"""

import argparse
import shutil
import sys
from pathlib import Path

import scanpy as sc
from scipy import io as spio
from scipy import sparse
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
            matrix_path = filt_dir / "matrix.mtx"
            barcodes_path = filt_dir / "barcodes.tsv"
            features_path = filt_dir / "features.tsv"
            if not matrix_path.exists():
                raise FileNotFoundError(f"Missing matrix file: {matrix_path}")
            if not barcodes_path.exists():
                raise FileNotFoundError(f"Missing barcode file: {barcodes_path}")
            if not features_path.exists():
                raise FileNotFoundError(f"Missing feature file: {features_path}")

            matrix = spio.mmread(matrix_path)
            if not sparse.issparse(matrix):
                matrix = sparse.csr_matrix(matrix)
            else:
                matrix = matrix.tocsr()

            barcodes = pd.read_csv(barcodes_path, sep="\t", header=None)
            features = pd.read_csv(features_path, sep="\t", header=None)
            if features.shape[1] < 2:
                raise ValueError(f"Unexpected features.tsv format in {features_path}")

            adata = sc.AnnData(X=matrix.T.tocsr())
            adata.obs_names = barcodes.iloc[:, 0].astype(str).to_numpy()
            adata.var_names = features.iloc[:, 1].astype(str).to_numpy()
            adata.var["gene_ids"] = features.iloc[:, 0].astype(str).to_numpy()
            if features.shape[1] > 2:
                adata.var["feature_type"] = features.iloc[:, 2].astype(str).to_numpy()
            adata.var_names_make_unique()
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


def cleanup_raw(osd_n: int) -> None:
    """Delete raw FASTQs and STARsolo BAM to free space after h5ad is saved."""
    osd_dir = SCRATCH / f"OSD-{osd_n}"

    # Delete FASTQ directory
    fastq_dir = osd_dir / "fastq"
    if fastq_dir.exists():
        size_gb = sum(f.stat().st_size for f in fastq_dir.glob("*.fastq.gz")) / 1e9
        shutil.rmtree(fastq_dir)
        print(f"  Deleted FASTQs: {fastq_dir} ({size_gb:.1f} GB freed)")

    # Delete BAM files (keep Solo.out/)
    starsolo_dir = osd_dir / "starsolo"
    if starsolo_dir.exists():
        for bam in starsolo_dir.glob("*.bam"):
            size_gb = bam.stat().st_size / 1e9
            bam.unlink()
            print(f"  Deleted BAM: {bam.name} ({size_gb:.1f} GB freed)")
        for bai in starsolo_dir.glob("*.bam.bai"):
            bai.unlink()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--osd", type=int, default=None,
                        help="Single OSD number to process (default: all)")
    parser.add_argument("--cleanup", action="store_true",
                        help="Delete raw FASTQs and BAM after successful h5ad write")
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
            if args.cleanup:
                cleanup_raw(osd_n)
            continue

        print(f"\n=== OSD-{osd_n} ({tissue}) ===")
        try:
            adata = load_starsolo(osd_n, tissue)
            adata = basic_filter(adata)
            adata.write_h5ad(out_path)
            h5ad_gb = out_path.stat().st_size / 1e9
            print(f"  Saved: {out_path} ({h5ad_gb:.2f} GB)")
            if args.cleanup:
                cleanup_raw(osd_n)
        except Exception as e:
            print(f"  ERROR: {e}")
            continue

    print("\nDone.")


if __name__ == "__main__":
    main()

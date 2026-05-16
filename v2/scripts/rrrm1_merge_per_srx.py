#!/usr/bin/env python3
"""
rrrm1_merge_per_srx.py
Merge per-SRX STARsolo matrices into per-tissue h5ads with FLT/GC condition labels.

Usage:
    python3 rrrm1_merge_per_srx.py --tissue blood
    python3 rrrm1_merge_per_srx.py --tissue eye
    python3 rrrm1_merge_per_srx.py --tissue muscle
    python3 rrrm1_merge_per_srx.py --tissue skin
    python3 rrrm1_merge_per_srx.py --all

Produces: $SCRATCH_DIR/rrrm1_scrna/OSD-{N}/OSD-{N}_{tissue}_labeled.h5ad
  - obs columns: condition (FLT/GC), age_months, source_name, srx, animal_id
  - Barcodes prefixed: {SRX}_{barcode}-1
"""

import argparse
import os
import sys
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path

SCRATCH = Path(os.environ.get("SCRATCH_DIR", "/path/to/scratch")) / "rrrm1_scrna"
MAP_CSV = Path(os.environ.get("HOME")) / "rrrm1_scrna" / "RRRM1_SRX_CONDITION_MAP.csv"

TISSUE_OSD = {
    "blood":  "OSD-918",
    "eye":    "OSD-920",
    "muscle": "OSD-924",
    "skin":   "OSD-934",
}


def load_starsolo_matrix(starsolo_dir: Path, srx: str) -> ad.AnnData:
    """Load STARsolo GeneFull filtered matrix and prefix barcodes with SRX."""
    import gzip
    import scipy.io as sio

    matrix_dir = starsolo_dir / "Solo.out" / "GeneFull" / "filtered"
    if not matrix_dir.exists():
        raise FileNotFoundError(f"No filtered matrix at {matrix_dir}")

    # STARsolo features.tsv has 2 columns (gene_id, gene_name) — not 3 like CellRanger.
    # Use manual loading instead of sc.read_10x_mtx to handle this.
    def read_lines(path):
        if path.with_suffix(path.suffix + ".gz").exists():
            with gzip.open(str(path) + ".gz", "rt") as f:
                return [line.strip().split("\t") for line in f]
        elif path.exists():
            with open(path) as f:
                return [line.strip().split("\t") for line in f]
        raise FileNotFoundError(f"Not found: {path} or {path}.gz")

    features = read_lines(matrix_dir / "features.tsv")
    barcodes = read_lines(matrix_dir / "barcodes.tsv")

    gene_ids = [f[0] for f in features]
    gene_names = [f[1] if len(f) > 1 else f[0] for f in features]
    bc_list = [b[0] for b in barcodes]

    # Read matrix (try .gz first, then uncompressed)
    mtx_gz = matrix_dir / "matrix.mtx.gz"
    mtx_plain = matrix_dir / "matrix.mtx"
    if mtx_gz.exists():
        with gzip.open(str(mtx_gz), "rb") as f:
            mat = sio.mmread(f).T.tocsr()  # genes × cells → cells × genes
    elif mtx_plain.exists():
        mat = sio.mmread(str(mtx_plain)).T.tocsr()
    else:
        raise FileNotFoundError(f"No matrix.mtx or matrix.mtx.gz in {matrix_dir}")

    adata = ad.AnnData(X=mat)
    adata.obs_names = [f"{srx}_{bc}" for bc in bc_list]
    adata.var_names = gene_ids
    adata.var["gene_symbols"] = gene_names

    return adata


def merge_tissue(tissue: str, srx_map: pd.DataFrame, overwrite: bool = False) -> Path:
    osd = TISSUE_OSD[tissue]
    osd_n = osd.replace("OSD-", "")
    out_path = SCRATCH / osd / f"{osd}_{tissue}_labeled.h5ad"

    if out_path.exists() and not overwrite:
        print(f"[{tissue}] Output already exists: {out_path}  (use --overwrite to redo)")
        return out_path

    tissue_rows = srx_map[srx_map["tissue"] == tissue]
    print(f"\n[{tissue}] Processing {len(tissue_rows)} SRX samples from {osd}")

    adatas = []
    for _, row in tissue_rows.iterrows():
        srx = row["srx"]
        starsolo_dir = SCRATCH / osd / "starsolo_per_srx" / srx

        if not starsolo_dir.exists():
            print(f"  MISSING: {starsolo_dir} — skipping {srx}")
            continue

        print(f"  Loading {srx} ({row['condition']} {row['age_months']}mo)...", end=" ", flush=True)
        try:
            adata = load_starsolo_matrix(starsolo_dir, srx)
        except FileNotFoundError as e:
            print(f"\n  ERROR: {e}")
            continue

        # Add per-animal metadata
        adata.obs["srx"] = srx
        adata.obs["animal_id"] = row["source_name"]
        adata.obs["condition"] = row["condition"]
        adata.obs["age_months"] = int(row["age_months"])
        adata.obs["osd"] = osd
        adata.obs["tissue"] = tissue

        print(f"{adata.n_obs} cells, {adata.n_vars} genes")
        adatas.append(adata)

    if not adatas:
        print(f"[{tissue}] ERROR: No valid SRX matrices found. Check STARsolo outputs.")
        sys.exit(1)

    print(f"[{tissue}] Concatenating {len(adatas)} animals...", flush=True)
    merged = ad.concat(adatas, join="outer", label=None, merge="first")

    # Fill any NaN that appeared from outer join
    for col in ["condition", "tissue", "osd", "srx", "animal_id"]:
        if col in merged.obs.columns:
            merged.obs[col] = merged.obs[col].fillna("unknown")
    merged.obs["age_months"] = merged.obs["age_months"].fillna(0).astype(int)

    print(f"[{tissue}] Total cells: {merged.n_obs}  genes: {merged.n_vars}")
    print(f"[{tissue}] Condition counts:\n{merged.obs['condition'].value_counts().to_string()}")

    # Basic QC metrics
    merged.var["mt"] = merged.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(
        merged, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    print(f"[{tissue}] Saving to {out_path}...", flush=True)
    merged.write_h5ad(out_path)
    print(f"[{tissue}] Done.")
    return out_path


def main():
    parser = argparse.ArgumentParser(description="Merge per-SRX STARsolo matrices")
    parser.add_argument("--tissue", choices=list(TISSUE_OSD.keys()),
                        help="Single tissue to process")
    parser.add_argument("--all", action="store_true", help="Process all tissues")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output")
    parser.add_argument("--map_csv", default=str(MAP_CSV),
                        help=f"SRX condition map CSV (default: {MAP_CSV})")
    args = parser.parse_args()

    if not args.tissue and not args.all:
        parser.error("Specify --tissue <name> or --all")

    srx_map = pd.read_csv(args.map_csv)
    print(f"Loaded {len(srx_map)} SRX entries from {args.map_csv}")

    tissues = list(TISSUE_OSD.keys()) if args.all else [args.tissue]
    for tissue in tissues:
        merge_tissue(tissue, srx_map, overwrite=args.overwrite)

    print("\nAll done.")


if __name__ == "__main__":
    main()

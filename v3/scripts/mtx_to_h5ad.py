#!/usr/bin/env python3
"""Convert MTX + metadata files (from R extraction) to h5ad for RRRM-2.

Usage: python mtx_to_h5ad.py <GLDS_NUMBER>
"""
import sys
import os
import numpy as np
import pandas as pd
import scipy.io
import anndata as ad

def main():
    if len(sys.argv) < 2:
        print("Usage: python mtx_to_h5ad.py <GLDS_NUMBER>")
        sys.exit(1)

    glds = sys.argv[1]
    data_dir = "/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/data/rrrm2"
    prefix = os.path.join(data_dir, f"GLDS_{glds}")

    print(f"=== Converting GLDS-{glds} to h5ad ===")

    # Load counts (genes x cells in MTX format)
    print("Loading counts matrix...")
    counts = scipy.io.mmread(f"{prefix}_counts.mtx")
    counts = counts.tocsc()  # genes x cells
    print(f"  Shape: {counts.shape}")

    # Load gene names and barcodes
    genes = [g.strip() for g in open(f"{prefix}_genes.txt")]
    barcodes = [b.strip() for b in open(f"{prefix}_barcodes.txt")]
    print(f"  Genes: {len(genes)}, Barcodes: {len(barcodes)}")

    # Transpose to cells x genes for AnnData
    counts_T = counts.T.tocsr()

    # Load metadata
    meta = pd.read_csv(f"{prefix}_metadata.csv", index_col=0)
    print(f"  Metadata: {meta.shape[0]} cells, {meta.shape[1]} columns")
    print(f"  Columns: {list(meta.columns)}")

    # Verify alignment
    assert len(barcodes) == counts_T.shape[0], f"Barcode count mismatch: {len(barcodes)} vs {counts_T.shape[0]}"
    assert len(genes) == counts_T.shape[1], f"Gene count mismatch: {len(genes)} vs {counts_T.shape[1]}"

    # Create AnnData
    adata = ad.AnnData(
        X=counts_T,
        obs=meta.loc[barcodes],  # Align metadata to barcode order
        var=pd.DataFrame(index=genes),
    )

    # Add embeddings if available
    for emb_name in ["umap", "pca", "harmony", "harmony_pca"]:
        emb_path = f"{prefix}_{emb_name}.csv"
        if os.path.exists(emb_path):
            print(f"  Adding embedding: {emb_name}")
            emb_df = pd.read_csv(emb_path, index_col=0)
            # Align to adata.obs_names
            emb_df = emb_df.loc[adata.obs_names]
            adata.obsm[f"X_{emb_name}"] = emb_df.values

    # Add dataset info
    adata.uns["glds"] = f"GLDS-{glds}"
    adata.uns["mission"] = "RRRM-2"

    # Summary
    print(f"\nAnnData summary:")
    print(f"  Shape: {adata.shape}")
    if "exp" in adata.obs.columns:
        print(f"  Conditions: {adata.obs['exp'].value_counts().to_dict()}")
    if "predicted.id" in adata.obs.columns:
        print(f"  Cell types: {adata.obs['predicted.id'].nunique()}")
        print(f"  {adata.obs['predicted.id'].value_counts().to_dict()}")

    # Save
    out_path = os.path.join(data_dir, f"GLDS_{glds}.h5ad")
    print(f"\nSaving to {out_path}...")
    adata.write_h5ad(out_path)
    print(f"  File size: {os.path.getsize(out_path) / 1e9:.2f} GB")
    print("Done!")


if __name__ == "__main__":
    main()

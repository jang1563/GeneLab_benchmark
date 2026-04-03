#!/usr/bin/env python3
"""
rrrm1_initial_scanpy.py
Initial exploratory Scanpy workflow for RRRM-1 scRNA-seq.

Outputs:
- global exploratory UMAP from the merged object
- tissue-specific processed h5ad objects
- tissue-specific marker tables
- summary CSV tables and PNG plots
"""

import os
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans

SCRATCH = Path(os.environ.get("SCRATCH_DIR", ".")) / "rrrm1_scrna"
INPUT_H5AD = SCRATCH / "RRRM1_merged_samples.h5ad"
OUTDIR = SCRATCH / "downstream_initial"
FIGDIR = OUTDIR / "figures"
TABLEDIR = OUTDIR / "tables"
OBJDIR = OUTDIR / "objects"
TISSUES = ["blood", "eye", "muscle", "skin"]


def ensure_dirs() -> None:
    for path in [OUTDIR, FIGDIR, TABLEDIR, OBJDIR]:
        path.mkdir(parents=True, exist_ok=True)


def preprocess(
    adata,
    batch_key: str | None,
    n_top_genes: int,
    resolution: float,
):
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    if batch_key is not None and batch_key in adata.obs:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            batch_key=batch_key,
            flavor="seurat",
            subset=True,
        )
    else:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor="seurat",
            subset=True,
        )

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    n_pcs = min(30, adata.obsm["X_pca"].shape[1])
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
    sc.tl.umap(adata)
    assign_clusters(adata, resolution=resolution)
    return adata


def assign_clusters(adata, resolution: float) -> None:
    """Use Leiden when available, otherwise fall back to PCA KMeans."""
    try:
        sc.tl.leiden(adata, resolution=resolution, key_added="cluster")
        adata.uns["cluster_method"] = "leiden"
    except ImportError:
        # Keep clustering dependency-light on Cayuga where igraph may be absent.
        n_clusters = max(4, min(12, round(resolution * 10)))
        km = KMeans(n_clusters=n_clusters, random_state=0, n_init=10)
        labels = km.fit_predict(adata.obsm["X_pca"][:, : min(20, adata.obsm["X_pca"].shape[1])])
        adata.obs["cluster"] = pd.Categorical(labels.astype(str))
        adata.uns["cluster_method"] = f"kmeans_{n_clusters}"


def save_global_outputs(adata) -> None:
    summary = (
        adata.obs.groupby(["tissue", "osd"], observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["tissue", "osd"])
    )
    summary.to_csv(TABLEDIR / "RRRM1_global_cell_counts.csv", index=False)

    qc_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    qc_summary = (
        adata.obs.groupby(["tissue", "osd"], observed=True)[qc_cols]
        .median()
        .reset_index()
        .sort_values(["tissue", "osd"])
    )
    qc_summary.to_csv(TABLEDIR / "RRRM1_global_qc_medians.csv", index=False)

    sc.pl.umap(
        adata,
        color=["tissue"],
        show=False,
        save=False,
    )
    plt.savefig(FIGDIR / "RRRM1_global_umap_by_tissue.png", dpi=200, bbox_inches="tight")
    plt.close()

    sc.pl.umap(
        adata,
        color=["osd"],
        show=False,
        save=False,
    )
    plt.savefig(FIGDIR / "RRRM1_global_umap_by_osd.png", dpi=200, bbox_inches="tight")
    plt.close()

    sc.pl.umap(
        adata,
        color=["cluster"],
        legend_loc="on data",
        show=False,
        save=False,
    )
    plt.savefig(FIGDIR / "RRRM1_global_umap_by_cluster.png", dpi=200, bbox_inches="tight")
    plt.close()

    adata.write_h5ad(OBJDIR / "RRRM1_global_exploratory.h5ad")


def marker_table(adata, groupby: str, top_n: int = 20) -> pd.DataFrame:
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    rows = []
    for group in groups:
        names = result["names"][group][:top_n]
        scores = result["scores"][group][:top_n]
        pvals_adj = result["pvals_adj"][group][:top_n]
        logfc = result["logfoldchanges"][group][:top_n]
        for gene, score, padj, lfc in zip(names, scores, pvals_adj, logfc):
            rows.append(
                {
                    groupby: group,
                    "gene": gene,
                    "score": float(score),
                    "pvals_adj": float(padj),
                    "logfoldchange": float(lfc),
                }
            )
    return pd.DataFrame(rows)


def save_tissue_outputs(adata) -> None:
    for tissue in TISSUES:
        tissue_adata = adata[adata.obs["tissue"] == tissue].copy()
        if tissue_adata.n_obs == 0:
            continue

        tissue_proc = preprocess(
            tissue_adata,
            batch_key=None,
            n_top_genes=min(3000, tissue_adata.n_vars),
            resolution=0.4,
        )

        sc.tl.rank_genes_groups(
            tissue_proc,
            groupby="cluster",
            method="wilcoxon",
            key_added="rank_genes_groups",
        )

        tissue_proc.write_h5ad(OBJDIR / f"RRRM1_{tissue}_processed.h5ad")

        counts = (
            tissue_proc.obs.groupby("cluster", observed=True)
            .size()
            .reset_index(name="n_cells")
            .sort_values("cluster")
        )
        counts.to_csv(TABLEDIR / f"RRRM1_{tissue}_cluster_counts.csv", index=False)

        markers = marker_table(tissue_proc, "cluster")
        markers.to_csv(TABLEDIR / f"RRRM1_{tissue}_cluster_markers_top20.csv", index=False)

        sc.pl.umap(
            tissue_proc,
            color=["cluster"],
            legend_loc="on data",
            show=False,
            save=False,
        )
        plt.savefig(FIGDIR / f"RRRM1_{tissue}_umap_by_cluster.png", dpi=200, bbox_inches="tight")
        plt.close()

        sc.pl.umap(
            tissue_proc,
            color=["total_counts", "n_genes_by_counts", "pct_counts_mt"],
            show=False,
            save=False,
        )
        plt.savefig(FIGDIR / f"RRRM1_{tissue}_umap_qc.png", dpi=200, bbox_inches="tight")
        plt.close()


def main() -> None:
    ensure_dirs()
    if not INPUT_H5AD.exists():
        raise FileNotFoundError(f"Missing merged input: {INPUT_H5AD}")

    adata = ad.read_h5ad(INPUT_H5AD)
    global_proc = preprocess(adata, batch_key="osd", n_top_genes=3000, resolution=0.6)
    save_global_outputs(global_proc)
    save_tissue_outputs(adata)

    print(f"Wrote outputs under: {OUTDIR}")


if __name__ == "__main__":
    main()

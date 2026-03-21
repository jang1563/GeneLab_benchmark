#!/usr/bin/env python3
"""Annotate GLDS-404 (RRRM-2 PBMC) cell types using marker genes.

GLDS-404 has seurat_clusters but no predicted.id (cell type labels).
This script assigns cell types based on marker gene expression scoring
per cluster, then saves the updated h5ad.

Mouse PBMC expected cell types + markers:
  T cell:       Cd3d, Cd3e, Cd3g
  CD4+ T:       Cd4
  CD8+ T:       Cd8a, Cd8b1
  B cell:       Cd79a, Ms4a1, Cd19
  NK cell:      Nkg7, Klrb1c, Gzma
  Monocyte:     Cd14, Lyz2, Csf1r
  Dendritic:    Itgax, Flt3, H2-Eb1
  Neutrophil:   S100a8, S100a9, Ly6g
  Erythroid:    Hbb-bt, Hba-a1, Hba-a2
  Platelet:     Pf4, Ppbp

Usage:
  python annotate_glds404_pbmc.py [--data-dir DIR] [--min-cells 10]
"""
import argparse
import json
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.filterwarnings("ignore", category=FutureWarning)

# Default HPC path
DEFAULT_DATA_DIR = Path(
    "/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/"
    "GeneLab_benchmark/v3/data/rrrm2"
)

# Mouse PBMC marker genes per cell type
MARKERS = {
    "T_cell": ["Cd3d", "Cd3e", "Cd3g"],
    "B_cell": ["Cd79a", "Ms4a1", "Cd19"],
    "NK_cell": ["Nkg7", "Klrb1c", "Gzma"],
    "monocyte": ["Cd14", "Lyz2", "Csf1r"],
    "dendritic": ["Itgax", "Flt3", "H2-Eb1"],
    "neutrophil": ["S100a8", "S100a9"],
    "erythroid": ["Hbb-bt", "Hba-a1", "Hba-a2"],
    "platelet": ["Pf4", "Ppbp"],
}

# Sub-type markers (applied to T_cell cluster only)
T_SUBTYPES = {
    "CD4_T_cell": ["Cd4"],
    "CD8_T_cell": ["Cd8a", "Cd8b1"],
}


def build_gene_index(adata) -> dict:
    """Build mapping from gene symbol → var index.

    Checks var_names and var["gene_symbols"/"gene_name"] columns.
    Returns: {symbol_lower: var_index}
    """
    gene_idx = {}

    # Primary: var_names as symbols
    for i, name in enumerate(adata.var_names):
        gene_idx[str(name).lower()] = i

    # Secondary: check for gene symbol columns
    for col in ["gene_symbols", "gene_name", "gene_short_name"]:
        if col in adata.var.columns:
            for i, sym in enumerate(adata.var[col].astype(str)):
                if sym and sym.lower() != "nan":
                    gene_idx[sym.lower()] = i
            break

    return gene_idx


def get_gene_vector(X, idx):
    """Extract a single gene column from X, return as 1D dense array."""
    col = X[:, idx]
    if sp.issparse(col):
        col = col.toarray().flatten()
    return np.asarray(col, dtype=np.float32).flatten()


def score_markers(adata, gene_idx: dict, markers: list) -> np.ndarray:
    """Score cells by mean log1p expression of marker genes.

    Returns: 1D array of scores (n_cells,).
    """
    X = adata.X
    found_indices = []
    for gene in markers:
        idx = gene_idx.get(gene.lower())
        if idx is not None:
            found_indices.append(idx)

    if not found_indices:
        return np.zeros(adata.n_obs)

    scores = np.zeros(X.shape[0], dtype=np.float32)
    for idx in found_indices:
        col = get_gene_vector(X, idx)
        scores += np.log1p(col)
    scores /= len(found_indices)

    return scores


def score_cluster_zscore(adata, gene_idx: dict, cluster_col: str) -> dict:
    """Score each cluster using z-scored log1p expression + detection rate.

    Steps:
    1. For each marker gene, compute per-cluster: mean log1p expression + detection rate
    2. Z-score across clusters per gene (so high-baseline genes don't dominate)
    3. For each cell type, score = mean z-scored expression × detection weight

    Returns: {cluster_id: {cell_type: combined_score}}
    """
    X = adata.X
    obs = adata.obs
    clusters = obs[cluster_col].astype(str)
    unique_cl = sorted(clusters.unique(), key=lambda x: int(x) if x.isdigit() else x)

    # Collect all unique marker genes
    all_markers = {}
    for ct, markers in MARKERS.items():
        for gene in markers:
            idx = gene_idx.get(gene.lower())
            if idx is not None:
                all_markers[gene] = idx

    if not all_markers:
        return {cl: {ct: 0.0 for ct in MARKERS} for cl in unique_cl}

    # Compute per-cluster stats for each marker gene
    # gene → {cluster: (mean_log1p, detection_rate)}
    gene_cluster_stats = {}
    for gene, idx in all_markers.items():
        col = get_gene_vector(X, idx)
        stats_per_cl = {}
        for cl in unique_cl:
            mask = (clusters == cl).values
            vals = col[mask]
            mean_log1p = float(np.log1p(vals).mean())
            det_rate = float((vals > 0).mean())
            stats_per_cl[cl] = (mean_log1p, det_rate)
        gene_cluster_stats[gene] = stats_per_cl

    # Z-score the mean_log1p across clusters per gene
    gene_cluster_zscore = {}
    for gene, stats_per_cl in gene_cluster_stats.items():
        means = np.array([stats_per_cl[cl][0] for cl in unique_cl])
        det_rates = np.array([stats_per_cl[cl][1] for cl in unique_cl])
        mu = means.mean()
        sigma = means.std()
        if sigma < 1e-10:
            zscores = np.zeros_like(means)
        else:
            zscores = (means - mu) / sigma

        gene_cluster_zscore[gene] = {
            cl: (float(z), float(d))
            for cl, z, d in zip(unique_cl, zscores, det_rates)
        }

    # Score each cluster for each cell type
    cluster_scores = {}
    for cl in unique_cl:
        ct_scores = {}
        for ct, markers in MARKERS.items():
            found_genes = [g for g in markers if g in gene_cluster_zscore]
            if not found_genes:
                ct_scores[ct] = 0.0
                continue
            # Combined score = mean(zscore × detection_rate) across markers
            total = 0.0
            for gene in found_genes:
                z, d = gene_cluster_zscore[gene][cl]
                total += z * d  # weight z-score by detection rate
            ct_scores[ct] = total / len(found_genes)
        cluster_scores[cl] = ct_scores

    return cluster_scores


def annotate_per_cell(adata, gene_idx: dict) -> pd.Series:
    """Score each cell for all marker sets, assign to highest-scoring type."""
    print("  Scoring all cells...")
    all_scores = {}
    for ct, markers in MARKERS.items():
        scores = score_markers(adata, gene_idx, markers)
        all_scores[ct] = scores
        found = sum(1 for g in markers if g.lower() in gene_idx)
        print(f"    {ct}: {found}/{len(markers)} markers found, "
              f"mean={scores.mean():.3f}, max={scores.max():.3f}")

    # Build score matrix (cells × cell_types)
    score_df = pd.DataFrame(all_scores, index=adata.obs_names)

    # Assign each cell to highest-scoring type
    # Cells with max score = 0 are "unassigned"
    max_scores = score_df.max(axis=1)
    assignments = score_df.idxmax(axis=1)
    assignments[max_scores == 0] = "unassigned"

    return assignments


def annotate_per_cluster(adata, gene_idx: dict, cluster_col: str) -> pd.Series:
    """Score each cluster using z-scored expression + detection rate.

    Uses score_cluster_zscore to avoid bias from high-baseline genes
    (e.g., hemoglobin in PBMCs).
    """
    clusters = adata.obs[cluster_col].astype(str)
    unique_cl = sorted(clusters.unique(), key=lambda x: int(x) if x.isdigit() else x)
    print(f"  Annotating {len(unique_cl)} clusters from '{cluster_col}'...")
    print(f"  Using z-scored log1p expression × detection rate")

    # Get z-scored cluster scores
    all_ct_scores = score_cluster_zscore(adata, gene_idx, cluster_col)

    cluster_info = {}
    for cl in unique_cl:
        n_cells = (clusters == cl).sum()
        ct_scores = all_ct_scores[cl]
        best_ct = max(ct_scores, key=ct_scores.get)
        best_score = ct_scores[best_ct]

        # Only assign if score is positive (above average for that marker set)
        if best_score <= 0:
            best_ct = "unassigned"

        cluster_info[cl] = {
            "assigned_type": best_ct,
            "best_score": best_score,
            "n_cells": int(n_cells),
            "all_scores": ct_scores,
        }

        # Show top 3 scores for each cluster
        sorted_scores = sorted(ct_scores.items(), key=lambda x: x[1], reverse=True)
        top3 = ", ".join(f"{ct}={s:.3f}" for ct, s in sorted_scores[:3])
        print(f"    Cluster {cl}: {best_ct} (n={n_cells}) [{top3}]")

    # Map back to cells
    cl_to_ct = {cl: info["assigned_type"] for cl, info in cluster_info.items()}
    assignments = adata.obs[cluster_col].astype(str).map(cl_to_ct)

    return assignments, cluster_info


def refine_t_subtypes(adata, assignments: pd.Series, gene_idx: dict) -> pd.Series:
    """Subdivide T_cell into CD4 and CD8 based on subtype markers."""
    t_mask = assignments == "T_cell"
    n_t = t_mask.sum()
    if n_t == 0:
        return assignments

    print(f"  Refining {n_t} T cells into CD4/CD8 subtypes...")
    t_sub = adata[t_mask]

    cd4_scores = score_markers(t_sub, gene_idx, T_SUBTYPES["CD4_T_cell"])
    cd8_scores = score_markers(t_sub, gene_idx, T_SUBTYPES["CD8_T_cell"])

    subtypes = []
    for cd4, cd8 in zip(cd4_scores, cd8_scores):
        if cd4 > cd8 and cd4 > 0:
            subtypes.append("CD4_T_cell")
        elif cd8 > cd4 and cd8 > 0:
            subtypes.append("CD8_T_cell")
        else:
            subtypes.append("T_cell")

    result = assignments.copy()
    result.loc[t_mask] = subtypes

    for st in ["CD4_T_cell", "CD8_T_cell", "T_cell"]:
        n = (result == st).sum()
        print(f"    {st}: {n}")

    return result


def main():
    parser = argparse.ArgumentParser(description="Annotate GLDS-404 PBMC cell types")
    parser.add_argument("--data-dir", type=str, default=str(DEFAULT_DATA_DIR))
    parser.add_argument("--method", choices=["per_cell", "per_cluster"], default="per_cluster")
    parser.add_argument("--cluster-col", default="seurat_clusters",
                        help="Cluster column to use for per-cluster annotation")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    h5ad_path = data_dir / "GLDS_404.h5ad"

    if not h5ad_path.exists():
        print(f"ERROR: {h5ad_path} not found")
        return

    print(f"=== GLDS-404 PBMC Cell Type Annotation ===")
    print(f"Method: {args.method}")
    print(f"Loading {h5ad_path}...")
    adata = ad.read_h5ad(h5ad_path)
    print(f"Shape: {adata.shape}")
    print(f"obs columns: {list(adata.obs.columns)}")

    # Build gene index
    gene_idx = build_gene_index(adata)
    print(f"Gene index: {len(gene_idx)} entries")

    # Check marker coverage
    total_markers = sum(len(v) for v in MARKERS.values())
    found_markers = sum(1 for ct, ms in MARKERS.items()
                        for m in ms if m.lower() in gene_idx)
    print(f"Marker coverage: {found_markers}/{total_markers}")

    # Annotate
    if args.method == "per_cluster":
        if args.cluster_col not in adata.obs.columns:
            print(f"ERROR: Cluster column '{args.cluster_col}' not found")
            print(f"Available: {list(adata.obs.columns)}")
            return

        assignments, cluster_info = annotate_per_cluster(
            adata, gene_idx, args.cluster_col)

        # Save cluster annotation summary
        summary_path = data_dir / "GLDS_404_cluster_annotation.json"
        with open(summary_path, "w") as f:
            json.dump(cluster_info, f, indent=2)
        print(f"Saved cluster summary: {summary_path}")

    else:
        assignments = annotate_per_cell(adata, gene_idx)

    # Refine T cell subtypes
    assignments = refine_t_subtypes(adata, assignments, gene_idx)

    # Add to adata
    adata.obs["predicted.id"] = assignments.values
    adata.obs["broad_celltype"] = assignments.values

    # Print summary
    print(f"\n=== Cell Type Summary ===")
    for ct in sorted(assignments.unique()):
        n = (assignments == ct).sum()
        pct = n / len(assignments) * 100
        print(f"  {ct}: {n:,} ({pct:.1f}%)")

    # Save updated h5ad
    output_path = data_dir / "GLDS_404.h5ad"
    print(f"\nSaving annotated h5ad to {output_path}...")
    adata.write_h5ad(output_path)
    print("Done!")


if __name__ == "__main__":
    main()

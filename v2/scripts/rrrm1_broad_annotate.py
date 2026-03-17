#!/usr/bin/env python3
"""
rrrm1_broad_annotate.py
First-pass broad annotation for RRRM-1 tissue-specific processed objects.

Approach:
- load each tissue-specific processed h5ad
- score curated marker sets with scanpy.tl.score_genes
- aggregate marker scores by cluster
- assign each cluster to the highest-scoring broad cell type
- save annotated h5ad, score tables, cluster annotations, and UMAP figures
"""

from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

BASE = Path("/athena/masonlab/scratch/users/jak4013/rrrm1_scrna/downstream_initial")
OBJDIR = BASE / "objects"
OUTDIR = BASE / "annotations"
OUT_OBJDIR = OUTDIR / "objects"
OUT_TABLEDIR = OUTDIR / "tables"
OUT_FIGDIR = OUTDIR / "figures"

MARKERS = {
    "blood": {
        "erythroid": ["Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Alas2", "Klf1"],
        "b_cell": ["Cd79a", "Cd79b", "Ms4a1", "Cd74", "H2-Aa"],
        "t_cell": ["Cd3d", "Cd3e", "Trac", "Ltb", "Il7r"],
        "nk_cell": ["Nkg7", "Klrd1", "Klrb1c", "Prf1", "Ccl5"],
        "monocyte_macrophage": ["Lyz2", "Lgals3", "Ctss", "Tyrobp", "Fcgr3"],
        "neutrophil": ["S100a8", "S100a9", "Retnlg", "Lcn2", "Mmp8"],
        "platelet_megakaryocyte": ["Ppbp", "Pf4", "Gng11", "Itga2b", "Gp9"],
        "dendritic_cell": ["Flt3", "Ccr7", "Itgax", "H2-Ab1", "Xcr1"],
    },
    "eye": {
        "epithelial": ["Krt8", "Krt18", "Krt19", "Epcam", "Krt14"],
        "corneal_conjunctival_epithelial": ["Krt12", "Krt13", "Krt14", "Krt15", "Tacstd2"],
        "lens_crystallin": ["Cryaa", "Cryab", "Cryba1", "Crybb2", "Mip"],
        "retinal_neuronal": ["Rbfox3", "Snap25", "Tubb3", "Elavl4", "Syt1"],
        "photoreceptor_like": ["Rho", "Gnat1", "Pde6a", "Rcvrn", "Crx"],
        "muller_glia": ["Rlbp1", "Slc1a3", "Glul", "Aqp4", "Vim"],
        "stromal_fibroblast": ["Col1a1", "Col1a2", "Dcn", "Lum", "Pdgfra"],
        "endothelial_perivascular": ["Kdr", "Klf2", "Pecam1", "Cdh5", "Rgs5"],
        "immune": ["Ptprc", "Lyz2", "Tyrobp", "Ctss", "H2-Ab1"],
    },
    "muscle": {
        "myonuclear_structural": ["Acta1", "Tnnt3", "Tpm1", "Myh1", "Ckm"],
        "satellite_myogenic_progenitor": ["Pax7", "Myf5", "Myod1", "Myog", "Vcam1"],
        "fibroadipogenic_progenitor": ["Pdgfra", "Col1a1", "Dcn", "Lum", "Pi16"],
        "fibroblast_matrix": ["Col1a1", "Col3a1", "Fn1", "Dcn", "Postn"],
        "endothelial": ["Pecam1", "Cdh5", "Kdr", "Emcn", "Klf2"],
        "pericyte_smooth_muscle": ["Rgs5", "Myl9", "Acta2", "Cspg4", "Pdgfrb"],
        "macrophage_myeloid": ["Lyz2", "Adgre1", "Tyrobp", "Ctss", "Lgals3"],
        "t_nk_lymphocyte": ["Cd3d", "Cd3e", "Trac", "Nkg7", "Ccl5"],
    },
    "skin": {
        "basal_keratinocyte": ["Krt5", "Krt14", "Krt15", "Dst", "Trp63"],
        "suprabasal_keratinocyte": ["Krt1", "Krt10", "Lor", "Flg", "Sprr1b"],
        "hair_follicle_epithelial": ["Krt17", "Krt79", "Lhx2", "Sox9", "Shh"],
        "fibroblast_dermal_stromal": ["Col1a1", "Col1a2", "Dcn", "Lum", "Pdgfra"],
        "endothelial": ["Pecam1", "Cdh5", "Kdr", "Emcn", "Klf2"],
        "pericyte_vascular_smooth_muscle": ["Rgs5", "Acta2", "Myl9", "Tagln", "Cspg4"],
        "immune_myeloid": ["Lyz2", "Tyrobp", "Ctss", "Fcgr3", "Adgre1"],
        "t_nk_lymphocyte": ["Cd3d", "Cd3e", "Trac", "Nkg7", "Ccl5"],
        "melanocyte_like": ["Pmel", "Dct", "Mlana", "Tyr", "Tyrp1"],
    },
}


def ensure_dirs() -> None:
    for path in [OUTDIR, OUT_OBJDIR, OUT_TABLEDIR, OUT_FIGDIR]:
        path.mkdir(parents=True, exist_ok=True)


def gene_universe(adata) -> set[str]:
    if adata.raw is not None:
        return set(map(str, adata.raw.var_names))
    return set(map(str, adata.var_names))


def score_marker_sets(adata, marker_sets: dict[str, list[str]]) -> list[str]:
    available = gene_universe(adata)
    score_cols = []
    for label, genes in marker_sets.items():
        present = [g for g in genes if g in available]
        if len(present) < 2:
            continue
        score_name = f"score_{label}"
        sc.tl.score_genes(
            adata,
            gene_list=present,
            score_name=score_name,
            use_raw=adata.raw is not None,
            random_state=0,
        )
        score_cols.append(score_name)
    return score_cols


def annotate_tissue(tissue: str) -> None:
    in_path = OBJDIR / f"RRRM1_{tissue}_processed.h5ad"
    if not in_path.exists():
        raise FileNotFoundError(f"Missing processed object: {in_path}")

    adata = ad.read_h5ad(in_path)
    score_cols = score_marker_sets(adata, MARKERS[tissue])
    if not score_cols:
        raise RuntimeError(f"No marker sets could be scored for tissue {tissue}")

    # ── Cluster-level scores (for diagnostics / high-confidence clusters) ──
    cluster_scores = (
        adata.obs.groupby("cluster", observed=True)[score_cols]
        .median()
        .sort_index()
    )
    renamed = {col: col.replace("score_", "") for col in score_cols}
    cluster_scores = cluster_scores.rename(columns=renamed)

    cluster_annotation = cluster_scores.idxmax(axis=1).rename("broad_celltype").to_frame()
    cluster_annotation["top_score"] = cluster_scores.max(axis=1)
    cluster_annotation["second_score"] = cluster_scores.apply(
        lambda row: row.nlargest(2).iloc[-1] if len(row) > 1 else row.iloc[0],
        axis=1,
    )
    cluster_annotation["score_margin"] = (
        cluster_annotation["top_score"] - cluster_annotation["second_score"]
    )
    cluster_annotation = cluster_annotation.reset_index()

    # ── Per-cell scoring (primary): robust to condition-driven clustering ──
    # Cluster-level argmax fails when FLT/GC form separate clusters due to batch
    # effects, causing entire FLT clusters to be misannotated (e.g., skin P0 bug).
    # Per-cell scoring assigns each cell based on its own marker expression.
    cell_scores = adata.obs[score_cols].copy()
    cell_scores.columns = [c.replace("score_", "") for c in cell_scores.columns]
    per_cell_label = cell_scores.idxmax(axis=1)

    # Use per-cell scoring as the primary annotation
    adata.obs["broad_celltype"] = per_cell_label.astype("category")

    # Keep cluster-level annotation for comparison / diagnostics
    cluster_map = dict(zip(cluster_annotation["cluster"], cluster_annotation["broad_celltype"]))
    adata.obs["broad_celltype_cluster"] = (
        adata.obs["cluster"].map(cluster_map).astype("category")
    )

    # Report agreement between per-cell and cluster-level
    agree = (adata.obs["broad_celltype"] == adata.obs["broad_celltype_cluster"]).sum()
    print(f"  {tissue}: per-cell vs cluster agreement: "
          f"{agree}/{len(adata)} ({agree/len(adata)*100:.1f}%)")

    celltype_counts = (
        adata.obs.groupby("broad_celltype", observed=True)
        .size()
        .reset_index(name="n_cells")
        .sort_values("n_cells", ascending=False)
    )

    score_table = cluster_scores.reset_index()
    score_table.to_csv(OUT_TABLEDIR / f"RRRM1_{tissue}_cluster_score_matrix.csv", index=False)
    cluster_annotation.to_csv(
        OUT_TABLEDIR / f"RRRM1_{tissue}_cluster_annotation.csv",
        index=False,
    )
    celltype_counts.to_csv(
        OUT_TABLEDIR / f"RRRM1_{tissue}_broad_celltype_counts.csv",
        index=False,
    )

    sc.pl.umap(
        adata,
        color=["broad_celltype"],
        legend_loc="right margin",
        show=False,
        save=False,
    )
    plt.savefig(
        OUT_FIGDIR / f"RRRM1_{tissue}_umap_broad_celltype.png",
        dpi=200,
        bbox_inches="tight",
    )
    plt.close()

    adata.write_h5ad(OUT_OBJDIR / f"RRRM1_{tissue}_annotated.h5ad")


def main() -> None:
    ensure_dirs()
    for tissue in MARKERS:
        annotate_tissue(tissue)
    print(f"Wrote broad annotation outputs under: {OUTDIR}")


if __name__ == "__main__":
    main()

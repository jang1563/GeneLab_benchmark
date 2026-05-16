#!/usr/bin/env python3
"""
rrrm1_singlecell_hardening.py

Benchmark-oriented single-cell hardening for RRRM-1:
- optional Scrublet doublet review on per-sample count-space h5ad
- conservative subtype refinement from curated marker sets
- ambient/leakage review from marker-set enrichment outside expected labels

Outputs are written under:
  $SCRATCH_DIR/rrrm1_scrna/downstream_initial/hardening
"""

import os
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

SCRATCH = Path(os.environ.get("SCRATCH_DIR", ".")) / "rrrm1_scrna"
ANNOT_BASE = SCRATCH / "downstream_initial" / "annotations"
OUTDIR = SCRATCH / "downstream_initial" / "hardening"
TABLEDIR = OUTDIR / "tables"
FIGDIR = OUTDIR / "figures"
OBJDIR = OUTDIR / "objects"

SAMPLE_H5AD = {
    "blood": SCRATCH / "OSD-918" / "OSD-918_blood_rrrm1.h5ad",
    "eye": SCRATCH / "OSD-920" / "OSD-920_eye_rrrm1.h5ad",
    "muscle": SCRATCH / "OSD-924" / "OSD-924_muscle_rrrm1.h5ad",
    "skin": SCRATCH / "OSD-934" / "OSD-934_skin_rrrm1.h5ad",
}

ANNOT_H5AD = {
    tissue: ANNOT_BASE / "objects" / f"RRRM1_{tissue}_annotated.h5ad"
    for tissue in ["blood", "eye", "muscle", "skin"]
}

SUBTYPE_MARKERS = {
    "blood": {
        "mature_erythroid": ["Hba-a1", "Hba-a2", "Hbb-bt", "Hbb-bs", "Klf1"],
        "b_cell_antigen_presentation": ["Cd79a", "Cd79b", "Ms4a1", "Cd74", "H2-Aa"],
        "t_cell": ["Cd3d", "Cd3e", "Trac", "Il7r", "Ltb"],
        "nk_cytotoxic": ["Nkg7", "Prf1", "Ccl5", "Klrd1"],
        "neutrophil_inflammatory": ["S100a8", "S100a9", "Retnlg", "Lcn2"],
        "monocyte_macrophage": ["Lyz2", "Ctss", "Tyrobp", "Lgals3"],
    },
    "eye": {
        "retinal_neuron": ["Rbfox3", "Snap25", "Tubb3", "Elavl4"],
        "muller_glia": ["Rlbp1", "Slc1a3", "Glul", "Aqp4"],
        "corneal_epithelium": ["Krt12", "Krt13", "Krt14", "Tacstd2"],
        "photoreceptor": ["Rho", "Gnat1", "Pde6a", "Rcvrn", "Crx"],
        "stromal": ["Col1a1", "Col1a2", "Dcn", "Lum"],
        "immune": ["Ptprc", "Tyrobp", "Ctss", "Lyz2"],
    },
    "muscle": {
        "macrophage_inflammatory": ["Lyz2", "Tyrobp", "Ctss", "Lgals3"],
        "t_nk_lymphocyte": ["Cd3d", "Cd3e", "Trac", "Nkg7", "Ccl5"],
        "fap": ["Pdgfra", "Col1a1", "Dcn", "Pi16"],
        "endothelial": ["Pecam1", "Cdh5", "Kdr", "Emcn"],
        "satellite_myogenic": ["Pax7", "Myf5", "Myod1", "Vcam1"],
        "myonuclear": ["Acta1", "Tnnt3", "Myh1", "Ckm"],
    },
    "skin": {
        "basal_keratinocyte": ["Krt5", "Krt14", "Krt15", "Trp63"],
        "suprabasal_keratinocyte": ["Krt1", "Krt10", "Lor", "Flg"],
        "hair_follicle": ["Krt17", "Krt79", "Lhx2", "Sox9"],
        "immune_myeloid": ["Lyz2", "Tyrobp", "Ctss", "Adgre1"],
        "t_nk_lymphocyte": ["Cd3d", "Cd3e", "Trac", "Nkg7"],
        "melanocyte": ["Pmel", "Dct", "Mlana", "Tyr"],
    },
}


def ensure_dirs() -> None:
    for path in [OUTDIR, TABLEDIR, FIGDIR, OBJDIR]:
        path.mkdir(parents=True, exist_ok=True)


def available_genes(adata) -> set[str]:
    return set(map(str, adata.var_names))


def run_scrublet_if_available(adata_counts):
    try:
        import scrublet as scr
    except Exception:
        return None

    X = adata_counts.X
    if sparse.issparse(X):
        X = X.tocsc()
    scrub = scr.Scrublet(X)
    scores, preds = scrub.scrub_doublets()
    return pd.DataFrame(
        {
            "cell_barcode": adata_counts.obs_names.astype(str),
            "scrublet_score": scores,
            "predicted_doublet": preds.astype(bool),
        }
    )


def subtype_and_leakage(tissue: str, adata_annot):
    marker_sets = SUBTYPE_MARKERS[tissue]
    genes = available_genes(adata_annot.raw.to_adata() if adata_annot.raw is not None else adata_annot)
    score_cols = []
    for label, markers in marker_sets.items():
        use_markers = [g for g in markers if g in genes]
        if len(use_markers) < 2:
            continue
        score_name = f"subtype_score_{label}"
        sc.tl.score_genes(
            adata_annot,
            gene_list=use_markers,
            score_name=score_name,
            use_raw=adata_annot.raw is not None,
            random_state=0,
        )
        score_cols.append(score_name)

    cluster_scores = (
        adata_annot.obs.groupby("cluster", observed=True)[score_cols]
        .median()
        .rename(columns=lambda c: c.replace("subtype_score_", ""))
        .reset_index()
    )

    dominant = cluster_scores.set_index("cluster")
    dominant_label = dominant.idxmax(axis=1).rename("refined_subtype")
    dominant_score = dominant.max(axis=1).rename("top_score")
    second_score = dominant.apply(
        lambda row: row.nlargest(2).iloc[-1] if len(row) > 1 else row.iloc[0], axis=1
    ).rename("second_score")
    subtype_table = pd.concat([dominant_label, dominant_score, second_score], axis=1).reset_index()
    subtype_table["score_margin"] = subtype_table["top_score"] - subtype_table["second_score"]

    mapping = dict(zip(subtype_table["cluster"], subtype_table["refined_subtype"]))
    adata_annot.obs["refined_subtype"] = adata_annot.obs["cluster"].map(mapping).astype("category")

    broad = adata_annot.obs["broad_celltype"].astype(str)
    leakage_rows = []
    for label in dominant.columns:
        cell_scores = adata_annot.obs[f"subtype_score_{label}"]
        by_broad = cell_scores.groupby(broad, observed=True).median().sort_values(ascending=False)
        top_broad = by_broad.index[0]
        leakage_rows.append(
            {
                "candidate_program": label,
                "highest_broad_celltype": top_broad,
                "highest_median_score": float(by_broad.iloc[0]),
                "second_broad_celltype": by_broad.index[1] if len(by_broad) > 1 else top_broad,
                "second_median_score": float(by_broad.iloc[1]) if len(by_broad) > 1 else float(by_broad.iloc[0]),
                "score_margin": float(by_broad.iloc[0] - (by_broad.iloc[1] if len(by_broad) > 1 else 0.0)),
            }
        )
    leakage = pd.DataFrame(leakage_rows).sort_values("score_margin")
    return adata_annot, cluster_scores, subtype_table, leakage


def process_tissue(tissue: str) -> None:
    sample = ad.read_h5ad(SAMPLE_H5AD[tissue])
    annot = ad.read_h5ad(ANNOT_H5AD[tissue])

    doublet_table = run_scrublet_if_available(sample)
    if doublet_table is not None:
        annot.obs = annot.obs.merge(
            doublet_table.set_index("cell_barcode"),
            left_index=True,
            right_index=True,
            how="left",
        )
        if "scrublet_score" in annot.obs:
            annot.obs["scrublet_score"] = pd.to_numeric(
                annot.obs["scrublet_score"], errors="coerce"
            )
        if "predicted_doublet" in annot.obs:
            annot.obs["predicted_doublet"] = pd.Categorical(
                annot.obs["predicted_doublet"]
                .map({True: "True", False: "False"})
                .fillna("NA")
            )
        doublet_summary = (
            doublet_table["predicted_doublet"].value_counts(dropna=False)
            .rename_axis("predicted_doublet")
            .reset_index(name="n_cells")
        )
        doublet_summary.to_csv(TABLEDIR / f"RRRM1_{tissue}_doublet_summary.csv", index=False)

    annot, score_matrix, subtype_table, leakage = subtype_and_leakage(tissue, annot)
    score_matrix.to_csv(TABLEDIR / f"RRRM1_{tissue}_refined_subtype_score_matrix.csv", index=False)
    subtype_table.to_csv(TABLEDIR / f"RRRM1_{tissue}_refined_subtype_by_cluster.csv", index=False)
    leakage.to_csv(TABLEDIR / f"RRRM1_{tissue}_ambient_leakage_review.csv", index=False)

    if "refined_subtype" in annot.obs:
        sc.pl.umap(
            annot,
            color=["refined_subtype"],
            legend_loc="right margin",
            show=False,
            save=False,
        )
        plt.savefig(FIGDIR / f"RRRM1_{tissue}_umap_refined_subtype.png", dpi=200, bbox_inches="tight")
        plt.close()

    annot.write_h5ad(OBJDIR / f"RRRM1_{tissue}_hardened.h5ad")


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tissue",
        choices=["blood", "eye", "muscle", "skin"],
        default=None,
        help="Run hardening for a single tissue only",
    )
    args = parser.parse_args()

    ensure_dirs()
    tissues = [args.tissue] if args.tissue else ["blood", "eye", "muscle", "skin"]
    for tissue in tissues:
        process_tissue(tissue)
    print(f"Wrote hardening outputs under: {OUTDIR}")


if __name__ == "__main__":
    main()

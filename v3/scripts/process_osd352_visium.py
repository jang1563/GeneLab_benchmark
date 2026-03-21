#!/usr/bin/env python3
"""
process_osd352_visium.py — Load OSD-352 RR-3 Brain Visium data
  - SpaceRanger outputs → AnnData per section
  - QC filtering
  - Pseudo-bulk aggregation (section-level and animal-level)
  - Spatial coordinate loading for Moran's I

Data: Masarapu et al. 2024 Nature Communications
  - 3 FLT + 3 GC animals, 12 Visium sections (4 per slide, 3 slides)
  - Sample names: Sample_{slide}_{capture_area}
"""

import sys
import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy import sparse


# ── Paths ──────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent  # v3/
DATA_DIR = BASE_DIR / "data" / "spatial" / "osd352"
VISIUM_DIR = DATA_DIR / "visium_data" / "visium_data"  # double nesting from zip
PROCESSED_DIR = DATA_DIR / "processed"

# ── Section → animal/condition mapping ──────────────────────────────────────
# From Masarapu et al. 2024: 3 slides (158, 159, 304), 4 capture areas each
# Mapping must be determined from Seurat metadata or paper supplementary
# Slide 158 = animals from one group, slide 159 = another, etc.
# We'll load this from seurat_metadata.csv if available, otherwise from
# a manually curated mapping.

# Fallback mapping (to be verified with Seurat metadata)
# Based on typical Visium experiment design:
#   Each slide has 4 capture areas (A1, B1, C1, D1)
#   With 6 animals (3 FLT + 3 GC), 2 sections per animal = 12 sections
#   3 slides × 4 areas = 12 sections
# This needs verification - see load_section_metadata()

def load_section_metadata(data_dir: Path) -> pd.DataFrame:
    """Load section → animal → condition mapping.

    Priority:
    1. seurat_metadata.csv (extracted from Seurat RDS)
    2. Manual mapping from paper/supplementary
    """
    seurat_meta_path = data_dir / "seurat_metadata.csv"
    if seurat_meta_path.exists():
        meta = pd.read_csv(seurat_meta_path, index_col=0)
        # Extract unique section info from Seurat barcode metadata
        # Seurat barcodes: {barcode}-1_{section_id} or {barcode}_{section_id}
        # orig.ident typically contains sample/section identifier
        if "orig.ident" in meta.columns:
            sections = meta["orig.ident"].unique()
            print(f"  Found {len(sections)} sections in Seurat metadata: {sections}")

        # Look for condition/animal columns
        condition_cols = [c for c in meta.columns
                         if c.lower() in ("condition", "group", "type",
                                           "sample_type", "flight_status")]
        animal_cols = [c for c in meta.columns
                       if c.lower() in ("animal", "animal_id", "mouse_id",
                                         "subject", "replicate")]
        cluster_cols = [c for c in meta.columns
                        if "cluster" in c.lower() or "seurat_clusters" == c.lower()]

        info = {"seurat_meta": meta}
        if condition_cols:
            info["condition_col"] = condition_cols[0]
            print(f"  Condition column: {condition_cols[0]}")
            print(f"  Values: {meta[condition_cols[0]].unique()}")
        if animal_cols:
            info["animal_col"] = animal_cols[0]
        if cluster_cols:
            info["cluster_col"] = cluster_cols[0]
            print(f"  Cluster column: {cluster_cols[0]}")
            n_clusters = meta[cluster_cols[0]].nunique()
            print(f"  N clusters: {n_clusters}")
        return info

    # Fallback: manual mapping (must be filled in after Seurat extraction)
    print("  WARNING: No seurat_metadata.csv found. Using fallback mapping.")
    print("  Run extract_seurat_metadata.R first!")
    return None


def load_visium_section(section_dir: Path) -> ad.AnnData:
    """Load a single Visium section from SpaceRanger output.

    Expected files in section_dir:
    - filtered_feature_bc_matrix.h5
    - tissue_positions_list.csv
    - scalefactors_json.json
    - tissue_hires_image.png (optional for analysis)
    """
    h5_path = section_dir / "filtered_feature_bc_matrix.h5"
    positions_path = section_dir / "tissue_positions_list.csv"
    scalefactors_path = section_dir / "scalefactors_json.json"

    if not h5_path.exists():
        raise FileNotFoundError(f"No h5 file in {section_dir}")

    # Load expression matrix
    adata = sc.read_10x_h5(str(h5_path))
    adata.var_names_make_unique()

    # Load spatial coordinates
    if positions_path.exists():
        # Visium tissue_positions_list.csv format:
        # barcode, in_tissue, array_row, array_col, pxl_row_in_fullres, pxl_col_in_fullres
        pos = pd.read_csv(positions_path, header=None,
                          names=["barcode", "in_tissue", "array_row", "array_col",
                                 "pxl_row", "pxl_col"])
        pos = pos.set_index("barcode")

        # Filter to barcodes present in expression matrix
        common_barcodes = adata.obs_names.intersection(pos.index)
        adata = adata[common_barcodes].copy()
        pos = pos.loc[common_barcodes]

        # Store spatial coordinates
        adata.obsm["spatial"] = pos[["pxl_row", "pxl_col"]].values.astype(float)
        adata.obs["in_tissue"] = pos["in_tissue"].values
        adata.obs["array_row"] = pos["array_row"].values
        adata.obs["array_col"] = pos["array_col"].values

        # Filter to in-tissue spots only
        in_tissue_mask = adata.obs["in_tissue"] == 1
        adata = adata[in_tissue_mask].copy()
    else:
        print(f"  WARNING: No tissue_positions_list.csv in {section_dir}")

    # Load scale factors
    if scalefactors_path.exists():
        with open(scalefactors_path) as f:
            adata.uns["spatial_scalefactors"] = json.load(f)

    return adata


def qc_filter(adata: ad.AnnData, section_name: str,
              min_genes: int = 200, min_cells: int = 3,
              max_mito_pct: float = 20.0) -> ad.AnnData:
    """Standard QC filtering for Visium spots."""
    n_before = adata.n_obs

    # Calculate QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("mt-") | \
                      adata.var_names.str.startswith("Mt-") | \
                      adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Filter spots
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Filter high-mito spots
    mito_mask = adata.obs["pct_counts_mt"] < max_mito_pct
    adata = adata[mito_mask].copy()

    n_after = adata.n_obs
    print(f"  {section_name}: {n_before} → {n_after} spots "
          f"(removed {n_before - n_after}, "
          f"median genes={adata.obs['n_genes_by_counts'].median():.0f})")
    return adata


def normalize_section(adata: ad.AnnData) -> ad.AnnData:
    """Normalize counts: total-count normalize + log1p."""
    # Store raw counts
    adata.layers["counts"] = adata.X.copy()
    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def pseudobulk_section(adata: ad.AnnData) -> np.ndarray:
    """Compute pseudo-bulk for a section: mean expression across all spots."""
    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()
    return np.mean(X, axis=0).astype(np.float32)


def load_all_sections(visium_dir: Path, verbose: bool = True) -> dict:
    """Load, QC, and normalize all Visium sections.

    Returns:
        dict: section_name → AnnData
    """
    section_dirs = sorted([d for d in visium_dir.iterdir()
                           if d.is_dir() and d.name.startswith("Sample_")])

    if not section_dirs:
        raise FileNotFoundError(f"No Sample_* directories found in {visium_dir}")

    sections = {}
    for sd in section_dirs:
        name = sd.name
        if verbose:
            print(f"Loading {name}...")
        try:
            adata = load_visium_section(sd)
            adata = qc_filter(adata, name)
            adata = normalize_section(adata)
            adata.obs["section"] = name
            sections[name] = adata
        except Exception as e:
            print(f"  ERROR loading {name}: {e}")
            continue

    print(f"\nLoaded {len(sections)} sections:")
    for name, adata in sections.items():
        print(f"  {name}: {adata.n_obs} spots × {adata.n_vars} genes")

    return sections


def build_pseudobulk_matrix(sections: dict, gene_names: np.ndarray = None
                            ) -> tuple:
    """Build pseudo-bulk matrix: sections × genes.

    Args:
        sections: dict of section_name → AnnData
        gene_names: common gene set (if None, use intersection)

    Returns:
        X: np.ndarray (n_sections × n_genes)
        section_names: list of section names
        gene_names: np.ndarray of gene names
    """
    if gene_names is None:
        # Find common genes across all sections
        gene_sets = [set(adata.var_names) for adata in sections.values()]
        common_genes = sorted(set.intersection(*gene_sets))
        gene_names = np.array(common_genes)
        print(f"Common genes across {len(sections)} sections: {len(gene_names)}")

    section_names = sorted(sections.keys())
    X = np.zeros((len(section_names), len(gene_names)), dtype=np.float32)

    for i, name in enumerate(section_names):
        adata = sections[name]
        # Subset to common genes
        adata_sub = adata[:, gene_names]
        X[i] = pseudobulk_section(adata_sub)

    return X, section_names, gene_names


def main():
    parser = argparse.ArgumentParser(
        description="Process OSD-352 RR-3 Brain Visium data")
    parser.add_argument("--data-dir", type=str,
                        default=str(DATA_DIR),
                        help="Path to OSD-352 data directory")
    parser.add_argument("--check-only", action="store_true",
                        help="Only check data structure, don't process")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    visium_dir = data_dir / "visium_data" / "visium_data"

    print("=" * 60)
    print("OSD-352 RR-3 Brain Visium — Data Processing")
    print("=" * 60)

    # Step 1: Check metadata
    print("\n[1] Loading section metadata...")
    meta_info = load_section_metadata(data_dir)

    # Step 2: List sections
    print("\n[2] Available sections...")
    section_dirs = sorted([d for d in visium_dir.iterdir()
                           if d.is_dir() and d.name.startswith("Sample_")])
    for sd in section_dirs:
        files = list(sd.iterdir())
        print(f"  {sd.name}: {len(files)} files — "
              f"{[f.name for f in files if not f.name.startswith('.')]}")

    if args.check_only:
        print("\n[check-only] Done.")
        return

    # Step 3: Load and process all sections
    print("\n[3] Loading sections...")
    sections = load_all_sections(visium_dir)

    # Step 4: Build pseudo-bulk matrix
    print("\n[4] Building pseudo-bulk matrix...")
    X, section_names, gene_names = build_pseudobulk_matrix(sections)
    print(f"  Pseudo-bulk matrix: {X.shape} (sections × genes)")

    # Step 5: Save processed data
    processed_dir = data_dir / "processed"
    processed_dir.mkdir(parents=True, exist_ok=True)

    np.save(processed_dir / "pseudobulk_matrix.npy", X)
    np.save(processed_dir / "gene_names.npy", gene_names)
    with open(processed_dir / "section_names.json", "w") as f:
        json.dump(section_names, f)

    # Save per-section spot counts for reference
    spot_info = {name: {"n_spots": adata.n_obs, "n_genes": adata.n_vars}
                 for name, adata in sections.items()}
    with open(processed_dir / "section_info.json", "w") as f:
        json.dump(spot_info, f, indent=2)

    # Save merged AnnData (all sections concatenated) for spatial analysis
    print("\n[5] Saving merged AnnData...")
    adatas = []
    for name in section_names:
        a = sections[name]
        a.obs["section"] = name
        adatas.append(a)
    merged = ad.concat(adatas, label="section", keys=section_names)
    merged.write_h5ad(processed_dir / "osd352_visium_merged.h5ad")
    print(f"  Saved: {merged.n_obs} total spots × {merged.n_vars} genes")

    print("\n[DONE] Processing complete.")
    print(f"  Outputs in: {processed_dir}")


if __name__ == "__main__":
    main()

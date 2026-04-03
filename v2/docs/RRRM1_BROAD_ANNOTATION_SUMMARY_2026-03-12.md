# RRRM-1 Broad Annotation Summary

Date: 2026-03-12
Scope: first-pass broad cell type annotation for RRRM-1 scRNA-seq downstream outputs

## Output location

- Base directory:
  - `$SCRATCH_DIR/rrrm1_scrna/downstream_initial/annotations`

## Generated files

- UMAP figures:
  - `RRRM1_blood_umap_broad_celltype.png`
  - `RRRM1_eye_umap_broad_celltype.png`
  - `RRRM1_muscle_umap_broad_celltype.png`
  - `RRRM1_skin_umap_broad_celltype.png`
- Annotated objects:
  - `RRRM1_blood_annotated.h5ad`
  - `RRRM1_eye_annotated.h5ad`
  - `RRRM1_muscle_annotated.h5ad`
  - `RRRM1_skin_annotated.h5ad`
- Tables per tissue:
  - `*_broad_celltype_counts.csv`
  - `*_cluster_annotation.csv`
  - `*_cluster_score_matrix.csv`
- Cross-tissue summary:
  - `summary/RRRM1_broad_celltype_counts_all_tissues.csv`
  - `summary/RRRM1_broad_celltype_proportions_all_tissues.csv`
  - `figures/RRRM1_broad_celltype_stacked_bar.png`

## First-pass cell type composition

### Blood (`OSD-918`)

- erythroid: `3652`
- b_cell: `290`
- t_cell: `237`
- neutrophil: `177`
- monocyte_macrophage: `25`
- nk_cell: `14`

Interpretation:
- Blood is strongly erythroid-heavy.
- Proportions:
  - erythroid `83.1%`
  - B cell `6.6%`
  - T cell `5.4%`

### Eye (`OSD-920`)

- retinal_neuronal: `1136`
- muller_glia: `405`
- corneal_conjunctival_epithelial: `225`
- stromal_fibroblast: `134`
- photoreceptor_like: `119`
- immune: `105`
- endothelial_perivascular: `56`
- lens_crystallin: `26`

Interpretation:
- Eye contains a mixed retinal/glial signature with smaller epithelial, stromal, and immune compartments.
- Proportions:
  - retinal neuronal `51.5%`
  - Muller glia `18.4%`
  - corneal/conjunctival epithelial `10.2%`

### Muscle (`OSD-924`)

- macrophage_myeloid: `6125`
- t_nk_lymphocyte: `3950`
- fibroadipogenic_progenitor: `2121`
- endothelial: `2035`
- satellite_myogenic_progenitor: `778`
- fibroblast_matrix: `516`
- pericyte_smooth_muscle: `263`
- myonuclear_structural: `230`

Interpretation:
- Muscle is immune-rich with stromal and endothelial support populations.
- Proportions:
  - macrophage/myeloid `38.2%`
  - T/NK `24.7%`
  - fibroadipogenic progenitor `13.2%`
  - endothelial `12.7%`

### Skin (`OSD-934`)

- immune_myeloid: `5881`
- basal_keratinocyte: `4745`
- suprabasal_keratinocyte: `4186`
- hair_follicle_epithelial: `790`
- t_nk_lymphocyte: `256`
- melanocyte_like: `18`

Interpretation:
- Skin shows expected keratinocyte structure with a large immune myeloid component.
- Proportions:
  - immune myeloid `37.0%`
  - basal keratinocyte `29.9%`
  - suprabasal keratinocyte `26.4%`

## Notes

- This is a first-pass marker-score annotation, not a final manuscript-level label set.
- Cluster labels should still be checked against `*_cluster_markers_top20.csv` and `*_cluster_annotation.csv`.
- The current output is sufficient to start figure drafting and tissue-specific deeper annotation.
- The stacked bar plot provides a quick cross-tissue composition panel for reports or slides.

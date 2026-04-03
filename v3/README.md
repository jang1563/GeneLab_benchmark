# GeneLabBench v3

**Multi-species, spatial, scRNA-seq, and foundation model evaluation**
**Date**: 2026-03-18 ~ 2026-03-20
**Status**: Complete

## Version History

| Version | Focus | Date | Status |
|---------|-------|------|--------|
| v1.0 | Bulk RNA-seq benchmark (6 tissues, LOMO, 25+ tasks) | 2025-12 ~ 2026-02 | Complete |
| v1.3 | +Temporal (T1-T3) +Foundation Models (scGPT, Geneformer) +LLM | 2026-02 ~ 2026-03-09 | Complete |
| v2.0 | +Cross-species (E1-E3) +PBMC (F1) +RRRM-1 scRNA-seq (F2) | 2026-03-07 ~ 2026-03-18 | Complete |
| **v3.0** | +Multi-species (E4) +Spatial (F3) +RRRM-2 (F5) +FM (UCE, scFoundation) | 2026-03-18 ~ 2026-03-20 | **Complete** |

## Directory Structure

```
v3/
  docs/           # Plans and runbooks (3 files)
  scripts/        # Analysis scripts: 20 Python + 8 Bash + 2 R (30 files)
  evaluation/     # Output JSONs with metrics (19 files)
  figures/        # D3.js interactive HTML figures (5 files)
  processed/      # Intermediate data (not tracked in git)
  logs/           # SLURM job logs (not tracked in git)
```

## Key Results

### Phase 1 — E4 Multi-species (Mouse + Drosophila)
- Mouse intra-species KEGG NES correlations: significant tissue pairs (eye-liver r=+0.32, kidney-skin r=-0.27)
- Drosophila cross-species: negative correlations only (skin r=-0.59, thymus r=-0.19)
- Arabidopsis: no cross-kingdom orthologs available

### Phase 2 — F3 Spatial Visium (Brain, OSD-352 RR-3)
- **Negative result**: Section-level AUROC=0.139, animal-level AUROC=0.444
- Bulk RNA-seq AUROC=0.000 (genuine overfitting at n=6, not a bug)
- PC1 (42.5%) = slide batch effect, not spaceflight condition
- Brain tissue has no detectable spaceflight signal

### Phase 3 — F5 RRRM-2 scRNA-seq (4 tissues)
- **PBMC**: NK cell 0.845***, T cell 0.752***, CD4 0.651***, neutrophil 0.629***, monocyte 0.575***
- **Spleen**: B cell 0.562***, macrophage 0.521***, dendritic 0.520*
- **Bone marrow**: All 14 cell types near chance (0.27-0.54) — no spaceflight signal
- **Cross-mission** (RRRM-1 vs RRRM-2): low concordance (r=+0.011)

### Phase 4 — Extended Tissues + Radiation
- A7 Lung: 0.608, A7b Colon: 0.900, A8 Skin: 0.669
- Radiation→ISS concordance: r=+0.14 ~ +0.44; HLU→ISS: r=-0.20 ~ -0.50

### Phase 5 — B_ext Cross-tissue Transfer
- 7x7 transfer matrix (42 pairs), Method A range 0.35-0.80, Method C range 0.43-0.87

### Foundation Model Comparison (UCE + scFoundation vs PCA-LR)
- **UCE** (seeded): thymus 0.632 (p=0.031) only significant; all others near chance
- **scFoundation**: liver 0.635 (p=0.001), gastro 0.691 (p=0.035) significant
- **Verdict**: Both FMs underperform PCA-LR baseline (liver 0.670, 6-tissue mean 0.758)
- Pre-trained cell atlas knowledge does not improve spaceflight detection

## Figures

| Figure | Content |
|--------|---------|
| Fig1 | Multi-species NES concordance heatmaps + phylogenetic tree |
| Fig2 | Spatial Visium overview (AUROC, SVGs, quality, Venn) |
| Fig3 | RRRM-2 scRNA-seq AUROC + fGSEA heatmaps |
| Fig4 | Extended tissues + radiation + cross-tissue transfer matrix |
| Fig5 | Foundation model comparison (5 models x 7 tissues) |

## Evaluation JSONs

| File | Phase | Description |
|------|-------|-------------|
| `E4_multispecies_nes.json` | E4 | Mouse + Drosophila NES correlations |
| `E5_phylogenetic.json` | E5 | Phylogenetic distance metrics |
| `F3a_visium_classification.json` | F3 | Section/animal AUROC |
| `F3b_visium_svg.json` | F3 | Spatially variable genes |
| `F3d_visium_cross_resolution.json` | F3 | Bulk vs spatial comparison |
| `F5A_rrrm2_composition.json` | F5 | Cell composition per animal |
| `F5B_rrrm2_fgsea.json` | F5 | Pathway enrichment per cell type |
| `F5C_rrrm2_loao.json` | F5 | Cell-type AUROC (LOAO-CV) |
| `F5D_rrrm2_cross_mission.json` | F5 | RRRM-1 vs RRRM-2 concordance |
| `F5E_rrrm2_bone_marrow.json` | F5 | Bone marrow benchmark |
| `A7_lung.json` | A7 | Lung tissue classification |
| `A7b_colon.json` | A7b | Colon tissue classification |
| `A8_extended_skin.json` | A8 | Skin tissue classification |
| `R1_radiation_classification.json` | R1 | Radiation classification |
| `R2_radiation_dge_concordance.json` | R2 | DGE concordance |
| `R3_combined_vs_single.json` | R3 | Combined vs single model |
| `B_ext_transfer_matrix.json` | B_ext | 7x7 cross-tissue transfer |
| `FM_uce.json` | FM | UCE 7-tissue evaluation |
| `FM_scfoundation.json` | FM | scFoundation 7-tissue evaluation |

## HPC Environment

- **Cluster**: Cornell Cayuga (Slurm)
- **Work directory**: `$GENELAB_ROOT/`
- **Conda**: `~/miniconda3/miniconda3/`
- **Environments**: `perturb_seq_new` (general), `uce_env` (UCE), `scfoundation` (scFoundation)
- **Partition**: `scu-gpu` (A40 48GB / A100 80GB)

## Documentation

- `docs/PLAN_V3.md` — Master execution plan with verification checkpoints
- `docs/DATA_CATALOG_V3.md` — Dataset inventory (accessions, sample sizes, tissues)
- `docs/HPC_FM_RUNBOOK.md` — Foundation model setup and execution guide

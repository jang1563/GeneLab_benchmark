# GeneLab Benchmark v2.0

**Status**: Complete (2026-03-18)

## Scope

v2.0 extends the v1.0 mouse bulk RNA-seq benchmark with temporal dynamics, cross-species conservation, and single-cell resolution:

| Category | Description | Data | Status |
|----------|-------------|------|--------|
| **T** | Temporal Dynamics | ISS-T/LAR timing, recovery, age×spaceflight | **Complete** |
| **E** | Cross-Species Conservation | Mouse bulk ↔ JAXA cfRNA (human), duration effects | **Complete** |
| **F1** | Human PBMC Single-Cell | I4 snRNA-seq cell-type fGSEA, cfRNA origin | **Complete** |
| **F2** | Mouse Multi-Tissue scRNA-seq | RRRM-1 (OSD-918/920/924/934), 4 tissues | **Complete** |

---

## Completed Analysis (Phases 1–5)

### T1: ISS-T vs LAR Preservation Artifact
- GC AUROC ≥ FLT AUROC → preservation method dominates, not spaceflight biology
- Cross-mission transfer confirms consistent artifact (RR-6↔RR-8 gene AUROC 0.84–0.99)

### T2: Post-Flight Recovery Signature
- PCA recovery ratio: RR-8 = 0.652 (stronger), RR-6 = 0.842 (partial)
- 25/27 Hallmark pathways show recovery in RR-8, with MYC/biosynthesis overshoot

### T3: Age × Spaceflight Interaction
- "Spaceflight amplifies aging" SUPPORTED: OLD AUROC=0.945 vs YNG=0.679 (Δ=+0.266)

### E1: Cross-Species NES Conservation
- Mouse liver bulk vs JAXA cfRNA Hallmark NES: Spearman r=0.352

### E2: Duration Effect
- Short-duration I4 PBMC r=−0.095 vs long-duration JAXA r=+0.352 (Δr=+0.446)

### E3: cfRNA Cellular Origin
- Anti-correlation between cfRNA and expected tissue-of-origin signatures

### F1: I4 PBMC Cell-Type Pathway Analysis
- 10 cell types × 50 Hallmark pathways fGSEA from GLDS-562 snRNA-seq
- Temporal analysis: delayed PBMC response mirrors mouse T2 recovery overshoot

### Publication Figures
- Fig1_temporal.html: T1/T2/T3 (3 panels)
- Fig2_crossspecies.html: E1/E2/E3 (3 panels)
- Fig3_pbmc_celltype.html: F1 NES heatmap + F1 temporal (2 panels)
- All: D3.js v7, Okabe-Ito palette, self-contained HTML+SVG

---

## RRRM-1 scRNA-seq (F2, In Progress)

### Pipeline Status

| Step | Script | Status |
|------|--------|--------|
| OSDR download | `rrrm1_osdr_download.sh` | Complete |
| STARsolo alignment | `rrrm1_starsolo_job.sh` | Complete |
| H5AD conversion | `rrrm1_h5ad_convert.py` | Complete |
| Per-SRX merge | `rrrm1_merge_per_srx.py` | Complete |
| Tissue merge | `rrrm1_merge_h5ad.py` | Complete |
| QC + processing | `rrrm1_initial_scanpy.py` | Complete |
| Broad annotation | `rrrm1_broad_annotate.py` | Complete |
| Doublet removal | `rrrm1_singlecell_hardening.py` | Complete |
| **F2 benchmark tasks** | `rrrm1_benchmark.py` | **Complete** |

### Cell Counts (post-hardening)

| Tissue (OSD) | Total Cells | Major Cell Types |
|--------------|-------------|------------------|
| Blood (918) | 4,377 | erythroid 83%, B cell 7%, T cell 5% |
| Eye (920) | 2,197 | retinal neuronal 52%, Muller glia 18% |
| Muscle (924) | 15,669 | macrophage/myeloid 38%, T/NK 25%, FAP 13% |
| Skin (934) | 15,838 | immune myeloid 37%, basal keratinocyte 30% |

### Planned Benchmark Tasks (F2)

| Task | Question | Method |
|------|----------|--------|
| **F2-A** | Cell-type composition shift (FLT vs GC) | Wilcoxon per cell type, LOAO |
| **F2-B** | Cell-type-specific pathway activity | Pseudo-bulk DESeq2 + fGSEA |
| **F2-C** | Cell-level spaceflight classifier | PCA-LR per cell type, LOAO |
| **F2-D** | Cross-species concordance (RRRM-1 ↔ I4 PBMC) | NES Spearman r |

**Blocker**: FLT/GC condition labels need to be mapped from OSDR metadata → h5ad obs. SRX→condition mapping available in `docs/RRRM1_SRX_CONDITION_MAP.csv`.

---

## Directory Structure

```
v2/
├── README.md                       ← This file
├── scripts/                        ← v2 analysis scripts (19 Python)
│   ├── temporal_analysis.py        ← T1/T2/T3 temporal dynamics
│   ├── cross_species_nes_comparison.py  ← E1/E2 cross-species NES
│   ├── e3_cfrna_celltype_origin.py      ← E3 cfRNA origin analysis
│   ├── mission_conservation_comparison.py ← Cross-mission NES
│   ├── generate_main_figures.py    ← v2 integrated figures (D3.js v7)
│   ├── i4_celltype_pathway_heatmap.py   ← F1 NES heatmap
│   ├── i4_temporal_heatmap.py      ← F1 temporal analysis
│   ├── llm_generate_prompts.py     ← Tier 3 prompt generation
│   ├── llm_call_api.py             ← Tier 3 API calls
│   ├── llm_parse_responses.py      ← Tier 3 response parsing
│   ├── pipeline_version_compare.py ← J1 pipeline comparison
│   ├── annotate_radiation.py       ← T4 radiation metadata
│   ├── rrrm1_h5ad_convert.py       ← RRRM-1 STARsolo → h5ad
│   ├── rrrm1_merge_per_srx.py      ← RRRM-1 per-SRX merge
│   ├── rrrm1_merge_h5ad.py         ← RRRM-1 tissue merge
│   ├── rrrm1_initial_scanpy.py     ← RRRM-1 QC + processing
│   ├── rrrm1_broad_annotate.py     ← RRRM-1 marker-based annotation
│   ├── rrrm1_singlecell_hardening.py ← RRRM-1 doublet removal
│   └── rrrm1_annotation_summary.py ← RRRM-1 cross-tissue summary
├── evaluation/                     ← v2 results and metrics
│   └── V2_RESULTS_SUMMARY.md      ← Comprehensive v2 results
├── docs/                           ← v2 design documents
│   ├── V2_PAPER_CONTENT.md         ← v2 manuscript draft
│   ├── DATA_CATALOG_V2.md          ← v2 data sources (JAXA, I4, RRRM-1)
│   ├── V2_SUPPLEMENT_TABLES.md     ← Supplementary tables
│   ├── RRRM1_SC_BENCHMARK_PLAN_V3.md ← F2 benchmark task definitions
│   ├── RRRM1_SRX_CONDITION_MAP.csv ← SRX → FLT/GC mapping (32 samples)
│   ├── RRRM1_BROAD_ANNOTATION_SUMMARY_2026-03-12.md
│   └── RRRM1_TISSUE_MARKERS_2026-03-12.md
├── data/                           ← v2 additional data
│   └── jaxa_cfrna/                 ← JAXA cfRNA NES data
├── processed/                      ← v2 intermediate outputs
│   ├── T_temporal/                 ← T1/T2/T3 JSON results
│   ├── T4_radiation/               ← Radiation metadata
│   ├── E_crossspecies/             ← E1/E2/E3 results
│   ├── F1_scrna/                   ← I4 PBMC cell-type fGSEA
│   └── llm_responses/             ← Tier 3 archived responses
└── figures/                        ← v2 integrated HTML figures
    ├── Fig1_temporal.html
    ├── Fig2_crossspecies.html
    └── Fig3_pbmc_celltype.html
```

## Relationship to v1

- v1 code lives in the project root (`scripts/`, `evaluation/`, etc.)
- Shared utilities from `scripts/utils.py` can be imported by v2 scripts
- v2 scripts should NOT modify v1 files
- v2 results complement v1 findings (temporal context, cell-type resolution, cross-species validation)

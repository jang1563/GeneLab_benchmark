# GeneLab Benchmark v2 — Data Catalog

*Created: 2026-03-09. Tracks all datasets for v2 analyses (E/F/G categories).*

---

## Already Available (No Download Required)

### Mouse Temporal Analysis (v2 Section 1: T1/T2/T3)

| Dataset | Location | Description | Analysis |
|---------|----------|-------------|---------|
| Mouse liver fGSEA (6 missions) | `processed/fgsea/liver/*.csv` | Per-mission Hallmark/KEGG/Reactome NES values | T1/T2/E1 |
| RR-8 liver temporal groups | `data/mouse/liver/` | ISS-T and LAR samples (BSL/FLT/GC × 3 groups) | T1/T2 |
| RR-8 liver aging cohort | `data/mouse/liver/` | OLD (32 week) + YNG (10 to 12 week) mice | T3 |
| T1/T2/T3 results | `v2/processed/T_temporal/*.json` | Computed results (recovery ratios, AUROCs, etc.) | Done |

---

## Category E: Cross-Species Conservation

### E1: Mouse Liver × Human cfRNA (JAXA CFE)

**GO decision: YES — all data local, E1 immediately executable**

| Dataset | Path | Format | Notes |
|---------|------|--------|-------|
| Human cfRNA DE (3-group) | `/SpaceOmicsBench/v2_public/data/processed/cfrna_3group_de_noleak.csv` | 26,845 genes × 30 cols | JAXA CFE OSD-530, pre/in-flight/post groups |
| Human cfRNA timecourse (11-group) | `/SpaceOmicsBench/v2_public/data/processed/cfrna_11group_timecourse.csv` | 26,845 genes × 23 cols | Pre1-3, Flight1-4, Post1-4 normalized means |
| Mouse liver fGSEA Hallmark NES | `processed/fgsea/liver/{mission}_fgsea_hallmark.csv` | 6 missions × ~50 pathways | Per-pathway NES + padj |
| Mouse-Human ortholog table | `data/mouse/ensembl_mouse_human_orthologs_with_symbols.tsv` | 25,166 pairs | ENSMUSG ↔ ENSG + human gene symbol |

**Data structure:**
- `cfrna_3group_de_noleak.csv`: columns include `gene` (human symbol), `edge_pre_vs_flight_diff` (mean_flight − mean_pre), `edge_pre_vs_flight_fc` (fold change)
- `cfrna_11group_timecourse.csv`: per-timepoint normalized means; groups: pre1-3, flight1-4, post1-4
- Mouse fGSEA CSVs: columns `pathway, pval, padj, NES, size, tissue, mission`

**Analysis approach (E1):**
1. Rank human cfRNA genes by `edge_pre_vs_flight_diff` (in-flight vs. pre-flight)
2. Run fGSEA (Python `gseapy`) against MSigDB Hallmark gene sets (Homo sapiens)
3. Extract NES per pathway → human cfRNA NES vector (50 pathways)
4. Compute mission-averaged mouse liver NES vector from 6 missions
5. Spearman r (human NES vs. mouse NES), 95% CI by bootstrap
6. Also: per-mission comparison (mouse RR-8 NES vs. human NES, RR-6 vs. human NES)

**No ortholog mapping needed** at the Hallmark pathway level — Hallmark gene set names are consistent across species.

**Script:** `v2/scripts/cross_species_nes_comparison.py`
**Output:** `v2/evaluation/E1_crossspecies_nes.json` + scatter plot HTML

### E2: Cross-Mission Stability (Short vs. Long Duration)

| Comparison | Mouse data | Human data | Metric |
|-----------|-----------|-----------|--------|
| RR-6 (35d) vs. human | `processed/fgsea/liver/RR-6_fgsea_hallmark.csv` | JAXA CFE cfRNA NES | Spearman r |
| RR-8 (35d) vs. human | `processed/fgsea/liver/RR-8_fgsea_hallmark.csv` | JAXA CFE cfRNA NES | Spearman r |
| All missions vs. human | Mission-averaged NES | JAXA CFE cfRNA NES | Spearman r |

---

## Category F: Mouse Single-Cell Analysis

### F1/F2: RRRM-1/RR-8 Multi-Tissue scRNA-seq

**GO decision: YES — RRRM-1 benchmark-aligned subset recovered on Cayuga**

| Dataset | OSDR ID | Expected size | Status |
|---------|---------|--------------|--------|
| RRRM-1 multi-tissue scRNA | Selected subset: OSD-918, OSD-920, OSD-924, OSD-934 | ~308 GB raw FASTQ for selected subset | Downloaded, processed, and benchmark-ready |
| RR-8 bone marrow scRNA | TBD from OSDR catalog | ~5-10 GB | Not downloaded |

**Current RRRM-1 state:** per-sample `Solo.out`, `.h5ad`, merged object, tissue-level processed objects, and broad cell-type annotations are available on Cayuga scratch. See `RRRM1_BENCHMARK_READY_MANIFEST_2026-03-12.csv`.

**HPC pipeline (Cayuga):**
```
1. Check OSDR API for processed data first
2. Download to /athena/masonlab/scratch/users/jak4013/GeneLab_v2/
3. Scanpy QC: min_genes=200, min_cells=3, mt_pct<20%
4. Leiden clustering + marker gene annotation
5. Per-cell-type FLT vs. GC: AUROC (PCA-LR, RSKF)
6. Composition analysis (F2): fraction FLT vs. GC per cell type
```

**Current script set:** `rrrm1_starsolo_job.sh`, `rrrm1_h5ad_convert.py`, `rrrm1_merge_h5ad.py`, `rrrm1_initial_scanpy.py`, `rrrm1_broad_annotate.py`, `rrrm1_annotation_summary.py`

**Expected cell types (based on tissue):**
- Spleen: T cells, B cells, NK cells, macrophages, dendritic cells
- Bone marrow: HSC, progenitors, erythroid, myeloid
- PBMC: CD4 T, CD8 T, B, NK, monocytes
- Liver (if available): hepatocytes, Kupffer cells, stellate cells

---

## Category G: Microbiome (Conditional)

**GO decision: TBD — need to verify mouse gut microbiome availability in OSDR**

| Data source | Location | Status | Notes |
|------------|---------|--------|-------|
| Human cfRNA microbiome | `/SpaceOmicsBench/v2_public/data/processed/microbiome_human_*.csv` | EXISTS | OSD-530 human |
| Human JAXA microbiome | `/SpaceOmicsBench/v2_public/data/processed/microbiome_*.csv` | EXISTS | Gut taxonomy + pathways |
| Mouse gut microbiome | OSDR GLDS-?? | NOT CHECKED | Need to search |

**SpaceOmicsBench microbiome files found:**
- `microbiome_gut_taxonomy_cpm.csv`
- `microbiome_gut_pathways_cpm.csv`
- `microbiome_human_taxonomy_cpm.csv`
- `microbiome_human_pathways_cpm.csv`
- `microbiome_env_taxonomy_cpm.csv`
- `microbiome_env_pathways_cpm.csv`

**Novel angle:** Human-mouse cross-species microbiome comparison during spaceflight.
**Blocker:** Mouse gut microbiome OSDR data needs verification.
**Decision:** Proceed only if mouse gut microbiome data available for ≥3 missions.

---

## Relationship to SpaceOmicsBench

SpaceOmicsBench (v2_public, v3, v3_dev) contains processed human astronaut data:
- **OSD-530** (JAXA CFE cfRNA): 6 astronauts, 11 timepoints = primary for E1
- Metabolomics, proteomics, microbiome, clinical data = potential for G category
- v3 includes AX-2 data in `axiom2_*/` directories

GeneLab_benchmark v2 adds:
- Mouse multi-tissue bulk RNA-seq × 6 missions (already in v1)
- Cross-species NES correlation (E1 = new)
- Mouse single-cell (F1 = new)

**Key SpaceOmicsBench file for E1:**
```
/Users/jak4013/Dropbox/Bioinformatics/Claude/SpaceOmicsBench/v2_public/data/processed/cfrna_3group_de_noleak.csv
```
- 26,845 human genes (gene symbols, e.g., MT-ND1, ACTB)
- `edge_pre_vs_flight_diff`: mean_flight_normalized − mean_pre_normalized (use as fGSEA rank)
- `edge_pre_vs_flight_fc`: fold change (alternative rank)
- No per-sample raw counts in CSV; these are group-level summary statistics

---

## Mouse fGSEA Hallmark NES — Quick Reference

Mouse liver NES inter-mission Spearman r (from `evaluation/NES_conservation_hallmark.json`):

| Mission pair | r |
|-------------|---|
| RR-6 ↔ RR-8 | 0.422 |
| RR-1 ↔ RR-8 | 0.271 |
| RR-1 ↔ RR-3 | 0.265 |
| MHU-2 ↔ RR-6 | 0.452 (highest) |
| Mean (all pairs) | ~0.15 |

→ Moderate within-mouse cross-mission NES conservation. E1 will compare mission-averaged mouse NES vs. human cfRNA NES.

---

## File Size and HPC Decision Matrix

| Dataset | Size | Location | HPC needed? |
|---------|------|---------|------------|
| cfRNA 3-group DE | 2 MB | Local (SpaceOmicsBench) | No |
| Mouse liver fGSEA CSVs | <1 MB each | Local (GeneLab) | No |
| Ortholog table | 5 MB | Local (GeneLab) | No |
| RRRM-1 scRNA counts | >10 GB | OSDR (not downloaded) | YES |
| RR-8 scRNA counts | ~5-10 GB | OSDR (not downloaded) | YES |

**Phase 2 (E1) can run entirely on MacBook** (< 5 min).
**Phase 3 (F1) requires Cayuga HPC** (large files + GPU for clustering at scale).

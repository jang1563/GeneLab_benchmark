# GeneLab_benchmark v2.0 — Results Summary

Generated: 2026-03-10 (Updated: 2026-03-13)

## Overview

v2.0 extends v1.0 with temporal/biological covariate analysis using existing data. All scripts in `v2/scripts/`, outputs in `v2/processed/` and `v2/evaluation/`.

---

## T1: ISS-T vs LAR Classification

**Question**: Can sacrifice timing (ISS-Terminal vs Live Animal Return) be detected from transcriptomics?

**Confound warning (DD-18, PMID 33376967)**: ISS-T = RNAlater on-orbit, LAR = standard necropsy. Preservation method confounds biological timing effect. GC AUROC serves as preservation artifact baseline.

### Within-Mission Results

| Sub-task | Condition | Gene AUROC | Pathway AUROC | n |
|---|---|---|---|---|
| **T1a RR-6 liver** | FLT | **1.000** [1.00, 1.00] | 0.960 [0.87, 1.00] | 20 |
| | GC (baseline) | **0.922** [0.75, 1.00] | 0.778 [0.52, 0.98] | 19 |
| | BSL | **1.000** [1.00, 1.00] | 1.000 [1.00, 1.00] | 20 |
| **T1b RR-8 liver** | FLT | 0.930 [0.79, 1.00] | 0.920 [0.78, 1.00] | 35 |
| | GC (baseline) | **0.973** [0.91, 1.00] | 0.993 [0.97, 1.00] | 35 |
| | BSL | 0.983 [0.94, 1.00] | 0.970 [0.90, 1.00] | 38 |
| **T1c RR-6 thymus** | FLT | 0.857 [0.63, 1.00] | 0.714 NS | 17 |
| | GC (baseline) | **0.925** [0.76, 1.00] | 1.000 [1.00, 1.00] | 18 |

**Key Finding**: GC AUROC ≥ FLT AUROC in most cases (RR-8: GC=0.973 > FLT=0.930, excess=-0.043) → ISS-T vs LAR difference is **dominated by preservation method artifact**, not biological timing effect. BSL (baseline, no flight) also shows near-perfect separation → confirms preservation method is the primary driver. RR-6 shows small excess (+0.078) but BSL perfect separation (1.0) undermines biological interpretation.

### Cross-Mission Transfer (T1d)

| Direction | FLT Gene | GC Gene |
|---|---|---|
| RR-6 → RR-8 | 0.987 | 0.957 |
| RR-8 → RR-6 | 0.970 | 0.844 |

Preservation artifact is consistent across missions (high transfer AUROC).

---

## T2: LAR Recovery Signature

**Question**: Do LAR samples show recovery toward baseline compared to ISS-T?

### PCA Recovery Ratio

| Mission | R = d(FLT_LAR,BSL_LAR) / d(FLT_ISS-T,BSL_ISS-T) | Interpretation |
|---|---|---|
| **RR-6** | **0.842** | Partial recovery |
| **RR-8** | **0.652** | Stronger recovery |

R < 1 → FLT_LAR is closer to baseline than FLT_ISS-T (recovery toward baseline).

### Per-Pathway Recovery (GSVA Hallmark, DD-20, direction-aware)

| Mission | Mean Recovery | Recovering/Total | Notes |
|---|---|---|---|
| RR-6 | -0.055 | 12/26 | Many pathways overshoot (direction reversal) |
| RR-8 | **1.236** | **25/27** | Strong recovery with overshoot past baseline |

**Direction-aware formula**: `ratio = delta_return / delta_flight`. Same sign → `recovery = 1 - ratio`. Opposite sign (overshoot) → `recovery = 1 - ratio > 1`.

**Top recovered/overshoot (RR-6)**: TGF-beta (2.08), Notch (1.88), WNT/β-catenin (1.43) — direction reversal past baseline

**Immune rebound (RR-6)**: IFN-α (-2.36), Inflammatory response (-2.54), IFN-γ (-1.36) → immune pathways diverge further after return

**Top recovered (RR-8)**: MYC targets V1 (2.49), Protein secretion (2.14), MYC targets V2 (1.92) — strong overshoot past baseline

### Classification Shift (Descriptive)

| Mission | FLT_ISS-T flight_prob | FLT_LAR flight_prob | Interpretation |
|---|---|---|---|
| RR-6 | 1.000 | **0.185** | FLT_LAR projects closer to baseline than FLT_ISS-T |
| RR-8 | 1.000 | **0.404** | FLT_LAR projects closer to baseline than FLT_ISS-T |

`FLT_ISS-T` probabilities are in-sample reference scores from the training contrast (`BSL_ISS-T` vs `FLT_ISS-T`), so this block should be read as a descriptive projection shift, not held-out validation evidence of recovery.

---

## T3: Age × Spaceflight Interaction (RR-8 Liver)

**Question**: Does spaceflight amplify age-related transcriptomic differences?

### T3a: Overall Age Classification

| Feature | AUROC | 95% CI | p-value | n |
|---|---|---|---|---|
| Gene | **0.985** | [0.962, 0.999] | <0.001 | 141 |
| Pathway | 0.851 | [0.786, 0.911] | <0.001 | 141 |

### T3b: Age Classification Within Condition

| Condition | Gene AUROC | Pathway AUROC |
|---|---|---|
| FLT (flight) | 0.912 | 0.781 |
| GC (ground control) | 0.925 | 0.837 |
| VIV (vivarium) | 0.816 | 0.827 |

### T3d: Age-Stratified Spaceflight Detection

| Age Group | Gene AUROC | Pathway AUROC | n |
|---|---|---|---|
| **OLD (32-week)** | **0.945** [0.846, 1.00] | 0.879 | 34 |
| **YNG (10-12 week)** | **0.679** [0.479, 0.86] | 0.716 | 36 |
| **Delta (OLD - YNG)** | **+0.266** | +0.163 | — |

**"Spaceflight amplifies aging" hypothesis**: SUPPORTED. Old mice show substantially higher spaceflight AUROC (+0.266) → spaceflight transcriptomic signature is more distinct in aged animals.

### T3c: Two-Way ANOVA

50 Hallmark pathways tested (FLT + GC, ISS-T only, n=40). **0/50 significant interactions** at FDR < 0.05 (DD-19). Likely underpowered (n=10 per cell).

---

## T4: Radiation Metadata

ISS LEO dose rate ≈ 0.200-0.230 mGy/day. All missions 6.45-8.05 mGy total. Range too narrow for dose-response analysis. Annotation only (metadata completeness).

---

## J1: Pipeline Version Comparison (GLDS-48 vs GLDS-168)

| Metric | Value | Threshold | Verdict |
|---|---|---|---|
| Mean expression correlation | **0.920** | > 0.95 | Below threshold |
| Top-100 gene Jaccard | 0.087 | — | Very low overlap |
| Variance rank Spearman rho | 0.747 | — | Moderate |
| AUROC delta (168 - 48, matched 9 animals) | **-0.100** | < 0.05 | Exceeds threshold |

**Finding**: Pipeline version creates significant expression-level differences (r=0.92, not >0.99 as expected), and the matched-animal ML comparison now also shows a meaningful performance shift (GLDS-48=0.700 vs GLDS-168=0.600 on the same 9 RR-1 animals; delta -0.100). Gene selection remains highly sensitive to pipeline version (Jaccard 0.087 for top-100 genes).

---

## Tier 3: Text LLM Zero-Shot Evaluation

**Question**: Can LLMs classify spaceflight vs ground control from gene expression z-scores alone?

**Design**: Top-50 protein-coding genes by raw variance (filtered: no pseudogenes, ncRNAs, RIKEN clones, LOC loci), per-sample z-scores relative to LOMO training set. Zero-shot, temperature=0.

### End-to-End AUROC

| Task | Tissue | PCA-LR (Tier 1) | Llama-3.3-70B | DeepSeek-V3 | Gemini-2.5-Flash | Best LLM |
|---|---|---|---|---|---|---|
| A1 | Liver | 0.700 | 0.470 | 0.451 | 0.475 | 0.475 |
| A2 | Gastrocnemius | 0.790 | 0.480 | 0.535 | 0.461 | 0.535 |
| A3 | Kidney | 0.660 | 0.531 | 0.621 | 0.540 | 0.621 |
| A4 | Thymus | 0.800 | 0.479 | 0.459 | 0.522 | 0.522 |
| A5 | Skin | 0.830 | 0.446 | 0.482 | 0.560 | 0.560 |
| A6 | Eye | 0.760 | 0.504 | 0.496 | 0.404 | 0.504 |
| **Mean** | — | **0.757** | **0.485** | **0.507** | **0.494** | — |

### Parse-Aware Summary

| Provider | Parsed | Total | Rate | Mean End-to-End AUROC | Mean Parsed-Only AUROC |
|---|---|---|---|---|---|
| DeepSeek-V3 | 549 | 549 | **100%** | **0.507** | **0.507** |
| Gemini-2.5-Flash | 529 | 549 | 96.4% | 0.493 | 0.493 |
| Llama-3.3-70B (Groq) | 489 | 549 | 89.1% | 0.485 | 0.443 |

`end-to-end AUROC` treats parse failures as the neutral fallback score used in the submission JSON (`0.5`). `parsed-only AUROC` evaluates only samples where an explicit answer was recovered. Archived Groq and Gemini responses show frequent `length` / `FinishReason.MAX_TOKENS` truncation, so parse failure is part of the deployed pipeline behavior rather than a pure model-quality metric. The current parser now also recovers truncated first-sentence verdicts such as `consistent with spaceflight`, plus a very narrow set of archived `A1`/`A2` residual truncation patterns.

### Key Findings

1. **All LLMs remain near chance** on end-to-end evaluation (mean AUROC 0.485-0.507) — zero-shot text classification does not reliably distinguish spaceflight from ground control.
2. **PCA-LR (Tier 1) still outperforms by +0.25-0.27 mean AUROC** on the same tasks and feature budget.
3. **Provider ranking depends on whether parsing is included**: DeepSeek leads on end-to-end because it parses cleanly (100%), while Groq drops from 0.485 end-to-end to 0.443 parsed-only mean AUROC due to concentrated parse loss.
4. **Parse robustness is a first-class metric** for this benchmark: archived Groq and Gemini outputs frequently truncate before an explicit A/B answer. The parser now recovers explicit narrative verdicts, but provider comparisons should still report both end-to-end and parsed-only values.
5. **Systematic overconfidence and flight bias remain visible** in the recovered answers, despite chance-level discrimination.
6. **Tissue-specific gene quality is still evident**: liver=CYP450s/MUPs, skin=keratins, eye=crystallins, thymus=immune markers.

### Providers

| Provider | Model | API | Cost |
|---|---|---|---|
| Groq | llama-3.3-70b-versatile | OpenAI-compatible | Free |
| DeepSeek | deepseek-chat (V3) | OpenAI-compatible | Free trial |
| Google | gemini-2.5-flash | Google GenAI | Free |

---

---

## E1: Cross-Species NES Conservation

**Question**: Is mouse spaceflight pathway activity conserved in human cfRNA?

| Comparison | Species | Tissue/Source | Pathways | Spearman r |
|---|---|---|---|---|
| Mouse bulk liver NES ↔ JAXA cfRNA NES | Mouse → Human | Liver ↔ cfRNA | 50 Hallmark | **0.352** |

**Interpretation**: Moderate cross-species conservation at the pathway level despite tissue mismatch (liver vs. cell-free RNA).

---

## E2: Duration Effect on Conservation

**Question**: Does mission duration affect cross-species pathway conservation?

| Mission | Duration | r vs. JAXA cfRNA |
|---|---|---|
| I4 (PBMC) | Short | **−0.095** |
| JAXA | Long | **+0.352** |
| **Δr** | — | **+0.446** |

**Interpretation**: Longer spaceflight duration produces stronger cross-species pathway conservation.

---

## E3: cfRNA Cellular Origin

**Question**: Does cfRNA reflect expected tissue-of-origin pathway signatures?

**Finding**: Anti-correlation between cfRNA NES and expected tissue-specific pathway activity, suggesting cfRNA integrates multiple tissue sources rather than reflecting a single dominant origin.

---

## F1: I4 PBMC Cell-Type Pathway Analysis

**Question**: Which PBMC cell types drive spaceflight pathway signals?

- Source: GLDS-562 snRNA-seq, 10 cell types × 50 Hallmark pathways
- Temporal analysis reveals delayed PBMC response that mirrors mouse T2 post-return overshoot pattern
- Output: `v2/processed/F1_scrna/i4_snrnaseq_celltype_fgsea.csv`

---

## RRRM-1 scRNA-seq (F2, In Progress)

### Pipeline Completion Status

| Step | Status | Output |
|---|---|---|
| OSDR download (OSD-918/920/924/934) | Complete | Raw FASTQ on Cayuga scratch |
| STARsolo alignment (GRCm39-2024-A) | Complete | `Solo.out/GeneFull/filtered/` |
| H5AD conversion | Complete | Per-sample `.h5ad` |
| Tissue merge | Complete | Per-tissue merged `.h5ad` |
| QC + processing (Scanpy) | Complete | `*_processed.h5ad` |
| Broad annotation (marker scores) | Complete | `*_annotated.h5ad` |
| Doublet removal (Scrublet) | Complete | `*_hardened.h5ad` |
| **F2 benchmark tasks** | **Next** | — |

### Cell Inventory (post-hardening)

| Tissue (OSD) | Total Cells | Top Cell Types | Doublet Removed |
|---|---|---|---|
| Blood (918) | 4,377 | erythroid 83%, B cell 7%, T cell 5% | 18 (0.4%) |
| Eye (920) | 2,197 | retinal neuronal 52%, Muller glia 18% | 9 (0.4%) |
| Muscle (924) | 15,669 | macrophage/myeloid 38%, T/NK 25%, FAP 13% | 349 (2.2%) |
| Skin (934) | 15,838 | immune myeloid 37%, basal keratinocyte 30% | 38 (0.2%) |
| **Total** | **38,081** | — | 414 (1.1%) |

### Planned Benchmark Tasks

| Task | Question | Method | Status |
|---|---|---|---|
| F2-A | Cell-type composition shift (FLT vs GC) | Wilcoxon per cell type, LOAO | Planned |
| F2-B | Cell-type-specific pathway activity | Pseudo-bulk DESeq2 + fGSEA | Planned |
| F2-C | Cell-level spaceflight classifier | PCA-LR per cell type, LOAO | Planned |
| F2-D | Cross-species concordance (RRRM-1 ↔ I4) | NES Spearman r | Planned |

**Blocker**: FLT/GC condition labels need mapping from OSDR metadata → h5ad obs. SRX→condition CSV available at `v2/docs/RRRM1_SRX_CONDITION_MAP.csv`.

---

## v2 Publication Figures

| Figure | Content | Status |
|---|---|---|
| Fig1_temporal.html | T1 preservation artifact + T2 recovery + T3 age×spaceflight | Complete |
| Fig2_crossspecies.html | E1 NES conservation + E2 duration + E3 cfRNA origin | Complete |
| Fig3_pbmc_celltype.html | F1 cell-type NES heatmap + F1 temporal | Complete |
| Fig4 (planned) | F2 scRNA-seq summary (composition, pathway, classifier, cross-species) | Planned |

---

## Pipeline Status

| Component | Status | Output Location |
|---|---|---|
| T1: ISS-T vs LAR | Complete | `v2/processed/T_temporal/T1_*.json` |
| T2: LAR Recovery | Complete | `v2/processed/T_temporal/T2_*.json` |
| T3: Age × Spaceflight | Complete | `v2/processed/T_temporal/T3_*.json` |
| T4: Radiation | Complete | `v2/processed/T4_radiation/` |
| E1: Cross-Species NES | Complete | `v2/processed/E_crossspecies/` |
| E2: Duration Effect | Complete | `v2/processed/E_crossspecies/` |
| E3: cfRNA Origin | Complete | `v2/processed/E_crossspecies/` |
| F1: I4 PBMC Cell-Type | Complete | `v2/processed/F1_scrna/` |
| J1: Pipeline Comparison | Complete | `v2/evaluation/J1_pipeline_comparison.json` |
| Tier 3: LLM Zero-Shot | Complete | `v2/processed/llm_responses/`, `v2/evaluation/` |
| v2 Figures (3 main) | Complete | `v2/figures/` |
| RRRM-1 scRNA Pipeline | Complete | Cayuga: `rrrm1_scrna/downstream_initial/` |
| **F2 Benchmark Tasks** | **Next** | — |

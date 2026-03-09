# GeneLab_benchmark v2.0 — Results Summary

Generated: 2026-03-07

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

### Classification Recovery

| Mission | FLT_ISS-T flight_prob | FLT_LAR flight_prob | Interpretation |
|---|---|---|---|
| RR-6 | 1.000 | **0.185** | Strong recovery |
| RR-8 | 1.000 | **0.404** | Moderate recovery |

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
| AUROC delta (168 - 48) | **+0.049** | < 0.05 | Marginal |

**Finding**: Pipeline version creates significant expression-level differences (r=0.92, not >0.99 as expected), but downstream ML performance is minimally affected (AUROC delta +0.049). Gene selection is highly sensitive to pipeline version (Jaccard 0.087 for top-100 genes).

---

## Tier 3: Text LLM Zero-Shot Evaluation

**Question**: Can LLMs classify spaceflight vs ground control from gene expression z-scores alone?

**Design**: Top-50 protein-coding genes by raw variance (filtered: no pseudogenes, ncRNAs, RIKEN clones, LOC loci), per-sample z-scores relative to LOMO training set. Zero-shot, temperature=0.

### Per-Task AUROC

| Task | Tissue | PCA-LR (Tier 1) | Llama-3.3-70B | DeepSeek-V3 | Gemini-2.5-Flash | Best LLM |
|---|---|---|---|---|---|---|
| A1 | Liver | 0.700 | 0.527 | 0.436 | 0.523 | 0.527 |
| A2 | Gastrocnemius | 0.790 | 0.544 | 0.514 | 0.438 | 0.544 |
| A3 | Kidney | 0.660 | 0.440 | 0.495 | 0.494 | 0.495 |
| A4 | Thymus | 0.800 | 0.533 | 0.421 | 0.602 | 0.602 |
| A5 | Skin | 0.830 | 0.451 | 0.467 | 0.580 | 0.580 |
| A6 | Eye | 0.760 | 0.407 | 0.492 | 0.393 | 0.492 |
| **Mean** | — | **0.757** | **0.484** | **0.471** | **0.505** | — |

### Parse Rates

| Provider | Parsed | Total | Rate | max_tokens |
|---|---|---|---|---|
| DeepSeek-V3 | 549 | 549 | **100%** | 1000 |
| Gemini-2.5-Flash | 512 | 549 | 93.3% | 1000 |
| Llama-3.3-70B (Groq) | 489 | 549 | 89.1% | 500 |

Parse failures mostly from truncated chain-of-thought reasoning (Groq at 500 tokens) or verbose responses without explicit A/B answer (Gemini).

### Key Findings

1. **All LLMs at chance level** (AUROC 0.47–0.51 mean) — zero-shot text classification cannot distinguish spaceflight from ground control
2. **PCA-LR (Tier 1) outperforms by +0.25–0.29** mean AUROC — numerical ML on same features dramatically superior
3. **No provider consistently best** — ranking varies by tissue (noise-level variation)
4. **Systematic overconfidence**: All models predict 0.80–0.95 confidence regardless of correctness
5. **Flight bias**: Models tend to predict Flight (A) more often than Ground (B), suggesting spaceflight-related training data bias
6. **Tissue-specific gene quality confirmed**: Liver=CYP450s/MUPs, Skin=Keratins, Eye=Crystallins, Thymus=immune markers

### Providers

| Provider | Model | API | Cost |
|---|---|---|---|
| Groq | llama-3.3-70b-versatile | OpenAI-compatible | Free |
| DeepSeek | deepseek-chat (V3) | OpenAI-compatible | Free trial |
| Google | gemini-2.5-flash | Google GenAI | Free |

---

## Pipeline Status

| Component | Status | Output Location |
|---|---|---|
| T1: ISS-T vs LAR | Complete | `v2/processed/T_temporal/T1_*.json` |
| T2: LAR Recovery | Complete | `v2/processed/T_temporal/T2_*.json` |
| T3: Age × Spaceflight | Complete | `v2/processed/T_temporal/T3_*.json` |
| T4: Radiation | Complete | `v2/processed/T4_radiation/` |
| J1: Pipeline Comparison | Complete | `v2/evaluation/J1_pipeline_comparison.json` |
| Tier 3: LLM Zero-Shot | **Complete** | `v2/processed/llm_responses/`, `v2/evaluation/` |
| fGSEA v2: Additional DBs | Pending | — |

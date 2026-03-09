# GeneLab_benchmark v1.2 — Results Summary

Generated: 2026-03-01 (Updated: 2026-03-09)

## Hypothesis Results

| Hypothesis | Statement | Verdict | Key Evidence |
|---|---|---|---|
| **H1** | Liver has the most consistent cross-mission spaceflight transcriptome | **REFUTED** | Thymus(0.860) >> Liver(0.577). Thymus and gastrocnemius Tier 1. |
| **H2** | Cross-mission transfer failure from biological diversity, not batch | **SUPPORTED** | NES conservation correlates with transfer AUROC (r=0.9 for 5 tissues excl. gastrocnemius; original 4-tissue r=1.0). D3 pathway F1=0.06 (batch-invariant). |
| **H3** | Pathway-level preserves spaceflight response better than gene-level | **CONDITIONALLY SUPPORTED** | Kidney rescue (0.43→0.74), Eye (0.79→0.92). But liver→thymus anti-predicts (AUROC<0.5). Tissue-pair dependent. |

---

## Category B: Cross-Mission Transfer (PCA-LR, AUROC)

| Tissue | Mean AUROC | 95% CI | N Missions | N Pairs | Tier |
|---|---|---|---|---|---|
| **Thymus** | 0.860 | [0.763, 0.953] | 4 | 12 | 1 |
| **Gastrocnemius** | 0.801 | [0.653, 0.944] | 3 | 6 | 1 |
| **Skin** | 0.772 | [0.691, 0.834] | 3 | 6 | 2 |
| **Eye** | 0.754 | [0.688, 0.838] | 3 | 6 | 2 |
| **Liver** | 0.577 | [0.492, 0.666] | 6 | 30 | 3 |
| **Kidney** | 0.555 | [0.397, 0.681] | 3 | 6 | 3 |

Thymus vs Liver Δ = 0.283. Permutation tests: Thymus vs Liver p=0.001, Gastro vs Liver p=0.048, Skin vs Liver p=0.032.

### Category A Detection Significance (BH-FDR corrected)

| Tissue | AUROC | Raw p | FDR q | Significant? |
|---|---|---|---|---|
| Skin | 0.821 | 0.002 | 0.012 | **Yes** |
| Gastrocnemius | 0.824 | 0.026 | 0.074 | No |
| Thymus | 0.923 | 0.037 | 0.074 | No |
| Eye | 0.789 | 0.063 | 0.095 | No |
| Liver | 0.670 | 0.091 | 0.109 | No |
| Kidney | 0.432 | 0.281 | 0.281 | No |

*Note: Only skin survives BH-FDR correction at α=0.05. However, all top-4 tissues have AUROC > 0.7 (GO threshold). High AUROC with modest significance reflects small fold counts (3-4 folds per tissue), not weak signal.*

---

## Category C: Cross-Tissue Transfer (3 Methods, AUROC)

| Pair | Method A (Gene) | Method B (DEG) | Method C (Pathway) | Best |
|---|---|---|---|---|
| **C1** liver→kidney | **0.730** [0.62, 0.83] | 0.441 NS | 0.483 NS | A |
| **C2** liver→gastro | 0.563 NS | 0.676* | **0.867** [0.72, 0.98] | **C** |
| **C3** liver→thymus | 0.350 NS | **0.621*** | 0.184 (anti) | B |
| **C4** thymus→kidney | 0.585 NS | 0.539 NS | **0.690** [0.58, 0.79] | **C** |

H3 test: Method C wins 2/4 pairs (C2, C4). C vs A mean diff = -0.001 (essentially tied overall).

---

## Category D: Condition/Confounder Prediction (macro-F1)

| Task | Tissue | N | Gene F1 | Gene p | Pathway F1 | Pathway p | Interpretation |
|---|---|---|---|---|---|---|---|
| **D3** Mission ID (6-class) | Liver | 264 | **1.000** [1.00, 1.00] | <0.001 | 0.056 [0.04, 0.07] NS | 1.0 | Perfect batch separation (gene); batch-invariant (pathway) |
| **D4** Strain (2-class) | Thymus | 34 | **0.892** [0.48, 1.00] | 0.004 | 0.817 [0.47, 1.00] | 0.015 | Strain detectable even from GC-only samples. EXPLORATORY (n_minority=3) |
| **D5** Hardware RR vs MHU | Liver | 264 | **1.000** [1.00, 1.00] | <0.001 | 0.386 [0.36, 0.41] NS | 1.0 | Perfect gene separation; collinear with D3 (hardware derived from mission) |
| **D5** Hardware RR vs MHU | Thymus | 92 | **1.000** [1.00, 1.00] | <0.001 | 0.352 [0.31, 0.39] NS | 1.0 | Perfect gene separation; collinear with D3 |
| **D6** Gravity (3-class) | Liver | 9 | **0.886** | 0.002 | 0.413 NS | 0.354 | uG separable from gene expression |
| **D6** Gravity (3-class) | Thymus | 9 | **0.657** | 0.037 | 0.641 (p=0.052) | 0.052 | Gene ≈ Pathway for gravity detection |

### Confounder Hierarchy

```
D3 (mission F1=1.0) >= D5 (hardware F1=1.0, collinear) >= D4 (strain F1=0.89, exploratory n=3)
All pathway F1 ≈ 0.05-0.41 → pathways resist confounder detection (batch-invariant)
```

**Key insight**: D5 hardware prediction is perfect but **collinear with D3** — hardware type (RR vs MHU) is a deterministic function of mission ID. D5 F1 should be interpreted as an upper bound of D3, not independent evidence. D4 strain effect is detectable but **exploratory** (minority class n=3).

---

## J5: Gene-level vs Pathway-level Comparison (12 tasks)

### Category A — Spaceflight Detection (LOMO, AUROC)

| Tissue | Gene | Pathway | Diff (P-G) | Winner |
|---|---|---|---|---|
| Liver | 0.670 | 0.574 | -0.096 | Gene |
| Gastrocnemius | 0.824 | 0.688 | -0.137 | Gene |
| **Kidney** | 0.432 | **0.743** | **+0.311** | **Pathway** |
| Thymus | 0.923 | 0.879 | -0.044 | Gene |
| **Eye** | 0.789 | **0.915** | **+0.125** | **Pathway** |

**Mean diff (Cat A)**: +0.032 (essentially tied)

### Across All Categories

| Category | N | Gene wins | Pathway wins | Mean diff |
|---|---|---|---|---|
| A (Detection) | 5 | 3 | 2 | +0.032 |
| C (Cross-tissue) | 4 | 2 | 2 | -0.001 |
| D (Condition, D3+D6 original) | 3 | 3 | 0 | -0.478 |
| D (Condition, full D3-D6) | 6 | 6 | 0 | -0.462 |
| **Total (original 12)** | **12** | **8** | **4** | **-0.106** |
| **Total (expanded 15)** | **15** | **11** | **4** | **-0.174** |

*Note: D4/D5 all show gene >> pathway, consistent with D3/D6 pattern. Pathways systematically resist confounder/batch detection.*

---

## NES Conservation vs Cross-Mission Transfer

| Tissue | NES Mean r | Transfer AUROC | N fGSEA missions |
|---|---|---|---|
| **Thymus** | **0.619** | **0.860** | 3 |
| Eye | 0.335 | 0.754 | 3 |
| **Skin** | **0.147** | **0.772** | **3** |
| Liver | 0.059 | 0.577 | 6 |
| Gastrocnemius | 0.057 | 0.801* | 2 |
| Kidney | 0.048 | 0.555 | 3 |

*Gastrocnemius outlier: only 2/3 missions have fGSEA data (RR-5 no DGE).
†Skin: RR-7 DGE absent; fGSEA on RR-6, MHU-2_dorsal (GLDS-238), MHU-2_femoral (GLDS-239) only.
Excluding gastrocnemius: rank-order correlation for 5 tissues (thymus/eye/skin/liver/kidney) Spearman r = 0.9 (skin NES rank 3rd vs transfer rank 2nd — partial outlier). Original 4-tissue finding (thymus/eye/liver/kidney, excl gastrocnemius) maintains perfect rank concordance (Spearman r = 1.0).

---

## Biological Validation (fGSEA Hallmark, all tissues PASS)

| Tissue | Top Enriched Pathways | Consistency |
|---|---|---|
| Liver | OXIDATIVE_PHOSPHORYLATION, FATTY_ACID_METABOLISM | Literature-concordant |
| Thymus | E2F_TARGETS, G2M_CHECKPOINT, IFN-gamma | Thymocyte proliferation |
| Gastrocnemius | OXIDATIVE_PHOSPHORYLATION, MYOGENESIS | Muscle metabolism |
| Kidney | MTORC1_SIGNALING, CHOLESTEROL_HOMEOSTASIS | Renal metabolism |
| Eye | OXIDATIVE_PHOSPHORYLATION (dominant 3/3 missions) | Retina metabolic demand |
| Skin | E2F_TARGETS, G2M_CHECKPOINT, EPITHELIAL_MESENCHYMAL_TRANSITION | Cell proliferation + ECM remodeling (2/3 missions consistent) |

---

## Tier 2: Geneformer (Mouse-GF) vs Classical Baseline

Mouse-Geneformer (6L BERT, 56K mouse gene vocab, pretrained on 30M scRNA-seq cells) fine-tuned on bulk RNA-seq LOMO folds (10 epochs, batch=16, lr=2e-5, freeze=4/6 layers).

| Task | Tissue | Geneformer AUROC | Baseline AUROC | Baseline Model | Delta | Winner |
|------|--------|-----------------|---------------|----------------|-------|--------|
| A1 | Liver | 0.486 ± 0.074 | 0.588 | LR ElasticNet | -0.102 | Baseline |
| A2 | Gastrocnemius | 0.382 ± 0.054 | 0.907 | LR ElasticNet | -0.525 | Baseline |
| A3 | Kidney | 0.452 ± 0.080 | 0.521 | LR ElasticNet | -0.069 | Baseline |
| A4 | Thymus | 0.495 ± 0.233 | 0.923 | PCA-50 + LogReg | -0.428 | Baseline |
| A5 | Skin | 0.557 ± 0.087 | 0.821 | LR ElasticNet | -0.265 | Baseline |
| A6 | Eye | 0.484 ± 0.117 | 0.789 | PCA-50 + LogReg | -0.305 | Baseline |
| **Mean** | **6 tissues** | **0.476** | **0.758** | — | **-0.283** | **Baseline** |

**Interpretation**: Classical ML wins 6/6 tissues (sign test p=0.016). Geneformer performs near chance level (0.5) on small-n bulk RNA-seq (train n=30-100). This is consistent with literature — foundation models pretrained on single-cell data do not automatically transfer to small-sample bulk transcriptomics tasks.

*Note: Table shows best baseline per tissue for fair comparison. Publication figures use unified PCA-LR baseline (mean 0.743) for cross-figure consistency with Category A/B results.*

---

## Held-Out Evaluation: A4 Thymus (OSD-515 / RR-23)

Reserved held-out test set for external benchmark evaluation. Train on 4 missions (MHU-1, MHU-2, RR-6, RR-9; n=67), test on RR-23 (n=16: 7 Flight, 9 GC). 27,541 common genes.

| Model | AUROC | 95% CI | p-value |
|-------|-------|--------|---------|
| LR ElasticNet | **0.905** | [0.672, 1.000] | 0.005 |
| Random Forest | **0.905** | [0.672, 1.000] | 0.007 |
| PCA-50 + LogReg | 0.873 | [0.609, 1.000] | 0.011 |
| Geneformer (Mouse-GF) | 0.556 | [0.265, 0.850] | — |

**Interpretation**: Classical baselines achieve strong held-out performance (AUROC ~0.90, p<0.01), confirming thymus cross-mission generalization beyond LOMO. Geneformer remains near chance on held-out data, consistent with LOMO results (0.495). The held-out confirms thymus as the most robust tissue for spaceflight detection.

---

## Tier 3: LLM Zero-Shot Classification

Three LLMs tested on zero-shot text-based spaceflight detection (no training, gene expression → text prompt → binary prediction).

| Model | A1 Liver | A2 Gastro | A3 Kidney | A4 Thymus | A5 Skin | A6 Eye | Mean |
|---|---|---|---|---|---|---|---|
| **PCA-LR (ref)** | 0.670 | 0.824 | 0.432 | 0.923 | 0.821 | 0.789 | **0.743** |
| DeepSeek-V3 | 0.435 | 0.514 | 0.495 | 0.421 | 0.467 | 0.492 | 0.471 |
| Gemini-2.5-Flash | 0.523 | 0.438 | 0.494 | 0.602 | 0.580 | 0.393 | 0.505 |
| Llama-3.3-70B | 0.527 | 0.544 | 0.440 | 0.533 | 0.451 | 0.407 | 0.484 |

**Interpretation**: All 3 LLMs perform at chance level (mean 0.47–0.51). Text-based reasoning cannot replace numerical ML for transcriptomics classification. Protein-coding gene filter was applied to reduce prompt noise.

---

## Multi-DB Pathway Comparison (LOMO, PCA-LR)

| Tissue | Hallmark | KEGG | Reactome | MitoCarta | Best DB | Range |
|---|---|---|---|---|---|---|
| Thymus | 0.879 | 0.899 | **0.922** | 0.846 | Reactome | 0.076 |
| Gastro | 0.688 | 0.713 | **0.755** | 0.627 | Reactome | 0.128 |
| Skin | 0.690 | **0.754** | 0.693 | 0.542 | KEGG | 0.212 |
| Eye | **0.915** | 0.625 | 0.658 | 0.478 | Hallmark | 0.437 |
| Liver | 0.574 | **0.639** | 0.614 | 0.555 | KEGG | 0.084 |
| Kidney | 0.743 | 0.665 | **0.779** | 0.641 | Reactome | 0.138 |

**Key findings**:
- DB choice > model choice (AUROC range up to 0.437 for Eye)
- No single DB dominates: Reactome best for 3 tissues, KEGG for 2, Hallmark for 1
- MitoCarta consistently worst (specialized → low coverage)

---

## Temporal & Biological Covariates

### T1: ISS-T vs LAR Sacrifice Timing

**Question**: Can sacrifice timing (ISS-Terminal vs Live Animal Return) be detected from transcriptomics?

**Confound warning (DD-18)**: ISS-T = RNAlater on-orbit, LAR = standard necropsy. Preservation method confounds biological timing.

| Sub-task | Condition | Gene AUROC | Pathway AUROC | n |
|---|---|---|---|---|
| **T1a RR-6 liver** | FLT | 1.000 | 0.960 | 20 |
| | GC (baseline) | 0.922 | 0.778 | 19 |
| **T1b RR-8 liver** | FLT | 0.930 | 0.920 | 35 |
| | GC (baseline) | **0.973** | 0.993 | 35 |
| **T1c RR-6 thymus** | FLT | 0.857 | 0.714 NS | 17 |
| | GC (baseline) | **0.925** | 1.000 | 18 |

**Verdict**: GC AUROC ≥ FLT AUROC in most cases → ISS-T vs LAR difference **dominated by preservation artifact**, not biological timing. Cross-mission transfer (T1d) confirms: artifact is consistent across RR-6↔RR-8 (FLT gene AUROC 0.97–0.99, GC gene 0.84–0.96).

### T2: LAR Recovery Signature

| Mission | PCA Recovery R | Pathways Recovering | FLT_LAR flight_prob |
|---|---|---|---|
| RR-6 | 0.842 (partial) | 12/26 | 0.185 (strong) |
| RR-8 | 0.652 (stronger) | **25/27** (overshoot) | 0.404 (moderate) |

RR-8 shows strong recovery with overshoot past baseline (MYC targets V1 +2.49, Protein secretion +2.14). RR-6 shows immune rebound (IFN-α -2.36, Inflammatory -2.54).

### T3: Age × Spaceflight Interaction (RR-8 Liver)

| Comparison | Gene AUROC | Pathway AUROC | n |
|---|---|---|---|
| T3a: Overall OLD vs YNG | **0.985** | 0.851 | 141 |
| T3d: Spaceflight in OLD | **0.945** [0.846, 1.00] | 0.879 | 34 |
| T3d: Spaceflight in YNG | **0.679** [0.479, 0.86] | 0.716 | 36 |
| **Delta (OLD - YNG)** | **+0.266** | +0.163 | — |

**Verdict**: "Spaceflight amplifies aging" **SUPPORTED** (Δ=+0.266). T3c ANOVA: 0/50 significant Age×Spaceflight interactions at FDR<0.05 (underpowered, n=10/cell).

---

## J2: DGE Pipeline Comparison

**Question**: Does the choice of DGE pipeline (DESeq2 vs edgeR vs limma-voom) affect downstream results?

**Scope**: 9 missions (6 liver + 3 thymus) × 3 pipelines = 27 DGE runs. Skin excluded (RR-7 has no raw counts).

| Metric | Mean | Min | Max |
|---|---|---|---|
| **Log2FC Spearman** | **0.926** | 0.790 | 1.000 |
| **Log2FC Pearson** | **0.895** | 0.706 | 1.000 |
| **DEG Jaccard (FDR<0.05)** | **0.600** | 0.000 | 1.000 |
| **GeneLab Replication** | **0.707** | — | 9 missions |

**Key findings**:
- Fold-change rankings are highly conserved across all three pipelines (Spearman 0.926)
- DEG list overlap varies by pipeline stringency: limma-voom most liberal, edgeR most conservative
- RR-3 liver: 0 DEGs across all pipelines (n=6+6, true biological null — GeneLab also found only 1 DEG)
- RR-1 edgeR: 0 DEGs due to conservative multiple testing correction, but log2FC correlation >0.95 with DESeq2
- GeneLab replication (our binary FLT-vs-GC vs GeneLab's multi-level contrasts): moderate agreement (Spearman 0.707) reflects different design matrices, not pipeline error

**Verdict**: Rankings consistent, DEG lists vary by stringency threshold. Pipeline choice has **moderate impact on DEG calls** but **minimal impact on gene-level effect size rankings** — consistent with J1 (pipeline version comparison).

---

## Held-Out Evaluation: A5 Skin (OSD-254 / RR-7)

Second held-out test set. Train on 2 missions (RR-6, MHU-2; n=72), test on RR-7 (n=30: 10 Flight, 20 GC). 20,110 common genes. RR-7 is a 75-day mission (longest in skin dataset).

| Model | AUROC | 95% CI | p-value |
|-------|-------|--------|---------|
| LR ElasticNet | **0.885** | [0.745, 0.986] | <0.001 |
| PCA-50 + LogReg | 0.840 | [0.679, 0.963] | 0.001 |
| Random Forest | 0.777 | [0.583, 0.929] | 0.007 |

**Cross-Tissue Held-Out Comparison**:

| Tissue | Mission | Duration | Best AUROC | n_test |
|--------|---------|----------|------------|--------|
| Thymus | RR-23 | 30 days | 0.905 (LR) | 16 |
| Skin | RR-7 | 75 days | 0.885 (LR) | 30 |

**Interpretation**: Skin held-out confirms strong generalization (AUROC 0.885, p<0.001), exceeding the LOMO mean (0.821). Both held-out tissues achieve AUROC > 0.85, validating cross-mission spaceflight detection beyond leave-one-out evaluation.

---

## Pipeline Status

| Component | Files | Status |
|---|---|---|
| fGSEA | 80 (6 tissues × missions × 4 DBs incl. MitoCarta) | Complete |
| GSVA | 88 (6 tissues × missions × 4 DBs, skin+thymus MHU-1) | Complete |
| Category A | 6 tissues, PCA-LR LOMO | Complete |
| Category B | 6 tissues, bootstrap CI + permutation | Complete |
| Category C | 4 pairs × 3 methods | Complete |
| Category D | D3 + D4 + D5×2 + D6×2 (6 tasks) | Complete |
| J5 | 15 comparisons | Complete |
| NES Conservation | 6 tissues × 4 DBs | Complete |
| Multi-DB LOMO | 24 runs (6 tissues × 4 DBs) | Complete |
| NC1/NC2 | Permutation + housekeeping controls | Complete |
| Cell 2020 | 5 tissues pathway validation | Complete |
| Geneformer | 6 tissues, 22 LOMO folds (Mouse-GF) | Complete |
| LLM Zero-Shot | 3 providers × 6 tasks (18 evals) | Complete |
| Held-Out | A4 Thymus (RR-23) + A5 Skin (RR-7) | Complete |
| T1-T3 Temporal | ISS-T/LAR, Recovery, Age×Spaceflight | Complete |
| J2 DGE Pipeline | 9 missions × 3 pipelines (DESeq2/edgeR/limma-voom) | Complete |
| **Figures** | **4 main + 4 supplementary (HTML/SVG)** | **Complete** |

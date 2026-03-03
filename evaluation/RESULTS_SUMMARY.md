# GeneLab_benchmark v1.0 — Results Summary

Generated: 2026-03-01

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

Permutation tests: Thymus vs Liver p=0.001, Gastro vs Liver p=0.048, Skin vs Liver p=0.032.

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

| Task | Gene | Pathway | p-value | Interpretation |
|---|---|---|---|---|
| **D3** Liver 6-class mission | **1.000** | 0.056 NS | <0.001 | Perfect batch separation (gene); batch-invariant (pathway) |
| **D6** Liver uG/AG/GC | **0.886** | 0.413 NS | 0.002 | uG separable from gene expression |
| **D6** Thymus uG/AG/GC | **0.657** | 0.641 (p=0.052) | 0.037 | Gene ≈ Pathway for gravity detection |

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
| D (Condition) | 3 | 3 | 0 | -0.478 |
| **Total** | **12** | **8** | **4** | **-0.106** |

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

## Pipeline Status

| Component | Files | Status |
|---|---|---|
| fGSEA | 60 (6 tissues × missions × 3 DBs) | Complete |
| GSVA | 54 (5 tissues × missions × 3 DBs) | Complete |
| Category B | 5 tissues, bootstrap CI + permutation | Complete |
| Category C | 4 pairs × 3 methods | Complete |
| Category D | D3 + D6×2 | Complete |
| J5 | 12 comparisons | Complete |
| NES Conservation | 6 tissues | Complete |
| Geneformer | Tokenized, HPC script ready | Awaiting HPC |

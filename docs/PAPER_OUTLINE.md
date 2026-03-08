# GeneLab Benchmark — Paper Outline (Draft)

**Target**: Genome Biology / Nature Methods
**Working Title**: "A Multi-Tissue Machine Learning Benchmark for Cross-Mission Generalization of Spaceflight Transcriptomic Signatures"

---

## Structure

### Title & Abstract (~250 words)
- Problem: No cross-mission ML benchmark for spaceflight transcriptomics; liver assumed as gold standard
- Approach: 6 tissues × 17 missions × ~450 samples, LOMO CV, 3-tier models (classical ML, Geneformer, text LLMs)
- Key results: (1) Thymus >> Liver (0.860 vs 0.577, p=0.001), (2) 17.8× pathway batch-resistance, (3) NES conservation predicts transfer (r=0.9), (4) Foundation models fail on bulk RNA-seq

---

### 1. Introduction (~1000 words)

1a. **Gap**: NASA OSDR has 400+ datasets analyzed in isolation. No benchmark for cross-mission generalization.
- Existing: Cell 2020 (Beheshti) descriptive atlas, SOMA (Nature 2024) multi-omics — neither tests generalization
- Liver-centric assumption untested by cross-mission ML

1b. **Why ML benchmarking matters**: DGE finds DEGs but can't measure generalization. LOMO directly tests reproducibility.

1c. **Pre-registered hypotheses**:
- H1: Liver most consistent (REFUTED)
- H2: Transfer failure from biology, not batch (SUPPORTED)
- H3: Pathways preserve spaceflight better (CONDITIONALLY SUPPORTED)

---

### 2. Results (~3000 words)

2a. **Benchmark design & data curation** (Fig 1)
- 6 tissues, 17 missions, 24 OSD studies, ~450 binary samples
- LOMO: mission = independence unit. Feature selection within loop (DD-03)
- 7 evaluation categories: A (detection), B (transfer), C (cross-tissue), D (confounder), J5 (gene vs pathway), NC (negative control), Validation

2b. **Thymus, not liver, has the most generalizable spaceflight signature** (Fig 2)
- Category A LOMO: Thymus 0.923 > Gastro 0.907 > Skin 0.821 > Eye 0.811 > Liver 0.653 > Kidney 0.593
- Category B transfer: Thymus 0.860 >> Liver 0.577 (p=0.001)
- H1 REFUTED. More data (liver 6 missions) reveals more heterogeneity, not more consistency

2c. **Pathway features are 17.8× more batch-resistant than gene features** (Fig 3)
- D3 mission identity: Gene F1=1.000, Pathway F1=0.056
- D5 hardware: Gene F1=1.000, Pathway F1=0.386
- Pathway rescue: Kidney gene 0.432 → pathway 0.743 (+0.311), Eye gene 0.789 → pathway 0.915 (+0.126)
- H2 SUPPORTED: batch effects are gene-driven, pathway-absorbed

2d. **NES conservation predicts cross-mission transfer success** (Fig 4)
- fGSEA NES pairwise Spearman r per tissue → rank-correlates with Category B AUROC
- 5-tissue r=0.9, 4-tissue r=1.0
- Practical: fGSEA on 2-3 missions predicts ML feasibility without training models

2e. **Cross-tissue transfer depends on biological conservation** (Fig 5, supplementary detail)
- Pathway transfer wins 2/4 pairs (C2 liver→gastro 0.867, C4 thymus→kidney 0.690)
- Anti-prediction: C3 liver→thymus 0.184 — pathway features actively mislead
- H3 CONDITIONALLY SUPPORTED

2f. **Foundation models and text LLMs fail on bulk RNA-seq** (Fig 6)
- Mouse-Geneformer: mean 0.476 vs Classical ML 0.758 (delta=-0.283), 6/6 tissues classical wins
- Text LLMs (3 providers, zero-shot): mean 0.47-0.51, all at chance level
- Held-out RR-23 thymus: Classical 0.905 vs Geneformer 0.556
- scRNA pretraining does not transfer to small-n bulk; zero-shot text reasoning cannot replace numerical ML

2g. **Validation** (supplementary)
- Cell 2020 concordance: 71.7% pathway direction match, 47× gene enrichment
- NC1 permutation: mean AUROC=0.50 ± 0.03 (PASS)
- NC2 housekeeping: 0.49-0.55 (PASS, thymus/eye elevated → batch leakage noted)
- Held-out RR-23: AUROC=0.905, p=0.005

---

### 3. Discussion (~1500 words)

3a. **Why thymus > liver**: Stereotyped response (E2F/G2M up, IFN-γ down) with large effect sizes. Liver heterogeneity from 6 different mission contexts. Thymic involution as conserved stress response.

3b. **Gene vs pathway complementarity**: Genes detect confounders (useful for QC), pathways resist batch (useful for transfer). Not substitutes — complementary feature levels.

3c. **NES conservation as a predictive biomarker**: Pathway-level biology drives gene-level generalization. Practical screening tool for new tissues/species.

3d. **Foundation model limitations**: scRNA→bulk domain gap. Small-n prevents effective fine-tuning. Rank-based tokenization loses quantitative info. Not a paradigm failure — domain mismatch. Text LLMs confirm: numerical ML required.

3e. **Limitations**: Small sample sizes (n=6-40). Per-mission normalization (not joint). Strain confounding in MHU-1 thymus. Single species (mouse). One FM tested (Geneformer). Temporal confounding (ISS-T vs LAR, explored in supplementary).

3f. **Future directions**: Cross-species (human, C. elegans). Domain-adapted FMs. Joint normalization. Radiation dose-response (requires ground irradiation).

---

### 4. Methods (~2500 words)

4a. Data curation: NASA OSDR sources, 24 OSD studies, preprocessing (DESeq2, log2, variance filter)
4b. Label assignment: Binary Flight/Ground, excluded BC/AG, keyword matching
4c. LOMO framework: Mission as independence unit, feature selection within loop
4d. Baseline models: LR (ElasticNet), RF (500 trees), PCA-50 + LR. StandardScaler per fold
4e. Geneformer: Architecture (6L BERT, 256 hidden, 56K vocab), fine-tuning (freeze 4/6, 10 epochs, lr=2e-5, A40 GPU)
4f. Text LLM: 3 providers (Llama-3.3-70B, DeepSeek-V3, Gemini-2.5-Flash), per-sample z-score prompts, temperature=0
4g. Pathway analysis: fGSEA (Hallmark/KEGG/Reactome), GSVA (sample-level), NES conservation
4h. Statistical framework: Bootstrap CI (N=2000), permutation test (N=1000), pseudocount correction
4i. Evaluation categories: A-D, J5, NC, validation criteria (Go/No-Go thresholds)
4j. External validation: Cell 2020 concordance, SHAP overlap, negative controls

---

### Figures (Main, 6 panels)

| Fig | Content | Key Message |
|-----|---------|-------------|
| 1 | Study design: tissues × missions matrix, LOMO diagram, 3-tier pipeline | Benchmark overview |
| 2 | Category A+B results: AUROC bars by tissue (CI), transfer heatmaps | Thymus >> Liver |
| 3 | D3 batch: Gene PCA (perfect separation) vs Pathway PCA (mixed) | 17.8× batch resistance |
| 4 | NES conservation vs transfer AUROC scatter (r=0.9 line) | Pathway biology predicts ML |
| 5 | Gene vs Pathway comparison: J5 summary + kidney/eye rescue | Complementary features |
| 6 | Geneformer + LLM vs Classical ML paired bars (6 tissues) | Classical ML wins |

### Supplementary

- S1: Cross-tissue transfer (Category C, 4 pairs × 3 methods)
- S2: fGSEA enrichment heatmaps (6 tissues × missions × 3 DBs)
- S3: Cell 2020 concordance detail
- S4: Negative controls (NC1, NC2)
- S5: RR-23 held-out detail
- S6: Confounder hierarchy (D3-D6)
- S7: Geneformer training curves
- S8: LLM prompt template + parse rates
- S9: v2 temporal analysis (T1-T4 summary — ISS-T/LAR, recovery, age interaction, radiation)
- S10: J1 pipeline version comparison (GLDS-48 vs GLDS-168)

---

### Data/Code Availability
- NASA OSDR (public): 24 OSD studies
- HuggingFace: jang1563/genelab-benchmark (processed data + task definitions)
- GitHub: analysis code, evaluation framework (MIT license)

---

## Word Budget (Genome Biology Research Article)

| Section | Target | Notes |
|---------|--------|-------|
| Abstract | 250 | Structured |
| Introduction | 1000 | 3 hypotheses |
| Results | 3000 | 7 subsections |
| Discussion | 1500 | 6 points |
| Methods | 2500 | 10 subsections |
| **Total** | **~8250** | Genome Biology limit ~8000-10000 |

---

## Novel Claims (reviewer-ready)

1. **First multi-tissue cross-mission ML benchmark** for spaceflight transcriptomics
2. **Thymus > Liver** with statistical significance (p=0.001) — overturns consensus
3. **17.8× pathway batch-resistance** — first quantitative measurement
4. **NES conservation predicts transfer** (r=0.9) — previously unreported relationship
5. **First systematic FM evaluation** on spaceflight bulk RNA-seq (Geneformer + 3 text LLMs)
6. **Pathway rescue** of kidney/eye — practical demonstration of feature-level complementarity

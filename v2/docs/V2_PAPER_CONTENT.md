# GeneLab Benchmark v2 — Paper Content

*Draft: 2026-03-09. English draft for manuscript preparation.*
*Analysis dates: 2026-03-07 (T1/T2/T3, RR-6/RR-8 liver). Data: NASA OSDR.*

---

## Section 1: Mouse Temporal Dynamics

### 1.1 Sample Preservation Artifacts in Multi-Timepoint Studies (T1)

The RR-8 mission collected liver tissue using two parallel preservation protocols: samples fixed in RNAlater immediately after in-flight euthanasia (ISS-T; n=20 flight, n=20 ground) and samples processed by vivarium necropsy after return to Earth (LAR; n=15 flight, n=15 ground). We asked whether ISS-T vs. LAR grouping reflects biological spaceflight response or preservation artifacts.

A classifier trained to distinguish ISS-T from LAR samples achieved AUROC = 0.930 (95% CI: 0.787–1.0) in FLT samples and AUROC = 0.973 (0.905–1.0) in GC (ground control) samples, yielding an excess signal of −0.043 (FLT − GC). The near-zero excess indicates that the ISS-T vs. LAR transcriptomic difference is present equally in both flight and ground conditions and therefore reflects the preservation method rather than spaceflight biology. This finding is consistent with published reports of RNAlater-induced transcriptomic perturbations in liver tissue (Galley et al. 2021; PMID 33376967).

**Key result:** ISS-T vs. LAR AUROC excess ≈ 0 → preservation artifact, not spaceflight signal.

Cross-mission transfer (T1d) confirmed that the FLT signature is highly consistent across missions: gene-level classifiers trained on RR-6 predicted RR-8 with AUROC = 0.987 (0.948–1.0), and RR-8→RR-6 achieved AUROC = 0.970 (0.879–1.0), demonstrating stable cross-mission transferability.

---

### 1.2 Post-Flight Recovery of Spaceflight-Induced Transcriptomic Perturbations (T2) ★★★

We leveraged the dual ISS-T/LAR sampling design in RR-8 as a natural experiment: ISS-T samples represent the in-flight transcriptomic state, while LAR samples collected after vivarium acclimation represent the post-return state (3–11 days on the ground). We characterized recovery at both the PCA trajectory level and individual pathway level.

#### PCA Trajectory Analysis

We computed PCA distances between six experimental groups (BSL_ISS-T, BSL_LAR, FLT_ISS-T, FLT_LAR, with the first three principal components explaining 12.8%, 8.2%, and 6.5% of variance, respectively). We defined a recovery ratio as:

> Recovery ratio = d(FLT_LAR, BSL_LAR) / d(FLT_ISS-T, BSL_ISS-T)

where the numerator captures the residual flight displacement after return under LAR processing and the denominator captures the in-flight displacement under ISS-T processing. In RR-8, this ratio was **0.652**, indicating that the post-return profile is closer to its baseline reference than the in-flight profile (compare: RR-6 recovery ratio = 0.842, reflecting the shorter duration RR-6 mission that may have induced less severe perturbations).

As a descriptive projection, a spaceflight detector trained on ISS-T samples assigned FLT_ISS-T a mean flight probability of 0.9997, while LAR samples were assigned 0.404 — a >2.4-fold reduction, close to the 0.5 random-assignment threshold. Because the FLT_ISS-T score is measured on the training contrast itself, we interpret this as supportive baselineward shift rather than held-out validation of recovery.

#### Pathway-Level Recovery

Of 27 Hallmark pathways with |Δ_flight| ≥ 0.10 NES units in RR-8, **25/27 (93%) showed recovery** (mean recovery fraction = 1.24, median = 1.27), with recovery fractions > 1.0 indicating overshoot beyond baseline. Two pathways did not recover: ALLOGRAFT_REJECTION and DNA_REPAIR.

Top overshoot pathways (recovery fraction, meaning post-return NES surpassed baseline in the opposite direction):

| Pathway | Δ_flight (NES) | Recovery fraction | Interpretation |
|---------|---------------|-------------------|----------------|
| MYC_TARGETS_V1 | −0.132 | 2.49× | Biosynthesis rebound |
| PROTEIN_SECRETION | −0.117 | 2.14× | Secretory pathway rebound |
| MYC_TARGETS_V2 | −0.275 | 1.92× | Biosynthesis rebound |
| UNFOLDED_PROTEIN_RESPONSE | −0.130 | 1.78× | ER stress rebound |
| OXIDATIVE_PHOSPHORYLATION | −0.294 | 1.75× | Mitochondrial recovery |
| GLYCOLYSIS | −0.131 | 1.75× | Metabolic recovery |
| MTORC1_SIGNALING | −0.236 | 1.67× | Anabolic rebound |
| ADIPOGENESIS | −0.212 | 1.79× | Lipid metabolism recovery |

The systematic overshoot in MYC targets and protein secretion pathways suggests compensatory transcriptional upregulation following in-flight suppression of biosynthetic programs, consistent with a refeeding-like response after stress removal.

**Key results:**
- PCA recovery ratio = 0.652 (65% trajectory recovery in RR-8)
- 25/27 pathways recover (93%), majority with overshoot
- Classifier flight probability drops from 0.9997 (ISS-T) to 0.404 (LAR) post-return
- MYC target overshoot (2.49×) implicates compensatory biosynthetic rebound

---

### 1.3 Age Amplifies Spaceflight Transcriptomic Response (T3) ★★★

The RR-8 mission included a unique aging cohort: both young adult (YNG; 10 to 12 week) and aged (OLD; 32 week) mice were flown simultaneously, enabling direct within-mission comparison of age-dependent spaceflight response.

#### Age Classification

Old and young mice were highly distinguishable transcriptomically (T3a): gene-level AUROC = 0.985 (95% CI: 0.962–0.999, permutation p < 0.0001), pathway-level AUROC = 0.851 (0.786–0.911). Within-condition age classification also remained highly accurate (T3b: FLT gene AUROC = 0.912, GC gene AUROC = 0.925), confirming that age imprints a robust and spaceflight-persistent transcriptomic signature.

#### Spaceflight Detection by Age Group

We trained spaceflight classifiers (FLT vs. GC) separately in each age group to test whether aging amplifies the transcriptomic response to spaceflight:

| Age group | Spaceflight AUROC (gene) | 95% CI | Permutation p |
|-----------|--------------------------|--------|---------------|
| OLD (32 week) | **0.945** | 0.846–1.0 | <0.0001 |
| YNG (10-12 week) | **0.679** | 0.479–0.860 | 0.033 |
| **Δ (OLD − YNG)** | **+0.266** | — | — |

The AUROC gap of +0.266 between OLD and YNG mice indicates that aged animals produce a substantially stronger and more detectable transcriptomic spaceflight response (T3d). Pathway-level results followed the same direction: OLD AUROC = 0.879 (0.744–0.978) vs. YNG = 0.716 (0.523–0.873).

#### Interaction Analysis

Two-way ANOVA for each pathway tested flight × age interaction across ISS-T samples (n=40; T3c). Zero of 50 Hallmark pathways showed a significant age × flight interaction after Benjamini-Hochberg FDR correction (minimum q = 0.365 for OXIDATIVE_PHOSPHORYLATION, nominal p = 0.007). This suggests that age amplifies spaceflight response in a global, non-pathway-specific manner rather than through specific interaction terms.

**Key results:**
- OLD mice: spaceflight AUROC = 0.945; YNG: 0.679; Δ = +0.266 (gene level)
- Pathway level: OLD 0.879 vs. YNG 0.716
- Age classification in FLT: gene AUROC = 0.912 (age signature maintained in space)
- No pathway-specific age × flight interaction (FDR > 0.05), suggesting global amplification
- Interpretation: aging broadly sensitizes hepatic transcriptome to spaceflight perturbation

---

## Section 2: Cross-Species Conservation (Category E) ★★

### E1: Mouse Liver vs. Human cfRNA NES Correlation

We compared Hallmark pathway NES values between mouse liver (mission-averaged across 6 NASA OSDR missions) and human plasma cfRNA from JAXA CFE astronauts (OSD-530; 6 astronauts, 11 timepoints). Human cfRNA NES was computed using fGSEA with gene rankings derived from `edge_pre_vs_flight_diff` (in-flight minus pre-flight normalized mean expression) against MSigDB Hallmark v7.5.1 human gene sets.

**Key result:** Spearman r = **0.352** (95% CI: 0.070–0.571), permutation p = **0.013**, n = 50 pathways.

The moderate but statistically significant cross-species conservation indicates that the direction of Hallmark pathway enrichment during spaceflight is partially shared between mouse liver bulk RNA-seq and human plasma cfRNA. Note that human cfRNA NES values were predominantly positive (46/50 pathways) due to the global increase in cfRNA abundance during spaceflight; the correlation captures relative differences in enrichment magnitude across pathways.

**Per-mission conservation:**

| Mouse mission | Spearman r | p-value | Duration |
|--------------|-----------|---------|----------|
| MHU-2 (Japan, ISS) | 0.007 | 0.962 | ~7 days |
| RR-1 (shuttle, 37d) | -0.096 | 0.507 | 37 days |
| RR-3 (shuttle, 16d) | 0.208 | 0.147 | 16 days |
| RR-6 (ISS, 35d) | **0.343** | **0.015** | 35 days |
| RR-8 (ISS, 35d) | 0.203 | 0.156 | 35 days |
| RR-9 (ISS, 33d) | **0.344** | **0.014** | 33 days |

RR-6 and RR-9 (both ISS missions, ~35 days) show the highest cross-species conservation (r ≈ 0.34, p < 0.02 each). MHU-2 and RR-1 show no significant conservation, suggesting mission-specific or hardware-dependent variability.

**Concordant pathways (|NES| > 0.5 in both species; upregulated in both):**

| Pathway | Human NES | Mouse NES |
|---------|-----------|-----------|
| HALLMARK_ANGIOGENESIS | +1.65 | +1.32 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | +1.44 | +1.41 |
| HALLMARK_COAGULATION | +1.61 | +0.97 |
| HALLMARK_ALLOGRAFT_REJECTION | +1.61 | +0.93 |
| HALLMARK_COMPLEMENT | +1.41 | +0.92 |
| HALLMARK_IL2_STAT5_SIGNALING | +1.21 | +0.90 |
| HALLMARK_APOPTOSIS | +1.25 | +0.82 |
| HALLMARK_INFLAMMATORY_RESPONSE | +1.22 | +0.78 |
| HALLMARK_INTERFERON_ALPHA_RESPONSE | +1.59 | +0.74 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | +1.50 | +0.72 |

**Discordant pathways (|NES| > 0.5, opposite signs):**

| Pathway | Human NES | Mouse NES | Interpretation |
|---------|-----------|-----------|----------------|
| HALLMARK_BILE_ACID_METABOLISM | -0.85 | +0.73 | Liver-specific (mouse) vs. blood (human) |
| HALLMARK_MITOTIC_SPINDLE | +1.34 | -0.90 | Cell cycle direction differs |
| HALLMARK_MYC_TARGETS_V1 | +0.75 | -0.69 | Biosynthesis direction differs |
| HALLMARK_ESTROGEN_RESPONSE_EARLY | +0.90 | -0.73 | Hormonal response differs |
| HALLMARK_HEME_METABOLISM | +0.82 | -0.57 | RBC turnover in cfRNA context |

**Interpretation:** Immune-related pathways (INFLAMMATORY_RESPONSE, INTERFERON responses, COMPLEMENT, ALLOGRAFT_REJECTION), tissue remodeling (EMT, ANGIOGENESIS), and coagulation are conserved across species and tissues, supporting a systemic spaceflight stress response. Discordant pathways reflect tissue-specific biology (liver bile acid metabolism) or opposite regulatory directions in blood vs. solid tissue contexts.

**Sign agreement:** 54% (27/50 pathways), marginally above 50% chance level. The combination of significant Spearman r (p = 0.013) and modest sign agreement (54%) suggests partial directionality conservation with substantial tissue-specific variation.

**Scripts and outputs:**
- R fGSEA: `v2/scripts/run_human_cfrna_fgsea.R`
- Python analysis: `v2/scripts/cross_species_nes_comparison.py`
- Results JSON: `v2/evaluation/E1_crossspecies_nes.json`
- Figure: `v2/figures/E1_crossspecies_scatter.html`

### E2: Mission Duration Drives Cross-Species NES Conservation ★★

We compared two human spaceflight cfRNA datasets with the mouse liver NES: JAXA CFE (OSD-530; 6 professional astronauts; 120-day ISS mission; n=26,845 genes) used in E1, and Inspiration4 (I4; 4 civilian crew; 3-day high-altitude mission; n=5,346 genes, log2FoldChange from DESeq2). Both datasets were compared against the same mouse liver mission-averaged NES (48 pathways in common with I4; 50 with JAXA).

| Mission | Duration | Spearman r | 95% CI | p (permutation) |
|---------|----------|-----------|--------|-----------------|
| I4 (Inspiration4) | 3 days | −0.095 | −0.379–0.211 | 0.516 (NS) |
| JAXA CFE | 120 days | +0.352 | +0.070–0.571 | 0.013 |
| **Δr (JAXA − I4)** | — | **+0.446** | — | — |

The absence of conservation in I4 (r ≈ 0, p = 0.52) contrasts sharply with the significant positive correlation in JAXA (r = 0.352, p = 0.013), yielding a Δr = +0.446. This large effect size suggests that mission duration — rather than civilian vs. professional crew composition — is the primary driver of cross-species NES conservation. The 3-day I4 mission likely induces transcriptional responses that are too transient and noisy to produce a coherent cross-species pathway signal, whereas the 120-day JAXA mission allows biologically meaningful pathway responses to develop.

**Caveat:** The I4 dataset has only 5,346 genes (subset common to I4 and Polaris Dawn cfRNA) vs. 26,845 for JAXA; reduced gene coverage may partially contribute to lower I4 NES reliability (only 48/50 Hallmark pathways computed with minSize=10).

**Key results:**
- I4 (3d): r = −0.095 (NS) — no cross-species conservation
- JAXA (120d): r = +0.352 (p = 0.013) — significant conservation
- Δr = +0.446 — large duration-dependent effect

**Scripts and outputs:**
- R fGSEA: `v2/scripts/run_i4_cfrna_fgsea.R`
- Python analysis: `v2/scripts/mission_conservation_comparison.py`
- Results JSON: `v2/evaluation/E2_mission_conservation.json`
- Figure: `v2/figures/E2_duration_conservation.html`

---

## Section 3: Cell-Type-Resolved Spaceflight Response (Category F)

### F1: I4 PBMC Cell-Type-Specific Pathway Response (snRNA-seq) ★★

**Data:** GLDS-562 (OSD-570), Inspiration4 PBMC snRNA-seq (10x Multiome), SpaceOmicsBench v2_public
**Comparison:** I4-FP1: R+1 vs. pre-flight (L-92, L-44, L-3); Seurat FindMarkers per cell type
**Pipeline:** Seurat FindMarkers avg_log2FC → fGSEA Hallmark per cell type (minSize=10)
**Cell types:** CD14+ Monocyte, CD16+ Monocyte, Dendritic Cell, Natural Killer Cell,
              B Cell, CD4+ T Cell, CD8+ T Cell, Other T Cell, Other, PBMC Pseudobulk (n=10)
**Pathways computed:** 13–30 per cell type (mediated by FindMarkers gene count; 34 total unique)
**Output:** `v2/processed/F1_scrna/i4_snrnaseq_celltype_fgsea.csv` (200 rows)
**Figure:** `v2/figures/F1_celltype_pathway_heatmap.html`

#### Key Results

**Universal response — all cell types:**

| Pathway | NES range | n cell types sig | Interpretation |
|---------|-----------|-----------------|----------------|
| MYC_TARGETS_V1 | −2.14 to −2.90 | 10/10 | Ribosome biogenesis suppressed in all PBMCs at R+1 |
| OXIDATIVE_PHOSPHORYLATION | −1.99 to −2.21 | 5/10 | Mitochondrial metabolism suppressed |

**Cell-type-specific responses (padj<0.05 per cell type):**

| Cell Type | Sig pathways | Top pathway | NES | Interpretation |
|-----------|-------------|-------------|-----|----------------|
| Dendritic Cell | **5** | MYC_TARGETS_V1 | −2.14 | Strongest innate immune response |
| Natural Killer Cell | **4** | MYC_TARGETS_V1 / TNFA_NFKB | −2.51 / −2.02 | NK cytotoxicity + inflammation |
| CD14+ Monocyte | **3** | MYC_TARGETS_V1 + P53_PATHWAY | −2.86 / −2.07 | Classical monocyte stress |
| CD16+ Monocyte | **3** | MYC_TARGETS_V1 + OXPHOS | −2.27 / −2.11 | Non-classical monocyte metabolism |
| CD4+ T Cell | **1** | MYC_TARGETS_V1 | −2.21 | Adaptive immune cells less responsive |

**Innate > Adaptive:** Innate immune cells (DC+NK+Mono) show 3–5 sig pathways vs. T cells (1–2).

**PBMC Pseudobulk vs. cell type:** MYC_TARGETS_V1 (NES=−2.44) and OXPHOS (−1.99) — consistent with dominant monocyte/NK signal driving bulk result.

#### Connection to E1/E2 (cfRNA)
The I4 cfRNA analysis (E2) used bulk plasma cfRNA (R+1: NES = −0.095 vs. mouse, NS). The cell-type-resolved snRNA-seq reveals that MYC_TARGETS_V1 suppression is universally present across PBMCs at R+1 — suggesting cfRNA dilution of cell-type-specific signals may explain the weak cross-species correlation for short missions.

#### F1-Temporal: Delayed and Persistent PBMC Response ★★

**Comparison across 5 timepoints:**

| Timepoint | Comparison | Total sig (sum across cell types) | Interpretation |
|-----------|-----------|-----------------------------------|----------------|
| FP1 | R+1 vs. pre | 26 | Immediate spaceflight response |
| LP1 | R+45 vs. pre | **45** | Peak divergence from baseline |
| LP3 | R+45+R+82 vs. pre | 13 | Partial normalization by R+82 |
| RP1 | R+45 vs. R+1 | **63** | Peak recovery dynamics |
| RP3 | R+45+R+82 vs. R+1 | 21 | Diminishing recovery signal |

**Key finding: LP1 > FP1** — At R+45 (45 days after landing), PBMCs are MORE transcriptionally distinct from pre-flight than at R+1. The spaceflight-induced state persists/amplifies during early return, paralleling the mouse liver T2 overshoot (93% pathway recovery with mean 1.75× overshoot fraction).

**RP1 peak**: The R+45 vs. R+1 comparison (recovery phase) shows 63 total significant pathway-cell type combinations — the highest of all timepoints — indicating active transcriptional reorganization during post-flight recovery.

Cell types most affected at R+45 (LP1):
- Other T Cell: 9 sig pathways (vs. 2 at FP1)
- Dendritic Cell: 8 (vs. 5)
- CD16+ Monocyte: 7 (vs. 3)
- B Cell: 6 (vs. 1)

**Cross-species parallel**: The human PBMC temporal pattern mirrors mouse liver T2 — peak transcriptional response occurs not at flight landing (R+1 ≡ ISS-T) but during post-return recovery (R+45 ≡ LAR with overshoot).

**Scripts and outputs:**
- R fGSEA: `v2/scripts/run_i4_snrnaseq_temporal_fgsea.R`
- Python heatmap: `v2/scripts/i4_temporal_heatmap.py`
- Data: `v2/processed/F1_scrna/i4_snrnaseq_temporal_fgsea.csv` (834 rows, 5 sheets × 10 cell types)
- Figure: `v2/figures/F1_temporal_heatmap.html`

#### RRRM-1 scRNA-seq Pipeline (benchmark-ready subset)
OSD-904~934 raw FASTQ availability was confirmed on Cayuga, but the benchmark-focused subset was narrowed to four tissues:
- **Selected 4 tissues** matching benchmark: OSD-918 (blood), OSD-920 (eye), OSD-924 (muscle), OSD-934 (skin)
- Total selected raw footprint: ~308 GB, 8 libraries per OSD, 10x Chromium 3' v3
- STARsolo recovery completed on Cayuga with STAR 2.7.3a-compatible settings
- Final outputs now available:
  - per-sample `Solo.out/GeneFull/filtered`
  - per-sample `.h5ad`
  - merged `.h5ad`
  - tissue-specific processed objects
  - first-pass broad cell-type annotation tables
- Broad composition highlights:
  - blood is erythroid-heavy
  - eye is retinal neuronal / Muller glia dominant
  - muscle is immune-rich
  - skin contains strong keratinocyte and immune-myeloid compartments

---

### E3: PBMC Cell-Type Transcriptome is Anti-Correlated with Plasma cfRNA at Flight Onset ★★★

**Question:** Which PBMC cell type's pathway response best predicts the plasma cfRNA NES pattern at R+1?

**Data:** E2 cfRNA NES (48 pathways) × F1 snRNA-seq cell-type NES (10 types × 34 common pathways)
**Method:** Spearman r per cell type, permutation p (n=10,000), bootstrap 95% CI (n=1,000)
**Output:** `v2/evaluation/E3_cfrna_origin.json`, `v2/figures/E3_cfrna_origin.html`

#### Key Results

| Cell Type | Spearman r | 95% CI | p | n pathways |
|-----------|-----------|--------|---|------------|
| Natural Killer Cell | −0.505 | — | 0.028* | 19 |
| CD4+ T Cell | −0.505 | — | 0.081 | 13 |
| PBMC Pseudobulk | −0.486 | — | 0.034* | 19 |
| CD8+ T Cell | −0.494 | — | 0.057 | 16 |
| CD14+ Monocyte | −0.475 | — | 0.056 | 17 |
| Other T Cell | −0.421 | — | 0.042* | 24 |
| Dendritic Cell | −0.395 | — | 0.038* | 28 |

**Critical finding: All cell-type NES values are NEGATIVELY correlated with cfRNA NES.**

The PBMC intracellular transcriptome response at R+1 is systematically anti-correlated with plasma cfRNA pathway signals. The most striking example: **MYC_TARGETS_V1** is strongly downregulated in all PBMC cell types (NES −2.1 to −2.9) but is among the most upregulated pathways in cfRNA (NES +2.0, MYC_TARGETS_V2). This inversion suggests that PBMC cells actively suppress MYC-driven biosynthesis at flight onset, while simultaneously releasing existing MYC-target transcripts into the plasma (via apoptotic bodies, extracellular vesicles, or stress-induced mRNA secretion).

#### Mechanistic Interpretation

The E1/E2/E3 results together form a coherent narrative:
1. **E2**: I4 cfRNA shows no conservation with mouse liver at R+1 (r = −0.095, NS)
2. **E3**: I4 PBMC transcriptome is anti-correlated with cfRNA at R+1 (r ≈ −0.5)
3. **Together**: At 3-day flight onset, plasma cfRNA signals are neither mouse-conserved (E2) nor PBMC-derived (E3) — suggesting cfRNA at R+1 represents **passive release of pre-existing cellular RNAs from stressed/apoptosing cells**, rather than active transcriptional programs observed in bulk tissues.

This contrasts with JAXA (120-day, E1 r = +0.352): longer duration allows genuine tissue-specific pathway changes to accumulate and be reflected in plasma cfRNA.

**Scripts and outputs:**
- `v2/scripts/e3_cfrna_celltype_origin.py`
- `v2/evaluation/E3_cfrna_origin.json`
- `v2/figures/E3_cfrna_origin.html`

---

## Methods Notes

### T1 Method Detail
- Classifier: Logistic Regression (PCA-reduced, 50 components)
- CV: RepeatedStratifiedKFold (5-fold × 10 repeats)
- AUROC 95% CI: bootstrap (1000 iterations)
- Permutation test: 10,000 label shuffles

### T2 Method Detail
- PCA: StandardScaler → PCA(n_components=3) on combined ISS-T + LAR samples
- Recovery ratio: d(FLT_ISS-T → FLT_LAR) / d(BSL_ISS-T → FLT_ISS-T), Euclidean distance in PC1-3 space
- Pathway NES: fGSEA Hallmark v7.5, pre-ranked by condition mean vs. BSL_ISS-T
- Recovery threshold: |Δ_flight| ≥ 0.10 NES for inclusion in recovery analysis
- Recovery fraction = Δ_return / |Δ_flight|; fraction > 1 = overshoot; fraction 0–1 = partial recovery; fraction < 0 = deepening

### T3 Method Detail
- Age groups: OLD = 32 week (n=17 FLT, n=17 GC, n=16 VIV); YNG = 10 to 12 week (n=18 each)
- Classifier: Logistic Regression, RepeatedStratifiedKFold (RSKF)
- ANOVA: Two-way (flight × age) per pathway, BH-FDR correction, α = 0.05
- Timing filter for T3c: ISS-T only (to exclude preservation-confounded LAR)

### E1 Method Detail
- Human data: JAXA CFE cfRNA (OSD-530), n=26,845 genes ranked by `edge_pre_vs_flight_diff`
- Mouse data: 6 missions × Hallmark NES, arithmetic mission average
- fGSEA: R `fgsea` package, Hallmark v7.5.1 GMT (Homo sapiens), minSize=15, nPermSimple=10,000
- Correlation: Spearman r, bootstrap 95% CI (n=1,000), permutation p (n=10,000)

### F1 Method Detail
- Data: GLDS-562 snRNA-seq, sheet I4-FP1, header row 7, `avg_log2FC` column as ranking statistic
- FindMarkers cutoff: `|avg_log2FC| ≥ 0.25` (Seurat default) → 615–1,070 genes per cell type
- fGSEA: R `fgsea` package, Hallmark v7.5.1 (Homo sapiens GMT), minSize=10, maxSize=500, nPermSimple=10,000, seed=42
- Deduplication: if gene appears twice, keep row with highest |avg_log2FC|
- Significance threshold: padj < 0.05 (Benjamini-Hochberg)

### E3 Method Detail
- cfRNA NES: from E2 output (I4 fGSEA, 48 pathways)
- Cell-type NES: from F1 output (I4 snRNA-seq FP1, per cell type)
- Common pathways: 34 (intersection of 48 cfRNA and 34 unique snRNA pathways)
- Effective n per cell type: 13–30 (varies due to FindMarkers gene count filter)
- Spearman r: scipy.stats.spearmanr on paired (non-NaN) pathway NES values
- Bootstrap 95% CI: n_boot=1,000, random resampling with replacement
- Permutation p: n_perm=10,000 label shuffles
- Significance: p < 0.05 (permutation)

### E2 Method Detail
- I4 data: Inspiration4 cfRNA DESeq2 results (log2FoldChange_I4), 5,346 genes (common with Polaris Dawn)
- Data source: SpaceOmicsBench `2025_01_08_v2.cfRNA/SpaceOmicsBench_v2.0/data/cfrna_crossmission_r1.csv`
- fGSEA: same GMT as E1, minSize=10 (to accommodate reduced gene count), 48/50 pathways computed
- Comparison: Spearman r between I4 NES and mouse liver NES (48 common pathways)
- JAXA result from E1 (50 pathways) compared with I4 result (48 pathways); common-pathway comparison used for Δr

---

## Supplementary Notes

### Tier 3 LLM parsing metrics
For the Tier 3 zero-shot LLM benchmark, the main text should report only the end-to-end AUROC because that is the deployed benchmark output submitted for evaluation. Parsed-only AUROC, provider-specific parse rates, and truncation-related parser failure analysis should be placed in the Supplementary Information or reviewer response.

Recommended table location: `v2/docs/V2_SUPPLEMENT_TABLES.md` as Supplementary Table S1, with end-to-end AUROC retained in the main text and parse-aware decomposition moved entirely to Supplementary Information.

Current supplement-level summary:
- DeepSeek-V3: 549/549 parsed, mean end-to-end AUROC = 0.507, mean parsed-only AUROC = 0.507
- Gemini-2.5-Flash: 529/549 parsed, mean end-to-end AUROC = 0.493, mean parsed-only AUROC = 0.493
- Llama-3.3-70B (Groq): 489/549 parsed, mean end-to-end AUROC = 0.485, mean parsed-only AUROC = 0.443

Interpretation:
- DeepSeek's end-to-end and parsed-only scores are identical because parse rate is 100%.
- Gemini improves after conservative parser recovery of truncated narrative verdicts, but provider comparison still depends partly on output-format robustness.
- Groq remains the clearest case where parsed-only and end-to-end performance should be shown together because truncation removes explicit answers in a substantial fraction of samples.

### Why T2 RR-8 is the primary recovery analysis
RR-8 is preferred over RR-6 for T2 because (1) larger cohort (ISS-T n=20, LAR n=15 per condition vs. n=10 each in RR-6), (2) longer mission duration allowing stronger in-flight perturbation to develop, and (3) the aging cohort design provides additional biological context. RR-6 is shown as comparison in Extended Data.

### T3 biological interpretation
The +0.266 AUROC delta between OLD and YNG mice is consistent with known age-related decline in homeostatic resilience. Aged animals may have reduced capacity to buffer transcriptomic responses to environmental stress, resulting in larger and more persistent gene expression changes. The absence of significant pathway-specific interaction terms (T3c) suggests that aging does not redirect the spaceflight response to novel pathways but rather amplifies the same pathways that respond in young animals.

---

*Next update: Section 2 E1/E2 results after Phase 2 execution.*
*Target journal: Cell Systems / npj Systems Biology.*

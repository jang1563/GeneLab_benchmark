# SpaceMed: A Causal Analysis Framework for Spaceflight Biology and Mars-Risk Hypotheses

## Integrated Nature Manuscript — One Coherent Narrative

### Opening Hook (Main Text, ~200 words)
**The Mars Problem Unmasks ISS Confounding**
- Artemis aims for Moon; Mars missions = 900 days
- ISS missions confound three stressors: µg + low GCR + isolation
- Mars trajectory: lower µg (0.38g on surface) BUT 5× higher GCR + extended isolation
- **Central question:** Where do Mars-like stressor regimes exceed what ISS-trained linear extrapolation can support?
- This is why NASA TRISH needs stressor-decomposed evidence before operational claims
- **Our solution:** SpaceMed — causal analysis framework linking mouse OSDR,
  human Inspiration4, and ground analog datasets to answer 3 questions:
  1. Which mouse pathway signals transfer to human spaceflight data? (Species transfer -> AUROC 0.888)
  2. Where do Mars-like stressor regimes break linear extrapolation? (Decompose ISS confounders -> factorial model -> uncertainty)
  3. Which perturbations reverse tissue signatures enough to justify
     pre-clinical validation? (LINCS + CRISPR orthogonal support)

### Results (Main Text, ~2000 words total for Nature)

#### Pillar 1: Species Transfer Supports Mouse Pathways as Human-Relevant Features
**Setting:** Inspiration4 crew (n=4 × 10 timepoints) + 6 mouse OSDR tissues (n=57 missions, 8 tissues)

**Key finding:** Mouse tissue-level pathway NES (Normalized Enrichment Score) improves supervised classification of human post-flight vs pre-flight compartment conservation compared with SpaceOmics features alone.
- Unsupervised: I4 PBMC mean NES Spearman r = 0.777 across tissues
- Supervised (RF on 8 SpaceOmics features + 6 mouse tissue NES): **AUROC 0.888 [0.854–0.918]** vs baseline 0.712
- **+0.175 AUROC lift** (95% CI [0.134, 0.219]) after the HPC full-MSigDB rerun

**Implication:** Despite gene-level sequence divergence, **pathway-level conservation is informative**. Mouse OSDR data can support human-spaceflight hypothesis generation, but it does not by itself validate Mars risk forecasting.

**Mechanistic insight:** Pathway conservation correlates with tissue function (immune genes thymus ↔ PBMC; metabolic genes liver ↔ plasma lipids). Species transfer strength depends on tissue homology, not just taxonomy.

#### Pillar 2: Factorial Analogs Decompose Confounded Spaceflight Signals
**Setting:** 4 ground-based 2×2 or 2×2×2 factorial analog studies (OSD-211 spleen, OSD-237 skin, OSD-202 brain, OSD-719 endocrine).

**Key finding #1 — Interactions Dominate:** Linear model: log2(CPM+1) ~ HLU + IR + HLU:IR + [Time]
- Interaction terms explain **44–61% of variance** in top-responsive genes
- Example: thymus 61% (HLU×IR dominates), skin 48%, brain 52%
- **NOT additive.** Synergy between stressors is mechanistically fundamental.
- ICP causal inference confirms interactions most stable across environments (ICP=0.52–0.54 vs main effects 0.40–0.41)

**Key finding #2 — Radiation Quality Flips Sign:**
- Low-LET γ-rays (OSD-211, 0.04 Gy ⁵⁷Co): r_flight = +0.36 (thymus signature matches flight)
- High-LET HZE mixed field (OSD-719, NSRL ⁵⁶Fe + ²⁸Si): r_flight = **–0.22** (opposite sign!)
- **Critical implication:** Dose equivalence (Gray) does not imply biological equivalence (LET-dependent). Mars GCR (high-LET) may produce transcriptomic responses that differ from ISS low-LET radiation and therefore needs separate validation.

#### Pillar 3: LINCS + CRISPR Prioritize Multi-Tissue Perturbation Hypotheses
**Setting:** L1000CDS2 chemical perturbation library (reverse tissue-specific flight signatures) + Enrichr CRISPR KO library (orthogonal support).

**Key finding #1 — Pareto-Optimal Multi-Tissue Signature Reversers:**
- 26 drugs reverse ≥2 tissues
- Pareto front (maximize min reversal, minimize signature aggravation): **2 lead hypotheses**
- Current Pareto front: CGP-60474 (3 tissues, mean_reversal=0.1457, min=0.0735) + quinacrine hydrochloride (2 tissues, mean_reversal=0.1364, min=0.0909)
- Dorsomorphin remains a gastrocnemius/liver tissue-level AMPK/BMP-axis hit, but not a current Pareto-front lead

**Key finding #2 — Chemical-CRISPR Convergence:**
- CDK inhibitors (AZD-5438, AT-7519) share CDK9 as one target
- CRISPR library: CDK9 KO reverses kidney UP genes, adj-p=0.1588, 7 overlapping genes
- **Target plausibility is suggestive, not validated.** CRISPR-library convergence
  strengthens the hypothesis but is not FDR-significant and does not establish
  safety or efficacy.

**Key finding #3 — Tissue-Specific Novel Targets:**
- CRISPR-only hits (no chemical match): FOXO3 (liver), IL21R/CXCL13 (thymus)
- These represent exploratory target axes for follow-up perturbation studies

#### Pillar Integration: Mars Extrapolation with CIs
**Setting:** Apply factorial coefficients to Mars mission stressor doses (µg=0.62×, HZE=12.9×, Time=7.5× vs ISS analog).

**Key finding #1 — Linear Extrapolation Breaks:**
- Exploratory Mars response: ΔExpr = β_stressor × dose_stressor
- Spearman r(projected, flight-observed) ≈ 0 across tissues — **not a point predictor**
- Interpretation: Interaction terms non-linear at >5× dose amplification. Gene response saturates or recruits additional pathways not captured by factorial model.
- **Conservative stance:** Mars projections flag regime change, not point estimates. Mechanistic constraints needed.

**Key finding #2 — Top Amplified Genes (Bootstrap CIs):**
- WNT10B (brain-to-eye analog): ~1052× amplification flag, with wide bootstrap uncertainty
- KRTAP19-2 (skin analog): ~414× amplification flag, with wide bootstrap uncertainty
- These genes are **extremely sensitive in the linear stress test** and require non-linear calibration before biological promotion.

**Key finding #3 — Causal Ranking (ICP):**
- Stressor effect sizes ranked by causal stability:
  - Time: ICP=0.540 (most causally invariant — mission duration is primary driver)
  - Interactions HLUxIRxT: ICP=0.522
  - Main effects HLU, IR: ICP=0.40 (least stable, most context-dependent)

### Discussion

#### What This Reveals About Spaceflight Biology
1. **Pathway vs. Gene Divergence:** Species transfer works at pathway level even though individual genes diverge. Implication: use functional readouts (pathway scores) not raw gene counts for predictive models.

2. **Synergy, Not Additivity:** Interaction dominance (44–61%) suggests spaceflight stressors activate qualitatively different cellular states when combined. Single-stressor studies (e.g., HLU alone) are incomplete for Mars hypothesis testing.

3. **Radiation Quality Matters Fundamentally:** High-LET GCR on Mars should not be assumed to follow ISS dose-response. Separate countermeasure-hypothesis testing is needed for Mars vs ISS.

4. **Time is Most Stable in the ICP Ranking:** ICP=0.540 — mission duration is more stable across available environments than individual dose metrics. Implication: 900-day Mars hypotheses should explicitly model temporal degeneration, not only acute stress response.

#### Pre-Clinical Hypotheses for Artemis/Mars
1. **Immune Reconstitution Axis:** Thymus signatures nominate CDK/lymphoid
   signaling as a validation target, but no crew protocol should be implied from
   transcriptomic reversal alone.

2. **Muscle Preservation Axis:** Gastrocnemius signatures nominate AMPK and
   related metabolic pathways for rodent testing alongside exercise or loading
   interventions.

3. **GCR-Specific Discovery:** Mars GCR differs from ISS radiation enough that
   radiation countermeasure hypotheses should be tested separately
   instead of extrapolated directly from ISS-like exposure.

#### Limitations
- **SOMA n=4:** Existence proof; larger cohorts (TRISH-funded missions 2025–2027) will validate
- **Linear Mars model breaks:** Indicates mechanistic complexity not fully captured. Future: dose-response threshold models + network biology constraints
- **ICP identifiability:** Two environments per stressor minimum; some combinations unidentifiable. Future: more ground analogs (e.g., Antarctic isolation data)
- **Chemical drug coverage:** LINCS ~10,000 compounds; natural products, biologics, gene therapy missing. Future: multi-modality screening

#### Future Directions: The Causal DAG
SpaceMed's long-term strength is a causal graph that can eventually support
counterfactual queries such as "if an intervention targets a nominated pathway,
what pathway response is projected under a Mars-like stressor regime?" This
requires:
- Ranking candidate edges via Invariant Causal Prediction-style stability scoring (current ICP time=0.540)
- Assigning quantitative causal effects (from factorial model βs)
- Validating counterfactual interventions in rodent RCTs before operational use

### Methods (Brief for Nature; Extended Methods online)

**Pillar 1 — Species Transfer:**
- fGSEA pathway NES (Hallmark, KEGG, Reactome, MitoCarta, C2CGP, C5BP: 8,428 unique mouse pathways in the full-MSigDB bridge)
- Spearman correlation (mouse NES vs human compartment NES)
- Supervised RF/LR on 14 features (8 SpaceOmics aggregate + 6 mouse tissue NES), 5-fold CV + 1000 bootstrap

**Pillar 2 — Factorial Decomposition:**
- OLS per gene: log2(CPM+1) ~ design matrix (2×2 or 2×2×2)
- ICP across 9 environments to rank causal stressors
- Comparison of low-LET (⁵⁷Co γ) vs high-LET (HZE ⁵⁶Fe + ²⁸Si) via Spearman r to flight

**Pillar 3 — Countermeasure-Hypothesis Generation:**
- L1000CDS2 API query with tissue up/dn signatures (top 150 genes each)
- Enrichr LINCS_L1000_CRISPR_KO library for target support
- Pareto filtering: ≥2 tissues reversed, dominant on efficiency frontier

**Mars Extrapolation:**
- Linear projection of factorial β × Mars dose (µg=0.62, HZE=12.9, T=7.5)
- Bootstrap CIs from 1000 resamples of β ± SE

### Figures (Main Paper: 4 figures; Extended Data: 8 figures)

**Figure 1: Species Transfer (BRIDGE)**
- A: NES heatmap (mouse tissue × I4 compartment)
- B: Supervised AUROC improvement (baseline 0.712 → 0.888)
- C: Feature importance / ablation

**Figure 2: Stressor Decomposition (DECOMPOSE)**
- A: Interaction dominance (variance fraction per tissue)
- B: ICP causal stability heatmap
- C: Radiation quality (low-LET vs high-LET opposite signs)

**Figure 3: Mars Extrapolation + Perturbation Hypotheses**
- A: Dose amplification factors (1.0× ISS → 0.62–12.9× Mars)
- B: Mars-projected vs flight-observed gene response (scatter, top movers labeled)
- C: Pareto front (drug × tissue reversal scores)

**Figure 4: Causal Integration & Future Interventional Model**
- Causal DAG: Stressors → Tissue responses → Human outcomes → Interventions
- Annotated with ICP scores + key species-transfer pathways
- Future counterfactual query example: "pathway intervention → projected tissue
  response under Mars-like stressor regime"

**Extended Data Figures:**
- ED1: Full tissue NES conservation (all mouse × human pairs)
- ED2–4: Factorial models per tissue (volcano plots, interaction heatmaps)
- ED5: CRISPR orthogonal-support matrix (chemical targets × KO hits)
- ED6: Multi-tissue drug matrix (26 multi-tissue drugs, Pareto filtering)
- ED7: Mars bootstrap CIs (top 50 genes with 95% credible intervals)
- ED8: ICP stability across stressor pairs and tissues

### Supplementary Tables
- **ST1:** Species transfer NES correlation matrix (all tissues × all I4 compartments)
- **ST2:** Factorial coefficients (all genes, all analogs, with SE and p-values)
- **ST3:** Mars extrapolation predictions (all genes, all analogs, with bootstrap CIs)
- **ST4:** CRISPR orthogonal support (chemical target genes × KO library matches)
- **ST5:** L1000CDS2 reversal scores (all 26 multi-tissue drugs, all tissues)
- **ST6:** Results summary table (consolidated metrics per tissue)

### Abstract (250 words)
[Narrative summarizing all three pillars; see v8/MANUSCRIPT_DRAFT.md]

---

## Publication Strategy

**Venue:** *Nature* (Integrated)
- High impact, accepts systems-level cross-species+multi-stressor papers
- Comparable prior: Overbey et al. 2024 *Nat* (Inspiration4 multi-omics); Kiffer et al. 2023 *Nat Commun* (NSRL GCR)
- **Positioning:** "Causal, multi-species spaceflight biology framework" —
  unifying previous separate efforts (v1–7 ML, v5 drug targets, SOMA, GCRsim,
  HLU analogs) while keeping clinical claims exploratory until validation.

**Timeline:**
- Figures & revisions: 2 weeks (bootstrap CIs finish this week)
- Internal review: 1 week
- Submission: Early May 2026
- Review cycle: 8–12 weeks (fast-track if Nature recognizes novelty)
- Revised manuscript: July 2026

**Data Availability:**
- OSDR: Public (https://genelab.nasa.gov)
- SOMA: Public (Overbey et al. 2024 supplementary)
- Code: GitHub repo (public at acceptance)
- Results tables: Zenodo (DOI at submission)

**Reviewers to Target:**
- Space biology experts: Duane Graveline (NASA), Thoru Pederson (multi-omics)
- Causal inference: Jonas Peters (Invariant Causal Prediction developer), Elias Bareinboim (do-calculus)
- Drug repurposing: Avi Ma'ayan (LINCS/Enrichr), Chris Tansey (CRISPR screens)
- Radiation biology: David Brenner, Susan Bailey (HZE vs gamma expertise)

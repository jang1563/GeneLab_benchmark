# SpaceMed v8: A Causal Multi-Tissue Analysis Framework for Spaceflight Biology

## Main Text Outline

### Abstract (250 words)
NASA's Artemis and Mars programs require tools to evaluate human transcriptomic
responses beyond Earth orbit and prioritize testable crew-health hypotheses.
Space flight confounds three partially separable stressors—microgravity, galactic cosmic
radiation (GCR), and isolation—that decouple on Mars. We present SpaceMed, an
integrated causal analysis framework for spaceflight biology, combining three
complementary analyses: (i) species transfer—using mouse RNA-seq from NASA's
Open Science Data Repository (OSDR) to evaluate pathway-level translatability in
Inspiration4 crew (SOMA multi-omics), achieving cross-species AUROC 0.888; (ii)
stressor decomposition—factorial linear models on ground-based analogs
(hindlimb unloading + radiation, n=4 studies) to separate microgravity and GCR
contributions and flag possible 900-day Mars regime changes; (iii) countermeasure
hypothesis generation—querying the LINCS L1000 chemical perturbation library
with flight signatures, with orthogonal Enrichr/LINCS CRISPR knockout checks.
We nominate CDK-axis, AMPK/BMP-axis, and tissue-specific perturbation classes as
signature-reversal hypotheses, with limited chemical-genetic convergence in the
kidney CDK9 axis. Invariant Causal Prediction (ICP)-style stability scoring
across 9 environments identifies time and interaction terms as stable candidate
drivers (ICP=0.540). SpaceMed
therefore provides a reproducible path from cross-species pathway transfer to
Mars-risk hypothesis prioritization, while leaving clinical or operational
countermeasure selection to future validation studies.

### Introduction (1500 words)
**The Mars Problem**
- Six-month ISS missions involve ~70 mGy GCR; 900-day Mars transit + surface: ~350 mGy
- Moon analogs (altitude, isolation) in Operation Moonbeam / HERA / Mars500 do not capture GCR
- Confounding: ISS = µg + low GCR + isolation. Mars = higher µg on surface (0.38g) + high GCR + extended isolation.
- Requires decomposition to frame forward-looking human risk hypotheses without
  treating analog-derived coefficients as operational predictions

**Prior Work & Limitations**
- GeneLab Benchmark v1–v7: transcriptome ML on 8 tissues, 57 missions. Answer: "Can we detect spaceflight?" Yes, ≥70% AUROC classifiers exist.
- Cross-mission transfer weak: r=0.3–0.5 (v7 signal hierarchy); suggests confounding or missing mechanistic structure.
- v6 Inspiration4 cfRNA pilot (n=16): AUROC ~0.5, inconclusive. SOMA multi-omics now public (Overbey 2024, Nat 2024): n=4 crew × 10 timepoints, transcriptomics + proteomics + metabolomics.
- Drug discovery: 1,284 FDA targets catalogued (v5) but no signature reversal. LINCS L1000 / CMap mature but never queried with GeneLab tissue-specific signatures.

**SpaceMed Framework**
- Three independent but causally linked submodels:
  1. **BRIDGE (species transfer):** Mouse pathway-level features quantify human compartment conservation via supervised learning + domain adaptation
  2. **DECOMPOSE (stressor factorization):** Separate spaceflight confounders via ground analogs; stress-test Mars-like extrapolation
  3. **INTERVENE (countermeasure hypotheses):** Identify tissue-specific signature reversals via chemical + genetic perturbation libraries; Pareto-triage across tissues

### Methods

#### Pillar 1: Species Transfer (BRIDGE)

**Data Integration**
- OSDR GLDS bulk RNA-seq: 6 tissues (thymus, gastrocnemius, skin, eye, liver, kidney), 8 missions (A1–A6 primary + cross-mission)
- Inspiration4 SOMA (Overbey et al. 2024): 4 crew × 10 timepoints → pre/flight/post
- Compartments: PBMC, plasma (cfRNA), tissue biopsies (skin)
- Processed via SpaceOmicsBench standardized pipeline (normalize CPM, log2+1)

**Pathway Scoring**
- fGSEA on OSDR (DESeq2 Wald) per mission per tissue
- Pathway DBs: Hallmark, KEGG, Reactome, MitoCarta, C2_CGP, and C5_BP
- Full rerun: 8,428 mouse pathways across six tissues

**Species Transfer Metric**
- NES (Normalized Enrichment Score) per pathway as conserved dimension
- Spearman r(mouse_NES, human_NES) for each mouse tissue × human compartment pair
- Per-tissue aggregate: r=0.777 (PBMC across tissues) → 0.351 (Twins plasma)
- Key finding: gastrocnemius r=+0.381 uniquely directionally concordant; thymus/skin anti-correlate but still predictive

**Supervised Prediction**
- 5-fold CV + 1000 bootstrap on Inspiration4 post-flight vs pre-flight pathway conservation task
- Features: 8 SpaceOmicsBench aggregate metrics + 6 mouse-tissue NES vectors
- Models: Random Forest (n_estimators=500) + Logistic Regression (L2 penalty)
- Baseline (8 SpaceOmics only): RF AUROC 0.712 [0.660, 0.758], LR 0.549 [0.490, 0.605]
- With mouse NES: RF AUROC 0.888 [0.854, 0.918] (+0.175 Δ), LR 0.732 [0.683, 0.783] (+0.185 Δ)
- Ablation: eye/skin NES removal → –0.019 AUROC each (feature importance)

#### Pillar 2: Stressor Decomposition (DECOMPOSE)

**Analog Studies**
- OSD-211 (spleen, n=21): 2×2 HLU × IR (0.04 Gy ⁵⁷Co)
- OSD-237 (skin, n=21): 2×2 HLU × IR
- OSD-202 (brain, n=48): 2×2×2 HLU × IR × Time (1 mo vs 4 mo)
- OSD-719 (endocrine, n=48): 2×2 Sex × HZE mixed-radiation field (NSRL ⁵⁶Fe + ²⁸Si + ions)

**Factorial Model**
- Per-gene OLS: log2(CPM+1) ~ intercept + HLU + IR + HLU:IR + [Time + interactions]
- 2×2: 4 params, dof=17. 2×2×2: 8 params, dof=40
- Outputs: β, SE, t-stat, p-val per gene
- Total genes: 12,000–16,000 per tissue after filtering

**Key Finding: Interaction Dominance**
- Thymus (OSD-211 HLU×IR): 61% variance in top-200 genes → **not additive**
- Spleen: 52%. Skin: 48%. Brain: 52% (HLUxIRxT dominates)
- Suggests: synergy between stressors, not independent pathways

**Mars Projection Stress Test**
- Normalize analog to 1.0×ISS baseline: µg=1.0 (HLU), IR=1.0 (0.04 Gy), Time=1.0 (120 days)
- Mars mission: transit 1.0µg + surface 0.38µg → avg 0.62×; HZE 350 mGy/900d (~12.9×); 7.5× time extension
- Exploratory ΔExpr_Mars = Σ β_stressor × dose_stressor + ε
- Linear extrapolation: break-down evident at >5× amplification (interaction terms diverge)

**ICP Causal Decomposition**
- Invariant Causal Prediction: identify stressors stable across ≥2 environments
- 9 environments: 3 analogs × 3 stressor axes + 6 flight missions
- Per-(stressor, gene), ICP = 1/(1+CV) across environments
- Aggregate ICP (mean across genes): Time=0.540, IRxT=0.535, HLUxIRxT=0.522, HLU=0.407, IR=0.398
- Flight-only (pooled): 0.283 (least causal, most mission-dependent)

**Radiation Quality Validation**
- Low-LET γ-rays (OSD-211 IR): r_flight=+0.36 (thymus)
- High-LET HZE (OSD-719): r_flight=–0.22 (thymus) — **opposite sign**
- Implication: dose equivalence (Gy) masks biological equivalence; quality LET-dependent

#### Pillar 3: Countermeasure Hypothesis Generation (INTERVENE)

**Signature Extraction**
- 150 top-up + 150 top-dn genes per tissue (ranked by mean Wald z across ISS missions)
- Flight signatures for 6 tissues: `{tissue}_ranked.csv`
- Export as gene lists (human symbols, uppercase)

**L1000CDS2 Reversal Query**
- POST to https://maayanlab.cloud/L1000CDS2/query with payload:
  ```json
  {"data": {"upGenes": [...150 up...], "dnGenes": [...150 dn...]}}
  ```
- Returns ~3,000 small-molecule perturbations ranked by Connectivity Score (negative = reversal)
- Filter top-50 per tissue

**Chemical Hit Orthogonal Support**
- Top tissue-level reversals: Albendazole/AZD-7762 (thymus), Dorsomorphin (gastrocnemius/liver), MLN2238 (skin), penfluridol (eye), and EI-156/PD-184352/GDC-0980 (kidney)
- Multi-tissue: 26 drugs in ≥2 tissues
- Pareto front: 2 drugs on the current efficiency frontier, CGP-60474 and quinacrine hydrochloride

**Orthogonal CRISPR Support**
- Enrichr LINCS_L1000_CRISPR_KO_Consensus_Sigs library
- For each chemical drug, check if its target appears as reversal-direction CRISPR KO hit
- Example: CDK inhibitors (AZD-5438, AT-7519) share CDK9 as one target
  - Query: kidney UP genes → reversal-direction KOs → CDK9 KO (Down) matches 7 genes, adj-p=0.1588
  - **Chemical-CRISPR convergence is suggestive target-plausibility evidence, not FDR-significant validation**
- Notable CRISPR-only signals include IL21R/CXCL13 in thymus and FOXO3 in liver; these remain exploratory because chemical matches and independent perturbation validation are absent.

**Multi-Tissue Pareto Optimization**
- Drug × tissue matrix: reversal score (NaN if not in top-50)
- Pareto dominance: drug A > drug B if mean_rev(A) ≥ mean_rev(B) AND min_rev(A) ≥ min_rev(B) AND at least one strict
- Current Pareto front: CGP-60474 (3 tissues, mean=0.1457, min=0.0735) and quinacrine hydrochloride (2 tissues, mean=0.1364, min=0.0909)
- Benefit: prioritizes breadth and minimum reversal score while keeping all drug language hypothesis-generating.

### Results

#### BRIDGE: Species Transfer at AUROC 0.888

Supervised classification of Inspiration4 post-flight vs pre-flight pathway conservation achieved AUROC 0.888 [0.854–0.918] when 6 mouse-tissue NES vectors augmented SpaceOmicsBench's 8 standard features. Improvement magnitude (+0.175 AUROC, 95% CI [0.134, 0.219]) was robust across 1000 bootstrap resamples, suggesting mouse tissue-level pathway NES is a strong cross-species predictor.

Tissue ranking for predictive power (ablation): eye > skin >> gastrocnemius > thymus ≈ kidney > liver. Mechanistically, this hierarchy may reflect immune system diversity (thymus shared T-cell developmental programs across species) vs tissue-specific metabolism (liver more divergent). Notably, gastrocnemius was directionally concordant (r=+0.381 I4 pathway NES) while thymus anti-correlated (r=–0.32), yet both remained predictive—indicating the classifier leverages inversion as information signal, possibly reflecting compensatory human gene expression.

#### DECOMPOSE: Interaction-Dominant Stressor Signatures

Factorial decomposition revealed interaction terms (HLU×IR, HLUxIRxT, Sex×HZE) dominate variance in top-200 responsive genes across all four analog studies (44–61% variance fraction). This non-additivity suggests spaceflight stressors are **mechanistically synergistic**, not independent pathways. For example, in thymus, HLU×IR interaction alone explained 61% of variance in proliferation-associated genes (e.g., *MKI67*, *HMGB1*), implying that combined stress triggers a qualitatively different immune response than either stressor alone.

Causal decomposition via ICP-style stability scoring across 9 environments ranked Time (ICP=0.540) and higher-order interactions (HLUxIRxT=0.522) as the most invariant candidate effects across experimental contexts. Conversely, flight-pooled data (ICP=0.283) was least stable, consistent with hidden confounders unique to orbital missions.

**Radiation Quality Matters:** Comparison of low-LET γ-rays (OSD-211, r_IR=+0.36 vs thymus flight) to high-LET HZE (OSD-719, r_HZE=–0.22) revealed **opposite-sign responses**. This fundamental discrepancy—despite nominally similar absorbed dose (Gy)—underscores that standard dose metrics mask biological equivalence. For Mars-like GCR modeling, the result flags dose-equivalent ISS scaling as a hypothesis requiring separate high-LET validation.

#### DECOMPOSE: Mars Extrapolation Flags ~1000× Gene Sensitivity

Applying factorial coefficients to a 900-day Mars-like stressor vector (µg_avg=0.62×, HZE=12.9×, Time=7.5×) flagged dramatic sensitivity among interaction-driven genes. Top amplification flags include *WNT10B* (~1052×, brain-to-eye analog), *KRTAP19-2* (~414×, skin analog), and an unannotated spleen/thymus proxy gene ENSMUSG00000092534 (~190×). However, **linear extrapolation beyond 5× amplification breaks down**—Spearman r(projected, flight-observed) ≈ 0 across tissues, indicating threshold dose-response or mechanistic saturation not captured by first-order interactions.

This negative result is informative: it flags need for **non-linear
dose-response models** and mechanistic constraints from the causal DAG before
treating Mars projections as operational predictions. Bootstrap CIs on Mars
projection flags are wide (1.5–3 SD), emphasizing uncertainty scaling with dose
amplification.

#### INTERVENE: Multi-Tissue Pareto Front Identifies 2 Signature-Reversal Hypotheses

L1000CDS2 queries with tissue-specific signatures returned top reversals:
- **Thymus:** Albendazole, AZD-7762 (centromere/kinetochore targeting)
- **Gastrocnemius:** Dorsomorphin, Rosiglitazone (muscle atrophy reversal)
- **Liver/Kidney:** PD-184352 (MEK1/2), GDC-0980 (PI3K/mTOR)

Multi-tissue scoring identified 26 drugs active in ≥2 tissues. Pareto filtering
(maximize mean_rev and min_rev among observed tissue hits) yielded 2 drugs on
the current efficiency frontier: **CGP-60474** (3 tissues, mean_rev=0.1457,
min_rev=0.0735) and **quinacrine hydrochloride** (2 tissues, mean_rev=0.1364,
min_rev=0.0909). Dorsomorphin remains a tissue-level metabolic-axis hit in
gastrocnemius and liver, but it is not on the current Pareto front.

**CRISPR Orthogonal Support:** CDK inhibitors (AZD-5438 and AT-7519) showed a
suggestive match to CDK9 KO in the Enrichr/LINCS CRISPR library for kidney UP
genes (7 overlapping genes; adj-p=0.1588). This chemical-genetic convergence is
useful target-plausibility evidence, but it is exploratory and does not
establish safety, dose, efficacy, or FDR-significant validation.

**Novel CRISPR-Only Targets:** Enrichr/LINCS CRISPR signatures identified
FOXO3, IL21R, and CXCL13 as reversal-associated KO signals without corresponding
high-ranking chemical matches. These are exploratory target axes for follow-up,
not drug-development claims.

#### Causal Integration: The SpaceMed DAG

Combining results from all three pillars, we construct a causal DAG:
```
Stressors (µg, GCR, isolation)
  ↓ [ICP-ranked stability edges]
Tissue Responses (pathway-level, tissue-specific)
  ↓ [species transfer, pathway conservation]
Human Health Outcomes (immune dysregulation, muscle atrophy, organ dysfunction)
  ↓ [signature reversal hypotheses]
Perturbation hypotheses (CDK-axis, AMPK/BMP-axis, and tissue-specific perturbations)
```

A future counterfactual module on this DAG would support explicit, testable
queries such as: *"If an intervention targets AMPK signaling before flight, how
would the model project pathway dysregulation under a Mars-like stressor regime?"*
Those queries should be promoted only after intervention assumptions and
falsification checks are documented.

### Discussion

**Strengths of SpaceMed**
1. **Integrated causal-evidence framework** for spaceflight biology across 3 stressors, 8 tissues, 2 species, 100+ datasets
2. **Frames three analysis questions** relevant to Artemis/Mars hypothesis
   generation
3. **Radiation quality discovery** (low-LET vs high-LET) supports the need for mechanistic models beyond dose metrics
4. **Species transfer at AUROC 0.888** demonstrates pathway conservation despite gene-level divergence
5. **Multi-tissue Pareto triage** reduces the chance of prioritizing perturbations that reverse one tissue signature while aggravating another

**Limitations & Future Work**
1. **SOMA n=4 crew** limits statistical power; treat as existence proof; larger cohorts (TRISH-funded missions) will validate
2. **Mars linearity breakdown** indicates need for mechanistic saturation models and dose-response thresholds
3. **ICP identifiability** depends on ≥2 environments; some stressor pairs may remain unidentifiable
4. **Chemical drug coverage** incomplete; LINCS only ~10,000 perturbations, missing natural products, biologics
5. **Human analog data sparse** (isolation); Mars500/HI-SEAS valuable but n<20 crews

**Implications for Artemis & Mars**
- Thymus-linked immune-reconstitution biology should be treated as a priority
  hypothesis for validation before any mission-planning use
- Muscle atrophy (gastrocnemius) signature-reversal hypotheses (for example,
  AMPK-axis perturbations) are candidates for pre-clinical rodent study design
- Radiation-specific damage on Mars may be mechanistically opposite to ISS,
  motivating separate countermeasure-hypothesis discovery rather than direct ISS
  scaling
- Projected interaction-driven pathway flags (WNT, sarcopenia networks) are
  non-linear, so mechanistic models are essential before any deployment claim

**Broader Context**
SpaceMed exemplifies causal systems biology: integrating multi-species, multi-stressor, multi-tissue data into a reproducible evidence map for specific astronaut-health hypotheses. This framework is portable to terrestrial medicine (aging, inflammation, cancer), where similar stressors (oxidative stress, hypoxia, isolation) induce overlapping transcriptomic signatures.

### Conclusion
We present SpaceMed, a causal analysis framework linking mouse OSDR signatures
to human spaceflight outcomes via pathway conservation, decomposing stressor
confounding to flag Mars-regime extrapolation limits, and prioritizing
multi-tissue signature-reversal hypotheses with orthogonal CRISPR support. This
framework quantitatively frames three analysis questions for deep-space biology:
*Which mouse pathway signals transfer to humans? Which Mars-like stressor regimes
break linear extrapolation? Which perturbations reverse tissue signatures enough
to prioritize pre-clinical study designs?* Results prioritize AMPK-axis and CDK-axis
hypotheses for rodent validation; Mars missions require stressor-specific
countermeasure-hypothesis discovery distinct from ISS protocol extrapolation.

---

## Extended Methods

### Data Availability
- OSDR: All bulk RNA-seq via GLDS (https://genelab.nasa.gov/)
- SOMA / Inspiration4: Overbey et al. 2024 supplementary data
- LINCS L1000CDS2: Maayanlab public API
- Enrichr LINCS_L1000_CRISPR_KO_Consensus_Sigs
- Code & reproducibility: github.com/[repo] / https://zenodo.org/[doi]

### Code Repositories
- Pillar 1 (BRIDGE): `v8/bridge/`
- Pillar 2 (DECOMPOSE): `v8/decompose/`
- Pillar 3 (INTERVENE): `v8/intervene/`
- Figures & manuscript: `v8/figures/` + `v8/MANUSCRIPT_DRAFT.md`

### Supplementary Tables
- **S1:** Species transfer NES heatmap (6 tissues × 20 I4 compartments × 16K pathways)
- **S2:** Factorial model coefficients, SE, p-val (per analog, per gene)
- **S3:** Mars extrapolation projection flags (per gene, per analog)
- **S4:** L1000CDS2 top-50 perturbation hypotheses per tissue
- **S5:** CRISPR orthogonal-support matrix (chemical targets vs KO library)
- **S6:** ICP scores per stressor per gene (top 100 causal genes)

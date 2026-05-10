# v8 — SpaceMed: Causal, Multi-Species Platform for Astronaut-Health Hypotheses

Status: **In progress (2026-04-17)** | Depends on: v1–v7

v8 extends the GeneLab benchmark beyond "can ML detect spaceflight?" to three
research questions for Artemis and Mars hypothesis generation:

1. **BRIDGE** — Which rodent pathway signals transfer to human spaceflight data? (Inspiration4 SOMA link)
2. **DECOMPOSE** — Where do Mars-like stressor regimes break linear extrapolation? (µg / GCR / isolation factorial)
3. **INTERVENE** — Which perturbations reverse tissue signatures enough to justify follow-up testing? (LINCS chemical + CRISPR-KO signature reversal)

Unified by an **invariant-causal-prediction evidence map** for future
counterfactual modeling. The current DAG ranks stability candidates; it does not
yet support operational do-calculus or crew countermeasure selection.

## Layout

```
v8/
├── bridge/       # Pillar 1 — rodent → human (SOMA)
├── decompose/    # Pillar 2 — Mars stressor factorial
├── intervene/    # Pillar 3 — signature-reversal hypotheses
├── causal/       # integrative causal evidence map
├── evaluation/   # result JSONs
├── figures/      # publication figures
└── provenance/   # run manifests linking claims to inputs, commands, outputs
```

See `ARTIFACTS.md` for the v8 tracking policy. Large GEO caches and raw
SpaceOmicsBench/SOMA inputs are intentionally kept out of Git.

Implementation planning and release gates are tracked in
`docs/V8_IMPLEMENTATION_RESEARCH.md`. HPC entrypoints live under `scripts/`:
`hpc_v8_bridge.sh`, `hpc_v8_decompose.sh`, `hpc_v8_intervene.sh`,
`hpc_v8_causal.sh`, and `hpc_v8_summary.sh`.

## External data source

Human spaceflight data is supplied by a parallel **SpaceOmicsBench** checkout,
configured with `SPACEOMICS_ROOT` or discovered as a sibling directory at
`../SpaceOmicsBench`. That project curates Inspiration4 (GLDS-530 cfRNA +
GLDS-563 proteomics/metabolomics), NASA Twins, JAXA CFE, and Axiom-2 data with
pre-computed fGSEA NES tables:

- `v2_public/data/processed/gt_conserved_pathways_GeneLab.csv` — mouse reference (292 pathways)
- `v2_public/data/processed/gt_conserved_pathways_i4_pbmc.csv` — I4 PBMC × 9 celltypes
- `v2_public/data/processed/gt_conserved_pathways_NASA_Twins.csv` — Twins × 12 compartments
- `v2_public/data/processed/cross_mission_pathway_features.csv` — labeled cross-mission features
- `v2_public/data/processed/conserved_pbmc_to_skin.csv` — gene-level PBMC↔skin mapping

## Reuses (do not rewrite)

- `scripts/run_fgsea.R`, `scripts/compute_pathway_scores.R` — pathway scoring
- `scripts/evaluate_submission.py` — AUROC + bootstrap + permutation
- `scripts/cross_tissue_transfer.py` — transfer scaffolding
- v5 DGIdb/ChEMBL druggability annotations
- v7 `unified/signal_hierarchy.py` — method comparison
- SpaceOmicsBench `missions/scripts/ingest_osdr.py` — OSDR→benchmark CSV
- SpaceOmicsBench `v2_public/evaluation/signature_query.py` — hypergeometric signature search

## Pillar 1 Results (2026-04-17)

### 1a. Aggregate mouse → human pathway-NES conservation (`bridge/link_spaceomicsbench.py`)

| Mission | n compartments | Mean Spearman r | Max r | Mean |NES|≥1 concordance |
|---|---|---|---|---|
| **Inspiration4 (I4) PBMC** | 9 celltypes | **0.777** | 0.842 | **95.5%** |
| NASA Twins | 21 | 0.351 | 0.646 | 72.3% |

Far above v6's raw-expression cfRNA baseline (AUROC ~0.5). Confirms H2
(batch effects are gene-driven, pathway-absorbed) extends cross-species.

### 1b. Per-tissue mouse NES refined (`bridge/tissue_nes_bridge.py`)

Breaking the aggregate into 6 mouse tissues reverses the expected ranking:

| Mouse tissue | Mean r vs I4 PBMC (8 celltypes) | Mean sign concordance |NES|≥1 |
|---|---|---|
| **gastrocnemius** | **+0.599** | **97.4%** |
| thymus | −0.211 | 5.9% |
| liver | −0.358 | 6.5% |
| skin | −0.360 | 17.5% |
| kidney | −0.531 | 4.8% |
| eye | −0.543 | 0.0% |

**Gastrocnemius is the uniquely positive cross-species bridge** to human peripheral
immune pathway response. All other mouse tissues run *opposite-sign* — consistent
with v1's Drosophila-mouse negative correlation (r=−0.19 to −0.59, README.md:150).
Biological read: spaceflight muscle catabolism/oxphos/inflammation pathways match
peripheral-blood immune activation direction; solid-tissue resident programs invert
in circulating leukocytes.

### 1c. Supervised I2 conservation with mouse prior (`bridge/supervised_conservation.py`)

Adding 6 mouse-tissue NES features to SpaceOmicsBench's 8 I4 baseline features on
I2 (n=452 pathways, conserved-in-Twins binary), 5-fold CV, bootstrap CI (N=1000):

**With Hallmark/KEGG/Reactome/MitoCarta only (7% I2 pathway coverage):**

| Model | Base AUROC | +Mouse AUROC | Δ (bootstrap 95% CI) |
|---|---|---|---|
| LR (balanced) | 0.549 | 0.571 | +0.020 [−0.006, +0.044] |
| **RF (500 trees)** | 0.710 | **0.723** | **+0.013 [+0.0005, +0.026]** ✓ |

### 1d. MSigDB expansion + full re-run (`bridge/tissue_nes_bridge.py`, `bridge/supervised_conservation.py`)

`scripts/run_fgsea.R` extended with `c2cgp` (3,537 chemical/genetic perturbations)
and `c5bp` (7,580 GO Biological Process) gene sets — total 8,428 unique pathways
across all 6 tissues. Re-run 2026-04-18.

**Tissue NES bridge with full MSigDB (8,428 pathways):**

| Mouse tissue | Mean r vs I4 PBMC (9 compartments) |
|---|---|
| **gastrocnemius** | **+0.381** |
| kidney | −0.088 |
| skin | −0.267 |
| liver | −0.289 |
| thymus | −0.371 |
| eye | −0.553 |

Gastrocnemius remains uniquely positive. Mean r dropped from +0.599 to +0.381 because C2/C5 dilute with tissue-biology pathways that don't transfer to blood — but directional ranking is preserved.

**Supervised with full MSigDB (+6 mouse-tissue NES, full I2 pathway overlap):**

| Model | Base AUROC | +Mouse AUROC | Δ (bootstrap 95% CI) |
|---|---|---|---|
| LR (balanced) | 0.549 | **0.732** | **+0.185 [+0.128, +0.246]** ✓✓ |
| **RF (500 trees)** | 0.712 | **0.888** | **+0.175 [+0.134, +0.219]** ✓✓ |

Leave-one-tissue-out ablation (LR): **eye** (Δ=−0.019) and **skin** (Δ=−0.019) are the largest contributors with full MSigDB coverage — suggesting these tissues carry the C2/C5 pathway information predictive of I4↔Twins conservation, even though their NES anti-correlates with human PBMC direction (directional inversion is informative to the classifier).

**The projected AUROC lift materialized: RF 0.712 → 0.888 (+0.175), entirely from adding mouse NES features.**

**Leakage audit:** `bridge_leakage_audit.json` records the supervised feature
contract on HPC: the `label` column is excluded from all 14 model features, the
I2 pathways and aggregated mouse NES pathways are unique under the `pathway`
merge key, deterministic 5-fold stratified splits are hashed, and no single
feature nearly perfectly reproduces the label (max single-feature oriented
AUROC=0.645). This audits the local v8 modeling contract; freezing the upstream
SpaceOmicsBench feature-builder command remains a beta requirement.

## Pillar 3 First Pass (2026-04-17)

### 3a. Tissue signature export (`intervene/export_signatures.py`)

Aggregates GeneLab DESeq2 Wald statistics across each tissue's missions via
per-mission z-score → mean. Emits top-150 up / top-150 dn in human-symbol
format (LINCS/CMap-ready). Striking biology surfaces without any ML:

| Tissue | Top up | Top down |
|---|---|---|
| Thymus | CLEC2D, DSE, MUC16, HPGD | H2AC14, RRM2, PCLAF, H2AZ1 (proliferation/histones — thymic involution) |
| Gastro | ENHO, PPP1R3C, SEC14L5 | **PER3, PER2**, TNNT1, TNNC1, MYLK4 (circadian + contractile) |
| Liver | CTRC, TRY4, **BMAL1**, CDKN1A | CLPS, MUC6, MUP-PS6 |
| Kidney | IFI202B (IFN), MYH4, CCL21D | GVIN1, IGHV1-78, PCBP2 |
| Skin | HOXA5, ART2B, FBP2 | DEFB8, TIMP1, WFDC12 |
| Eye | SST, PTPRN, PAPPA2 | NME2, NR4A1, UPK3B, **PER2** |

Circadian-clock genes (PER2, PER3, BMAL1) appear as flight-responsive across
multiple tissues — an orthogonal confirmation of a known space-biology
signal without explicit targeting.

### 3b. L1000CDS2 signature reversal (`intervene/lincs_query.py`)

First genome-wide LINCS query using GeneLab spaceflight signatures (2026-04-18).
L1000CDS2 (aggravate=False) queries top-50 up ∪ top-50 down per tissue against
the full L1000 chemical perturbation database; scores by gene-set overlap reversal.

| Tissue | Top reversal hit | Mechanism | Reversed genes (dn→up) |
|---|---|---|---|
| Thymus | **ALBENDAZOLE** | Microtubule/anti-parasitic | CENPA, CENPE, MKI67 (anti-involution) |
| Thymus | **AZD-7762** | CHK1/CHK2 checkpoint inhibitor | CENPE, CENPF, MKI67 |
| Gastrocnemius | **Dorsomorphin** | AMPK/BMP inhibitor | TNNT1↑, PCK1↑, TNNC1↑ (anti-atrophy) |
| Gastrocnemius | Geldanamycin | HSP90 inhibitor | DNAJB9, TNNC1, TNNT1 |
| Liver | **Dorsomorphin** | AMPK/BMP inhibitor | KRT23, NPC1 |
| Kidney | **PD-184352** | MEK1/2 inhibitor | NR4A1↓, PCBP2↑ |
| Kidney | **GDC-0980** | PI3K/mTOR inhibitor | TXNIP↑, TSPAN4↑ |
| Eye | Penfluridol | Antipsychotic / Ca²⁺ channel | cell cycle/proliferation axis |

**Key biological insights:**
1. **Dorsomorphin (AMPK/BMP inhibitor)** appears as a top reversal hit in both gastrocnemius AND liver, nominating AMPK/BMP signaling as a shared metabolic-tissue follow-up axis. Transcriptomic reversal alone does not establish anti-atrophy efficacy.
2. **Thymus involution hypothesis**: Albendazole + AZD-7762 (CHK1 inhibitor) both induce centromere/kinetochore proliferation genes that are suppressed in spaceflight-involuated thymus. This is a biologically coherent pro-lymphopoiesis hypothesis for controlled validation.
3. **MEK/mTOR signaling** (kidney): PD-184352 + GDC-0980 target innate immune/metabolic stress genes (NR4A1, TXNIP) that are perturbed in renal spaceflight response.

All top-50 per-tissue CSVs: `v8/intervene/evaluation/lincs_{tissue}_top50.csv`

### 3c. Multi-tissue Pareto front (`intervene/pareto_multi_tissue.py`)

215 unique drugs identified; 26 appear in ≥2 tissue top-50 lists. Pareto-dominated
drugs removed (a drug A dominates B if mean AND min reversal both ≥).

**Drugs hitting 3–4 tissues (broadest cross-tissue signature-reversal hypotheses):**

| Drug | N tissues | Tissues | Mechanism | Mean rev score |
|---|---|---|---|---|
| **AZD-5438** | **4** | skin, eye, liver, kidney | CDK1/2/9 inhibitor | 0.114 |
| **Geldanamycin** | **4** | gastrocnemius, skin, liver, kidney | HSP90 inhibitor | 0.096 |
| **CGP-60474** | 3 | gastrocnemius, skin, eye | CDK1/2/5 inhibitor | 0.146 |
| Mitoxantrone | 3 | skin, eye, kidney | TOP2/anthracenedione | 0.118 |
| AT-7519 | 3 | gastrocnemius, skin, eye | CDK1/2/4/6/9 | 0.103 |
| KU-0063794 | 3 | eye, liver, kidney | mTORC1/2 inhibitor | 0.070 |

**Pharmacological convergence:**
1. **CDK inhibitors dominate** (AZD-5438, CGP-60474, AT-7519): Spaceflight-driven cell cycle dysregulation is a recurrent reversible signature across solid tissues. CDK biology is therefore a target axis for follow-up, not a treatment recommendation.
2. **HSP90 inhibitors** (geldanamycin, NVP-AUY922): Heat shock response is spaceflight-upregulated; HSP90 perturbation reverses parts of the downstream client-protein signature and needs independent toxicity-aware validation.
3. **mTOR pathway** (GDC-0980, KU-0063794): mTOR integrates mechanical sensing and metabolic adaptation, making it a plausible microgravity-response axis for controlled perturbation studies.

**Pre-clinical triage:** AZD-5438 (4 tissues, CDK) and Geldanamycin (4 tissues, HSP90) are
high-priority signatures for NASA Rodent Research-style validation because they
span multiple tissues and have existing pharmacology. Safety, dose, and efficacy
remain untested for this use case.

**Current Pareto front:** after maximizing both mean reversal and the minimum
observed tissue reversal score, the validated front contains **CGP-60474** (3
tissues; mean=0.146, min=0.074) and **quinacrine hydrochloride** (2 tissues;
mean=0.136, min=0.091). The broad CDK/HSP90 hits remain useful mechanistic axes,
but they are not identical to the present Pareto frontier.

Outputs: `v8/intervene/evaluation/pareto_front.csv`, `multi_tissue_drug_matrix.csv`,
`multi_tissue_drug_scores.json`.

**API snapshot reproducibility:** `api_snapshot_manifest.json` records
deterministic hashes for the L1000CDS2 query payloads, Enrichr CRISPR query
gene lists, and all tracked parsed API outputs. The manifest does not re-call
external APIs and does not store raw response dumps; concrete upstream
db-version pinning remains a beta-freeze item.

**Safety-aware triage:** `safety_triage.csv` joins current multi-tissue reversal
scores with class-level safety liabilities. The table keeps CGP-60474 and
quinacrine as Pareto-front pathway hypotheses, while broad CDK/HSP90/TOP2 and
tool-compound hits remain mechanism flags only, not countermeasure
recommendations.

### 3d. Enrichr/LINCS CRISPR KO orthogonal support (`intervene/perturb_seq_orthog.py`)

Queries the Enrichr `LINCS_L1000_CRISPR_KO_Consensus_Sigs` library (3,742 gene
perturbations) with each tissue's UP/DN signatures. Checks whether Pillar 3b
drug hits' known targets surface as CRISPR-KO reversers, providing orthogonal
support where chemical and genetic perturbation signatures converge.

**Top CRISPR-KO reversers per tissue:**

| Tissue | UP-reversal KOs (KO→Down genes overlap spaceflight UP) | DN-reversal KOs (KO→Up genes overlap spaceflight DN) |
|---|---|---|
| Thymus | IL21R, HAGH, CXCL13 (lymphoid signaling) | ATP1B4, DOCK4, PIP4K2C |
| Gastrocnemius | ARSG, ETV1 (skeletal muscle TF), SUCLA2 | MAP4K1/ACTN4, ICAM2, ETS1 |
| Skin | CTSS, FOLR3, POU1F1 | SMTNL1, NEXN, NECTIN4 |
| Eye | NDUFA13 (complex I), HDAC1, PIGL | RNF212, PSMD1, ACHE |
| Liver | **FOXO3** (metabolism master TF), UCK1, TSPAN13 | PTPN2, NNT, SAMD4A |
| Kidney | ASCL2, SAMD4B, GNG4 | ODC1, NEU2, SRMS |

**Orthogonal drug-target convergence (kidney):**

| Chemical drug (Pillar 3b) | CRISPR KO target with matching support | Overlap genes | adj-p |
|---|---|---|---|
| **AZD-5438** (CDK1/2/9) | CDK9 KO↓ reverses kidney UP genes | CFD, EGR1, NR4A1, ASNS, PLIN1, MEOX2, CLDN1 | 0.16 |
| **AT-7519** (CDK1/2/4/6/9) | CDK9 KO↓ reverses kidney UP genes | CFD, EGR1, NR4A1, ASNS, PLIN1, MEOX2, CLDN1 | 0.16 |

**This is orthogonal support, not efficacy validation**: two independent chemical
CDK inhibitors and the CRISPR knockout of their shared target (CDK9) point in
the same reversal direction for 7 genes spanning innate immunity (CFD),
immediate early stress response (EGR1, NR4A1), amino acid metabolism (ASNS),
adipocyte biology (PLIN1), vascular TF (MEOX2), and epithelial junctions
(CLDN1). The adjusted p-value is exploratory, so this remains a kidney
spaceflight-signature hypothesis.

**Biologically notable CRISPR hits with no chemical precedent yet:**
- **FOXO3 KO** reverses liver spaceflight UP — FoxO3 is the master regulator of hepatic autophagy/gluconeogenesis; its induction in spaceflight liver has long been hypothesized but never orthogonally confirmed.
- **IL21R/CXCL13 KO** reverses thymus UP — the lymphoid survival/chemotaxis axis directly targets thymic-involution biology.
- **ETV1 KO** reverses gastrocnemius UP — ETV1 is a skeletal-muscle-specific transcription factor; its KO would downregulate muscle genes, which matches the atrophy direction in spaceflight.

All CRISPR hit tables: `v8/intervene/evaluation/crispr_orthog_{tissue}_{up,dn}.csv`

### 3e. Offline drug-reversal scorer (`intervene/offline_reversal_scorer.py`)

Hypergeometric enrichment of DGIdb drug-target sets (from v5 consensus panel)
against each tissue's top-150 up∪down gene set. LINCS-independent first pass
(bg=1 per drug via panel source; proper scaling needs full DGIdb query).

Biologically coherent top hits surface with FDR-pending (p<0.05 uncorrected):

| Tissue | Top hit | Target | Note |
|---|---|---|---|
| Gastrocnemius | **BEPRIDIL, LEVOSIMENDAN** | TNNC1↓ | Cardiac/muscle troponin C sensitizers — muscle atrophy axis |
| Liver | **MAVACAMTEN** | MYH7↑ | Cardiac myosin inhibitor; ectopic MYH7 up in liver is novel |
| Kidney | **DANICOPAN** | CFD↑ | Complement factor D inhibitor (PNH); innate immunity axis |
| Thymus | TOP2 poisons (doxorubicin etc.) | TOP2A↓ | Mimic, not reverse — directional caveat |
| Skin | LIOTHYRONINE, DEXAMETHASONE | THRSP↓ | Thyroid/lipid dermal axis |

Pipeline smoke-tested; next step is LINCS/L2S2 query with full drug×gene matrix
for properly-powered reversal scores.

## Pillar 2 First Result (2026-04-17)

### 2a. Factorial HLU × radiation decomposition (`decompose/factorial_analog.py`)

OSDR 2×2 factorial analog studies discovered cached in `data/mouse/`:

| Study | Tissue | n | Design |
|---|---|---|---|
| OSD-211 | spleen | 21 | HLU × 0.04 Gy 57Co γ |
| OSD-237 | skin | 21 | HLU × 0.04 Gy 57Co γ |
| OSD-202 | retina+brain | 84 | HLU × 0.04 Gy 57Co × time (not yet fit) |

Per-gene OLS: `log2(norm+1) ~ HLU + IR + HLU:IR`. Compared β coefficients
to flight Wald-z signatures via Spearman on flight-responsive genes (|z|>1):

| Analog → Flight | r_HLU | r_IR | r_interaction | Top-200 var attribution |
|---|---|---|---|---|
| **Spleen → Thymus** | −0.193 | **+0.361** | −0.323 | HLU 27% / IR 11% / **Int 62%** |
| Skin → Skin | −0.088 | −0.188 | +0.050 | HLU 46% / IR 14% / Int 40% |

**Key finding**: thymus flight response aligns best with radiation (r=+0.36),
with strong antagonistic µg×rad interaction dominating variance attribution.
On Mars (higher GCR than ISS, same-order µg), thymic consequences could be
**disproportionately worse** than linear extrapolation would predict.

Skin factorial explains little of flight skin signal — the 0.04 Gy low-LET
57Co analog may be too weak. HZE (Fe-56, Si-28 at NSRL; OSD-73 human, OSD-109
cardiomyocyte) extension needed for Mars-realistic GCR.

Significant-gene counts per factor (p<0.05 uncorrected):
- Spleen: HLU=1742, IR=821, HLU×IR=1993
- Skin:   HLU=5621, IR=659, HLU×IR=1511  (HLU-dominant in numerosity)
- Brain:  HLU=3436, IR=3414, **T=6028 (time largest)**, HLUxIR=3559, HLUxT=1850, IRxT=1459, HLUxIRxT=1503

### 2b. Full analog × flight-tissue matrix (2026-04-17)

Expanded factorial to **OSD-202 brain (2×2×2 HLU × IR × Time)**. Comparing each
analog's coefficients against every flight tissue signature reveals cross-tissue
spillover patterns not visible in same-tissue comparison:

Spearman r (flight-responsive |z|>1) between analog coefficients and flight Wald-z:

| Analog | flight tissue | r_HLU | r_IR | r_int | r_T | Note |
|---|---|---|---|---|---|---|
| Spleen | **thymus** | −0.20 | **+0.36** | −0.32 | — | Radiation-axis match |
| Spleen | **eye** | **−0.43** | −0.23 | +0.33 | — | Eye has strong anti-µg pattern |
| Spleen | liver | −0.30 | −0.37 | +0.35 | — | Interaction-dominant |
| Skin | **thymus** | **−0.51** | +0.31 | **+0.52** | — | Large sign-flip + interaction |
| Brain | **liver** | **+0.28** | +0.18 | — | +0.25 | Brain→liver axis (HLU × Time) |
| Brain | eye | +0.16 | +0.04 | — | +0.05 | Same-region weak |

**Novel biology surfaces**:
1. **Cross-tissue spillover**: ground-stressor responses in spleen/skin/brain
   analogs sometimes predict flight response in *different* tissues better than
   in own tissue — consistent with systemic inflammatory/metabolic axes dominating
   on flight.
2. **Time factor is the biggest single contributor to brain analog** (6028 sig
   genes vs 3436 HLU + 3414 IR) — ground-based aging/isolation at 4 vs 1 month
   moves more genes than either stressor alone.
3. **Skin analog µg response is anti-correlated with thymus flight** — interesting
   candidate for biomarker cross-talk (dermal HLU bystander signal).

All coefficient CSVs written to `v8/decompose/evaluation/factorial_{analog}_betas.csv`.

### 2c. HZE analog factorial (OSD-719 Sex × mixed_radiation_field, 2026-04-19)

First high-LET analog added: OSD-719 (Burke et al. 2023, PMID 38409284),
mouse adrenal/endocrine tissue, n=12, 2×2 Sex × HZE (NSRL mixed radiation field
= Fe-56 + Si-28 + p + C + O in flight-relevant fluence).

| Factor | n sig p<0.05 |
|---|---|
| Sex | 464 |
| HZE | 373 |
| Sex × HZE | 510 |

**HZE adrenal analog vs flight tissues** (Spearman r on flight-responsive |z|>1 genes):

| Flight tissue | r_Sex | r_HZE | r_SexxHZE | n |
|---|---|---|---|---|
| **Thymus** | −0.240 | **−0.219** | +0.160 | 780 |
| **Eye** | −0.058 | **−0.210** | +0.123 | 1,925 |
| Liver | +0.124 | −0.086 | +0.070 | 445 |
| Gastrocnemius | −0.051 | +0.053 | +0.002 | 2,282 |
| Skin | +0.114 | +0.101 | −0.090 | 800 |
| Kidney | +0.075 | +0.016 | −0.081 | 856 |

**Variance attribution (thymus, top-200 flight genes):**

| Factor | Var fraction |
|---|---|
| Sex | 32.4% |
| HZE | 23.6% |
| **Sex × HZE interaction** | **44.0%** |

**Novel findings:**
1. **Lymphoid + eye are the HZE-responsive flight tissues** — direction is anti-correlated (mirror of the gastrocnemius paradox; HZE adrenal response is opposite-sign from flight lymphoid/ocular).
2. **Sex × HZE interaction is the dominant factor** (44% of top-200 variance in thymus) — larger than either main effect. Implication for Mars studies: mixed-sex cohorts may show divergent HZE responses not captured by averaging.
3. Compared to 0.04 Gy ⁵⁷Co low-LET analog (OSD-211, thymus r_IR=+0.36), **high-LET HZE produces OPPOSITE-sign response** (r=−0.22). This confirms radiation quality matters fundamentally — dose equivalence does not equal biological equivalence at the transcriptome level.

### 2d. Mars 900-day extrapolation (`decompose/mars_extrapolate.py`)

Runs an exploratory linear projection of per-gene analog coefficients to
Mars-like mission conditions:

| Parameter | Analog reference | Mars setpoint | Scale factor |
|---|---|---|---|
| µg_dose | HLU full-unload = 1 | 0.62 (transit 1.0 + 0.38 surface) | 0.62× |
| HZE_dose | 0.04 Gy ⁵⁷Co analog | 350 mGy over 900 d | 12.9× |
| Time | 4 mo vs 1 mo | 900 d | 7.5× |

**Key result**: Linear extrapolation fails dramatically beyond the calibrated
analog dose range. Spearman r(Mars-projected, flight) is ≈0 across all three
tissues (thymus: r=0.0593; skin: r=−0.1263; eye: r=−0.0164), confirming Mars is
**not "more ISS"** — the transcriptomic response is qualitatively different
when stressors are scaled 12.9×.

**Top 5 Mars-amplified genes per tissue** (|mars_Δ / iss_Δ|, requiring |iss|>0.1):

| Tissue | Mars-amplified genes (fold vs ISS) |
|---|---|
| Spleen → thymus | ENSMUSG00000092534 (+182×), GM6192 (+134×), YBX2 (+129×) |
| Skin | KRTAP19-2 (+414×), CRISP1 (+248×), KRTAP5-5 (+196×) |
| Brain → eye | **WNT10B (+1052×)**, and many retinal-specific genes ≥1000× |

**Research implication for Mars modeling**: The factorial model flags possible
qualitative regime changes rather than operational point predictions.
Interaction terms (µg×IR×T) compound multiplicatively at Mars-like doses and
dominate over main effects in this linear stress test. WNT10B in brain/eye
response is particularly striking because Wnt signaling governs neuronal stem
cell renewal, making it a high-priority target for non-linear and tissue-specific
follow-up rather than a direct SANS prediction.

**Caveat**: linear extrapolation beyond calibration range is unreliable. Next
step requires either (a) non-linear dose-response models calibrated with HZE
data (OSD-73 bronchial, OSD-109 cardiomyocyte) or (b) mechanistic pathway
constraints from the ICP causal DAG.

**Bounded-dose sensitivity**: `mars_saturation_summary.json` applies cap5,
sqrt-after-5x, and log-after-5x guardrail scenarios to ask whether >5x Mars
amplification flags survive conservative damping. It finds spleen -> thymus
has 12,961 linear >5x flags, with 7,749 robust and 5,212 saturation-sensitive;
skin has 11,428 linear flags, with 8,254 robust and 3,174
saturation-sensitive; brain -> eye has 22,634 linear flags, but only 3,813
robust and 18,821 saturation-sensitive. Thus the brain/eye ~1000x linear
flags are the most sensitive to extrapolation assumptions and should be
treated as high-priority hypothesis flags, not point predictions.

Outputs: `v8/decompose/evaluation/mars_extrapolation_{analog}.csv`,
`mars_summary.json`, `mars_saturation_sensitivity_{analog}.csv`, and
`mars_saturation_summary.json`.

**Raw-cache readiness:** `raw_cache_audit.json` now records a complete HPC
bundle for the DECOMPOSE rerun: OSD-719, OSD-211, OSD-237, and OSD-202
count/sample-table caches are present and checksummed. The full factorial,
Mars, and bootstrap outputs were regenerated from that cache on HPC.

## Causal DAG — ICP across Stressor Environments (2026-04-18)

### 4a. Invariant Causal Prediction (`causal/build_dag.py`)

477,001-row environment matrix spanning 9 environments (3 analog tissues × 4 stressor
axes + 6 flight signatures), 35,559 genes. Per-(stressor, gene) ICP stability score =
1/(1+CV) across environments where that gene is observed in ≥2 settings.

**Stressor aggregate ICP (causally stable = consistent across environments):**

| Stressor | Mean ICP | Max ICP | N genes |
|---|---|---|---|
| **T (Time/Isolation)** | **0.540** | 0.961 | 82 |
| IRxT | 0.535 | 0.988 | 82 |
| HLUxIRxT | 0.522 | 0.962 | 82 |
| HLUxT | 0.488 | 0.934 | 82 |
| HLU (µg) | 0.407 | 1.000 | 27,029 |
| IR (radiation) | 0.398 | 0.999 | 27,029 |
| HLUxIR (interaction) | 0.390 | 1.000 | 27,029 |
| **combined_flight (ISS)** | **0.283** | 0.999 | 21,607 |

**Key interpretation:** Time/Isolation is the most causally stable stressor across
all analog environments (highest mean ICP). Combined ISS flight response has the
*lowest* ICP — confirming it conflates multiple stressors in non-invariant ways.
This is the computational analog of the ICP identifiability argument: ISS observations
alone cannot separate µg from GCR from isolation, while ground analogs each provide a
pure-stressor interventional environment.

**Top causal stressor → tissue edges (mean ICP over top-50 flight-responsive genes):**

| Stressor | Tissue | ICP score | N genes |
|---|---|---|---|
| IR | skin | 0.440 | 90 |
| HLU | thymus | 0.439 | 99 |
| HLU | skin | 0.424 | 90 |
| IR | kidney | 0.410 | 89 |
| HLU | gastrocnemius | 0.405 | 99 |
| IR | eye | 0.397 | 93 |

Skin radiation response and thymic µg response are the most causally invariant
tissue-specific edges — matching the factorial decomposition finding that
thymus is radiation-dominant and skin is HLU-dominant in ground analogs.

Outputs: `v8/causal/evaluation/icp_stressor_pathway_scores.csv`, `dag_edges.csv`,
`icp_dag_summary.json`.

## Status as of 2026-04-18

All three pillars have first-pass results. Summary of what's complete vs pending:

| Task | Status | Key result |
|---|---|---|
| Pillar 1a: Aggregate mouse→human NES (all tissues) | ✅ Done | I4 mean r=0.777 |
| Pillar 1b: Per-tissue NES bridge (4 DBs) | ✅ Done | Gastrocnemius r=+0.599 |
| Pillar 1c: Supervised conservation (4 DBs) | ✅ Done | RF AUROC +0.013 |
| Pillar 1d: MSigDB C2/C5 expansion fGSEA | ✅ Done | 8,428 pathways, all 6 tissues |
| Pillar 1d: Re-run bridge + supervised (full DBs) | ✅ Done | RF AUROC **0.712→0.888** (+0.175) |
| Pillar 2a: HLU×IR factorial (spleen, skin analog) | ✅ Done | Thymus radiation-dominant, Int=61% |
| Pillar 2b: HLU×IR×T brain 3-way + cross-tissue matrix | ✅ Done | T largest factor in brain |
| Pillar 3a: Tissue signature export (6 tissues) | ✅ Done | Circadian clock genes confirmed |
| Pillar 3b: L1000CDS2 reversal query (6 tissues) | ✅ Done | Dorsomorphin (gastro+liver), CGP-60474 (3 tissues), AZD-5438 (4 tissues) |
| Pillar 3c: Multi-tissue Pareto front | ✅ Done | CGP-60474 + quinacrine hydrochloride on current Pareto front |
| Pillar 3d: CRISPR KO orthogonal support | ✅ Done | CDK9 KO suggestively supports AZD-5438/AT-7519 target plausibility (kidney, 7 shared genes; adj-p=0.16) |
| Pillar 2c: HZE analog factorial (OSD-719 adrenal/endocrine) | ✅ Done | Sex×HZE dominates thymus variance (44%); HZE≠γ-rays sign |
| Pillar 2d: Mars 900-day extrapolation | ✅ First pass + bounded sensitivity | Linear model breaks at 12.9× dose; brain-to-eye flags are highly saturation-sensitive |
| Causal DAG: ICP stressor stability | ✅ Done | T/Isolation most causally invariant |
| Manuscript integration | ⏳ Pending | — |

## Milestones

| Sprint | Pillar | Deliverable | Status |
|---|---|---|---|
| M1–M3 | Bridge | Species-transfer AUROC matrix + RF 0.888 with mouse NES features | ✅ Complete |
| M4–M6 | Decompose | Factorial stressor model + ICP causal DAG | ✅ First pass complete |
| M7–M9 | Intervene + Integrate | LINCS reversal + Pareto front + CRISPR-KO orthogonal support + manuscript | 🔄 In progress |

## Next priorities

1. **Fetch HZE analog data** (OSD-73 human bronchial, OSD-109 cardiomyocyte) from OSDR for Mars-realistic GCR factorial extension and non-linear dose-response calibration.
2. **Manuscript figure set**: Species-transfer matrix, supervised AUROC waterfall (0.712→0.888), stressor ICP heatmap, multi-tissue drug Pareto front, chemical-CRISPR convergence figure.
3. **Mars calibration upgrade**: combine bootstrap CIs with bounded-dose sensitivity to select robust genes for HZE-calibrated non-linear follow-up.

## Target venues

*Nature* / *Nature Medicine* / *Cell* (integrated) or *Nat Methods* + *Nat Comms* + *Cell Rep Med* (split).

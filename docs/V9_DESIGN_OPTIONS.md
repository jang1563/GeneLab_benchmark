# v9 Design Options

Status: planning document
Depends on: v8 beta freeze and `docs/V9_EXTERNAL_DEEP_RESEARCH_2026_05_21.md`
Companion matrix: `docs/V9_SOURCE_AND_COMPETITOR_MATRIX.md`
Long-horizon execution plan: `docs/V9_LONG_HORIZON_EXECUTION_PLAN.md`
Organoid/species review: `docs/V9_ORGANOID_AND_SPECIES_EXTENSION_REVIEW_2026_05_21.md`

## North star

v9 should turn GeneLab Benchmark from a strong project-specific benchmark into a
community-usable benchmark platform for space biology AI.

Working name:

> SpaceBio-Bench

Working claim:

> SpaceBio-Bench evaluates whether biological AI models generalize across
> spaceflight missions, tissues, species, modalities, perturbations, and
> stressor regimes.

## Design constraints inherited from v8

- v8 beta must be frozen before v9 results are promoted.
- Bulk LOMO tasks keep macro-F1, balanced accuracy, AUROC, and calibration as
  primary metrics.
- cell-eval-style metrics are used only for scRNA-seq/AnnData tasks where DE
  recovery and perturbation discrimination are biologically meaningful.
- Mars and countermeasure outputs are hypothesis-generation artifacts unless
  independent validation exists.
- Every promoted v9 claim must link to source, code, output, metric, and
  validation manifests.

## Option matrix

| Option | Description | Novelty | Tractability | Risk | Best output |
|---|---|---:|---:|---:|---|
| A. Platform spine | Package, task registry, metric profiles, manifests, dataset bundle | High | High | Medium maintenance | Resource paper + release |
| B. Mission-held-out virtual-cell track | RRRM scRNA-seq/AnnData tasks aligned with cell-eval ideas | Very high | Medium | Model adapters and GPU friction | Methods/resource paper |
| C. Radiation-quality stressor track | Low-LET vs high-LET, HLU, time/isolation, nonlinear Mars-regime tests | Very high | Medium | Data sparsity and causal overclaim risk | Science paper |
| D. Intervention hypothesis track | LINCS/Enrichr/CRISPR signature triage with safety gates | Medium | High | Countermeasure overclaim risk | Supplementary module |
| E. Gated human data protocol | EXPAND/SOMA restricted-track schema and governance | High | Medium-low | Access and privacy | Protocol / future-ready track |
| F. Human organoid and multispecies extension | Public human iPSC organoid RNA-seq/proteomics plus non-mouse OSDR species | High | Medium | Small-N, modality mixing, orthology complexity | Public extension track |

Recommended v9 shape:

1. Build Option A as the required spine.
2. Choose Option B as the first flagship public benchmark.
3. Keep Option C as the second flagship or a v9.1 paper.
4. Keep Option D as hypothesis-only.
5. Prepare Option E as a schema/protocol, not a public dependency.
6. Add Option F only after schema support prevents mouse tissue, human
   organoid, and non-mouse species conflation.

## Recommended v9 architecture

```text
spacebio_bench/
  tasks/
    bulk_lomo/
    bridge_cross_species/
    sc_spaceflight/
    stressor_radiation_quality/
    intervention_hypothesis/
  metrics/
    genelab_minimal.py
    genelab_full.py
    genelab_sc.py
    stressor_regime.py
  data/
    loaders.py
    manifests.py
  submissions/
    validate.py
  reports/
    cards.py
```

This package does not need to replace the current repo layout immediately. It
can start as a thin API over existing processed tasks and v8 outputs.

## Task taxonomy

| Task family | Input | Output | Split | Primary metrics | Public-ready |
|---|---|---|---|---|---|
| Bulk LOMO | Bulk RNA-seq features | Flight vs Ground prediction | Leave-one-mission-out | macro-F1, balanced accuracy, AUROC | Yes |
| Tissue transfer | Tissue/task matrix | Held-out tissue prediction | Leave-one-tissue or tissue-pair | transfer F1, entropy, calibration | Yes |
| Cross-species BRIDGE | Mouse pathway NES + human summaries | Human alignment / post-flight state | Species/domain split | AUROC, rank correlation, calibration | Yes with summaries |
| Human organoid spaceflight | Human iPSC-derived neural organoid RNA-seq | LEO vs ground or response-signature recovery | Donor/type/microglia-blocked pilot | balanced accuracy, AUROC, DE direction match, rank correlation | Yes, pilot only |
| Multispecies spaceflight | Fly/plant/other OSDR transcriptomics | Species-specific or ortholog/pathway transfer | Species or study held out | species-specific F1, rank correlation, pathway concordance | Yes after schema expansion |
| scRNA-seq spaceflight | AnnData per cell type | DE recovery or state prediction | Held-out mission/cell type | overlap_at_N, direction match, mission discrimination | Yes |
| Radiation quality | Analog stressor tables | Regime classification / response projection | Held-out analog or stressor | sign stability, saturation sensitivity, uncertainty | Yes |
| Intervention hypothesis | Tissue signatures | Reversal/target triage | Tissue-held-out | reversal score, orthogonal support, safety flags | Yes, hypothesis-only |
| Human gated track | Restricted human multi-modal data | Hidden evaluation outputs | Access-controlled | pre-registered metrics only | Protocol first |

## Metric profiles

### `genelab_minimal`

- macro-F1
- balanced accuracy
- AUROC
- calibration error
- mission discrimination

### `genelab_full`

- all `genelab_minimal` metrics
- per-tissue F1
- cross-mission transfer entropy
- tissue-held-out degradation
- species/domain shift delta
- confidence interval by bootstrap

### `genelab_sc`

- DE overlap_at_N
- DE direction match
- perturbation or mission discrimination
- AUROC/AUPRC for state labels when applicable
- expression MAE only when a model predicts expression values directly

### `stressor_regime`

- low-LET/high-LET sign consistency
- stressor attribution stability
- nonlinear saturation sensitivity
- held-out analog correlation
- uncertainty width and calibration

## Flagship v9 paper options

### Paper 1: SpaceBio-Bench platform

Question:

Can a public benchmark test biological AI models under real spaceflight domain
shift rather than ordinary held-out samples?

Core contributions:

- task registry
- frozen public release
- metric profiles
- baseline model suite
- manifest-validated provenance
- living benchmark protocol

Why it is strong:

- It protects the v1-v8 work as infrastructure.
- It is less likely to become stale than a single-model leaderboard.
- It invites community submissions.

### Paper 2: Mission-held-out virtual-cell evaluation

Question:

Do single-cell and virtual-cell models recover spaceflight state changes under
held-out mission and held-out cell-type conditions?

Core contributions:

- RRRM-1/RRRM-2 AnnData task definitions
- cell-eval-aligned DE metrics where valid
- mission discrimination metric
- baseline adapters for simple models and selected foundation models
- explicit comparison to bulk LOMO behavior

Why it is strong:

- It connects GeneLab to the fastest-moving AI biology evaluation field.
- It preserves GeneLab's unique mission-held-out axis.

### Paper 3: Radiation-quality nonlinear stressor benchmark

Question:

Do transcriptomic models distinguish biologically different radiation and
combined-stressor regimes that have similar simplified dose summaries?

Core contributions:

- low-LET versus high-LET task construction
- HLU/radiation/time interaction tasks
- nonlinear saturation and uncertainty tests
- Mars-regime falsification framing

Why it is strong:

- It is scientifically distinctive.
- It builds directly on the v8 radiation-quality signal.
- It avoids overclaiming Mars point prediction.

## First 30 days

1. Freeze v8 beta or explicitly label v9 work as pre-v8-beta design only.
2. Create `spacebio_bench` skeleton with no new model claims.
3. Write a task manifest schema covering source, split, metric, and output.
4. Implement `genelab_minimal` and `mission_discrimination`.
5. Convert one RRRM scRNA-seq task into an AnnData-style v9 task card.
6. Create a public dataset inventory table with source URLs, access status,
   checksum status, privacy class, and release target.
7. Run only baseline methods first: logistic regression, random forest,
   PCA-LR, simple nearest-centroid, and DE baselines.

## First 60 days

1. Add `genelab_sc` metrics.
2. Add RRRM-1/RRRM-2 single-cell task cards.
3. Add method adapters for selected single-cell foundation or virtual-cell
   models only when installation and licensing are reproducible.
4. Add radiation-quality task definitions from existing v8 DECOMPOSE sources.
5. Generate first v9 baseline report.
6. Build Hugging Face dataset card draft for v9.

## First 90 days

1. Release v9-alpha benchmark snapshot.
2. Add RO-Crate-compatible metadata export.
3. Add Zenodo-ready release manifest.
4. Draft the SpaceBio-Bench platform manuscript.
5. Decide whether the first science flagship is virtual-cell evaluation or
   radiation-quality modeling.

## Acceptance criteria for v9-alpha

- At least three task families are runnable from one command.
- Every task has a manifest with source, split, metric, and output schema.
- Every metric has a documented biological interpretation.
- Baseline methods run without proprietary data.
- Public outputs can be packaged for GitHub plus Hugging Face plus Zenodo.
- No intervention output is phrased as a clinical or crew-health
  recommendation.

## Risks and guardrails

| Risk | Guardrail |
|---|---|
| v9 becomes too broad | Platform spine plus one flagship track only |
| Virtual-cell adapters consume all effort | Require simple baselines first |
| Human data access blocks release | Public core cannot depend on gated data |
| Countermeasure language overreaches | Keep INTERVENE as hypothesis-only |
| Mars projection overclaims | Use regime/falsification language |
| Data drift breaks reproducibility | Frozen snapshots plus separate living track |
| Benchmark becomes stale | Task API and metric profiles outlive model names |

## Immediate next implementation choices

Preferred first code path:

1. Add a v9 planning directory or package skeleton.
2. Define task manifest and metric profile schema.
3. Implement `mission_discrimination` on existing v1-v8 embeddings or summary
   tables where available.
4. Convert one RRRM single-cell task into a v9 task card.

Preferred first manuscript path:

1. Draft the SpaceBio-Bench platform outline.
2. Use v8 as the translational incubator evidence.
3. Use v9-alpha as the benchmark/platform contribution.

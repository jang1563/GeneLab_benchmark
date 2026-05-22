# v9 External Deep Research Memo

Date: 2026-05-21
Scope: External ecosystem scan before committing GeneLab Benchmark v9 design.
Status: planning reference; no v9 results are claimed here.
Companion docs: `docs/V9_DESIGN_OPTIONS.md`,
`docs/V9_SOURCE_AND_COMPETITOR_MATRIX.md`

## Executive conclusion

The strongest v9 opportunity is not a larger v8 intervention paper. It is a
public, living, provenance-first benchmark for space biology AI:

> SpaceBio-Bench: a mission-held-out, cross-species, perturbation-aware
> benchmark for spaceflight domain shift.

The external ecosystem has two fast-moving halves that are not yet well joined:

- Space biology data infrastructure is expanding through NASA OSDR, ALSDA,
  SOMA, and TRISH EXPAND.
- AI biology evaluation is converging around continuous benchmark frameworks
  such as OpenProblems and perturbation-aware virtual-cell evaluation suites
  such as Arc cell-eval.

GeneLab Benchmark already owns a complementary niche: leave-one-mission-out
generalization, cross-tissue/cross-species transfer, and spaceflight-specific
domain shift. v9 should make that niche explicit and reusable.

## Non-goals for v9

- Do not turn v9 into a clinical countermeasure recommendation project.
- Do not replace bulk LOMO metrics with cell-eval perturbation metrics.
- Do not build a generic "latest model leaderboard" that becomes stale as soon
  as the next virtual-cell model appears.
- Do not require controlled-access human datasets for the public benchmark
  core.
- Do not merge v8 beta freeze work with v9 design work.

## Source landscape

### 1. NASA OSDR is the public data spine

NASA OSDR combines GeneLab omics and ALSDA non-omics life-science data. The
public OSDR overview describes a repository spanning omics, physiology,
phenotype, behavior, imaging, and environmental data. The OSDR Biological Data
API exposes query and REST interfaces over metadata and data tables.

Implication for v9:

- Use OSDR accessions and API-derived metadata as first-class source records.
- Treat API query payloads, timestamps, filters, returned file URLs, and output
  checksums as part of the v9 task manifest.
- Keep public tasks reproducible without depending on local OSDR cache paths.

Primary sources:

- NASA OSDR: https://science.nasa.gov/biological-physical/data/osdr/
- OSDR Biological Data API: https://visualization.osdr.nasa.gov/biodata/api/

### 2. SOMA is the human-spaceflight anchor, not a large-N benchmark

The 2024 Space Omics and Medical Atlas paper presents an astronaut multi-omics
resource spanning public human spaceflight studies, including Inspiration4 and
other crew datasets. SOMA is strategically important because it connects rodent
and analog GeneLab signals to human spaceflight molecular measurements.

Limitations for v9:

- Human crew counts remain small.
- Data access, privacy, consent, and reidentification risk are central.
- SOMA should anchor human validation and bridge tasks, not become the only
  ground truth for broad model ranking.

Implication for v9:

- Use SOMA as a human validation layer for BRIDGE-like tasks.
- Separate public SOMA-derived lightweight summaries from controlled or large
  raw inputs.
- Preserve claims as cross-species/human-alignment evidence unless prospective
  validation data exist.

Primary sources:

- SOMA Nature 2024: https://www.nature.com/articles/s41586-024-07639-y
- SOMA Browser: https://soma.weill.cornell.edu/apps/SOMA_Browser/
- SOMA data and code collection: https://www.nature.com/collections/ebdbcahdgc/data-code

### 3. TRISH EXPAND is a future gated-data bridge

TRISH EXPAND is opening commercial spaceflight health and performance data to
researchers through a proposal and access process. The public program materials
describe multi-modal data such as omics, wearables, cognitive testing,
questionnaires, and other health/performance measures.

Implication for v9:

- v9 should be "EXPAND-ready" without requiring EXPAND access for the open
  benchmark.
- Define schema contracts for restricted human tracks now: subject-level
  metadata, timepoint harmonization, assay modality, privacy class, and allowed
  outputs.
- Keep a gated evaluation track separate from the public leaderboard.

Primary sources:

- TRISH EXPAND: https://cdn.bcm.edu/academic-centers/space-medicine/expand
- BCM release on EXPAND opening to researchers:
  https://www.bcm.edu/news/trish-opens-commercial-spaceflight-database-to-scientific-community

### 4. OpenProblems is the closest platform pattern

OpenProblems provides a community benchmark framework for single-cell analysis:
tasks, datasets, metrics, methods, reproducible runs, and leaderboard-style
evaluation. Its value for GeneLab is architectural rather than biological.

Implication for v9:

- Model v9 around a task registry plus metric profiles.
- Separate data loaders, method adapters, and metrics.
- Support frozen public snapshots and a living track that can refresh when new
  OSDR data appear.

Primary sources:

- OpenProblems Nature Biotechnology 2025:
  https://www.nature.com/articles/s41587-025-02694-w
- OpenProblems GitHub: https://github.com/openproblems-bio
- OpenProblems docs: https://openproblems.bio/

### 5. cell-eval is relevant only where the biology matches

Arc cell-eval and the Virtual Cell Challenge formalize perturbation-prediction
evaluation using metrics such as differential expression score, perturbation
discrimination score, and expression error. These metrics are valuable for
single-cell perturbation or DE-recovery tasks, but they are not replacements for
bulk mission-level classification metrics.

Implication for v9:

- Add a `genelab_sc` metric profile for RRRM-1/RRRM-2 and future AnnData tasks.
- Use DE overlap, direction match, and perturbation/mission discrimination for
  scRNA-seq tracks.
- Keep bulk LOMO metrics as primary for bulk transcriptomics.

Primary sources:

- Arc cell-eval: https://github.com/ArcInstitute/cell-eval
- Arc Virtual Cell Challenge wrap-up:
  https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up
- VCC metrics article:
  https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge

### 6. Virtual-cell model evaluation is moving quickly and remains unstable

Arc STATE, Xaira X-Cell, and related virtual-cell efforts show that perturbation
prediction is now a major frontier. However, model availability, licensing,
input requirements, and benchmark sensitivity vary substantially. scArchon also
highlights that virtual-cell model ranking can be metric-dependent and difficult
to reproduce in the wild.

Implication for v9:

- Prepare adapters for model classes, not one-off fragile integrations.
- Use strong baselines and calibration checks, not just foundation-model names.
- Make the unique claim "spaceflight domain shift" rather than "we ranked the
  newest model this month."

Primary sources:

- Arc STATE: https://arcinstitute.org/news/virtual-cell-model-state
- X-Cell model card: https://huggingface.co/Xaira-Therapeutics/X-Cell
- scArchon Genome Biology 2026:
  https://link.springer.com/article/10.1186/s13059-026-04104-z

### 7. LINCS/L1000 resources remain hypothesis-generation tools

L1000CDS2 supports LINCS signature reversal and mimicking queries. Enrichr
supports enrichment workflows including gene-list submission and library
queries. These are useful for v8 INTERVENE-style triage but do not validate
spaceflight countermeasures.

Implication for v9:

- Keep intervention outputs in a "hypothesis triage" layer.
- Require payload hashes, raw-response checksums, query timestamps, library
  versions, and parsed-output hashes.
- Pair any drug or target claims with safety triage and independent validation
  requirements.

Primary sources:

- L1000CDS2 help: https://maayanlab.cloud/L1000CDS2/help/
- Enrichr: https://maayanlab.cloud/Enrichr/
- CLUE/CMap: https://clue.io/

### 8. Radiation quality is a v9-grade scientific opening

NASA human-spaceflight risk materials continue to emphasize radiation,
isolation, distance from Earth, gravity fields, and hostile/closed environments.
For GeneLab, the important v9 opportunity is not total dose alone. It is the
biological non-equivalence of stressor regimes: low-LET versus high-LET
radiation, radiation plus unloading, time/isolation, sex effects, tissue
specificity, and nonlinear saturation.

Implication for v9:

- Build a radiation-quality and nonlinear-stressor subtrack.
- Treat Mars outputs as falsification/regime-change tests, not point forecasts.
- Use uncertainty, saturation sensitivity, and external validation plans as
  acceptance criteria.

Primary sources:

- NASA five hazards: https://www.nasa.gov/hrp/hazards/
- NASA space radiation overview:
  https://www.nasa.gov/directorates/esdmd/hhp/space-radiation/
- NASA NTRS human-system risk resources:
  https://ntrs.nasa.gov/citations/20210023466

### 9. Provenance should become a first-class benchmark feature

v8 already moved toward manifest-backed claims. v9 should formalize this as a
platform feature: frozen task manifests, source checksums, method cards,
submission validation, RO-Crate-compatible metadata, Hugging Face dataset
bundles, and Zenodo DOI snapshots.

Implication for v9:

- Create one manifest schema shared by tasks, sources, metrics, and runs.
- Export an optional RO-Crate metadata bundle for each frozen release.
- Keep Git for code and small summaries; put public bundles on Hugging Face and
  immutable snapshots on Zenodo.

Primary sources:

- RO-Crate: https://www.researchobject.org/ro-crate/
- FAIR principles: https://www.nature.com/articles/sdata201618
- Hugging Face dataset cards:
  https://huggingface.co/docs/hub/datasets-cards
- Zenodo GitHub integration:
  https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content

## Strategic gap map

| Gap | External state | GeneLab v9 opening |
|---|---|---|
| Spaceflight OOD benchmark | Space data exist, but public AI benchmarks are fragmented | Mission-held-out benchmark with frozen splits |
| Cross-species human bridge | SOMA anchors human spaceflight, but small-N limits model ranking | Human-alignment validation layer |
| Virtual-cell evaluation | cell-eval/OpenProblems cover general single-cell tasks | Spaceflight-specific AnnData subtrack |
| Radiation quality | NASA recognizes radiation risk, but benchmark tasks rarely separate LET/regime | Low-LET vs high-LET and nonlinear stressor tasks |
| Countermeasure claims | LINCS/Enrichr support signature triage | Hypothesis-only intervention module with safety gates |
| Provenance | FAIR/RO-Crate patterns exist, often not enforced in benchmarks | Manifest-validated benchmark releases |

## Recommended v9 thesis

GeneLab Benchmark v9 should claim:

> Existing virtual-cell and omics benchmarks do not test whether models
> generalize across spaceflight missions, tissues, species, and stressor
> regimes. SpaceBio-Bench fills this gap with frozen mission-held-out tasks,
> spaceflight-specific single-cell evaluation, radiation-quality stress tests,
> and manifest-backed provenance.

## Research tasks before implementation

1. Build a source inventory table for all candidate v9 datasets.
2. Build a competitor benchmark matrix: OpenProblems, cell-eval, VCC, scArchon,
   CMap/LINCS, and existing space-omics resources.
3. Define v9 task contracts before writing any model code.
4. Choose one flagship scientific track for the first v9 paper:
   mission-held-out virtual-cell evaluation or radiation-quality stressor
   modeling.
5. Keep gated human data readiness as a protocol track, not a public data
   dependency.

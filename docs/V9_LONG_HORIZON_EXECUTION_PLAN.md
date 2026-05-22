# v9 Long-Horizon Execution Plan

Status: working execution plan
Date: 2026-05-21
Program name: SpaceBio-Bench v9
Companion docs:

- `docs/V9_EXTERNAL_DEEP_RESEARCH_2026_05_21.md`
- `docs/V9_DESIGN_OPTIONS.md`
- `docs/V9_SOURCE_AND_COMPETITOR_MATRIX.md`
- `docs/V9_LONG_RUN_OPERATING_PROTOCOL.md`
- `docs/V9_PUBLIC_BULK_PACKAGE_DESIGN.md`
- `docs/V9_ORGANOID_AND_SPECIES_EXTENSION_REVIEW_2026_05_21.md`
- `docs/V9_MULTISPECIES_FEATURE_STRATEGY.md`
- `docs/v9_hf_dataset_card.md`
- `v9/OPERATING_BACKLOG.md`

## Purpose

This document turns the v9 concept into a long-running work program. It is meant
to let future sessions continue without rediscovering the strategy, rearguing
scope, or accidentally mixing v8 release work with v9 development.

The v9 program should be run like a benchmark platform build, not like a single
analysis sprint.

## North star

SpaceBio-Bench v9 is a public, living, provenance-first benchmark for biological
AI under spaceflight domain shift.

The durable claim:

> Existing omics and virtual-cell benchmarks do not directly test whether models
> generalize across spaceflight missions, tissues, species, modalities,
> perturbations, and stressor regimes. SpaceBio-Bench fills that gap with frozen
> mission-held-out tasks, spaceflight-specific single-cell evaluation,
> radiation-quality stress tests, and manifest-backed provenance.

## Current state at plan creation

Completed v9 scaffold:

- External ecosystem research docs are present.
- `spacebio_bench` package skeleton exists.
- Task-manifest and metric-profile schemas exist.
- `mission_discrimination` metric exists and is tested.
- Legacy bulk LOMO task exporter exists.
- Eight bulk LOMO v9 manifests exist under `v9/task_manifests/`.
- Task registry and index builder exist.
- `v9/task_manifest_index.csv` and `.json` exist.
- Submission validation, local CSV evaluation, and evaluation run-manifest
  writing exist.
- Read-only bulk task loading and fold-level data indexing exist.
- Nearest-centroid baseline runner exists and has been run across all eight
  generated bulk LOMO manifests.
- `v9/reports/nearest_centroid/` contains per-task predictions, metrics, run
  manifests, and a cross-task summary.
- PCA-LR and L2 logistic-regression sklearn baseline runners exist and have
  been run across all eight generated bulk LOMO manifests.
- `v9/reports/bulk_lomo_baseline_summary.csv` contains 24 evaluated
  task/baseline rows across three simple baselines.
- `v9/source_inventory.csv` contains 22 deduplicated public bulk source rows
  with OSDR accession, GLDS prefix, mission, tissue, task ids, privacy class,
  and current checksum status.
- `v9/source_checksum_audit.csv` contains OSDR API file-list evidence and
  checksum-manifest evidence for all 22 public bulk source rows: 22 API-ok
  sources, 39 checksum manifests, 8,439 parsed MD5 entries, and 8,275 payload
  name matches against OSDR listings.
- `v9/datapackage.draft.json` contains a draft Frictionless Data Package
  descriptor for current v9 public bulk metadata, provenance evidence, and
  baseline outputs.
- `docs/v9_hf_dataset_card.md` contains draft Hugging Face-style
  release-facing language with explicit not-frozen and not-payload-verified
  status.
- `docs/V9_ORGANOID_AND_SPECIES_EXTENSION_REVIEW_2026_05_21.md` confirms public
  human neural organoid RNA-seq and non-mouse GeneLab/OSDR resources as credible
  v9 extension candidates, with separate task-family boundaries from current
  mouse bulk LOMO.
- Current generated task manifests and `v9/source_inventory.csv` include
  organism, material, model-system, assay-modality, and feature-namespace
  metadata fields.
- Draft extension source inventories exist for human organoids
  (`v9/human_organoid/source_inventory.draft.csv`) and multispecies pilots
  (`v9/multispecies/source_inventory.draft.csv`).
- A draft `human_organoid_spaceflight` task manifest exists at
  `v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json`, with
  donor, organoid-type, microglia-condition, and single-mission limitations
  recorded before any baseline claims.
- `v9/human_organoid/source_checksum_audit.draft.csv` contains API-ok and
  checksum-manifest-parsed evidence for OSD-863 and OSD-871; payload hashing
  remains pending.
- `v9/human_organoid/sample_table_audit.draft.csv` parses one SampleTable per
  human organoid source, with 19 OSD-863 rows and 23 OSD-871 rows.
- `v9/human_organoid/sample_factors.draft.csv` parses all 42 public organoid
  samples into disease context, Ground/LEO label, microglia condition, and
  organoid type. The draft task manifest now carries sample-count-backed folds,
  while baselines and payload hashes remain pending.
- `v9/human_organoid/expression_matrix_audit.draft.csv` records canonical
  OSDR normalized-count matrix download and sample-column alignment for both
  human organoid sources. OSD-863 has 30,408 feature rows with 19/19 samples
  matched, and OSD-871 has 30,269 feature rows with 23/23 samples matched.
- `v9/human_organoid/reports/nearest_centroid/` contains the first draft-only
  human organoid pilot baseline. It evaluates four sample-count-backed folds
  over 27,986 common human gene features and is explicitly marked
  `pilot_baseline_only_not_leaderboard`.
- `v9/human_organoid/reports/sensitivity/` contains a 20-variant preprocessing
  sensitivity grid, and
  `v9/human_organoid/reports/ORGANOID_BASELINE_ROBUSTNESS_REVIEW.md` records
  the first confounding review. The conservative default remains
  log1p/z-score/2,000 variable genes until donor and library-size diagnostics
  are complete.
- `v9/human_organoid/sample_scale_diagnostics.draft.csv`,
  `v9/human_organoid/group_scale_diagnostics.draft.csv`, and
  `v9/human_organoid/donor_metadata_audit.draft.csv` record the first local
  donor/library-size audit. Current OSDR SampleTables expose only `condition`,
  but GEO GSE259421 series metadata now recovers Subject and iPSC-line
  identifiers for 42/42 public samples. Donor holdouts remain diagnostic-only
  because donor, source, organoid fate, and disease context are not
  independently crossed. Label-associated sample-total/sparsity differences
  keep raw-count sensitivity variants interpretively unsafe.
- `v9/human_organoid/geo_sample_metadata.draft.csv` records GEO-derived sample
  titles, Subject ids, iPSC-line ids, BioSample/SRA accessions, cell type, and
  treatment metadata; these fields are merged into
  `v9/human_organoid/sample_factors.draft.csv`.
- `v9/human_organoid/reports/donor_diagnostics/` evaluates GEO-derived donor
  holdouts as a diagnostic-only fold family. The decision note at
  `v9/human_organoid/reports/ORGANOID_DONOR_AWARE_SPLIT_DECISION.md` keeps
  donor holdouts out of the default benchmark split because donor, source,
  organoid fate, and disease context are coupled.
- `v9/human_organoid/signature_reference_audit.draft.csv` records the first
  DE/signature reference audit. OSDR lists public differential-expression
  tables and contrast-definition CSVs for both OSD-863 and OSD-871; each source
  has 56 parsed contrast pairs, including four direct matched Ground Control
  versus Space Flight contrasts and four reversed matches. The decision note at
  `v9/human_organoid/reports/ORGANOID_SIGNATURE_METRIC_REFERENCE_DECISION.md`
  keeps DE/signature metrics reference-backed but non-primary until a frozen
  contrast/signature contract is added.
- `docs/V9_MULTISPECIES_FEATURE_STRATEGY.md` fixes the initial rule that
  species-native tasks may use species-local gene ids, while cross-species
  bridge tasks should start with pathway/NES-style feature spaces.
- Full repository unit tests pass with v9 tests included.

Important boundary:

- Existing staged v8 beta work remains separate and must not be overwritten.
- Current v9 manifests are alpha scaffolds, not frozen benchmark release
  claims. They preserve legacy source identity with
  `checksum_status=legacy_task_source_unfrozen`.

## Guiding principles

1. Freeze claims before expanding scope.
2. Make every task runnable from a manifest.
3. Make every metric biologically interpretable.
4. Prefer strong baselines before foundation-model adapters.
5. Keep public benchmark tasks independent of gated human data.
6. Keep intervention and Mars-regime outputs hypothesis-only unless separately
   validated.
7. Maintain frozen snapshots and living refreshes as separate tracks.
8. Add one flagship scientific track at a time.

## Program phases

### Phase 0: v8/v9 boundary and repo hygiene

Goal:

Keep v8 beta release work stable while v9 grows in a separate namespace.

Exit criteria:

- v8 beta staged work is either committed, branched, or intentionally preserved.
- v9 files are clearly separated under `spacebio_bench/`, `v9/`, `docs/V9_*`,
  `scripts/*v9*`, and v9-specific tests.
- Full tests pass.
- No v9 artifact is described as a frozen release unless source checksums and
  release metadata are finalized.

Current status:

- Mostly satisfied, except v8 beta is still staged separately.

### Phase 1: Platform spine

Goal:

Turn current v9 scaffolding into a usable local benchmark API.

Deliverables:

- Task registry with query helpers.
- Task manifest validation and index generation.
- Data-loader interface for public bulk LOMO tasks.
- Submission validator for prediction CSVs.
- Evaluation runner for `genelab_minimal`.
- Baseline runner for simple public methods.

Exit criteria:

- One command can evaluate a baseline on at least one bulk LOMO task.
- One command can evaluate all generated bulk LOMO manifests.
- Outputs include predictions, metrics, run manifest, and indexable summary.
- Tests cover successful and failing submission cases.

Current status:

- Satisfied for nearest-centroid, PCA-LR, and L2 logistic-regression bulk
  baselines.
- Source inventory and checksum-manifest audit are complete for the current
  public bulk source set.
- Random forest and later model adapters remain in Phase 2/WS4 after public data
  package boundaries are clear.

### Phase 2: Frozen public bulk benchmark alpha

Goal:

Promote the legacy bulk LOMO tasks into a v9-alpha frozen public benchmark
snapshot.

Deliverables:

- Source checksum audit for all public bulk LOMO inputs.
- Freeze status promoted from `legacy_task_source_unfrozen` to an explicit v9
  alpha source-freeze state.
- Dataset inventory table with OSDR accession, GLDS prefix, tissue, mission,
  privacy class, checksum status, and release target.
- Baseline report for PCA-LR, logistic regression, random forest, and nearest
  centroid.
- Hugging Face dataset-card draft for v9 alpha.
- Zenodo-ready release manifest draft.

Exit criteria:

- All public bulk tasks regenerate from documented inputs.
- All baseline outputs can be reproduced from a clean checkout plus public data
  bundle.
- `v9-alpha` result language is limited to benchmark infrastructure and bulk
  baseline reproduction.

Current status:

- Source checksum-manifest evidence exists for all 22 source rows.
- Payload download and local hash verification is not complete; therefore
  current audit rows deliberately keep `freeze_ready=false`.
- Public bulk package design is complete as a draft boundary and descriptor.
- Dataset-card draft language is complete and stays explicit about draft status
  and payload verification gaps.
- The next step is RO-Crate export design, so release metadata can link sources,
  tasks, package resources, workflow runs, and checksum evidence.

### Phase 3: Mission-held-out single-cell track

Goal:

Add the first flagship scientific extension: spaceflight scRNA-seq evaluation
aligned with cell-eval ideas where biologically valid.

Deliverables:

- RRRM-1/RRRM-2 task inventory.
- AnnData task-card schema extension.
- `genelab_sc` evaluator with DE overlap, direction match, mission
  discrimination, and state-label metrics where applicable.
- Baselines: per-cell-type DE, logistic regression, nearest centroid, simple
  embedding baselines.
- Optional model adapters only when installation, licensing, and inference are
  reproducible.

Exit criteria:

- At least one scRNA-seq task runs end to end with public data.
- Bulk metrics and single-cell metrics remain clearly separated.
- The first v9 science draft can claim a spaceflight OOD single-cell evaluation
  track without claiming broad virtual-cell superiority.

### Phase 4: Radiation-quality and nonlinear-stressor track

Goal:

Add the second flagship scientific extension: radiation quality and nonlinear
combined-stressor benchmark tasks.

Deliverables:

- v8 DECOMPOSE source inventory promoted into v9 stressor task manifests.
- Low-LET versus high-LET task definitions.
- HLU/radiation/time interaction task definitions.
- Saturation and uncertainty metrics under `stressor_regime`.
- Falsification language for Mars-regime projections.

Exit criteria:

- At least one radiation-quality task runs end to end.
- Outputs distinguish sign stability, saturation sensitivity, and uncertainty.
- No result is described as a Mars point prediction.

### Phase 4A: Human organoid and multispecies extension

Goal:

Add public non-mouse extension tracks after the v9 alpha public bulk spine is
stable.

Deliverables:

- Manifest/source schema fields for organism, material type, model system,
  assay modality, feature namespace, and orthology strategy.
- Human neural organoid RNA-seq source inventory for OSD-863/OSD-871 and GEO
  GSE259421.
- One draft `human_organoid_spaceflight` task card with donor, organoid type,
  and microglia-condition blocking explicitly handled.
- Canonical normalized-count matrices aligned to parsed organoid sample factors
  before any baseline is reported.
- One transparent draft-only human organoid baseline with conservative
  single-mission limitations.
- Non-mouse source inventory for existing Drosophila and Arabidopsis pilot
  assets.
- Ortholog/pathway feature strategy for cross-species tasks.

Exit criteria:

- Existing mouse bulk manifests still validate unchanged.
- Organoid and non-mouse sources have task-family boundaries and no claims are
  made as current v9-alpha results.
- A first organoid or multispecies baseline is run only after split rules,
  source integrity checks, and matrix/sample alignment evidence are written.

### Phase 5: Human bridge and gated-data readiness

Goal:

Prepare human multi-omics and commercial astronaut data integration without
making restricted data a public benchmark dependency.

Deliverables:

- SOMA/SpaceOmicsBench bridge task contracts.
- Human validation summary schema.
- EXPAND-ready gated-track protocol.
- Privacy and allowed-output classes.
- Hidden-evaluation protocol draft.

Exit criteria:

- Public tasks run without restricted data.
- Gated track can be documented for future access-controlled evaluation.
- Human bridge claims are framed as validation/alignment evidence, not
  large-N model ranking.

### Phase 6: Release, manuscript, and community layer

Goal:

Package v9 as a durable benchmark resource.

Deliverables:

- `spacebio_bench` API docs.
- Benchmark paper outline and first full draft.
- Dataset card.
- Model/submission card template.
- RO-Crate-compatible metadata export.
- Hugging Face dataset bundle.
- Zenodo DOI snapshot.
- Community submission guide.

Exit criteria:

- External users can understand task scope, download data, run baselines,
  validate submissions, and cite the release.
- Frozen and living tracks are clearly separated.
- Manuscript claims are backed by manifests, source records, outputs, and tests.

## Workstreams

### WS1: Source and provenance

Owner scope:

- OSDR accessions, GLDS prefixes, source checksums, API payloads, task manifests,
  and freeze metadata.

Near-term tasks:

- Treat the v9 source inventory table as complete for generated bulk LOMO
  manifests.
- Treat `v9/source_checksum_audit.csv/json` as complete for OSDR API and
  checksum-manifest evidence.
- Treat `v9/datapackage.draft.json` as complete for draft descriptor planning.
- Treat `docs/v9_hf_dataset_card.md` as complete for draft dataset-card
  language.
- Keep payload download and local checksum verification as a separate freeze
  step before promoting manifest checksum metadata.
- Decide how public bulk matrices are packaged outside Git.

Definition of done:

- Every task source has accession, URL, checksum state, privacy class, and
  release target.

### WS2: Benchmark package API

Owner scope:

- `spacebio_bench` loaders, registry, metrics, evaluators, submissions, reports.

Near-term tasks:

- Harden the bulk loader abstraction beyond legacy path adapters.
- Keep prediction/submission validation aligned with baseline output metadata.
- Add baseline-level report aggregation helpers as additional methods arrive.

Definition of done:

- A user can run `load_task`, `evaluate`, and `write_report` without reading
  internal repo scripts.

### WS3: Metrics

Owner scope:

- Bulk, single-cell, stressor, and calibration metrics.

Near-term tasks:

- Implement macro-F1, balanced accuracy, AUROC, and calibration wrappers.
- Integrate `mission_discrimination` into evaluator outputs.
- Define exact `genelab_sc` metric formulas before implementing them.

Definition of done:

- Each metric has a code implementation, a test, and a biological
  interpretation.

### WS4: Baselines

Owner scope:

- Simple reproducible methods before foundation models.

Near-term tasks:

- Treat nearest-centroid baseline as complete for bulk LOMO alpha scaffold.
- Treat logistic regression/PCA-LR adapter as complete for the current bulk LOMO
  alpha scaffold.
- Add random forest only after source inventory/checksum work, so stronger
  baseline expansion does not outrun provenance hardening.

Definition of done:

- Every public v9-alpha task has at least two simple baselines.

### WS5: Single-cell flagship

Owner scope:

- RRRM-1/RRRM-2 task cards, AnnData loader, DE metrics, optional adapters.

Near-term tasks:

- Inventory existing scRNA-seq assets and scripts.
- Draft one RRRM task manifest.
- Decide whether cell-eval is a dependency, optional extra, or reference-only.

Definition of done:

- One scRNA-seq mission-held-out task runs with `genelab_sc`.

### WS6: Radiation/stressor flagship

Owner scope:

- Radiation quality, HLU, time/isolation, nonlinear response, uncertainty.

Near-term tasks:

- Map v8 DECOMPOSE outputs to v9 stressor task manifests.
- Define `stressor_regime` metrics precisely.
- Add held-out analog validation checks.

Definition of done:

- One stressor task runs with saturation/uncertainty outputs and conservative
  language.

### WS7: Packaging and release

Owner scope:

- Hugging Face, Zenodo, RO-Crate, dataset cards, release notes.

Near-term tasks:

- Draft v9 dataset card.
- Add release manifest schema.
- Add RO-Crate exporter design.

Definition of done:

- A frozen v9-alpha snapshot can be archived and cited.

### WS8: Manuscript and positioning

Owner scope:

- Platform manuscript, science flagship manuscript, claims table, related work.

Near-term tasks:

- Draft SpaceBio-Bench platform outline.
- Add claim-to-artifact map.
- Add competitor matrix to manuscript notes.

Definition of done:

- Every manuscript claim maps to task, metric, output, and source manifest.

### WS9: Organoid and multispecies extension

Owner scope:

- Public human organoid data, non-mouse OSDR/GeneLab species, ortholog/pathway
  feature strategies, and task-family boundaries beyond mouse bulk LOMO.

Near-term tasks:

- Treat the organoid/species extension review as complete planning evidence.
- Add source/task schema support before importing OSD-863, OSD-871, OSD-207,
  OSD-37, or OSD-120 into v9 manifests.
- Treat OSD-863/OSD-871 sample factors, canonical normalized-count matrices, the
  first nearest-centroid pilot baseline, and sensitivity review as draft-only
  extension artifacts, while keeping payload freeze and leaderboard language out
  of scope.
- Keep organoid RNA-seq, organoid proteomics, species-specific flight/ground,
  and cross-species transfer as separate task families.

Definition of done:

- Public non-mouse task candidates have explicit organism, material, modality,
  feature namespace, split strategy, and claim boundary metadata.

## Decision gates

| Gate | Question | Required evidence | Decision |
|---|---|---|---|
| G0 | Is v8 safe from v9 churn? | Clean v8/v9 file separation and tests | Continue v9 only after boundary is clear |
| G1 | Is the platform spine runnable? | One task end-to-end with baseline metrics | Move from planning scaffold to alpha implementation |
| G2 | Are public bulk tasks frozen? | Source checksums and clean regeneration | Publish v9-alpha bulk snapshot |
| G3 | Which flagship comes first? | scRNA-seq data readiness vs radiation task readiness | Choose single-cell or radiation first |
| G4 | Can external users reproduce it? | docs, package API, dataset bundle, tests | Prepare public release |
| G5 | Are manuscript claims safe? | Claim-to-artifact map and no overclaim language | Draft submission package |

## Cadence

Recommended work rhythm:

- Use `docs/V9_LONG_RUN_OPERATING_PROTOCOL.md` as the default operating mode
  when the user asks to continue v9 work.
- Treat a continue request as a request for a 30-60 minute uninterrupted block
  unless the user explicitly asks for a short answer or status only.
- Every long session starts by reading `v9/OPERATING_BACKLOG.md` and
  `git status --short --branch`.
- Every long block should combine local build work, external primary-source
  research, and verification/doc continuity where possible.
- Every implementation session should end with tests and a short status note.
- Every major result must add or update a manifest before being described in
  manuscript language.
- Every new task family gets one minimal happy-path test and one failure-mode
  test before broad expansion.
- Every week or equivalent work block should update the backlog statuses.

## Standard validation commands

```bash
python scripts/export_v9_task_manifests.py
python scripts/build_v9_task_index.py
python scripts/build_v9_task_data_index.py
python scripts/run_v9_nearest_centroid.py --output-dir v9/reports/nearest_centroid
python scripts/run_v9_sklearn_baselines.py --output-dir v9/reports/sklearn_baselines
python scripts/build_v9_baseline_summary.py
python scripts/build_v9_source_inventory.py
python scripts/audit_v9_source_checksums.py
python scripts/build_v9_datapackage_draft.py
python scripts/evaluate_v9_submission.py --task-manifest v9/task_manifests/A2_gastrocnemius_bulk_lomo.json --submission path/to/predictions.csv --output-dir path/to/report
python -m unittest tests/test_v9_spacebio_bench.py
python -m unittest discover -s tests
```

Use heavier data-generation or HPC commands only after the manifest-level
contracts are stable.

## Long-term risks

| Risk | Why it matters | Guardrail |
|---|---|---|
| Scope creep | v9 spans platform, single-cell, stressor, human, intervention | Platform spine plus one flagship at a time |
| v8/v9 contamination | v8 beta release still has staged work | Keep namespaces and commits separate |
| Data drift | OSDR and external APIs can change | Freeze snapshots and keep living track separate |
| Weak baselines | Foundation models can distract from benchmark validity | Require simple baselines first |
| Metric mismatch | cell-eval metrics do not fit bulk LOMO | Metric profiles by task family |
| Human privacy | SOMA/EXPAND data can be sensitive | Public core cannot depend on gated data |
| Countermeasure overclaim | LINCS hits are not clinical validation | Hypothesis-only intervention language |
| Mars overclaim | Linear projections are stress tests | Regime/falsification language only |
| Species conflation | Mouse tissue, human organoids, plants, insects, and microbes do not share one raw feature space | Separate task families plus ortholog/pathway feature namespaces |
| Organoid overclaim | Human organoid datasets are small and often single-mission | Pilot/task-card language until source checks, splits, and baselines are complete |

## Preferred next build sequence

1. Add submission validator.
2. Add evaluator for prediction CSVs.
3. Add minimal bulk data loader or adapter over existing processed task outputs.
4. Add checksum audit over the generated source inventory.
5. Add public bulk package design and future `datapackage.json` resource shape.
6. Add dataset-card draft with conservative draft/freeze language.
7. Add RO-Crate export design.
8. Add release manifest scaffold.
9. Add schema design for organism/material/model-system fields before non-mouse
   expansion.
10. Recover donor/iPSC-line metadata from external GEO or publication sources
   before any manuscript-facing organoid interpretation.
11. Decide random forest versus payload checksum verification as the next bulk
   step.
12. Decide single-cell, radiation, or multispecies as the next flagship branch
   after the organoid pilot baseline reaches a clean draft checkpoint.

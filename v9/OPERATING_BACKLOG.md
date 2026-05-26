# v9 Operating Backlog

Status: living backlog for SpaceBio-Bench v9
Last updated: 2026-05-26

Use this file at the start of each long v9 work session. Keep tasks small enough
to finish or leave in a clear handoff state, but organize them into long-running
workstreams.

Long-run rule:

- When the user asks to continue, follow
  `docs/V9_LONG_RUN_OPERATING_PROTOCOL.md`.
- Default to one 30-60 minute uninterrupted block, not a 1-3 minute checkpoint.
- Pair local implementation with external primary-source research whenever the
  next decision depends on standards, APIs, or release conventions.
- Send final answers only at real artifact checkpoints unless the user asks for
  a short status update.

## Current checkpoint

Done:

- External v9 research docs.
- v9 design options and source/competitor matrix.
- `spacebio_bench` skeleton.
- Minimal task-manifest validator.
- Metric profiles.
- `mission_discrimination` metric.
- Legacy bulk LOMO manifest exporter.
- Eight generated bulk LOMO task manifests.
- Task registry and manifest index.
- Submission validator for prediction CSVs.
- `genelab_minimal` prediction evaluator.
- Evaluation report writer with `metrics.json` and `run_manifest.json`.
- CLI evaluation entrypoint: `scripts/evaluate_v9_submission.py`.
- Read-only bulk task loader and fold-level data index builder.
- Nearest-centroid bulk baseline package runner and CLI.
- All eight generated bulk LOMO manifests evaluated with nearest-centroid under
  `v9/reports/nearest_centroid/`.
- PCA-LR and L2 logistic-regression sklearn baseline runner and CLI.
- All eight generated bulk LOMO manifests evaluated with PCA-LR and L2 logistic
  regression under `v9/reports/sklearn_baselines/`.
- Cross-baseline bulk report generated at
  `v9/reports/bulk_lomo_baseline_summary.csv`.
- Source inventory builder and CLI.
- Source-level v9 inventory generated at `v9/source_inventory.csv` with 22
  public bulk source rows.
- Source checksum audit generated at `v9/source_checksum_audit.csv` with OSDR
  API evidence and parsed checksum-manifest evidence for all 22 public bulk
  source rows.
- Public bulk package design drafted and `v9/datapackage.draft.json` generated.
- Draft Hugging Face-style dataset card written at
  `docs/v9_hf_dataset_card.md`.
- Organoid and non-mouse species extension review written at
  `docs/V9_ORGANOID_AND_SPECIES_EXTENSION_REVIEW_2026_05_21.md`.
- Organism/material/model-system schema fields added to generated mouse bulk
  manifests, task index, and source inventory.
- Draft human organoid source inventory generated at
  `v9/human_organoid/source_inventory.draft.csv`.
- Draft human organoid task manifest generated at
  `v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json`.
- Draft human organoid checksum audit generated at
  `v9/human_organoid/source_checksum_audit.draft.csv`.
- Draft human organoid sample-table audit generated at
  `v9/human_organoid/sample_table_audit.draft.csv`.
- Draft human organoid sample-factor table generated at
  `v9/human_organoid/sample_factors.draft.csv`.
- Draft human organoid normalized-count expression matrix audit generated at
  `v9/human_organoid/expression_matrix_audit.draft.csv`, with canonical
  OSD-863 and OSD-871 matrices downloaded under `v9/human_organoid/matrices/`.
- Draft human organoid nearest-centroid pilot baseline generated under
  `v9/human_organoid/reports/nearest_centroid/`.
- Draft human organoid sensitivity grid and robustness review generated under
  `v9/human_organoid/reports/sensitivity/` and
  `v9/human_organoid/reports/ORGANOID_BASELINE_ROBUSTNESS_REVIEW.md`.
- Draft human organoid donor metadata and library-size diagnostics generated at
  `v9/human_organoid/donor_metadata_audit.draft.csv`,
  `v9/human_organoid/sample_scale_diagnostics.draft.csv`, and
  `v9/human_organoid/group_scale_diagnostics.draft.csv`.
- GEO GSE259421 sample metadata parsed at
  `v9/human_organoid/geo_sample_metadata.draft.csv`, recovering Subject and
  iPSC-line identifiers for 42/42 public human organoid samples and merging
  them into `v9/human_organoid/sample_factors.draft.csv`.
- Donor-aware holdout evaluation generated under
  `v9/human_organoid/reports/donor_diagnostics/`, with donor folds kept
  diagnostic-only in
  `v9/human_organoid/reports/ORGANOID_DONOR_AWARE_SPLIT_DECISION.md`.
- Human organoid DE/signature reference audit generated at
  `v9/human_organoid/signature_reference_audit.draft.csv`, showing public OSDR
  DE tables and contrast-definition CSVs for OSD-863/OSD-871, with the metric
  policy recorded in
  `v9/human_organoid/reports/ORGANOID_SIGNATURE_METRIC_REFERENCE_DECISION.md`.
- Human organoid derived DE reference generated at
  `v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz`
  with 2 sources, 8 direct contrasts, 242,708 gene/contrast rows, 2,368
  FDR<=0.05 rows, and a response-signature contract recorded in
  `v9/human_organoid/reports/ORGANOID_DE_REFERENCE_CONTRACT.md`.
- Human organoid response-signature diagnostic scorer implemented in
  `spacebio_bench/signature_metrics.py` and wired into the evaluator. It
  validates `response_signature.csv`, joins to the derived DE reference, and
  computes `de_direction_match` plus `signature_rank_correlation` when a valid
  artifact is supplied. Scorer policy is recorded in
  `v9/human_organoid/reports/ORGANOID_RESPONSE_SIGNATURE_SCORER.md`.
- Human organoid response-signature smoke report generated at
  `v9/human_organoid/reports/response_signature_smoke/`, exercising the scorer
  against the real derived DE reference with a 40-row mirrored fixture and
  explicit non-model claim boundary.
- Human organoid model-produced response-signature adapter design written at
  `v9/human_organoid/reports/ORGANOID_RESPONSE_SIGNATURE_ADAPTER_DESIGN.md`,
  selecting a non-leaky source-transfer empirical signature baseline as the
  first real adapter path.
- Human organoid source-transfer response-signature adapter implemented at
  `v9/human_organoid/reports/source_transfer_signature/`, producing 223,888
  response-signature rows and diagnostic scores
  `de_direction_match=0.7706734867860188` and
  `signature_rank_correlation=0.1760078660242601`.
- Human organoid source-transfer diagnostic review written at
  `v9/human_organoid/reports/ORGANOID_SOURCE_TRANSFER_SIGNATURE_REVIEW.md`,
  keeping the source-transfer signature as the first conservative diagnostic
  baseline, showing direction agreement above trivial sign baselines, and
  selecting a per-condition source-transfer adapter design as the next block.
- Human organoid per-condition source-transfer adapter design written at
  `v9/human_organoid/reports/ORGANOID_PER_CONDITION_SIGNATURE_ADAPTER_DESIGN.md`,
  selecting a microglia-matched source-transfer empirical signature as the next
  comparative diagnostic and deferring disease-matched signatures to a partial
  coverage follow-up.
- Human organoid microglia-matched source-transfer response-signature adapter
  implemented at
  `v9/human_organoid/reports/microglia_source_transfer_signature/`, producing
  223,888 response-signature rows and diagnostic scores
  `de_direction_match=0.7889125799573561` and
  `signature_rank_correlation=0.1500722829461316`. Compared with the global
  source-transfer diagnostic, direction agreement improves by 0.018239 while
  rank correlation decreases by 0.025936.
- Human organoid source-transfer adapter comparison review written at
  `v9/human_organoid/reports/ORGANOID_SOURCE_TRANSFER_ADAPTER_COMPARISON_REVIEW.md`,
  keeping the global adapter as the first conservative diagnostic, keeping the
  microglia-matched adapter as a secondary condition-sensitivity diagnostic, and
  selecting a partial shared-control disease+microglia adapter design as the
  next block.
- Human organoid shared-control disease+microglia adapter design written at
  `v9/human_organoid/reports/ORGANOID_SHARED_CONTROL_SIGNATURE_ADAPTER_DESIGN.md`,
  limiting disease+microglia matching to four shared `no_known_diseases`
  contrasts and requiring explicit skipped-contrast metadata for PPMS and
  sporadic Parkinson disease contrasts.
- Human organoid partial shared-control source-transfer response-signature
  adapter implemented at
  `v9/human_organoid/reports/shared_control_source_transfer_signature/`,
  producing 111,944 response-signature rows for four shared-control contrasts,
  skipping four disease-specific contrasts, and reporting partial-coverage
  `de_direction_match=0.5899419729206963` plus
  `signature_rank_correlation=0.03566701395309356`.
- Human organoid classifier feature-effect contract written at
  `v9/human_organoid/reports/ORGANOID_CLASSIFIER_FEATURE_EFFECT_CONTRACT.md`,
  separating classifier-derived model effects from log2FC
  `response_signature.csv` artifacts and selecting an L2 logistic gene-space
  feature-effect pilot as the next implementation path.
- Human organoid L2 logistic feature-effect pilot implemented at
  `v9/human_organoid/reports/logistic_feature_effect/`, producing 16,000
  feature-effect rows and diagnostic scores
  `feature_effect_direction_match=0.6078431372549019` and
  `feature_effect_rank_correlation=0.08672800238082004`, with top-k DE overlap
  1/50, 5/100, 10/250, and 14/500.
- Human organoid feature-effect diagnostic review written at
  `v9/human_organoid/reports/ORGANOID_FEATURE_EFFECT_DIAGNOSTIC_REVIEW.md`,
  keeping `feature_effect.csv` as an optional secondary diagnostic and selecting
  top-k/null calibration before PCA-LR reconstructed-weight work.
- Human organoid feature-effect top-k null calibration implemented and reviewed
  at
  `v9/human_organoid/reports/ORGANOID_FEATURE_EFFECT_NULL_CALIBRATION_REVIEW.md`,
  adding expected overlap, enrichment, and exact hypergeometric upper-tail
  p-values to `feature_effect_top_k_de_overlap`. Regenerated L2 logistic
  feature-effect top-k diagnostics show aggregate enrichment at K=100, K=250,
  and K=500, but heterogeneous per-contrast behavior.
- Human organoid PCA-LR reconstructed feature-effect design written at
  `v9/human_organoid/reports/ORGANOID_PCA_LR_FEATURE_EFFECT_DESIGN.md`,
  defining the `pca.components_.T @ logistic.coef_[0]` reconstruction in
  train-standardized gene space and approving one constrained secondary
  diagnostic pilot.
- Human organoid PCA-LR reconstructed feature-effect pilot implemented at
  `v9/human_organoid/reports/pca_lr_feature_effect/` and reviewed at
  `v9/human_organoid/reports/ORGANOID_PCA_LR_FEATURE_EFFECT_REVIEW.md`.
  The report emits 16,000 feature-effect rows, matches L2 logistic aggregate
  top-k enrichment exactly, and is slightly weaker in sign/rank diagnostics, so
  it is kept as an optional negative comparison.
- Human organoid diagnostic family consolidated at
  `v9/human_organoid/reports/ORGANOID_DIAGNOSTIC_CONSOLIDATION_AND_RELEASE_BOUNDARY.md`,
  marking task/provenance, derived DE reference, global source-transfer response
  signatures, and calibrated L2 logistic feature effects as the main draft-alpha
  organoid surface while keeping the other diagnostics secondary, exploratory,
  or plumbing-only.
- Draft multispecies source inventory generated at
  `v9/multispecies/source_inventory.draft.csv`.
- Multispecies feature namespace strategy written at
  `docs/V9_MULTISPECIES_FEATURE_STRATEGY.md`.
- Multispecies candidate source deep audit written at
  `v9/multispecies/reports/MULTISPECIES_CANDIDATE_SOURCE_DEEP_AUDIT.md`,
  selecting OSD-37 and OSD-207 as go-after-source-audit species-native task
  candidates and deferring OSD-120 to a light/genotype interaction task.
- Multispecies sample-factor table generated at
  `v9/multispecies/sample_factors.draft.csv`, parsing 124 local samples across
  OSD-207, OSD-37, and OSD-120 into Ground/LEO labels, condition strata,
  genotype/ecotype factors, and OSD-120 light-treatment factors.
- Multispecies expression-matrix audit generated at
  `v9/multispecies/expression_matrix_audit.draft.csv`, confirming local
  normalized-count matrix columns align with parsed sample factors for all
  three draft multispecies candidates.
- Multispecies OSDR checksum audit generated at
  `v9/multispecies/source_checksum_audit.draft.csv`, with all three draft
  sources returning API-ok and checksum-manifest-parsed evidence.
- Multispecies checksum/local payload review written at
  `v9/multispecies/reports/MULTISPECIES_CHECKSUM_AND_LOCAL_PAYLOAD_AUDIT.md`,
  confirming the six local SampleTable and normalized-count matrix files used
  by the current scaffold match OSDR processed MD5 manifest entries.
- Draft species-native multispecies task manifests generated under
  `v9/multispecies/task_manifests/` for OSD-37 Arabidopsis and OSD-207
  Drosophila, with an index at
  `v9/multispecies/task_manifest_index.draft.csv`.
- Multispecies species-native manifest design note written at
  `v9/multispecies/reports/MULTISPECIES_SPECIES_NATIVE_TASK_MANIFEST_DESIGN.md`,
  keeping OSD-120 deferred to a separate light/genotype interaction-task
  design.
- Multispecies read-only loader and nearest-centroid draft baseline implemented
  for OSD-37 and OSD-207 species-native manifests.
- Multispecies nearest-centroid baseline outputs generated under
  `v9/multispecies/reports/nearest_centroid/`, with summary metrics at
  `v9/multispecies/reports/nearest_centroid/multispecies_baseline_summary.csv`.
- Multispecies baseline feasibility review written at
  `v9/multispecies/reports/MULTISPECIES_BASELINE_FEASIBILITY_REVIEW.md`,
  marking OSD-37 as the cleaner first plant feasibility example and OSD-207 as
  valid but condition-stratum heterogeneous.
- Multispecies baseline sensitivity grid generated under
  `v9/multispecies/reports/sensitivity/`, covering 40 evaluated rows across
  OSD-37 and OSD-207.
- Multispecies baseline sensitivity review written at
  `v9/multispecies/reports/MULTISPECIES_BASELINE_SENSITIVITY_REVIEW.md`,
  keeping the conservative default at `log1p` + train-fold `zscore` + top
  2,000 variable genes and confirming `w1118_KCNQ370` as OSD-207's recurring
  fragile stratum.
- OSD-120 interaction-task design written at
  `v9/multispecies/reports/OSD120_INTERACTION_TASK_DESIGN.md`, selecting a
  separate `multispecies_light_interaction_spaceflight` task family rather than
  merging OSD-120 into the simpler OSD-37 plant task.
- OSD-120 interaction-task manifest generated at
  `v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json`,
  with a draft interaction index at
  `v9/multispecies/interaction_task_manifest_index.draft.csv`.
- OSD-120 interaction nearest-centroid baseline generated under
  `v9/multispecies/reports/interaction_nearest_centroid/`, with separate
  primary genotype/ecotype, secondary light-treatment, and diagnostic
  condition-stratum fold-family reports.
- OSD-120 interaction baseline feasibility review written at
  `v9/multispecies/reports/OSD120_INTERACTION_BASELINE_FEASIBILITY_REVIEW.md`,
  keeping the baseline draft-only and selecting an OSD-120 sensitivity grid as
  the next multispecies block.
- OSD-120 interaction sensitivity grid generated under
  `v9/multispecies/reports/interaction_sensitivity/`, covering 60 evaluated
  rows across 20 preprocessing variants and three fold families.
- OSD-120 interaction sensitivity review written at
  `v9/multispecies/reports/OSD120_INTERACTION_SENSITIVITY_REVIEW.md`, keeping
  `log1p` + train-fold `zscore` + top 2,000 genes as the conservative default
  while recording repeated fragile OSD-120 strata.
- OSD-120 interaction fold-detail aggregation generated at
  `v9/multispecies/reports/interaction_sensitivity/fold_detail_summary.csv`
  and `.json`, exposing 220 fold-level sensitivity rows with variant id,
  held-out value, balanced accuracy, within-variant rank, and fragility flags.
- OSD-120 interaction L2 logistic baseline design and implementation completed,
  with design/review notes at
  `v9/multispecies/reports/OSD120_INTERACTION_LOGISTIC_BASELINE_DESIGN.md` and
  `v9/multispecies/reports/OSD120_INTERACTION_LOGISTIC_BASELINE_REVIEW.md`.
- OSD-120 interaction L2 logistic reports generated under
  `v9/multispecies/reports/interaction_logistic_l2/`, showing higher aggregate
  balanced accuracy than nearest centroid but larger diagnostic
  condition-stratum delta.
- OSD-120 interaction L2 logistic fold-detail comparison generated at
  `v9/multispecies/reports/interaction_logistic_l2/fold_detail_summary.csv`
  and `fold_detail_comparison_vs_nearest_centroid.csv`; across 11 held-out
  folds, logistic improves 8, ties 2, and worsens only
  `Col.0.PhyD|Dark.Treatment`.
- OSD-120 interaction L2 logistic sensitivity generated under
  `v9/multispecies/reports/interaction_logistic_l2_sensitivity/`, covering six
  top-gene/C variants and showing that the `Col.0.PhyD|Dark.Treatment` failure
  is feature-count sensitive rather than regularization-sensitive.
- OSD-120 interaction logistic feature-set audit generated under
  `v9/multispecies/reports/interaction_logistic_feature_audit/`, showing that
  top-2,000 extra features strongly change coefficient rankings across the
  fragile OSD-120 folds.
- OSD-120 sparse feature-stability branch design written at
  `v9/multispecies/reports/OSD120_INTERACTION_SPARSE_BRANCH_DESIGN.md`,
  selecting a draft-only L1 logistic pilot before elastic-net or opaque models.
- OSD-120 sparse L1 logistic pilot generated under
  `v9/multispecies/reports/interaction_logistic_sparse_l1/`, with
  `tvg2000_log1p_zscore_l1_c1` clearing the fragile-fold gates and worsening
  zero nearest-centroid fold rows.
- OSD-120 sparse L1 stability audit generated under
  `v9/multispecies/reports/interaction_logistic_sparse_l1_stability/`, keeping
  `C=1.0` as the performance-leading sparse comparator and `C=0.3` as a compact
  light-treatment stability comparator.
- OSD-120 interaction baseline ladder, diagnostic candidate package,
  figure/table package, paired comparator, diagnostic artifact manifest, and
  release-readiness gap audit generated under
  `v9/multispecies/reports/`. The release-readiness audit records 5 pass
  items, 3 public-alpha blockers, 3 broader-release needs-work items, and 1
  acceptable draft limitation.
- OSD-120 diagnostic-required payload freeze manifest generated under
  `v9/multispecies/reports/interaction_payload_freeze_manifest/`. It parses 533
  OSDR processed MD5 entries, confirms the two diagnostic-required payloads
  match local MD5 hashes, and keeps the remaining 531 processed payloads outside
  the current diagnostic freeze scope.
- OSD-120 diagnostic public-alpha card draft generated under
  `v9/multispecies/reports/interaction_public_alpha_card/`. It records source
  scope, payload-freeze boundary, diagnostic result surface, seven allowed
  claims, five disallowed claims, inspectable files, external context links,
  and remaining release work.
- OSD-120 rebuild gate and environment lock generated under
  `v9/multispecies/reports/interaction_rebuild_gate/`. It records one
  preflight command, 8 script-backed packaging steps, 40 hashed outputs, Python
  and package versions, and no missing rebuild inputs.
- OSD-120 public metadata package skeleton generated under
  `v9/multispecies/reports/interaction_public_metadata_package/`. It separates
  one public-now diagnostic metadata draft from three not-public-now release
  targets, records 20 metadata fields, and keeps 5 DOI/creator/license-style
  placeholders unresolved by design.
- OSD-120 RO-Crate/citation freeze scaffold generated under
  `v9/multispecies/reports/interaction_ro_crate_citation_scaffold/`. It emits
  `ro-crate-metadata.draft.json`, `datapackage.draft.json`, 13 validation
  checks, and 11 citation-freeze checklist items while keeping archive blockers
  explicit.
- OSD-120 archive identifier/license decision gate generated under
  `v9/multispecies/reports/interaction_archive_decision_gate/`. It records six
  archive identifier options, five license/rights components, six
  creator/contributor components, six external guidance anchors, and selects
  no local archive identifier for the current diagnostic draft.
- OSD-120 release-owner citation metadata fill scaffold generated under
  `v9/multispecies/reports/interaction_citation_metadata_fill/`. It records 16
  owner-fill fields, keeps four current-draft values retained, leaves 12
  release blockers explicit, and emits only a descriptor patch preview with no
  RO-Crate/Data Package mutation.
- OSD-120 archive-release deferral/application guard generated under
  `v9/multispecies/reports/interaction_archive_release_deferral_guard/`. It
  records 11 guard checks, 9 blockers, 2 pass checks, 9 deferral/action rows,
  and blocks descriptor mutation while owner metadata remains absent.
- OSD-120 diagnostic metadata release-note closeout generated under
  `v9/multispecies/reports/interaction_diagnostic_metadata_release_note/`. It
  records 3 allowed diagnostic metadata claims, 3 prohibited release claims, 6
  owner metadata retry items, and a no-archive/no-leaderboard public note.
- V9 recenter decision generated under `v9/reports/recenter_decision/` and
  `docs/V9_REFOCUS_001_POST_OSD120_RECENTER_DECISION.md`. It selects public
  bulk alpha readiness as the next active lane and defers the first
  single-cell flagship until after a bulk alpha freeze-gap matrix.
- Public bulk alpha freeze-gap matrix generated under
  `v9/reports/public_bulk_alpha_gap_matrix/` and
  `docs/V9_PUBLIC_BULK_ALPHA_FREEZE_GAP_MATRIX_REVIEW.md`. It records 6 pass
  rows, 2 blockers, 2 needs-update rows, 3 allowed current claims, and 3
  prohibited release claims.
- Purpose-drift audit written at
  `docs/V9_PURPOSE_DRIFT_AUDIT_2026_05_26.md`. It concludes that v9 remains on
  mission but has a moderate sequencing risk: OSD-120 is aligned as a
  Phase 4A/provenance case study, but the branch should close after
  V9-MULTI-035 unless owner metadata appears. V9-MULTI-035 has now closed that
  branch and hands off to V9-REFOCUS-001.
- v9 unit tests integrated into full repository tests.

Current test baseline:

```bash
python -m unittest discover -s tests
```

Expected status after V9-BULK-ALPHA-001 public bulk alpha freeze-gap matrix:

- 220 tests passing.

Known boundary:

- v8 beta work is staged separately.
- `submissions/` is unrelated untracked work and should not be touched unless
  explicitly requested.

## Now: platform spine hardening

### V9-NOW-001: Submission validator

Status: complete as of 2026-05-21.

Goal:

Validate prediction CSVs against a task manifest.

Files likely touched:

- `spacebio_bench/submissions/validate.py`
- `tests/test_v9_spacebio_bench.py`

Requirements:

- Check required columns from manifest output schema.
- Validate non-empty rows.
- Validate labels are known when label domain is provided.
- Allow optional `flight_probability` and `embedding_*` columns.
- Return structured validation report, not just boolean.

Done when:

- Tests cover valid submission, missing required column, empty file, and unknown
  task id or manifest mismatch.

### V9-NOW-002: Evaluator for prediction CSVs

Status: complete as of 2026-05-21 for lightweight local CSV evaluation.

Goal:

Compute `genelab_minimal` metrics from validated predictions.

Files likely touched:

- `spacebio_bench/evaluate.py`
- `spacebio_bench/metrics/`
- `tests/test_v9_spacebio_bench.py`

Requirements:

- Macro-F1.
- Balanced accuracy.
- AUROC when probabilities exist.
- Calibration error when probabilities exist.
- Mission discrimination when embedding columns exist.
- Graceful metric skip with reason when required columns are missing.

Done when:

- One synthetic CSV can be evaluated end to end.

### V9-NOW-003: Run manifest for evaluations

Status: complete as of 2026-05-21 for local evaluation reports.

Goal:

Every evaluation output records task id, input file hashes, metrics, command,
timestamp, and package version placeholder.

Files likely touched:

- `spacebio_bench/reports/`
- `spacebio_bench/manifests.py`
- `spacebio_bench/schemas/`

Done when:

- Evaluator writes `metrics.json` and `run_manifest.json`.

## Active next task

Continue with
`V9-BULK-ALPHA-002: metadata-only alpha snapshot decision`.
Treat the OSD-120 metadata branch as closed unless explicit owner-supplied
release metadata appears. The next block should decide whether a metadata-only
public bulk alpha snapshot is acceptable with explicit payload-hash blockers,
or whether payload mirroring and local hash verification must precede any alpha
wording.
`V9-THEN-005: RO-Crate export design` remains open for the release-packaging
lane, but the active scientific lane now returns to multispecies expansion
because the organoid diagnostic family has a coherent draft-alpha boundary.

## Next: public bulk benchmark alpha

### V9-NEXT-001: Bulk data loader adapter

Status: complete as of 2026-05-21 for read-only legacy fold-file loading.

Goal:

Expose legacy processed task data through a v9 loader.

Requirements:

- Start read-only.
- Prefer existing processed task files and task metadata.
- No new data generation in first pass.
- Raise clear error when required processed files are missing.

Done when:

- `load_task("A2_gastrocnemius_bulk_lomo")` returns fold metadata and paths or
  tables needed for baseline evaluation.

### V9-NEXT-002: Nearest-centroid baseline

Status: complete as of 2026-05-21 for all eight generated bulk LOMO manifests.

Goal:

Add a transparent baseline that is easy to debug and useful for mission
discrimination.

Implemented files:

- `spacebio_bench/baselines/nearest_centroid.py`
- `scripts/run_v9_nearest_centroid.py`
- `v9/reports/nearest_centroid/`
- `v9/reports/README.md`

Done when:

- Runs on one task.
- Writes predictions and metrics.
- Has deterministic tests.

Actual output:

- Per-task `predictions.csv`, `metrics.json`, and `run_manifest.json`.
- Cross-task summary:
  `v9/reports/nearest_centroid/bulk_lomo_summary.csv`.
- All eight tasks completed with `status=evaluated`.

### V9-NEXT-003: Logistic/PCA-LR baseline adapter

Status: complete as of 2026-05-21 for all eight generated bulk LOMO manifests.

Goal:

Connect v9 runner to existing baseline patterns without rewriting all legacy
code.

Implemented files:

- `spacebio_bench/baselines/sklearn_classifiers.py`
- `scripts/run_v9_sklearn_baselines.py`
- `v9/reports/sklearn_baselines/`

Done when:

- At least one canonical bulk LOMO task can be evaluated from v9 CLI.
- Baseline output paths are manifest-linked.

Actual output:

- PCA-LR per-task `predictions.csv`, `metrics.json`, and `run_manifest.json`.
- L2 logistic-regression per-task `predictions.csv`, `metrics.json`, and
  `run_manifest.json`.
- Cross-sklearn summary:
  `v9/reports/sklearn_baselines/bulk_lomo_summary.csv`.
- All sixteen sklearn task/baseline rows completed with `status=evaluated`.

### V9-NEXT-004: All-bulk-task baseline report

Status: complete as of 2026-05-21 for nearest-centroid, PCA-LR, and L2 logistic
regression.

Goal:

Generate a v9-alpha bulk report across all 8 current manifests and all selected
simple baselines.

Output:

- `v9/reports/nearest_centroid/bulk_lomo_summary.csv`
- `v9/reports/sklearn_baselines/bulk_lomo_summary.csv`
- `v9/reports/bulk_lomo_baseline_summary.csv`
- `v9/reports/bulk_lomo_baseline_summary.json`
- `v9/reports/README.md`

Done when:

- Report is regenerated by one command and checked by tests.

## Then: source freeze and packaging

### V9-THEN-001: Source inventory

Status: complete as of 2026-05-21 for generated bulk LOMO manifests.

Goal:

Create a public source inventory table independent of generated task manifests.

Output:

- `v9/source_inventory.csv`
- `v9/source_inventory.json`

Implemented files:

- `spacebio_bench/sources.py`
- `scripts/build_v9_source_inventory.py`

Fields:

- source_id
- glds_prefix
- osd_url
- tissue
- mission
- task_ids
- access_status
- privacy_class
- checksum_status
- release_target
- notes

Actual output:

- 22 unique public source rows deduplicated from the eight generated bulk LOMO
  manifests.
- Source rows aggregate `task_ids`, `tissue`, `variant`, access/privacy status,
  OSDR URL, GLDS prefix, mission label, and current checksum status.
- Current checksum status remains `legacy_task_source_unfrozen` for all rows.

### V9-THEN-002: Checksum audit

Status: complete as of 2026-05-21 for OSDR API and checksum-manifest evidence.

Goal:

Promote `legacy_task_source_unfrozen` toward a v9 alpha source freeze.

Implemented files:

- `spacebio_bench/source_audit.py`
- `scripts/audit_v9_source_checksums.py`
- `v9/source_checksum_audit.csv`
- `v9/source_checksum_audit.json`

Actual output:

- 22 of 22 source rows returned `api_status=ok`.
- 22 of 22 source rows returned `audit_status=checksum_manifest_parsed`.
- 39 checksum manifest-like files were discovered through the OSDR file-list
  API.
- 8,439 MD5 checksum entries were parsed from OSDR checksum manifests.
- 8,275 parsed entries matched OSDR file-list payload names by exact, basename,
  or suffix matching.
- `freeze_ready=false` for all rows because payload files were not downloaded
  and hashed in this audit.

Done when:

- Each public bulk source has file-level checksum evidence or an explicit
  pending reason.

### V9-THEN-003: Public data package design

Status: complete as of 2026-05-21 for draft descriptor and package boundary.

Goal:

Design how v9 public bulk artifacts should be packaged outside Git before
dataset-card or RO-Crate release language is written.

Output:

- `docs/V9_PUBLIC_BULK_PACKAGE_DESIGN.md`
- `spacebio_bench/datapackage.py`
- `scripts/build_v9_datapackage_draft.py`
- `v9/datapackage.draft.json`

Actual output:

- Draft package boundary separates metadata spine, public bulk payload bundle,
  benchmark output bundle, and deferred/excluded artifacts.
- `v9/datapackage.draft.json` describes 11 resources.
- The descriptor includes 8 task manifests, 24 prediction CSVs, 24 metrics JSON
  files, and 24 run-manifest JSON files as file collections.
- Large public bulk fold payload files remain indexed by `v9/task_data_index.csv`
  but are not payload-hashed yet.

Done when:

- Git-tracked artifacts, downloadable data bundle artifacts, checksum audit
  outputs, and future `datapackage.json` resources are separated clearly.
- The design explains that checksum-manifest evidence is present, while payload
  verification remains a separate freeze step.

### V9-THEN-004: Dataset card draft

Status: complete as of 2026-05-21 for draft release-facing language.

Goal:

Prepare a v9 Hugging Face dataset card before upload.

Output:

- `docs/v9_hf_dataset_card.md`

Actual output:

- Draft YAML metadata for a Hugging Face-style dataset card.
- Dataset summary, intended uses, out-of-scope uses, structure, task table,
  source table, provenance/integrity status, baseline summary, privacy/access
  boundaries, limitations, licensing/citation notes, and maintenance checklist.
- Explicit draft language:
  `draft_not_frozen` and
  `checksum_manifests_parsed_payloads_not_hashed`.
- OSDR citation/access guidance included from NASA OSDR FAQ and terms pages.

Done when:

- Includes task table, source policy, citation, limitations, privacy notes,
  and artifact split.

### V9-THEN-005: RO-Crate export design

Goal:

Define metadata export before implementation.

Output:

- `docs/V9_RO_CRATE_EXPORT_DESIGN.md`

Done when:

- Lists entities, files, authorship, source datasets, workflows, and checksums.

### V9-THEN-006: Organism/material schema design

Status: complete as of 2026-05-21 for source/task metadata scaffolding and
draft extension inventories.

Goal:

Extend v9 source/task schema planning for non-mouse and non-tissue systems
before importing organoid or multispecies data.

Output:

- Schema/code changes covering organism, taxon id, material type, model system,
  assay modality, feature namespace, and orthology strategy.
- Draft human organoid and multispecies source inventories.

Actual output:

- Existing mouse bulk manifests still validate.
- `v9/source_inventory.csv` now carries organism/material/model-system fields
  for the current 22 mouse bulk source rows.
- `v9/human_organoid/source_inventory.draft.csv` contains OSD-863 and OSD-871.
- `v9/human_organoid/source_checksum_audit.draft.csv` contains API-ok and
  checksum-manifest-parsed rows for OSD-863 and OSD-871.
- `v9/human_organoid/sample_table_audit.draft.csv` contains parsed sample-table
  evidence: 19 rows for OSD-863 and 23 rows for OSD-871.
- `v9/human_organoid/sample_factors.draft.csv` parses all 42 sample conditions
  into disease context, Ground/LEO label, microglia condition, and organoid type.
- `v9/multispecies/source_inventory.draft.csv` contains OSD-207, OSD-37, and
  OSD-120.
- Draft source rows for human organoids and non-mouse species no longer need to
  overload the current `tissue` field.

Done when:

- Source/task schema metadata exists and v9 tests pass.

## Flagship branch A: single-cell track

### V9-SC-001: RRRM asset inventory

Goal:

Find every RRRM-1/RRRM-2 scRNA-seq input, script, result, and doc in the repo.

Output:

- `v9/sc_spaceflight/asset_inventory.md`

### V9-SC-002: AnnData task manifest draft

Goal:

Draft one `sc_spaceflight` manifest without requiring full evaluator support.

Output:

- `v9/sc_spaceflight/task_manifests/*.json`

### V9-SC-003: `genelab_sc` metric specification

Goal:

Write formulas and required inputs before code.

Output:

- `docs/V9_SC_METRIC_SPEC.md`

## Flagship branch B: radiation/stressor track

### V9-STRESS-001: DECOMPOSE asset inventory

Goal:

Map v8 DECOMPOSE outputs into candidate v9 stressor tasks.

Output:

- `v9/stressor_radiation_quality/asset_inventory.md`

### V9-STRESS-002: Stressor task manifest draft

Goal:

Draft one `stressor_radiation_quality` manifest.

Output:

- `v9/stressor_radiation_quality/task_manifests/*.json`

### V9-STRESS-003: `stressor_regime` metric specification

Goal:

Define sign consistency, saturation sensitivity, uncertainty, and held-out
analog metrics.

Output:

- `docs/V9_STRESSOR_METRIC_SPEC.md`

## Flagship branch C: human organoid and multispecies extension

### V9-ORG-001: Human organoid schema extension

Status: initial scaffold complete as of 2026-05-21; task manifest support still
needs a concrete organoid task card.

Goal:

Extend manifest/source schemas for public human organoid spaceflight studies.

Output:

- Schema fields for organism, material type, model system, organoid type,
  microglia condition, donor block, disease context, and assay modality.

### V9-ORG-002: OSD-863/OSD-871 source inventory

Status: draft inventory complete as of 2026-05-21.

Goal:

Build a draft source inventory for public human neural organoid RNA-seq.

Output:

- `v9/human_organoid/source_inventory.draft.csv`
- `v9/human_organoid/source_inventory.draft.json`

### V9-ORG-003: Human organoid task-card draft

Status: draft manifest complete as of 2026-05-21.

Goal:

Draft one `human_organoid_spaceflight` task manifest for GSE259421/OSD-863/
OSD-871 without claiming a frozen benchmark result.

Output:

- `v9/human_organoid/task_manifests/*.json`
- `v9/human_organoid/task_manifest_index.draft.csv`
- `v9/human_organoid/task_manifest_index.draft.json`

Actual output:

- `draft_human_organoid_spaceflight` manifest with OSD-863 and OSD-871 source
  records.
- Draft blocked LEO/ISS-versus-ground split design.
- Explicit donor, organoid type, and microglia condition blocking factors.
- Explicit single-mission limitation and pending payload-audit status.
- Updated after sample-table audit to record 19 OSD-863 rows and 23 OSD-871
  rows, with manual condition-factor mapping still pending.
- Updated after sample-factor parsing to record 42 parsed samples and four
  sample-count-backed draft folds.

### V9-ORG-004: OSD-863/OSD-871 checksum-manifest audit

Status: complete as of 2026-05-21 for API and checksum-manifest evidence.

Goal:

Run the existing OSDR file-list/checksum audit over the human organoid draft
sources.

Output:

- `v9/human_organoid/source_checksum_audit.draft.csv`
- `v9/human_organoid/source_checksum_audit.draft.json`

Actual output:

- OSD-863: API ok, 338 listed files, 2 checksum manifests, 335 parsed MD5
  entries, 296 payload-name matches.
- OSD-871: API ok, 402 listed files, 2 checksum manifests, 399 parsed MD5
  entries, 352 payload-name matches.
- `freeze_ready=false` remains correct because payload files were not
  downloaded and hashed.

### V9-ORG-005: OSD-863/OSD-871 sample-table audit

Status: complete as of 2026-05-21 for sample-table discovery and row counts.

Goal:

Find and parse OSDR SampleTable files for the human organoid draft sources.

Output:

- `v9/human_organoid/sample_table_audit.draft.csv`
- `v9/human_organoid/sample_table_audit.draft.json`

Actual output:

- OSD-863: API ok, one SampleTable file parsed, 19 rows.
- OSD-871: API ok, one SampleTable file parsed, 23 rows.
- Combined sample-table rows match the 42 public GSE259421 samples.
- OSDR encodes disease context, ground/space condition, and microglia state in
  a compact `condition` field; this field is parsed in `V9-ORG-006`.

### V9-ORG-006: Human organoid condition-factor parser

Status: complete as of 2026-05-21 for OSD-863/OSD-871 draft factors.

Goal:

Parse the compact OSDR `condition` field into explicit sample-level factors.

Output:

- `v9/human_organoid/sample_factors.draft.csv`
- `v9/human_organoid/sample_factors.draft.json`

Actual output:

- 42 parsed sample rows.
- Labels: 22 Ground, 20 LEO_or_ISS.
- Organoid types: 19 cortical, 23 dopaminergic.
- Microglia conditions: 20 with_microglia, 22 without_microglia.
- Disease contexts: 21 no_known_diseases, 8 primary_progressive_multiple_sclerosis,
  13 sporadic_parkinson_disease.
- Draft manifest now includes four sample-count-backed folds:
  holdout cortical, holdout dopaminergic, holdout with_microglia, and holdout
  without_microglia.

### V9-ORG-007: Human organoid expression-matrix audit

Status: complete as of 2026-05-21 for canonical normalized-count matrix
download and sample alignment.

Goal:

Find the processed OSDR normalized-count matrices for OSD-863 and OSD-871,
download only the canonical non-rRNArm normalized matrices, hash them, and
verify their sample columns against `sample_factors.draft.csv`.

Output:

- `v9/human_organoid/expression_matrix_audit.draft.csv`
- `v9/human_organoid/expression_matrix_audit.draft.json`
- `v9/human_organoid/matrices/GLDS-716_rna_seq_Normalized_Counts_GLbulkRNAseq.csv`
- `v9/human_organoid/matrices/GLDS-720_rna_seq_Normalized_Counts_GLbulkRNAseq.csv`

Actual output:

- OSD-863: one canonical normalized matrix downloaded, SHA-256 recorded,
  30,408 feature rows, 19 sample columns, 19/19 parsed sample-factor rows
  matched.
- OSD-871: one canonical normalized matrix downloaded, SHA-256 recorded,
  30,269 feature rows, 23 sample columns, 23/23 parsed sample-factor rows
  matched.
- Draft task manifest now records `expression_matrix_status` as
  `matrix_downloaded_sample_aligned` in both `split` and `provenance`.

### V9-ORG-008: Human organoid baseline loader and first baseline

Status: complete as of 2026-05-21 for a draft nearest-centroid pilot baseline.

Goal:

Add a small, deterministic loader over the aligned human organoid matrices and
run the first no-claim pilot baseline over the sample-count-backed folds.

Likely files:

- `spacebio_bench/data/organoid.py`
- `spacebio_bench/baselines/organoid.py`
- `scripts/run_v9_human_organoid_baseline.py`
- `v9/human_organoid/reports/`
- `tests/test_v9_spacebio_bench.py`

Done when:

- The loader returns feature rows, sample IDs, labels, and fold definitions from
  the draft manifest plus matrix audit.
- One transparent baseline writes predictions, metrics, and run-manifest files
  under `v9/human_organoid/reports/`.
- The output language remains `draft_not_frozen` and does not claim a stable
  leaderboard.

Actual output:

- `load_human_organoid_task()` returns 42 samples, 27,986 common human gene
  features across OSD-863/OSD-871, two source matrix paths, and four
  sample-count-backed draft folds.
- `scripts/run_v9_human_organoid_baseline.py` writes:
  `v9/human_organoid/reports/nearest_centroid/draft_human_organoid_spaceflight/predictions.csv`,
  `metrics.json`, `run_manifest.json`, and
  `human_organoid_baseline_summary.csv/json`.
- Current draft-only summary row:
  `n_predictions=84`, `balanced_accuracy=0.5295454545`,
  `auroc=0.6147727273`, `calibration_error=0.02537185971`.
- `mission_discrimination` is skipped because the pilot predictions do not emit
  embeddings and the data are single-mission.

### V9-ORG-009: Human organoid baseline robustness review

Status: complete as of 2026-05-21 for first-pass sensitivity and confounding review.

Goal:

Stress-test the draft organoid baseline before using it in any manuscript or
release-facing language.

Likely work:

- Run sensitivity variants for `top_variable_genes`, transform, and scaling.
- Add a small baseline summary comparison table under
  `v9/human_organoid/reports/`.
- Review fold leakage and confounding risk across organoid type, microglia
  condition, disease context, and source accession.
- Decide whether the next organoid metric should be signature correlation, DE
  direction agreement, or only classifier metrics until DE references are
  audited.

Done when:

- A short robustness note names which baseline configuration remains the default
  and which claims are still prohibited.

Actual output:

- `scripts/run_v9_human_organoid_sensitivity.py` runs a 20-variant grid:
  `transform in {log1p, none}`, `scaling in {zscore, none}`, and
  `top_variable_genes in {100, 500, 2000, 5000, 27986}`.
- `v9/human_organoid/reports/sensitivity/human_organoid_baseline_summary.csv`
  records all variant metrics and per-variant report paths.
- `v9/human_organoid/reports/ORGANOID_BASELINE_ROBUSTNESS_REVIEW.md` keeps
  `log1p/zscore/top_variable_genes=2000` as the conservative default.
- Sensitivity range across 20 variants:
  balanced accuracy 0.5159090909 to 0.7693181818, AUROC 0.5676136364 to
  0.9295454545, calibration error 0.0242551690 to 0.2445559995.
- The top balanced-accuracy variant is `tvg2000_none_zscore`, but it is not
  promoted to default because raw-scale variants show worse calibration and need
  library-size and confounding audits.

### V9-ORG-010: Human organoid donor and library-size audit

Status: complete as of 2026-05-21 for local matrix-scale diagnostics; updated
after V9-ORG-011 with GEO-derived donor/iPSC-line metadata.

Goal:

Resolve the biggest remaining interpretability blockers for the organoid pilot:
donor/iPSC-line structure and expression-scale diagnostics.

Likely work:

- Extract donor or iPSC-line identifiers from GEO/OSDR metadata or publication
  tables.
- Add sample-total/library-size diagnostics for both normalized matrices.
- Add a compact fold-level metric table for default and sensitivity-best
  variants.
- Decide whether donor-blocked folds are possible or whether donor remains a
  documented limitation.

Done when:

- The robustness review can say whether high raw-count sensitivity scores are
  likely biological signal, scale artifact, or unresolved.

Actual output:

- `scripts/audit_v9_human_organoid_diagnostics.py` generated per-sample scale,
  grouped scale, and donor metadata availability diagnostics.
- Current OSDR SampleTable files expose only `condition`; donor/iPSC-line fields
  are not available in the local OSDR sample table audit.
- GEO GSE259421 series metadata now recovers Subject and iPSC-line identifiers
  for 42/42 samples, merged into `v9/human_organoid/sample_factors.draft.csv`.
- Median sample sum differs by label: Ground is about 32.0M and LEO_or_ISS is
  about 28.1M across the 27,986 common features.
- Median zero fraction also differs by label: Ground is about 0.1074 and
  LEO_or_ISS is about 0.08345.
- `v9/human_organoid/reports/ORGANOID_DONOR_LIBRARY_SIZE_AUDIT.md` keeps the
  raw-count sensitivity variants unresolved and preserves the log1p/z-score
  default.

### V9-ORG-011: Human organoid external donor metadata recovery

Status: complete as of 2026-05-21.

Goal:

Search GEO series/sample metadata, supplementary tables, and the publication for
donor or iPSC-line identifiers that are absent from the current OSDR SampleTable.

Likely work:

- Fetch or inspect GEO series metadata for GSE259421.
- Map GSM samples to donor/iPSC-line if fields exist.
- Extend `sample_factors.draft.csv` with `donor_or_line_id` only if evidence is
  explicit.
- If donor metadata cannot be recovered, update the task manifest to mark donor
  blocking as unresolved rather than merely required.

Done when:

- Donor blocking is either implemented with source-backed sample IDs or closed
  as unavailable for this draft task.

Actual output:

- `spacebio_bench/organoid_geo.py` parses GEO GSE259421 series matrix metadata.
- `scripts/build_v9_human_organoid_geo_metadata.py` writes
  `v9/human_organoid/geo_sample_metadata.draft.csv/json` and merges donor fields
  into `v9/human_organoid/sample_factors.draft.csv/json`.
- Recovered donor/line distribution:
  Subject1/`051121-01-MR-017` n=11, Subject2/`AK003-01-MR-008` n=8,
  Subject3/`UEC741iPS517` n=10, Subject4/`HDF410iPS504` n=13.
- The task manifest records
  `sample_factor_status=condition_factors_and_geo_donor_line_metadata_parsed`
  and `donor_metadata_status=parsed_from_geo_series_matrix`.
- Default benchmark folds remain organoid-type and microglia-condition folds.
  Donor holdouts are recorded as `donor_diagnostic_folds` only because donor,
  source, organoid fate, and disease context are not independently crossed.

### V9-ORG-012: Human organoid donor-aware split decision

Status: complete as of 2026-05-22.

Goal:

Decide whether GEO-derived donor holdouts should remain diagnostic-only or
become a separate conservatively worded draft task family.

Likely work:

- Compare default four folds against donor diagnostic folds without changing the
  default baseline contract.
- Write a short split-design note explaining the confounding structure.
- Decide whether a new `human_organoid_donor_diagnostic` task id is useful or
  whether donor-aware metrics should remain manifest metadata only.

Done when:

- Donor-aware split language is explicit enough that no reader can mistake the
  current organoid pilot for a clean donor-generalization benchmark.

Actual output:

- `HumanOrganoidTaskData` now exposes `diagnostic_folds` from manifest
  `donor_diagnostic_folds`.
- `run_human_organoid_donor_diagnostics` and
  `scripts/run_v9_human_organoid_donor_diagnostics.py` evaluate donor holdouts
  without changing the default fold family.
- Diagnostic outputs:
  `v9/human_organoid/reports/donor_diagnostics/human_organoid_baseline_summary.csv`,
  per-task `predictions.csv`, `metrics.json`, and `run_manifest.json`.
- Donor diagnostic default metrics:
  balanced accuracy 0.5318181818, AUROC 0.6954545455, calibration error
  0.02808852089, macro-F1 0.5138888889.
- `v9/human_organoid/reports/ORGANOID_DONOR_AWARE_SPLIT_DECISION.md` records
  the decision to keep donor holdouts diagnostic-only and not create a separate
  leaderboard task yet.

### V9-ORG-013: Human organoid DE/signature metric reference design

Status: complete as of 2026-05-21.

Goal:

Move beyond classifier-only organoid metrics by defining the first DE/signature
reference metric that matches the biological claims in GSE259421 and Marotta et
al. 2024.

Output:

- `spacebio_bench/organoid_signature_audit.py`
- `scripts/audit_v9_human_organoid_signature_references.py`
- `v9/human_organoid/signature_reference_audit.draft.csv`
- `v9/human_organoid/signature_reference_audit.draft.json`
- `v9/human_organoid/reports/ORGANOID_SIGNATURE_METRIC_REFERENCE_DECISION.md`
- `v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json`

Decision:

- OSDR lists public differential-expression tables and contrast-definition CSVs
  for both OSD-863/GLDS-716 and OSD-871/GLDS-720.
- Each source has 56 parsed contrast pairs across 8 groups, including four
  direct matched Ground Control versus Space Flight contrasts and four reversed
  matches.
- `de_direction_match` and `signature_rank_correlation` remain declared and are
  now reference-backed, but non-primary until the response-signature scorer is
  promoted.

### V9-ORG-014: Human organoid frozen DE contrast extraction and signature metric contract

Status: complete as of 2026-05-21 for draft input contract and skip-aware
evaluator plumbing.

Goal:

Turn the audited public OSDR DE references into a compact, checksummed reference
artifact and an evaluator contract for gene-level response signatures.

Implemented work:

- Download the canonical non-rRNArm OSDR DE tables for OSD-863 and
  OSD-871.
- Extract the four direct matched Ground Control versus Space Flight contrasts
  per source.
- Normalize log2 fold-change orientation to `LEO_or_ISS - Ground`.
- Write a compressed derived reference table and JSON manifest.
- Define a `response_signature.csv` submission artifact and skip-aware
  evaluator behavior for submissions that do not provide gene-level signatures.

Done when:

- A derived DE reference table exists with source, contrast, gene id,
  log2 fold-change, adjusted p-value, and orientation metadata.
- The evaluator can report DE/signature metrics as skipped with a precise reason
  unless a valid response-signature artifact is supplied.

Actual output:

- `spacebio_bench/organoid_de_reference.py`
- `scripts/build_v9_human_organoid_de_reference.py`
- `v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz`
- `v9/human_organoid/de_references/human_organoid_de_reference_manifest.draft.json`
- `v9/human_organoid/reports/ORGANOID_DE_REFERENCE_CONTRACT.md`
- Evaluator skip entries for `de_direction_match` and
  `signature_rank_correlation` when `response_signature.csv` is missing.

### V9-ORG-015: Human organoid response-signature scorer

Status: complete as of 2026-05-22 for diagnostic scorer plumbing.

Goal:

Compute diagnostic DE/signature metrics when a valid `response_signature.csv`
artifact is supplied.

Implemented work:

- Add `response_signature.csv` validation.
- Load `.csv.gz` reference tables.
- Join by `source_id`, `contrast_id`, and `gene_symbol`, with Ensembl fallback.
- Compute `de_direction_match` on significant reference genes.
- Compute `signature_rank_correlation` over joined all-gene signatures.
- Report per-contrast details and aggregate status in `metrics.json`.

Done when:

- A synthetic response-signature artifact produces computed
  `de_direction_match` and `signature_rank_correlation` values.
- Missing, invalid, or underpowered response signatures still skip with clear
  reasons.

Actual output:

- `spacebio_bench/signature_metrics.py`
- `spacebio_bench/evaluate.py`
- `scripts/evaluate_v9_submission.py`
- `v9/human_organoid/reports/ORGANOID_RESPONSE_SIGNATURE_SCORER.md`
- Tests covering computed, missing, and invalid response-signature artifacts.

### V9-ORG-016: Human organoid response-signature scorer dry run

Status: complete as of 2026-05-22 for smoke-test report generation.

Goal:

Create a reproducible end-to-end dry run for the response-signature scorer
against the real OSD-863/OSD-871 derived DE reference.

Implemented work:

- Generate a small example `response_signature.csv` from selected real
  reference rows, clearly marked as a scorer smoke test rather than a model
  result.
- Evaluate a synthetic prediction CSV plus example response-signature artifact.
- Write a report under `v9/human_organoid/reports/response_signature_smoke/`.
- Keep all language diagnostic-only and non-leaderboard.

Done when:

- The scorer runs against the real `.csv.gz` reference table.
- Output `metrics.json` contains computed DE/signature diagnostics and
  per-contrast details.
- The report states that the response signature is an example artifact, not a
  biological model claim.

Actual output:

- `scripts/run_v9_human_organoid_response_signature_smoke.py`
- `v9/human_organoid/reports/response_signature_smoke/predictions.csv`
- `v9/human_organoid/reports/response_signature_smoke/response_signature.csv`
- `v9/human_organoid/reports/response_signature_smoke/metrics.json`
- `v9/human_organoid/reports/response_signature_smoke/run_manifest.json`
- `v9/human_organoid/reports/response_signature_smoke/README.md`

### V9-ORG-017: Human organoid model-produced response-signature adapter design

Status: complete as of 2026-05-22.

Goal:

Design how a real baseline or model adapter can emit `response_signature.csv`
without using the derived DE reference as an input.

Implemented work:

- Review whether a fold-level training-only expression-difference baseline can
  produce per-contrast predicted signatures without reference leakage.
- Define how source/disease/microglia contrast IDs map to train/test sample
  factors.
- Decide whether the first real adapter should be train-fold-only empirical
  signatures, a simple classifier-derived feature effect, or a separate model
  output contract only.
- Keep outputs diagnostic-only until confounding and fold leakage risks are
  reviewed.

Done when:

- A design note names the first non-leaky adapter path and its risks.
- The next implementation block has files, tests, and artifact boundaries
  specified.

Actual decision:

- First adapter:
  `organoid_source_transfer_empirical_signature`.
- Predict OSD-863 target contrasts from OSD-871 training samples.
- Predict OSD-871 target contrasts from OSD-863 training samples.
- Do not use target-source expression or DE references for signature
  generation.
- Keep output diagnostic-only.

Actual output:

- `v9/human_organoid/reports/ORGANOID_RESPONSE_SIGNATURE_ADAPTER_DESIGN.md`

### V9-ORG-018: Human organoid source-transfer response-signature adapter

Status: complete as of 2026-05-22.

Goal:

Implement the source-transfer empirical response-signature adapter and generate
a diagnostic-only report.

Implemented work:

- Add `spacebio_bench/response_signature_adapters.py`.
- Add `scripts/run_v9_human_organoid_source_transfer_signature.py`.
- Use organoid-type/source holdout folds to generate one train-source global
  LEO/ISS-minus-Ground signature per target source.
- Emit `response_signature.csv` with required scorer columns plus provenance
  columns such as `fold_id`, `training_source_id`, `target_source_id`,
  `n_train_ground`, `n_train_leo_or_iss`, and
  `reference_usage_policy=reference_not_used_for_signature_generation`.
- Evaluate against the derived DE reference and write outputs under
  `v9/human_organoid/reports/source_transfer_signature/`.

Done when:

- Target-source samples are excluded from signature generation.
- The evaluator computes diagnostic DE/signature metrics from the generated
  artifact.
- The report clearly states that the baseline is source-transfer, weak, and
  non-leaderboard.

Actual output:

- `spacebio_bench/response_signature_adapters.py`
- `scripts/run_v9_human_organoid_source_transfer_signature.py`
- `v9/human_organoid/reports/source_transfer_signature/predictions.csv`
- `v9/human_organoid/reports/source_transfer_signature/response_signature.csv.gz`
- `v9/human_organoid/reports/source_transfer_signature/response_signature_metadata.json`
- `v9/human_organoid/reports/source_transfer_signature/metrics.json`
- `v9/human_organoid/reports/source_transfer_signature/run_manifest.json`
- `v9/human_organoid/reports/source_transfer_signature/README.md`

Observed diagnostic result:

- response rows: 223,888.
- joined rows: 223,888.
- significant reference rows scored for direction: 2,346.
- direction matches: 1,808.
- `de_direction_match`: 0.7706734867860188.
- `signature_rank_correlation`: 0.1760078660242601.

### V9-ORG-019: Human organoid source-transfer diagnostic review

Status: complete as of 2026-05-22.

Goal:

Interpret the first model-produced response-signature diagnostic report and
decide whether it should become the default organoid signature baseline.

Likely work:

- Inspect per-contrast `de_direction_match` and `signature_rank_correlation`
  details in `metrics.json`.
- Compare direction-match behavior between OSD-863 and OSD-871 targets.
- Review whether global source-transfer signatures are too crude for
  disease/microglia-specific contrasts.
- Decide if the next implementation should add per-condition training
  signatures, classifier-coefficient signatures, or keep the current baseline
  as the conservative first diagnostic.

Done when:

- A review note records what the source-transfer diagnostic does and does not
  mean.
- The next adapter/scoring improvement is selected with explicit leakage
  boundaries.

Actual output:

- `v9/human_organoid/reports/ORGANOID_SOURCE_TRANSFER_SIGNATURE_REVIEW.md`

Decision:

- Keep `organoid_source_transfer_empirical_signature` as the first conservative
  response-signature diagnostic baseline.
- Do not promote it to a primary metric or leaderboard result.
- Design a per-condition source-transfer adapter next, while preserving
  train-source-only and reference-not-used leakage boundaries.

### V9-ORG-020: Human organoid per-condition source-transfer adapter design

Status: complete as of 2026-05-22.

Goal:

Design a second response-signature adapter that compares against the current
global source-transfer signature by matching disease context, microglia
condition, or both on the training-source side when feasible.

Likely work:

- Inventory which target contrasts have matching training-source strata.
- Decide skip versus labeled fallback behavior when a target stratum is absent
  from the training source.
- Preserve the rule that target-source expression and DE-reference tables are
  not used for signature generation.
- Define output metadata fields so global and per-condition signatures can be
  compared without ambiguity.

Done when:

- A design note records per-condition adapter rules, skip/fallback behavior,
  leakage boundaries, files, and tests for implementation.

Actual output:

- `v9/human_organoid/reports/ORGANOID_PER_CONDITION_SIGNATURE_ADAPTER_DESIGN.md`

Decision:

- Implement a microglia-matched source-transfer empirical signature next.
- Defer disease-matched and disease+microglia matched signatures because disease
  context is not fully crossed between OSD-863 and OSD-871.
- Keep all signatures train-source-only and reference-not-used during
  generation.

### V9-ORG-021: Human organoid microglia-matched source-transfer adapter

Status: complete as of 2026-05-22.

Goal:

Implement the microglia-matched source-transfer empirical response-signature
adapter and compare it with the current global source-transfer baseline.

Likely work:

- Extend `spacebio_bench/response_signature_adapters.py` with a condition-
  filtered source-transfer helper.
- Add a runner under
  `scripts/run_v9_human_organoid_microglia_source_transfer_signature.py`.
- Generate a report under
  `v9/human_organoid/reports/microglia_source_transfer_signature/`.
- Evaluate `de_direction_match` and `signature_rank_correlation`.
- Add a README comparison against
  `v9/human_organoid/reports/source_transfer_signature/metrics.json`.

Done when:

- The report writes a valid `response_signature.csv.gz`.
- All emitted target contrasts are conditioned on the target contrast's
  microglia state using training-source samples only.
- The report records diagnostic metrics and deltas versus the global
  source-transfer baseline.
- Targeted unit tests pass.

Actual output:

- `spacebio_bench/response_signature_adapters.py`
- `scripts/run_v9_human_organoid_microglia_source_transfer_signature.py`
- `v9/human_organoid/reports/microglia_source_transfer_signature/predictions.csv`
- `v9/human_organoid/reports/microglia_source_transfer_signature/response_signature.csv.gz`
- `v9/human_organoid/reports/microglia_source_transfer_signature/response_signature_metadata.json`
- `v9/human_organoid/reports/microglia_source_transfer_signature/metrics.json`
- `v9/human_organoid/reports/microglia_source_transfer_signature/run_manifest.json`
- `v9/human_organoid/reports/microglia_source_transfer_signature/README.md`

Observed diagnostic result:

- response rows: 223,888.
- joined rows: 223,888.
- significant reference rows scored for direction: 2,345.
- direction matches: 1,850.
- `de_direction_match`: 0.7889125799573561.
- `signature_rank_correlation`: 0.1500722829461316.
- direction delta versus global source-transfer: +0.018239093171337317.
- rank-correlation delta versus global source-transfer:
  -0.025935583078128516.

### V9-ORG-022: Human organoid source-transfer adapter comparison review

Status: complete as of 2026-05-22.

Goal:

Compare the global and microglia-matched source-transfer response-signature
diagnostics and decide which adapter family should come next.

Likely work:

- Review per-contrast deltas between global and microglia-matched reports.
- Decide whether direction improvement is enough to keep microglia matching as
  a standard secondary diagnostic despite weaker rank correlation.
- Decide whether to add a partial control-only disease-matched diagnostic or
  move to classifier-coefficient signatures.
- Keep all claims diagnostic-only and non-leaderboard.

Done when:

- A short comparison note records the adapter interpretation and next adapter
  choice.

Actual output:

- `v9/human_organoid/reports/ORGANOID_SOURCE_TRANSFER_ADAPTER_COMPARISON_REVIEW.md`

Decision:

- Keep the global source-transfer adapter as the first conservative diagnostic.
- Keep the microglia-matched adapter as a secondary condition-sensitivity
  diagnostic.
- Design a partial shared-control disease+microglia source-transfer diagnostic
  next.

### V9-ORG-023: Human organoid shared-control disease+microglia adapter design

Status: complete as of 2026-05-22.

Goal:

Design a partial-coverage source-transfer diagnostic for the shared
`no_known_diseases` strata, matching both disease context and microglia
condition when the opposite-source training stratum exists.

Likely work:

- Identify eligible target contrasts where disease and microglia are both
  available in the opposite-source training fold.
- Define skipped-contrast metadata for PPMS and sporadic Parkinson disease
  target contrasts.
- Keep the adapter partial-coverage and non-leaderboard.
- Decide whether to implement the partial report immediately or leave it as a
  design-only diagnostic path.

Done when:

- A design note records eligible contrasts, skipped contrasts, leakage
  boundaries, output schema, and tests for an optional implementation.

Actual output:

- `v9/human_organoid/reports/ORGANOID_SHARED_CONTROL_SIGNATURE_ADAPTER_DESIGN.md`

Decision:

- Emit only the four target contrasts with shared `no_known_diseases` plus
  matched microglia training strata.
- Explicitly skip the four disease-specific target contrasts.
- Keep aggregate metrics partial-coverage and non-leaderboard.

### V9-ORG-024: Human organoid shared-control source-transfer adapter

Status: complete as of 2026-05-22.

Goal:

Implement the partial shared-control disease+microglia source-transfer
response-signature adapter and generate a partial-coverage diagnostic report.

Likely work:

- Extend `spacebio_bench/response_signature_adapters.py` with disease+microglia
  condition filtering and skipped-contrast metadata.
- Add a runner under
  `scripts/run_v9_human_organoid_shared_control_source_transfer_signature.py`.
- Generate a report under
  `v9/human_organoid/reports/shared_control_source_transfer_signature/`.
- Evaluate emitted rows against the derived DE reference.
- Include README coverage and per-contrast comparison tables.

Done when:

- The report writes a valid partial `response_signature.csv.gz`.
- Four shared-control contrasts are emitted and four disease-specific contrasts
  are skipped in metadata.
- Targeted tests pass.

Actual output:

- `spacebio_bench/response_signature_adapters.py`
- `scripts/run_v9_human_organoid_shared_control_source_transfer_signature.py`
- `v9/human_organoid/reports/shared_control_source_transfer_signature/predictions.csv`
- `v9/human_organoid/reports/shared_control_source_transfer_signature/response_signature.csv.gz`
- `v9/human_organoid/reports/shared_control_source_transfer_signature/response_signature_metadata.json`
- `v9/human_organoid/reports/shared_control_source_transfer_signature/metrics.json`
- `v9/human_organoid/reports/shared_control_source_transfer_signature/run_manifest.json`
- `v9/human_organoid/reports/shared_control_source_transfer_signature/README.md`

Observed diagnostic result:

- emitted target contrasts: 4.
- skipped target contrasts: 4.
- response rows: 111,944.
- joined rows: 111,944.
- significant reference rows scored for direction: 517.
- direction matches: 305.
- `de_direction_match`: 0.5899419729206963.
- `signature_rank_correlation`: 0.03566701395309356.

Decision:

- Do not promote disease+microglia stratification further as a stronger
  empirical diagnostic.
- Treat this as evidence that the small shared-control strata are too noisy for
  better all-gene response-signature agreement.
- Move next to classifier-derived feature-effect contract design.

### V9-ORG-025: Human organoid classifier-coefficient signature contract design

Status: complete as of 2026-05-23.

Goal:

Design how classifier-derived feature effects can be reported and scored
without labeling them as log2FC response signatures.

Likely work:

- Review whether PCA-LR or L2 logistic coefficients can be mapped to a separate
  feature-effect artifact.
- Define whether coefficient diagnostics should be rank-only, top-k overlap, or
  sign-only against DE references.
- Keep the current `response_signature.csv` log2FC contract separate from any
  coefficient artifact.
- Decide files, tests, and report boundary for an optional implementation.

Done when:

- A contract note names the coefficient artifact schema, metric policy, leakage
  boundary, and whether implementation should proceed.

Actual output:

- `v9/human_organoid/reports/ORGANOID_CLASSIFIER_FEATURE_EFFECT_CONTRACT.md`

Decision:

- Do not write classifier coefficients to `response_signature.csv`.
- Add a separate optional `feature_effect.csv` contract.
- Keep rank/sign/top-k feature-effect diagnostics non-primary.
- Implement L2 logistic gene-space coefficients before PCA-LR reconstructed
  coefficients.

### V9-ORG-026: Human organoid L2 logistic feature-effect pilot

Status: complete as of 2026-05-23.

Goal:

Implement a source-transfer L2 logistic feature-effect pilot that writes
`feature_effect.csv` and scores rank/sign/top-k diagnostics against the derived
DE reference without treating coefficients as log2FC signatures.

Likely work:

- Add `spacebio_bench/feature_effects.py`.
- Add `scripts/run_v9_human_organoid_logistic_feature_effect.py`.
- Generate `v9/human_organoid/reports/logistic_feature_effect/`.
- Add tests for artifact validation, no target-source leakage, and diagnostic
  metric computation.

Done when:

- The pilot writes a valid `feature_effect.csv`.
- The report states coefficients are not log2FC response signatures.
- Targeted tests pass.

Actual output:

- `spacebio_bench/feature_effects.py`
- `scripts/run_v9_human_organoid_logistic_feature_effect.py`
- `v9/human_organoid/reports/logistic_feature_effect/predictions.csv`
- `v9/human_organoid/reports/logistic_feature_effect/feature_effect.csv.gz`
- `v9/human_organoid/reports/logistic_feature_effect/feature_effect_metadata.json`
- `v9/human_organoid/reports/logistic_feature_effect/metrics.json`
- `v9/human_organoid/reports/logistic_feature_effect/run_manifest.json`
- `v9/human_organoid/reports/logistic_feature_effect/README.md`

Observed diagnostic result:

- feature-effect rows: 16,000.
- joined rows: 16,000.
- significant reference rows scored for direction: 204.
- direction matches: 124.
- `feature_effect_direction_match`: 0.6078431372549019.
- `feature_effect_rank_correlation`: 0.08672800238082004.
- top-k DE overlap: 1/50, 5/100, 10/250, 14/500.

Decision:

- Keep the feature-effect artifact separate from `response_signature.csv`.
- Treat L2 logistic coefficients as a weak/moderate secondary diagnostic, not a
  response-signature replacement.
- Review whether null/top-k calibration is needed before PCA-LR reconstructed
  weights.

### V9-ORG-027: Human organoid feature-effect diagnostic review

Status: complete as of 2026-05-23.

Goal:

Interpret the L2 logistic feature-effect pilot and decide whether coefficient
diagnostics should become a standard secondary report.

Likely work:

- Compare feature-effect direction/rank/top-k diagnostics with source-transfer
  response-signature diagnostics.
- Decide whether top-k overlap needs null calibration.
- Decide whether PCA-LR reconstructed gene weights should be attempted or
  deferred.

Done when:

- A review note records the feature-effect interpretation and next adapter or
  calibration choice.

Actual output:

- `v9/human_organoid/reports/ORGANOID_FEATURE_EFFECT_DIAGNOSTIC_REVIEW.md`

Decision:

- Keep `feature_effect.csv` as an optional secondary diagnostic.
- Do not promote L2 logistic coefficients to response-signature status.
- Add top-k/null calibration before PCA-LR reconstructed gene weights.

### V9-ORG-028: Human organoid feature-effect top-k null calibration

Status: complete as of 2026-05-23.

Goal:

Add expected-overlap and null-calibration details to
`feature_effect_top_k_de_overlap`.

Likely work:

- Extend `spacebio_bench/feature_effects.py` top-k overlap output with expected
  overlap, enrichment, and hypergeometric or permutation-style p-value.
- Regenerate `v9/human_organoid/reports/logistic_feature_effect/`.
- Update tests to cover calibrated top-k fields.

Done when:

- The top-k metric reports observed overlap, expected overlap, enrichment, and a
  null-calibration field.
- Targeted and full tests pass.

Actual output:

- `spacebio_bench/feature_effects.py`
- `scripts/run_v9_human_organoid_logistic_feature_effect.py`
- `v9/human_organoid/reports/logistic_feature_effect/metrics.json`
- `v9/human_organoid/reports/logistic_feature_effect/README.md`
- `v9/human_organoid/reports/ORGANOID_FEATURE_EFFECT_NULL_CALIBRATION_REVIEW.md`

Decision:

- Keep calibrated top-k fields in `feature_effect_top_k_de_overlap`.
- Interpret aggregate top-k enrichment as a secondary feature-effect signal.
- Do not promote feature effects to primary benchmark metrics.
- Move next to PCA-LR reconstructed feature-effect design, not immediate
  implementation.

### V9-ORG-029: Human organoid PCA-LR reconstructed feature-effect design

Status: complete as of 2026-05-23.

Goal:

Design a PCA-LR reconstructed gene-space feature-effect adapter, then decide
whether it is worth implementing.

Likely work:

- Review the existing PCA-LR baseline pipeline and identify where train-fold
  scaling, PCA fitting, and logistic fitting occur.
- Define the reconstruction formula from PCA components and logistic
  coefficients back to gene-space effects.
- Preserve the `feature_effect.csv` schema and calibrated top-k diagnostics.
- Require direct comparison against the calibrated L2 logistic feature-effect
  report.
- Define a go/no-go rule before writing a runner.

Done when:

- A design note records leakage boundaries, reconstruction math, output schema,
  and implementation decision.

Actual output:

- `v9/human_organoid/reports/ORGANOID_PCA_LR_FEATURE_EFFECT_DESIGN.md`

Decision:

- Implement one constrained PCA-LR reconstructed feature-effect pilot.
- Reuse `feature_effect.csv`, not `response_signature.csv`.
- Require calibrated comparison against the L2 logistic feature-effect report.

### V9-ORG-030: Human organoid PCA-LR reconstructed feature-effect pilot

Status: complete as of 2026-05-23.

Goal:

Implement the PCA-LR reconstructed gene-space feature-effect adapter described
in V9-ORG-029.

Likely work:

- Add `PCA_LR_FEATURE_EFFECT_ID` and
  `build_pca_lr_reconstructed_feature_effect` to
  `spacebio_bench/feature_effects.py`.
- Add `scripts/run_v9_human_organoid_pca_lr_feature_effect.py`.
- Generate `v9/human_organoid/reports/pca_lr_feature_effect/`.
- Add unit and CLI tests covering reconstruction math and report output.
- Compare calibrated diagnostics against
  `v9/human_organoid/reports/logistic_feature_effect/metrics.json`.

Done when:

- The PCA-LR feature-effect report is generated.
- A comparison note records whether PCA-LR reconstruction improves, matches, or
  weakens the L2 logistic feature-effect diagnostic.
- Targeted and full tests pass.

Actual output:

- `spacebio_bench/feature_effects.py`
- `scripts/run_v9_human_organoid_pca_lr_feature_effect.py`
- `v9/human_organoid/reports/pca_lr_feature_effect/predictions.csv`
- `v9/human_organoid/reports/pca_lr_feature_effect/feature_effect.csv.gz`
- `v9/human_organoid/reports/pca_lr_feature_effect/feature_effect_metadata.json`
- `v9/human_organoid/reports/pca_lr_feature_effect/metrics.json`
- `v9/human_organoid/reports/pca_lr_feature_effect/run_manifest.json`
- `v9/human_organoid/reports/pca_lr_feature_effect/README.md`
- `v9/human_organoid/reports/ORGANOID_PCA_LR_FEATURE_EFFECT_REVIEW.md`

Observed diagnostic result:

- feature-effect rows: 16,000.
- `feature_effect_direction_match`: 0.5980392156862745.
- `feature_effect_rank_correlation`: 0.08668664748698189.
- calibrated top-k overlap: 1/50, 5/100, 10/250, 14/500.

Decision:

- Keep PCA-LR reconstructed feature effects as an optional negative comparison.
- Do not promote PCA-LR over the simpler L2 logistic feature-effect report.
- Do not pursue more PCA-LR reconstruction variants unless a future task uses a
  genuinely low-rank PCA setting.

### V9-ORG-031: Human organoid diagnostic consolidation and release boundary

Status: complete as of 2026-05-23.

Goal:

Consolidate the completed organoid diagnostic family and decide the next active
lane.

Likely work:

- Summarize response-signature diagnostics: global source-transfer,
  microglia-matched, shared-control partial.
- Summarize feature-effect diagnostics: L2 logistic, calibrated top-k, PCA-LR
  reconstructed negative comparison.
- Mark which artifacts are part of the draft v9 alpha surface and which are
  exploratory-only.
- Update `v9/README.md`, this backlog, and the long-run protocol.
- Decide whether to return to multispecies expansion or release packaging next.

Done when:

- A consolidation note exists with the organoid diagnostic release boundary and
  the next active lane decision.

Actual output:

- `v9/human_organoid/reports/ORGANOID_DIAGNOSTIC_CONSOLIDATION_AND_RELEASE_BOUNDARY.md`

Decision:

- Main draft-alpha organoid surface: task/provenance artifacts, derived DE
  reference, global source-transfer response signatures, and calibrated L2
  logistic feature effects.
- Secondary/exploratory: microglia-matched, shared-control, PCA-LR
  reconstructed, donor-holdout, and smoke-fixture reports.
- Return the active scientific lane to multispecies expansion.

### V9-MULTI-001: Non-mouse source inventory

Status: draft inventory complete as of 2026-05-21 for OSD-207, OSD-37, and
OSD-120.

Goal:

Inventory Drosophila and Arabidopsis pilot sources already present in the repo,
then identify additional OSDR species candidates.

Output:

- `v9/multispecies/source_inventory.draft.csv`
- `v9/multispecies/source_inventory.draft.json`

### V9-MULTI-002: Ortholog/pathway feature strategy

Status: initial strategy complete as of 2026-05-21.

Goal:

Decide whether each multispecies task uses raw gene ids, ortholog groups,
pathways, or NES-style signatures.

Output:

- `docs/V9_MULTISPECIES_FEATURE_STRATEGY.md`

Actual decision:

- Species-native fly and plant tasks may use species-local gene ids.
- Cross-species bridge tasks should start with pathway/NES-style features.
- Ortholog-group features are a later sensitivity layer after mapping rules are
  explicit.

### V9-MULTI-003: Multispecies candidate source deep audit

Status: complete as of 2026-05-23.

Goal:

Turn the draft non-mouse source inventory into source-specific task candidates.

Likely work:

- Inspect `v9/multispecies/source_inventory.draft.csv` and `.json`.
- Summarize per-source evidence for OSD-207, OSD-37, and OSD-120.
- Classify each source as species-native classification, response-signature,
  pathway/NES bridge, or defer.
- Identify missing sample-factor, matrix, and checksum evidence before task
  manifests can be generated.

Done when:

- A multispecies source deep-audit note records go/defer decisions and the next
  task-manifest block.

Actual output:

- `v9/multispecies/reports/MULTISPECIES_CANDIDATE_SOURCE_DEEP_AUDIT.md`

Decision:

- OSD-37: first plant go-after-source-audit candidate.
- OSD-207: first fly go-after-source-audit candidate.
- OSD-120: defer to Arabidopsis light/genotype interaction design.
- Do not start with raw-gene cross-species or ortholog bridge tasks.

### V9-MULTI-004: Multispecies sample-factor and matrix audit scaffold

Status: complete as of 2026-05-23.

Goal:

Create the first parsed multispecies sample-factor and matrix-alignment audit
artifacts for OSD-37, OSD-207, and OSD-120.

Likely work:

- Add parser logic for the simple two-column OSD sample tables.
- Extract `source_id`, `sample_id`, `true_label`, organism, genotype/ecotype,
  light treatment where present, and condition stratum.
- Write `v9/multispecies/sample_factors.draft.csv` and `.json`.
- Write `v9/multispecies/expression_matrix_audit.draft.csv` and `.json`.
- Verify count-matrix columns match parsed sample ids.
- Keep task manifests deferred until source checksums and parser output are
  reviewed.

Done when:

- Parsed sample factors and expression matrix audit artifacts exist for all
  three multispecies candidates.

Actual output:

- `v9/multispecies/sample_factors.draft.csv`
- `v9/multispecies/sample_factors.draft.json`
- `v9/multispecies/expression_matrix_audit.draft.csv`
- `v9/multispecies/expression_matrix_audit.draft.json`

Actual result:

- OSD-207: 32 samples, 15,999 feature rows, 32/32 matrix/sample-factor columns
  aligned, 16 Ground and 16 LEO/ISS labels across four Drosophila
  background/variant strata.
- OSD-37: 56 samples, 28,067 feature rows, 56/56 columns aligned, 28 Ground
  and 28 LEO/ISS labels across four Arabidopsis ecotypes.
- OSD-120: 36 samples, 24,740 feature rows, 36/36 columns aligned, 18 Ground
  and 18 LEO/ISS labels across genotype/ecotype by light-treatment strata.
- Task manifests remain deferred until multispecies OSDR checksum evidence is
  generated and reviewed.

### V9-MULTI-005: Multispecies OSDR checksum audit

Status: complete as of 2026-05-23.

Goal:

Generate OSDR file-list and checksum-manifest evidence for OSD-207, OSD-37,
and OSD-120 before drafting species-native task manifests.

Likely work:

- Run `scripts/audit_v9_source_checksums.py` against
  `v9/multispecies/source_inventory.draft.csv`.
- Write `v9/multispecies/source_checksum_audit.draft.csv` and `.json`.
- Compare API/checksum evidence against the already-present local processed
  matrices and SampleTables.
- Record whether each source is task-manifest-ready, checksum-evidence-ready
  but payload-hash-pending, or blocked.
- Keep OSD-120 as interaction-task deferred unless checksum/sample evidence
  changes the decision.

Done when:

- All three multispecies source rows have API/checksum audit status and a
  concise readiness decision.

Actual output:

- `v9/multispecies/source_checksum_audit.draft.csv`
- `v9/multispecies/source_checksum_audit.draft.json`
- `v9/multispecies/reports/MULTISPECIES_CHECKSUM_AND_LOCAL_PAYLOAD_AUDIT.md`

Actual result:

- OSD-207: API OK, 534 OSDR files, 1 checksum manifest, 370 parsed MD5 entries,
  366 payload-name matches; local SampleTable and normalized-count matrix MD5s
  match the OSDR processed manifest.
- OSD-37: API OK, 986 OSDR files, 2 checksum manifests, 927 parsed MD5 entries,
  926 payload-name matches; local SampleTable and normalized-count matrix MD5s
  match the OSDR processed manifest.
- OSD-120: API OK, 851 OSDR files, 1 checksum manifest, 533 parsed MD5 entries,
  532 payload-name matches; local SampleTable and normalized-count matrix MD5s
  match the OSDR processed manifest.
- OSD-37 and OSD-207 can proceed to species-native draft task-manifest design.
  OSD-120 remains deferred to genotype/ecotype by light-treatment interaction
  design.

### V9-MULTI-006: Multispecies species-native task-manifest design

Status: complete as of 2026-05-23.

Goal:

Draft task manifests for the first species-native plant and fly tasks without
mixing species into a raw-gene cross-species benchmark.

Likely work:

- Design one OSD-37 Arabidopsis species-native task with ecotype-aware
  stratification.
- Design one OSD-207 Drosophila species-native task with genotype/background
  strata made explicit.
- Decide split contract: within-source blocked cross-validation versus
  condition-stratum holdout, with a preference for conservative blocked splits
  that do not claim mission generalization.
- Reuse the V9-MULTI-004 sample-factor and matrix audits plus V9-MULTI-005
  checksum evidence in each manifest source/provenance record.
- Keep OSD-120 out of the first manifest set and open a separate interaction
  design note.

Done when:

- Draft JSON manifests exist for OSD-37 and OSD-207.
- A task-manifest index exists for `v9/multispecies/task_manifests/`.
- Tests validate both manifests and their source/sample-count boundaries.

Actual output:

- `v9/multispecies/task_manifests/draft_osd37_arabidopsis_seedling_spaceflight.json`
- `v9/multispecies/task_manifests/draft_osd207_drosophila_whole_body_spaceflight.json`
- `v9/multispecies/task_manifest_index.draft.csv`
- `v9/multispecies/task_manifest_index.draft.json`
- `v9/multispecies/reports/MULTISPECIES_SPECIES_NATIVE_TASK_MANIFEST_DESIGN.md`

Actual result:

- OSD-37 manifest records 56 samples, four ecotype condition-stratum candidate
  folds, aligned matrix evidence, checksum-manifest evidence, and local
  SampleTable/matrix MD5 match status.
- OSD-207 manifest records 32 samples, four Drosophila background/variant
  condition-stratum candidate folds, aligned matrix evidence,
  checksum-manifest evidence, and local SampleTable/matrix MD5 match status.
- OSD-120 remains explicitly deferred to a separate genotype/ecotype by
  light-treatment interaction-task design.

### V9-MULTI-007: Multispecies species-native baseline feasibility

Status: complete as of 2026-05-23.

Goal:

Run conservative draft-only baseline feasibility checks for the two
species-native multispecies manifests.

Likely work:

- Add a multispecies task loader for the audited local normalized-count
  matrices and sample-factor rows.
- Align matrix columns to manifest sample ids for OSD-37 and OSD-207.
- Reuse simple train-only preprocessing patterns from organoid/bulk baselines
  where appropriate, without creating a leaderboard claim.
- Run at least one baseline that respects condition-stratum candidate folds.
- Write draft reports under `v9/multispecies/reports/`.

Done when:

- OSD-37 and OSD-207 have draft baseline reports or a documented blocker.
- Tests cover the loader and at least one smoke baseline path.

Actual output:

- `spacebio_bench/data/multispecies.py`
- `spacebio_bench/baselines/multispecies.py`
- `scripts/run_v9_multispecies_baseline.py`
- `v9/multispecies/reports/nearest_centroid/`
- `v9/multispecies/reports/MULTISPECIES_BASELINE_FEASIBILITY_REVIEW.md`

Actual result:

- OSD-37 nearest-centroid feasibility: 56 predictions,
  `balanced_accuracy=0.8214285714`, `auroc=0.9196428571`,
  `condition_stratum_holdout_delta=0.125`.
- OSD-207 nearest-centroid feasibility: 32 predictions,
  `balanced_accuracy=0.875`, `auroc=0.9296875`,
  `condition_stratum_holdout_delta=0.375`.
- OSD-37 is the cleaner first plant feasibility example. OSD-207 remains valid
  but needs sensitivity review because the `w1118_KCNQ370` heldout stratum is
  weaker than the other fly folds.
- `mission_discrimination` is skipped by design because these are
  single-source/single-mission-context species-native pilots.

### V9-MULTI-008: Multispecies baseline sensitivity grid

Status: complete as of 2026-05-23.

Goal:

Run a small preprocessing sensitivity grid for the two draft multispecies
species-native baselines.

Likely work:

- Add a multispecies sensitivity runner mirroring the organoid sensitivity grid.
- Evaluate `log1p` versus no transform, `zscore` versus no scaling, and top
  100/500/2,000/5,000/all variable-gene settings where feasible.
- Summarize whether OSD-37 stays stable across preprocessing choices.
- Check whether OSD-207's `w1118_KCNQ370` fold remains the fragile stratum.
- Keep all outputs under a draft-only, not-leaderboard claim boundary.

Done when:

- A sensitivity summary exists under `v9/multispecies/reports/sensitivity/`.
- A short review decides whether the default nearest-centroid setting remains
  acceptable for the next multispecies checkpoint.

Actual output:

- `scripts/run_v9_multispecies_sensitivity.py`
- `v9/multispecies/reports/sensitivity/multispecies_baseline_summary.csv`
- `v9/multispecies/reports/sensitivity/multispecies_baseline_summary.json`
- `v9/multispecies/reports/MULTISPECIES_BASELINE_SENSITIVITY_REVIEW.md`

Actual result:

- OSD-37 sensitivity range: balanced accuracy 0.5357142857 to 0.8571428571,
  median 0.7678571429. The conservative default remains acceptable because it
  balances score and fold stability.
- OSD-207 sensitivity range: balanced accuracy 0.75 to 0.9375, median 0.875.
  `w1118_KCNQ370` is the weakest heldout stratum in 18/20 variants.
- Keep `log1p`, train-fold `zscore`, and top 2,000 train-fold variable genes as
  the default baseline setting for the next checkpoint.

### V9-MULTI-009: OSD-120 interaction-task design

Status: complete as of 2026-05-23.

Goal:

Design the deferred OSD-120 Arabidopsis root task as a genotype/ecotype by
light-treatment interaction benchmark, without merging it into the simpler
OSD-37 species-native plant task.

Likely work:

- Review OSD-120 sample factors, matrix audit, and checksum evidence.
- Decide whether the primary split should hold out genotype/ecotype,
  light-treatment, or genotype/ecotype-by-light strata.
- Define which metric should represent interaction robustness rather than
  simple binary LEO/Ground classification.
- Draft a design note before generating a manifest.
- Generate a manifest only if each proposed fold has enough Ground and LEO/ISS
  samples to evaluate.

Done when:

- A design note records whether OSD-120 should become a draft manifest now or
  remain deferred.

Actual output:

- `v9/multispecies/reports/OSD120_INTERACTION_TASK_DESIGN.md`

Actual result:

- OSD-120 has a complete 3 genotype/ecotype x 2 light-treatment x 2 label x 3
  replicate structure.
- Genotype/ecotype holdout, light-treatment holdout, and genotype/ecotype by
  light condition-stratum holdout all have balanced train/test labels.
- Decision: OSD-120 should become a separate draft interaction task family,
  tentatively `multispecies_light_interaction_spaceflight`, instead of being
  merged into the OSD-37 species-native task.

### V9-MULTI-010: OSD-120 interaction manifest implementation

Status: complete as of 2026-05-23.

Goal:

Implement the OSD-120 interaction manifest in a separate manifest directory and
keep it separate from the OSD-37/OSD-207 species-native manifests.

Likely work:

- Add an OSD-120 interaction manifest builder.
- Write
  `v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json`.
- Write an interaction task-manifest index.
- Add tests validating primary genotype/ecotype folds, secondary light folds,
  and tertiary condition-stratum diagnostic folds.
- Decide whether to extend the existing multispecies loader or create a
  separate interaction loader in the next baseline block.

Done when:

- The OSD-120 manifest validates and has sample-count-backed fold metadata for
  genotype/ecotype, light-treatment, and condition-stratum holdouts.

Actual output:

- `spacebio_bench.extension_tasks.build_osd120_interaction_task_manifest`
  builds the separate interaction-task manifest from the multispecies source
  inventory, sample factors, expression audit, and checksum audit.
- `scripts/build_v9_osd120_interaction_task_manifest.py` writes the manifest and
  draft interaction index.
- `v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json`
  records 36 samples, three genotype/ecotype primary folds, two light-treatment
  secondary folds, and six condition-stratum diagnostic folds.

### V9-MULTI-011: OSD-120 interaction loader and baseline feasibility

Status: complete as of 2026-05-23.

Goal:

Load the OSD-120 interaction manifest against the audited local normalized-count
matrix and run a draft-only feasibility baseline without mixing it into the
OSD-37/OSD-207 species-native baseline report.

Likely work:

- Decide whether the existing multispecies loader can support the interaction
  manifest directly or needs a small interaction-specific adapter.
- Align OSD-120 expression columns to the new manifest fold metadata.
- Run nearest-centroid feasibility across genotype/ecotype primary folds first,
  then light-treatment and condition-stratum diagnostics as secondary outputs.
- Write a review note under `v9/multispecies/reports/` that interprets primary
  versus diagnostic folds and preserves the draft-only claim boundary.

Done when:

- A reproducible OSD-120 interaction baseline report exists with run manifests,
  predictions, metrics, and a conservative review note.

Actual output:

- `spacebio_bench.data.load_multispecies_interaction_task` loads the OSD-120
  interaction manifest, aligns 36 samples to 24,740 expression features, and
  exposes three fold families: 3 primary genotype/ecotype folds, 2 secondary
  light-treatment folds, and 6 diagnostic condition-stratum folds.
- `scripts/run_v9_osd120_interaction_baseline.py` runs the draft interaction
  nearest-centroid baseline and writes one report per fold family.
- `v9/multispecies/reports/interaction_nearest_centroid/multispecies_baseline_summary.csv`
  records three evaluated fold-family rows. Balanced accuracy is 0.6667 for all
  three families; AUROC is 0.7346 for primary genotype/ecotype, 0.7840 for
  secondary light-treatment, and 0.7654 for diagnostic condition-stratum.
- `v9/multispecies/reports/OSD120_INTERACTION_BASELINE_FEASIBILITY_REVIEW.md`
  keeps the report draft-only and notes fold-level fragility in
  `Col.0.PhyD`, `Light.Treatment`, and selected dark condition strata.

### V9-MULTI-012: OSD-120 interaction sensitivity grid

Status: complete as of 2026-05-23.

Goal:

Check whether the OSD-120 interaction baseline's fold-family and stratum-level
behavior is stable under preprocessing choices before any stronger baseline or
claim language is considered.

Likely work:

- Add an OSD-120 interaction sensitivity runner mirroring the species-native
  multispecies sensitivity grid.
- Evaluate log1p/no-transform, z-score/no-scaling, and multiple train-fold
  variable-gene settings across the three interaction fold families.
- Summarize which fold family and held-out strata are repeatedly fragile.
- Write a sensitivity review under `v9/multispecies/reports/`.

Done when:

- A sensitivity summary exists under
  `v9/multispecies/reports/interaction_sensitivity/`, and a review states
  whether the current `log1p` + train-fold `zscore` + top 2,000 genes default
  remains acceptable for OSD-120.

Actual output:

- `scripts/run_v9_osd120_interaction_sensitivity.py` runs the OSD-120
  interaction sensitivity grid.
- `v9/multispecies/reports/interaction_sensitivity/multispecies_baseline_summary.csv`
  contains 60 evaluated rows: 20 preprocessing variants for each of the three
  interaction fold families.
- `v9/multispecies/reports/OSD120_INTERACTION_SENSITIVITY_REVIEW.md` records
  mean balanced accuracy ranges by fold family, repeated fragile strata, and the
  decision to keep the conservative default for now.

### V9-MULTI-013: OSD-120 interaction fold-detail aggregation

Status: complete as of 2026-05-23.

Goal:

Make the OSD-120 sensitivity grid's repeated fragile strata machine-readable
without requiring readers or scripts to open 60 individual `metrics.json` files.

Likely work:

- Add a small aggregation script for OSD-120 interaction sensitivity fold
  details.
- Read the interaction sensitivity summary plus each metrics file.
- Write a fold-detail table with fold family, held-out factor/value, variant id,
  balanced accuracy, rank within variant, and fragile-stratum flags.
- Write a compact review/update that points to the machine-readable table.

Done when:

- `v9/multispecies/reports/interaction_sensitivity/fold_detail_summary.csv`
  and `.json` exist and tests validate row counts and key repeated fragile
  strata.

Actual output:

- `scripts/build_v9_osd120_interaction_fold_details.py` builds the fold-detail
  aggregation table from the interaction sensitivity summary and nested
  metrics reports.
- `v9/multispecies/reports/interaction_sensitivity/fold_detail_summary.csv`
  contains 220 rows: 60 primary genotype/ecotype rows, 40 secondary
  light-treatment rows, and 120 diagnostic condition-stratum rows.
- The generated table confirms the most repeated lowest-performing held-out
  values: `Light.Treatment` in 19/20 secondary variants,
  `Wassilewskija.ecotype|Dark.Treatment` in 16/20 diagnostic variants, and
  `Wassilewskija.ecotype` in 12/20 primary variants.

### V9-MULTI-014: OSD-120 interaction L2 logistic baseline design

Status: complete as of 2026-05-25.

Goal:

Design the first stronger but still transparent OSD-120 interaction baseline
without erasing the fold-family separation or draft-only claim boundary.

Likely work:

- Review whether the existing sklearn baseline machinery can be adapted to the
  OSD-120 interaction task.
- Specify the L2 logistic preprocessing, feature selection, solver, and output
  schema.
- Require fold-family separated reports and reuse
  `fold_detail_summary.csv` for fragile-stratum comparison.
- Write a design note before implementation if the adaptation is not trivial.

Done when:

- A short design/review note states whether to implement L2 logistic directly
  next, which files would change, and how the nearest-centroid sensitivity
  table will be used as the comparison baseline.

Actual output:

- `v9/multispecies/reports/OSD120_INTERACTION_LOGISTIC_BASELINE_DESIGN.md`
  records the implementation decision and model contract.
- `scripts/run_v9_osd120_interaction_logistic.py` runs the draft OSD-120 L2
  logistic baseline over all three interaction fold families.
- `v9/multispecies/reports/interaction_logistic_l2/multispecies_baseline_summary.csv`
  records three evaluated fold-family rows. Balanced accuracy is 0.7778 for
  primary genotype/ecotype, 0.7500 for secondary light-treatment, and 0.8611
  for diagnostic condition-stratum.
- `v9/multispecies/reports/OSD120_INTERACTION_LOGISTIC_BASELINE_REVIEW.md`
  compares logistic against nearest centroid and keeps the report diagnostic.

### V9-MULTI-015: OSD-120 logistic fold-detail comparison

Status: complete as of 2026-05-25.

Goal:

Compare L2 logistic fold-level behavior against the nearest-centroid
fold-detail table so stronger aggregate scores do not obscure repeated fragile
strata.

Likely work:

- Add a fold-detail aggregation path for `interaction_logistic_l2`.
- Write logistic fold-detail CSV/JSON tables.
- Join nearest-centroid sensitivity/default fold detail with logistic fold
  detail by fold family and held-out value.
- Write a comparison review focused on `Col.0.PhyD|Dark.Treatment`,
  `Light.Treatment`, and primary genotype/ecotype fragility.

Done when:

- A machine-readable logistic fold-detail table and side-by-side comparison note
  exist, with tests validating row counts and key weak strata.

Actual output:

- `scripts/build_v9_osd120_logistic_fold_comparison.py` rebuilds the logistic
  fold-detail table and side-by-side comparison from the generated logistic
  summary and nearest-centroid fold-detail summary.
- `v9/multispecies/reports/interaction_logistic_l2/fold_detail_summary.csv`
  and `.json` contain 11 logistic held-out fold rows.
- `v9/multispecies/reports/interaction_logistic_l2/fold_detail_comparison_vs_nearest_centroid.csv`
  and `.json` compare those rows against the default nearest-centroid
  `tvg2000_log1p_zscore` rows.
- `v9/multispecies/reports/OSD120_INTERACTION_LOGISTIC_BASELINE_REVIEW.md`
  records the V9-MULTI-015 conclusion: logistic improves 8/11 held-out folds,
  ties 2/11, and worsens `Col.0.PhyD|Dark.Treatment`.
- Tests cover the API, CLI, and generated output row counts/key weak strata.

### V9-MULTI-016: OSD-120 logistic regularization sensitivity

Status: complete as of 2026-05-25.

Goal:

Check whether the V9-MULTI-015 fragile fold is a stable biological/modeling
warning or an artifact of a single logistic configuration.

Likely work:

- Add a small OSD-120 logistic sensitivity grid over `C` and top-variable-gene
  counts while preserving train-fold-only preprocessing.
- Keep fold families separated and reuse the V9-MULTI-015 comparison table
  shape for side-by-side fragile-fold review.
- Focus first on `Col.0.PhyD|Dark.Treatment`, `Col.0.PhyD`, and
  `Light.Treatment`.
- Keep the claim boundary draft-only and do not promote logistic outputs into a
  benchmark leaderboard.

Done when:

- A compact logistic regularization/feature-count sensitivity table exists with
  tests and a short review identifying whether the worsened fold persists.

Actual output:

- `scripts/run_v9_osd120_interaction_logistic_sensitivity.py` runs the compact
  top-gene/C grid and writes summary, fold-detail, and nearest-centroid
  comparison tables.
- `v9/multispecies/reports/interaction_logistic_l2_sensitivity/multispecies_baseline_summary.csv`
  contains 18 evaluated rows: six logistic variants by three fold families.
- `v9/multispecies/reports/interaction_logistic_l2_sensitivity/fold_detail_summary.csv`
  contains 66 logistic held-out fold rows.
- `v9/multispecies/reports/interaction_logistic_l2_sensitivity/fold_detail_comparison_vs_nearest_centroid.csv`
  contains 66 nearest-centroid comparison rows.
- `v9/multispecies/reports/OSD120_INTERACTION_LOGISTIC_SENSITIVITY_REVIEW.md`
  records the conclusion: top 2,000 genes preserves stronger light-treatment
  behavior but keeps `Col.0.PhyD|Dark.Treatment` at 0.3333, while top 500 genes
  restores that dark stratum to 0.6667 but weakens `Light.Treatment`.

### V9-MULTI-017: OSD-120 logistic feature-set audit

Status: complete as of 2026-05-25.

Goal:

Explain the top-500 versus top-2,000 tradeoff from V9-MULTI-016 before adding a
more complex classifier.

Likely work:

- Compare train-fold selected feature sets for top 500 and top 2,000 logistic
  variants, with special attention to `Col.0.PhyD|Dark.Treatment`,
  `Light.Treatment`, and `Col.0.PhyD`.
- Extract compact fold-level coefficient/effect summaries for the two variant
  families without treating coefficients as biological DE claims.
- Identify whether extra 1,500 genes mainly change the problematic dark stratum
  or introduce broader light-treatment instability.
- Keep outputs diagnostic and draft-only.

Done when:

- A feature-set/coefficient audit table and review note explain the observed
  V9-MULTI-016 tradeoff well enough to decide whether to try a non-linear model
  or adjust the transparent baseline family.

Actual output:

- `scripts/audit_v9_osd120_logistic_feature_sets.py` writes selected-feature
  summary and coefficient audit tables for the focused OSD-120 folds.
- `v9/multispecies/reports/interaction_logistic_feature_audit/feature_set_audit_summary.csv`
  contains six summary rows: three focus folds by top-500/top-2,000 logistic
  variants.
- `v9/multispecies/reports/interaction_logistic_feature_audit/feature_coefficient_audit.csv`
  contains 7,500 selected-feature coefficient rows.
- `v9/multispecies/reports/OSD120_INTERACTION_LOGISTIC_FEATURE_SET_AUDIT_REVIEW.md`
  records that top-10 absolute coefficient overlap between top-500 and
  top-2,000 variants is only 1/10 to 3/10 across the focus folds.

### V9-MULTI-018: OSD-120 sparse feature-stability branch design

Status: complete as of 2026-05-25.

Goal:

Design the next OSD-120 modeling branch using the V9-MULTI-017 feature-set
audit, without jumping straight to opaque model claims.

Likely work:

- Decide whether the next branch should be sparse logistic, stability-selected
  logistic, or a small non-linear diagnostic.
- Define gates that any branch must satisfy: improve `Light.Treatment`, avoid
  recreating the `Col.0.PhyD|Dark.Treatment` failure, and keep fold-family
  separated reporting.
- Preserve the current nearest-centroid and L2 logistic outputs as comparator
  baselines.
- Keep the work draft-only until all fragile folds are reviewed side by side.

Done when:

- A short branch design note exists with exact candidate configs, outputs,
  decision gates, and tests to implement before running the next model.

Actual output:

- `v9/multispecies/reports/OSD120_INTERACTION_SPARSE_BRANCH_DESIGN.md`
  records the external method context, local evidence, candidate L1 configs,
  fragile-fold gates, output contract, and test plan.
- The selected next branch is a draft-only L1 logistic pilot over top 2,000
  genes with `C` values 0.01, 0.03, 0.1, 0.3, and 1.0.
- Elastic-net/stability-selection remain fallback paths only after the simpler
  sparse pilot passes or fails the fragile-fold gates.

### V9-MULTI-019: OSD-120 sparse L1 logistic pilot

Status: complete as of 2026-05-25.

Goal:

Implement and run the sparse L1 logistic pilot designed in V9-MULTI-018.

Likely work:

- Add sparse logistic config support without breaking existing L2 summaries.
- Run top-2,000 L1 logistic over `C` values 0.01, 0.03, 0.1, 0.3, and 1.0 for
  all three OSD-120 fold families.
- Write fold-detail and nearest-centroid comparison tables.
- Write a sparse coefficient audit with non-zero feature counts per fold.
- Review the candidate variants against the fragile-fold gates.

Done when:

- Sparse L1 reports, fold-detail comparison, coefficient audit, tests, and a
  review note determine whether any L1 variant passes the V9-MULTI-018 gates.

Actual output:

- `scripts/run_v9_osd120_interaction_sparse_l1.py` runs the sparse L1 pilot,
  fold-detail comparison, and sparse coefficient audit.
- `v9/multispecies/reports/interaction_logistic_sparse_l1/multispecies_baseline_summary.csv`
  contains 15 evaluated rows: five L1 `C` values by three fold families.
- `v9/multispecies/reports/interaction_logistic_sparse_l1/fold_detail_summary.csv`
  contains 55 sparse L1 held-out fold rows.
- `v9/multispecies/reports/interaction_logistic_sparse_l1/feature_set_audit_summary.csv`
  contains 15 focus-fold coefficient-audit summary rows, and
  `feature_coefficient_audit.csv` contains 30,000 selected-feature rows.
- `v9/multispecies/reports/OSD120_INTERACTION_SPARSE_L1_REVIEW.md` selects
  `tvg2000_log1p_zscore_l1_c1` as the current best transparent OSD-120
  diagnostic candidate: focus folds pass, 9/11 nearest-centroid fold rows
  improve, 2/11 tie, and 0/11 worsen.

### V9-MULTI-020: OSD-120 sparse L1 stability audit

Status: complete as of 2026-05-26.

Goal:

Test whether the useful sparse L1 region from V9-MULTI-019 has stable selected
features under train-fold subsampling.

Likely work:

- Focus on `C=0.3` and `C=1.0`.
- Refit sparse L1 models on deterministic balanced subsamples within each
  training fold.
- Record feature selection frequencies and non-zero coefficient sign stability
  for `Col.0.PhyD|Dark.Treatment`, `Light.Treatment`, and `Col.0.PhyD`.
- Keep performance evaluation on held-out folds separate from subsampling-based
  feature stability.

Done when:

- A stability audit table and review note determine whether the sparse L1
  candidate is stable enough to remain the leading transparent OSD-120
  comparator.

Actual output:

- `scripts/audit_v9_osd120_sparse_l1_stability.py` writes deterministic
  train-fold subsampling stability tables for the sparse L1 focus variants.
- `v9/multispecies/reports/interaction_logistic_sparse_l1_stability/stability_summary.csv`
  contains six rows: two L1 `C` values by three focus folds.
- `v9/multispecies/reports/interaction_logistic_sparse_l1_stability/stability_feature_frequency.csv`
  contains 12,000 feature-frequency rows.
- `v9/multispecies/reports/OSD120_INTERACTION_SPARSE_L1_STABILITY_REVIEW.md`
  records that `C=1.0` is stronger for `Col.0.PhyD|Dark.Treatment` and
  `Col.0.PhyD`, while `C=0.3` is the compact light-treatment stability
  comparator.

### V9-MULTI-021: OSD-120 interaction baseline ladder consolidation

Status: complete as of 2026-05-26.

Goal:

Consolidate the OSD-120 interaction modeling ladder into one comparison note
before adding any new model family.

Likely work:

- Compare nearest-centroid default, L2 logistic default, L2 top-500 sensitivity,
  sparse L1 `C=0.3`, and sparse L1 `C=1.0`.
- Summarize fold-family aggregate metrics, fragile-focus fold behavior, and
  feature-stability evidence.
- Decide which variant becomes the leading transparent OSD-120 diagnostic
  comparator and which variants remain sensitivity controls.
- Keep all language draft-only and avoid leaderboard promotion.

Done when:

- A single baseline-ladder review table exists and points the next block toward
  either finalizing the OSD-120 draft diagnostic surface or testing a new model
  family with clear gates.

Actual output:

- `scripts/build_v9_osd120_baseline_ladder.py` writes the OSD-120 ladder
  consolidation tables.
- `v9/multispecies/reports/interaction_baseline_ladder/baseline_ladder_summary.csv`
  contains five candidate rows: nearest centroid, dense L2 default, top-500 L2,
  sparse L1 `C=0.3`, and sparse L1 `C=1.0`.
- `v9/multispecies/reports/interaction_baseline_ladder/baseline_ladder_focus_folds.csv`
  contains 15 focus-fold rows for `Col.0.PhyD|Dark.Treatment`,
  `Light.Treatment`, and `Col.0.PhyD` across those five candidates.
- `v9/multispecies/reports/OSD120_INTERACTION_BASELINE_LADDER_REVIEW.md`
  selects `sparse_l1_c1` as the primary draft transparent diagnostic candidate
  and keeps `sparse_l1_c0p3` as the compact stability comparator.

### V9-MULTI-022: OSD-120 diagnostic surface candidate package

Status: complete as of 2026-05-26.

Goal:

Package the selected OSD-120 sparse diagnostic surface into figure-ready and
claim-boundary-ready artifacts without adding a new model family.

Likely work:

- Build a concise figure-ready candidate table from the ladder, sparse L1
  coefficient audit, and stability audit.
- Add a candidate-package review note that separates evidence, limitations,
  and draft-only claims.
- Define the exact artifact set that would travel into a manuscript or poster
  appendix.
- Keep nearest centroid, dense L2, and sparse L1 `C=0.3` as controls.

Done when:

- A candidate package links each OSD-120 diagnostic statement to the underlying
  manifest, summary table, fold-detail comparison, coefficient audit, stability
  audit, and test evidence.

Actual output:

- `scripts/build_v9_osd120_diagnostic_candidate_package.py` writes the
  OSD-120 sparse diagnostic candidate package.
- `v9/multispecies/reports/interaction_diagnostic_candidate_package/candidate_package_summary.csv`
  contains one package row for `sparse_l1_c1`.
- `v9/multispecies/reports/interaction_diagnostic_candidate_package/candidate_focus_evidence.csv`
  contains three fragile-focus evidence rows.
- `v9/multispecies/reports/interaction_diagnostic_candidate_package/candidate_stable_feature_evidence.csv`
  contains 19 stable-feature evidence rows selected in at least half of
  deterministic balanced train-fold subsamples.
- `v9/multispecies/reports/interaction_diagnostic_candidate_package/candidate_claim_map.csv`
  contains five claim-boundary rows, including an external-context row for
  OSD-120 CARA light/genotype biology.
- `v9/multispecies/reports/OSD120_INTERACTION_DIAGNOSTIC_CANDIDATE_PACKAGE_REVIEW.md`
  records supported and unsupported claims for using the package in a
  manuscript or poster draft.

### V9-MULTI-023: OSD-120 figure/table draft package

Status: complete as of 2026-05-26.

Goal:

Turn the candidate package into a human-facing figure/table draft suitable for
poster, manuscript supplement, or ASGSR-style presentation material.

Likely work:

- Generate a compact markdown or CSV figure table from
  `candidate_package_summary.csv` and `candidate_focus_evidence.csv`.
- Add a stable-feature appendix table with conservative labels.
- Add caption text that states the draft-only boundary and within-source split.
- Decide whether the compact `sparse_l1_c0p3` comparator needs a paired panel.

Done when:

- A figure/table draft exists with caption, source artifact links, and a clear
  list of claims that are allowed versus not allowed.

Actual output:

- `scripts/build_v9_osd120_figure_table_package.py` writes the human-facing
  OSD-120 figure/table draft package.
- `v9/multispecies/reports/interaction_figure_table_package/figure_main_focus_table.csv`
  contains three display-ready focus rows for `sparse_l1_c1`.
- `v9/multispecies/reports/interaction_figure_table_package/figure_stable_feature_appendix.csv`
  contains 19 stable sparse-feature appendix rows.
- `v9/multispecies/reports/interaction_figure_table_package/figure_caption.md`
  provides draft caption text with allowed and disallowed claims.
- `v9/multispecies/reports/interaction_figure_table_package/figure_claim_boundary_box.md`
  provides the short claim guardrail box.
- `v9/multispecies/reports/OSD120_INTERACTION_FIGURE_TABLE_DRAFT_REVIEW.md`
  records that the main focus table is ready as the primary draft table surface.

### V9-MULTI-024: OSD-120 compact comparator paired figure table

Status: complete as of 2026-05-26.

Goal:

Decide whether the compact sparse L1 `C=0.3` comparator should be shown beside
the primary `C=1.0` candidate, and if useful, emit a paired table that makes the
performance-versus-stability tradeoff explicit.

Likely work:

- Reuse the ladder, candidate package, and stability outputs.
- Build a paired focus table for `sparse_l1_c1` versus `sparse_l1_c0p3`.
- Keep the paired table optional and clearly marked as a comparator/control.
- Decide whether the next public-facing table needs the paired comparator or
  should stay with the simpler primary-candidate table.

Done when:

- A paired comparator review states whether `sparse_l1_c0p3` should appear in
  the figure/table set or remain only in the ladder and stability review.

Actual output:

- `scripts/build_v9_osd120_paired_comparator_table.py` writes the optional
  sparse L1 `C=1.0` versus `C=0.3` comparator tables.
- `v9/multispecies/reports/interaction_paired_comparator_table/paired_comparator_summary.csv`
  contains one paired-summary row.
- `v9/multispecies/reports/interaction_paired_comparator_table/paired_focus_comparator_table.csv`
  contains three focus-fold paired-comparator rows.
- `v9/multispecies/reports/interaction_paired_comparator_table/paired_comparator_decision.md`
  recommends appendix/supplement placement for the compact comparator.
- `v9/multispecies/reports/OSD120_INTERACTION_PAIRED_COMPARATOR_REVIEW.md`
  records the decision: keep `sparse_l1_c1` as the main figure candidate and
  use `sparse_l1_c0p3` only as an optional compactness comparator.

### V9-MULTI-025: OSD-120 diagnostic artifact manifest

Status: complete as of 2026-05-26.

Goal:

Index the OSD-120 diagnostic evidence set as one traceable artifact manifest
before moving to another model family or source.

Likely work:

- Build a machine-readable manifest for the task manifest, ladder, sparse L1
  audit, stability audit, candidate package, figure/table package, paired
  comparator, review notes, and tests.
- Include claim boundary, source URLs, generated file counts, and validation
  commands.
- Keep draft-only language and avoid release/freeze claims.

Done when:

- A single OSD-120 diagnostic artifact manifest can answer which files support
  each local claim and which tests validate the current package.

Actual output:

- `scripts/build_v9_osd120_diagnostic_artifact_manifest.py` writes the
  OSD-120 diagnostic artifact manifest and claim map.
- `v9/multispecies/reports/interaction_diagnostic_artifact_manifest/diagnostic_artifact_manifest.csv`
  contains 26 artifact rows with existence status, byte size, line count, row
  count where applicable, SHA-256, claim ids, validation scope, and generator.
- `v9/multispecies/reports/interaction_diagnostic_artifact_manifest/diagnostic_claim_artifact_map.csv`
  contains seven claim rows linking local claims to artifact ids, paths,
  validation tests, external context URLs, and limitations.
- `v9/multispecies/reports/OSD120_INTERACTION_DIAGNOSTIC_ARTIFACT_MANIFEST_REVIEW.md`
  records that the OSD-120 diagnostic package is internally traceable but still
  draft-only.

### V9-MULTI-026: OSD-120 release-readiness gap audit

Status: complete as of 2026-05-26.

Goal:

Check which gaps remain before the OSD-120 diagnostic package could move from a
draft evidence set toward a cleaner public alpha artifact.

Likely work:

- Use the artifact manifest, task manifest, checksum audit, source inventory,
  and local payload audit.
- Identify missing freeze, source checksum, payload hash, data-card,
  reproducibility, and claim-language gaps.
- Separate "must fix before public alpha" from "acceptable draft limitation".

Done when:

- A release-readiness gap table and review note state the blockers, acceptable
  limitations, and next concrete action for OSD-120.

Actual output:

- `scripts/build_v9_osd120_release_readiness_gap_audit.py` writes the
  OSD-120 release-readiness audit tables.
- `v9/multispecies/reports/interaction_release_readiness_gap_audit/release_readiness_summary.csv`
  contains one summary row with
  `public_alpha_ready=False`, `blocker_count=3`, `needs_work_count=3`,
  `pass_count=5`, and `acceptable_draft_limitation_count=1`.
- `v9/multispecies/reports/interaction_release_readiness_gap_audit/release_readiness_gap_table.csv`
  contains 12 gap rows. The public-alpha blockers are full OSDR payload
  freeze, source release-target promotion, and an OSD-120-specific public
  card/citation draft.
- `v9/multispecies/reports/interaction_release_readiness_gap_audit/release_readiness_external_references.csv`
  records seven external release-reference anchors: NASA OSDR, OSDR FAQ, FAIR,
  RO-Crate, DataCite, Hugging Face dataset cards, and GitHub/Zenodo citation
  guidance.
- `v9/multispecies/reports/OSD120_INTERACTION_RELEASE_READINESS_GAP_AUDIT.md`
  records the decision: internally traceable, but not public-alpha ready.

### V9-MULTI-027: OSD-120 payload freeze manifest

Status: complete as of 2026-05-26.

Goal:

Promote OSD-120 from checksum-manifest evidence plus local two-file spot-check
to an explicit payload freeze manifest suitable for public-alpha promotion.

Likely work:

- Start from `v9/multispecies/source_checksum_audit.draft.csv` and the OSD-120
  processed checksum manifest URL.
- Enumerate expected OSD-120 processed payload rows from the parsed MD5
  manifest.
- Hash local OSD-120 payload files that are actually used by the diagnostic
  package.
- Record missing/not-downloaded payload rows separately from hash matches.
- Keep `release_status=draft_not_frozen` unless the freeze manifest proves the
  required payload scope.

Done when:

- A payload freeze table separates matched local files, missing payloads, and
  out-of-scope OSDR payloads, with a review note saying whether the source
  freeze blocker is resolved or still open.

Actual output:

- `scripts/build_v9_osd120_payload_freeze_manifest.py` writes the
  OSD-120 diagnostic payload freeze summary and manifest.
- `v9/multispecies/reports/interaction_payload_freeze_manifest/payload_freeze_summary.csv`
  contains one summary row with 533 parsed OSDR processed checksum entries, 2
  diagnostic-required payloads, 2 MD5 matches, 0 missing required payloads, 0
  required checksum mismatches, and 531 out-of-scope processed payloads.
- `v9/multispecies/reports/interaction_payload_freeze_manifest/payload_freeze_manifest.csv`
  contains one row per parsed OSDR processed checksum entry. The local
  `SampleTable_GLbulkRNAseq.csv` and `Normalized_Counts_GLbulkRNAseq.csv`
  entries are marked `required_payload_md5_matched`.
- `v9/multispecies/reports/OSD120_INTERACTION_PAYLOAD_FREEZE_MANIFEST_REVIEW.md`
  records the decision: the diagnostic-required payload scope is frozen, but a
  full OSDR processed payload mirror is still not claimed.

### V9-MULTI-028: OSD-120 diagnostic public-alpha card draft

Status: complete as of 2026-05-26.

Goal:

Draft an OSD-120-specific public-alpha card that external readers can use
without confusing the diagnostic package with a frozen full benchmark release.

Likely work:

- Start from the claim map, artifact manifest, release-readiness gap table, and
  payload freeze manifest.
- Include OSDR credit/citation language, source/task scope, input payload
  freeze boundary, artifact structure, metrics, allowed claims, disallowed
  claims, and limitations.
- Keep the wording explicit that only the two diagnostic-required processed
  payloads are MD5 matched locally.

Done when:

- A card draft can be used as the source-specific public alpha README/card for
  OSD-120, while preserving draft/not-frozen and diagnostic-only boundaries.

Actual output:

- `scripts/build_v9_osd120_public_alpha_card.py` writes the OSD-120
  diagnostic public-alpha card draft and summary.
- `v9/multispecies/reports/interaction_public_alpha_card/public_alpha_card.md`
  records source scope, payload-freeze boundary, diagnostic result surface,
  allowed claims, disallowed claims, inspectable files, external context links,
  remaining release work, and OSDR/GeneLab credit language.
- `v9/multispecies/reports/interaction_public_alpha_card/public_alpha_card_summary.csv`
  contains one summary row with card status, candidate ids, payload freeze
  decision, artifact/claim counts, diagnostic metrics, and the next required
  block.
- `v9/multispecies/reports/OSD120_INTERACTION_PUBLIC_ALPHA_CARD_REVIEW.md`
  records the decision: this is a source-specific diagnostic alpha card, not a
  frozen benchmark release.

### V9-MULTI-029: OSD-120 rebuild gate and environment lock

Status: complete as of 2026-05-26.

Goal:

Make the OSD-120 diagnostic package reproducible through a single rebuild or
preflight gate with explicit runtime/package context.

Likely work:

- Start from all OSD-120 builder scripts from V9-MULTI-021 through
  V9-MULTI-028.
- Add a small manifest or shell/Python preflight that runs the non-model
  packaging steps in order and records Python/package versions.
- Avoid rerunning expensive sparse L1 model grids unless explicitly requested;
  the first gate should verify and rebuild packaging artifacts from existing
  model outputs.
- Update the public-alpha card or review note with the rebuild command.

Done when:

- One command can rebuild the OSD-120 diagnostic packaging layer and emit a
  run/rebuild manifest with command, status, runtime, package version context,
  and output paths.

Outputs:

- `scripts/rebuild_v9_osd120_diagnostic_package.py`
- `v9/multispecies/reports/interaction_rebuild_gate/rebuild_gate_summary.csv`
- `v9/multispecies/reports/interaction_rebuild_gate/rebuild_gate_steps.csv`
- `v9/multispecies/reports/interaction_rebuild_gate/rebuild_gate_environment.csv`
- `v9/multispecies/reports/OSD120_INTERACTION_REBUILD_GATE_REVIEW.md`

Result:

- The preflight gate reports `ready_existing_outputs_present`.
- It covers 8 packaging steps from V9-MULTI-021 through V9-MULTI-028.
- It hashes 40 packaging outputs and records Python, NumPy, scikit-learn,
  SciPy, and pandas versions.
- It keeps model-grid reruns, full OSDR payload mirroring, and frozen
  benchmark-release language out of scope.

### V9-MULTI-030: source release target and public metadata package

Status: complete as of 2026-05-26.

Goal:

Turn the OSD-120 diagnostic alpha surface into a clearer public metadata
package path without overstating release readiness.

Likely work:

- Start from `interaction_public_alpha_card/`,
  `interaction_payload_freeze_manifest/`, `interaction_rebuild_gate/`, and the
  diagnostic artifact manifest.
- Draft a source release-target decision table separating local diagnostic
  alpha, source package, and future frozen benchmark release.
- Add a machine-readable metadata skeleton that can later become RO-Crate or
  Data Package output.
- Keep full OSDR payload mirror, LOMO/cross-study, biomarker, and operational
  plant-growth claims explicitly disallowed.

Done when:

- The next artifact states what can be public now, what remains local-only or
  draft-only, and what metadata fields are still blocking a citable frozen
  release.

Outputs:

- `scripts/build_v9_osd120_public_metadata_package.py`
- `v9/multispecies/reports/interaction_public_metadata_package/public_metadata_summary.csv`
- `v9/multispecies/reports/interaction_public_metadata_package/source_release_target_decision.csv`
- `v9/multispecies/reports/interaction_public_metadata_package/public_metadata_field_table.csv`
- `v9/multispecies/reports/interaction_public_metadata_package/public_metadata_external_references.csv`
- `v9/multispecies/reports/interaction_public_metadata_package/public_metadata_skeleton.json`
- `v9/multispecies/reports/OSD120_INTERACTION_PUBLIC_METADATA_PACKAGE_REVIEW.md`

Result:

- `diagnostic_alpha_metadata_draft` is the only public-now target and remains
  draft-limited.
- `source_specific_public_alpha_package`, `full_osdr_payload_mirror_release`,
  and `frozen_v9_benchmark_release` are not-public-now targets.
- The metadata field table has 20 fields: 12 present, 3 partial, and 5
  placeholders.
- The skeleton records DataCite 4.7, RO-Crate 1.2, Hugging Face dataset-card,
  and OSDR citation/credit anchors without claiming DOI, license, or frozen
  benchmark release status.

### V9-MULTI-031: RO-Crate and citation freeze scaffold

Status: complete as of 2026-05-26.

Goal:

Promote the V9-MULTI-030 metadata skeleton into a stricter export scaffold that
can later become a citable archive package once placeholder fields are resolved.

Likely work:

- Start from `interaction_public_metadata_package/public_metadata_skeleton.json`
  and `public_metadata_field_table.csv`.
- Emit a draft `ro-crate-metadata.json`-style artifact or validation table
  without pretending the placeholder DOI/creator/license fields are solved.
- Add a citation freeze checklist for OSDR source citation, SpaceBio-Bench
  package citation, version string, rights/license, and related identifiers.
- Keep the current diagnostic-alpha claim boundary intact.

Done when:

- There is one machine-readable export scaffold and one citation-freeze
  checklist showing exactly which fields block an archival release.

Outputs:

- `scripts/build_v9_osd120_ro_crate_citation_scaffold.py`
- `v9/multispecies/reports/interaction_ro_crate_citation_scaffold/ro_crate_export_summary.csv`
- `v9/multispecies/reports/interaction_ro_crate_citation_scaffold/ro_crate_validation_table.csv`
- `v9/multispecies/reports/interaction_ro_crate_citation_scaffold/citation_freeze_checklist.csv`
- `v9/multispecies/reports/interaction_ro_crate_citation_scaffold/ro-crate-metadata.draft.json`
- `v9/multispecies/reports/interaction_ro_crate_citation_scaffold/datapackage.draft.json`
- `v9/multispecies/reports/OSD120_INTERACTION_RO_CRATE_CITATION_SCAFFOLD_REVIEW.md`

Result:

- The scaffold status is
  `draft_scaffold_ready_archive_blocked_by_citation_placeholders`.
- RO-Crate graph has 30 entities, including 26 file Data Entities.
- Data Package descriptor has 26 resources.
- Validation has 13 checks: 9 pass, 1 needs review, and 3 archive blockers.
- Citation checklist has 11 items: 5 pass, 2 needs review, and 4 blockers.
- The next blockers are version/archive identifier, creator/contributor list,
  and license/rights decisions.

### V9-MULTI-032: archive identifier and license decision gate

Status: complete as of 2026-05-26.

Goal:

Resolve or explicitly defer the archive identifier, release version,
creator/contributor, and license/rights decisions that keep the RO-Crate and
Data Package scaffold from becoming citable.

Likely work:

- Start from `interaction_ro_crate_citation_scaffold/citation_freeze_checklist.csv`
  and `ro_crate_validation_table.csv`.
- Draft an archive decision matrix for DOI/Zenodo/SWHID/RAiD/no-archive options.
- Draft a rights/license decision table separating upstream OSDR data credit,
  local metadata package rights, and generated diagnostic tables.
- Keep unresolved fields as blockers unless the user explicitly supplies a
  release owner/creator/license decision.

Done when:

- The archive path, license path, and creator/contributor path are each marked
  as selected, deferred, or blocked with clear evidence and next action.

Outputs:

- `scripts/build_v9_osd120_archive_decision_gate.py`
- `v9/multispecies/reports/interaction_archive_decision_gate/archive_decision_summary.csv`
- `v9/multispecies/reports/interaction_archive_decision_gate/archive_identifier_option_matrix.csv`
- `v9/multispecies/reports/interaction_archive_decision_gate/license_rights_decision_table.csv`
- `v9/multispecies/reports/interaction_archive_decision_gate/creator_contributor_decision_table.csv`
- `v9/multispecies/reports/interaction_archive_decision_gate/archive_decision_external_references.csv`
- `v9/multispecies/reports/OSD120_INTERACTION_ARCHIVE_IDENTIFIER_LICENSE_DECISION_REVIEW.md`

Result:

- The current diagnostic draft selects `current_no_archive_diagnostic_draft`
  and `osdr_source_citation_only`, but does not mint a DOI or local archive
  identifier.
- Zenodo/GitHub DOI, CITATION.cff, and SWHID paths are deferred until release
  owner, version, creator/contributor, and license decisions are supplied.
- Full OSDR payload mirror archiving is blocked because the current package
  intentionally freezes only diagnostic-required payloads.
- Local code license, generated metadata table license, payload-mirror rights,
  local package creators, local package contributors, publisher/maintainer, and
  release version remain blockers for a citable archive.

### V9-MULTI-033: release owner supplied citation metadata fill

Status: complete as of 2026-05-26.

Goal:

Fill only the citation metadata that a release owner explicitly supplies, while
keeping every unsupplied DOI, creator/contributor, version, publisher, and
license field blocked rather than invented.

Likely work:

- Start from
  `interaction_archive_decision_gate/archive_decision_summary.csv`,
  `archive_identifier_option_matrix.csv`,
  `license_rights_decision_table.csv`, and
  `creator_contributor_decision_table.csv`.
- Add a release-owner metadata intake template for archive route, version tag,
  package title, creators, contributors, publisher/maintainer, local license,
  generated metadata license, and OSDR source citation text.
- Add validation that separates `supplied`, `not_supplied`, `blocked`, and
  `not_applicable_current_draft` states without making default legal or
  authorship decisions.
- If no owner-supplied fields are available, emit a no-op fill report that
  preserves the V9-MULTI-032 no-archive diagnostic boundary.

Done when:

- The repository has a script-backed citation metadata intake/fill scaffold
  that can safely update the RO-Crate/Data Package descriptors only when
  explicit release-owner values are present.

Outputs:

- `scripts/build_v9_osd120_citation_metadata_fill.py`
- `v9/multispecies/reports/interaction_citation_metadata_fill/citation_metadata_fill_summary.csv`
- `v9/multispecies/reports/interaction_citation_metadata_fill/release_owner_metadata_intake_template.csv`
- `v9/multispecies/reports/interaction_citation_metadata_fill/citation_metadata_fill_status.csv`
- `v9/multispecies/reports/interaction_citation_metadata_fill/citation_descriptor_patch_preview.json`
- `v9/multispecies/reports/OSD120_INTERACTION_CITATION_METADATA_FILL_REVIEW.md`

Result:

- No release-owner metadata file was supplied, so the fill status is
  `no_owner_metadata_supplied_no_descriptor_changes`.
- The scaffold records 16 owner-fill fields: 4 retained current-draft values
  and 12 release blockers.
- The descriptor preview explicitly reports `mutates_ro_crate=false` and
  `mutates_datapackage=false`.
- If an owner metadata CSV/JSON is supplied later, the CLI emits a patch preview
  only; descriptor mutation remains gated behind a later application block.

### V9-MULTI-034: owner metadata application guard or archive release deferral

Status: complete as of 2026-05-26.

Goal:

Decide whether release-owner metadata has arrived and either apply it through a
strict descriptor mutation guard or formally defer archive release while keeping
the diagnostic metadata package usable.

Likely work:

- Start from
  `interaction_citation_metadata_fill/citation_metadata_fill_summary.csv`,
  `release_owner_metadata_intake_template.csv`,
  `citation_metadata_fill_status.csv`, and
  `citation_descriptor_patch_preview.json`.
- If no owner metadata is supplied, emit a release-deferral memo and keep the
  public output as diagnostic metadata only.
- If owner metadata is supplied, validate every supplied field before mutating
  RO-Crate/Data Package/CITATION.cff draft outputs.
- Keep legal, authorship, publisher, archive route, and license decisions
  non-inferred.

Done when:

- The next artifact either records a clean deferral path or applies only
  validated owner-supplied metadata into regenerated descriptors with tests.

Outputs:

- `scripts/build_v9_osd120_archive_release_deferral_guard.py`
- `v9/multispecies/reports/interaction_archive_release_deferral_guard/archive_release_deferral_summary.csv`
- `v9/multispecies/reports/interaction_archive_release_deferral_guard/owner_metadata_application_guard.csv`
- `v9/multispecies/reports/interaction_archive_release_deferral_guard/archive_release_deferral_actions.csv`
- `v9/multispecies/reports/interaction_archive_release_deferral_guard/descriptor_mutation_guard.json`
- `v9/multispecies/reports/OSD120_INTERACTION_ARCHIVE_RELEASE_DEFERRAL_GUARD_REVIEW.md`

Result:

- Archive release is deferred with status
  `archive_release_deferred_no_owner_metadata`.
- The guard records 11 checks: 9 blockers and 2 pass checks.
- The 2 pass checks are descriptor mutation prevention and current diagnostic
  metadata surface retention.
- The action table records 6 release-owner actions plus 3 selected
  maintainer-side deferral/guard actions.
- Descriptor mutation is not allowed; RO-Crate, Data Package, and CITATION.cff
  generation remain blocked until owner metadata is supplied and validated.

### V9-MULTI-035: diagnostic metadata release note or owner metadata intake retry

Status: complete as of 2026-05-26.

Goal:

Use the V9-MULTI-034 deferral guard to choose the next public-facing path:
either write a diagnostic metadata release note that avoids archive-release
claims, or prepare a retry protocol for owner-supplied citation metadata intake.

Likely work:

- Start from
  `interaction_archive_release_deferral_guard/archive_release_deferral_summary.csv`,
  `owner_metadata_application_guard.csv`,
  `archive_release_deferral_actions.csv`, and
  `descriptor_mutation_guard.json`.
- Draft a diagnostic metadata release note that says what is inspectable now and
  what is explicitly not released as an archive.
- Keep DOI, CITATION.cff, release tag, creator/publisher, license, and exact
  OSDR source citation out of the public release note unless owner metadata is
  supplied.
- Optionally emit a small owner-metadata retry checklist for a later release
  owner review.

Done when:

- The next artifact provides a clean user-facing diagnostic metadata note or a
  clearly bounded owner-intake retry packet without changing archive metadata.

Outputs:

- `scripts/build_v9_osd120_diagnostic_metadata_release_note.py`
- `v9/multispecies/reports/interaction_diagnostic_metadata_release_note/diagnostic_metadata_release_summary.csv`
- `v9/multispecies/reports/interaction_diagnostic_metadata_release_note/diagnostic_metadata_release_note_sections.csv`
- `v9/multispecies/reports/interaction_diagnostic_metadata_release_note/diagnostic_metadata_public_claims.csv`
- `v9/multispecies/reports/interaction_diagnostic_metadata_release_note/owner_metadata_retry_checklist.csv`
- `v9/multispecies/reports/interaction_diagnostic_metadata_release_note/OSD120_DIAGNOSTIC_METADATA_NOTE.md`
- `v9/multispecies/reports/OSD120_INTERACTION_DIAGNOSTIC_METADATA_RELEASE_NOTE_REVIEW.md`

Result:

- The closeout note status is
  `diagnostic_metadata_note_ready_no_archive_release`.
- The current OSD-120 public surface is `diagnostic_metadata_only`.
- The claim table records 3 allowed diagnostic metadata claims and 3 prohibited
  release claims.
- The owner retry checklist records 6 missing owner-supplied metadata items:
  owner metadata file, archive route/identifier, version/date/year,
  creator/contributor/publisher, license scope, and exact OSDR study citation.
- No DOI, CITATION.cff, release tag, creator/publisher, license, archive
  identifier, or leaderboard promotion is claimed.
- The next required block is
  `V9-REFOCUS-001: post-OSD-120 recenter decision`.

Drift guard:

- This is a closeout block, not the start of another OSD-120 release-metadata
  sequence.
- If no owner metadata is supplied, do not create another archive/citation
  guard. Close the branch and return to public bulk alpha or single-cell
  flagship planning.

### V9-REFOCUS-001: post-OSD-120 recenter decision

Status: complete as of 2026-05-26.

Goal:

Choose the next v9 lane after the OSD-120 metadata closeout so the project
returns to its original benchmark-platform sequence.

Likely work:

- Start from `docs/V9_PURPOSE_DRIFT_AUDIT_2026_05_26.md`.
- Compare public bulk alpha readiness versus first single-cell flagship
  readiness.
- Pick one next 30-60 minute implementation block with files, tests, and
  expected outputs.

Done when:

- The backlog names either a public bulk alpha recenter block or a single-cell
  flagship inventory block as the next active lane.

Outputs:

- `scripts/build_v9_recenter_decision.py`
- `v9/reports/recenter_decision/recenter_decision_summary.csv`
- `v9/reports/recenter_decision/recenter_candidate_matrix.csv`
- `v9/reports/recenter_decision/recenter_next_block_actions.csv`
- `v9/reports/recenter_decision/single_cell_asset_paths.json`
- `docs/V9_REFOCUS_001_POST_OSD120_RECENTER_DECISION.md`

Result:

- Selected next lane: `public_bulk_alpha`.
- Selected next block:
  `V9-BULK-ALPHA-001: public bulk alpha freeze-gap matrix`.
- Public bulk readiness snapshot: 8 task manifests, 33 folds, 22 source rows,
  22/22 API-ok source rows, 22/22 checksum-manifest-parsed rows, 24 evaluated
  baseline rows, and 0/22 `freeze_ready=true` source rows.
- Single-cell readiness snapshot: 54 legacy RRRM/single-cell-related paths,
  `genelab_sc` metric profile present, but 0 local h5ad/loom/mtx files found
  in the repo scan and 0 v9 `sc_spaceflight` task manifests.
- Decision boundary:
  `recenter_decision_only_no_new_benchmark_or_release_claim`.

### V9-BULK-ALPHA-001: public bulk alpha freeze-gap matrix

Status: complete as of 2026-05-27.

Goal:

Create a machine-readable pass/blocker matrix for the public bulk alpha
scaffold before returning to single-cell flagship implementation.

Likely work:

- Start from `v9/source_inventory.csv`, `v9/source_checksum_audit.csv`,
  `v9/task_manifest_index.csv`, `v9/reports/bulk_lomo_baseline_summary.csv`,
  `v9/datapackage.draft.json`, and `docs/v9_hf_dataset_card.md`.
- Separate API/checksum-manifest evidence from locally verified payload hash
  evidence.
- Decide whether a metadata-only alpha snapshot is acceptable before full
  payload mirroring.
- Keep public bulk core separate from organoid and multispecies draft tracks.

Done when:

- Every public bulk alpha blocker has a status, evidence path, owner action,
  and allowed/disallowed claim boundary.

Outputs:

- `scripts/build_v9_public_bulk_alpha_gap_matrix.py`
- `v9/reports/public_bulk_alpha_gap_matrix/public_bulk_alpha_gap_summary.csv`
- `v9/reports/public_bulk_alpha_gap_matrix/public_bulk_alpha_gap_matrix.csv`
- `v9/reports/public_bulk_alpha_gap_matrix/payload_hash_boundary.csv`
- `v9/reports/public_bulk_alpha_gap_matrix/public_bulk_alpha_claim_boundary.csv`
- `v9/reports/public_bulk_alpha_gap_matrix/package_update_plan.csv`
- `docs/V9_PUBLIC_BULK_ALPHA_FREEZE_GAP_MATRIX_REVIEW.md`

Result:

- The public bulk scaffold remains the active lane, but frozen release language
  is blocked.
- Gap matrix: 6 pass rows, 2 blocker rows, and 2 needs-update rows.
- Payload boundary: 22/22 sources have checksum-manifest evidence, but 0/22
  sources are locally payload-hash verified with `freeze_ready=true`.
- Claim boundary: `public_bulk_alpha_gap_matrix_no_release_approval`.
- Next required block:
  `V9-BULK-ALPHA-002: metadata-only alpha snapshot decision`.

### V9-BULK-ALPHA-002: metadata-only alpha snapshot decision

Status: active next.

Goal:

Decide whether public bulk v9 can publish a metadata-only alpha snapshot with
explicit payload-hash blockers, or whether payload mirroring and local hash
verification must precede any alpha wording.

Likely work:

- Start from
  `v9/reports/public_bulk_alpha_gap_matrix/public_bulk_alpha_gap_summary.csv`,
  `public_bulk_alpha_gap_matrix.csv`, `payload_hash_boundary.csv`,
  `public_bulk_alpha_claim_boundary.csv`, and
  `docs/V9_PUBLIC_BULK_ALPHA_FREEZE_GAP_MATRIX_REVIEW.md`.
- Compare two paths: metadata-only alpha snapshot versus payload-mirror-first.
- Update dataset-card/package language only if the chosen wording is explicit
  about checksum-manifest evidence versus local payload hash verification.

Done when:

- The backlog and user-facing draft language name either metadata-only alpha
  with blockers or payload-mirror-first as the next path.

## Manuscript and community layer

### V9-MANU-001: Platform paper outline

Output:

- `docs/V9_SPACEBIO_BENCH_MANUSCRIPT_OUTLINE.md`

Required sections:

- motivation
- benchmark gap
- task families
- metric profiles
- baseline results
- provenance/release model
- limitations
- community submission path

### V9-MANU-002: Claim-to-artifact map

Output:

- `v9/claim_to_artifact_map.csv`

Done when:

- Every claim maps to task manifest, source records, code, output, and test or
  validation evidence.

## Operating rules

- Do not expand to single-cell or radiation implementation until platform spine
  can evaluate at least one bulk task end to end.
- Do not add foundation-model adapters before simple baselines.
- Do not make public release claims while `checksum_status` remains
  `legacy_task_source_unfrozen`.
- Do not treat gated human data as required for public v9.
- Do not use intervention language that implies clinical or crew-health
  recommendation.
- Do not describe Mars outputs as point predictions.

## Session handoff template

At the end of a long v9 session, update this block or the relevant task status:

```text
Last task worked:
Files changed:
Commands run:
Tests passed:
Known blockers:
Recommended next action:
```

Latest handoff, 2026-05-21:

```text
Last task worked:
V9-ORG-003 and V9-ORG-004 after initial V9-ORG/V9-MULTI scaffolding.
Files changed:
spacebio_bench source/schema/registry/legacy helpers; extension inventory CLI;
human organoid task-manifest CLI; generated task/source indexes; draft
human_organoid and multispecies inventories; human organoid checksum audit; v9
docs/backlog/README.
Commands run:
scripts/export_v9_task_manifests.py
scripts/build_v9_task_index.py
scripts/build_v9_source_inventory.py
scripts/build_v9_extension_source_inventories.py
scripts/build_v9_human_organoid_task_manifest.py
scripts/audit_v9_source_checksums.py --source-inventory v9/human_organoid/source_inventory.draft.csv --csv v9/human_organoid/source_checksum_audit.draft.csv --json v9/human_organoid/source_checksum_audit.draft.json
scripts/audit_v9_sample_tables.py --source-inventory v9/human_organoid/source_inventory.draft.csv --csv v9/human_organoid/sample_table_audit.draft.csv --json v9/human_organoid/sample_table_audit.draft.json
scripts/build_v9_human_organoid_sample_factors.py --source-inventory v9/human_organoid/source_inventory.draft.csv --csv v9/human_organoid/sample_factors.draft.csv --json v9/human_organoid/sample_factors.draft.json
scripts/build_v9_datapackage_draft.py
Tests passed:
/usr/local/bin/python3 -m unittest tests/test_v9_spacebio_bench.py
/usr/local/bin/python3 -m unittest discover -s tests
Known blockers:
Human organoid API/checksum-manifest audit is done, but payload files have not
been downloaded and hashed. Sample tables and condition factors are parsed, but
payload download/hash, baseline-ready expression matrices, and baselines are not
done for extension tracks.
Recommended next action:
Audit/download the processed normalized-count matrices for OSD-863/OSD-871 and
build the first human-organoid baseline loader over the 42 parsed samples.
```

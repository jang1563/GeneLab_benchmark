# v9 Operating Backlog

Status: living backlog for SpaceBio-Bench v9
Last updated: 2026-05-21

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
- Draft multispecies source inventory generated at
  `v9/multispecies/source_inventory.draft.csv`.
- Multispecies feature namespace strategy written at
  `docs/V9_MULTISPECIES_FEATURE_STRATEGY.md`.
- v9 unit tests integrated into full repository tests.

Current test baseline:

```bash
python -m unittest discover -s tests
```

Expected status after DE/signature reference audit:

- 110 tests passing.

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

Continue with `V9-ORG-014: Human organoid frozen DE contrast extraction and
signature metric contract`.
`V9-THEN-005: RO-Crate export design` remains open for the release-packaging
lane, but the active implementation lane is now the organoid extension because
the source, sample, matrix, GEO donor metadata, donor diagnostics, pilot
baseline, robustness, DE reference audit, and local diagnostic artifacts are
aligned.

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
  now reference-backed, but non-primary until a frozen contrast subset, sign
  orientation, and signature-output schema are added.

### V9-ORG-014: Human organoid frozen DE contrast extraction and signature metric contract

Status: active next task.

Goal:

Turn the audited public OSDR DE references into a compact, checksummed reference
artifact and an evaluator contract for gene-level response signatures.

Likely work:

- Download or stream the canonical non-rRNArm OSDR DE tables for OSD-863 and
  OSD-871.
- Extract the four direct matched Ground Control versus Space Flight contrasts
  per source.
- Normalize log2 fold-change orientation to `LEO_or_ISS - Ground`.
- Write a compact reference table and checksum manifest.
- Define a `response_signature.csv` submission artifact and skip-aware
  evaluator behavior for submissions that do not provide gene-level signatures.

Done when:

- A small derived DE reference table exists with source, contrast, gene id,
  log2 fold-change, adjusted p-value, and orientation metadata.
- The evaluator can report DE/signature metrics as skipped with a precise reason
  unless a valid response-signature artifact is supplied.

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

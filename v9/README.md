# SpaceBio-Bench v9 Planning Artifacts

Status: planning/alpha scaffold. This directory does not claim v9 benchmark
results yet.

## Current contents

- `task_manifests/*.json`: v9 task manifests exported from the existing
  `tasks/*/task_info.json` bulk LOMO tasks.
- `task_manifest_index.csv` and `task_manifest_index.json`: compact registry
  summaries generated from the task manifests.
- `task_data_index.csv` and `task_data_index.json`: fold-level data path/count
  summaries for the public bulk LOMO scaffold.
- `reports/nearest_centroid/`: first runnable v9-alpha baseline report across
  all eight generated bulk LOMO manifests.
- `reports/sklearn_baselines/`: PCA-LR and L2 logistic-regression reports
  across all eight generated bulk LOMO manifests.
- `reports/bulk_lomo_baseline_summary.csv` and `.json`: normalized
  cross-baseline summary for nearest-centroid plus sklearn baselines.
- `source_inventory.csv` and `source_inventory.json`: source-level inventory
  deduplicated from all generated bulk LOMO task manifests.
- `source_checksum_audit.csv` and `source_checksum_audit.json`: OSDR API
  file-list and checksum-manifest evidence for the 22 public bulk source rows.
- `datapackage.draft.json`: draft Frictionless Data Package descriptor for
  v9 public bulk metadata, provenance evidence, and baseline outputs.
- `human_organoid/source_inventory.draft.csv` and `.json`: draft source
  inventory for public human neural organoid RNA-seq extension candidates
  OSD-863 and OSD-871.
- `human_organoid/task_manifests/*.json`: draft human organoid task manifest.
- `human_organoid/task_manifest_index.draft.csv` and `.json`: draft index for
  human organoid extension task manifests.
- `human_organoid/source_checksum_audit.draft.csv` and `.json`: OSDR API
  file-list and checksum-manifest evidence for OSD-863 and OSD-871.
- `human_organoid/sample_table_audit.draft.csv` and `.json`: OSDR sample-table
  discovery and row-count evidence for OSD-863 and OSD-871.
- `human_organoid/sample_factors.draft.csv` and `.json`: parsed sample-level
  disease, Ground/LEO, microglia, organoid-type, GEO Subject, and iPSC-line
  factors for the 42 public GSE259421 samples.
- `human_organoid/geo_sample_metadata.draft.csv` and `.json`: GEO GSE259421
  series-matrix-derived sample titles, Subject identifiers, iPSC-line ids,
  BioSample/SRA accessions, cell type, and treatment metadata.
- `human_organoid/expression_matrix_audit.draft.csv` and `.json`: OSDR
  normalized-count matrix discovery, download, SHA-256, row-count, and
  sample-column alignment evidence for OSD-863 and OSD-871.
- `human_organoid/signature_reference_audit.draft.csv` and `.json`: OSDR
  DE/signature reference audit showing public differential-expression tables,
  contrast-definition files, direct/reversed spaceflight contrast counts, and
  the current metric policy.
- `human_organoid/matrices/`: downloaded human organoid normalized-count
  matrices referenced by the draft expression-matrix audit.
- `human_organoid/reports/nearest_centroid/`: first draft human organoid pilot
  baseline outputs, including predictions, metrics, run manifest, and summary.
- `human_organoid/reports/sensitivity/`: draft sensitivity grid for the human
  organoid pilot baseline across transform, scaling, and variable-gene settings.
- `human_organoid/reports/ORGANOID_BASELINE_ROBUSTNESS_REVIEW.md`: conservative
  robustness and confounding review for the draft organoid pilot baseline.
- `human_organoid/sample_scale_diagnostics.draft.csv` and `.json`: per-sample
  expression scale diagnostics over common human gene features.
- `human_organoid/group_scale_diagnostics.draft.csv` and `.json`: grouped
  sample-total and sparsity summaries by source, label, organoid type,
  microglia condition, disease context, donor/Subject id, and iPSC-line id.
- `human_organoid/donor_metadata_audit.draft.csv` and `.json`: donor/iPSC-line
  metadata availability audit showing OSDR SampleTables lack donor fields but
  GEO series metadata recovers Subject/iPSC-line ids for all 42 samples.
- `human_organoid/reports/ORGANOID_DONOR_LIBRARY_SIZE_AUDIT.md`: interpretation
  note for donor metadata gaps and library-size diagnostics.
- `human_organoid/reports/donor_diagnostics/`: diagnostic-only GEO
  donor/Subject holdout baseline outputs.
- `human_organoid/reports/ORGANOID_DONOR_AWARE_SPLIT_DECISION.md`: decision
  note keeping donor holdouts diagnostic-only rather than promoting them to the
  default benchmark split.
- `human_organoid/reports/ORGANOID_SIGNATURE_METRIC_REFERENCE_DECISION.md`:
  decision note keeping DE/signature metrics reference-backed but non-primary
  until a frozen contrast/signature contract is added.
- `multispecies/source_inventory.draft.csv` and `.json`: draft source
  inventory for non-mouse pilot sources OSD-207, OSD-37, and OSD-120.
- `OPERATING_BACKLOG.md`: living long-run backlog for continuing v9 work across
  sessions.

## Planning docs

- `docs/V9_LONG_HORIZON_EXECUTION_PLAN.md`: long-term execution plan,
  workstreams, phases, decision gates, and risks.
- `docs/V9_PUBLIC_BULK_PACKAGE_DESIGN.md`: public bulk package boundary and
  draft Data Package design.
- `docs/v9_hf_dataset_card.md`: draft Hugging Face-style dataset card for the
  public bulk scaffold.
- `docs/V9_DESIGN_OPTIONS.md`: v9 option matrix and initial 30/60/90-day
  design path.
- `docs/V9_EXTERNAL_DEEP_RESEARCH_2026_05_21.md`: external ecosystem research.
- `docs/V9_SOURCE_AND_COMPETITOR_MATRIX.md`: source and competitor matrix.
- `docs/V9_ORGANOID_AND_SPECIES_EXTENSION_REVIEW_2026_05_21.md`: evidence
  review for human organoid and non-mouse species extension candidates.
- `docs/V9_MULTISPECIES_FEATURE_STRATEGY.md`: feature-namespace strategy for
  species-native and cross-species bridge tasks.

## Regeneration

```bash
python scripts/export_v9_task_manifests.py
python scripts/build_v9_task_index.py
python scripts/build_v9_task_data_index.py
python scripts/run_v9_nearest_centroid.py \
  --output-dir v9/reports/nearest_centroid
python scripts/run_v9_sklearn_baselines.py \
  --output-dir v9/reports/sklearn_baselines
python scripts/build_v9_baseline_summary.py
python scripts/build_v9_source_inventory.py
python scripts/build_v9_extension_source_inventories.py
python scripts/audit_v9_source_checksums.py \
  --source-inventory v9/human_organoid/source_inventory.draft.csv \
  --csv v9/human_organoid/source_checksum_audit.draft.csv \
  --json v9/human_organoid/source_checksum_audit.draft.json
python scripts/audit_v9_sample_tables.py \
  --source-inventory v9/human_organoid/source_inventory.draft.csv \
  --csv v9/human_organoid/sample_table_audit.draft.csv \
  --json v9/human_organoid/sample_table_audit.draft.json
python scripts/build_v9_human_organoid_sample_factors.py \
  --source-inventory v9/human_organoid/source_inventory.draft.csv \
  --csv v9/human_organoid/sample_factors.draft.csv \
  --json v9/human_organoid/sample_factors.draft.json
python scripts/build_v9_human_organoid_geo_metadata.py \
  --series-matrix-gz path/to/GSE259421_series_matrix.txt.gz
python scripts/audit_v9_expression_matrices.py \
  --source-inventory v9/human_organoid/source_inventory.draft.csv \
  --sample-factors v9/human_organoid/sample_factors.draft.csv \
  --matrix-dir v9/human_organoid/matrices \
  --csv v9/human_organoid/expression_matrix_audit.draft.csv \
  --json v9/human_organoid/expression_matrix_audit.draft.json
python scripts/audit_v9_human_organoid_signature_references.py
python scripts/build_v9_human_organoid_task_manifest.py
python scripts/run_v9_human_organoid_baseline.py
python scripts/run_v9_human_organoid_sensitivity.py
python scripts/audit_v9_human_organoid_diagnostics.py
python scripts/run_v9_human_organoid_donor_diagnostics.py
python scripts/audit_v9_source_checksums.py
python scripts/build_v9_datapackage_draft.py
python scripts/evaluate_v9_submission.py \
  --task-manifest v9/task_manifests/A2_gastrocnemius_bulk_lomo.json \
  --submission path/to/predictions.csv \
  --output-dir path/to/report
```

The exporter reads legacy task metadata from `tasks/`, attaches public OSDR
source records from `scripts/fetch_osdr.py`, and validates each generated JSON
against the minimal `spacebio_bench` manifest contract.

The index builder reads the generated manifests and writes compact CSV/JSON
summaries for quick inspection and downstream loaders.

The data-index builder validates legacy fold files for each bulk LOMO manifest
and writes fold-level path/count summaries.

The evaluator validates a prediction CSV, computes available `genelab_minimal`
metrics, and writes `metrics.json` plus `run_manifest.json`.

The nearest-centroid runner reads the legacy processed fold matrices through
the v9 loader, writes one `predictions.csv`, `metrics.json`, and
`run_manifest.json` per task, then writes a cross-task summary at
`v9/reports/nearest_centroid/bulk_lomo_summary.csv`.

The sklearn runner currently evaluates PCA-LR and L2 logistic regression using
train-fold-fitted preprocessing only. The cross-baseline summary builder
normalizes nearest-centroid and sklearn outputs into
`v9/reports/bulk_lomo_baseline_summary.csv`.

The source-inventory builder deduplicates `source_records` across all task
manifests and writes one row per OSDR source accession. The current inventory
covers 22 public bulk source rows and now includes organism, material,
model-system, assay-modality, and feature-namespace fields. It keeps
`checksum_status` at
`legacy_task_source_unfrozen` because the generated manifests have not yet been
promoted to a frozen v9 release state.

The extension-inventory builder writes draft inventories for public human
organoid and multispecies candidates. These files are planning artifacts, not
frozen benchmark tasks. They keep organoids, plants, and fly data outside the
mouse bulk LOMO leaderboard until split rules, feature namespaces, checksum
audits, and baselines are separately defined.

The human-organoid task-manifest builder converts the draft OSD-863/OSD-871
source inventory into a draft `human_organoid_spaceflight` manifest and a
separate extension task index. It records a blocked LEO/ISS-versus-ground pilot
design, donor/organoid/microglia blocking requirements, the single-mission
limitation, sample-factor-backed draft folds, and normalized expression matrix
alignment status. It also records that public OSDR DE tables and contrast
definitions are available for future DE/signature metrics, while keeping those
metrics non-primary until a frozen contrast/signature contract is added. It
still treats payload checksum verification as pending and keeps organoid
baseline language draft-only.

The human-organoid checksum audit uses the same OSDR file-list/checksum parser
as the public mouse bulk audit. It currently records API-ok file listings and
parsed raw/processed MD5 manifests for OSD-863 and OSD-871, while keeping
`freeze_ready=false` because payload files have not been downloaded and hashed.

The human-organoid sample-table audit finds and parses one OSDR SampleTable per
source. It currently records 19 sample-table rows for OSD-863 and 23 for
OSD-871, matching the 42 public samples reported for GSE259421. Manual condition
factor mapping is still required before fold generation.

The human-organoid sample-factor builder parses the compact OSDR `condition`
field into `disease_context`, `spaceflight_condition`, `true_label`, and
`microglia_condition`. The GEO metadata builder then augments those rows from
GSE259421 series-matrix metadata with `donor_or_line_id`, `ipsc_line_id`, sample
title, cell type, and treatment. The draft task manifest carries four default
sample-count-backed draft folds: hold out cortical organoids, hold out
dopaminergic organoids, hold out with-microglia samples, and hold out
without-microglia samples. GEO-derived donor/Subject holdouts are recorded as
diagnostic folds only because donor, source, organoid fate, and disease context
are not independently crossed.

The human-organoid expression-matrix audit discovers and downloads the canonical
OSDR `*_rna_seq_Normalized_Counts_GLbulkRNAseq.csv` matrices, excluding rRNArm
and unnormalized files. The current audit records one aligned matrix per source:
OSD-863 has 30,408 feature rows and 19/19 matched sample columns, while OSD-871
has 30,269 feature rows and 23/23 matched sample columns.

The human-organoid nearest-centroid pilot baseline loads both aligned matrices,
uses their 27,986 common human gene features, applies train-fold-only log1p,
variable-gene selection, and z-scoring, then evaluates the four
sample-count-backed draft folds. Current summary output is
`v9/human_organoid/reports/nearest_centroid/human_organoid_baseline_summary.csv`
with `release_status=draft_not_frozen` and
`claim_boundary=pilot_baseline_only_not_leaderboard`.

The human-organoid sensitivity runner evaluates 20 preprocessing variants across
log1p/no-transform, z-score/no-scaling, and five train-fold variable-gene
settings. The robustness review keeps the log1p/z-score/2,000-gene baseline as
the conservative default, despite higher raw-count sensitivity scores, because
raw-scale variants show worse calibration and need library-size/confounding
audits before interpretation.

The human-organoid diagnostics audit records that current OSDR SampleTable files
expose only `condition`, while GEO GSE259421 series metadata recovers
Subject/iPSC-line ids for all 42 public samples. Donor-aware folds remain
diagnostic-only because the design is small and coupled. The audit also shows
label-associated sample-total and sparsity differences, which keeps raw-count
sensitivity variants out of the default baseline.

The donor diagnostic runner evaluates the GEO-derived Subject/iPSC-line holdouts
without changing the default benchmark fold family. Current donor diagnostic
outputs stay under `v9/human_organoid/reports/donor_diagnostics/` with
`claim_boundary=donor_diagnostic_only_not_leaderboard`.

The human-organoid signature-reference audit records that OSDR lists public
differential-expression tables and contrast-definition CSVs for both organoid
sources. Each source has 56 parsed contrast pairs, including four direct
matched Ground Control versus Space Flight contrasts and four reversed matches.
The draft task manifest therefore marks `de_direction_match` and
`signature_rank_correlation` as reference-backed future metrics, but keeps
classification metrics primary until a frozen contrast/signature contract is
implemented.

The checksum-audit script queries the OSDR biological data API file-list
endpoint for each source, records API response hashes, discovers checksum
manifest-like files, downloads the manifest text, parses MD5 entries in the
OSDR formats observed so far, and records match counts against the OSDR file
listing. The current audit covers all 22 public bulk source rows, finds 39
checksum manifests, parses 8,439 checksum entries, and matches 8,275 entries
back to listed OSDR payload names.

The Data Package draft builder writes `v9/datapackage.draft.json`. It describes
11 current resources: metadata indexes, source/provenance tables, baseline
summary tables, eight task manifests, 24 prediction CSVs, 24 metrics JSON files,
and 24 run-manifest JSON files. It deliberately excludes large fold payload
hashing until a dedicated public payload manifest block.

## Interpretation

- `checksum_status=legacy_task_source_unfrozen` means these manifests preserve
  the current legacy source identity but are not yet a v9 frozen release.
- `freeze_ready=false` in `source_checksum_audit.csv` means checksum-manifest
  evidence exists, but payload files have not yet been downloaded and hashed in
  this audit.
- The public benchmark core does not depend on gated human datasets.
- Intervention or Mars-regime claims are intentionally out of scope for these
  bulk LOMO task manifests.

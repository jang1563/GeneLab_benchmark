# SpaceBio-Bench v9 Public Bulk Metadata Catalog

This directory contains the public metadata catalog for the SpaceBio-Bench v9
bulk RNA-seq task surface. It is organized for direct review of task manifests,
fold indexes, public OSDR source coverage, checksum-audit summaries, and
reference baseline outputs.

For the reader-facing card, see
[`docs/v9_hf_dataset_card.md`](../docs/v9_hf_dataset_card.md).

## Included

- `task_manifests/*.json`: eight public bulk LOMO task manifests.
- `task_manifest_index.csv` and `.json`: task registry summary.
- `task_data_index.csv` and `.json`: fold-level row-count and path registry.
- `source_inventory.csv` and `.json`: public OSDR source rows.
- `source_checksum_audit.csv` and `.json`: OSDR API and checksum-manifest
  summary rows.
- `datapackage.draft.json`: Frictionless Data Package descriptor for the
  metadata catalog.
- `reports/bulk_lomo_baseline_summary.csv` and `.json`: normalized summary of
  reference baselines.
- `reports/nearest_centroid/bulk_lomo_summary.csv` and `.json`:
  nearest-centroid baseline summary.
- `reports/sklearn_baselines/bulk_lomo_summary.csv` and `.json`: PCA-LR and
  L2 logistic-regression baseline summary.
- Per-task baseline `predictions.csv`, `metrics.json`, and `run_manifest.json`
  files referenced by `datapackage.draft.json`.

## Scope Notes

- The v9 metadata catalog is separate from archived fold-matrix payload
  bundles.
- Controlled-access human sequence data is not part of this public bulk catalog.
- Clinical, crew-health, intervention, countermeasure, operational-readiness,
  and state-of-the-art leaderboard claims are outside the catalog scope.

## Interpretation

The files here support review of task definitions, source coverage, fold
indexes, checksum-audit summaries, and baseline workflow outputs. They are best
read as a public metadata catalog for method-development and release review.

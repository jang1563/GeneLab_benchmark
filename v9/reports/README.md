# v9 Public Bulk Reports

This directory contains compact report tables for the v9 public bulk metadata
catalog.

- `bulk_lomo_baseline_summary.csv` and `.json`: cross-baseline summary across
  eight public bulk LOMO task manifests and three baseline families.
- `nearest_centroid/bulk_lomo_summary.csv` and `.json`: nearest-centroid
  baseline summary.
- `sklearn_baselines/bulk_lomo_summary.csv` and `.json`: PCA-LR and L2
  logistic-regression baseline summary.
- Per-task baseline `predictions.csv`, `metrics.json`, and `run_manifest.json`
  files for the baseline families above.

The baseline rows are reproducible reference anchors for the public bulk task
catalog. Use per-task rows together with pooled means when comparing methods.

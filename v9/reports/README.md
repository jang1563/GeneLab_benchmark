# v9 Public Bulk Metadata-Alpha Reports

Status: scaffold evidence, not a frozen public benchmark release.

This directory contains compact summary reports for the v9 public bulk
metadata-alpha subset:

- `bulk_lomo_baseline_summary.csv` and `.json`: cross-baseline summary across
  eight public bulk LOMO task manifests and three scaffold baselines.
- `nearest_centroid/bulk_lomo_summary.csv` and `.json`: nearest-centroid
  scaffold baseline summary.
- `sklearn_baselines/bulk_lomo_summary.csv` and `.json`: PCA-LR and L2
  logistic-regression scaffold baseline summary.
- Per-task scaffold baseline `predictions.csv`, `metrics.json`, and
  `run_manifest.json` files for the baseline families above.
- `public_bulk_alpha_gap_matrix/`: payload-boundary and package-readiness gap
  tables.
- `public_bulk_alpha_snapshot_decision/`: selected metadata-alpha release path,
  allowed language, blocked language, and next actions.

The baseline rows validate evaluation plumbing and provide comparison anchors.
They are not tuned leaderboard endpoints and should not be described as model
superiority evidence.

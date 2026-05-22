# v9 Reports

Status: alpha scaffold outputs, not a frozen public benchmark release.

## Cross-Baseline Bulk LOMO Summary

Regeneration commands:

```bash
python scripts/run_v9_nearest_centroid.py --output-dir v9/reports/nearest_centroid
python scripts/run_v9_sklearn_baselines.py --output-dir v9/reports/sklearn_baselines
python scripts/build_v9_baseline_summary.py
```

Outputs:

- `bulk_lomo_baseline_summary.csv`
- `bulk_lomo_baseline_summary.json`

Current summary scope:

- 8 generated bulk LOMO task manifests.
- 3 simple baselines: nearest-centroid, PCA-LR, and L2 logistic regression.
- 24 evaluated task/baseline rows.

Best AUROC per task in the current scaffold:

| Task | Best Baseline | AUROC | Balanced Accuracy |
|---|---|---:|---:|
| A1_liver_bulk_lomo | logistic_regression_l2 | 0.6468 | 0.5511 |
| A1_liver_bulk_lomo_combat | logistic_regression_l2 | 0.6236 | 0.6027 |
| A1_liver_bulk_lomo_iss_only | logistic_regression_l2 | 0.6559 | 0.5405 |
| A2_gastrocnemius_bulk_lomo | logistic_regression_l2 | 0.7773 | 0.5000 |
| A3_kidney_bulk_lomo | logistic_regression_l2 | 0.6722 | 0.6294 |
| A4_thymus_bulk_lomo | nearest_centroid | 0.6578 | 0.5557 |
| A5_skin_bulk_lomo | logistic_regression_l2 | 0.7847 | 0.7326 |
| A6_eye_bulk_lomo | pca_logistic_regression | 0.7206 | 0.5471 |

The `flight_probability` column in all current baselines is a model score used
for ranking/AUROC and calibration diagnostics. The scaffold should not describe
these values as validated posterior probabilities.

## Nearest-Centroid Bulk LOMO Baseline

Regeneration command:

```bash
python scripts/run_v9_nearest_centroid.py --output-dir v9/reports/nearest_centroid
```

Outputs:

- `nearest_centroid/<task_id>/predictions.csv`
- `nearest_centroid/<task_id>/metrics.json`
- `nearest_centroid/<task_id>/run_manifest.json`
- `nearest_centroid/bulk_lomo_summary.csv`
- `nearest_centroid/bulk_lomo_summary.json`

The `flight_probability` column is a deterministic distance-ratio score derived
from Euclidean distances to the Ground and Flight centroids. It is useful for
ranking/AUROC checks, but it should not be interpreted as a calibrated posterior
probability.

Current all-task summary:

| Task | Macro-F1 | Balanced Accuracy | AUROC | Mission Discrimination |
|---|---:|---:|---:|---:|
| A1_liver_bulk_lomo | 0.4873 | 0.4936 | 0.4942 | 0.8798 |
| A1_liver_bulk_lomo_combat | 0.4929 | 0.5053 | 0.5279 | 0.8684 |
| A1_liver_bulk_lomo_iss_only | 0.5391 | 0.5402 | 0.5633 | 0.8997 |
| A2_gastrocnemius_bulk_lomo | 0.3333 | 0.5000 | 0.7031 | 1.0000 |
| A3_kidney_bulk_lomo | 0.6062 | 0.6043 | 0.6519 | 0.6864 |
| A4_thymus_bulk_lomo | 0.5286 | 0.5557 | 0.6578 | 0.9950 |
| A5_skin_bulk_lomo | 0.7294 | 0.7418 | 0.7530 | 0.6569 |
| A6_eye_bulk_lomo | 0.5898 | 0.6074 | 0.7059 | 1.0000 |

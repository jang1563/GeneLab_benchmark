# Human Organoid PCA-LR Reconstructed Feature-Effect Pilot

Status: diagnostic feature-effect pilot; not a leaderboard result.

This report generates `feature_effect.csv.gz` from source-transfer
training samples only. Feature effects are PCA-LR coefficients
reconstructed into train-standardized gene space using
`pca.components_.T @ logistic.coef_[0]`. They are not log2FC
response signatures and must not be interpreted as expression fold
changes.

The DE reference table is not used during feature-effect generation;
it is used only afterward for rank, sign, and calibrated top-k
diagnostic scoring.

The accompanying `predictions.csv` is a two-row evaluator plumbing
fixture. Classification metrics in this report are not interpretable
as organoid model performance; only the feature-effect diagnostics are
reviewed.

- classifier model: `organoid_pca_lr_reconstructed_gene_space_feature_effect`
- feature-effect rows: 16000
- reference usage policy: `reference_not_used_for_effect_generation`
- claim boundary: feature-effect diagnostic only, not a biological
  generalization claim and not a leaderboard claim.

## Feature-Effect Metrics

| Metric | PCA-LR Reconstructed | L2 Logistic Reference |
|---|---:|---:|
| `feature_effect_direction_match` | 0.598039 | 0.607843 |
| `feature_effect_rank_correlation` | 0.086687 | 0.086728 |

## Top-K DE Overlap

| K | Overlap | Fraction | Expected | Enrichment | Hypergeometric p>=x |
|---:|---:|---:|---:|---:|---:|
| 50 | 1 | 0.020000 | 0.637500 | 1.568627 | 0.474071 |
| 100 | 5 | 0.050000 | 1.275000 | 3.921569 | 0.009090 |
| 250 | 10 | 0.040000 | 3.187500 | 3.137255 | 0.001414 |
| 500 | 14 | 0.028000 | 6.375000 | 2.196078 | 0.004966 |

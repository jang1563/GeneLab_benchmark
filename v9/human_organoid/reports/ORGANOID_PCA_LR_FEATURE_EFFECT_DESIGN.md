# Human Organoid PCA-LR Reconstructed Feature-Effect Design

Status: design complete
Date: 2026-05-23
Task: `draft_human_organoid_spaceflight`
Block: `V9-ORG-029`

## Scope

This note decides how to reconstruct PCA-LR classifier coefficients back into a
gene-space `feature_effect.csv` artifact for the draft human organoid task.

It does not implement the runner. The goal is to remove ambiguity before
adding a second classifier-derived feature-effect adapter.

## External API Checks

The design was checked against official scikit-learn documentation:

- [`StandardScaler`](https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html):
  standardizes each feature using training-set mean and standard deviation.
- [`PCA`](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html):
  `components_` has shape `(n_components, n_features)` and stores principal
  axes in feature space.
- [`LogisticRegression`](https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html):
  binary `coef_` has shape `(1, n_features)` and represents coefficients in the
  classifier decision function.

## Local Pipeline Boundary

The existing v9 PCA-LR bulk baseline in
`spacebio_bench/baselines/sklearn_classifiers.py` uses:

```text
StandardScaler()
PCA(n_components=adaptive_n, random_state=42)
LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000)
```

For the human organoid feature-effect adapter, the PCA-LR path should mirror
the current L2 logistic feature-effect leakage boundary rather than the bulk
LOMO file loader:

1. Use only source-transfer folds.
2. For an OSD-863 target contrast, train only on OSD-871 training samples.
3. For an OSD-871 target contrast, train only on OSD-863 training samples.
4. Do not use target-source expression rows.
5. Do not use DE references during feature-effect generation.
6. Use the DE reference only after artifact generation for diagnostic scoring.

## Proposed Training Pipeline

For each source-transfer fold:

1. Load train-source expression rows and binary labels.
2. Apply the same transform as the L2 logistic feature-effect pilot:
   `log1p` by default.
3. Select train-fold top variable genes, default `top_variable_genes=2000`.
4. Fit train-fold `StandardScaler` on selected genes.
5. Fit train-fold PCA on scaled selected genes with:
   `n_components=min(requested_components, n_train - 1, n_selected_features)`.
6. Fit binary logistic regression on train-fold PC scores.
7. Reconstruct gene-space coefficients from fitted PCA and logistic weights.
8. Emit reconstructed coefficients for each target-source DE contrast.

The first implementation should use `whiten=False`. Whitened PCA changes the
component scaling and should be a separate explicit design choice.

## Reconstruction Math

Let:

- `X` be the transformed selected-gene expression matrix.
- `Z = (X - scaler.mean_) / scaler.scale_` be train-standardized gene space.
- `C = pca.components_`, shape `(q, p)`, with `q` PCs and `p` selected genes.
- `m = pca.mean_`, shape `(p,)`.
- `T = (Z - m) @ C.T` be PC scores.
- `b_pc = logistic.coef_[0]`, shape `(q,)`, with positive class
  `LEO_or_ISS`.

The PC-space decision function is:

```text
logit = intercept + T @ b_pc
```

Substitute the PCA transform:

```text
logit = intercept + (Z - m) @ C.T @ b_pc
```

The standardized gene-space reconstructed coefficient is therefore:

```text
b_gene_standardized = C.T @ b_pc
```

The PCA mean only changes the intercept:

```text
intercept_adjusted = intercept - m @ b_gene_standardized
```

The feature-effect artifact should emit `b_gene_standardized`, not the adjusted
intercept. This keeps the scale parallel to the L2 logistic gene-space pilot:
classifier effect per train-standardized transformed gene.

## Output Contract

The adapter should reuse the existing `feature_effect.csv` required columns.

Recommended IDs and fields:

| Field | Value |
|---|---|
| `classifier_model_id` | `organoid_pca_lr_reconstructed_gene_space_feature_effect` |
| `effect_scale` | `pca_lr_reconstructed_standardized_logistic_coefficient` |
| `model_space` | `reconstructed_gene_space_from_pca` |
| `training_scope` | `source_transfer_organoid_type_holdout_train_samples` |
| `target_scope` | `target_source_contrast` |
| `reference_usage_policy` | `reference_not_used_for_effect_generation` |
| `artifact_claim_boundary` | `diagnostic_only_not_leaderboard` |

Optional metadata should include:

- `training_source_id`
- `target_source_id`
- `transform`
- `scaling`
- `feature_selection`
- `pca_components_requested`
- `pca_components_used`
- `pca_explained_variance_ratio_sum`
- `pca_whiten`
- `regularization`
- `classifier_solver`
- `random_state`
- `reconstruction_formula`

## Validation Requirements

The implementation must include unit tests that check:

1. Target-source samples are excluded from training.
2. Reconstructed coefficients have one row per selected gene per target
   contrast.
3. `b_gene_standardized` equals `pca.components_.T @ logistic.coef_[0]` on a
   small synthetic fit.
4. The artifact validates through `validate_feature_effect_artifact`.
5. `compute_feature_effect_metrics` reports calibrated top-k fields for the
   PCA-LR artifact.
6. The runner writes `feature_effect.csv.gz`, `feature_effect_metadata.json`,
   `metrics.json`, `run_manifest.json`, and `README.md`.

## Interpretation Rules

PCA-LR reconstructed coefficients are not:

- log2FC response signatures;
- DE estimates;
- causal gene effects;
- primary benchmark metrics.

They are:

- a linear reconstruction of a PC-space classifier;
- model-effect diagnostics in train-standardized transformed gene space;
- comparable to the L2 logistic feature-effect artifact only through the
  same post hoc diagnostic metrics.

PCA sign flips should not change the reconstructed gene-space vector because
the paired logistic coefficient flips with the PC direction. The product
`C.T @ b_pc` is the stable quantity to compare.

## Comparison Plan

The PCA-LR report should be compared against:

```text
v9/human_organoid/reports/logistic_feature_effect/metrics.json
```

Required comparison fields:

- `feature_effect_direction_match`
- `feature_effect_rank_correlation`
- top-k observed overlap
- top-k expected overlap
- top-k enrichment
- top-k hypergeometric p>=x
- per-contrast calibrated top-k behavior

The report should explicitly state whether PCA-LR reconstruction improves,
matches, or weakens the L2 logistic feature-effect diagnostic.

## Go/No-Go Decision

Go for one constrained pilot implementation.

Reason:

- the reconstruction math is linear and auditable;
- the current feature-effect schema already separates classifier effects from
  response signatures;
- calibrated top-k diagnostics can judge whether reconstruction adds signal;
- PCA-LR is an important baseline family elsewhere in the project.

Boundary:

- implement only as an optional secondary diagnostic;
- do not promote to primary metrics;
- do not write reconstructed weights into `response_signature.csv`;
- do not continue beyond one pilot if calibrated diagnostics are weaker than
  the L2 logistic feature-effect report.

## Next Block

`V9-ORG-030: Human organoid PCA-LR reconstructed feature-effect pilot`

Expected implementation:

- add `build_pca_lr_reconstructed_feature_effect` to
  `spacebio_bench/feature_effects.py`;
- add a CLI runner under `scripts/`;
- generate
  `v9/human_organoid/reports/pca_lr_feature_effect/`;
- compare calibrated diagnostics against
  `v9/human_organoid/reports/logistic_feature_effect/`.

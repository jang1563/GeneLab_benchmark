# Human Organoid Classifier Feature-Effect Contract

Status: design decision
Date: 2026-05-23
Task: `draft_human_organoid_spaceflight`
Block: `V9-ORG-025`

## Decision

Classifier-derived coefficients should not be written to
`response_signature.csv`.

They should use a separate artifact:

```text
feature_effect.csv
```

Reason:

- `response_signature.csv` currently means predicted
  `LEO_or_ISS - Ground` log2 fold-change.
- logistic-regression coefficients are discriminative model weights, not
  expression fold changes.
- PCA-LR coefficients may live in PCA component space unless they are explicitly
  projected back to gene space.
- mixing coefficients and log2FC values under one scorer would make the metric
  biologically ambiguous.

## External Checks

The relevant scikit-learn references are:

- LogisticRegression documents `coef_` as model coefficients whose shape depends
  on binary versus multiclass use:
  https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html
- The permutation-importance guide describes feature importance as a fitted
  model inspection result, not an intrinsic property of a feature:
  https://scikit-learn.org/stable/modules/permutation_importance.html
- The common-pitfalls guide emphasizes that preprocessing and model-building
  must be fit on training data only:
  https://scikit-learn.org/stable/common_pitfalls.html#data-leakage

These references support the v9 rule that coefficients may be diagnostic model
effects, but must not be presented as expression-response signatures.

## Artifact Name

```text
feature_effect.csv
```

This artifact is optional and diagnostic-only. It should be accepted only by a
separate scorer or report path, never by the existing response-signature log2FC
scorer.

## Required Columns

```text
task_id
fold_id
source_id
contrast_id
feature_namespace
gene_symbol
ensembl_id
feature_id
effect_value
effect_direction_positive_class
effect_scale
model_space
classifier_model_id
training_scope
target_scope
positive_label
reference_usage_policy
```

Column meaning:

| Column | Meaning |
|---|---|
| `effect_value` | signed model effect in the declared `effect_scale`. |
| `effect_direction_positive_class` | class whose odds/score increase when effect is positive, usually `LEO_or_ISS`. |
| `effect_scale` | e.g. `standardized_logistic_coefficient`, `pca_reconstructed_linear_weight`, or `permutation_importance_delta_metric`. |
| `model_space` | e.g. `gene_space`, `pca_component_space`, or `pca_reconstructed_gene_space`. |
| `reference_usage_policy` | must be `reference_not_used_for_effect_generation`. |

## Recommended Optional Columns

```text
transform
scaling
feature_selection
n_train_ground
n_train_leo_or_iss
n_features_model_input
n_features_emitted
regularization
classifier_solver
random_state
artifact_claim_boundary
```

## Allowed Initial Producers

### L2 Logistic Regression In Gene Space

This is the cleanest first coefficient producer.

Model:

```text
train-only transform -> train-only scaling -> L2 logistic regression
```

Feature-effect interpretation:

```text
standardized_logistic_coefficient
```

Positive values mean that higher train-standardized feature values increase the
model score or odds for `LEO_or_ISS`.

Boundary:

- not log2FC;
- not a DE estimate;
- valid only in the fitted train-fold feature/scaling context.

### PCA-LR Reconstructed Gene-Space Weight

Allowed later, but not as the first implementation.

If PCA-LR is used, the native classifier coefficients are in PCA component
space. A gene-space effect can be reconstructed only if the report records the
projection rule:

```text
gene_weight = PCA_components.T @ logistic_coefficients
```

If StandardScaler is used before PCA, the report must state whether weights are
still in standardized-gene units or transformed back to original feature units.

Boundary:

- reconstructed linear weight, not log2FC;
- sensitive to PCA components and train-fold preprocessing;
- rank-only interpretation is safer than magnitude interpretation.

### Permutation Importance

Allowed as a model-inspection follow-up, not the first coefficient artifact.

Permutation importance is model-agnostic, but it depends on the validation set,
the metric, repeats, and feature correlations. For the tiny human organoid task,
it should be treated as exploratory and reported with uncertainty.

## Scoring Policy

Do not run these metrics on `feature_effect.csv`:

- `de_direction_match` as currently defined for predicted log2FC;
- log2FC magnitude error;
- any metric whose name implies expression fold-change prediction.

Allowed diagnostic metrics:

1. `feature_effect_rank_correlation`
   - Spearman correlation between signed model effects and reference log2FC.
   - Uses joined gene rows.
   - Interpreted only as rank alignment.
2. `feature_effect_direction_match`
   - Sign match between model effect and reference log2FC.
   - Only valid when `effect_direction_positive_class=LEO_or_ISS`.
   - Must be labeled as model-effect sign agreement, not DE direction recovery.
3. `feature_effect_top_k_de_overlap`
   - Overlap between top absolute-effect genes and significant DE genes.
   - Recommended K values: 50, 100, 250, 500.
   - Useful when coefficient magnitudes are noisy but ranking is informative.

All three metrics remain non-primary and non-leaderboard for the current
organoid task.

## Leakage Boundary

Feature-effect generation must use:

- train-fold expression only;
- train-fold preprocessing only;
- train-fold labels only;
- no OSDR DE tables;
- no derived DE-reference table;
- no target-source expression rows for source-transfer reports.

The DE reference may be used only after artifact generation for diagnostic
scoring.

## Relationship To Current Response-Signature Diagnostics

Current log2FC response-signature diagnostics:

| Adapter | Role |
|---|---|
| `organoid_source_transfer_empirical_signature` | first conservative global diagnostic |
| `organoid_microglia_matched_source_transfer_empirical_signature` | secondary condition-sensitivity diagnostic |
| `organoid_shared_control_disease_microglia_source_transfer_empirical_signature` | partial-coverage negative/weak stratification diagnostic |

The classifier feature-effect artifact should answer a different question:

> Do discriminative model weights trained without reference leakage rank or sign
> genes similarly to the post hoc DE reference?

It should not answer:

> What log2 fold-change does the model predict for each gene?

## Recommended First Implementation

`V9-ORG-026: Human organoid L2 logistic feature-effect pilot`

Implement a source-transfer L2 logistic feature-effect report in gene space.

Recommended path:

```text
v9/human_organoid/reports/logistic_feature_effect/
```

Likely files:

- `spacebio_bench/feature_effects.py`
- `scripts/run_v9_human_organoid_logistic_feature_effect.py`
- `tests/test_v9_spacebio_bench.py`
- `v9/human_organoid/reports/logistic_feature_effect/README.md`

Recommended model:

```text
log1p normalized counts
train-fold z-score scaling
optional train-fold top variable genes
L2 LogisticRegression(class_weight="balanced")
```

The first implementation should use source-transfer folds:

- train OSD-871, emit feature effects for OSD-863 target source;
- train OSD-863, emit feature effects for OSD-871 target source.

This mirrors the non-leaky global source-transfer response-signature boundary.

## Required Tests

- `feature_effect.csv` validates required columns.
- target-source samples are excluded from model fitting.
- `response_signature.csv` scorer does not consume `feature_effect.csv`.
- feature-effect scorer computes rank/sign/top-k diagnostics on a synthetic
  reference fixture.
- report README states that coefficients are not log2FC signatures.

## Final Decision

Proceed with a separate `feature_effect.csv` contract. Do not overload
`response_signature.csv`. The next implementation should start with L2 logistic
gene-space coefficients because they have the clearest model-space semantics
among the available classifier options.

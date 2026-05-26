# OSD-120 Interaction L2 Logistic Baseline Design

Status: draft-only V9-MULTI-014 checkpoint  
Date: 2026-05-25

## Goal

Design the first stronger transparent baseline for the OSD-120 Arabidopsis root
interaction task without losing the separation between genotype/ecotype,
light-treatment, and condition-stratum fold families.

This baseline is still draft-only. It is not a leaderboard result and should be
compared against the nearest-centroid reports through fold-family and fold-detail
tables, not only aggregate metrics.

## External Implementation Context

The design follows the existing v9 sklearn bulk baseline pattern and current
scikit-learn behavior:

- `sklearn.linear_model.LogisticRegression` implements regularized logistic
  regression and supports L2 regularization with solvers including `liblinear`
  and `lbfgs`.
- `sklearn.preprocessing.StandardScaler` documents why scale differences can
  dominate estimator objectives, which supports keeping train-fold-only
  standardization explicit for high-dimensional expression features.

References:

- https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html
- https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html

## Inputs

- Task manifest:
  `v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json`
- Loader:
  `spacebio_bench.data.load_multispecies_interaction_task`
- Nearest-centroid baseline:
  `v9/multispecies/reports/interaction_nearest_centroid/`
- Sensitivity grid:
  `v9/multispecies/reports/interaction_sensitivity/multispecies_baseline_summary.csv`
- Fold-detail comparison table:
  `v9/multispecies/reports/interaction_sensitivity/fold_detail_summary.csv`

## Design Decision

Implement a dedicated OSD-120 interaction L2 logistic runner rather than routing
through the bulk `sklearn_classifiers` module.

Reason:

- The bulk sklearn runner assumes `BulkTaskData` and legacy fold files.
- OSD-120 already loads as `MultispeciesTaskData` with in-memory aligned
  expression features and explicit fold families.
- The interaction task needs fold-family separated output directories and
  delta metrics for genotype, light, and condition-stratum holdouts.

## Model Contract

Baseline id:

- `multispecies_interaction_logistic_regression_l2`

Default preprocessing:

- transform: `log1p`
- scaling: train-fold `zscore`
- feature selection: top 2,000 train-fold variable genes
- classifier: `LogisticRegression(penalty="l2", solver="liblinear", C=1.0,
  class_weight="balanced", max_iter=5000, random_state=42)`

The `liblinear` solver is selected to match the existing v9 bulk L2 logistic
baseline and because the task is binary with far more candidate features than
training samples.

## Output Contract

The runner should write one report per fold family:

- `primary_genotype_or_ecotype_holdout`
- `secondary_light_treatment_holdout`
- `diagnostic_condition_stratum_holdout`

Each report should include:

- `predictions.csv`
- `metrics.json`
- `run_manifest.json`

The summary should be:

- `v9/multispecies/reports/interaction_logistic_l2/multispecies_baseline_summary.csv`
- `v9/multispecies/reports/interaction_logistic_l2/multispecies_baseline_summary.json`

Rows should reuse the existing multispecies baseline summary schema so the
nearest-centroid and logistic reports are directly comparable.

## Comparison Rule

The logistic baseline must be compared against the nearest-centroid baseline
using:

1. Fold-family aggregate metrics.
2. Delta metrics for the applicable held-out factor.
3. `fold_detail_summary.csv` for repeated fragile strata.

Aggregate balanced accuracy alone is insufficient because V9-MULTI-012 showed
that condition-specific fragility can hide behind similar aggregate scores.

## Decision

Proceed with direct implementation. The adaptation is local to
`spacebio_bench/baselines/multispecies.py`, a small CLI, generated reports, and
tests. No task-manifest schema change is required.

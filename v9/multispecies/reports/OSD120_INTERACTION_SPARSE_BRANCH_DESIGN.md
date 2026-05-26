# OSD-120 Interaction Sparse Feature-Stability Branch Design

Status: draft-only V9-MULTI-018 design checkpoint  
Date: 2026-05-25

## Scope

This note designs the next OSD-120 modeling branch after the L2 logistic
feature-set audit. It does not run a new model. It defines the smallest
transparent sparse branch worth implementing and the gates it must pass before
any more complex classifier is considered.

## Local Evidence

The current diagnostic surface is:

- nearest-centroid default `tvg2000_log1p_zscore`
- L2 logistic default `tvg2000_log1p_zscore_c1`
- L2 logistic sensitivity `tvg500_log1p_zscore_c1`
- feature-set audit comparing top 500 and top 2,000 selected genes

Key facts from V9-MULTI-016/017:

- top 2,000 L2 logistic preserves `Light.Treatment` better but keeps
  `Col.0.PhyD|Dark.Treatment` at balanced accuracy 0.3333.
- top 500 L2 logistic restores `Col.0.PhyD|Dark.Treatment` to 0.6667 but drops
  `Light.Treatment` to 0.5000.
- top 500 is nested inside top 2,000 by train-fold variance rank, yet top-10
  coefficient overlap is only 1/10 to 3/10 across focus folds.

This argues for a sparse branch that can admit the broader top-2,000 candidate
space while limiting coefficient turnover from extra features.

## External Method Context

scikit-learn `LogisticRegression` supports L1 and L2 regularization with the
`liblinear` solver for binary problems, while elastic-net requires `saga`.
Smaller `C` means stronger regularization. The documentation also notes that
`sag`/`saga` convergence depends on similarly scaled features, matching the
existing train-fold `zscore` preprocessing rule.

Stability selection is a general high-dimensional variable-selection wrapper
based on subsampling plus a selection algorithm. It is attractive for this
setting because the scientific risk is unstable feature admission, but it
should be a second step after a small sparse pilot establishes useful candidate
regularization strengths.

References:

- https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html
- https://arxiv.org/abs/0809.2932

## Selected Branch

Implement a first sparse L1 logistic pilot before elastic-net or non-linear
models.

Rationale:

- OSD-120 is binary, so `liblinear` can reuse the current logistic dependency
  path with minimal new solver risk.
- L1 directly tests whether sparse coefficients can keep the useful top-2,000
  candidate space while suppressing unstable extra-feature dominance.
- Elastic-net and stability-selection are still valid, but they should be
  gated by the simpler sparse pilot outcome.

## Candidate Configs

Fixed settings:

- transform: `log1p`
- scaling: train-fold `zscore`
- top variable genes: 2,000
- classifier: logistic regression with L1 semantics, `solver="liblinear"`,
  `class_weight="balanced"`, `max_iter=5000`, `random_state=42`
- implementation note: use `l1_ratio=1` with current scikit-learn rather than
  the deprecated `penalty="l1"` argument

Grid:

- `C`: 0.01, 0.03, 0.1, 0.3, 1.0
- fold families: primary genotype/ecotype, secondary light-treatment,
  diagnostic condition-stratum

Optional fallback only if L1 is empty or unstable:

- elastic-net logistic with `solver="saga"`, `l1_ratio`: 0.25, 0.5, 0.75,
  `C`: 0.03, 0.1, 0.3

## Required Outputs

The L1 pilot should mirror the existing logistic report shape:

- `multispecies_baseline_summary.csv` and `.json`
- per-variant `predictions.csv`, `metrics.json`, and `run_manifest.json`
- `fold_detail_summary.csv` and `.json`
- `fold_detail_comparison_vs_nearest_centroid.csv` and `.json`
- sparse coefficient audit with non-zero feature counts per fold
- short review note

## Gates

The branch is useful only if at least one sparse variant satisfies all of:

- `Col.0.PhyD|Dark.Treatment` balanced accuracy >= 0.5000, ideally 0.6667.
- `Light.Treatment` balanced accuracy > 0.5556, ideally >= 0.6667.
- `Col.0.PhyD` balanced accuracy >= 0.5833.
- diagnostic condition-stratum delta <= 0.3333.
- secondary light-treatment delta <= 0.1111.
- no focus fold has an empty fitted model.
- non-zero coefficient counts are recorded and bounded enough to interpret.

If no variant passes these gates, the branch should be recorded as a negative
diagnostic and the next step should move to stability selection or pathway-level
features rather than more scalar `C` tuning.

## Test Plan

Implementation should add tests for:

- L1 config validation and variant ids.
- a small temp-output CLI/API smoke run over one fold family.
- generated summary row counts.
- generated fold-detail comparison row counts.
- generated sparse coefficient audit row counts and non-zero coefficient counts.
- gate summary values for `Col.0.PhyD|Dark.Treatment`, `Light.Treatment`, and
  `Col.0.PhyD`.

## Decision

Proceed next with a draft-only sparse L1 logistic pilot. Do not implement a
random forest, SVM, neural model, or promoted benchmark claim until the sparse
transparent branch either passes or fails the fragile-fold gates.

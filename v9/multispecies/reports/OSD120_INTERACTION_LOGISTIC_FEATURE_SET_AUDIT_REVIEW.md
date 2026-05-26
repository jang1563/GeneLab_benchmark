# OSD-120 Interaction Logistic Feature-Set Audit Review

Status: draft-only V9-MULTI-017 checkpoint  
Date: 2026-05-25

## Scope

This review audits why the V9-MULTI-016 logistic sensitivity grid trades off
two surfaces: top 500 genes improves `Col.0.PhyD|Dark.Treatment`, while top
2,000 genes improves `Light.Treatment`. The audit is model-diagnostic only.
Classifier coefficients are not treated as biological differential-expression
claims.

Focus folds:

- `Col.0.PhyD|Dark.Treatment`
- `Col.0.PhyD`
- `Light.Treatment`

Compared variants:

- `tvg500_log1p_zscore_c1`
- `tvg2000_log1p_zscore_c1`

Outputs:

- `v9/multispecies/reports/interaction_logistic_feature_audit/feature_set_audit_summary.csv`
- `v9/multispecies/reports/interaction_logistic_feature_audit/feature_coefficient_audit.csv`

## Summary

| Held-out fold | Top 500 BA | Top 2,000 BA | Top 2,000 extra features | Top-10 coefficient overlap |
|---|---:|---:|---:|---:|
| `Col.0.PhyD|Dark.Treatment` | 0.6667 | 0.3333 | 1,500 | 2/10 |
| `Col.0.PhyD` | 0.6667 | 0.5833 | 1,500 | 1/10 |
| `Light.Treatment` | 0.5000 | 0.7222 | 1,500 | 3/10 |

The top-500 feature set is nested inside top 2,000 by train-fold variance rank,
but the fitted coefficient surface changes strongly after the extra 1,500 genes
enter. Top-10 absolute coefficient overlap is only 1-3 genes across the three
focus folds.

Examples of top-2,000 high-magnitude extra-feature coefficients:

| Held-out fold | Top extra coefficient features |
|---|---|
| `Col.0.PhyD|Dark.Treatment` | `AT5G48010`, `AT1G71380`, `AT2G26400` |
| `Col.0.PhyD` | `AT5G04700`, `AT5G48010`, `AT5G11180` |
| `Light.Treatment` | `AT5G04700`, `AT3G01345`, `AT5G42600` |

## Interpretation

The V9-MULTI-016 tradeoff is not explained by L2 regularization strength. It is
explained by feature admission. The extra 1,500 top-variance genes are not just
small perturbations; they can dominate the fitted coefficient ranking.

This gives a concrete reason not to promote either setting as a resolved model:

- top 500 genes gives the cleaner dark condition-stratum behavior but weakens
  light-treatment holdout behavior.
- top 2,000 genes gives the cleaner light-treatment behavior but reintroduces
  the `Col.0.PhyD|Dark.Treatment` failure.

The current transparent logistic family is therefore useful as a diagnostic,
but it does not yet provide a single robust OSD-120 interaction baseline.

## Decision

Keep `tvg2000_log1p_zscore_c1` as the direct comparator to the established
nearest-centroid default, and keep `tvg500_log1p_zscore_c1` as a documented
feature-count sensitivity comparator.

The next block should design the next modeling branch rather than tune `C`
again. A reasonable path is a small, draft-only sparse/feature-stability
experiment, with explicit gates: it must improve `Light.Treatment` without
recreating the `Col.0.PhyD|Dark.Treatment` failure, and it must preserve
fold-family separated reporting.

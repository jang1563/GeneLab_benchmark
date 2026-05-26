# OSD-120 Interaction Sparse L1 Stability Review

Status: draft-only V9-MULTI-020 checkpoint  
Date: 2026-05-26

## Scope

This review audits selected-feature stability for the useful sparse L1 logistic
region from V9-MULTI-019. It does not change held-out performance estimates.
For each focus fold, the audit fixes the full train-fold top-2,000 candidate
feature set, then refits sparse L1 logistic models on 20 deterministic balanced
subsamples of the training fold.

Focus variants:

- `tvg2000_log1p_zscore_l1_c0p3`
- `tvg2000_log1p_zscore_l1_c1`

Focus folds:

- `Col.0.PhyD|Dark.Treatment`
- `Col.0.PhyD`
- `Light.Treatment`

Outputs:

- `v9/multispecies/reports/interaction_logistic_sparse_l1_stability/stability_summary.csv`
- `v9/multispecies/reports/interaction_logistic_sparse_l1_stability/stability_feature_frequency.csv`

## Stability Summary

| Variant | Fold | Ref BA | Ref nonzero | Median subsample nonzero | Stable >=0.5 | Stable >=0.8 | Mean Jaccard |
|---|---|---:|---:|---:|---:|---:|---:|
| `c0p3` | `Col.0.PhyD|Dark.Treatment` | 0.6667 | 9 | 7.0 | 5 | 4 | 0.4089 |
| `c1` | `Col.0.PhyD|Dark.Treatment` | 0.6667 | 10 | 12.5 | 10 | 7 | 0.4625 |
| `c0p3` | `Col.0.PhyD` | 0.9167 | 7 | 5.0 | 2 | 1 | 0.2762 |
| `c1` | `Col.0.PhyD` | 0.9167 | 12 | 9.0 | 7 | 2 | 0.3239 |
| `c0p3` | `Light.Treatment` | 0.8333 | 3 | 3.0 | 1 | 1 | 0.2823 |
| `c1` | `Light.Treatment` | 0.8333 | 10 | 6.0 | 2 | 1 | 0.1907 |

Top stable features by fold:

- `Col.0.PhyD|Dark.Treatment`, `c1`: `AT1G54970` at 1.00 selection frequency,
  plus `AT1G16440`, `AT1G71380`, `AT2G30660`, and `AT5G52700` at 0.95.
- `Col.0.PhyD`, `c1`: `AT2G30660` at 1.00 and `AT1G16440` at 0.90.
- `Light.Treatment`, both variants: `AT5G04700` at 1.00.

These are classifier-stability signals only. They are not differential
expression or biological mechanism claims.

## Interpretation

The sparse L1 candidate is not perfectly stable, but it is much more coherent
than the earlier L2 top-500/top-2,000 tradeoff. The useful feature-selection
signal is fold-specific:

- `C=1.0` is clearly stronger for the dark diagnostic fold and the primary
  `Col.0.PhyD` fold. It has more stable selected features and higher pairwise
  Jaccard than `C=0.3` on both.
- `C=0.3` is more compact for `Light.Treatment`: only three reference nonzero
  coefficients, with `AT5G04700` selected in every subsample.
- `C=1.0` keeps the best held-out gate profile from V9-MULTI-019, but its
  light-treatment stability is more diffuse than `C=0.3`.

The right conclusion is not that the sparse branch is final. It is that sparse
L1 now deserves to be the leading transparent OSD-120 diagnostic comparator,
with `C=1.0` as the performance-leading setting and `C=0.3` as a compact
stability comparator.

## Decision

Keep `tvg2000_log1p_zscore_l1_c1` as the leading transparent diagnostic
candidate because it clears the fragile-fold gates and worsens zero default
nearest-centroid fold rows.

Keep `tvg2000_log1p_zscore_l1_c0p3` as the compact sparse comparator,
especially for light-treatment stability.

The next block should consolidate the OSD-120 interaction baseline ladder into
a single comparison note: nearest-centroid default, L2 logistic default, L2
top-500 sensitivity, sparse L1 `c0p3`, and sparse L1 `c1`.

# OSD-120 Interaction Sparse L1 Logistic Review

Status: draft-only V9-MULTI-019 checkpoint  
Date: 2026-05-25

## Scope

This review evaluates the sparse L1 logistic pilot selected in
`OSD120_INTERACTION_SPARSE_BRANCH_DESIGN.md`. The pilot is still diagnostic and
draft-only. It is compared against the default nearest-centroid fold-detail
baseline and the earlier L2 logistic sensitivity results.

## Outputs

- `v9/multispecies/reports/interaction_logistic_sparse_l1/multispecies_baseline_summary.csv`
- `v9/multispecies/reports/interaction_logistic_sparse_l1/fold_detail_summary.csv`
- `v9/multispecies/reports/interaction_logistic_sparse_l1/fold_detail_comparison_vs_nearest_centroid.csv`
- `v9/multispecies/reports/interaction_logistic_sparse_l1/feature_set_audit_summary.csv`
- `v9/multispecies/reports/interaction_logistic_sparse_l1/feature_coefficient_audit.csv`

## Aggregate Results

| Variant | Primary BA | Secondary BA | Diagnostic BA | Light delta | Diagnostic delta |
|---|---:|---:|---:|---:|---:|
| `tvg2000_log1p_zscore_l1_c0p01` | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.0000 |
| `tvg2000_log1p_zscore_l1_c0p03` | 0.5000 | 0.5000 | 0.5000 | 0.0000 | 0.0000 |
| `tvg2000_log1p_zscore_l1_c0p1` | 0.6111 | 0.5000 | 0.7500 | 0.0000 | 0.3333 |
| `tvg2000_log1p_zscore_l1_c0p3` | 0.8611 | 0.7778 | 0.9167 | 0.1111 | 0.3333 |
| `tvg2000_log1p_zscore_l1_c1` | 0.9167 | 0.8333 | 0.8889 | 0.0000 | 0.3333 |

## Fragile-Fold Gates

| Variant | `Col.0.PhyD|Dark.Treatment` | `Light.Treatment` | `Col.0.PhyD` | Non-zero focus coefficients |
|---|---:|---:|---:|---|
| `c0p01` | 0.5000 | 0.5000 | 0.5000 | 0 / 0 / 0 |
| `c0p03` | 0.5000 | 0.5000 | 0.5000 | 0 / 0 / 0 |
| `c0p1` | 0.5000 | 0.5000 | 0.6667 | 3 / 0 / 1 |
| `c0p3` | 0.6667 | 0.8333 | 0.9167 | 9 / 3 / 7 |
| `c1` | 0.6667 | 0.8333 | 0.9167 | 10 / 10 / 12 |

The non-zero coefficient counts are ordered as
`Col.0.PhyD|Dark.Treatment / Light.Treatment / Col.0.PhyD`.

## Nearest-Centroid Comparison

Against the default nearest-centroid `tvg2000_log1p_zscore` fold-detail rows:

| Variant | Improved folds | Tied folds | Worsened folds |
|---|---:|---:|---:|
| `c0p01` | 0 | 3 | 8 |
| `c0p03` | 0 | 3 | 8 |
| `c0p1` | 3 | 4 | 4 |
| `c0p3` | 8 | 2 | 1 |
| `c1` | 9 | 2 | 0 |

## Interpretation

Very strong L1 regularization (`C=0.01` and `C=0.03`) collapses to empty models
on the focus folds and should be treated as negative diagnostics. `C=0.1` starts
to recover primary/diagnostic behavior but still fails the light-treatment gate.

The two useful variants are `C=0.3` and `C=1.0`. Both restore
`Col.0.PhyD|Dark.Treatment` to 0.6667 and keep `Light.Treatment` at 0.8333.
`C=1.0` is the strongest candidate because it has the best primary and
secondary aggregate balanced accuracy, zero worsened folds versus
nearest-centroid, and non-empty sparse coefficient sets on all focus folds.

This is the first OSD-120 logistic variant that clears the V9-MULTI-018
fragile-fold gates without the top-500/top-2,000 tradeoff seen in L2.

## Decision

Keep `tvg2000_log1p_zscore_l1_c1` as the current best transparent OSD-120
diagnostic candidate. Do not promote it to a benchmark claim yet.

V9-MULTI-020 follows this checkpoint with a train-fold subsampling stability
audit around `C=0.3` and `C=1.0`. That stability result should be used before
consolidating the final OSD-120 diagnostic ladder.

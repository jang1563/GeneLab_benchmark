# OSD-120 Interaction L2 Logistic Sensitivity Review

Status: draft-only V9-MULTI-016 checkpoint  
Date: 2026-05-25

## Scope

This review probes whether the V9-MULTI-015 logistic failure on
`Col.0.PhyD|Dark.Treatment` is stable under a compact L2 logistic sensitivity
grid. It keeps the OSD-120 interaction task draft-only and compares every
logistic variant against the default nearest-centroid fold-detail baseline,
`tvg2000_log1p_zscore`.

## Grid

Fixed settings:

- transform: `log1p`
- scaling: train-fold `zscore`
- classifier: L2 logistic regression, `solver="liblinear"`,
  `class_weight="balanced"`, `max_iter=5000`, `random_state=42`

Varied settings:

- top variable genes: 500, 2,000
- inverse regularization strength `C`: 0.1, 1.0, 10.0
- fold families: primary genotype/ecotype, secondary light-treatment,
  diagnostic condition-stratum

Outputs:

- `v9/multispecies/reports/interaction_logistic_l2_sensitivity/multispecies_baseline_summary.csv`
- `v9/multispecies/reports/interaction_logistic_l2_sensitivity/fold_detail_summary.csv`
- `v9/multispecies/reports/interaction_logistic_l2_sensitivity/fold_detail_comparison_vs_nearest_centroid.csv`

## Aggregate Results

| Variant | Primary BA | Secondary BA | Diagnostic BA | Diagnostic delta |
|---|---:|---:|---:|---:|
| `tvg500_log1p_zscore_c0p1` | 0.8333 | 0.6667 | 0.8611 | 0.3333 |
| `tvg500_log1p_zscore_c1` | 0.8333 | 0.6667 | 0.8611 | 0.3333 |
| `tvg500_log1p_zscore_c10` | 0.8333 | 0.6944 | 0.8611 | 0.3333 |
| `tvg2000_log1p_zscore_c0p1` | 0.7778 | 0.7222 | 0.8611 | 0.6667 |
| `tvg2000_log1p_zscore_c1` | 0.7778 | 0.7500 | 0.8611 | 0.6667 |
| `tvg2000_log1p_zscore_c10` | 0.7778 | 0.7500 | 0.8611 | 0.6667 |

## Fragile-Fold Results

| Variant group | `Col.0.PhyD|Dark.Treatment` BA | `Light.Treatment` BA | `Col.0.PhyD` BA | Readout |
|---|---:|---:|---:|---|
| top 500 genes, C 0.1/1.0 | 0.6667 | 0.5000 | 0.6667 | dark stratum improves, light holdout weakens |
| top 500 genes, C 10.0 | 0.6667 | 0.5556 | 0.6667 | dark stratum improves, light holdout only matches nearest centroid |
| top 2,000 genes, C 0.1 | 0.3333 | 0.6667 | 0.5833 | dark stratum fails, light remains acceptable |
| top 2,000 genes, C 1.0/10.0 | 0.3333 | 0.7222 | 0.5833 | dark stratum fails, light is strongest |

Across all six variants, logistic improves 8/11 nearest-centroid fold rows for
five variants and 8/11 with three ties for `tvg500_log1p_zscore_c10`. The
important change is not regularization strength but feature-count choice.

## Interpretation

The `Col.0.PhyD|Dark.Treatment` failure is not an invariant property of L2
logistic regression. It persists across all top-2,000 variants, but disappears
when the model is restricted to 500 train-fold variable genes. That fix is not
free: the top-500 variants weaken the secondary `Light.Treatment` holdout,
including two variants at balanced accuracy 0.5000.

This creates a genuine tradeoff. Top 2,000 genes gives the cleaner light
treatment behavior, while top 500 genes gives the cleaner diagnostic
condition-stratum spread. No tested logistic variant resolves both surfaces at
once.

## Decision

Keep the original top-2,000 L2 logistic baseline as the direct comparator to the
nearest-centroid default, but mark the top-500 sensitivity result as important
evidence that the dark-stratum failure is feature-set sensitive.

The next block should audit selected feature sets and fold-level coefficient
behavior for the top-500 versus top-2,000 logistic variants before trying a more
complex classifier.

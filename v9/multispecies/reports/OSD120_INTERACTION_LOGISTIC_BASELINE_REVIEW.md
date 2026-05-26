# OSD-120 Interaction L2 Logistic Baseline Review

Status: draft-only V9-MULTI-014/015 implementation checkpoint  
Date: 2026-05-25

## Scope

This review covers the first L2 logistic regression baseline for the OSD-120
Arabidopsis root interaction task. It follows the design recorded in
`OSD120_INTERACTION_LOGISTIC_BASELINE_DESIGN.md` and preserves the three
fold-family surfaces:

- primary genotype/ecotype holdout
- secondary light-treatment holdout
- diagnostic condition-stratum holdout

The result is diagnostic only. It is not a leaderboard result and should not be
merged into the OSD-37/OSD-207 species-native reports.

## Implementation

Baseline id:

- `multispecies_interaction_logistic_regression_l2`

Default configuration:

- transform: `log1p`
- scaling: train-fold `zscore`
- feature selection: top 2,000 train-fold variable genes
- classifier: `LogisticRegression` with L2 penalty, `solver="liblinear"`,
  `C=1.0`, `class_weight="balanced"`, `max_iter=5000`, and `random_state=42`

Outputs:

- `v9/multispecies/reports/interaction_logistic_l2/multispecies_baseline_summary.csv`
- `v9/multispecies/reports/interaction_logistic_l2/fold_detail_summary.csv`
- `v9/multispecies/reports/interaction_logistic_l2/fold_detail_comparison_vs_nearest_centroid.csv`
- per-fold-family `predictions.csv`, `metrics.json`, and `run_manifest.json`

## Summary

| Fold family | Logistic BA | Nearest-centroid BA | Delta | Logistic AUROC | Nearest-centroid AUROC | Logistic holdout delta | Nearest-centroid holdout delta |
|---|---:|---:|---:|---:|---:|---:|---:|
| Primary genotype/ecotype | 0.7778 | 0.6667 | +0.1111 | 0.8765 | 0.7346 | 0.3333 | 0.2500 |
| Secondary light-treatment | 0.7500 | 0.6667 | +0.0833 | 0.8488 | 0.7840 | 0.0556 | 0.2222 |
| Diagnostic condition-stratum | 0.8611 | 0.6667 | +0.1944 | 0.8981 | 0.7654 | 0.6667 | 0.3333 |

The L2 logistic baseline improves aggregate balanced accuracy and AUROC across
all three fold families. This is expected for a supervised linear classifier
relative to the nearest-centroid rule.

The improvement is not uniformly reassuring. The diagnostic
condition-stratum holdout delta doubles from 0.3333 to 0.6667 because one
stratum remains very weak while several others become perfect:

- `Col.0.PhyD|Dark.Treatment`: balanced accuracy 0.3333
- `Col.0|Dark.Treatment`: balanced accuracy 0.8333
- all other diagnostic strata: balanced accuracy 1.0000

Primary genotype/ecotype also still has a weak held-out genotype:

- `Col.0.PhyD`: balanced accuracy 0.5833
- `Col.0`: balanced accuracy 0.9167
- `Wassilewskija.ecotype`: balanced accuracy 0.8333

Secondary light-treatment improves the most cleanly: both light folds are above
0.70 and the light holdout delta drops to 0.0556.

## Fold-Detail Comparison

V9-MULTI-015 compares the logistic fold-detail table against the default
nearest-centroid sensitivity variant, `tvg2000_log1p_zscore`. The comparison
covers 11 held-out folds: logistic improves 8, ties 2, and worsens 1.

| Fold | Nearest-centroid BA | Logistic BA | Delta | Interpretation |
|---|---:|---:|---:|---|
| `Col.0.PhyD` | 0.5000 | 0.5833 | +0.0833 | still the weakest primary genotype fold |
| `Light.Treatment` | 0.5556 | 0.7222 | +0.1667 | secondary light holdout becomes more stable |
| `Col.0.PhyD|Dark.Treatment` | 0.5000 | 0.3333 | -0.1667 | only fold that gets worse |
| `Wassilewskija.ecotype|Dark.Treatment` | 0.5000 | 1.0000 | +0.5000 | nearest-centroid weak fold is resolved by logistic |

The side-by-side table is intentionally machine-readable so future OSD-120
models can be checked against repeated fragile folds instead of only aggregate
balanced accuracy.

## Interpretation

The logistic baseline is a useful stronger transparent comparator. It suggests
that linear discriminative training captures more of the OSD-120 signal than
the nearest-centroid baseline. But it also sharpens the same scientific warning
from V9-MULTI-012: aggregate performance can hide condition-stratum fragility.

The `Col.0.PhyD|Dark.Treatment` diagnostic fold should be tracked in every
future OSD-120 model comparison. A model that improves aggregate diagnostic
balanced accuracy while failing this fold should not be described as resolving
the interaction structure.

## Decision

Keep the L2 logistic report as the first stronger transparent diagnostic
baseline. Do not promote it to a primary benchmark claim.

V9-MULTI-016 follows this checkpoint with a compact regularization and
feature-count sensitivity grid. Its conclusion should be consulted before
moving to a more complex model class.

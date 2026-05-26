# OSD-120 Interaction Sensitivity Review

Status: draft-only V9-MULTI-012 checkpoint  
Date: 2026-05-23

## Scope

This review checks whether the first OSD-120 interaction nearest-centroid
baseline is stable under simple preprocessing choices. It does not promote
OSD-120 to a leaderboard task. The task remains a small, single-source
Arabidopsis root interaction pilot with explicit genotype/ecotype and
light-treatment structure.

External context still supports treating OSD-120 separately from the
species-native OSD-37/OSD-207 pilots. NASA OSDR describes GeneLab/OSDR as a
standardized repository for space-relevant omics across organisms including
plants, and the 2024 CARA Arabidopsis transcriptomics paper treats light and
genotype/ecotype as central response axes in OSD-120.

- NASA OSDR overview: https://science.nasa.gov/biological-physical/data/osdr/
- Zhou, Ferl, and Paul 2024, npj Microgravity:
  https://www.nature.com/articles/s41526-024-00417-0

## Inputs

- Baseline feasibility report:
  `v9/multispecies/reports/interaction_nearest_centroid/`
- Sensitivity summary:
  `v9/multispecies/reports/interaction_sensitivity/multispecies_baseline_summary.csv`
- Fold-detail summary:
  `v9/multispecies/reports/interaction_sensitivity/fold_detail_summary.csv`
- Task manifest:
  `v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json`
- Review predecessor:
  `v9/multispecies/reports/OSD120_INTERACTION_BASELINE_FEASIBILITY_REVIEW.md`

The grid covers 20 preprocessing variants across all three fold families, for
60 evaluated report rows:

- transform: `log1p`, `none`
- scaling: train-fold `zscore`, `none`
- top train-fold variable genes: 100, 500, 2,000, 5,000, 24,740
- fold families: primary genotype/ecotype, secondary light-treatment,
  diagnostic condition-stratum

## Fold-Family Summary

| Fold family | Mean balanced accuracy | Range | Mean AUROC | Range | Mean delta | Delta range |
|---|---:|---:|---:|---:|---:|---:|
| Primary genotype/ecotype | 0.6792 | 0.5833-0.8056 | 0.7093 | 0.6512-0.8025 | 0.1750 | 0.0833-0.3333 |
| Secondary light-treatment | 0.6250 | 0.5833-0.6667 | 0.6912 | 0.5772-0.7840 | 0.1278 | 0.0556-0.2222 |
| Diagnostic condition-stratum | 0.6611 | 0.5278-0.7500 | 0.7085 | 0.6049-0.8117 | 0.3167 | 0.1667-0.5000 |

The default setting from V9-MULTI-011, `log1p` + train-fold `zscore` + top
2,000 genes, is not the best-performing sensitivity variant:

| Fold family | Default balanced accuracy | Best variant | Best balanced accuracy |
|---|---:|---|---:|
| Primary genotype/ecotype | 0.6667 | `tvg100_log1p_zscore` | 0.8056 |
| Secondary light-treatment | 0.6667 | `tvg100_log1p_zscore` | 0.6667 |
| Diagnostic condition-stratum | 0.6667 | `tvg500_log1p_zscore` | 0.7500 |

However, the better scores come from smaller feature subsets and do not remove
the underlying interaction fragility. They should be treated as a baseline
tuning observation, not as a new benchmark claim.

## Fold-Detail Table

V9-MULTI-013 adds `fold_detail_summary.csv` and `.json` so the fold-level
fragility below is machine-readable. The table contains 220 rows: 60 primary
genotype/ecotype rows, 40 secondary light-treatment rows, and 120 diagnostic
condition-stratum rows. Each row records variant id, fold family, delta metric,
held-out factor/value, fold balanced accuracy, within-variant low-to-high rank,
and flags for lowest fold and balanced accuracy <= 0.50.

This follows the same design pressure as Frictionless tabular data resources:
CSV outputs should have stable fields and metadata-ready structure, rather than
requiring consumers to parse nested report files ad hoc.

## Repeated Fragile Strata

Across 20 variants, the most repeatedly weak held-out strata were:

| Fold family | Held-out value | Times weakest | Times <= 0.50 BA | Mean BA | Range |
|---|---|---:|---:|---:|---:|
| Primary genotype/ecotype | `Wassilewskija.ecotype` | 12/20 | 4/20 | 0.6417 | 0.5000-0.7500 |
| Primary genotype/ecotype | `Col.0.PhyD` | 9/20 | 2/20 | 0.6333 | 0.5000-0.8333 |
| Secondary light-treatment | `Light.Treatment` | 19/20 | 0/20 | 0.5639 | 0.5556-0.6111 |
| Diagnostic condition-stratum | `Wassilewskija.ecotype|Dark.Treatment` | 16/20 | 17/20 | 0.5417 | 0.5000-0.8333 |
| Diagnostic condition-stratum | `Col.0.PhyD|Dark.Treatment` | 9/20 | 10/20 | 0.5833 | 0.5000-0.6667 |

This confirms the V9-MULTI-011 concern: aggregate balanced accuracy hides
condition-specific weakness. The repeated low performance of
`Light.Treatment`, `Wassilewskija.ecotype|Dark.Treatment`, and
`Col.0.PhyD|Dark.Treatment` suggests the model is sensitive to the interaction
structure, not merely to one unlucky preprocessing choice.

## Decision

Keep the conservative default `log1p` + train-fold `zscore` + top 2,000 genes
for now, even though smaller 100- or 500-gene subsets often score higher. The
default remains consistent with the OSD-37/OSD-207 multispecies baseline and is
less likely to overfit this small 36-sample interaction pilot. The sensitivity
grid should be cited as diagnostic evidence that feature-count tuning can change
scores substantially.

The OSD-120 task should keep three separate evaluation surfaces:

1. Primary genotype/ecotype holdout as the main pilot split.
2. Secondary light-treatment holdout as a condition-transfer diagnostic.
3. Condition-stratum holdout as a small-sample stress diagnostic only.

## Next Step

The next useful block can now consider a simple L2 logistic baseline with the
same fold-family separation and the same draft-only claim boundary. It should
reuse `fold_detail_summary.csv` as the comparison table for repeated fragile
strata rather than relying only on aggregate metrics.

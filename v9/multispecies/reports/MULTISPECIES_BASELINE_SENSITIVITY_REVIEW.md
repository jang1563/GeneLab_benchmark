# Multispecies Baseline Sensitivity Review

Status: V9-MULTI-008 review  
Date: 2026-05-23

## Scope

This block tests whether the V9-MULTI-007 nearest-centroid feasibility results
are stable across simple preprocessing choices.

Tasks:

- `draft_osd37_arabidopsis_seedling_spaceflight`
- `draft_osd207_drosophila_whole_body_spaceflight`

Grid:

- transform: `log1p`, `none`
- scaling: train-fold `zscore`, `none`
- top variable genes: 100, 500, 2,000, 5,000, and all available features per
  task

All outputs remain draft-only and not leaderboard claims.

## Generated Artifacts

- `v9/multispecies/reports/sensitivity/multispecies_baseline_summary.csv`
- `v9/multispecies/reports/sensitivity/multispecies_baseline_summary.json`
- per-variant prediction, metrics, and run-manifest directories under
  `v9/multispecies/reports/sensitivity/`

Command:

```bash
/usr/local/bin/python3 scripts/run_v9_multispecies_sensitivity.py
```

The grid produced 40 evaluated rows: 20 variants for OSD-37 and 20 variants for
OSD-207.

## Aggregate Readout

| task_id | balanced_accuracy min | median | max | AUROC min | median | max | holdout_delta min | median | max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| draft_osd37_arabidopsis_seedling_spaceflight | 0.5357142857 | 0.7678571429 | 0.8571428571 | 0.6058673469 | 0.8781887755 | 0.9234693878 | 0.0833333333 | 0.1666666667 | 0.5 |
| draft_osd207_drosophila_whole_body_spaceflight | 0.75 | 0.875 | 0.9375 | 0.83984375 | 0.919921875 | 0.96875 | 0.125 | 0.375 | 0.5 |

## OSD-37 Interpretation

OSD-37 remains a valid and relatively stable first plant feasibility task, but
the preprocessing choice matters.

Mean balanced accuracy by transform/scaling group:

- `none` + `zscore`: 0.8107142857
- `log1p` + `none`: 0.7428571428
- `log1p` + `zscore`: 0.7321428571
- `none` + `none`: 0.5357142857

The default V9-MULTI-007 configuration (`log1p`, `zscore`, top 2,000 genes)
scores `balanced_accuracy=0.8214285714` with
`condition_stratum_holdout_delta=0.125`. It is not the highest-scoring variant,
but it is a conservative, stable setting. The best balanced-accuracy variant is
`tvg500_none_zscore` at 0.8571428571, but its fold range is larger
(`condition_stratum_holdout_delta=0.3125`).

Weakest heldout stratum frequency across the 20 variants:

- Col.0: 16 variants
- Ler.0: 4 variants

Decision: keep OSD-37 as the cleaner first plant example. Keep `log1p` +
train-fold `zscore` + top 2,000 genes as the conservative default for now
because it balances performance and fold stability without relying on raw
untransformed scale.

## OSD-207 Interpretation

OSD-207 remains a feasible fly task, but its condition-stratum heterogeneity is
real.

Mean balanced accuracy by transform/scaling group:

- `none` + `none`: 0.9375
- `none` + `zscore`: 0.86875
- `log1p` + `zscore`: 0.85625
- `log1p` + `none`: 0.825

The best balanced-accuracy variant is `tvg15999_log1p_zscore`
(`balanced_accuracy=0.9375`, `AUROC=0.9609375`,
`condition_stratum_holdout_delta=0.25`). Raw untransformed/unscaled variants
also score high, but they are not promoted as defaults because they may be more
exposed to count-scale effects.

Weakest heldout stratum frequency across the 20 variants:

- w1118_KCNQ370: 18 variants
- Canton.S_Sei.ts1: 2 variants

Decision: keep OSD-207 as a valid fly feasibility task, but mark
`w1118_KCNQ370` as the recurring fragile stratum. The default configuration is
acceptable as a conservative draft baseline, but OSD-207 should not be used as
the cleaner headline multispecies example without a fold-level caveat.

## Default Baseline Decision

Do not chase the best grid score. For both tasks, keep the V9-MULTI-007 default
for now:

- `log1p`
- train-fold `zscore`
- top 2,000 train-fold variable genes

This default remains easy to explain, uses train-only preprocessing, and avoids
over-promoting raw untransformed count scale. The sensitivity grid should be
reported alongside the default result for transparency.

## Next Action

Move to OSD-120 interaction-task design. OSD-120 has adequate source, matrix,
sample-factor, and checksum evidence, but should be modeled as a
genotype/ecotype by light-treatment interaction task rather than as a simple
species-native replicate of OSD-37.

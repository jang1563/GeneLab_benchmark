# Multispecies Baseline Feasibility Review

Status: V9-MULTI-007 review  
Date: 2026-05-23

## Scope

This block tests whether the two species-native multispecies manifests can be
loaded, aligned to local normalized-count matrices, and evaluated with a
conservative draft-only baseline.

Task manifests:

- `draft_osd37_arabidopsis_seedling_spaceflight`
- `draft_osd207_drosophila_whole_body_spaceflight`

Baseline:

- `multispecies_nearest_centroid`
- transform: `log1p`
- scaling: train-fold `zscore`
- feature selection: top 2,000 train-fold variable genes
- fold family: `condition_stratum_candidate_folds`
- claim boundary: `draft_feasibility_only_not_leaderboard`

## Generated Artifacts

- `v9/multispecies/reports/nearest_centroid/multispecies_baseline_summary.csv`
- `v9/multispecies/reports/nearest_centroid/multispecies_baseline_summary.json`
- `v9/multispecies/reports/nearest_centroid/draft_osd37_arabidopsis_seedling_spaceflight/`
- `v9/multispecies/reports/nearest_centroid/draft_osd207_drosophila_whole_body_spaceflight/`

Command:

```bash
/usr/local/bin/python3 scripts/run_v9_multispecies_baseline.py
```

## Summary Metrics

| task_id | n_predictions | balanced_accuracy | AUROC | calibration_error | condition_stratum_holdout_delta |
|---|---:|---:|---:|---:|---:|
| draft_osd37_arabidopsis_seedling_spaceflight | 56 | 0.8214285714 | 0.9196428571 | 0.3066498840 | 0.125 |
| draft_osd207_drosophila_whole_body_spaceflight | 32 | 0.8750000000 | 0.9296875000 | 0.3376823871 | 0.375 |

`mission_discrimination` is skipped for both tasks because each manifest is a
single-source/single-mission-context species-native pilot. This is expected and
is part of the claim boundary.

## Fold-Level Read

OSD-37 is the cleaner first plant feasibility task. All four ecotype holdouts
remain above 0.75 balanced accuracy, and the fold range is 0.125:

- Col.0: 0.8125
- Cvi.0: 0.8333333333
- Ler.0: 0.75
- Ws.2: 0.875

OSD-207 is feasible but more heterogeneous. The weak fold is
`w1118_KCNQ370`, producing a larger fold range of 0.375:

- Canton.S_Sei.ts1: 1.0
- Canton.S_Wild.Type: 0.875
- w1118_KCNQ370: 0.625
- w1118_Wild.Type: 1.0

## Decision

The multispecies species-native manifests are baseline-runnable.

Promote OSD-37 as the cleaner first plant draft feasibility example. Keep
OSD-207 as a valid fly feasibility task, but mark it as needing sensitivity
review because the `w1118_KCNQ370` heldout stratum is weaker than the others.

Do not promote either result as a leaderboard or release claim. The task family
is still draft-only, and the fold design is within-source condition-stratum
robustness rather than mission-held-out generalization.

## Next Action

Run a small sensitivity grid for both species-native tasks:

- `log1p` versus no transform;
- `zscore` versus no scaling;
- top 100, 500, 2,000, 5,000, and all available variable genes where feasible.

The purpose is not to maximize a score. The purpose is to learn whether OSD-37
is stable across preprocessing choices and whether OSD-207's `w1118_KCNQ370`
fold remains the main fragile stratum.

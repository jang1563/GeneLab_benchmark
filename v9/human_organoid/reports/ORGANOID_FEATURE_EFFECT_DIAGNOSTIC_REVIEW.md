# Human Organoid Feature-Effect Diagnostic Review

Status: diagnostic review complete
Date: 2026-05-23
Task: `draft_human_organoid_spaceflight`
Block: `V9-ORG-027`

## Scope

This note reviews the L2 logistic feature-effect pilot:

```text
v9/human_organoid/reports/logistic_feature_effect/
```

The goal is to decide whether classifier-derived feature effects should become
a standard secondary organoid diagnostic, and what must be added before trying
PCA-LR reconstructed weights.

## Reviewed Artifacts

- `v9/human_organoid/reports/logistic_feature_effect/feature_effect.csv.gz`
- `v9/human_organoid/reports/logistic_feature_effect/feature_effect_metadata.json`
- `v9/human_organoid/reports/logistic_feature_effect/metrics.json`
- `v9/human_organoid/reports/logistic_feature_effect/README.md`
- `v9/human_organoid/reports/ORGANOID_CLASSIFIER_FEATURE_EFFECT_CONTRACT.md`

The feature-effect artifact is separate from `response_signature.csv`.
Coefficients are discriminative model effects, not log2FC response signatures.

## Aggregate Comparison

| Report | Direction/Sign Agreement | Rank Correlation | Joined Rows | Scored Significant Rows |
|---|---:|---:|---:|---:|
| Global source-transfer response signature | 0.770673 | 0.176008 | 223,888 | 2,346 |
| Microglia-matched response signature | 0.788913 | 0.150072 | 223,888 | 2,345 |
| Shared-control partial response signature | 0.589942 | 0.035667 | 111,944 | 517 |
| L2 logistic feature effect | 0.607843 | 0.086728 | 16,000 | 204 |

The logistic feature-effect pilot is weaker than the two full-coverage
empirical response-signature diagnostics. It is stronger than the
shared-control partial diagnostic in rank correlation and slightly stronger in
sign agreement, but it is not a replacement for response signatures.

## Per-Contrast Review

| Contrast | Effect Sign | Effect Rank | Top 50 | Top 100 | Top 250 | Top 500 | Global Sign | Global Rank | Microglia Sign | Microglia Rank |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| OSD-863 no known disease, with microglia | 0.636364 | 0.075216 | 1 | 1 | 2 | 6 | 0.673469 | 0.114822 | 0.663265 | 0.089727 |
| OSD-863 no known disease, without microglia | 0.857143 | 0.137029 | 1 | 2 | 2 | 4 | 0.803681 | 0.215367 | 0.785276 | 0.199685 |
| OSD-863 PPMS, with microglia | 0.500000 | 0.082665 | 1 | 1 | 1 | 1 | 0.812500 | 0.207273 | 0.875000 | 0.169350 |
| OSD-863 PPMS, without microglia | 0.666667 | 0.112294 | 1 | 2 | 2 | 3 | 0.868687 | 0.206861 | 0.888889 | 0.183374 |
| OSD-871 no known disease, with microglia | 0.500000 | -0.018114 | 0 | 0 | 0 | 0 | 0.736842 | 0.134653 | 0.710526 | 0.076384 |
| OSD-871 no known disease, without microglia | 0.629630 | 0.084167 | 0 | 0 | 1 | 2 | 0.651584 | 0.071176 | 0.663636 | 0.036763 |
| OSD-871 sporadic Parkinson disease, with microglia | 0.500000 | 0.127598 | 2 | 4 | 5 | 6 | 0.738462 | 0.204739 | 0.726923 | 0.158517 |
| OSD-871 sporadic Parkinson disease, without microglia | 0.601852 | 0.133048 | 6 | 8 | 16 | 29 | 0.791178 | 0.267495 | 0.822192 | 0.291061 |

The best effect-sign contrast is OSD-863 no-known-disease without microglia
(`0.857143`), but its rank correlation remains below the global empirical
response signature. OSD-871 no-known-disease with microglia is essentially a
failure case for the feature-effect pilot.

## Top-K Interpretation

The aggregate top-k DE overlap is low in absolute terms:

| K | Observed Overlap | Random Expectation | Enrichment |
|---:|---:|---:|---:|
| 50 | 1 | 0.6375 | 1.57 |
| 100 | 5 | 1.2750 | 3.92 |
| 250 | 10 | 3.1875 | 3.14 |
| 500 | 14 | 6.3750 | 2.20 |

This suggested that classifier effects may weakly enrich for DE genes near the
top of the ranked list. At the time of this review, the report did not yet
include permutation or hypergeometric calibration, so this remained a
qualitative observation.

Follow-up calibration was added in
`v9/human_organoid/reports/ORGANOID_FEATURE_EFFECT_NULL_CALIBRATION_REVIEW.md`.
That review keeps K=50 as null-like, but supports aggregate top-k enrichment at
K=100, K=250, and K=500 under a hypergeometric random-top-k null.

## Decision

Keep `feature_effect.csv` as an optional secondary diagnostic artifact.

Do not make it part of the primary organoid task or leaderboard. Do not treat it
as a response-signature replacement.

Recommended interpretation:

> L2 logistic gene-space coefficients provide a weak model-effect diagnostic.
> They preserve some DE sign/rank signal and show modest top-k enrichment, but
> empirical source-transfer response signatures remain the stronger diagnostic
> surface for this small organoid task.

## Next Step

Before PCA-LR reconstructed weights, add top-k/null calibration for
feature-effect diagnostics.

Reason:

- top-k overlap is the one place where the feature-effect pilot hints at useful
  signal;
- raw overlap fractions are hard to interpret when only 204 of 16,000 joined
  feature/contrast rows are significant;
- PCA-LR reconstructed weights will add more modeling assumptions, so the
  scorer should first learn to report enrichment relative to a clear null.

## Recommended Next Block

`V9-ORG-028: Human organoid feature-effect top-k null calibration`

Expected work:

- add hypergeometric or permutation-style expected-overlap diagnostics for
  `feature_effect_top_k_de_overlap`;
- report expected overlap, enrichment, and a lightweight p-value or empirical
  percentile where feasible;
- keep calibration diagnostic-only and non-leaderboard;
- regenerate `v9/human_organoid/reports/logistic_feature_effect/`.

## Final Decision

Feature effects are worth keeping as a secondary diagnostic, but only with
explicit artifact separation and better top-k calibration. PCA-LR reconstructed
weights should remain deferred until calibrated L2 logistic feature-effect
diagnostics are stable.

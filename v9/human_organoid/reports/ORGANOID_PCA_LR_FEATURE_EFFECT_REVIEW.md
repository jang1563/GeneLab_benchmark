# Human Organoid PCA-LR Feature-Effect Review

Status: pilot review complete
Date: 2026-05-23
Task: `draft_human_organoid_spaceflight`
Block: `V9-ORG-030`

## Scope

This note reviews the PCA-LR reconstructed feature-effect pilot:

```text
v9/human_organoid/reports/pca_lr_feature_effect/
```

The pilot implements the V9-ORG-029 design. It reconstructs PCA-LR PC-space
classifier coefficients into train-standardized gene space with:

```text
pca.components_.T @ logistic.coef_[0]
```

The artifact remains a diagnostic `feature_effect.csv`, not a
`response_signature.csv`.

## Generated Artifacts

- `v9/human_organoid/reports/pca_lr_feature_effect/predictions.csv`
- `v9/human_organoid/reports/pca_lr_feature_effect/feature_effect.csv.gz`
- `v9/human_organoid/reports/pca_lr_feature_effect/feature_effect_metadata.json`
- `v9/human_organoid/reports/pca_lr_feature_effect/metrics.json`
- `v9/human_organoid/reports/pca_lr_feature_effect/run_manifest.json`
- `v9/human_organoid/reports/pca_lr_feature_effect/README.md`

The report emits 16,000 feature-effect rows. PCA retains all train-fold
available rank in the default run:

| Target Source | Training Source | PCA Components Used | Explained Variance Ratio Sum |
|---|---|---:|---:|
| OSD-863 | OSD-871 | 22 | 1.000000 |
| OSD-871 | OSD-863 | 18 | 1.000000 |

## Aggregate Comparison

| Metric | PCA-LR Reconstructed | L2 Logistic | Delta |
|---|---:|---:|---:|
| `feature_effect_direction_match` | 0.5980392156862745 | 0.6078431372549019 | -0.009804 |
| `feature_effect_rank_correlation` | 0.08668664748698189 | 0.08672800238082004 | -0.000041 |

The PCA-LR reconstruction is essentially tied with L2 logistic in rank
correlation and slightly weaker in direction agreement.

## Calibrated Top-K Comparison

| K | PCA-LR Observed | L2 Observed | Expected | Enrichment | Hypergeometric p>=x |
|---:|---:|---:|---:|---:|---:|
| 50 | 1 | 1 | 0.6375 | 1.5686 | 0.474071 |
| 100 | 5 | 5 | 1.2750 | 3.9216 | 0.009090 |
| 250 | 10 | 10 | 3.1875 | 3.1373 | 0.001414 |
| 500 | 14 | 14 | 6.3750 | 2.1961 | 0.004966 |

The aggregate top-k diagnostic is identical to the L2 logistic report. PCA-LR
therefore does not add a new top-ranked DE enrichment signal.

## Per-Contrast Comparison

| Contrast | PCA Direction | L2 Direction | Direction Delta | PCA Rank | L2 Rank | Rank Delta | PCA Top500 | L2 Top500 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| OSD-863 no known disease, with microglia | 0.636364 | 0.636364 | +0.000000 | 0.075315 | 0.075216 | +0.000098 | 6 | 6 |
| OSD-863 no known disease, without microglia | 0.857143 | 0.857143 | +0.000000 | 0.137565 | 0.137029 | +0.000536 | 3 | 4 |
| OSD-863 PPMS, with microglia | 0.500000 | 0.500000 | +0.000000 | 0.082737 | 0.082665 | +0.000072 | 1 | 1 |
| OSD-863 PPMS, without microglia | 0.666667 | 0.666667 | +0.000000 | 0.112742 | 0.112294 | +0.000448 | 3 | 3 |
| OSD-871 no known disease, with microglia | 0.500000 | 0.500000 | +0.000000 | -0.017486 | -0.018114 | +0.000627 | 0 | 0 |
| OSD-871 no known disease, without microglia | 0.629630 | 0.629630 | +0.000000 | 0.083006 | 0.084167 | -0.001161 | 2 | 2 |
| OSD-871 Parkinson disease, with microglia | 0.466667 | 0.500000 | -0.033333 | 0.126830 | 0.127598 | -0.000769 | 6 | 6 |
| OSD-871 Parkinson disease, without microglia | 0.592593 | 0.601852 | -0.009259 | 0.132900 | 0.133048 | -0.000147 | 29 | 29 |

The only meaningful direction differences are negative and occur in the OSD-871
Parkinson disease contrasts. Small rank gains in OSD-863 contrasts are too
small to justify a stronger claim.

## Interpretation

The reconstruction is mathematically valid and implementation-ready, but it
does not improve the feature-effect diagnostic surface.

Because PCA retains all available train-fold rank in this small setting, the
PCA-LR reconstruction behaves like a rotated full-rank linear classifier rather
than a genuinely compressed representation. The added interpretation layer
therefore increases complexity without adding diagnostic signal.

## Decision

`V9-ORG-030` is complete.

Keep the PCA-LR reconstructed report as an optional secondary diagnostic and
negative comparison.

Do not promote PCA-LR reconstructed coefficients over the simpler L2 logistic
gene-space feature-effect report.

Do not pursue more PCA-LR reconstructed-weight variants unless a future task
uses a genuinely low-rank PCA setting or a larger training sample regime where
PCA compression is substantive.

## Next Block

`V9-ORG-031: Human organoid diagnostic consolidation and release boundary`

Expected work:

- consolidate response-signature and feature-effect diagnostic decisions;
- update v9 README/report pointers for the completed organoid diagnostic family;
- define which organoid artifacts are part of the draft v9 alpha surface;
- decide whether the active lane should return to multispecies expansion or
  release packaging.

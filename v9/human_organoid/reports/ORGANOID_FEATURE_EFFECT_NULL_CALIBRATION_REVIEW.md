# Human Organoid Feature-Effect Null Calibration Review

Status: implementation review complete
Date: 2026-05-23
Task: `draft_human_organoid_spaceflight`
Block: `V9-ORG-028`

## Scope

This note records the top-k null calibration added to the L2 logistic
feature-effect diagnostic:

```text
v9/human_organoid/reports/logistic_feature_effect/
```

The calibration applies only to `feature_effect_top_k_de_overlap`. It does not
make feature effects a primary benchmark metric, and it does not make
classifier coefficients into `response_signature.csv` rows.

## Implemented Null

For each top-k entry, the scorer now reports:

- `n_feature_universe`: joined feature/contrast rows available for scoring;
- `n_selected`: top absolute feature-effect rows selected;
- `n_significant_reference_rows`: FDR <= 0.05 DE reference rows in that
  universe;
- `n_overlap`: selected rows that are significant in the DE reference;
- `expected_overlap`: `n_selected * n_significant_reference_rows /
  n_feature_universe`;
- `enrichment`: `n_overlap / expected_overlap`;
- `hypergeometric_p_value_greater_equal`: exact upper-tail probability under a
  random top-k draw without replacement;
- `null_model`: `hypergeometric_random_top_k_without_replacement`.

The null is intentionally simple. It asks whether a top-k list contains more DE
genes than expected if the same number of rows were selected randomly from the
same scored feature universe.

## Regenerated Report

The L2 logistic feature-effect report was regenerated after adding calibrated
top-k fields:

- `v9/human_organoid/reports/logistic_feature_effect/feature_effect.csv.gz`
- `v9/human_organoid/reports/logistic_feature_effect/metrics.json`
- `v9/human_organoid/reports/logistic_feature_effect/README.md`
- `v9/human_organoid/reports/logistic_feature_effect/run_manifest.json`

The direction and rank diagnostics are unchanged:

| Metric | Value |
|---|---:|
| `feature_effect_direction_match` | 0.6078431372549019 |
| `feature_effect_rank_correlation` | 0.08672800238082004 |

## Aggregate Calibrated Top-K

| K | Observed | Expected | Enrichment | Hypergeometric p>=x |
|---:|---:|---:|---:|---:|
| 50 | 1 | 0.6375 | 1.5686 | 0.474071 |
| 100 | 5 | 1.2750 | 3.9216 | 0.009090 |
| 250 | 10 | 3.1875 | 3.1373 | 0.001414 |
| 500 | 14 | 6.3750 | 2.1961 | 0.004966 |

The calibrated aggregate result supports a narrow claim: the top-ranked L2
logistic coefficients are enriched for DE genes above a random top-k null at
K=100, K=250, and K=500. K=50 is not distinguishable from the random-top-k null.

## Per-Contrast Pattern

The aggregate signal is heterogeneous.

Contrasts with nominal calibrated top-k evidence include:

| Contrast | K | Observed | Expected | Enrichment | Hypergeometric p>=x |
|---|---:|---:|---:|---:|---:|
| OSD-863 no known disease, with microglia | 500 | 6 | 2.7500 | 2.1818 | 0.033926 |
| OSD-863 PPMS, with microglia | 50 | 1 | 0.0500 | 20.0000 | 0.049387 |
| OSD-863 PPMS, without microglia | 100 | 2 | 0.3000 | 6.6667 | 0.032544 |
| OSD-871 Parkinson disease, without microglia | 50 | 6 | 2.7000 | 2.2222 | 0.049354 |

The OSD-871 no-known-disease contrasts remain weak or negative under
calibration. OSD-871 no-known with microglia has zero overlap through K=500,
and OSD-871 no-known without microglia is below random expectation at K=250 and
K=500.

## Decision

`V9-ORG-028` is complete.

Keep calibrated top-k fields in `feature_effect_top_k_de_overlap`.

Interpret the L2 logistic feature-effect artifact as:

> a secondary, model-effect diagnostic with calibrated aggregate top-k DE
> enrichment, but weak and heterogeneous per-contrast behavior.

Do not promote feature-effect diagnostics to primary benchmark metrics. Do not
move classifier coefficients into response-signature artifacts.

## Next Block

`V9-ORG-029: Human organoid PCA-LR reconstructed feature-effect design`

Before implementing PCA-LR reconstructed gene weights, write a design review
that fixes:

- the exact reconstruction formula from PCA components and logistic
  coefficients;
- the train-fold-only fit boundary for scaling, PCA, and logistic regression;
- the `feature_effect.csv` schema reuse rules;
- comparison requirements against the calibrated L2 logistic feature-effect
  report;
- a go/no-go rule for whether reconstructed weights are worth implementing.

# Human Organoid Model-Produced Response-Signature Adapter Design

Status: design decision
Date: 2026-05-22
Task: `draft_human_organoid_spaceflight`
Block: `V9-ORG-017`

## Decision

The first real `response_signature.csv` adapter should be a conservative
train-fold-only source-transfer empirical signature baseline.

It should not use:

- the derived DE reference table;
- the OSDR differential-expression table;
- any target-source expression rows for the target contrast;
- all-sample source-local LEO/Ground mean differences.

It should use only normalized expression matrices and sample-factor labels from
the training side of the existing organoid-type/source holdout folds.

## Why This Is The First Safe Adapter

The current draft task has 42 samples:

- OSD-863: 19 cortical neural organoid samples.
- OSD-871: 23 dopaminergic neural organoid samples.

Source and organoid type are perfectly coupled:

| Source | Organoid Type | Samples |
|---|---|---:|
| OSD-863 | cortical_neural_organoid | 19 |
| OSD-871 | dopaminergic_neural_organoid | 23 |

Disease context and donor/iPSC-line are also coupled with source:

| Source | Disease Contexts | Donor/iPSC-line Blocks |
|---|---|---|
| OSD-863 | no_known_diseases; primary_progressive_multiple_sclerosis | Subject1; Subject2 |
| OSD-871 | no_known_diseases; sporadic_parkinson_disease | Subject3; Subject4 |

This means a source-local empirical signature for an OSD-863 contrast would be
almost the same information source as the OSD-863 DE reference being scored. It
would be a scorer sanity check, not a model-produced benchmark prediction.

The safe first adapter should therefore predict each target source from the
other source:

| Target Contrasts | Training Fold | Training Source |
|---|---|---|
| all OSD-863 contrast ids | `holdout_organoid_type_cortical_neural_organoid` | OSD-871 |
| all OSD-871 contrast ids | `holdout_organoid_type_dopaminergic_neural_organoid` | OSD-863 |

This is biologically crude, but it is honest. It tests whether a generic
spaceflight response learned from one organoid fate/source transfers to the
other source-specific reference contrasts.

## Current Contrast Feasibility

All eight DE-reference contrasts have both Ground and LEO/ISS samples in the
source-local sample table:

| Source | Disease | Microglia | Ground | LEO/ISS |
|---|---|---|---:|---:|
| OSD-863 | no_known_diseases | with_microglia | 3 | 3 |
| OSD-863 | no_known_diseases | without_microglia | 3 | 2 |
| OSD-863 | primary_progressive_multiple_sclerosis | with_microglia | 2 | 2 |
| OSD-863 | primary_progressive_multiple_sclerosis | without_microglia | 2 | 2 |
| OSD-871 | no_known_diseases | with_microglia | 3 | 2 |
| OSD-871 | no_known_diseases | without_microglia | 2 | 3 |
| OSD-871 | sporadic_parkinson_disease | with_microglia | 3 | 2 |
| OSD-871 | sporadic_parkinson_disease | without_microglia | 4 | 4 |

These counts are enough to define source-local DE references, but too small and
too confounded to make source-local empirical signatures a legitimate model
prediction.

## Rejected Adapter Options

### All-Sample Empirical Signature

Rejected for benchmark use.

Computing `mean(LEO/ISS) - mean(Ground)` from all 42 samples and scoring it
against DE references built from the same public data would mix target and test
information into the prediction.

### Target-Source Empirical Signature

Rejected for benchmark use.

Computing source-local signatures from OSD-863 samples and scoring against the
OSD-863 DE reference is too close to reconstructing the reference from the
scored samples.

### Mirrored Reference Fixture

Kept only as a smoke test.

`v9/human_organoid/reports/response_signature_smoke/response_signature.csv`
mirrors selected reference rows. It validates scorer plumbing, but must never be
reported as model performance.

### Logistic Coefficient Signature

Deferred.

A linear classifier can produce feature coefficients, but those coefficients are
not log2 fold changes. They may become a useful ranking diagnostic after the
empirical source-transfer baseline exists, but they should not be the first
adapter.

## Proposed Adapter Contract

Name:

```text
organoid_source_transfer_empirical_signature
```

Input:

- `load_human_organoid_task()` feature matrix.
- `sample_factors.draft.csv`.
- organoid-type/source holdout folds from the task manifest.

Training rule:

1. For each target source, select the fold that holds out that source's
   organoid type.
2. Use only the training samples from the other source.
3. Compute one global training response signature per gene:

```text
log2(mean_normalized_expression_train_LEO_or_ISS + 1)
- log2(mean_normalized_expression_train_Ground + 1)
```

4. Emit that same training-derived signature for every target-source
   `contrast_id`.

Output:

```text
response_signature.csv
```

Required columns from the scorer contract:

```text
task_id
source_id
contrast_id
gene_symbol
ensembl_id
predicted_log2fc_leo_or_iss_minus_ground
```

Recommended optional columns:

```text
fold_id
signature_model_id
training_source_id
target_source_id
training_scope
n_train_ground
n_train_leo_or_iss
signature_generation_method
reference_usage_policy
```

Reference usage policy should be:

```text
reference_not_used_for_signature_generation
```

## Expected Limitations

- The predicted signature is global within a training source, not disease- or
  microglia-specific.
- OSD-863 and OSD-871 differ by organoid fate, disease blocks, donor/iPSC-line,
  and source-specific processing; poor signature agreement would be
  interpretable but not surprising.
- The metric remains diagnostic-only because the task is not mission-held-out
  and donor/source/disease factors are coupled.
- The output should be reported as a baseline response-signature diagnostic, not
  as a leaderboard result.

## V9-ORG-018 Implementation Plan

Files:

- `spacebio_bench/response_signature_adapters.py`
- `scripts/run_v9_human_organoid_source_transfer_signature.py`
- `tests/test_v9_spacebio_bench.py`
- `v9/human_organoid/reports/source_transfer_signature/`

Required tests:

- the adapter writes one response row per target contrast/gene pair;
- target-source samples are excluded from signature generation;
- `reference_not_used_for_signature_generation` appears in output metadata;
- evaluator computes DE/signature diagnostics from the generated artifact;
- run manifest records `predictions.csv`, `response_signature.csv`,
  `task_manifest`, and `reference_signature_table`.

Done when:

- A source-transfer empirical response-signature report exists under
  `v9/human_organoid/reports/source_transfer_signature/`.
- The report includes computed diagnostic signature metrics and an explicit
  non-leaderboard claim boundary.

## V9-ORG-018 Output

The source-transfer adapter is implemented and the diagnostic report exists at:

```text
v9/human_organoid/reports/source_transfer_signature/
```

The report writes compressed `response_signature.csv.gz` with 223,888 rows:
27,986 common Ensembl features times eight target-source contrasts. It records
`reference_not_used_for_signature_generation`, uses OSD-871 training samples for
OSD-863 target contrasts and OSD-863 training samples for OSD-871 target
contrasts, and keeps the output non-leaderboard.

Observed diagnostic scores in this draft run:

- `de_direction_match`: 0.7706734867860188
- `signature_rank_correlation`: 0.1760078660242601

These values are source-transfer diagnostics only; they are not leaderboard
claims.

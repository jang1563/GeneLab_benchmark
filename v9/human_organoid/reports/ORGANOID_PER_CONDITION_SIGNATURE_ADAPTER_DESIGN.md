# Human Organoid Per-Condition Source-Transfer Adapter Design

Status: design decision
Date: 2026-05-22
Task: `draft_human_organoid_spaceflight`
Block: `V9-ORG-020`

## Decision

The next response-signature adapter should be:

```text
organoid_microglia_matched_source_transfer_empirical_signature
```

It should match `microglia_condition` within the training source, while
preserving the same source-transfer leakage boundary as the current global
source-transfer baseline.

This should be implemented as a comparative diagnostic against:

```text
organoid_source_transfer_empirical_signature
```

It should not replace the global source-transfer baseline as the first
conservative diagnostic. It should test whether conditioning on microglia status
improves target-source DE direction and rank agreement.

## Why Microglia-Matched Comes First

The current source-transfer baseline computes one global training-source
LEO/ISS-minus-Ground signature and emits that same signature for every target
contrast in the held-out source. This is safe but crude because each target
contrast is disease- and microglia-specific.

The public sample-factor table supports a stricter condition-aware adapter for
microglia status:

| Source | Microglia | Ground | LEO/ISS |
|---|---|---:|---:|
| OSD-863 | with_microglia | 5 | 5 |
| OSD-863 | without_microglia | 5 | 4 |
| OSD-871 | with_microglia | 6 | 4 |
| OSD-871 | without_microglia | 6 | 7 |

Every target contrast has an opposite-source training stratum with both Ground
and LEO/ISS samples when matching only `microglia_condition`.

Disease matching is not fully crossed:

| Source | Disease Context | Ground | LEO/ISS |
|---|---|---:|---:|
| OSD-863 | no_known_diseases | 6 | 5 |
| OSD-863 | primary_progressive_multiple_sclerosis | 4 | 4 |
| OSD-871 | no_known_diseases | 5 | 5 |
| OSD-871 | sporadic_parkinson_disease | 7 | 6 |

Only `no_known_diseases` is shared across both sources. PPMS exists only in
OSD-863 and sporadic Parkinson disease exists only in OSD-871. A disease-matched
adapter would therefore score only the healthy/control contrasts unless it used
a fallback. That partial coverage is useful later, but it should not be the next
primary comparative adapter.

## Opposite-Source Availability

| Target Source | Training Source | Target Disease | Target Microglia | Microglia-Matched Available | Disease-Matched Available | Disease+Microglia Available |
|---|---|---|---|---|---|---|
| OSD-863 | OSD-871 | no_known_diseases | with_microglia | yes, 6/4 | yes, 5/5 | yes, 3/2 |
| OSD-863 | OSD-871 | no_known_diseases | without_microglia | yes, 6/7 | yes, 5/5 | yes, 2/3 |
| OSD-863 | OSD-871 | primary_progressive_multiple_sclerosis | with_microglia | yes, 6/4 | no | no |
| OSD-863 | OSD-871 | primary_progressive_multiple_sclerosis | without_microglia | yes, 6/7 | no | no |
| OSD-871 | OSD-863 | no_known_diseases | with_microglia | yes, 5/5 | yes, 6/5 | yes, 3/3 |
| OSD-871 | OSD-863 | no_known_diseases | without_microglia | yes, 5/4 | yes, 6/5 | yes, 3/2 |
| OSD-871 | OSD-863 | sporadic_parkinson_disease | with_microglia | yes, 5/5 | no | no |
| OSD-871 | OSD-863 | sporadic_parkinson_disease | without_microglia | yes, 5/4 | no | no |

Counts are shown as Ground/LEO_or_ISS in the opposite-source training stratum.

## Adapter Contract

Input:

- `load_human_organoid_task()` feature matrix.
- `sample_factors.draft.csv`.
- derived DE-reference manifest for source-specific target contrast ids.
- organoid-type/source holdout folds from the task manifest.

Training rule for each target contrast:

1. Select the source-transfer fold where the target source is held out.
2. Restrict training samples to the opposite source.
3. Restrict those training samples to the target contrast's
   `microglia_condition`.
4. Require at least one Ground and one LEO/ISS training sample in that stratum.
5. Compute one stratum-specific response signature per gene:

```text
log2(mean_normalized_expression_train_LEO_or_ISS_in_microglia_stratum + 1)
- log2(mean_normalized_expression_train_Ground_in_microglia_stratum + 1)
```

6. Emit that signature for the target contrast only.

Expected full-report coverage:

- 8 target contrasts.
- 27,986 common Ensembl features.
- 223,888 response rows, matching the global source-transfer report coverage.

## Required Output Columns

Keep the scorer-required columns:

```text
task_id
source_id
contrast_id
gene_symbol
ensembl_id
predicted_log2fc_leo_or_iss_minus_ground
```

Recommended provenance columns:

```text
fold_id
signature_model_id
training_source_id
target_source_id
training_scope
conditioning_strategy
conditioning_factor
conditioning_value
n_train_ground
n_train_leo_or_iss
n_condition_train_ground
n_condition_train_leo_or_iss
signature_generation_method
reference_usage_policy
```

The `reference_usage_policy` must remain:

```text
reference_not_used_for_signature_generation
```

## Rejected Primary Designs

### Disease-Matched Primary Adapter

Rejected as the next primary comparative adapter.

The disease context is not independently crossed between OSD-863 and OSD-871.
Only `no_known_diseases` exists in both sources. A disease-matched report would
either cover only four target contrasts or silently fall back for disease-only
contexts. That makes aggregate metrics harder to compare with the current
global source-transfer report.

### Disease+Microglia Matched Primary Adapter

Rejected as the next primary comparative adapter for the same reason.

It is biologically attractive but coverage is limited to the shared
`no_known_diseases` target contrasts.

### Implicit Fallback Adapter

Rejected for the first condition-aware report.

If a disease stratum is unavailable, silently falling back to global or
microglia-only signatures would mix adapter strategies inside one aggregate
score. Fallbacks are acceptable only if every row records the fallback reason
and the report separates metrics by `conditioning_strategy`.

### Target-Source Condition Signature

Rejected for benchmark scoring.

Computing a condition-specific signature from the target source would use the
same samples that define the scored reference contrasts and would collapse back
into reference reconstruction.

## Evaluation Plan

The microglia-matched report should be evaluated by the existing
`de_direction_match` and `signature_rank_correlation` scorer.

The report should include a small comparison table against the global
source-transfer baseline:

| Metric | Global Source-Transfer | Microglia-Matched Source-Transfer | Delta |
|---|---:|---:|---:|
| `de_direction_match` | from global report | from new report | new - global |
| `signature_rank_correlation` | from global report | from new report | new - global |

It should also compare per-contrast deltas, because an aggregate improvement
could hide failure in a specific disease or microglia block.

## Implementation Plan

Recommended block:

```text
V9-ORG-021: Human organoid microglia-matched source-transfer adapter
```

Likely files:

- `spacebio_bench/response_signature_adapters.py`
- `scripts/run_v9_human_organoid_microglia_source_transfer_signature.py`
- `tests/test_v9_spacebio_bench.py`
- `v9/human_organoid/reports/microglia_source_transfer_signature/`

Implementation notes:

- Reuse the existing source-transfer fold discovery.
- Add a general condition-filtered signature helper rather than duplicating the
  whole global adapter.
- Keep gene IDs Ensembl-first, matching the existing scorer fallback path.
- Add metadata for condition coverage and no-reference generation policy.
- Keep the classification `predictions.csv` fixture explicitly non-performance
  if the evaluator still requires a prediction CSV for report writing.

Required tests:

- target-source samples are excluded from signature generation;
- every emitted row has the target contrast's microglia condition in metadata;
- each emitted contrast has both Ground and LEO/ISS condition-stratum training
  counts;
- metadata records `reference_not_used_for_signature_generation`;
- the evaluator computes DE/signature diagnostics against a small reference
  fixture;
- the report README states that classification metrics are plumbing-only and
  response-signature metrics are diagnostic-only.

## Deferred Follow-Up

After the microglia-matched adapter is implemented and compared, revisit a
partial disease/disease+microglia diagnostic. That report should be explicitly
partial coverage and should write skipped-contrast metadata rather than using an
implicit fallback.

## Done Criteria

This design is complete when the backlog points to V9-ORG-021 as the active
implementation block and the next adapter has a clear leakage boundary, output
contract, tests, and report destination.

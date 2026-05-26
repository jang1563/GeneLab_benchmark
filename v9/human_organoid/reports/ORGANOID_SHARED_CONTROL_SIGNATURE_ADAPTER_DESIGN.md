# Human Organoid Shared-Control Source-Transfer Adapter Design

Status: design decision
Date: 2026-05-22
Task: `draft_human_organoid_spaceflight`
Block: `V9-ORG-023`

## Decision

The next optional response-signature adapter should be a partial-coverage
shared-control disease+microglia matched source-transfer diagnostic:

```text
organoid_shared_control_disease_microglia_source_transfer_empirical_signature
```

It should emit signatures only for target contrasts where the opposite-source
training fold contains matching:

- `disease_context = no_known_diseases`
- `microglia_condition = target microglia condition`
- at least one Ground and one LEO/ISS sample

It must explicitly skip disease-specific target contrasts whose disease context
is not present in the opposite source.

## Why Partial Coverage Is Required

Disease context is not fully crossed across the two public organoid sources:

| Source | Disease Contexts |
|---|---|
| OSD-863 | `no_known_diseases`; `primary_progressive_multiple_sclerosis` |
| OSD-871 | `no_known_diseases`; `sporadic_parkinson_disease` |

Only `no_known_diseases` exists in both sources. Therefore a disease-matched
adapter cannot cover all eight target contrasts without using an implicit
fallback. This design forbids silent fallback because mixed-strategy aggregate
metrics would be hard to interpret.

## Eligible Target Contrasts

The shared-control adapter should emit rows for these four target contrasts:

| Target Contrast | Target Source | Training Source | Training Stratum | Ground | LEO/ISS |
|---|---|---|---|---:|---:|
| `osd_863_no_known_diseases_with_microglia_leo_or_iss_vs_ground` | OSD-863 | OSD-871 | no_known_diseases + with_microglia | 3 | 2 |
| `osd_863_no_known_diseases_without_microglia_leo_or_iss_vs_ground` | OSD-863 | OSD-871 | no_known_diseases + without_microglia | 2 | 3 |
| `osd_871_no_known_diseases_with_microglia_leo_or_iss_vs_ground` | OSD-871 | OSD-863 | no_known_diseases + with_microglia | 3 | 3 |
| `osd_871_no_known_diseases_without_microglia_leo_or_iss_vs_ground` | OSD-871 | OSD-863 | no_known_diseases + without_microglia | 3 | 2 |

Expected emitted rows:

- 4 target contrasts.
- 27,986 common Ensembl features.
- 111,944 response-signature rows.

## Skipped Target Contrasts

The adapter should not emit rows for these target contrasts:

| Target Contrast | Target Source | Missing In Training Source | Skip Reason |
|---|---|---|---|
| `osd_863_primary_progressive_multiple_sclerosis_with_microglia_leo_or_iss_vs_ground` | OSD-863 | OSD-871 | PPMS absent in OSD-871 |
| `osd_863_primary_progressive_multiple_sclerosis_without_microglia_leo_or_iss_vs_ground` | OSD-863 | OSD-871 | PPMS absent in OSD-871 |
| `osd_871_sporadic_parkinson_disease_with_microglia_leo_or_iss_vs_ground` | OSD-871 | OSD-863 | sporadic Parkinson disease absent in OSD-863 |
| `osd_871_sporadic_parkinson_disease_without_microglia_leo_or_iss_vs_ground` | OSD-871 | OSD-863 | sporadic Parkinson disease absent in OSD-863 |

Skipped contrasts should be recorded in `response_signature_metadata.json` with
`status=skipped_missing_shared_disease_context`.

## Adapter Contract

Training rule for each target contrast:

1. Select the source-transfer fold where the target source is held out.
2. Restrict training samples to the opposite source.
3. Restrict training samples to the target disease context and target microglia
   condition.
4. If both Ground and LEO/ISS samples exist, compute:

```text
log2(mean_normalized_expression_train_LEO_or_ISS_matched_stratum + 1)
- log2(mean_normalized_expression_train_Ground_matched_stratum + 1)
```

5. If either label is missing, emit no response rows for that contrast and
   record an explicit skipped-contrast metadata row.

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
conditioning_factors
target_disease_context
target_microglia_condition
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

## Metric Interpretation

This adapter must be scored as partial coverage:

- aggregate `de_direction_match` is over emitted shared-control contrasts only;
- aggregate `signature_rank_correlation` is over emitted shared-control
  contrasts only;
- these values are not directly comparable to the all-eight-contrast global or
  microglia-matched aggregate scores;
- per-contrast deltas against global and microglia-matched reports are the main
  interpretation surface.

The README should include a coverage table:

| Coverage | Count |
|---|---:|
| emitted target contrasts | 4 |
| skipped target contrasts | 4 |
| emitted response rows | 111,944 |

## Rejected Options

### Silent Fallback

Rejected.

Falling back from disease+microglia matching to microglia-only or global
signatures would create mixed semantics inside one aggregate score.

### Full-Coverage Disease-Matched Claim

Rejected.

The public organoid disease contexts are source-specific except for
`no_known_diseases`, so a full-coverage disease-matched claim is not supported
by the current data.

### Classifier Coefficients Before Partial Shared-Control Test

Deferred.

Classifier coefficients are useful, but they are not log2 fold changes. The
shared-control adapter stays within the current response-signature log2FC
contract and is the cleaner next diagnostic.

## Implementation Plan

Recommended block:

```text
V9-ORG-024: Human organoid shared-control source-transfer adapter
```

Likely files:

- `spacebio_bench/response_signature_adapters.py`
- `scripts/run_v9_human_organoid_shared_control_source_transfer_signature.py`
- `tests/test_v9_spacebio_bench.py`
- `v9/human_organoid/reports/shared_control_source_transfer_signature/`

Implementation notes:

- Reuse the condition-filtered helper added for the microglia-matched adapter.
- Add disease+microglia filtering without changing scorer contracts.
- Preserve skipped-contrast metadata.
- Keep the report README explicit that classification metrics are plumbing-only
  and signature metrics are partial-coverage diagnostics.

Required tests:

- emits rows only for shared `no_known_diseases` target contrasts;
- records four skipped disease-specific contrasts;
- excludes target-source samples from signature generation;
- records `reference_not_used_for_signature_generation`;
- evaluator computes diagnostics on a small partial reference fixture;
- README contains partial-coverage language.

## Done Criteria

The design is complete when the backlog points to V9-ORG-024 as the active
implementation block and the partial-coverage interpretation is explicit.

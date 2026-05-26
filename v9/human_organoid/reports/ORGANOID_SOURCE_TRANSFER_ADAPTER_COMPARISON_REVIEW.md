# Human Organoid Source-Transfer Adapter Comparison Review

Status: diagnostic comparison complete
Date: 2026-05-22
Task: `draft_human_organoid_spaceflight`
Block: `V9-ORG-022`

## Scope

This note compares the two current human organoid response-signature
diagnostics:

- global source-transfer:
  `v9/human_organoid/reports/source_transfer_signature/`
- microglia-matched source-transfer:
  `v9/human_organoid/reports/microglia_source_transfer_signature/`

Both reports keep the same leakage boundary:

- target-source expression rows are not used for signature generation;
- derived DE-reference rows are not used for signature generation;
- DE-reference scoring is post hoc only;
- outputs are diagnostic-only and non-leaderboard.

## Aggregate Comparison

| Metric | Global Source-Transfer | Microglia-Matched | Delta |
|---|---:|---:|---:|
| `de_direction_match` | 0.770673486786 | 0.788912579957 | +0.018239093171 |
| `signature_rank_correlation` | 0.176007866024 | 0.150072282946 | -0.025935583078 |

The microglia-matched adapter improves sign-level direction agreement but
worsens all-gene rank correlation. This is a useful condition-sensitivity
result, not evidence that the microglia-matched adapter should replace the
global source-transfer diagnostic.

## Per-Contrast Delta

| Contrast | Significant Rows | Direction Delta | Rank Delta | Microglia Direction | Microglia Rank |
|---|---:|---:|---:|---:|---:|
| OSD-863 no known disease, with microglia | 98 | -0.010204 | -0.025095 | 0.663265 | 0.089727 |
| OSD-863 no known disease, without microglia | 163 | -0.018405 | -0.015682 | 0.785276 | 0.199685 |
| OSD-863 PPMS, with microglia | 16 | +0.062500 | -0.037923 | 0.875000 | 0.169350 |
| OSD-863 PPMS, without microglia | 99 | +0.020202 | -0.023487 | 0.888889 | 0.183374 |
| OSD-871 no known disease, with microglia | 38 | -0.026316 | -0.058269 | 0.710526 | 0.076384 |
| OSD-871 no known disease, without microglia | 220 | +0.012053 | -0.034413 | 0.663636 | 0.036763 |
| OSD-871 sporadic Parkinson disease, with microglia | 260 | -0.011538 | -0.046222 | 0.726923 | 0.158517 |
| OSD-871 sporadic Parkinson disease, without microglia | 1,451 | +0.031013 | +0.023566 | 0.822192 | 0.291061 |

Counts:

- direction improved in 4 of 8 contrasts;
- rank correlation improved in 1 of 8 contrasts;
- both direction and rank improved in 1 of 8 contrasts;
- both direction and rank worsened in 4 of 8 contrasts.

## Group-Level Interpretation

| Group | Significant Rows | Weighted Direction Delta | Mean Rank Delta |
|---|---:|---:|---:|
| OSD-863 targets | 376 | -0.002660 | -0.025546 |
| OSD-871 targets | 1,969 | +0.022169 | -0.028835 |
| with microglia | 412 | -0.009709 | -0.041877 |
| without microglia | 1,933 | +0.024134 | -0.012504 |
| no known diseases | 519 | -0.004525 | -0.033365 |
| PPMS | 115 | +0.026087 | -0.030705 |
| sporadic Parkinson disease | 1,711 | +0.024547 | -0.011328 |

The aggregate direction improvement is not uniform. It is driven mainly by
OSD-871 targets, especially the large sporadic Parkinson disease without
microglia contrast. The microglia-matched adapter is therefore best treated as a
condition-sensitivity diagnostic rather than a universally stronger baseline.

## Decision

Keep both adapters:

1. `organoid_source_transfer_empirical_signature`
   - role: first conservative global response-signature diagnostic;
   - best current all-gene rank behavior;
   - easier to interpret as a coarse cross-source spaceflight response.
2. `organoid_microglia_matched_source_transfer_empirical_signature`
   - role: secondary condition-sensitivity diagnostic;
   - slightly better aggregate direction agreement;
   - weaker rank behavior and mixed per-contrast deltas.

Do not promote either adapter to a primary metric or leaderboard result.

## Next Adapter Choice

The next block should design a partial control-only disease+microglia matched
diagnostic before moving to classifier coefficients.

Reason:

- disease context is not fully crossed across OSD-863 and OSD-871;
- `no_known_diseases` is the only disease context shared by both sources;
- the shared healthy/control block can test whether matching both disease
  context and microglia condition helps when the stratum actually exists;
- classifier coefficients are still not log2FC values and need a separate
  score contract before they can be mixed with response-signature log2FC
  metrics.

The next diagnostic must be explicitly partial coverage:

- score only target contrasts whose opposite-source training stratum exists;
- write skipped-contrast metadata for PPMS and sporadic Parkinson disease
  target contrasts;
- do not silently fall back to global or microglia-only signatures;
- report aggregate metrics only over emitted contrasts and label them as
  partial coverage.

## Recommended Next Block

`V9-ORG-023: Human organoid shared-control disease+microglia adapter design`

Expected outputs:

- `v9/human_organoid/reports/ORGANOID_SHARED_CONTROL_SIGNATURE_ADAPTER_DESIGN.md`
- implementation plan for an optional partial report under
  `v9/human_organoid/reports/shared_control_source_transfer_signature/`

Done criteria:

- the design names eligible target contrasts;
- skipped disease-specific contrasts are explicit;
- leakage boundaries match the current source-transfer adapters;
- scorer interpretation is partial-coverage and non-leaderboard.

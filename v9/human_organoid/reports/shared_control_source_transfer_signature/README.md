# Human Organoid Shared-Control Source-Transfer Response-Signature Baseline

Status: partial-coverage diagnostic baseline; not a leaderboard result.

This report generates `response_signature.csv.gz` from source-transfer
training samples matched to shared `no_known_diseases` plus each target
contrast's microglia condition. Source-specific disease contexts are
skipped without fallback. The DE reference table is not used during
signature generation; it is used only afterward for diagnostic scoring
by the evaluator.

The accompanying `predictions.csv` is a two-row evaluator plumbing
fixture. Classification metrics in this report are not interpretable
as organoid model performance; only the partial response-signature
diagnostics are reviewed.

- signature model: `organoid_shared_control_disease_microglia_source_transfer_empirical_signature`
- response rows: 111944
- conditioning strategy: `shared_control_disease_microglia_source_transfer`
- reference usage policy: `reference_not_used_for_signature_generation`
- claim boundary: shared-control partial diagnostic only, not a
  biological generalization claim and not a leaderboard claim.

## Coverage

| Coverage | Count |
|---|---:|
| emitted target contrasts | 4 |
| skipped target contrasts | 4 |
| emitted response rows | 111944 |

## Partial-Coverage Metrics

| Metric | Value |
|---|---:|
| `de_direction_match` | 0.589942 |
| `signature_rank_correlation` | 0.035667 |

These aggregate values cover emitted shared-control contrasts only and
are not directly comparable to all-eight-contrast aggregate reports.

## Emitted Contrast Comparison

| Contrast | Global Direction | Microglia Direction | Shared Direction | Global Rank | Microglia Rank | Shared Rank |
|---|---:|---:|---:|---:|---:|---:|
| `osd_863_no_known_diseases_with_microglia_leo_or_iss_vs_ground` | 0.673469 | 0.663265 | 0.577320 | 0.114822 | 0.089727 | 0.045019 |
| `osd_863_no_known_diseases_without_microglia_leo_or_iss_vs_ground` | 0.803681 | 0.785276 | 0.595092 | 0.215367 | 0.199685 | 0.033068 |
| `osd_871_no_known_diseases_with_microglia_leo_or_iss_vs_ground` | 0.736842 | 0.710526 | 0.605263 | 0.134653 | 0.076384 | 0.042236 |
| `osd_871_no_known_diseases_without_microglia_leo_or_iss_vs_ground` | 0.651584 | 0.663636 | 0.589041 | 0.071176 | 0.036763 | 0.037082 |

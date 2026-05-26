# Human Organoid Microglia-Matched Source-Transfer Response-Signature Baseline

Status: diagnostic baseline; not a leaderboard result.

This report generates `response_signature.csv.gz` from source-transfer
training samples matched to each target contrast's microglia condition.
The DE reference table is not used during signature generation; it is
used only afterward for diagnostic scoring by the evaluator.

The accompanying `predictions.csv` is a two-row evaluator plumbing
fixture. Classification metrics in this report are not interpretable
as organoid model performance; only the response-signature diagnostics
are reviewed.

- signature model: `organoid_microglia_matched_source_transfer_empirical_signature`
- response rows: 223888
- conditioning strategy: `microglia_matched_source_transfer`
- reference usage policy: `reference_not_used_for_signature_generation`
- claim boundary: microglia-matched source-transfer diagnostic only,
  not a biological generalization claim and not a leaderboard claim.

## Comparison To Global Source-Transfer

| Metric | Global Source-Transfer | Microglia-Matched | Delta |
|---|---:|---:|---:|
| `de_direction_match` | 0.770673 | 0.788913 | 0.018239 |
| `signature_rank_correlation` | 0.176008 | 0.150072 | -0.025936 |

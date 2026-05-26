# Human Organoid Source-Transfer Response-Signature Baseline

Status: diagnostic baseline; not a leaderboard result.

This report generates `response_signature.csv.gz` from source-transfer
training samples only. The DE reference table is not used during
signature generation; it is used only afterward for diagnostic
scoring by the evaluator.

The accompanying `predictions.csv` is a two-row evaluator plumbing
fixture. Classification metrics in this report are not interpretable
as organoid model performance; only the response-signature diagnostics
are reviewed.

- signature model: `organoid_source_transfer_empirical_signature`
- response rows: 223888
- reference usage policy: `reference_not_used_for_signature_generation`
- claim boundary: source-transfer diagnostic only, not a biological
  generalization claim and not a leaderboard claim.

# Human Organoid Response-Signature Smoke Test

Status: scorer smoke test; not a biological model result.

This directory exercises the response-signature evaluator against the
real derived human organoid DE reference. The response signature is
constructed directly from selected reference rows, so perfect metric
values only mean that the evaluator join/scoring path is functioning.

- reference: `v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz`
- response rows: 40
- claim boundary: diagnostic scorer plumbing only, not leaderboard
  evidence and not a model performance claim.

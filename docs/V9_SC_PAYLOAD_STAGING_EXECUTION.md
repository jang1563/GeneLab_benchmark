# V9-SC-006 Payload Staging Execution

Status: `payload_staging_execution_blocked_no_candidate_payload`

Task id: `draft_rrrm1_blood_single_cell_spaceflight`

Claim boundary: `payload_staging_execution_no_payload_or_score_claim`

## Decision

The canonical RRRM-1 blood AnnData payload is still not staged at
`v9/sc_spaceflight/payloads/rrrm1_blood/OSD-918_blood_rrrm1_bench.h5ad`. The preferred annotated legacy h5ad is the correct
first route if it can be found, but this execution package does not claim that
any h5ad was copied, downloaded, regenerated, or audited.

Selected route: `prepare_regeneration_from_starsolo_per_srx`

## Current Evidence

- Canonical payload status: `planned_repo_path_absent`
- Existing payload manifest status: `missing`
- Existing obs/var audit status: `obs_var_audit_skipped_no_local_payload`
- Candidate rows checked: 4
- Regeneration steps fixed: 5

## Regeneration Gate

If the annotated h5ad cannot be restaged, regeneration must start by confirming
per-SRX STARsolo matrices for `OSD-918` and then rerun the legacy RRRM-1 merge
and broad annotation scripts before applying the v9 obs/var/uns contract.

## Not Claimed

- No local h5ad payload is claimed.
- No payload checksum or shape is claimed beyond the existing missing-payload
  manifest.
- No obs/var pass is claimed while the canonical payload remains absent.
- No evaluator, leaderboard result, or legacy RRRM score promotion is claimed.

## Next Block

Run `V9-SC-006b: locate STARsolo per-SRX matrices or restage annotated h5ad`.

# V9-SC-004 AnnData Payload Staging And Obs/Var Audit Plan

Status: `payload_staging_plan_ready_no_local_payload`

Task id: `draft_rrrm1_blood_single_cell_spaceflight`

Claim boundary: `payload_staging_plan_only_no_local_payload_or_score_claim`

## Decision

The first `sc_spaceflight` payload target is `v9/sc_spaceflight/payloads/rrrm1_blood/OSD-918_blood_rrrm1_bench.h5ad` for
RRRM-1 blood (`OSD-918`). This block does not stage, copy,
download, or regenerate the payload. It fixes the gate that must be satisfied
before the manifest can become runnable.

## Required Gate

- Required obs fields: 9
- Required var fields: 2
- Required uns fields: 4
- Audit requirement rows: 27
- Staging actions: 7
- Local payload present now: `false`

## Source Preference

Use the annotated legacy RRRM-1 blood h5ad first if it can be staged locally,
because it is the candidate with broad cell-type labels. If it is unavailable,
regenerate from per-SRX STARsolo matrices and then reapply the v9 obs/var/uns
contract. A labeled but unannotated h5ad can help recover sample metadata but
is not enough for the `genelab_sc` diagnostic surface.

## Not Claimed

- No local h5ad payload is claimed.
- No payload checksum is claimed.
- No obs/var audit pass is claimed.
- No evaluator, leaderboard result, or legacy RRRM score promotion is claimed.

## Inputs Used

- `v9/sc_spaceflight/task_manifests/draft_rrrm1_blood_single_cell_spaceflight.json`
- `v9/sc_spaceflight/metric_spec_summary.csv`
- `v2/docs/RRRM1_BENCHMARK_READY_MANIFEST_2026-03-12.csv`
- `v2/docs/RRRM1_BROAD_ANNOTATION_SUMMARY_2026-03-12.md`

## Next Block

Run `V9-SC-005: AnnData obs/var audit implementation`. That block should add
the actual h5ad reader/auditor, but still skip cleanly until a local payload is
available.

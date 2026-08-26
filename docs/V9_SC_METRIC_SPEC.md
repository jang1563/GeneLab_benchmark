# V9-SC-003 `genelab_sc` Metric Specification

Status: `genelab_sc_metric_spec_ready_no_evaluator`

Profile: `genelab_sc`

Claim boundary: `genelab_sc_metric_spec_only_no_evaluator_or_score_claim`

## Scope

This is a metric specification only. It defines formulas, required inputs, and
skip policy for the first non-runnable RRRM-1 blood `sc_spaceflight` manifest:
`draft_rrrm1_blood_single_cell_spaceflight`.

No evaluator, local AnnData payload, benchmark score, or legacy RRRM score
promotion is claimed by this block.

## Metric Roles

- Primary after payload freeze: `state_label_auroc`, `state_label_auprc`
- Diagnostic representation: `mission_discrimination`
- Diagnostic DE recovery: `de_overlap_at_n`, `de_direction_match`
- Optional reconstruction: `expression_mae_when_applicable`

## Required Boundaries

- State-label metrics require `predictions.csv` with cell labels and
  `flight_probability`.
- Embedding diagnostics require `embedding_*` columns and are skippable.
- DE recovery metrics require both `response_signature.csv` and a frozen
  v9 single-cell DE reference table.
- Expression MAE is optional and must skip unless a reconstruction artifact and
  normalization layer are declared.
- Every skipped metric must report `NA` plus a machine-readable skip reason.

## Next Block

Run `V9-SC-004: AnnData payload staging and obs/var audit plan`. The metric
spec is now fixed enough to decide exactly which local AnnData object and
metadata fields must be staged or regenerated before evaluator implementation.

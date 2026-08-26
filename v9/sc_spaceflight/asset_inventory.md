# V9-SC-001 RRRM Asset Inventory

Status: `legacy_rrrm_assets_indexed_no_local_anndata_payload`

Claim boundary: `single_cell_asset_inventory_only_no_v9_sc_task_or_payload_claim`

## Inventory Summary

- Total legacy RRRM/single-cell asset paths: 54
- RRRM-1 asset paths: 41
- RRRM-2 asset paths: 13
- Documentation rows: 8
- Script rows: 31
- Evaluation/config rows: 5
- Figure rows: 2
- Processed-summary rows: 2
- Generated-cache rows excluded from promotion: 6
- Local `.h5ad` files found: 0
- Local `.loom` files found: 0
- Local `.mtx` files found: 0
- v9 `sc_spaceflight` task manifests found: 0
- Metric profile status: `genelab_sc_profile_present`

## Decision

The single-cell lane has enough legacy RRRM documentation, scripts, figures,
and evaluation surfaces to justify a v9 flagship manifest draft, but not enough
local matrix payload evidence to claim a runnable v9 single-cell task. This
inventory is therefore evidence only: it does not create a v9 `sc_spaceflight`
task manifest, does not claim local AnnData payload availability, and does not
promote legacy RRRM scores as v9 benchmark results.

## Strongest Candidate Surface

- RRRM-1 has compact sample/QC/benchmark-ready manifests in `v2/docs/`.
- RRRM-2 has legacy benchmark scripts and evaluation JSON surfaces under
  `v3/scripts/` and `v3/evaluation/`.
- The next manifest draft should choose one candidate task and explicitly list
  required AnnData `obs`, `var`, label, split, and provenance fields before any
  evaluator implementation.

## Next Block

Run `V9-SC-002: AnnData task manifest draft`. Start from this inventory and
draft one non-runnable `sc_spaceflight` manifest with explicit payload blockers
instead of trying to regenerate RRRM matrices in the same block.

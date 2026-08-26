# V9-SC-002 AnnData Task Manifest Draft

Status: `draft_anndata_manifest_contract_ready_payload_blocked`

Task id: `draft_rrrm1_blood_single_cell_spaceflight`

Claim boundary: `anndata_manifest_contract_only_no_local_payload_or_score_claim`

## Draft Decision

Draft one non-runnable `sc_spaceflight` manifest for RRRM-1 blood
(`OSD-918`). This is the most compact first single-cell flagship
contract because the repo already has RRRM-1 sample, QC, annotation, and
benchmark-ready evidence for blood, while local AnnData payload files are still
absent.

## Evidence

- Samples in condition map: 8
- Flight samples: 4
- Ground samples: 4
- QC cell count: 4395
- QC gene count: 19064
- Local `.h5ad` payload count: 0
- Minimal manifest validator: pass

## Not Claimed

- This is not a runnable v9 single-cell benchmark task.
- It does not claim local AnnData payload availability.
- It does not promote legacy RRRM scripts, figures, or scores as v9 results.
- It does not define final `genelab_sc` formulas.

## Next Block

Run `V9-SC-003: genelab_sc metric specification`. Define metric formulas,
required inputs, skip policy, and which metrics remain diagnostic before any
single-cell evaluator implementation.

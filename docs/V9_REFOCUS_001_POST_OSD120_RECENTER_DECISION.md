# V9-REFOCUS-001 Post-OSD-120 Recenter Decision

Status: `selected_public_bulk_alpha_recenter`

Selected next lane: `public_bulk_alpha`

Selected next block: `V9-BULK-ALPHA-001: public bulk alpha freeze-gap matrix`

## Decision

Choose public bulk alpha readiness before the first single-cell flagship
implementation.

This is not a retreat from single-cell work. It is a sequencing correction:
the public bulk scaffold already has task manifests, source inventory,
checksum-manifest evidence, baseline outputs, a draft Data Package, and a
dataset-card draft. It still lacks a clean alpha gap matrix and payload hash
boundary. The single-cell lane has valuable RRRM legacy assets and a
`genelab_sc` metric profile, but no v9 `sc_spaceflight` task manifests and no
local h5ad/loom/mtx files found by the current repo scan.

## Readiness Snapshot

| Lane | Score | Decision | Main gap |
|---|---:|---|---|
| Public bulk alpha | 90 | selected next | payload hash/release claim boundary |
| Single-cell flagship | 35 | deferred | RRRM asset inventory before AnnData task cards |

## Selected Next Block

`V9-BULK-ALPHA-001: public bulk alpha freeze-gap matrix` should produce a machine-readable gap matrix for the
public bulk alpha scaffold, with explicit pass/blocker rows for source
inventory, checksum evidence, payload hash verification, package metadata,
dataset-card language, and baseline output evidence.

## External Method Anchors

- AnnData stores expression matrix `X` with observation metadata in `obs`,
  variable metadata in `var`, and h5ad-backed on-disk storage. A single-cell
  flagship should therefore start from explicit AnnData/obs/var inventory.
- cell-eval expects predicted and real AnnData inputs for perturbation-response
  evaluation, so v9 should not add model adapters before the real/predicted
  AnnData contract is defined.
- OpenProblems emphasizes task APIs and benchmark contracts; v9 should follow
  that discipline for `sc_spaceflight` after the asset inventory.

## Claim Boundary

This is a recenter decision only. It does not create a new benchmark release,
does not promote OSD-120 to the v9-alpha core, and does not claim a
single-cell flagship result.

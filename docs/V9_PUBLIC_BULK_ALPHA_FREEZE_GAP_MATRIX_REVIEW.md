# V9-BULK-ALPHA-001 Public Bulk Alpha Freeze-Gap Matrix Review

Status: `metadata_alpha_gap_matrix_ready_payload_hash_blocked`

Claim boundary: `public_bulk_alpha_gap_matrix_no_release_approval`

## Decision

The public bulk alpha scaffold is close enough to stay as the active lane, but
not ready for frozen release language. The strongest blocker is payload hash
verification: `0/22` public bulk sources are
`freeze_ready=true`, even though `22/22` have
parsed OSDR checksum-manifest evidence.

## Matrix Summary

- Pass rows: 6
- Blocker rows: 2
- Needs-update rows: 2
- Allowed current claims: 3
- Prohibited release claims: 3

## Current Allowed Language

It is safe to say that the public bulk scaffold has indexed task manifests,
source inventory, OSDR API/checksum-manifest evidence, draft Data Package
metadata, dataset-card draft language, and reproduced simple baselines.

## Current Prohibited Language

Do not claim a frozen public benchmark release, locally verified payload mirror,
complete release metadata descriptor, extension-track leaderboard, or approved
alpha release.

## Next Block

Run `V9-BULK-ALPHA-002: metadata-only alpha snapshot decision`. That block
should decide whether the project can publish a metadata-only alpha snapshot
with explicit payload-hash blockers, or whether payload mirroring and local hash
verification must precede any alpha wording.

## External Guidance Anchors

- OSDR API evidence can establish source/file-list traceability, but local
  payload verification is a separate benchmark-release claim.
- Frictionless Data Package resources can describe current tables and reports,
  but package metadata must not imply release completeness while payload and
  license/citation boundaries remain unresolved.

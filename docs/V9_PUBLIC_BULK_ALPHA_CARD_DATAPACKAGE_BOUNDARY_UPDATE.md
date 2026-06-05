# V9-BULK-ALPHA-003 Dataset Card And Data Package Boundary Update

Status: `metadata_alpha_boundary_applied`

Claim boundary: `metadata_only_public_bulk_alpha_no_payload_release`

## Decision Applied

The V9-BULK-ALPHA-002 metadata-only alpha decision is now reflected in both
public-facing metadata surfaces:

- `docs/v9_hf_dataset_card.md`
- `v9/datapackage.draft.json`

Both artifacts say the public bulk lane is a metadata-only alpha snapshot, not
a frozen payload release or locally hash-verified payload bundle.

## Data Package Update

The draft Data Package now records:

- `spacebio_bench:release_status = metadata_alpha_not_frozen`
- `spacebio_bench:alpha_snapshot_status = metadata_only_alpha_snapshot`
- `spacebio_bench:claim_boundary = metadata_only_public_bulk_alpha_no_payload_release`
- `spacebio_bench:payload_release_allowed = false`
- `spacebio_bench:payload_verification_status = checksum_manifests_parsed_payloads_not_hashed`

It contains 21 resources. Ten resources are explicit
`alpha_boundary_metadata` tables from the public bulk gap matrix and snapshot
decision packages.

## Dataset Card Update

The dataset card now uses `SpaceBio-Bench v9 Public Bulk Metadata Alpha` as the
public label and states:

- Release status: metadata-only alpha snapshot, not frozen.
- The card is not frozen release language.
- The package does not include a local payload mirror.
- Payload-level hash verification remains pending.
- Organoid and multispecies draft tracks are not public bulk alpha core tasks.

## Next

The public bulk alpha claim/payload boundary is now explicit enough to return
to `V9-SC-001: RRRM asset inventory`, unless the user explicitly chooses to
start the deferred payload-mirror freeze lane first.

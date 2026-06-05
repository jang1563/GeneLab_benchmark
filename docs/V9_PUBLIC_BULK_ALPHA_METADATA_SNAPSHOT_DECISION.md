# V9-BULK-ALPHA-002 Metadata-Only Alpha Snapshot Decision

Status: `metadata_only_alpha_snapshot_allowed_with_payload_blockers`

Selected path: `metadata_only_alpha_snapshot`

Deferred path: `payload_mirror_first`

Claim boundary: `metadata_only_public_bulk_alpha_no_payload_release`

## Decision

Proceed with a metadata-only alpha snapshot for the public bulk lane, with
explicit payload-hash blockers. This is not a frozen payload release. It is a
bounded review surface for task manifests, source inventory, OSDR API and
checksum-manifest evidence, baseline summaries, and provenance reports.

The payload-release path remains blocked because `0/22`
public bulk sources are locally payload-hash verified, while
`22/22` sources have parsed checksum-manifest
evidence.

## Option Comparison

| Path | Decision | Status |
| --- | --- | --- |
| metadata-only alpha snapshot | selected | allowed with explicit blockers |
| payload mirror first | deferred | valid for future payload release, not required before metadata alpha |
| no alpha until payload frozen | rejected | too conservative for metadata scaffold |

## Allowed Language

- `SpaceBio-Bench v9 public bulk metadata alpha`
- The snapshot documents public mouse bulk LOMO task/source/provenance metadata.
- OSDR file-list and checksum-manifest evidence has been parsed for all 22
  public bulk source rows.
- Payload mirroring and local payload-hash verification remain pending.

## Prohibited Language

- Frozen public benchmark release.
- Frozen payload mirror.
- Locally hash-verified data bundle.
- DOI/archive release, complete release Data Package, or leaderboard claim.
- Organoid or multispecies draft tracks as public bulk alpha core tasks.

## External Guidance Anchors

- Hugging Face dataset cards are README/metadata surfaces meant to help users
  understand dataset contents, context, and responsible use:
  https://huggingface.co/docs/hub/datasets-cards
- Frictionless Data Package descriptors separate package metadata from resource
  entries and can describe metadata resources without implying a local payload
  mirror:
  https://specs.frictionlessdata.io/data-package/
- OSDR API file-list and metadata endpoints support source/file traceability,
  while local benchmark payload hashing remains a separate project claim:
  https://visualization.osdr.nasa.gov/biodata/api/
- NASA OSDR should remain the credited upstream source for space biology data:
  https://science.nasa.gov/reference/osdr-faq/

## Next Block

Run `V9-BULK-ALPHA-003: dataset card and Data Package alpha boundary update`. That block should update
`docs/v9_hf_dataset_card.md` and `v9/datapackage.draft.json` using the claim
boundary in this decision package.

# V9-SC-006c OSDR Processed Payload Discovery

Status: `osdr_processed_payload_unavailable_owner_scratch_request_required`

Task id: `draft_rrrm1_blood_single_cell_spaceflight`

Claim boundary: `osdr_processed_payload_discovery_only_no_payload_copy_or_score_claim`

## Decision

The OSDR Biodata file-list endpoint for `OSD-918` was queried at
`https://visualization.osdr.nasa.gov/biodata/api/v2/dataset/OSD-918/files/`. The current file list records 19 files: 1
metadata file, 1 raw checksum manifest,
1 raw MultiQC report, and
16 raw FASTQ files. It does not list a processed h5ad,
processed STARsolo matrix bundle, or processed checksum manifest for the
canonical RRRM-1 blood AnnData target.

## Coverage

- Expected OSD-918 blood SRX rows: 8
- Complete expected raw FASTQ pairs listed by OSDR: 8/8
- Processed h5ad files listed: 0
- Processed STARsolo matrix rows listed: 0
- Canonical copy allowed: `false`
- Regeneration allowed: `false`
- Owner scratch request required: `true`

## Not Claimed

- No OSDR payload was downloaded.
- No source payload hash, AnnData shape, or obs/var pass is claimed.
- No canonical v9 h5ad copy was made.
- No raw FASTQ-to-STARsolo regeneration was started.
- No evaluator, leaderboard result, or legacy RRRM score promotion is claimed.

## Next Block

Run `V9-SC-006d: owner scratch path intake or raw FASTQ regeneration feasibility decision`.

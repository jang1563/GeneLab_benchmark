# OSD-120 Payload Freeze Manifest Review

Status: V9-MULTI-027 complete as of 2026-05-26.

This block addresses the first public-alpha blocker from the
release-readiness gap audit: OSD-120 payload freeze. It does not download the
full OSDR processed payload set. Instead, it makes the diagnostic package scope
explicit by comparing the live OSDR processed MD5 manifest against the two local
OSD-120 payload files required by the current diagnostic task.

Primary outputs:

- `interaction_payload_freeze_manifest/payload_freeze_summary.csv`
- `interaction_payload_freeze_manifest/payload_freeze_summary.json`
- `interaction_payload_freeze_manifest/payload_freeze_manifest.csv`
- `interaction_payload_freeze_manifest/payload_freeze_manifest.json`

## Source Manifest

The source checksum manifest is the OSD-120 processed GLbulkRNAseq MD5 manifest:

`https://osdr.nasa.gov/geode-py/ws/studies/OSD-120/download?source=datamanager&file=GLDS-120_rna_seq_processed_md5sum_GLbulkRNAseq.tsv`

The fetched manifest SHA-256 is
`2035d1a07e862f16044678ace3c0e4747c4d3b0734872cacac29e3dfd86b9523`, matching
the earlier multispecies source checksum audit.

## Result

| field | value |
|---|---:|
| parsed OSDR checksum entries | 533 |
| diagnostic-required payloads | 2 |
| diagnostic-required payloads matched | 2 |
| diagnostic-required payloads missing | 0 |
| diagnostic-required checksum mismatches | 0 |
| OSDR processed payloads outside current diagnostic scope | 531 |
| diagnostic required payload freeze ready | True |
| full OSDR payload freeze ready | False |

The two diagnostic-required OSDR payloads both match local MD5 hashes:

| OSDR manifest payload | local file | status |
|---|---|---|
| `SampleTable_GLbulkRNAseq.csv` | `data/multispecies/arabidopsis/GLDS-120_rna_seq_SampleTable_GLbulkRNAseq.csv` | `required_payload_md5_matched` |
| `Normalized_Counts_GLbulkRNAseq.csv` | `data/multispecies/arabidopsis/GLDS-120_rna_seq_Normalized_Counts_GLbulkRNAseq.csv` | `required_payload_md5_matched` |

The rRNA-removed normalized matrix, DE tables, trimming reports, and other OSDR
processed outputs remain listed as
`osdr_processed_payload_not_required_for_diagnostic`. They are not missing
diagnostic inputs; they are simply outside the current OSD-120 diagnostic
freeze scope.

## Decision

The OSD-120 diagnostic-required payload scope is frozen for the current draft
diagnostic package. This is narrower than a full OSDR payload freeze.

Safe wording:

- The current OSD-120 diagnostic package uses two OSDR processed payloads.
- Both required local payloads match the OSDR processed MD5 manifest.
- The full OSDR processed payload set is not downloaded or frozen.

Not safe:

- Do not claim a complete OSD-120 OSDR mirror.
- Do not claim all processed OSDR payloads are locally hash-verified.
- Do not promote the diagnostic package to a full benchmark release solely from
  this two-file freeze.

## Next Block

Proceed to `V9-MULTI-028: OSD-120 diagnostic public-alpha card draft`. The card
should state this narrower freeze boundary clearly: diagnostic-required payloads
are MD5 matched; the broader OSDR processed payload set remains out of scope.


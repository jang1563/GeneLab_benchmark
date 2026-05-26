# Multispecies Checksum And Local Payload Audit

Status: V9-MULTI-005 review  
Date: 2026-05-23

## Scope

This review audits OSDR file-list/checksum evidence for the three draft
multispecies candidates:

- OSD-207 / GLDS-207: Drosophila whole-body bulk RNA-seq
- OSD-37 / GLDS-37: Arabidopsis seedling-pool bulk RNA-seq
- OSD-120 / GLDS-120: Arabidopsis root bulk RNA-seq with light/genotype strata

The local sample/matrix audit from V9-MULTI-004 already confirmed that the
processed count matrices align with parsed sample factors. This block checks
whether those local files are backed by OSDR checksum manifest evidence.

## External Source Check

Primary-source anchors checked during this block:

- NASA OSDR Biological Data API documentation:
  https://visualization.osdr.nasa.gov/biodata/api/
- NASA OSDR API announcement and access page:
  https://www.nasa.gov/osdr-api/

Operating implication: source evidence should be derived from OSDR API file
listings and checksum manifests, not only from hand-maintained local file names.

## Generated Artifacts

- `v9/multispecies/source_checksum_audit.draft.csv`
- `v9/multispecies/source_checksum_audit.draft.json`

Command:

```bash
/usr/local/bin/python3 scripts/audit_v9_source_checksums.py \
  --source-inventory v9/multispecies/source_inventory.draft.csv \
  --csv v9/multispecies/source_checksum_audit.draft.csv \
  --json v9/multispecies/source_checksum_audit.draft.json
```

## Source-Level Result

| source_id | API status | OSDR files | checksum manifests | parsed checksum entries | payload-name matches | status |
|---|---:|---:|---:|---:|---:|---|
| OSD-207 | ok | 534 | 1 | 370 | 366 | checksum_manifest_parsed |
| OSD-37 | ok | 986 | 2 | 927 | 926 | checksum_manifest_parsed |
| OSD-120 | ok | 851 | 1 | 533 | 532 | checksum_manifest_parsed |

All three draft sources have live OSDR API file-list evidence and parseable MD5
checksum manifests. The generic source checksum audit keeps
`freeze_ready=false` because that script records manifest evidence but does not
download and hash every payload listed by each source.

## Local Payload Spot-Check

The six local files used by the V9-MULTI-004 input scaffold were hashed locally
and compared against the processed OSDR MD5 manifests. OSDR manifests list the
payload basenames without the local `GLDS-*` prefix, but the hashes match.

| source_id | local file | manifest payload | expected md5 | observed md5 | result |
|---|---|---|---|---|---|
| OSD-207 | `data/multispecies/drosophila/GLDS-207_rna_seq_Normalized_Counts_GLbulkRNAseq.csv` | `Normalized_Counts_GLbulkRNAseq.csv` | `a556e54773a6b3a6bd3297aac556bd6e` | `a556e54773a6b3a6bd3297aac556bd6e` | match |
| OSD-207 | `data/multispecies/drosophila/GLDS-207_rna_seq_SampleTable_GLbulkRNAseq.csv` | `SampleTable_GLbulkRNAseq.csv` | `624e6c55be02361b7311e665fa27a46f` | `624e6c55be02361b7311e665fa27a46f` | match |
| OSD-37 | `data/multispecies/arabidopsis/GLDS-37_rna_seq_Normalized_Counts_GLbulkRNAseq.csv` | `Normalized_Counts_GLbulkRNAseq.csv` | `538c2d107eb466d766a6b6bd9e0361ae` | `538c2d107eb466d766a6b6bd9e0361ae` | match |
| OSD-37 | `data/multispecies/arabidopsis/GLDS-37_rna_seq_SampleTable_GLbulkRNAseq.csv` | `SampleTable_GLbulkRNAseq.csv` | `2bb3c6dfe39c61933eec5daf95159947` | `2bb3c6dfe39c61933eec5daf95159947` | match |
| OSD-120 | `data/multispecies/arabidopsis/GLDS-120_rna_seq_Normalized_Counts_GLbulkRNAseq.csv` | `Normalized_Counts_GLbulkRNAseq.csv` | `d2e9e835f60ef0752a47ecdcbca1b9af` | `d2e9e835f60ef0752a47ecdcbca1b9af` | match |
| OSD-120 | `data/multispecies/arabidopsis/GLDS-120_rna_seq_SampleTable_GLbulkRNAseq.csv` | `SampleTable_GLbulkRNAseq.csv` | `7e4a87241e4272909c8d9a47b3d45e1f` | `7e4a87241e4272909c8d9a47b3d45e1f` | match |

## Readiness Decision

OSD-37 is ready to proceed to a species-native Arabidopsis task-manifest design
block, subject to the existing draft/not-frozen release boundary.

OSD-207 is ready to proceed to a species-native Drosophila task-manifest design
block, subject to the same draft/not-frozen boundary.

OSD-120 has adequate source, sample, matrix, and local payload checksum
evidence, but remains deferred from the first species-native task because its
scientific design requires explicit genotype/ecotype by light-treatment
stratification. It should be handled as an interaction-design task, not merged
into the first OSD-37 plant task.

## Remaining Release Boundary

This audit is enough for draft task-manifest design. It is not yet a complete
public release freeze because:

- the generic source audit does not hash every OSDR-listed payload;
- only the six local files used by the current multispecies scaffold were
  locally hash-checked;
- the task manifests have not yet declared split units, baseline outputs, or
  release-card wording for these species-native tasks.

Next action: draft a multispecies task-manifest design for OSD-37 and OSD-207,
keeping OSD-120 in an interaction-task queue.

# OSD-120 Release-Readiness Gap Audit

Status: V9-MULTI-026 complete as of 2026-05-26.

This block reviews whether the OSD-120 diagnostic package can move from local
draft evidence to a cleaner public alpha artifact. It does not rerun models and
does not expand the scientific claim boundary. The audit ties the existing
artifact manifest, claim map, task manifest, checksum audit, sample/matrix
audits, and external release conventions into one gap table.

Primary outputs:

- `interaction_release_readiness_gap_audit/release_readiness_summary.csv`
- `interaction_release_readiness_gap_audit/release_readiness_summary.json`
- `interaction_release_readiness_gap_audit/release_readiness_gap_table.csv`
- `interaction_release_readiness_gap_audit/release_readiness_gap_table.json`
- `interaction_release_readiness_gap_audit/release_readiness_external_references.csv`
- `interaction_release_readiness_gap_audit/release_readiness_external_references.json`

## Decision

The current OSD-120 diagnostic package is internally traceable, but it is not
yet public-alpha ready.

Summary:

| field | value |
|---|---:|
| public alpha ready | False |
| blockers | 3 |
| needs-work items | 3 |
| pass items | 5 |
| acceptable draft limitations | 1 |
| indexed diagnostic artifacts | 26 |
| missing artifacts | 0 |
| unhashed artifacts | 0 |
| mapped claims | 7 |

The package can be used as a draft diagnostic evidence bundle. It should not be
called a frozen v9 benchmark, release snapshot, leaderboard baseline, or
leave-one-mission-out result.

## Blockers Before Public Alpha

| blocker | evidence | required fix |
|---|---|---|
| Full OSDR payload freeze missing | `source_checksum_audit.draft.csv` has `freeze_ready=false` for OSD-120 | Build an OSD-120 payload freeze manifest with expected MD5 values, observed local hashes, and missing payload rows. |
| Source remains draft/not frozen | task manifest and source inventory still use draft/not-frozen status | Define a public-alpha release target only after payload freeze and release metadata are complete. |
| OSD-120-specific public card missing | no source-specific card exists for this diagnostic package | Draft a card with OSDR credit, source scope, intended use, artifact structure, limitations, citation/version language, and claim boundaries. |

## Passed Draft-Readiness Checks

- The diagnostic artifact manifest indexes 26 files and all exist with SHA-256
  hashes.
- The claim map links seven claims to artifacts, tests, limitations, and
  external context URLs.
- OSD-120 expression-matrix columns align to 36 parsed sample-factor rows.
- The local OSD-120 SampleTable and normalized-count matrix have MD5 spot-check
  evidence against the processed OSDR checksum manifest.
- External OSD-120 context is limited to source/task framing, not performance
  claims.

## Needs Work Before Broader Release

- Add machine-readable release metadata for the OSD-120 diagnostic package,
  ideally as a Data Package or RO-Crate-linked bundle.
- Add a release-specific rebuild gate and environment/runtime lock for the
  diagnostic outputs.
- Define a versioned archive/DOI plan only after the package is promoted beyond
  internal draft.

## External Standards Used

The audit records seven external reference anchors:

- NASA OSDR repository overview:
  https://science.nasa.gov/biological-physical/data/osdr/
- NASA OSDR FAQ and credit/citation guidance:
  https://science.nasa.gov/reference/osdr-faq/
- FAIR principles:
  https://www.go-fair.org/fair-principles
- RO-Crate 1.2 specification:
  https://www.researchobject.org/ro-crate/specification.html
- DataCite Metadata Schema:
  https://schema.datacite.org/
- Hugging Face dataset-card documentation:
  https://huggingface.co/docs/hub/datasets-cards
- GitHub/Zenodo citation guidance:
  https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content

Operating implication: public alpha needs source identifiers, access and reuse
metadata, artifact and workflow provenance, clear dataset-card language, and a
version/archive plan. The current local artifact manifest is a strong start but
does not replace payload freeze or release metadata.

## Next Blocks

Recommended sequence:

1. `V9-MULTI-027`: OSD-120 payload freeze manifest.
2. `V9-MULTI-028`: OSD-120 diagnostic public-alpha card draft.
3. `V9-MULTI-029`: OSD-120 diagnostic rebuild gate and environment lock.
4. `V9-THEN-005`: RO-Crate/Data Package export design integration.


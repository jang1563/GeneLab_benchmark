# OSD-120 Public Metadata Package Review

Block: V9-MULTI-030

Package id: `osd120_diagnostic_public_metadata_package`

Status: `draft_metadata_skeleton_ready_not_release_frozen`

This block separates the current public diagnostic metadata surface from future
source-package or frozen-benchmark release claims. It creates a machine-readable
metadata skeleton aligned with DataCite 4.7, RO-Crate 1.2, Hugging Face dataset
card metadata, and OSDR citation/credit expectations.

## Release Target Decision

- `diagnostic_alpha_metadata_draft`: ready_now_with_draft_limitations; public_now=True.
- `source_specific_public_alpha_package`: partially_ready_source_release_target_pending; public_now=False.
- `full_osdr_payload_mirror_release`: blocked_full_payload_freeze_missing; public_now=False.
- `frozen_v9_benchmark_release`: blocked_version_doi_and_broader_validation_missing; public_now=False.

Only `diagnostic_alpha_metadata_draft` is public-now, and only with draft
limitations. Source-specific alpha packaging, full OSDR payload mirroring, and
frozen benchmark release claims remain blocked.

## Metadata Field Status

- present: 12
- partial: 3
- placeholder: 5

The placeholders are intentional: they prevent false DOI, creator, publisher,
rights, and license claims before a formal release target exists.

## Standards Anchors

- DataCite Metadata Schema 4.7: `https://schema.datacite.org/`
- RO-Crate 1.2 metadata: `https://www.researchobject.org/ro-crate/specification/1.2/metadata.html`
- Hugging Face dataset cards: `https://huggingface.co/docs/hub/datasets-cards`
- NASA OSDR FAQ/citation guidance: `https://science.nasa.gov/reference/osdr-faq/`

## Claim Boundary

This package is a diagnostic public metadata skeleton. It does not claim a
complete local OSDR payload mirror, leave-one-mission-out generalization,
validated biomarkers, operational plant-growth recommendations, or a frozen v9
benchmark release.

## Next Block

V9-MULTI-032 should make explicit archive identifier, release version,
creator/contributor, and license/rights decisions before any citable archive or
DOI path is attempted.

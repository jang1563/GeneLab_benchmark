# OSD-120 Citation Metadata Fill Review

Block: V9-MULTI-033

Fill id: `osd120_release_owner_citation_metadata_fill`

Status: `no_owner_metadata_supplied_no_descriptor_changes`

This block creates a release-owner metadata intake/fill scaffold. It does not
mutate the draft RO-Crate or Data Package descriptors, does not mint a DOI, does
not infer creator order, and does not choose a license.

## Owner-Supplied Fields

- No owner-supplied citation metadata fields were provided.

## Remaining Release Blockers

- `future_archive_identifier`: versioned_archive_plan.
- `release_version_string`: versioned_archive_plan.
- `release_date`: versioned_archive_plan.
- `publication_year`: versioned_archive_plan.
- `local_package_creators`: formal_license_and_creator_metadata.
- `local_package_contributors`: formal_license_and_creator_metadata.
- `publisher_maintainer_identity`: formal_license_and_creator_metadata.
- `local_code_license`: formal_license_and_creator_metadata.
- `local_metadata_tables_license`: formal_license_and_creator_metadata.
- `upstream_osdr_dataset_citation`: source_inventory_release_target.
- `repository_url`: repository_public_and_owner_access_unconfirmed.
- `citation_cff_type`: formal_license_and_creator_metadata.

## Descriptor Patch Policy

Patch previews are emitted only for fields explicitly supplied by a release
owner. The current generated descriptors remain diagnostic drafts until owner
metadata, archive route, version, creator/contributor, publisher, and license
decisions are complete.

## External Guidance Anchors

- GitHub/Zenodo DOI release paths require a public repository, owner access,
  GitHub release, and license.
- GitHub CITATION.cff can describe software or dataset citation metadata, but
  author, version, DOI or URL, and release-date fields must be supplied.
- OSDR source citation should come from the OSDR study page citation button in
  BibTeX or RIS form.
- DataCite 4.7 citable metadata requires identifier, creator, title,
  publisher, publication year, and resource type.
- Local license metadata should use an SPDX identifier or explicit custom terms
  after review.

## Claim Boundary

The current output is an intake and preview scaffold only. No archive
identifier, no DOI, no frozen release version, no creator list, no publisher,
no local package license, and no descriptor mutation are claimed.

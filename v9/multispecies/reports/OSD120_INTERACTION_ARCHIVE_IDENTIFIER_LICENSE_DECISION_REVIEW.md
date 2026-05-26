# OSD-120 Archive Identifier And License Decision Gate Review

Block: V9-MULTI-032

Decision id: `osd120_archive_identifier_license_decision_gate`

Status: `current_draft_no_archive_selected_release_decisions_blocked`

This block does not mint a DOI, choose a license, or invent a creator list. It
selects a safe current-draft path with no local archive identifier, records OSDR
source citation as upstream credit, and keeps release-owner decisions blocked
until supplied explicitly.

## Archive Identifier Options

- `current_no_archive_diagnostic_draft`: selected_for_current_draft_only; current_draft_selected=True.
- `zenodo_github_release_doi`: deferred_release_owner_required; current_draft_selected=False.
- `citation_cff_without_archive`: deferred_after_creator_version_decision; current_draft_selected=False.
- `osdr_source_citation_only`: selected_for_upstream_credit_not_local_archive; current_draft_selected=True.
- `software_heritage_swhid_related_identifier`: deferred_code_archive_related_identifier_only; current_draft_selected=False.
- `full_osdr_payload_mirror_archive`: blocked_full_payload_freeze_missing; current_draft_selected=False.

## License And Rights Decisions

- `upstream_osdr_data_credit`: needs_review.
- `local_code_license`: blocked.
- `local_metadata_tables_license`: blocked.
- `diagnostic_payload_mirror_rights`: blocked.
- `third_party_literature_links`: selected_citation_only.

## Creator And Contributor Decisions

- `upstream_osdr_dataset_citation`: needs_review.
- `local_package_creators`: blocked.
- `local_package_contributors`: blocked.
- `publisher_maintainer_identity`: blocked.
- `package_title`: selected_for_current_draft.
- `release_version_string`: blocked.

## External Guidance Anchors

- GitHub/Zenodo DOI requires a public repository, repository-owner access,
  GitHub release, and license.
- GitHub CITATION.cff can describe dataset/software citation metadata, but
  needs author, version, and DOI or URL fields.
- OSDR dataset citation should be taken from the OSDR study page citation
  button in BibTeX or RIS form.
- DataCite 4.7 citable metadata requires identifier, creators, title,
  publisher, publication year, and resource type.
- License decisions should use a clear SPDX identifier or explicit custom
  terms after review.

## Claim Boundary

The current diagnostic package remains draft metadata only. No DOI, no archive
identifier, no frozen release version, no local package license, no full OSDR
payload mirror, and no frozen benchmark release are claimed.

## Next Block

V9-MULTI-033 should fill citation metadata only after the release owner supplies
creator/contributor, license/rights, archive route, and version decisions.

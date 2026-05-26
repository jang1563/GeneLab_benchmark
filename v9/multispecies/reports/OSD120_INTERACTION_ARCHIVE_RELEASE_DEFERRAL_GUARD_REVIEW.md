# OSD-120 Archive Release Deferral And Application Guard Review

Block: V9-MULTI-034

Guard id: `osd120_archive_release_deferral_application_guard`

Status: `archive_release_deferred_no_owner_metadata`

The V9-MULTI-033 owner metadata fill has no owner-supplied citation metadata.
This guard therefore defers archive release, blocks descriptor mutation, and
keeps the current OSD-120 public surface as diagnostic metadata only.

## Guard Blockers

- `owner_metadata_file_present`: defer_archive_release_and_keep_diagnostic_metadata_only.
- `owner_supplied_patch_preview_present`: do_not_apply_empty_patch_preview.
- `all_release_blockers_resolved`: defer_archive_release_until_blockers_are_zero.
- `archive_identifier_supplied`: keep_no_archive_identifier_for_current_draft.
- `release_version_date_supplied`: do_not_use_build_date_or_draft_string_as_release_version.
- `creator_publisher_metadata_supplied`: do_not_infer_authorship_or_publisher_from_repository_state.
- `license_rights_metadata_supplied`: do_not_choose_spdx_or_custom_license_without_owner_review.
- `osdr_source_citation_supplied`: retain_generic_osdr_credit_but_defer_source_alpha_citation.
- `citation_release_surface_supplied`: do_not_emit_citation_cff_or_github_release_metadata.

## Deferral Actions

- `supply_release_owner_metadata_file`: required_owner_action.
- `decide_archive_route_identifier`: required_owner_action.
- `freeze_version_date_year`: required_owner_action.
- `freeze_creator_contributor_publisher`: required_owner_action.
- `choose_license_scope`: required_owner_action.
- `capture_exact_osdr_study_citation`: required_owner_action.
- `retain_diagnostic_metadata_surface`: selected_deferral_action.
- `prevent_descriptor_mutation`: selected_guard_action.
- `defer_archive_release`: selected_deferral_action.

## Mutation Policy

No RO-Crate, Data Package, CITATION.cff, DOI, release tag, license, creator,
publisher, or source-citation freeze is applied by this block. A future
application block may mutate descriptors only if owner-supplied metadata exists,
all release blockers validate, and the patch preview is non-empty and reviewed.

## External Guidance Anchors

- GitHub/Zenodo DOI release paths depend on public repository access, owner
  approval where needed, a GitHub release, and explicit license/reuse terms.
- GitHub releases require a tag and release metadata; draft/pre-release states
  are separate from a citable archive.
- Zenodo release ingestion can fail on release metadata, so metadata should be
  validated before release.
- DataCite citable metadata requires identifier, creator, title, publisher,
  publication year, and resource type.

## Claim Boundary

Archive release is deferred. The current outputs remain diagnostic metadata
only: no archive identifier, no DOI, no frozen release version, no local license,
no creator/publisher metadata, no full OSDR payload mirror, and no descriptor
mutation are claimed.

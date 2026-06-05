---
title: SpaceBio-Bench Release Readiness Card
page_type: release_readiness_card
status: public_review_ready
last_reviewed: 2026-06-04
claim_boundary: benchmark_release_readiness_card_draft_no_release_approval
---

# SpaceBio-Bench Release Readiness Card

## Card Purpose

This card defines release tiers, readiness gates, blocked claims, and required
evidence for SpaceBio-Bench artifacts. It is designed to keep documentation,
metadata-alpha scaffolds, diagnostic lanes, frozen payload releases, and
archive-ready releases separate without making a stronger release claim than
the evidence supports.

This card does not approve a new release. It records the conditions a surface
must satisfy before its public wording can become stronger.

Branch note: on the default `main` branch, v9-specific evidence paths such as
`v9/...` and `docs/V9_*` refer to the curated public bulk metadata-alpha subset
included in this repository. Payload matrices and draft extension lanes are
excluded from the current public-review path.

## Release Tiers

| Tier | Meaning | Current examples | Public claim status |
|---|---|---|---|
| `historical_result_surface` | Existing result surface with canonical documentation | v1-v7 GeneLab Benchmark; v7.1 canonical result surface; v7.1.2 documentation/card/evidence-visibility patch | Allowed when using canonical v7.1 result language and v7.1.2 patch qualifiers |
| `metadata_alpha` | Public metadata, provenance, task, and baseline scaffold without frozen payload mirror | v9 public bulk alpha | Allowed only with metadata-only and no-payload qualifiers |
| `diagnostic_alpha` | Draft lane used for feasibility, metric, payload, or asset diagnostics | v9 single-cell, organoid, multispecies draft lanes | Draft-only; no public leaderboard or frozen release claim |
| `frozen_payload_release` | Future payload bundle with local mirror and verified payload hashes | Not yet active for v9 public bulk | Blocked until payload gates pass |
| `doi_archive_release` | Future citable archive release with final metadata, license, citations, and integrity manifest | Not yet active | Blocked until archive gates pass |
| `excluded` | Artifact or claim that should not enter public benchmark release | controlled human sequence data; unsupported countermeasure claims; Mars-regime predictions | Not releasable under current scope |

## Current Readiness State

| Area | Status | Evidence |
|---|---|---|
| v7.1 canonical result surface | Active as public result source of truth | `docs/CANONICAL_RESULTS_V7_1.md` |
| v9 public bulk task registry | Pass for metadata alpha | `v9/task_manifest_index.csv` |
| v9 public bulk fold registry | Pass for metadata alpha | `v9/task_data_index.csv` |
| v9 public bulk source inventory | Pass for metadata alpha | `v9/source_inventory.csv` |
| OSDR checksum-manifest evidence | Pass for metadata alpha | `v9/source_checksum_audit.csv` |
| Payload-level local hash verification | Blocker for frozen payload release | `freeze_ready_source_count=0` in snapshot decision |
| Baseline output evidence | Pass for metadata alpha | `v9/reports/bulk_lomo_baseline_summary.csv` |
| Dataset-card alpha boundary | Applied | `docs/v9_hf_dataset_card.md`; `docs/V9_PUBLIC_BULK_ALPHA_CARD_DATAPACKAGE_BOUNDARY_UPDATE.md` |
| Draft Data Package alpha boundary | Applied | `v9/datapackage.draft.json`; `docs/V9_PUBLIC_BULK_ALPHA_CARD_DATAPACKAGE_BOUNDARY_UPDATE.md` |
| Paper archive metadata | Release-candidate ready | `CITATION.cff`; `.zenodo.json`; `docs/RELEASE_ARCHIVE_CARD.md`; `docs/RELEASE_ARCHIVE_MANIFEST.md` |

## Metadata Alpha Requirements

A `metadata_alpha` release surface may be described publicly only when:

- Task manifests are generated and indexed.
- Fold definitions are indexed.
- Source inventory is present.
- Source provenance or checksum-manifest evidence is present.
- Baseline outputs, if reported, have predictions, metrics, and run manifests.
- Dataset card states metadata-only alpha status.
- Draft Data Package or equivalent descriptor states release status and payload
  boundary.
- Payload-release blockers are explicit.
- Draft extension lanes are excluded from the core alpha scope unless their
  own claim boundary is named.

The current v9 public bulk lane satisfies this tier after
`V9-BULK-ALPHA-003`, with the boundary
`metadata_only_public_bulk_alpha_no_payload_release`.

## Frozen Payload Release Requirements

A future `frozen_payload_release` requires:

- Local payload mirror for release-target fold matrices.
- Payload-level SHA-256 manifest for every distributed payload file.
- Verification report showing all expected payload hashes pass.
- Release `datapackage.json`, not only `datapackage.draft.json`.
- Dataset card language updated from metadata alpha to frozen payload release.
- License and reuse status reviewed.
- Dataset-specific OSDR citations filled from OSDR study pages.
- Regression checks for row counts, fold definitions, selected-gene counts,
  source ids, and no local path leaks.

Until these gates pass, avoid all frozen-payload wording.

## DOI Or Archive Release Requirements

A future `doi_archive_release` requires everything in
`frozen_payload_release`, plus:

- Final release manifest.
- DataCite-aligned metadata: title, creators, contributors, version,
  identifiers, related identifiers, resource type, license, and descriptions.
- RO-Crate or equivalent research-object metadata linking data, code, task
  manifests, reports, and provenance.
- Archive deposit target selected, such as Zenodo or another approved
  repository.
- Citation and acknowledgement text reviewed for OSDR and individual upstream
  datasets.
- Changelog entry and maintenance policy.

## Diagnostic Alpha Requirements

`diagnostic_alpha` lanes may be useful internally or for reviewer discussion,
but each must keep its own boundary. A diagnostic lane should include:

- A named claim boundary.
- Source inventory or asset inventory.
- Metric or payload contract if scores are discussed.
- Explicit statement that the lane is not a leaderboard.
- Explicit statement that the lane is not a frozen payload release.
- Owner action or blocker list for promotion.

Examples include single-cell asset inventory, organoid diagnostic surfaces, and
multispecies feasibility or interaction-task scaffolds.

## Excluded Or Blocked Claims

Do not promote these claims under the current readiness state:

| Claim | Reason |
|---|---|
| Frozen v9 public benchmark release | Payload verification has not passed |
| Frozen v9 payload mirror | Local payload mirror and hash manifest are pending |
| DOI/archive-ready release | Final metadata, license, citation, and archive manifest are pending |
| State-of-the-art leaderboard | Baselines are scaffold anchors, not tuned endpoints |
| Uniform foundation-model leaderboard | Current FM rows are mixed-surface or adapter-dependent |
| Countermeasure recommendation | Current evidence is not intervention validation |
| Crew-health or clinical decision support | Current scope excludes clinical recommendations |
| Mars-regime prediction | Current tasks are not validated for Mars-regime extrapolation |

## Readiness Rules

- Move a surface to a stronger tier only when its evidence files pass the tier
  requirements.
- Preserve the old tier label in changelogs when a surface moves tiers.
- Do not use a stronger tier label in slide, paper, README, or dataset-card
  text before the readiness card and claim register are updated.
- When a count changes, update the source artifact, dataset card, system card,
  and claim register in the same edit batch.
- If an artifact is diagnostic-only, its filename or report header should say
  so.

## Correction Notes

If a release-facing artifact overstates the current boundary:

1. Correct the artifact.
2. Add a note in the relevant decision or readiness document.
3. Update the claim register if the allowed or blocked wording changed.
4. If external users may have seen the artifact, add a public correction note
   in the release manifest or changelog for that surface.

## Companion Documents

- `docs/SPACEBIOBENCH_SYSTEM_CARD.md`
- `docs/SPACEBIOBENCH_EVALUATION_CARD.md`
- `docs/SPACEBIOBENCH_CLAIM_REGISTER.md`
- `docs/RELEASE_ARCHIVE_CARD.md`
- `docs/RELEASE_ARCHIVE_MANIFEST.md`
- `docs/RELEASE_ARCHIVE_CHECKLIST.md`
- `docs/CANONICAL_RESULTS_V7_1.md`
- `docs/v9_hf_dataset_card.md`
- `docs/V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md`
- `docs/V9_PUBLIC_BULK_ALPHA_CARD_DATAPACKAGE_BOUNDARY_UPDATE.md`
- `v9/reports/public_bulk_alpha_gap_matrix/public_bulk_alpha_gap_summary.csv`
- `v9/reports/public_bulk_alpha_snapshot_decision/snapshot_decision_summary.csv`
- `v9/datapackage.draft.json`

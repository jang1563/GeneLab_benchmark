---
title: GeneLab Benchmark Release Archive Checklist
page_type: release_archive_checklist
status: release_candidate_checklist
last_reviewed: 2026-06-04
claim_boundary: checklist_no_new_release_or_doi_claim
---

# GeneLab Benchmark Release Archive Checklist

## Current Status

This checklist separates what is ready for a paper-supporting archive from what
must wait until the final manuscript, tag, and DOI are available.

## Ready Now

- Root `README.md` presents the project, canonical result doc, transparency
  card pack, portfolio brief, and HF dataset link.
- `CITATION.cff` exists and is machine-readable citation metadata.
- `.zenodo.json` exists as candidate Zenodo metadata.
- `LICENSE` is present and declares MIT for code.
- `docs/hf_dataset_card.md` is synchronized with the public HF card surface.
- SpaceBio-Bench system, evaluation, readiness, claim, and transparency cards
  are `public_review_ready`.
- `docs/RELEASE_ARCHIVE_CARD.md` documents archive scope and exclusions.
- `docs/RELEASE_ARCHIVE_MANIFEST.md` documents archive contents and final tag
  procedure.
- CI has passed for the recent public-card merge PRs.

## Needs Final Author Or Manuscript Review

- Confirm author order and affiliations in `CITATION.cff`.
- Confirm creator metadata in `.zenodo.json`.
- Confirm final paper title and whether the archive title should match exactly.
- Confirm final version string and release date.
- Add DOI after Zenodo or the selected archive mints it.
- Add dataset-specific OSDR citations for the manuscript subset.

## Needs Final Release Action

- Create annotated Git tag.
- Create GitHub release from the tag.
- Verify Zenodo deposition metadata.
- Download source archive and record SHA-256 checksum.
- Store release DOI, tag, release URL, and checksum in a release-specific
  manifest.
- Update README, HF card, and citation metadata with DOI if appropriate.

## Still Blocked For Stronger Claims

- Frozen v9 public bulk payload release.
- DOI/archive-ready v9 metadata-alpha release.
- Complete RO-Crate research-object release.
- Clinical, crew-health, intervention, countermeasure, or Mars-regime claims.
- Foundation-model leaderboard claims without matched adapter/run manifests.

## Red-Flag Checks To Re-Run

```bash
git diff --check
git ls-files | rg '(^output/|^outputs/|__pycache__|\\.DS_Store$|\\.log$|\\.tmp$|\\.bak$|~$|\\.zip$|\\.tar$|\\.tar\\.gz$|\\.gz$|\\.pdf$|\\.pptx$|\\.png$|\\.jpg$|\\.jpeg$|\\.gif$)'
rg --glob '!docs/RELEASE_ARCHIVE_CHECKLIST.md' -n '(gh[pousr]_[A-Za-z0-9_]{20,}|sk-[A-Za-z0-9]{20,}|hf_[A-Za-z0-9]{20,}|/Users/jak4013|~/.claude)' README.md docs scripts tests .github
rg --glob '!docs/RELEASE_ARCHIVE_CHECKLIST.md' -n 'status: draft|draft_public_ready|blob/codex|codex/spacebiobench' README.md docs
```

Expected interpretation:

- Generated or binary artifacts should not appear unless intentionally archived.
- Secret scans should not show literal credentials.
- Local personal paths should appear only in tests that explicitly forbid them.
- Public cards should remain `public_review_ready`.

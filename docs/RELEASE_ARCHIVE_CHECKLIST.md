---
title: SpaceBio-Bench Release Archive Checklist
page_type: release_archive_checklist
status: public_ready
last_reviewed: 2026-06-16
---

# SpaceBio-Bench Release Archive Checklist

## Current Status

This checklist separates repository materials that are ready now from final
metadata actions that should happen near DOI or archive deposition.

## Ready Now

- Root `README.md` presents the project, canonical result document, public
  documentation map, and HF dataset link.
- `CITATION.cff` is machine-readable citation metadata using the current public
  repository identity.
- `.zenodo.json` provides archive metadata using the current public repository
  identity.
- `LICENSE` declares MIT for code.
- `docs/hf_dataset_card.md` is synchronized with the public HF card surface.
- SpaceBio-Bench public cards are marked `public_ready`.
- `docs/RELEASE_ARCHIVE_CARD.md` summarizes archive scope.
- `docs/RELEASE_ARCHIVE_MANIFEST.md` lists archive file groups and final tag
  procedure.
- CI and public release QA pass for the current public docs.

## Final Metadata Review

- Confirm manuscript author list and affiliations in `CITATION.cff`.
- Confirm creator metadata in `.zenodo.json`.
- Confirm final paper title and archive title.
- Confirm final version string and release date.
- Add DOI after Zenodo or the selected archive service mints it.
- Add dataset-specific OSDR citations for the manuscript subset.

## Final Release Actions

- Confirm Git tag `v7.1.2` or a successor tag from the intended public commit.
- Create or update the GitHub release from the tag.
- Verify Zenodo deposition metadata.
- Download source archive and record SHA-256 checksum.
- Store release DOI, tag, release URL, and checksum in a release-specific
  manifest.
- Update README, HF card, and citation metadata with DOI if appropriate.

## Public QA Commands

```bash
git diff --check
make release-qa
make hpc-public-qa
```

Optional archive hygiene checks:

```bash
git ls-files | rg '(^output/|^outputs/|__pycache__|\\.DS_Store$|\\.log$|\\.tmp$|\\.bak$|~$|\\.zip$|\\.tar$|\\.tar\\.gz$|\\.gz$|\\.pdf$|\\.pptx$|\\.jpg$|\\.jpeg$|\\.gif$)'
rg -n '(gh[pousr]_[A-Za-z0-9_]{20,}|sk-[A-Za-z0-9]{20,}|hf_[A-Za-z0-9]{20,})' README.md docs scripts tests .github
```

Expected interpretation:

- Generated or binary artifacts should appear only when intentionally tracked.
- Secret scans are expected to show no literal credentials.
- Public cards should remain `public_ready`.

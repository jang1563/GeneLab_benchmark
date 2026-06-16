---
title: SpaceBio-Bench Release Archive Manifest
page_type: release_archive_manifest
status: public_ready
last_reviewed: 2026-06-16
---

# SpaceBio-Bench Release Archive Manifest

## Scope

This manifest lists the repository materials that support the current public
archive surface for SpaceBio-Bench / GeneLab Benchmark.

## Current Repository Inventory

| Field | Value |
|---|---|
| Public release tag | `v7.1.2` |
| Repository branch | `main` |
| Root citation file | `CITATION.cff` |
| Zenodo metadata file | `.zenodo.json` |
| License | `LICENSE` |
| HF card source | `docs/hf_dataset_card.md` |
| Release archive status | `docs/RELEASE_ARCHIVE_CARD.md` |
| Release archive checklist | `docs/RELEASE_ARCHIVE_CHECKLIST.md` |
| v9 metadata catalog | `v9/`; `docs/v9_hf_dataset_card.md` |

## Archive File Groups

| Class | Files |
|---|---|
| Entry point | `README.md` |
| Citation and license | `CITATION.cff`; `.zenodo.json`; `LICENSE` |
| Canonical results | `docs/CANONICAL_RESULTS_V7_1.md`; `docs/hf_dataset_card.md` |
| Public documentation map | `docs/SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md` |
| Public cards | `docs/SPACEBIOBENCH_SYSTEM_CARD.md`; `docs/SPACEBIOBENCH_EVALUATION_CARD.md`; `docs/SPACEBIOBENCH_RELEASE_READINESS_CARD.md`; `docs/SPACEBIOBENCH_CLAIM_REGISTER.md` |
| Archive status docs | `docs/RELEASE_ARCHIVE_CARD.md`; `docs/RELEASE_ARCHIVE_MANIFEST.md`; `docs/RELEASE_ARCHIVE_CHECKLIST.md` |
| v9 metadata catalog | `v9/`; `docs/v9_hf_dataset_card.md` |
| Reproducibility scripts | `scripts/`; `.github/workflows/` |
| Public task metadata | `tasks/` |
| Evaluation outputs | `evaluation/`; `v2/`; `v3/`; `v4/`; `v5/`; `v6/`; `v7/` |
| Data catalog helpers | `DATA_CATALOG.md`; `GLDS_verified.json`; `data/DOWNLOAD_LOG.json` |

## Data Package Notes

The GitHub repository is not the complete raw-data archive. Large source data
and upstream biological payloads remain with NASA OSDR and the individual study
pages.

Hugging Face hosts the processed public fold package for selected benchmark
tasks. The HF dataset card should remain synchronized with
`docs/hf_dataset_card.md`.

The v9 public bulk surface is a metadata catalog. It records task manifests,
source rows, fold indexes, checksum-audit summaries, and reference baseline
outputs.

## Suggested Final Archive Procedure

1. Confirm manuscript title, author list, affiliations, and release version.
2. Update `CITATION.cff` and `.zenodo.json` with final metadata.
3. Confirm Git tag `v7.1.2` or a successor tag from the intended public commit.
4. Create or update the GitHub release from that tag.
5. Let Zenodo mint the DOI from the GitHub release, or deposit manually with
   `.zenodo.json` metadata.
6. Download the generated source archive and compute SHA-256.
7. Record DOI, tag, release URL, and archive checksum in a release-specific
   manifest.
8. Re-run public release QA before citing the archive.

## Checksum Commands

Use these commands after the final release tag exists:

```bash
git archive --format=tar.gz --prefix=GeneLab_benchmark-v7.1.2/ \
  -o GeneLab_benchmark-v7.1.2.tar.gz v7.1.2
shasum -a 256 GeneLab_benchmark-v7.1.2.tar.gz
```

For a source checkout:

```bash
git rev-parse HEAD
git status --short
git ls-files | wc -l
du -sh .
```

## Notes For Final Review

- Review final manuscript author-list metadata before deposition.
- Add dataset-specific OSDR citations for the manuscript subset.
- Keep the v9 metadata catalog separate from v7.1 result statements unless the
  manuscript explicitly discusses it as a metadata resource.

---
title: GeneLab Benchmark Release Archive Manifest
page_type: release_archive_manifest
status: release_candidate_manifest
last_reviewed: 2026-06-05
claim_boundary: manifest_no_new_release_or_payload_claim
---

# GeneLab Benchmark Release Archive Manifest

## Scope

This manifest describes the paper-supporting release-candidate archive surface
for GeneLab Benchmark on the default `main` branch.

It is a manifest for repository archive readiness, not a DOI certificate and
not a frozen v9 payload manifest.

## Current Repository Inventory

| Field | Value |
|---|---:|
| Audited repository scope | Default `main` release tree plus curated v9 metadata-alpha subset |
| Curated v9 subset size | ~2 MB; metadata/provenance/baseline evidence only |
| Root citation file | `CITATION.cff` |
| Zenodo metadata file | `.zenodo.json` |
| License | `LICENSE` |
| HF card source | `docs/hf_dataset_card.md` |
| Release archive card | `docs/RELEASE_ARCHIVE_CARD.md` |
| Release archive checklist | `docs/RELEASE_ARCHIVE_CHECKLIST.md` |

## Critical Archive Files

| Class | Files |
|---|---|
| Entry point | `README.md` |
| Citation and license | `CITATION.cff`; `.zenodo.json`; `LICENSE` |
| Canonical results | `docs/CANONICAL_RESULTS_V7_1.md`; `docs/hf_dataset_card.md` |
| Transparency entry | `docs/SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md` |
| System/evaluation cards | `docs/SPACEBIOBENCH_SYSTEM_CARD.md`; `docs/SPACEBIOBENCH_EVALUATION_CARD.md` |
| Readiness/claim cards | `docs/SPACEBIOBENCH_RELEASE_READINESS_CARD.md`; `docs/SPACEBIOBENCH_CLAIM_REGISTER.md` |
| Archive controls | `docs/RELEASE_ARCHIVE_CARD.md`; `docs/RELEASE_ARCHIVE_MANIFEST.md`; `docs/RELEASE_ARCHIVE_CHECKLIST.md` |
| v9 metadata-alpha evidence | `v9/`; `docs/V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md`; `docs/V9_PUBLIC_BULK_ALPHA_CARD_DATAPACKAGE_BOUNDARY_UPDATE.md` |
| Reproducibility scripts | `scripts/`; `.github/workflows/` |
| Public task metadata | `tasks/` |
| Evaluation outputs | `evaluation/`; `v2/`; `v3/`; `v4/`; `v5/`; `v6/`; `v7/` |
| Data download/provenance helpers | `DATA_CATALOG.md`; `GLDS_verified.json`; `data/DOWNLOAD_LOG.json` |

## Payload Boundary

The GitHub repository is not the complete raw-data archive. Large source data
and upstream biological payloads remain governed by NASA OSDR and individual
study pages.

Hugging Face hosts the public feature-matrix package for selected benchmark
tasks. The HF dataset card should remain synchronized with `docs/hf_dataset_card.md`.

The curated v9 public bulk surface remains metadata-alpha only. It should not
be cited as a frozen payload archive until payload-level hashes and release
manifests are complete.

## Suggested Final Archive Procedure

1. Confirm final manuscript title, final author list, affiliations, and release
   version.
2. Update `CITATION.cff` and `.zenodo.json` with final metadata.
3. Create or confirm annotated Git tag `v7.1.2` from the final
   public-card/archive-metadata commit. The older `v7.1` GitHub release
   predates the card-pack and archive-metadata files.
4. Create a GitHub release from that tag.
5. Let Zenodo mint the DOI from the GitHub release, or deposit manually with
   `.zenodo.json` metadata.
6. Download the generated source archive and compute SHA-256.
7. Record the DOI, tag, release URL, and archive checksum in this manifest or a
   release-specific manifest.
8. Re-run CI and the public-card red-flag checks before citing the archive.

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

## Known Archive Caveats

- The paper archive DOI has not yet been minted.
- The older `v7.1` GitHub release predates the current card-pack and
  archive-metadata polish; cite `v7.1.2` or a successor tag for this card
  surface once minted.
- Final manuscript author-list metadata should be reviewed before final
  deposition.
- NASA OSDR source datasets need dataset-specific citation review for the
  manuscript subset.
- v9 metadata-alpha evidence should remain separate from v1-v7 final result
  claims unless the manuscript explicitly scopes it as future or supplemental
  work.

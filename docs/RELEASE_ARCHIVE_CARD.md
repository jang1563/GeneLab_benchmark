---
title: SpaceBio-Bench Release Archive Status
page_type: release_archive_status
status: public_ready
last_reviewed: 2026-06-16
---

# SpaceBio-Bench Release Archive Status

## Purpose

This card summarizes the repository materials that support a paper or archive
release for SpaceBio-Bench / GeneLab Benchmark. It is a planning aid for
repository packaging, citation metadata, and dataset links.

## Current Archive Scope

| Field | Current value |
|---|---|
| Public release line | `v7.1.2` documentation, public-card, citation, and metadata patch |
| Result surface | v7.1 canonical historical benchmark results |
| Repository entry point | `main` branch and `v7.1.2` tag |
| Dataset access | Hugging Face public fold package plus GitHub task metadata |
| v9 public bulk | Metadata catalog for task/source/fold/audit/baseline records |
| DOI status | Ready for final metadata review before archive deposition |

## Archive Contents

The repository archive should include:

- source code and evaluation scripts in `scripts/`;
- public task metadata, labels, and fold definitions in `tasks/`;
- result JSON and summary artifacts in `evaluation/` and versioned evaluation
  directories;
- public documentation in `README.md`, `docs/`, and `CITATION.cff`;
- release metadata in `.zenodo.json`;
- Hugging Face dataset card source in `docs/hf_dataset_card.md`;
- SpaceBio-Bench public cards in `docs/SPACEBIOBENCH_*.md`;
- v9 public bulk metadata catalog files under `v9/`.

## Scope Notes

The repository archive does not bundle raw NASA OSDR source payloads or
controlled-access human sequence data. Source biological data remain governed
by NASA OSDR and the individual OSDR study pages.

The v9 public bulk surface is documented as a metadata catalog. Larger payload
bundles, DOI-specific archive records, or research-object packaging can be
added as separate release work when the relevant metadata and checks are ready.

## Before DOI Deposition

- Confirm manuscript title, author list, affiliations, and release date.
- Review `CITATION.cff` and `.zenodo.json` against the final manuscript.
- Confirm the intended GitHub tag and release URL.
- Confirm the Hugging Face dataset card links to the intended repository tag.
- Record dataset-specific OSDR citations for the manuscript subset.
- Generate a source archive checksum after selecting the final release archive.

## Current Status

| Area | Status | Public file |
|---|---|---|
| README entry point | Ready | `README.md` |
| Citation metadata | Ready for final manuscript review | `CITATION.cff` |
| Zenodo metadata | Ready for final manuscript review | `.zenodo.json` |
| Code license | Ready | `LICENSE` |
| HF dataset card | Ready | `docs/hf_dataset_card.md` |
| Public documentation map | Ready | `docs/SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md` |
| Release manifest | Ready | `release/release_manifest.json` |
| Archive manifest | Ready | `docs/RELEASE_ARCHIVE_MANIFEST.md` |
| Source archive checksum | To be generated after archive selection | release-specific checksum record |
| DOI | To be minted by Zenodo or the selected archive service | archive record |

## References

- GitHub citation files:
  https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files
- Citation File Format:
  https://citation-file-format.github.io/
- Zenodo GitHub integration:
  https://help.zenodo.org/docs/deposit/github/
- DataCite Metadata Schema:
  https://schema.datacite.org/
- RO-Crate:
  https://www.researchobject.org/ro-crate/

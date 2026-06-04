---
title: GeneLab Benchmark Release Archive Card
page_type: release_archive_card
status: release_candidate_metadata_ready
last_reviewed: 2026-06-04
claim_boundary: archive_metadata_no_new_result_or_payload_claim
---

# GeneLab Benchmark Release Archive Card

## Card Purpose

This card defines what would be archived for a paper-supporting release of
GeneLab Benchmark, what is intentionally excluded, and which gates remain before
minting a final DOI or claiming a frozen research-object release.

It does not create a DOI, approve a new result release, or freeze v9 payloads.

## Archive Candidate Summary

| Field | Current value |
|---|---|
| Release candidate | v7.1 public documentation and card-pack surface |
| Repository branch | `main` as public entry point; `v3` as detailed v9 evidence branch |
| Current result boundary | v1-v7 canonical historical benchmark surface |
| Documentation patch | v7.1 consistency patch, no new result generation |
| Dataset payload boundary | HF feature-matrix package plus GitHub task/evaluation metadata |
| DOI status | Not minted in this card |
| Archive status | Release-candidate metadata ready; final tag/DOI still pending |

## Archive Contents

The archive candidate should include:

- Source code and evaluation scripts in `scripts/`.
- Public task metadata, labels, and fold definitions in `tasks/`.
- Public result JSON and summary artifacts in `evaluation/` and versioned
  evaluation directories.
- Human-facing release documentation in `README.md`, `docs/`, and
  `CITATION.cff`.
- Release metadata in `.zenodo.json`.
- Hugging Face dataset card text in `docs/hf_dataset_card.md`.
- SpaceBio-Bench transparency cards in `docs/SPACEBIOBENCH_*.md`.

## Excluded From The Archive Claim

The archive candidate should not claim to include:

- Raw NASA OSDR source payloads.
- Controlled-access human sequence data.
- Frozen v9 public bulk payload files.
- Locally hash-verified v9 payload bundle.
- DOI-ready RO-Crate research object.
- Clinical, crew-health, countermeasure, intervention, or Mars-regime guidance.

## Required Before DOI Or Final Release Tag

Before a DOI-oriented release:

- Confirm author order, affiliations, and manuscript title in `CITATION.cff`
  and `.zenodo.json`.
- Create an annotated GitHub release tag from the final commit.
- Generate and store a checksum for the release source archive.
- Confirm the Hugging Face dataset card links to the intended canonical branch
  and release tag.
- Record individual OSDR dataset citations for any analysis subset used in the
  manuscript.
- Decide whether v9 metadata-alpha artifacts are cited as future work, a
  supplemental surface, or excluded from the paper archive.
- If claiming research-object completeness, add RO-Crate or equivalent
  provenance metadata.

## Current Readiness

| Gate | Status | Note |
|---|---|---|
| Public README entry point | Pass | Links card pack and portfolio brief |
| Citation metadata | Partial pass | `CITATION.cff` exists; author/title review still needed |
| Zenodo metadata | Candidate ready | `.zenodo.json` added for DOI deposition |
| License | Pass for code | MIT license present |
| HF dataset card | Pass | Links transparency cards on canonical `v3` |
| Claim-boundary cards | Pass | System, evaluation, readiness, and claim cards are public-review ready |
| Release archive manifest | Pass | `docs/RELEASE_ARCHIVE_MANIFEST.md` added |
| Source archive checksum | Pending final tag | Generate after final GitHub release tag |
| DOI | Pending | Mint through Zenodo or equivalent archive after metadata review |
| v9 frozen payload | Blocked | Payload-level hash verification remains pending |

## External Standards Used

- GitHub citation files:
  https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files
- Citation File Format:
  https://citation-file-format.github.io/
- Zenodo GitHub integration:
  https://help.zenodo.org/docs/deposit/github/
- DataCite Metadata Schema 4.7:
  https://schema.datacite.org/meta/kernel-4.7/
- RO-Crate:
  https://www.researchobject.org/ro-crate/1.2/

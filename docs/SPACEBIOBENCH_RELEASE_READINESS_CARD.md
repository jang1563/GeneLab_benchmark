---
title: SpaceBio-Bench Release Status Card
page_type: release_status_card
status: public_ready
last_reviewed: 2026-06-16
---

# SpaceBio-Bench Release Status Card

## Purpose

This card summarizes the current public status of the SpaceBio-Bench release
surfaces and the additional work needed for larger dataset or archive releases.

## Current Public Status

| Surface | Public status | Current wording |
|---|---|---|
| v7.1 GeneLab Benchmark | Canonical historical result surface | Use for v1-v7 results and citation context |
| v7.1.2 public patch | Active documentation and metadata patch | Use for README, HF card, citation, Zenodo, and release metadata updates |
| Hugging Face public fold package | Active processed dataset package | Use for selected LOMO fold downloads and public task access |
| v9 public bulk | Metadata catalog | Use for task manifests, source records, fold indexes, audit summaries, and reference baseline rows |
| Extension workspaces | Development workspaces | Keep separate from the v7.1 public result surface |

## Active Release Checks

| Area | Public file |
|---|---|
| v7.1 result source | `docs/CANONICAL_RESULTS_V7_1.md` |
| v7.1 dataset card source | `docs/hf_dataset_card.md` |
| Citation metadata | `CITATION.cff`; `.zenodo.json` |
| Release manifest | `release/release_manifest.json` |
| v9 task registry | `v9/task_manifest_index.csv` |
| v9 fold registry | `v9/task_data_index.csv` |
| v9 source inventory | `v9/source_inventory.csv` |
| v9 checksum audit | `v9/source_checksum_audit.csv` |
| v9 baseline summaries | `v9/reports/bulk_lomo_baseline_summary.csv` |
| v9 metadata catalog card | `docs/v9_hf_dataset_card.md` |

## For A Larger v9 Payload Release

A larger v9 payload release would need:

- release-target fold matrices packaged as a public payload bundle;
- payload-level SHA-256 manifest for distributed files;
- verification report for expected payload hashes;
- release `datapackage.json` or equivalent descriptor;
- reviewed license and source reuse language;
- dataset-specific OSDR citations;
- regression checks for row counts, fold definitions, selected-gene counts,
  source IDs, and local-path hygiene.

## For A DOI Or Archive Release

A DOI or archive-oriented release would additionally need:

- final release manifest;
- DataCite-aligned metadata;
- related identifiers for repository, dataset, and archive records;
- final resource type, license, contributors, and descriptions;
- research-object metadata such as RO-Crate when appropriate;
- archive deposit target such as Zenodo;
- reviewed citation and acknowledgment text.

## Public Reporting Rules

- Keep v7.1 result statements tied to the canonical v7.1 result document.
- Describe v7.1.2 as a documentation, card, citation, and metadata patch.
- Describe v9 public bulk as a metadata catalog.
- Pair pooled metrics with per-task or per-fold context.
- Update README, HF card, citation, release manifest, and public docs together
  when release labels or counts change.

## Companion Documents

- [System card](SPACEBIOBENCH_SYSTEM_CARD.md)
- [Evaluation card](SPACEBIOBENCH_EVALUATION_CARD.md)
- [Public statement guide](SPACEBIOBENCH_CLAIM_REGISTER.md)
- [Canonical v7.1 results](CANONICAL_RESULTS_V7_1.md)
- [v9 metadata catalog card](v9_hf_dataset_card.md)

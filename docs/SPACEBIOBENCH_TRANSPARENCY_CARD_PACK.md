---
title: SpaceBio-Bench Public Documentation Map
page_type: public_documentation_map
status: public_ready
last_reviewed: 2026-06-16
---

# SpaceBio-Bench Public Documentation Map

This map points readers to the main public documentation for SpaceBio-Bench:
system scope, evaluation interpretation, release status, and recommended public
wording. It is a navigation layer, not a new benchmark result surface.

## Current Public Summary

- **v7.1 GeneLab Benchmark** is the canonical historical result surface for
  v1-v7 results.
- **v7.1.2** is the public-card and metadata patch over v7.1. It updates
  documentation, citation metadata, release metadata, and access guidance
  without adding new result generation.
- **Hugging Face dataset** provides processed public fold packages for selected
  mission-held-out tasks.
- **v9 public bulk** is a metadata catalog for public bulk RNA-seq task
  manifests, source records, fold indexes, checksum-audit summaries, and
  reference baselines.
- Extension workspaces for single-cell, organoid, multispecies, and
  translational analyses remain separate from the v7.1 public result surface.

## Three-Minute Review Map

| Step | Open this | What it gives you |
|---|---|---|
| 1 | [System card](SPACEBIOBENCH_SYSTEM_CARD.md) | Benchmark surfaces, components, intended use, and scope notes |
| 2 | [Evaluation card](SPACEBIOBENCH_EVALUATION_CARD.md) | How to read tasks, folds, metrics, baselines, and pooled summaries |
| 3 | [Benchmark integrity note](BENCHMARK_INTEGRITY.md) | Which labels are public and what can or cannot be called blind evaluation |
| 4 | [Release status card](SPACEBIOBENCH_RELEASE_READINESS_CARD.md) | Current public status for v7.1, v7.1.2, HF, and v9 catalog surfaces |
| 5 | [Public statement guide](SPACEBIOBENCH_CLAIM_REGISTER.md) | Preferred wording for common public statements |
| Cross-check | [Canonical v7.1 results](CANONICAL_RESULTS_V7_1.md), [v7.1 HF card](hf_dataset_card.md), and [v9 metadata catalog card](v9_hf_dataset_card.md) | The source document for a result, dataset, or catalog statement |

## Card Index

| Card | Question | File |
|---|---|---|
| System Card | What is SpaceBio-Bench and what surfaces exist? | [docs/SPACEBIOBENCH_SYSTEM_CARD.md](SPACEBIOBENCH_SYSTEM_CARD.md) |
| Evaluation Card | How should task, fold, baseline, metric, and pooled results be interpreted? | [docs/SPACEBIOBENCH_EVALUATION_CARD.md](SPACEBIOBENCH_EVALUATION_CARD.md) |
| Benchmark Integrity | Which evaluation labels are public, and what validation claims are supported? | [docs/BENCHMARK_INTEGRITY.md](BENCHMARK_INTEGRITY.md) |
| Release Status Card | What is active now, and what would be needed for a larger release? | [docs/SPACEBIOBENCH_RELEASE_READINESS_CARD.md](SPACEBIOBENCH_RELEASE_READINESS_CARD.md) |
| Public Statement Guide | What wording should be used for public descriptions? | [docs/SPACEBIOBENCH_CLAIM_REGISTER.md](SPACEBIOBENCH_CLAIM_REGISTER.md) |
| v7.1 Canonical Results | What is the locked public result and scope source for v1-v7? | [docs/CANONICAL_RESULTS_V7_1.md](CANONICAL_RESULTS_V7_1.md) |
| v7.1 HF Dataset Card | What is in the public processed fold package? | [docs/hf_dataset_card.md](hf_dataset_card.md) |
| v9 Metadata Catalog Card | What does the v9 public bulk metadata catalog contain? | [docs/v9_hf_dataset_card.md](v9_hf_dataset_card.md) |

## Public Update Checklist

- Keep the root README, Hugging Face card, citation metadata, and release
  manifest aligned on release labels and titles.
- Keep v7.1 result statements tied to
  [CANONICAL_RESULTS_V7_1.md](CANONICAL_RESULTS_V7_1.md).
- Keep v9 catalog statements tied to
  [v9_hf_dataset_card.md](v9_hf_dataset_card.md) and the files under `v9/`.
- Report per-task and per-fold metrics alongside pooled summaries.
- Cite NASA OSDR and the individual OSDR datasets used in downstream analyses.

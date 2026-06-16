---
title: SpaceBio-Bench System Card
page_type: system_card
status: public_ready
last_reviewed: 2026-06-16
---

# SpaceBio-Bench System Card

## Purpose

This card describes SpaceBio-Bench as a benchmark system. It summarizes the
public result surfaces, task structure, data sources, evaluation logic, and
scope notes that help readers interpret the repository.

## System Summary

SpaceBio-Bench evaluates whether AI/ML and foundation-model methods generalize
spaceflight transcriptomic signatures across missions. The core design is
mission-held-out validation: train on one set of missions, hold out a mission,
and test whether the method recognizes the flight-vs-ground signal in the
held-out context.

## Public Surfaces

| Surface | Current public role | Entry point |
|---|---|---|
| v7.1 GeneLab Benchmark | Canonical historical result surface for v1-v7 | [CANONICAL_RESULTS_V7_1.md](CANONICAL_RESULTS_V7_1.md) |
| v7.1.2 public patch | Documentation, citation, metadata, and access update over v7.1 | [README.md](../README.md) |
| Hugging Face dataset | Processed public fold package for selected LOMO tasks | [hf_dataset_card.md](hf_dataset_card.md) |
| v9 public bulk | Metadata catalog for task manifests, source rows, fold indexes, audit summaries, and baseline outputs | [v9_hf_dataset_card.md](v9_hf_dataset_card.md) |
| Extension workspaces | Development areas for additional biological settings and modalities | versioned workspace docs under `docs/` and `v*/` |

## System Components

| Component | Role | Public files |
|---|---|---|
| Task manifests | Define tissues, missions, labels, feature namespaces, and metrics | `v9/task_manifests/*.json`; `v9/task_manifest_index.csv` |
| Fold data index | Tracks held-out missions, row counts, selected-gene counts, and fold paths | `v9/task_data_index.csv` |
| OSDR source inventory | Maps OSDR accessions, GLDS prefixes, tissues, missions, and access status | `v9/source_inventory.csv` |
| Checksum audit | Summarizes OSDR API and checksum-manifest parsing | `v9/source_checksum_audit.csv` |
| Data package descriptor | Lists metadata and output resources for the v9 catalog | `v9/datapackage.draft.json` |
| Baseline outputs | Provide reference workflow rows for the v9 catalog | `v9/reports/bulk_lomo_baseline_summary.csv` |
| Dataset cards | Explain public dataset and metadata-catalog surfaces | `docs/hf_dataset_card.md`; `docs/v9_hf_dataset_card.md` |
| Canonical result document | Summarizes v7.1 public result scope and headline results | `docs/CANONICAL_RESULTS_V7_1.md` |

## Intended Uses

SpaceBio-Bench is intended for:

- evaluating generalization under mission-held-out transcriptomics shift;
- comparing methods on fixed task, fold, and feature definitions;
- auditing which public OSDR accessions support each task;
- reporting task-level and fold-level behavior before pooled summaries;
- developing public benchmark packaging for OSDR-derived biological data.

## Scope Notes

The public benchmark surfaces focus on task definitions, processed public fold
packages, baseline outputs, and result interpretation. They do not provide
clinical recommendations, crew-health guidance, countermeasure selection,
dosing guidance, or Mars-regime prediction.

The v9 public bulk surface is a metadata catalog. It is separate from archived
v9 fold-matrix payload bundles and from development workspaces for other
modalities.

## Evaluation Philosophy

Mission shift is the main evaluation pressure. Pooled averages are useful for
orientation, but the benchmark should be read at the tissue, task, fold, and
mission level before drawing conclusions from aggregate metrics.

Foundation-model and adapter comparisons require matched task inputs,
adapters, tissues, and evaluation surfaces. The v7.1 result document records
the current public summary for these comparisons.

## Known Limitations

- Mouse bulk RNA-seq is not a complete representation of space biology.
- Mission labels can combine spaceflight exposure, vehicle, hardware, age,
  protocol, tissue handling, processing, and batch effects.
- Some task labels include analog or special mission labels such as MHU-2 and
  OSD-397.
- Bulk RNA-seq does not resolve cell-type-specific effects.
- Legacy processed fold matrices can reflect earlier preprocessing choices.
- Baseline rows are reference workflow comparisons rather than tuned leaderboard
  endpoints.

## Recommended Reporting

- Cite NASA OSDR and the individual upstream OSDR datasets used in an analysis.
- Report per-task and per-fold metrics alongside pooled summaries.
- State the release surface used for each result.
- Use the public statement guide for concise release, dataset, and result
  wording: [SPACEBIOBENCH_CLAIM_REGISTER.md](SPACEBIOBENCH_CLAIM_REGISTER.md).

## Companion Documents

- [Public documentation map](SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md)
- [Evaluation card](SPACEBIOBENCH_EVALUATION_CARD.md)
- [Release status card](SPACEBIOBENCH_RELEASE_READINESS_CARD.md)
- [Public statement guide](SPACEBIOBENCH_CLAIM_REGISTER.md)
- [Canonical v7.1 results](CANONICAL_RESULTS_V7_1.md)
- [v7.1 Hugging Face dataset card](hf_dataset_card.md)
- [v9 metadata catalog card](v9_hf_dataset_card.md)

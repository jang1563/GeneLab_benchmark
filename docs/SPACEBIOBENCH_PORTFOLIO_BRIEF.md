---
title: SpaceBio-Bench Portfolio Brief
page_type: portfolio_brief
status: public_ready
last_reviewed: 2026-06-16
---

# SpaceBio-Bench Portfolio Brief

## One-Sentence Summary

SpaceBio-Bench is a mission-held-out transcriptomics benchmark for public NASA
OSDR space-biology data, built to evaluate whether AI/ML and foundation-model
methods generalize biological signatures across missions.

## Public Artifact Links

| Artifact | Link | Review purpose |
|---|---|---|
| GitHub repository | [jang1563/GeneLab_benchmark](https://github.com/jang1563/GeneLab_benchmark) | Source tree, documentation, tests, and release history |
| Hugging Face dataset | [jang1563/genelab-benchmark](https://huggingface.co/datasets/jang1563/genelab-benchmark) | Public fold package and dataset card |
| Public documentation map | [SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md](SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md) | System, evaluation, release-status, and statement-guide navigation |

## What I Built

- A cross-mission benchmark surface for public space-biology transcriptomics,
  with canonical v1-v7 result documentation and a v7.1.2 public-card/metadata
  patch.
- A live Hugging Face dataset card and processed public fold package for
  selected leave-one-mission-out tasks.
- A v9 public bulk metadata catalog with 8 task manifests, 33 fold definitions,
  22 public OSDR source rows, 24 reference baseline runs, and a Frictionless
  Data Package descriptor.
- A public documentation map that connects system scope, evaluation
  interpretation, release status, and recommended public wording.

## Why It Matters

Space-biology transcriptomics is small-sample and mission-shifted. A model can
look strong in a pooled summary while failing on a specific tissue, mission, or
task variant. SpaceBio-Bench treats that as an evaluation-design problem: every
result is tied back to a task, fold, tissue, feature surface, and release
surface.

The project demonstrates how to build a benchmark that is technically useful
and scientifically legible: tasks are explicit, source rows are traceable,
baseline runs are documented, and public language stays aligned with the
underlying release surface.

## Technical Depth

| Area | Evidence |
|---|---|
| Benchmark design | Mission-held-out LOMO task manifests and fold definitions |
| Bioinformatics scope | Public OSDR-derived mouse bulk RNA-seq tasks |
| Evaluation engineering | Baseline predictions, metrics, run records, and pooled-summary interpretation |
| Data access | Hugging Face public fold package plus GitHub documentation |
| Metadata packaging | v9 source inventory, checksum audit, and Data Package descriptor |
| Scientific communication | Public system card, evaluation card, release status card, and statement guide |

## Public Release Discipline

The public release surface is intentionally separated into:

- v7.1 canonical historical results;
- v7.1.2 documentation, card, citation, and metadata updates;
- the Hugging Face processed fold package;
- the v9 public bulk metadata catalog;
- development workspaces for future benchmark extensions.

That separation keeps benchmark scores, dataset access, metadata catalogs, and
extension work easy to inspect without flattening them into one headline.

## Role-Relevant Signal

This project is most relevant to roles involving:

- AI evaluation and benchmark design;
- bioinformatics and space-biology data infrastructure;
- scientific ML and biological foundation-model evaluation;
- AI safety evaluation for biological-domain claims;
- research engineering for reproducible scientific artifacts;
- technical communication for biological AI benchmarks.

## Fast Review Path

For a quick review, read these in order:

1. [Public Documentation Map](SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md)
2. [System Card](SPACEBIOBENCH_SYSTEM_CARD.md)
3. [Evaluation Card](SPACEBIOBENCH_EVALUATION_CARD.md)
4. [Release Status Card](SPACEBIOBENCH_RELEASE_READINESS_CARD.md)
5. [Public Statement Guide](SPACEBIOBENCH_CLAIM_REGISTER.md)
6. [v9 Metadata Catalog Card](v9_hf_dataset_card.md)
7. [Canonical v7.1 Results](CANONICAL_RESULTS_V7_1.md)
8. [Live Hugging Face Dataset](https://huggingface.co/datasets/jang1563/genelab-benchmark)

## Resume-Ready Summary

Built SpaceBio-Bench, a mission-held-out space-biology transcriptomics
benchmark using public NASA OSDR data. The project defines task/fold manifests,
processed public fold packages, reference baselines, checksum-audit summaries,
Hugging Face dataset-card documentation, release metadata, and public
interpretation guides for evaluating biological AI methods under mission shift.

## Cover-Letter Variant

SpaceBio-Bench reflects the kind of evaluation work I want to keep doing:
building biological AI benchmarks that are technically useful, scientifically
careful, and easy for other researchers to inspect. Beyond model scores, the
project documents task definitions, source coverage, baseline interpretation,
dataset access, and release status so that small datasets, mission confounding,
and foundation-model comparisons are not compressed into a single headline
number.

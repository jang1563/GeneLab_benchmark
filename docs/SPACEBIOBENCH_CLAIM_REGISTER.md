---
title: SpaceBio-Bench Public Statement Guide
page_type: public_statement_guide
status: public_ready
last_reviewed: 2026-06-16
---

# SpaceBio-Bench Public Statement Guide

## Purpose

This guide gives concise, public-ready wording for common SpaceBio-Bench
statements. Use it when writing README text, dataset cards, release notes,
abstracts, talks, or issue responses.

## Preferred Wording

| Topic | Use this wording | Context |
|---|---|---|
| Project name | SpaceBio-Bench / GeneLab Benchmark | SpaceBio-Bench is the forward-looking platform name; GeneLab Benchmark is the historical public name |
| Core task | Mission-held-out spaceflight transcriptomics benchmark | The independence unit is the held-out mission |
| v7.1 result surface | Canonical historical result surface for v1-v7 | Use for published result summaries and citation context |
| v7.1.2 patch | Documentation, public-card, citation, and metadata patch over canonical v7.1 results | The patch does not add new benchmark result generation |
| Hugging Face package | Processed public fold package for selected LOMO tasks | Use for direct dataset downloads |
| v9 public bulk | Public bulk metadata catalog | Use for task manifests, source records, fold indexes, audit summaries, and reference baselines |
| Foundation-model result summary | Tested gene-expression foundation models underperform tuned classical baselines on small-n bulk RNA-seq mission shift | Keep tied to the documented v7.1 result surface |
| Classical baseline summary | PCA-LR is the strongest 8-tissue gene-level baseline in v4, with mean AUROC 0.776 | Use with v7.1 canonical result context |
| Held-out validation | Thymus RR-23 AUROC 0.905; skin RR-7 AUROC 0.885 | Use with the v7.1 held-out validation context |
| OSDR source data | Source data are derived from public NASA OSDR studies | Cite the relevant OSDR study pages |

## Result Interpretation

- Report task, fold, tissue, method, feature surface, and release surface.
- Pair pooled summaries with per-task or per-fold rows.
- Treat benchmark scores as evidence for the stated benchmark task.
- Keep dataset, model, and release labels aligned across README, HF card,
  citation metadata, and release notes.

## Scope Language

Use scope language that is precise and compact:

- "The public benchmark evaluates mission-held-out transcriptomics
  generalization."
- "The v9 public bulk surface is a metadata catalog."
- "The Hugging Face dataset provides processed public fold packages."
- "Clinical, crew-health, countermeasure, and operational-readiness use cases
  are outside the current benchmark scope."

## Citation And Acknowledgment

Use the GitHub `CITATION.cff` metadata for the software and benchmark citation.
For source data, cite NASA OSDR and the individual OSDR studies used in the
analysis.

Recommended OSDR acknowledgment:

> Data are courtesy of the NASA Open Science Data Repository.

OSDR resource citation:

Gebre S G, Scott R T, Saravia-Butler A M, Lopez D K, Sanders L M, and Costes S
V. 2024. NASA Open Science Data Repository: Open Science for Life in Space.
Nucleic Acids Research 53(D1): D1697-D1710.
https://doi.org/10.1093/nar/gkae1116

## Companion Documents

- [Public documentation map](SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md)
- [System card](SPACEBIOBENCH_SYSTEM_CARD.md)
- [Evaluation card](SPACEBIOBENCH_EVALUATION_CARD.md)
- [Release status card](SPACEBIOBENCH_RELEASE_READINESS_CARD.md)
- [Canonical v7.1 results](CANONICAL_RESULTS_V7_1.md)
- [v7.1 Hugging Face dataset card](hf_dataset_card.md)
- [v9 metadata catalog card](v9_hf_dataset_card.md)

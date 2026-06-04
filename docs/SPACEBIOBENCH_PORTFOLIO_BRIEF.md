---
title: SpaceBio-Bench Portfolio Brief
page_type: portfolio_brief
status: draft_public_ready
last_reviewed: 2026-06-04
claim_boundary: portfolio_brief_no_new_release_claim
---

# SpaceBio-Bench Portfolio Brief

## One-Sentence Summary

SpaceBio-Bench is a mission-held-out transcriptomics benchmark and transparency
package for public NASA OSDR space-biology data, designed to evaluate model
generalization under mission shift while keeping provenance, evaluation scope,
release readiness, and claim boundaries explicit.

## What I Built

- A cross-mission benchmark surface for public space-biology transcriptomics,
  with canonical v1-v7 result documentation and a v7.1 consistency patch.
- A v9 public bulk metadata-alpha scaffold with 8 task manifests, 33 fold
  definitions, 22 public OSDR source rows, 24 scaffold baseline runs, and a
  draft Frictionless Data Package descriptor.
- A transparency card pack that separates system scope, evaluation
  interpretation, release readiness, and allowed/blocked claims.
- A Hugging Face-style dataset card that states the current metadata-alpha
  boundary without implying a frozen payload release.
- A claim register that prevents mixed-surface overclaiming across v1-v7
  results, the v8 translational extension, and v9 draft surfaces.

## Why It Matters

Space-biology transcriptomics is small-sample, mission-confounded, and easy to
overinterpret. A model can appear strong under a pooled score while failing on a
specific tissue, mission, or evaluation surface. SpaceBio-Bench treats that as
an evaluation-design problem, not just a modeling problem.

The project demonstrates how to build a benchmark that is technically useful
and release-disciplined: tasks are explicit, source rows are traceable, baseline
runs are documented, and unsupported claims are blocked before they become
public release language.

## Technical Depth

| Area | Evidence |
|---|---|
| Benchmark design | Mission-held-out LOMO task manifests and fold definitions |
| Bioinformatics scope | Public OSDR-derived mouse bulk RNA-seq task scaffold |
| Evaluation engineering | Baseline predictions, metrics, run manifests, and pooled-summary caveats |
| Provenance | Source inventory plus OSDR API and checksum-manifest evidence |
| Packaging | Draft Data Package descriptor and Hugging Face-style dataset card |
| Release discipline | Metadata-alpha boundary, payload blocker tracking, and readiness tiers |
| Scientific communication | System card, evaluation card, release readiness card, and claim register |

## Responsible Release Discipline

Current v9 public bulk status is intentionally limited:

- It is a metadata-only alpha, not a frozen payload release.
- Payload-level hash verification remains pending.
- Baselines are scaffold anchors, not a state-of-the-art leaderboard.
- Current evidence does not support clinical, crew-health, countermeasure,
  intervention, Mars-regime, or biological-mechanism claims.

This boundary is not a weakness of the project. It is part of the contribution:
the repository shows which claims are supported now, which claims are blocked,
and what evidence would be needed before stronger wording is appropriate.

## Role-Relevant Signal

This project is most directly relevant to roles involving:

- AI evaluation and benchmark design.
- Bioinformatics and space-biology data infrastructure.
- Scientific ML / foundation-model evaluation in biology.
- AI safety evaluation, especially claim-boundary and release-readiness work.
- Research engineering for reproducible, auditable scientific artifacts.
- Technical communication for responsible biological AI evaluation systems.

## Fast Review Path

For a quick review, read these in order:

1. `docs/SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md`
2. `docs/SPACEBIOBENCH_SYSTEM_CARD.md`
3. `docs/SPACEBIOBENCH_EVALUATION_CARD.md`
4. `docs/SPACEBIOBENCH_RELEASE_READINESS_CARD.md`
5. `docs/SPACEBIOBENCH_CLAIM_REGISTER.md`
6. `docs/v9_hf_dataset_card.md`
7. `docs/CANONICAL_RESULTS_V7_1.md`

## Resume-Ready Summary

Built SpaceBio-Bench, a mission-held-out space-biology transcriptomics benchmark
and transparency package using public NASA OSDR data. The project defines
task/fold manifests, scaffold baselines, provenance and checksum evidence,
Hugging Face-style dataset documentation, release-readiness tiers, and a claim
register to prevent unsupported model, clinical, or payload-release claims.

## Cover-Letter Variant

SpaceBio-Bench reflects the kind of evaluation work I want to keep doing:
building biological AI benchmarks that are technically useful, scientifically
careful, and release-disciplined. Beyond model scores, the project documents
task boundaries, source provenance, baseline interpretation, payload-readiness
blockers, and allowed versus unsupported claims. That combination is especially
important for biological AI evaluation, where overclaiming can happen quickly
when small datasets, mission confounding, and foundation-model comparisons are
compressed into a single headline number.

---
title: SpaceBio-Bench System Card
page_type: system_card
status: draft
last_reviewed: 2026-06-04
claim_boundary: benchmark_system_card_draft_no_new_release_claim
---

# SpaceBio-Bench System Card

## Card Purpose

This card documents SpaceBio-Bench as a benchmark system, not as a trained
model. It summarizes the project's data surfaces, task definitions, evaluation
logic, provenance controls, release boundaries, and known limitations.

This card does not replace the Hugging Face dataset cards, task manifests,
release manifests, or model-specific documentation. Individual trained models
or adapters should receive separate model cards when they become release
artifacts.

## System Summary

SpaceBio-Bench is a mission-held-out transcriptomics benchmark scaffold for
public space-biology data. It organizes public OSDR-derived expression tasks
into task manifests, fold definitions, baseline outputs, provenance metadata,
and release-boundary documents.

The project currently has multiple surfaces with different maturity levels:

| Surface | Current status | Primary boundary |
|---|---|---|
| v1-v7 GeneLab Benchmark | Canonical historical result surface | Cross-mission mouse transcriptomics benchmark and result synthesis |
| v7.1 documentation patch | Documentation consistency patch | No new benchmark result generation |
| v8 SpaceMed | Incubating translational extension | Hypothesis and diagnostic work only; do not mix into v7.1 claims |
| v9 public bulk | Metadata-only alpha snapshot | Task/source/provenance scaffold; not a frozen payload release |
| v9 single-cell, organoid, multispecies | Draft diagnostic lanes | Asset, payload, metric, or feasibility scaffolds depending on lane |

## Intended Uses

SpaceBio-Bench is intended for:

- Evaluating method generalization under mission-held-out spaceflight or analog
  transcriptomics shift.
- Comparing classical ML, gene-expression foundation-model adapters, and other
  biological AI methods when inputs, tasks, and result surfaces are explicitly
  matched.
- Stress-testing whether benchmark claims remain grounded in source accessions,
  task manifests, fold definitions, and run manifests.
- Developing provenance-aware benchmark packaging for public OSDR-derived
  biological data.
- Supporting reviewer-facing analysis where per-task results, claim boundaries,
  and upstream data citations are visible.

## What This Enables

The card pack is meant to make the benchmark easier to inspect before a reader
trusts a headline number. It exposes what is being evaluated, which result
surface a claim belongs to, which files support that claim, and which release
conditions remain unfinished.

## Out Of Scope

SpaceBio-Bench should not be used to claim:

- Astronaut health-risk prediction.
- Clinical or crew-health recommendations.
- Operational readiness for space missions.
- Countermeasure, intervention, or dosing recommendations.
- Mars-regime point predictions.
- A frozen v9 payload release before local payload mirroring and payload-level
  hash verification are complete.
- A uniform foundation-model leaderboard when compared rows come from different
  task subsets, tissues, adapters, or evaluation surfaces.

## System Components

| Component | Role | Evidence |
|---|---|---|
| Task manifests | Define benchmark tasks, folds, labels, sources, and metric ids | `v9/task_manifests/*.json`; `v9/task_manifest_index.csv` |
| Fold data index | Tracks train/test files, held-out missions, row counts, and gene counts | `v9/task_data_index.csv` |
| Source inventory | Maps OSDR accessions, GLDS prefixes, tissues, missions, and access status | `v9/source_inventory.csv` |
| Source checksum audit | Records OSDR API and checksum-manifest evidence | `v9/source_checksum_audit.csv`; `v9/source_checksum_audit.json` |
| Draft Data Package | Machine-readable descriptor for metadata and output resources | `v9/datapackage.draft.json` |
| Baseline outputs | Validate task/evaluation workflow, not model superiority | `v9/reports/bulk_lomo_baseline_summary.csv` and per-baseline reports |
| Dataset cards | Human-facing dataset and release-boundary summaries | `docs/hf_dataset_card.md`; `docs/v9_hf_dataset_card.md` |
| Canonical result doc | Public-facing v7.1 result and scope source of truth | `docs/CANONICAL_RESULTS_V7_1.md` |
| Evaluation card | Task, metric, baseline, and result-interpretation controls | `docs/SPACEBIOBENCH_EVALUATION_CARD.md` |
| Release readiness card | Release tiers, readiness gates, and blocked release claims | `docs/SPACEBIOBENCH_RELEASE_READINESS_CARD.md` |
| Claim register | Claim, support, confidence, and blocked wording control | `docs/SPACEBIOBENCH_CLAIM_REGISTER.md` |

## Evaluation Philosophy

SpaceBio-Bench treats mission shift as the main evaluation pressure. The core
public bulk task shape is leave-one-mission-out classification within a tissue
context: train on some missions, hold out one mission, and evaluate whether the
method generalizes to the held-out mission.

Evaluation should be read at task and fold level before pooled summaries.
Pooled averages are useful for navigation, but they can hide tissue-specific,
mission-specific, or label-source-specific failure modes.

Current v9 public bulk baselines are scaffold baselines. They validate the
benchmark workflow and provide comparison anchors, but they are not tuned
leaderboard endpoints and should not be used for biological superiority claims.

## Data And Provenance Boundary

The v9 public bulk lane is currently a metadata-only alpha snapshot. It includes
task manifests, source inventory, OSDR checksum-manifest evidence,
alpha-boundary decision tables, baseline outputs, and a draft Frictionless Data
Package descriptor.

The current boundary is:

- OSDR API and checksum-manifest evidence exists for all 22 public bulk source
  rows.
- The package does not include a frozen local payload mirror.
- Payload-level hash verification for every distributed fold matrix remains
  pending.
- The descriptor is `v9/datapackage.draft.json`, not a release
  `datapackage.json`.

The v9 dataset card records:

- `spacebio_bench:release_status = metadata_alpha_not_frozen`
- `spacebio_bench:alpha_snapshot_status = metadata_only_alpha_snapshot`
- `spacebio_bench:claim_boundary = metadata_only_public_bulk_alpha_no_payload_release`
- `spacebio_bench:payload_release_allowed = false`
- `spacebio_bench:payload_verification_status = checksum_manifests_parsed_payloads_not_hashed`

## Results Boundary

The v7.1 canonical result surface is the current source of truth for public
v1-v7 scope accounting, headline result tables, foundation-model comparison
language, held-out validation, and submission-safe claim wording.

The v9 public bulk baseline results are a separate alpha scaffold surface. They
should be described as baseline workflow evidence, not as final model ranking or
biological mechanism evidence.

## Known Limitations

- Mouse bulk RNA-seq is not a complete representation of space biology.
- Mission labels can conflate spaceflight exposure, vehicle, hardware, age,
  protocol, tissue handling, processing, and batch effects.
- Some task labels include analog or special mission labels that should not be
  treated as interchangeable with ISS or deep-space exposure.
- Bulk RNA-seq does not resolve cell-type-specific effects.
- Legacy processed fold matrices may reflect earlier preprocessing choices.
- Current v9 public bulk packaging is not payload-frozen.
- Foundation-model comparisons need adapter-specific validation and matched
  evaluation surfaces.
- Strong benchmark scores do not prove biological mechanism, translational
  validity, or operational deployment readiness.

## Responsible Use

Users should:

- Cite OSDR and the individual upstream OSDR datasets used in a downstream
  analysis.
- Report per-task and per-fold metrics alongside any pooled summary.
- Keep v7.1 result claims separate from v8 translational and v9 metadata-alpha
  claims.
- State whether a result is from a canonical result surface, a metadata alpha,
  a diagnostic lane, or a draft feasibility lane.
- Treat intervention, countermeasure, clinical, and crew-health claims as out
  of scope unless a future release explicitly validates them.

## Release And Integrity Controls

Current controls:

- Explicit claim-boundary tables in v9 public bulk alpha reports.
- Source inventory and source checksum audit artifacts.
- Draft Data Package metadata with release-status and payload-boundary fields.
- Canonical v7.1 result document to prevent mixed-surface claims.
- Human-facing dataset cards for v7.1 and v9 public bulk alpha.

Readiness controls:

- Maintain the release readiness card with release tiers such as
  `metadata_alpha`, `diagnostic_alpha`, `frozen_payload_release`, and
  `doi_release`.
- Maintain the evaluation card that separates task validity, leakage controls,
  baseline status, and result interpretation.
- Add a payload-level SHA-256 manifest before any frozen v9 payload language.
- Add RO-Crate metadata for research-object citation and provenance.
- Add BagIt-style payload manifest checks if a distributable payload bundle is
  created.

## Companion Documents

Local evidence:

- `docs/CANONICAL_RESULTS_V7_1.md`
- `docs/hf_dataset_card.md`
- `docs/v9_hf_dataset_card.md`
- `docs/V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md`
- `docs/V9_PUBLIC_BULK_ALPHA_CARD_DATAPACKAGE_BOUNDARY_UPDATE.md`
- `docs/SPACEBIOBENCH_EVALUATION_CARD.md`
- `docs/SPACEBIOBENCH_RELEASE_READINESS_CARD.md`
- `docs/SPACEBIOBENCH_CLAIM_REGISTER.md`
- `v9/datapackage.draft.json`
- `v9/source_inventory.csv`
- `v9/source_checksum_audit.csv`
- `v9/task_manifest_index.csv`
- `v9/task_data_index.csv`

External best-practice anchors:

- Model Cards for Model Reporting: https://arxiv.org/abs/1810.03993
- Hugging Face model cards: https://huggingface.co/docs/hub/main/model-cards
- Hugging Face dataset cards: https://huggingface.co/docs/datasets/v2.7.0/en/dataset_card
- Datasheets for Datasets: https://www.microsoft.com/en-us/research/uploads/prod/2019/01/1803.09010.pdf
- NIST AI RMF: https://www.nist.gov/itl/ai-risk-management-framework
- NASA OSDR FAQ: https://science.nasa.gov/reference/osdr-faq/
- NASA OSDR Biological Data API: https://visualization.osdr.nasa.gov/biodata/api/
- FAIR principles: https://www.nature.com/articles/sdata201618
- RO-Crate: https://www.researchobject.org/ro-crate/technical_overview
- BagIt RFC 8493: https://www.rfc-editor.org/info/rfc8493/
- DataCite Metadata Schema: https://schema.datacite.org/

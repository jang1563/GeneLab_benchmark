---
title: SpaceBio-Bench Claim Register
page_type: evidence_register
status: draft
last_reviewed: 2026-06-04
claim_boundary: benchmark_claim_register_draft_no_new_release_claim
---

# SpaceBio-Bench Claim Register

## Purpose

This register ties major SpaceBio-Bench claims to support level, source files,
confidence, allowed language, and blocked language. It is designed to prevent
mixed-surface claims across v1-v7, v8, and v9.

The register is a documentation control, not a new release artifact and not a
source of new benchmark results.

## Evidence Vocabulary

- `primary`: local source file, generated artifact, official external page, or
  primary paper.
- `project-synthesis`: synthesis from multiple local primary sources.
- `inference`: reasonable interpretation from primary sources; should be
  phrased cautiously.
- `blocked`: not currently supported by the evidence boundary.
- `future`: plausible future claim after named blockers are resolved.

## Reader Summary

The current public-safe summary is:

- v1-v7 is the canonical historical benchmark result surface.
- v7.1 is a documentation consistency patch, not a new result release.
- v8 is an incubating translational extension and should not be mixed into
  v7.1 claims.
- v9 public bulk is a metadata-only alpha with explicit payload blockers.
- Current v9 baselines are scaffold anchors, not final leaderboard entries.
- Clinical, crew-health, intervention, countermeasure, and Mars-regime claims
  are out of scope for the current evidence boundary.

## Current Supported Claims

| ID | Claim | Support | Primary source | Confidence | Allowed wording | Blocked wording |
|---|---|---|---|---|---|---|
| SBB-C001 | v1-v7 GeneLab Benchmark evaluates cross-mission generalization of mouse spaceflight transcriptomic signatures. | primary | `docs/CANONICAL_RESULTS_V7_1.md`; `docs/hf_dataset_card.md` | high | "cross-mission mouse spaceflight transcriptomics benchmark" | "clinical astronaut health predictor" |
| SBB-C002 | v7.1 is a documentation consistency patch, not a new result-generation release. | primary | `docs/CANONICAL_RESULTS_V7_1.md` | high | "v7.1 documentation patch; no new benchmark result generation" | "v7.1 adds new benchmark results" |
| SBB-C003 | The v4 multi-method surface covers 8 tissues, 8 classifiers, and 4 feature types for 256 evaluations. | primary | `docs/CANONICAL_RESULTS_V7_1.md` | high | "v4 multi-method evaluation: 8 tissues x 8 classifiers x 4 feature types" | "all methods were evaluated on every later v8/v9 task" |
| SBB-C004 | Foundation-model and text-LLM rows are mixed-surface snapshots and should not be described as a single uniform 8-tissue FM leaderboard. | primary | `docs/CANONICAL_RESULTS_V7_1.md` | high | "benchmark-surface summary with subset notes" | "single uniform 8-tissue foundation-model leaderboard" |
| SBB-C005 | Current gene-expression FMs in the canonical v7.1 snapshot do not automatically outperform tuned classical baselines under small-n bulk RNA-seq shift. | primary | `docs/CANONICAL_RESULTS_V7_1.md` | high | "do not automatically outperform tuned classical baselines" | "foundation models fail at space biology" |
| SBB-C006 | v8 SpaceMed is an incubating translational extension and should not be mixed into v7.1 benchmark claims. | primary | `docs/CANONICAL_RESULTS_V7_1.md`; `docs/V8_BETA_RELEASE_PLAN_2026_05_10.md` | high | "incubating translational extension" | "v8 proves countermeasure efficacy" |
| SBB-C007 | v9 public bulk is a metadata-only alpha snapshot, not a frozen payload release. | primary | `docs/v9_hf_dataset_card.md`; `docs/V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md`; `docs/V9_PUBLIC_BULK_ALPHA_CARD_DATAPACKAGE_BOUNDARY_UPDATE.md` | high | "SpaceBio-Bench v9 public bulk metadata alpha" | "frozen v9 public benchmark release" |
| SBB-C008 | v9 public bulk currently includes 8 generated public bulk LOMO task manifests, 6 tissue contexts, 22 deduplicated public OSDR source rows, 33 fold definitions, 24 baseline runs, and 21 draft Data Package resources. | primary | `docs/v9_hf_dataset_card.md`; `v9/datapackage.draft.json` | high | Use counts with "current public bulk draft" qualifier | Use counts as a frozen DOI release inventory |
| SBB-C009 | OSDR API and checksum-manifest evidence has been parsed for all 22 public bulk source rows, but local payload-level hash verification is still pending. | primary | `docs/v9_hf_dataset_card.md`; `v9/source_checksum_audit.csv`; `docs/V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md` | high | "checksum-manifest evidence parsed; payload hashing pending" | "locally hash-verified payload bundle" |
| SBB-C010 | v9 public bulk baselines validate the scaffold workflow and provide anchors, but are not tuned leaderboard endpoints. | primary | `docs/v9_hf_dataset_card.md`; `v9/reports/bulk_lomo_baseline_summary.csv` | high | "scaffold baselines" | "final model rankings" |
| SBB-C011 | Per-task and per-fold reporting should accompany pooled summaries because pooled averages can hide mission or tissue failures. | primary/project-synthesis | `docs/v9_hf_dataset_card.md`; `docs/CANONICAL_RESULTS_V7_1.md` | high | "report per-task results, not only pooled averages" | "single pooled score fully characterizes the method" |
| SBB-C012 | Mission labels can conflate biological spaceflight signal with vehicle, hardware, protocol, tissue handling, time, and processing effects. | primary/project-synthesis | `docs/v9_hf_dataset_card.md`; `docs/SPACEBIOBENCH_SYSTEM_CARD.md` | medium-high | "mission-shift benchmark with known confounding risks" | "pure microgravity effect estimator" |
| SBB-C013 | Public bulk tasks do not support clinical, crew-health, countermeasure, intervention, or Mars-regime claims. | primary | `docs/v9_hf_dataset_card.md`; `docs/CANONICAL_RESULTS_V7_1.md` | high | "benchmark evidence, not biological mechanism or operational recommendation" | "astronaut health-risk or countermeasure recommendation" |
| SBB-C014 | OSDR and individual OSDR datasets should be credited and cited for downstream analyses. | primary | `docs/v9_hf_dataset_card.md`; NASA OSDR FAQ | high | "Data are courtesy of the NASA Open Science Data Repository" plus dataset-specific citations | Hand-written substitute citations without checking OSDR study pages |
| SBB-C015 | Dataset cards should document contents, context, intended use, creation, responsible use, and limitations. | primary | Hugging Face dataset-card docs; Datasheets for Datasets | high | "dataset card as human-facing responsible-use surface" | "README with only download commands is sufficient" |
| SBB-C016 | Model cards are appropriate for individual trained models or adapters, but SpaceBio-Bench itself should be documented as a benchmark/system card. | primary/inference | Model Cards paper; Hugging Face model-card docs; `docs/SPACEBIOBENCH_SYSTEM_CARD.md` | high | "benchmark/system card for the project; model cards for individual baselines" | "the full benchmark is a model card" |

## Release Readiness And Future Claims

| ID | Claim | Support | Primary source | Confidence | Allowed wording | Blocked wording |
|---|---|---|---|---|---|---|
| SBB-C017 | A future frozen payload release should add machine-readable payload manifests and verification reports. | primary/inference | BagIt RFC 8493; `docs/v9_hf_dataset_card.md`; `docs/V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md` | high | "payload-level SHA-256 manifest before frozen payload language" | "payload freeze without payload-level hashes" |
| SBB-C018 | A future citable research-object release should add RO-Crate or equivalent provenance metadata. | primary/inference | RO-Crate technical overview; `docs/V9_LONG_RUN_OPERATING_PROTOCOL.md`; `docs/v9_hf_dataset_card.md` | medium-high | "future RO-Crate export for research-object provenance" | "current v9 alpha is already a complete citable RO-Crate release" |
| SBB-C019 | DataCite-style metadata is useful for future DOI-oriented release planning. | primary/inference | DataCite Metadata Schema; `docs/V9_LONG_RUN_OPERATING_PROTOCOL.md` | medium-high | "align release metadata with DataCite fields before DOI/archive release" | "DOI release ready without creator, version, related identifier, license, and resource type review" |
| SBB-C020 | NIST AI RMF can inform documentation structure, but SpaceBio-Bench is not a deployed AI product. | primary/inference | NIST AI RMF; `docs/SPACEBIOBENCH_SYSTEM_CARD.md` | medium | "use Govern/Map/Measure/Manage as documentation lenses" | "NIST compliance claim" |
| SBB-C021 | Evaluation should be interpreted through task, fold, source, payload-boundary, and run-manifest evidence before pooled summaries. | primary/project-synthesis | `docs/SPACEBIOBENCH_EVALUATION_CARD.md`; `docs/v9_hf_dataset_card.md`; `v9/task_data_index.csv` | high | "read per-task and per-fold metrics before pooled means" | "single pooled mean fully establishes benchmark performance" |
| SBB-C022 | v9 public bulk currently satisfies a metadata-alpha tier, not a frozen-payload or DOI/archive tier. | primary | `docs/SPACEBIOBENCH_RELEASE_READINESS_CARD.md`; `docs/V9_PUBLIC_BULK_ALPHA_CARD_DATAPACKAGE_BOUNDARY_UPDATE.md`; `v9/reports/public_bulk_alpha_snapshot_decision/snapshot_decision_summary.csv` | high | "metadata alpha with explicit payload blockers" | "frozen-payload or archive-ready release" |

## Blocked Or Future Claims

| Future claim | Current blocker | Evidence needed |
|---|---|---|
| Frozen v9 public bulk payload release | Payload mirroring and payload-level hash verification pending | Local payload mirror, SHA-256 manifest, verification report, release `datapackage.json` |
| DOI or archive-ready release | Metadata alpha status and license/citation review pending | DataCite-aligned metadata, final license, dataset-specific OSDR citations, archive manifest |
| Complete RO-Crate research object | RO-Crate export not yet created | `ro-crate-metadata.json`, entity graph, workflow and provenance links |
| Foundation-model leaderboard | Adapter validation and matched surfaces incomplete | Matched task inputs, adapter cards, run manifests, per-task metrics, leakage checks |
| Biological mechanism proof | Benchmark scores are not mechanistic validation | Independent biological validation, mechanistic assays, matched causal analysis |
| Countermeasure or intervention recommendation | v8/v9 diagnostic claims do not validate interventions | Controlled intervention evidence, safety review, translational validation |
| Crew-health or clinical decision support | Public benchmark scope excludes clinical recommendations | Controlled human-data review, clinical validation, institutional review |

## Maintenance Rules

- Update this register whenever public-facing result counts, release status,
  or allowed language changes.
- Do not promote a claim from `inference` to `primary` without adding the local
  or external primary source.
- Keep v7.1 result claims, v8 translational hypotheses, and v9 alpha scaffold
  claims in separate rows.
- When in doubt, prefer allowed wording that names the release surface and
  status explicitly.

## External Sources

- Model Cards for Model Reporting: https://arxiv.org/abs/1810.03993
- Hugging Face model cards: https://huggingface.co/docs/hub/main/model-cards
- Hugging Face dataset cards: https://huggingface.co/docs/datasets/v2.7.0/en/dataset_card
- Datasheets for Datasets: https://www.microsoft.com/en-us/research/uploads/prod/2019/01/1803.09010.pdf
- NIST AI RMF: https://www.nist.gov/itl/ai-risk-management-framework
- NASA OSDR FAQ: https://science.nasa.gov/reference/osdr-faq/
- RO-Crate: https://www.researchobject.org/ro-crate/technical_overview
- BagIt RFC 8493: https://www.rfc-editor.org/info/rfc8493/
- DataCite Metadata Schema: https://schema.datacite.org/

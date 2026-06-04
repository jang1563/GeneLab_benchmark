---
title: SpaceBio-Bench Claim Register
page_type: evidence_register
status: public_review_ready
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

Branch note: on the default `main` branch, v9-specific evidence paths such as
`v9/...` and `docs/V9_*` refer to artifacts maintained on the canonical `v3`
branch.

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

## Current Supported Claim Cards

### SBB-C001 - Cross-Mission Benchmark Surface

- Claim: v1-v7 GeneLab Benchmark evaluates cross-mission generalization of
  mouse spaceflight transcriptomic signatures.
- Support: primary.
- Sources: `docs/CANONICAL_RESULTS_V7_1.md`; `docs/hf_dataset_card.md`.
- Confidence: high.
- Use: "cross-mission mouse spaceflight transcriptomics benchmark".
- Avoid: "clinical astronaut health predictor".

### SBB-C002 - v7.1 Patch Boundary

- Claim: v7.1 is a documentation consistency patch, not a new result-generation
  release.
- Support: primary.
- Sources: `docs/CANONICAL_RESULTS_V7_1.md`.
- Confidence: high.
- Use: "v7.1 documentation patch; no new benchmark result generation".
- Avoid: "v7.1 adds new benchmark results".

### SBB-C003 - v4 Multi-Method Scope

- Claim: the v4 multi-method surface covers 8 tissues, 8 classifiers, and 4
  feature types for 256 evaluations.
- Support: primary.
- Sources: `docs/CANONICAL_RESULTS_V7_1.md`.
- Confidence: high.
- Use: "v4 multi-method evaluation: 8 tissues x 8 classifiers x 4 feature
  types".
- Avoid: "all methods were evaluated on every later v8/v9 task".

### SBB-C004 - Foundation-Model Snapshot Boundary

- Claim: foundation-model and text-LLM rows are mixed-surface snapshots, not a
  single uniform 8-tissue FM leaderboard.
- Support: primary.
- Sources: `docs/CANONICAL_RESULTS_V7_1.md`.
- Confidence: high.
- Use: "benchmark-surface summary with subset notes".
- Avoid: "single uniform 8-tissue foundation-model leaderboard".

### SBB-C005 - Classical Baseline Comparison

- Claim: current gene-expression FMs in the canonical v7.1 snapshot do not
  automatically outperform tuned classical baselines under small-n bulk RNA-seq
  shift.
- Support: primary.
- Sources: `docs/CANONICAL_RESULTS_V7_1.md`.
- Confidence: high.
- Use: "do not automatically outperform tuned classical baselines".
- Avoid: "foundation models fail at space biology".

### SBB-C006 - v8 Boundary

- Claim: v8 SpaceMed is an incubating translational extension and should not be
  mixed into v7.1 benchmark claims.
- Support: primary.
- Sources: `docs/CANONICAL_RESULTS_V7_1.md`;
  `docs/V8_BETA_RELEASE_PLAN_2026_05_10.md`.
- Confidence: high.
- Use: "incubating translational extension".
- Avoid: "v8 proves countermeasure efficacy".

### SBB-C007 - v9 Metadata-Alpha Boundary

- Claim: v9 public bulk is a metadata-only alpha snapshot, not a frozen payload
  release.
- Support: primary.
- Sources: `docs/v9_hf_dataset_card.md`;
  `docs/V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md`;
  `docs/V9_PUBLIC_BULK_ALPHA_CARD_DATAPACKAGE_BOUNDARY_UPDATE.md`.
- Confidence: high.
- Use: "SpaceBio-Bench v9 public bulk metadata alpha".
- Avoid: "frozen v9 public benchmark release".

### SBB-C008 - v9 Public Bulk Inventory

- Claim: v9 public bulk currently includes 8 generated public bulk LOMO task
  manifests, 6 tissue contexts, 22 deduplicated public OSDR source rows, 33
  fold definitions, 24 baseline runs, and 21 draft Data Package resources.
- Support: primary.
- Sources: `docs/v9_hf_dataset_card.md`; `v9/datapackage.draft.json`.
- Confidence: high.
- Use: counts with a "current public bulk draft" qualifier.
- Avoid: using counts as a frozen DOI release inventory.

### SBB-C009 - Checksum Evidence Boundary

- Claim: OSDR API and checksum-manifest evidence has been parsed for all 22
  public bulk source rows, but local payload-level hash verification is still
  pending.
- Support: primary.
- Sources: `docs/v9_hf_dataset_card.md`; `v9/source_checksum_audit.csv`;
  `docs/V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md`.
- Confidence: high.
- Use: "checksum-manifest evidence parsed; payload hashing pending".
- Avoid: "locally hash-verified payload bundle".

### SBB-C010 - Baseline Status

- Claim: v9 public bulk baselines validate the scaffold workflow and provide
  anchors, but are not tuned leaderboard endpoints.
- Support: primary.
- Sources: `docs/v9_hf_dataset_card.md`;
  `v9/reports/bulk_lomo_baseline_summary.csv`.
- Confidence: high.
- Use: "scaffold baselines".
- Avoid: "final model rankings".

### SBB-C011 - Per-Task Reporting

- Claim: per-task and per-fold reporting should accompany pooled summaries
  because pooled averages can hide mission or tissue failures.
- Support: primary/project-synthesis.
- Sources: `docs/v9_hf_dataset_card.md`; `docs/CANONICAL_RESULTS_V7_1.md`.
- Confidence: high.
- Use: "report per-task results, not only pooled averages".
- Avoid: "single pooled score fully characterizes the method".

### SBB-C012 - Mission Confounding

- Claim: mission labels can conflate biological spaceflight signal with
  vehicle, hardware, protocol, tissue handling, time, and processing effects.
- Support: primary/project-synthesis.
- Sources: `docs/v9_hf_dataset_card.md`; `docs/SPACEBIOBENCH_SYSTEM_CARD.md`.
- Confidence: medium-high.
- Use: "mission-shift benchmark with known confounding risks".
- Avoid: "pure microgravity effect estimator".

### SBB-C013 - Unsupported Operational Claims

- Claim: public bulk tasks do not support clinical, crew-health,
  countermeasure, intervention, or Mars-regime claims.
- Support: primary.
- Sources: `docs/v9_hf_dataset_card.md`; `docs/CANONICAL_RESULTS_V7_1.md`.
- Confidence: high.
- Use: "benchmark evidence, not biological mechanism or operational
  recommendation".
- Avoid: "astronaut health-risk or countermeasure recommendation".

### SBB-C014 - OSDR Credit

- Claim: OSDR and individual OSDR datasets should be credited and cited for
  downstream analyses.
- Support: primary.
- Sources: `docs/v9_hf_dataset_card.md`; NASA OSDR FAQ.
- Confidence: high.
- Use: "Data are courtesy of the NASA Open Science Data Repository" plus
  dataset-specific citations.
- Avoid: hand-written substitute citations without checking OSDR study pages.

### SBB-C015 - Dataset-Card Role

- Claim: dataset cards should document contents, context, intended use,
  creation, responsible use, and limitations.
- Support: primary.
- Sources: Hugging Face dataset-card docs; Datasheets for Datasets.
- Confidence: high.
- Use: "dataset card as human-facing responsible-use surface".
- Avoid: "README with only download commands is sufficient".

### SBB-C016 - System Card Versus Model Card

- Claim: model cards are appropriate for individual trained models or adapters,
  but SpaceBio-Bench itself should be documented as a benchmark/system card.
- Support: primary/inference.
- Sources: Model Cards paper; Hugging Face model-card docs;
  `docs/SPACEBIOBENCH_SYSTEM_CARD.md`.
- Confidence: high.
- Use: "benchmark/system card for the project; model cards for individual
  baselines".
- Avoid: "the full benchmark is a model card".

## Release Readiness And Future Claim Cards

### SBB-C017 - Frozen Payload Requirements

- Claim: a future frozen payload release should add machine-readable payload
  manifests and verification reports.
- Support: primary/inference.
- Sources: BagIt RFC 8493; `docs/v9_hf_dataset_card.md`;
  `docs/V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md`.
- Confidence: high.
- Use: "payload-level SHA-256 manifest before frozen payload language".
- Avoid: "payload freeze without payload-level hashes".

### SBB-C018 - Research-Object Provenance

- Claim: a future citable research-object release should add RO-Crate or
  equivalent provenance metadata.
- Support: primary/inference.
- Sources: RO-Crate technical overview; `docs/V9_LONG_RUN_OPERATING_PROTOCOL.md`;
  `docs/v9_hf_dataset_card.md`.
- Confidence: medium-high.
- Use: "future RO-Crate export for research-object provenance".
- Avoid: "current v9 alpha is already a complete citable RO-Crate release".

### SBB-C019 - DOI-Oriented Metadata

- Claim: DataCite-style metadata is useful for future DOI-oriented release
  planning.
- Support: primary/inference.
- Sources: DataCite Metadata Schema; `docs/V9_LONG_RUN_OPERATING_PROTOCOL.md`.
- Confidence: medium-high.
- Use: "align release metadata with DataCite fields before DOI/archive release".
- Avoid: "DOI release ready without creator, version, related identifier,
  license, and resource type review".

### SBB-C020 - NIST AI RMF Lens

- Claim: NIST AI RMF can inform documentation structure, but SpaceBio-Bench is
  not a deployed AI product.
- Support: primary/inference.
- Sources: NIST AI RMF; `docs/SPACEBIOBENCH_SYSTEM_CARD.md`.
- Confidence: medium.
- Use: "use Govern/Map/Measure/Manage as documentation lenses".
- Avoid: "NIST compliance claim".

### SBB-C021 - Evaluation Reading Order

- Claim: evaluation should be interpreted through task, fold, source,
  payload-boundary, and run-manifest evidence before pooled summaries.
- Support: primary/project-synthesis.
- Sources: `docs/SPACEBIOBENCH_EVALUATION_CARD.md`;
  `docs/v9_hf_dataset_card.md`; `v9/task_data_index.csv`.
- Confidence: high.
- Use: "read per-task and per-fold metrics before pooled means".
- Avoid: "single pooled mean fully establishes benchmark performance".

### SBB-C022 - v9 Release Tier

- Claim: v9 public bulk currently satisfies a metadata-alpha tier, not a
  frozen-payload or DOI/archive tier.
- Support: primary.
- Sources: `docs/SPACEBIOBENCH_RELEASE_READINESS_CARD.md`;
  `docs/V9_PUBLIC_BULK_ALPHA_CARD_DATAPACKAGE_BOUNDARY_UPDATE.md`;
  `v9/reports/public_bulk_alpha_snapshot_decision/snapshot_decision_summary.csv`.
- Confidence: high.
- Use: "metadata alpha with explicit payload blockers".
- Avoid: "frozen-payload or archive-ready release".

## Blocked Or Future Claims

### Frozen v9 Public Bulk Payload Release

- Current blocker: payload mirroring and payload-level hash verification are
  pending.
- Evidence needed: local payload mirror, SHA-256 manifest, verification report,
  and release `datapackage.json`.

### DOI Or Archive-Ready Release

- Current blocker: metadata alpha status and license/citation review are
  pending.
- Evidence needed: DataCite-aligned metadata, final license, dataset-specific
  OSDR citations, and archive manifest.

### Complete RO-Crate Research Object

- Current blocker: RO-Crate export has not yet been created.
- Evidence needed: `ro-crate-metadata.json`, entity graph, workflow links, and
  provenance links.

### Foundation-Model Leaderboard

- Current blocker: adapter validation and matched evaluation surfaces are
  incomplete.
- Evidence needed: matched task inputs, adapter cards, run manifests, per-task
  metrics, and leakage checks.

### Biological Mechanism Proof

- Current blocker: benchmark scores are not mechanistic validation.
- Evidence needed: independent biological validation, mechanistic assays, and
  matched causal analysis.

### Countermeasure Or Intervention Recommendation

- Current blocker: v8/v9 diagnostic claims do not validate interventions.
- Evidence needed: controlled intervention evidence, safety review, and
  translational validation.

### Crew-Health Or Clinical Decision Support

- Current blocker: public benchmark scope excludes clinical recommendations.
- Evidence needed: controlled human-data review, clinical validation, and
  institutional review.

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

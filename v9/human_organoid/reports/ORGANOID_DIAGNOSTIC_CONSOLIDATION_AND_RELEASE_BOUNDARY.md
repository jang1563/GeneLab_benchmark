# Human Organoid Diagnostic Consolidation And Release Boundary

Status: consolidation complete
Date: 2026-05-23
Task: `draft_human_organoid_spaceflight`
Block: `V9-ORG-031`

## Scope

This note consolidates the completed human organoid diagnostic family and marks
which artifacts belong to the draft v9 alpha surface versus exploratory or
negative-comparison work.

It closes the current organoid diagnostic loop:

- response-signature scorer and adapters;
- derived DE reference;
- source-transfer response signatures;
- classifier-derived feature-effect diagnostics;
- calibrated top-k feature-effect nulls;
- PCA-LR reconstructed feature-effect comparison.

## Release Boundary

The draft v9 alpha organoid surface should include three layers.

### Layer 1: Task And Provenance Surface

These are part of the draft-alpha organoid task definition:

| Artifact | Role |
|---|---|
| `v9/human_organoid/source_inventory.draft.csv` | source inventory |
| `v9/human_organoid/source_checksum_audit.draft.csv` | source-level OSDR evidence |
| `v9/human_organoid/sample_table_audit.draft.csv` | OSDR sample-table evidence |
| `v9/human_organoid/sample_factors.draft.csv` | parsed sample factors |
| `v9/human_organoid/geo_sample_metadata.draft.csv` | donor/Subject and iPSC-line recovery |
| `v9/human_organoid/expression_matrix_audit.draft.csv` | normalized matrix evidence |
| `v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json` | draft task manifest |
| `v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz` | derived DE reference |
| `v9/human_organoid/de_references/human_organoid_de_reference_manifest.draft.json` | DE reference manifest |

Boundary:

- still draft, not frozen;
- payload-level checksum freeze remains pending;
- classification metrics remain primary for the draft task;
- DE/signature and feature-effect metrics remain diagnostic.

### Layer 2: Default Diagnostic Surface

These are the recommended default diagnostic reports to show with the organoid
task:

| Artifact | Decision |
|---|---|
| `v9/human_organoid/reports/source_transfer_signature/` | default response-signature diagnostic |
| `v9/human_organoid/reports/ORGANOID_SOURCE_TRANSFER_SIGNATURE_REVIEW.md` | interpretation boundary |
| `v9/human_organoid/reports/logistic_feature_effect/` | default feature-effect diagnostic |
| `v9/human_organoid/reports/ORGANOID_FEATURE_EFFECT_NULL_CALIBRATION_REVIEW.md` | calibrated top-k interpretation |

Default diagnostic interpretation:

- global source-transfer response signatures are the first conservative
  response-signature baseline;
- L2 logistic gene-space feature effects are the first conservative
  classifier-effect baseline;
- both are non-primary and post hoc reference-scored;
- neither should be described as a leaderboard metric.

### Layer 3: Retained Secondary Or Exploratory Reports

These reports should stay in the repository but not be presented as default
diagnostic baselines:

| Artifact | Classification | Reason |
|---|---|---|
| `v9/human_organoid/reports/microglia_source_transfer_signature/` | retained secondary diagnostic | direction improves, rank worsens |
| `v9/human_organoid/reports/shared_control_source_transfer_signature/` | exploratory partial/negative diagnostic | partial coverage and weaker than broader signatures |
| `v9/human_organoid/reports/pca_lr_feature_effect/` | exploratory negative comparison | no signal gain over L2 logistic |
| `v9/human_organoid/reports/response_signature_smoke/` | scorer plumbing only | mirrored fixture, not model evidence |
| `v9/human_organoid/reports/donor_diagnostics/` | diagnostic-only split stress test | donor/source/fate/disease are coupled |

## Metric Boundary

Primary draft task metrics:

- `macro_f1`
- `balanced_accuracy`
- `auroc`
- `calibration_error`

Diagnostic-only metrics:

- `de_direction_match`
- `signature_rank_correlation`
- `feature_effect_direction_match`
- `feature_effect_rank_correlation`
- `feature_effect_top_k_de_overlap`

Reason:

- the organoid task has only 42 samples;
- disease, source, organoid fate, microglia condition, and donor are coupled;
- DE references are public and useful, but they should score optional artifacts
  after model generation rather than shape primary task optimization.

## Completed Adapter Decisions

| Diagnostic Family | Best Current Artifact | Decision |
|---|---|---|
| Response signature | global source-transfer | keep as default diagnostic |
| Condition-aware response signature | microglia-matched source-transfer | retain secondary, not replacement |
| Disease+microglia response signature | shared-control partial | keep as negative partial diagnostic |
| Feature effect | L2 logistic gene-space | keep as default feature-effect diagnostic |
| Feature-effect calibration | hypergeometric top-k null | keep in scorer output |
| PCA-LR reconstructed feature effect | PCA-LR pilot | keep as optional negative comparison |

## Draft Alpha Claim Language

Safe claim:

> The human organoid extension includes a draft, provenance-backed pilot task
> over OSD-863/OSD-871 neural organoid RNA-seq, with classification baselines
> and non-primary response-signature and feature-effect diagnostics.

Unsafe claim:

> The organoid task provides a frozen leaderboard-grade biological mechanism
> benchmark.

## Next Active Lane

Return to multispecies expansion.

Reason:

- the organoid diagnostic family now has a coherent draft-alpha boundary;
- the user's original expansion question included both organoids and
  non-mouse species;
- `V9-MULTI-001` and `V9-MULTI-002` already created a source inventory and
  feature-namespace strategy;
- the next unfinished scientific breadth gain is to turn the multispecies
  inventory into source-specific task candidates.

## Next Block

`V9-MULTI-003: Multispecies candidate source deep audit`

Expected work:

- inspect the existing multispecies source inventory;
- add per-source evidence summaries for OSD-207, OSD-37, and OSD-120;
- classify each candidate as species-native classification, response-signature,
  pathway/NES bridge, or defer;
- identify missing data needed before task manifests can be generated.

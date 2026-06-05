---
title: SpaceBio-Bench Evaluation Card
page_type: evaluation_card
status: public_review_ready
last_reviewed: 2026-06-04
claim_boundary: benchmark_evaluation_card_draft_no_new_result_claim
---

# SpaceBio-Bench Evaluation Card

## Evaluation Purpose

This card documents how SpaceBio-Bench evaluations should be interpreted. It
separates task validity, fold structure, metric reporting, baseline status, and
claim boundaries so that benchmark scores are not overread as biological
mechanism, translational readiness, or model superiority claims.

This card does not introduce new results. It summarizes evaluation evidence
already recorded in the v7.1 canonical result surface and the v9 public bulk
metadata-alpha scaffold.

Branch note: on the default `main` branch, v9-specific evidence paths such as
`v9/...` and `docs/V9_*` refer to artifacts maintained on the canonical `v3`
branch.

## Evaluation Surfaces

| Surface | Current status | Evaluation unit | Result boundary |
|---|---|---|---|
| v1-v7 / v7.1 canonical surface | Canonical historical result surface | Tissue, method, feature type, held-out validation, and FM snapshot rows | Cross-mission transcriptomics benchmark summary |
| v9 public bulk alpha | Metadata-only alpha scaffold | Task manifest, LOMO fold, baseline run, prediction row, run manifest | Workflow and provenance evidence, not leaderboard ranking |
| v9 draft extension lanes | Diagnostic or feasibility scaffolds | Asset inventory, metric spec, payload audit, draft task manifest, or diagnostic run | Draft-only; no public benchmark score claim |

## Evaluation Flow

```mermaid
flowchart LR
  A["Source inventory"] --> B["Task manifest"]
  B --> C["Held-out mission fold"]
  C --> D["Baseline or submitted predictions"]
  D --> E["Metrics with task/fold ids"]
  E --> F["Per-task interpretation"]
  F --> G["Pooled summary with caveats"]
  G --> H["Claim register language"]
```

The evaluation flow is intentionally claim-aware. A score is first interpreted
at the task and fold level, then summarized only with caveats about mission,
tissue, payload, baseline, and release-surface boundaries.

## Current v9 Public Bulk Evaluation Unit

The current v9 public bulk lane evaluates mission-held-out classification tasks
for public mouse bulk RNA-seq sources. The unit hierarchy is:

| Unit | Meaning | Evidence |
|---|---|---|
| Task | Tissue-specific bulk LOMO task with source ids, missions, feature namespace, and metric ids | `v9/task_manifest_index.csv`; `v9/task_manifests/*.json` |
| Fold | One held-out mission with train/test row counts and selected-gene counts | `v9/task_data_index.csv` |
| Baseline run | One baseline family evaluated on one task | `v9/reports/bulk_lomo_baseline_summary.csv` |
| Prediction file | Per-sample prediction output for a task/baseline run | per-baseline `predictions.csv` |
| Metrics file | Task-level metric output for a baseline run | per-baseline `metrics.json` |
| Run manifest | Provenance for a baseline execution | per-baseline `run_manifest.json` |

## Current Task Inventory

The v9 public bulk metadata-alpha scaffold currently indexes eight generated
bulk LOMO task manifests and 33 fold definitions:

| Task | Tissue | Variant | Missions | Folds | Sources |
|---|---|---:|---:|---:|---:|
| `A1_liver_bulk_lomo` | liver | canonical | 6 | 6 | 6 |
| `A1_liver_bulk_lomo_combat` | liver | combat | 6 | 6 | 6 |
| `A1_liver_bulk_lomo_iss_only` | liver | iss_only | 5 | 5 | 5 |
| `A2_gastrocnemius_bulk_lomo` | gastrocnemius | canonical | 3 | 3 | 3 |
| `A3_kidney_bulk_lomo` | kidney | canonical | 3 | 3 | 3 |
| `A4_thymus_bulk_lomo` | thymus | canonical | 4 | 4 | 3 |
| `A5_skin_bulk_lomo` | skin | canonical | 3 | 3 | 4 |
| `A6_eye_bulk_lomo` | eye | canonical | 3 | 3 | 3 |

These tasks are part of a metadata-only alpha snapshot. They do not imply a
frozen payload release.

## Metrics

Current v9 public bulk task manifests list these metric ids:

| Metric | Interpretation | Reporting note |
|---|---|---|
| `macro_f1` | Class-balanced F1 summary across labels | Report per task; small folds can be unstable |
| `balanced_accuracy` | Mean sensitivity across classes | Useful under class imbalance |
| `auroc` | Rank-based discrimination between flight/control labels | Can look strong even when calibration or fold reliability is weak |
| `calibration_error` | Probability calibration diagnostic | Lower is better; compare only across matched output formats |
| `mission_discrimination` | Diagnostic for mission-correlated structure where embeddings are available | Not emitted for every baseline family |

Metric reporting should include the task id, fold family, baseline id, variant,
and whether the result comes from a canonical result surface, metadata alpha, or
draft diagnostic lane.

## Current v9 Baseline Status

The current v9 public bulk scaffold includes 24 evaluated baseline rows across
three simple baseline families:

| Baseline family | Rows | Current interpretation |
|---|---:|---|
| `logistic_regression_l2` | 8 | Simple linear workflow baseline |
| `nearest_centroid` | 8 | Simple centroid-based diagnostic baseline |
| `pca_logistic_regression` | 8 | PCA feature-compression plus logistic baseline |

Mean metrics across the eight task rows are recorded in the v9 dataset card:

| Baseline | Macro-F1 | Balanced accuracy | AUROC | Calibration error | Mission discrimination |
|---|---:|---:|---:|---:|---:|
| Logistic regression L2 | 0.5532 | 0.5808 | 0.6870 | 0.3439 | NA |
| Nearest centroid | 0.5383 | 0.5685 | 0.6321 | 0.1132 | 0.8733 |
| PCA logistic regression | 0.5353 | 0.5619 | 0.6447 | 0.3747 | 0.9190 |

These are scaffold baselines. They validate the task/evaluation plumbing and
provide anchors for future methods, but they are not tuned leaderboard
endpoints.

## What This Enables

The evaluation card gives readers a compact way to check whether a score is
being read at the right level: task, fold, baseline run, or pooled summary. It
also makes clear when a number is a workflow sanity check rather than a model
ranking.

## Reading Order For Results

Read results in this order:

1. Release surface and claim boundary.
2. Task and fold definition.
3. Source provenance and payload verification status.
4. Per-fold and per-task metrics.
5. Baseline run manifest and prediction count.
6. Pooled summary, only after the above checks.

Do not start with a pooled mean alone. Pooled summaries can hide fold-level
failure, tissue-specific instability, mission-label confounding, or variant
differences.

## Leakage And Confounding Controls

Current controls:

- Mission-held-out fold definitions are explicit in `v9/task_data_index.csv`.
- Task-source relationships are explicit in `v9/task_manifest_index.csv`.
- v7.1 public wording requires subset notes when comparing mixed FM surfaces.
- v9 public bulk alpha separates public bulk tasks from draft organoid,
  single-cell, and multispecies lanes.
- v9 alpha documents that payload-level hash verification is pending.

Controls that remain future or lane-specific:

- Payload-level hash verification for every distributed fold matrix.
- Release-time revalidation of feature-selection and preprocessing leakage
  controls for any frozen v9 payload bundle.
- Adapter-specific validation for foundation models and virtual-cell models.
- Human or biological validation for mechanistic interpretation claims.

## Failure Modes To Watch

| Failure mode | Why it matters | Mitigation |
|---|---|---|
| Pooled-score overclaim | A high mean can hide poor tissue or mission performance | Always report per-task metrics |
| Mixed-surface FM comparison | Rows may come from different tissues, adapters, or versions | Label the reported surface for every row |
| Mission confounding | Mission labels can encode hardware, vehicle, age, protocol, or processing | Treat mission shift as benchmark pressure, not pure biology |
| Payload-boundary drift | Metadata alpha can be mistaken for frozen payload release | Keep payload hash status in every release-facing card |
| Extension-lane leakage | Draft organoid, single-cell, or multispecies lanes can be pulled into public bulk claims | Keep release target and claim boundary explicit |
| Biological interpretation overreach | Classification performance is not mechanistic proof | Separate benchmark evidence from mechanistic validation |

## Boundary Summary

The current evaluation surfaces should not be read as evidence for:

- A frozen v9 payload release.
- A state-of-the-art model leaderboard.
- A uniform 8-tissue foundation-model comparison.
- Clinical, crew-health, countermeasure, intervention, or Mars-regime validity.
- Biological mechanism proof from classifier performance.
- Payload-level integrity beyond the currently documented checksum-manifest
  evidence.

## Reproducibility Evidence

Local evidence:

- `docs/CANONICAL_RESULTS_V7_1.md`
- `docs/v9_hf_dataset_card.md`
- `v9/task_manifest_index.csv`
- `v9/task_data_index.csv`
- `v9/source_inventory.csv`
- `v9/source_checksum_audit.csv`
- `v9/reports/bulk_lomo_baseline_summary.csv`
- `v9/reports/nearest_centroid/bulk_lomo_summary.csv`
- `v9/reports/sklearn_baselines/bulk_lomo_summary.csv`
- `v9/datapackage.draft.json`

Companion cards:

- `docs/SPACEBIOBENCH_SYSTEM_CARD.md`
- `docs/SPACEBIOBENCH_CLAIM_REGISTER.md`
- `docs/SPACEBIOBENCH_RELEASE_READINESS_CARD.md`

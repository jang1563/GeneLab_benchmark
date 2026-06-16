---
title: SpaceBio-Bench Evaluation Card
page_type: evaluation_card
status: public_ready
last_reviewed: 2026-06-16
---

# SpaceBio-Bench Evaluation Card

## Purpose

This card explains how to read SpaceBio-Bench evaluations. It separates task
definition, fold structure, metric reporting, baseline rows, and pooled
summaries so readers can understand what a score means.

## Evaluation Surfaces

| Surface | Evaluation unit | Public interpretation |
|---|---|---|
| v1-v7 / v7.1 canonical surface | Tissue, method, feature type, held-out validation, and foundation-model snapshot rows | Historical cross-mission transcriptomics benchmark summary |
| Hugging Face public fold package | Public LOMO fold with matrices, labels, metadata, and selected genes | Processed data package for reproducible task-level evaluation |
| v9 public bulk catalog | Task manifest, LOMO fold, baseline run, prediction row, metric file, and run record | Metadata catalog and reference baseline surface |

## Evaluation Flow

| Stage | Evidence to inspect | Why it matters |
|---|---|---|
| 1. Source rows | OSDR accessions, tissue labels, mission labels, and access status | Identifies the public data behind a task |
| 2. Task manifest | Task ID, tissue, feature namespace, label map, and metric IDs | Defines what the evaluation is testing |
| 3. Held-out mission fold | Train/test mission split, row counts, and selected-gene counts | Separates mission-held-out validation from random-split performance |
| 4. Prediction and metric files | Predictions, task/fold IDs, AUROC, macro-F1, balanced accuracy, and calibration | Ties every score to a concrete task and fold |
| 5. Per-task interpretation | Tissue-specific and fold-specific behavior | Preserves mission and tissue variability |
| 6. Pooled summary | Aggregate metrics after task and fold checks | Gives a navigation-level summary |

## v9 Public Bulk Evaluation Unit

The v9 public bulk catalog organizes mission-held-out classification tasks for
public mouse bulk RNA-seq sources.

| Unit | Meaning | Public files |
|---|---|---|
| Task | Tissue-specific bulk LOMO task with missions, labels, feature namespace, and metric IDs | `v9/task_manifest_index.csv`; `v9/task_manifests/*.json` |
| Fold | One held-out mission with train/test row counts and selected-gene counts | `v9/task_data_index.csv` |
| Baseline run | One baseline family evaluated on one task | `v9/reports/bulk_lomo_baseline_summary.csv` |
| Prediction file | Per-sample prediction output for a task and baseline run | per-baseline `predictions.csv` |
| Metrics file | Task-level metric output for a baseline run | per-baseline `metrics.json` |
| Run record | Execution metadata for a baseline run | per-baseline `run_manifest.json` |

## Current Task Inventory

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

## Metrics

| Metric | Interpretation | Reporting note |
|---|---|---|
| `macro_f1` | Class-balanced F1 summary across labels | Report per task; small folds can be unstable |
| `balanced_accuracy` | Mean sensitivity across classes | Useful under class imbalance |
| `auroc` | Rank-based discrimination between flight/control labels | Compare with calibration and per-fold behavior |
| `calibration_error` | Probability calibration diagnostic | Lower is better; compare matched output formats |
| `mission_discrimination` | Diagnostic for mission-correlated structure where embeddings are available | Not emitted for every baseline family |

Metric reports should include task ID, fold ID, baseline or method ID, variant,
and release surface.

## Current v9 Baseline Rows

| Baseline family | Rows | Interpretation |
|---|---:|---|
| `logistic_regression_l2` | 8 | Simple linear workflow baseline |
| `nearest_centroid` | 8 | Simple centroid-based baseline |
| `pca_logistic_regression` | 8 | PCA feature-compression plus logistic baseline |

Mean metrics across the eight task rows:

| Baseline | Macro-F1 | Balanced accuracy | AUROC | Calibration error | Mission discrimination |
|---|---:|---:|---:|---:|---:|
| Logistic regression L2 | 0.5532 | 0.5808 | 0.6870 | 0.3439 | NA |
| Nearest centroid | 0.5383 | 0.5685 | 0.6321 | 0.1132 | 0.8733 |
| PCA logistic regression | 0.5353 | 0.5619 | 0.6447 | 0.3747 | 0.9190 |

These rows provide reproducible reference anchors. Use per-task rows together
with pooled means when comparing methods.

## Reading Order For Results

1. Release surface.
2. Task and fold definition.
3. Source rows and access status.
4. Per-fold and per-task metrics.
5. Baseline run record and prediction count.
6. Pooled summary.

Starting from the task and fold keeps tissue-specific behavior, mission-label
structure, and variant differences visible.

## Leakage And Confounding Controls

- Mission-held-out fold definitions are explicit in `v9/task_data_index.csv`.
- Task-source relationships are explicit in `v9/task_manifest_index.csv`.
- v7.1 public wording includes subset notes for mixed foundation-model
  comparisons.
- v9 public bulk tasks are separated from extension workspaces for other
  modalities.
- Feature selection and preprocessing should be checked within each training
  split for any newly packaged fold.

## Failure Modes To Watch

| Failure mode | Why it matters | Mitigation |
|---|---|---|
| Pooled-score overread | A high mean can hide weak tissue or mission rows | Report per-task metrics |
| Mixed-surface comparison | Rows can come from different tissues, adapters, or versions | Label the release surface for every row |
| Mission confounding | Mission labels can encode hardware, vehicle, age, protocol, or processing | Interpret mission shift as the benchmark pressure |
| Biological overinterpretation | Classification performance is not mechanistic validation | Separate benchmark scores from mechanistic follow-up |

## Reproducibility Files

- `docs/CANONICAL_RESULTS_V7_1.md`
- `docs/hf_dataset_card.md`
- `docs/v9_hf_dataset_card.md`
- `v9/task_manifest_index.csv`
- `v9/task_data_index.csv`
- `v9/source_inventory.csv`
- `v9/source_checksum_audit.csv`
- `v9/reports/bulk_lomo_baseline_summary.csv`
- `v9/reports/nearest_centroid/bulk_lomo_summary.csv`
- `v9/reports/sklearn_baselines/bulk_lomo_summary.csv`
- `v9/datapackage.draft.json`

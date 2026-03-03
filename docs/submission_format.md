# GeneLab Benchmark — Submission Format Specification

**Version**: 1.1
**Date**: 2026-03-01
**Status**: Draft

---

## Overview

External models submit predictions as a **JSON file**. Each submission covers one task (e.g., `A4` = Spaceflight Detection – Thymus) and includes predicted Flight probabilities for every test sample in every fold.

The evaluation server computes AUROC, permutation p-value, and CI automatically using `scripts/evaluate_submission.py`.
If a task ID maps to multiple task directories (for example `A1`), pass `--task-dir` explicitly when validating/evaluating.

---

## File Format

### Top-Level Structure

```json
{
  "task_id": "A4",
  "model_name": "MyModel_v1",
  "model_description": "Optional: brief description of model architecture",
  "submission_date": "2026-03-01",
  "predictions": {
    "<fold_id>": {
      "<sample_id>": <flight_probability>
    }
  }
}
```

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| `task_id` | string | ✓ | Task identifier: `A2`, `A4`, `B2`, `B4`, etc. |
| `model_name` | string | ✓ | Unique model identifier (alphanumeric + underscore, max 50 chars) |
| `model_description` | string | — | Optional free-text description |
| `submission_date` | string | — | ISO date format (YYYY-MM-DD) |
| `predictions` | dict | ✓ | Nested dict: fold_id → {sample_id → probability} |

### Predictions Structure

- **fold_id**: Must match the fold/pair directory name — format depends on category:
  - **Category A (LOMO)**: `fold_{MISSION}_test` (e.g., `fold_RR-9_test`)
  - **Category B (Cross-Mission Transfer)**: `pair_{TRAIN}_{TEST}` (e.g., `pair_MHU-2_RR-6`)
- **sample_id**: Must match the row index in `test_X.csv` exactly (from the corresponding directory)
- **flight_probability**: Float in `[0.0, 1.0]` — probability of being a Flight sample

> **Held-out folds** (`fold_*_holdout`): Labels are private. These are *optional* in Tier 1 submissions — the evaluator will note if holdout predictions are present/absent but does not require them.

---

## Example: Task A4 (Thymus LOMO, 4 folds)

```json
{
  "task_id": "A4",
  "model_name": "Geneformer_finetuned_v1",
  "model_description": "Geneformer (6-layer) fine-tuned on A4 train splits",
  "submission_date": "2026-03-01",
  "predictions": {
    "fold_MHU-1_test": {
      "MHU1_TMS_FLT_M1": 0.92,
      "MHU1_TMS_FLT_M2": 0.88,
      "MHU1_TMS_FLT_M3": 0.95,
      "MHU1_TMS_GC_M1":  0.07,
      "MHU1_TMS_GC_M2":  0.12,
      "MHU1_TMS_GC_M3":  0.03
    },
    "fold_MHU-2_test": {
      "MHU2_TMS_FLT_M1": 0.81,
      "MHU2_TMS_FLT_M2": 0.79,
      "MHU2_TMS_FLT_M3": 0.93,
      "MHU2_TMS_GC_M1":  0.18,
      "MHU2_TMS_GC_M2":  0.22,
      "MHU2_TMS_GC_M3":  0.09
    },
    "fold_RR-6_test": {
      "RR6_TMS_FLT_F1":  0.75,
      "...": "...",
      "RR6_TMS_GC_F17":  0.31
    },
    "fold_RR-9_test": {
      "RR9_TMS_FLT_F1":  0.98,
      "...": "...",
      "RR9_TMS_GC_F10":  0.04
    }
  }
}
```

---

## Example: Task B4 (Thymus Cross-Mission, 12 pairs)

```json
{
  "task_id": "B4",
  "model_name": "PCA-LR_baseline_B4",
  "model_description": "StandardScaler + PCA(50) + LogReg(L2, lbfgs), train-variance gene filter",
  "submission_date": "2026-03-01",
  "predictions": {
    "pair_MHU-1_MHU-2": {
      "MHU2_TMS_FLT_M1": 0.83,
      "MHU2_TMS_FLT_M2": 0.76,
      "MHU2_TMS_GC_M1":  0.21,
      "MHU2_TMS_GC_M2":  0.18
    },
    "pair_MHU-1_RR-6": {
      "RR6_TMS_FLT_F1":  0.91,
      "...": "...",
      "RR6_TMS_GC_F17":  0.09
    },
    "pair_MHU-2_MHU-1": {
      "MHU1_TMS_FLT_M1": 0.87,
      "...": "..."
    },
    "...": "... (12 pairs total)"
  }
}
```

Key difference from Category A:
- `fold_id` uses `pair_{TRAIN}_{TEST}` format (not `fold_{MISSION}_test`)
- All N×(N-1) directed pairs must be included (12 for B4, 6 for B2)

---

## Available Tasks (v1.0)

### Category A — Spaceflight Detection (LOMO)

| Task ID | Tissue | Folds | Test Samples per Fold | Input Features |
|---------|--------|-------|----------------------|----------------|
| `A2` | Gastrocnemius | 3 | 8–12 | ~21k genes (log2 normalized) |
| `A4` | Thymus | 4 | 6–35 | ~21k genes (log2 normalized) |

### Category B — Cross-Mission Transfer

| Task ID | Tissue | Missions | Pairs | Task Directory |
|---------|--------|---------|-------|----------------|
| `B4` | Thymus | MHU-1, MHU-2, RR-6, RR-9 | 12 | `tasks/B4_thymus_cross_mission/` |
| `B2` | Gastrocnemius | RR-1, RR-5, RR-9 | 6 | `tasks/B2_gastrocnemius_cross_mission/` |

Each pair directory is named `pair_{TRAIN}_{TEST}` (e.g., `pair_MHU-2_RR-6`).

Category B evaluation does **not** report a single GO/NO-GO. Instead, Transfer Pattern Summary (DD-17) is reported:
- AUROC ≥ 0.70 pairs count
- perm_p < 0.05 pairs count
- Large pairs (n≥10) mean AUROC + significance

*(Full AUROC matrix: see `processed/B_cross_mission/` and PHASE1_RESULTS.md)*

---

## Data Access

All input features and public labels are located in the `tasks/` directory:

```
tasks/
  A4_thymus_lomo/                     ← Category A task
    task_info.json                    ← task metadata
    fold_MHU-1_test/
      train_X.csv                     ← training features (samples × genes)
      train_y.csv                     ← training labels (1=Flight, 0=Ground)
      train_meta.csv                  ← training sample metadata
      test_X.csv                      ← test features (PUBLIC)
      test_y.csv                      ← test labels (PUBLIC for LOMO; PRIVATE for held-out)
      test_meta.csv                   ← test sample metadata
      selected_genes.txt              ← gene list (top 75th percentile variance)
    fold_MHU-2_test/
      ...
  B4_thymus_cross_mission/            ← Category B task
    pair_MHU-1_MHU-2/
      test_y.csv                      ← test labels (PUBLIC)
      ← test_X comes from the source tissue CSV (see PHASE1_RESULTS.md)
    pair_MHU-1_RR-6/
      ...
    pair_MHU-2_MHU-1/
      ...
    ... (12 pairs for B4)
  B2_gastrocnemius_cross_mission/     ← Category B task
    pair_RR-1_RR-5/
      test_y.csv
    ... (6 pairs for B2)
```

> **Category B feature access**: Train/test feature splits for each pair are computed on-the-fly from `processed/B_cross_mission/{tissue}/{tissue}_all_missions_log2_norm.csv` using the mission column. See `scripts/cross_mission_transfer.py` for the exact split logic and variance filter (train-only).

**Gene format**: Columns are Ensembl gene IDs (e.g., `ENSMUSG00000019773`).
**Label encoding**: `1.0` = Flight, `0.0` = Ground/Vivarium Control. Basal Control (BC) and Aseptic Ground Control (AG) samples are excluded.

---

## Evaluation Metrics

The evaluator (`scripts/evaluate_submission.py`) computes:

| Metric | Description |
|--------|-------------|
| **AUROC** | Area under ROC curve (per fold + mean ± SD) |
| **95% CI** | Bootstrap CI (N=2000) on mean AUROC |
| **perm_p** | Permutation p-value (N=1000, pseudocount applied) |
| **Go/No-Go** | AUROC > 0.700 AND CI lower > 0.500 |

---

## Submission Rules

1. **One submission per task**: Submit separate JSON files for A2, A4, B-tasks.
2. **All folds required**: Every fold in a task must have predictions for all test samples.
3. **No training label leakage**: Models may only use features from `train_X.csv` and `train_y.csv` within each fold (LOMO design — do not train on other folds' test labels).
4. **Probability required**: Submit calibrated probabilities, not binary labels. Binary 0/1 will be accepted but AUROC will be limited.
5. **Model name unique**: Each submission must have a unique `model_name`. Same name = overwrite.
6. **Text LLM submissions (Track B)**: Use the text-based input format defined in DD-16 (`docs/text_llm_format.md`). Same JSON output format applies.
7. **Ambiguous task IDs**: If multiple task directories exist for one task ID (e.g., `A1_liver_lomo`, `A1_liver_lomo_combat`, `A1_liver_lomo_iss_only`), evaluation requires `--task-dir`.

---

## Leaderboard Structure (Planned)

```
Rank | Model Name            | Track  | A2 AUROC | A4 AUROC | Mean  | perm_p
-----|-----------------------|--------|----------|----------|-------|-------
1    | Geneformer_ft_v1      | Tier 2 | 0.941    | 0.961    | 0.951 | 0.001
2    | PCA-LR_baseline       | Tier 1 | 0.907    | 0.923    | 0.915 | 0.001
3    | GPT4o_zeroshot        | Tier 3 | 0.712    | 0.698    | 0.705 | 0.018
```

Track classification:
- **Tier 1**: Classical ML (LR, RF, XGBoost, PCA-LR)
- **Tier 2**: Gene Expression Foundation Models (Geneformer, scGPT)
- **Tier 3**: Text-based LLMs (GPT-4o, Claude, Llama 3)

---

## Validation

Before submitting, validate your JSON with the provided schema checker:

```bash
python scripts/evaluate_submission.py \
    --submission path/to/your_submission.json \
    --task A4 \
    --validate-only

# A1 example (variant must be explicit)
python scripts/evaluate_submission.py \
    --submission path/to/your_A1_submission.json \
    --task A1 \
    --task-dir A1_liver_lomo \
    --validate-only
```

---

## Contact / Issues

This benchmark is part of the GeneLab_benchmark project based on NASA OSDR data.
GitHub Issues: [TBD when public]
Data source: [NASA OSDR](https://osdr.nasa.gov/bio/repo/)

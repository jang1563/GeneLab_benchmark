# tasks/ — GeneLab Benchmark Task Inputs

This directory contains the **public task inputs** for the GeneLab Benchmark (v1.0).
Each subdirectory is a self-contained task with all features and labels needed to train and evaluate a model.

---

## Directory Overview

```
tasks/
├── A2_gastrocnemius_lomo/          ← Category A, Task A2 (Gastrocnemius LOMO — GO)
├── A4_thymus_lomo/                 ← Category A, Task A4 (Thymus LOMO — GO)
├── A5_skin_lomo/                   ← Category A, Task A5 (Skin LOMO — GO)
├── A6_eye_lomo/                    ← Category A, Task A6 (Eye LOMO — GO, pathway-level)
├── A1_liver_lomo/                  ← Category A, Task A1 (Liver LOMO — NO-GO)
├── A1_liver_lomo_combat/           ← A1 sensitivity: ComBat-seq batch correction (Category J)
├── A1_liver_lomo_iss_only/         ← A1 sensitivity: ISS missions only subset
├── A3_kidney_lomo/                 ← Category A, Task A3 (Kidney LOMO — NO-GO)
├── B2_gastrocnemius_cross_mission/ ← Category B, Task B2 (Gastrocnemius Cross-Mission)
├── B4_thymus_cross_mission/        ← Category B, Task B4 (Thymus Cross-Mission)
├── B5_skin_cross_mission/          ← Category B, Task B5 (Skin Cross-Mission)
├── B6_eye_cross_mission/           ← Category B, Task B6 (Eye Cross-Mission)
├── B1_liver_cross_mission/         ← Category B, Task B1 (Liver Cross-Mission)
└── B3_kidney_cross_mission/        ← Category B, Task B3 (Kidney Cross-Mission)
```

---

## Category A — Spaceflight Detection (LOMO)

**Goal**: Binary classification (Flight vs. Ground) using Leave-One-Mission-Out CV.

### Folder Structure

```
A4_thymus_lomo/
  task_info.json           ← Task metadata (missions, fold sizes, gene counts)
  fold_RR-9_test/          ← One LOMO fold: test mission = RR-9
    train_X.csv            ← Training features (samples × genes), log2 normalized
    train_y.csv            ← Training labels (1.0=Flight, 0.0=Ground)
    train_meta.csv         ← Training sample metadata (mission, group, etc.)
    test_X.csv             ← Test features (PUBLIC)
    test_y.csv             ← Test labels (PUBLIC for LOMO folds)
    test_meta.csv          ← Test sample metadata
    selected_genes.txt     ← Gene list: top 75th percentile variance (train-only)
    fold_info.json         ← Per-fold statistics
    geneformer_tokens/     ← Tokenized inputs for Geneformer (Tier 2)
  fold_MHU-1_test/
    ...
  fold_RR-23_holdout/      ← Held-out fold: labels PRIVATE (not for public eval)
    train_X.csv
    test_X.csv
    (test_y.csv absent — private)
```

### Available A Tasks

| Task | Tissue | Folds | n samples (binary) | Status | Notes |
|------|--------|-------|-------------------|--------|-------|
| `A2` | Gastrocnemius | 3 | 32 | ✓ GO | LOMO: RR-1, RR-5, RR-9 |
| `A4` | Thymus | 4 | 67 | ✓ GO | LOMO: MHU-1†, MHU-2, RR-6, RR-9 |
| `A5` | Skin | 3 | 102 | ✓ GO | LOMO: MHU-2§, RR-6, RR-7 |
| `A6` | Eye | 3 | 37 | ✓ GO (pathway) | Gene-level CI fails; GSVA Hallmark AUROC=0.915 |
| `A1` | Liver | 6 | 193 | ✗ NO-GO | Pipeline heterogeneity |
| `A3` | Kidney | 3 | 118 | ✗ NO-GO | Low AUROC (0.593); pathway also fails |

†MHU-1 = Track 2b (GC strain = C57BL/6CR, not C57BL/6J). See [PHASE1_RESULTS.md](../docs/development_history/PHASE1_RESULTS.md).
§MHU-2 = OSD-238 (dorsal) + OSD-239 (femoral) merged; RR-7 = OSD-254 C57BL/6J non-BSL subset (n=30).

> **A1 variants (sensitivity analysis only, not in primary benchmark):**
> `A1_liver_lomo_combat/` — same as A1 but with ComBat-seq batch correction applied (Category J analysis).
> `A1_liver_lomo_iss_only/` — A1 restricted to ISS missions only (excludes ground-only missions).
> Both variants remain NO-GO and are provided for reproducibility. Not included in HuggingFace upload.

### Reading a Fold

```python
import pandas as pd

task = "A4_thymus_lomo"
fold = "fold_RR-9_test"
base = f"tasks/{task}/{fold}"

train_X = pd.read_csv(f"{base}/train_X.csv", index_col=0)  # (n_train, n_genes)
train_y = pd.read_csv(f"{base}/train_y.csv", index_col=0)  # (n_train, 1)
test_X  = pd.read_csv(f"{base}/test_X.csv",  index_col=0)  # (n_test, n_genes)
test_y  = pd.read_csv(f"{base}/test_y.csv",  index_col=0)  # (n_test, 1)

# Columns: Ensembl gene IDs (e.g., ENSMUSG00000019773)
# Labels:  1.0 = Flight,  0.0 = Ground/Vivarium Control
print(train_X.shape, train_y.iloc[:,0].value_counts().to_dict())
```

---

## Category B — Cross-Mission Transfer

**Goal**: Train on one mission, evaluate generalization to another mission it has never seen.

### Folder Structure

```
B4_thymus_cross_mission/
  pair_MHU-2_RR-6/         ← Train=MHU-2, Test=RR-6
    test_y.csv             ← Test labels (PUBLIC)
  pair_MHU-2_RR-9/
    test_y.csv
  pair_RR-9_MHU-2/         ← Train=RR-9, Test=MHU-2 (directed — different from above)
    test_y.csv
  ... (12 pairs total for B4)
```

> **Note**: Category B task directories contain only `test_y.csv`. Train and test features are
> generated on-the-fly from the processed log2-normalized expression matrix
> (`processed/B_cross_mission/{tissue}/{tissue}_all_missions_log2_norm.csv`)
> by splitting on the `mission` column. See `scripts/cross_mission_transfer.py` for the split
> logic and variance filter (computed on train only, matching DD-03).

### Available B Tasks

| Task | Tissue | Missions | Pairs | Mean AUROC (PCA-LR) |
|------|--------|---------|-------|---------------------|
| `B4` | Thymus | MHU-1, MHU-2, RR-6, RR-9 | 12 | 0.860 [0.763, 0.953] |
| `B5` | Skin | MHU-2, RR-6, RR-7 | 6 | 0.772 [0.691, 0.834] |
| `B6` | Eye | RR-1, RR-3, TBD | 6 | 0.754 [0.688, 0.838] |
| `B2` | Gastrocnemius | RR-1, RR-5, RR-9 | 6 | 0.801 [0.653, 0.944] |
| `B1` | Liver | MHU-2, RR-1, RR-3, RR-6, RR-8, RR-9 | 30 | 0.577 [0.492, 0.666] |
| `B3` | Kidney | RR-1, RR-3, RR-7 | 6 | 0.555 [0.397, 0.681] |

Pair naming: `pair_{TRAIN}_{TEST}` — directed (train→test order matters).

> **B3 kidney note**: Some kidney samples in the source data (`processed/B_cross_mission/kidney/`) carry
> the "BAL-TAL" strain label — a pre-existing annotation in NASA OSDR metadata, not introduced by this
> pipeline. A3 is already NO-GO; B3 AUROC (0.555) reflects this mixed-strain confound. Use B3 results
> with caution for cross-mission transfer conclusions.

### Evaluation (DD-17)

Category B does **not** produce a single GO/NO-GO. Use Transfer Pattern Summary:
- Count of AUROC ≥ 0.70 pairs
- Count of perm_p < 0.05 pairs
- Large pairs (n≥10) mean AUROC + significance (primary statistical unit)

```bash
python scripts/evaluate_submission.py \
    --submission evaluation/submission_PCALR_baseline_B4.json \
    --task B4

python scripts/evaluate_submission.py \
    --submission evaluation/submission_PCALR_baseline_B5.json \
    --task B5
```

---

## Labels

| Value | Meaning |
|-------|---------|
| `1.0` | Flight (spaceflight / microgravity condition) |
| `0.0` | Ground (vivarium control) |
| — | Basal Control (BC) and Aseptic Ground Control (AG) samples are excluded from all tasks |

---

## Gene Features

- **Format**: Log2(DESeq2 size-factor normalized counts + 1), per sample
- **Columns**: Ensembl mouse gene IDs (e.g., `ENSMUSG00000019773`)
- **Selection**: Top 75th percentile by variance, computed on training set only per fold (DD-03)
- **Typical count**: ~21,000 genes after variance filter

---

## Submission

Prepare predictions following [docs/submission_format.md](../docs/submission_format.md) and evaluate with:

```bash
python scripts/evaluate_submission.py --submission <your_file.json> --task <A1|A2|A3|A4|A5|A6|B1|B2|B3|B4|B5|B6>
```

Baseline submissions are in `evaluation/`.

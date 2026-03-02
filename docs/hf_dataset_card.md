---
license: cc-by-4.0
task_categories:
  - tabular-classification
tags:
  - genomics
  - transcriptomics
  - spaceflight
  - benchmarking
  - nasa-osdr
  - bulk-rna-seq
  - lomo
  - cross-mission-transfer
  - mouse
  - rna-seq
size_categories:
  - 1GB<n<10GB
language:
  - en
pretty_name: "GeneLab Spaceflight Transcriptomics Benchmark"
---

# GeneLab Spaceflight Transcriptomics Benchmark

**A public benchmark for evaluating AI/ML and Foundation Models on NASA OSDR spaceflight transcriptomics data.**

Version: v1.0-alpha | Code: [GitHub](https://github.com/jak4013/GeneLab_benchmark)

---

## Dataset Summary

GeneLab Benchmark provides standardized train/test splits for evaluating how well machine learning models — from classical baselines to gene expression foundation models (Geneformer, scGPT) and text-based LLMs — generalize **spaceflight transcriptomic signatures across missions**.

**Core challenge**: Train a model on one spaceflight mission's RNA-seq data. Can it classify samples from a different mission it has never seen?

**Data source**: NASA Open Science Data Repository (OSDR) — mouse multi-tissue bulk RNA-seq from ISS and ground control missions (C57BL/6J strain, Track 2a).

---

## Dataset Structure

This repository contains the **feature matrices** (train_X.csv, test_X.csv) for GO-status tasks. Labels, metadata, and fold structure files (train_y.csv, test_y.csv, fold_info.json, selected_genes.txt) are in the [GitHub repository](https://github.com/jak4013/GeneLab_benchmark).

```
genelab-benchmark/
├── A2_gastrocnemius_lomo/        ← Gastrocnemius: GO (AUROC=0.907)
│   ├── fold_RR-1_test/
│   │   ├── train_X.csv           ← Training features (samples × genes, log2-normalized)
│   │   └── test_X.csv            ← Test features
│   ├── fold_RR-5_test/
│   └── fold_RR-9_test/
│
├── A4_thymus_lomo/               ← Thymus: GO (AUROC=0.923)
│   ├── fold_MHU-1_test/
│   ├── fold_MHU-2_test/
│   ├── fold_RR-6_test/
│   └── fold_RR-9_test/
│
├── A5_skin_lomo/                 ← Skin: GO (AUROC=0.821)
│   ├── fold_MHU-2_test/
│   ├── fold_RR-6_test/
│   └── fold_RR-7_test/
│
└── A6_eye_lomo/                  ← Eye: GO pathway-level (GSVA Hallmark AUROC=0.915)
    ├── fold_RR-1_test/
    ├── fold_RR-3_test/
    └── fold_TBD_test/
```

---

## File Format

### Feature matrix (train_X.csv, test_X.csv)

- **Rows**: Sample IDs (e.g., `Mmus_C57-6J_SKN_FLT_25days_Rep1_F1`)
- **Columns**: Ensembl mouse gene IDs (e.g., `ENSMUSG00000021969`)
- **Values**: Log2(DESeq2 size-factor normalized counts + 1)
- **Gene selection**: Top 75th percentile variance, computed on **training missions only** (no test leakage — DD-03)
- **Typical shape**: (n_train × ~20,000 genes), (n_test × ~20,000 genes)

### Labels (in GitHub repo)

| Value | Meaning |
|-------|---------|
| `1` | Flight (spaceflight / microgravity condition) |
| `0` | Ground (vivarium control / ground control) |

Basal Control (BC) and Artificial Gravity (AG) samples are excluded from binary tasks.

---

## Task Statistics

### Category A — Spaceflight Detection (LOMO)

| Task | Tissue | Missions | Samples (binary) | Mean AUROC | Status |
|------|--------|---------|-----------------|-----------|--------|
| A2 | Gastrocnemius | RR-1, RR-5, RR-9 | 32 | 0.907 | ✓ GO |
| A4 | Thymus | MHU-1, MHU-2, RR-6, RR-9 | 67 | 0.923 | ✓ GO |
| A5 | Skin | MHU-2, RR-6, RR-7 | 102 | 0.821 | ✓ GO |
| A6 | Eye | RR-1, RR-3, TBD | 37 | 0.915† | ✓ GO (pathway) |

†A6 gene-level AUROC=0.811 (CI lower fails); pathway-level (GSVA Hallmark) AUROC=0.915 (all conditions pass).

---

## Downloading

### Option A: Python API (recommended)

```python
from huggingface_hub import hf_hub_download
import pandas as pd

# Download one fold's feature matrix
train_X = pd.read_csv(
    hf_hub_download(
        repo_id="jang1563/genelab-benchmark",
        filename="A5_skin_lomo/fold_RR-7_test/train_X.csv",
        repo_type="dataset",
    ),
    index_col=0
)
print(train_X.shape)  # (72, 20110)
```

### Option B: Download full task

```python
from huggingface_hub import snapshot_download

# Download all A5 skin files
snapshot_download(
    repo_id="jang1563/genelab-benchmark",
    repo_type="dataset",
    allow_patterns="A5_skin_lomo/**",
    local_dir="./data/benchmark",
)
```

### Option C: Clone GitHub + link HF data

```bash
git clone https://github.com/jak4013/GeneLab_benchmark
cd GeneLab_benchmark
python scripts/download_from_hf.py --task A5
```

---

## Evaluation

Evaluate your model predictions using the included evaluation script (GitHub):

```python
# Prepare submission JSON
submission = {
    "task_id": "A5",
    "model_name": "MyModel_v1",
    "predictions": {
        "fold_MHU-2_test": {"sample_id_1": 0.92, ...},
        "fold_RR-6_test":  {...},
        "fold_RR-7_test":  {...}
    }
}

# Evaluate (requires cloning GitHub repo)
# python scripts/evaluate_submission.py --submission my_submission.json --task A5
```

Metrics: Mean AUROC + 95% bootstrap CI + permutation p-value.
GO criteria: AUROC > 0.70 AND CI lower > 0.50 AND perm_p < 0.05.

---

## Source Data

All data derived from publicly available NASA OSDR datasets:

| OSD ID | Tissue | Mission | n (binary) |
|--------|--------|---------|-----------|
| OSD-101 | Gastrocnemius | RR-1 | 12 |
| OSD-401 | Gastrocnemius | RR-5 | 12 |
| OSD-326 | Gastrocnemius | RR-9 | 8 |
| OSD-289 | Thymus | MHU-1/MHU-2 | 12 |
| OSD-244 | Thymus | RR-6 | 35 |
| OSD-421 | Thymus | RR-9 | 20 |
| OSD-238 | Skin (dorsal) | MHU-2 | 18 |
| OSD-239 | Skin (femoral) | MHU-2 | 17 |
| OSD-243 | Skin | RR-6 | 37 |
| OSD-254 | Skin | RR-7 | 30 |
| OSD-100 | Eye | RR-1 | 12 |
| OSD-194 | Eye | RR-3 | 9 |
| OSD-397 | Eye | TBD | 16 |

---

## Preprocessing

1. DESeq2 size-factor normalization (per-mission)
2. log2(counts + 1) transformation
3. Global low-expression filter (≥20% samples with count > 1)
4. Top 75th percentile variance gene selection **per fold, train missions only** (DD-03)

---

## Citation

*(Manuscript in preparation)*

```bibtex
@dataset{kang2026genelab,
  title   = {GeneLab Spaceflight Transcriptomics Benchmark},
  author  = {Kang, Jaeyoung},
  year    = {2026},
  url     = {https://huggingface.co/datasets/jang1563/genelab-benchmark},
  note    = {v1.0-alpha}
}
```

Data source: NASA Open Science Data Repository (OSDR) — https://osdr.nasa.gov/bio/repo/

---

## License

Dataset: CC-BY-4.0
Code: MIT (see GitHub repository)
Source data: NASA OSDR public data

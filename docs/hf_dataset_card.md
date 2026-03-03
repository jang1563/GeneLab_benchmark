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
  - geneformer
  - foundation-model
  - pathway-analysis
size_categories:
  - 1GB<n<10GB
language:
  - en
pretty_name: "GeneLab Spaceflight Transcriptomics Benchmark"
---

# GeneLab Spaceflight Transcriptomics Benchmark

**A public benchmark for evaluating AI/ML and Foundation Models on NASA OSDR spaceflight transcriptomics data.**

Version: v1.0-alpha | Dataset freeze: 2026-03-01 | Code: [GitHub](https://github.com/jang1563/GeneLab_benchmark)

---

## Dataset Summary

GeneLab Benchmark provides standardized train/test splits for evaluating how well machine learning models — from classical baselines to gene expression foundation models (Geneformer, scGPT) and text-based LLMs — generalize **spaceflight transcriptomic signatures across missions**.

**Core challenge**: Train a model on one spaceflight mission's RNA-seq data. Can it classify samples from a different mission it has never seen?

**Data source**: NASA Open Science Data Repository (OSDR) — mouse multi-tissue bulk RNA-seq from ISS and ground control missions (C57BL/6J strain, Track 2a).

### Benchmark Scope

| Dimension | Coverage |
|---|---|
| Tissues | 6 (Liver, Gastrocnemius, Kidney, Thymus, Skin, Eye) |
| ISS Missions | 17 (RR-1 through RR-9, MHU-1, MHU-2, etc.) |
| Verified Studies | 24 OSD accessions |
| Samples | ~450 (binary Flight/Ground) |
| Task Categories | 7 (A–D, J, NC, Validation) |
| Evaluation Tasks | 25+ |

---

## Dataset Structure

This repository contains the **feature matrices** (train_X.csv, test_X.csv) for GO-status tasks. Labels, metadata, and fold structure files (train_y.csv, test_y.csv, fold_info.json, selected_genes.txt) are in the [GitHub repository](https://github.com/jang1563/GeneLab_benchmark).

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

## Phase 1 Results

### Category A — Spaceflight Detection (LOMO)

| Task | Tissue | Missions | Samples | Mean AUROC | Status |
|------|--------|---------|---------|-----------|--------|
| A4 | Thymus | MHU-1, MHU-2, RR-6, RR-9 | 67 | **0.923** | ✓ GO |
| A2 | Gastrocnemius | RR-1, RR-5, RR-9 | 32 | **0.907** | ✓ GO |
| A5 | Skin | MHU-2, RR-6, RR-7 | 102 | **0.821** | ✓ GO |
| A6 | Eye (pathway) | RR-1, RR-3, TBD | 37 | **0.915**† | ✓ GO |

†A6 gene-level AUROC=0.811 (CI lower fails); pathway-level (GSVA Hallmark) AUROC=0.915 (all conditions pass).

### Category B — Cross-Mission Transfer

| Tissue | Mean AUROC | 95% CI | N pairs | Tier |
|---|---|---|---|---|
| Thymus | **0.860** | [0.763, 0.953] | 12 | 1 |
| Gastrocnemius | **0.801** | [0.653, 0.944] | 6 | 1 |
| Skin | **0.772** | [0.691, 0.834] | 6 | 2 |
| Eye | 0.754 | [0.688, 0.838] | 6 | 2 |
| Liver | 0.577 | [0.492, 0.666] | 30 | 3 |
| Kidney | 0.555 | [0.397, 0.681] | 6 | 3 |

### Key Scientific Findings

| Hypothesis | Verdict | Evidence |
|---|---|---|
| H1: Liver most consistent cross-mission | **REFUTED** | Thymus (0.860) >> Liver (0.577) |
| H2: Transfer failure = biology, not batch | **SUPPORTED** | NES conservation r=0.9, D3 pathway F1=0.06 |
| H3: Pathway > gene conservation | **CONDITIONALLY SUPPORTED** | Kidney rescue (+0.31), Eye (+0.13), but tissue-pair dependent |

### Validation

- **Cell 2020 concordance**: 71.7% pathway direction match across 5 tissues (vs Beheshti et al., PMID 33242417)
- **Gene SHAP overlap**: 10.7% top-50 overlap (47× above random chance)
- **Negative controls**: NC1 permutation (0.50 ± 0.03) and NC2 housekeeping (0.49–0.55) both PASS

### Biological Enrichment (fGSEA Hallmark)

| Tissue | Top Pathways |
|---|---|
| Thymus | E2F_TARGETS, G2M_CHECKPOINT, IFN-gamma |
| Gastrocnemius | OXIDATIVE_PHOSPHORYLATION, MYOGENESIS |
| Eye | OXIDATIVE_PHOSPHORYLATION (dominant 3/3 missions) |
| Skin | E2F_TARGETS, EPITHELIAL_MESENCHYMAL_TRANSITION |
| Liver | OXIDATIVE_PHOSPHORYLATION, FATTY_ACID_METABOLISM |
| Kidney | MTORC1_SIGNALING, CHOLESTEROL_HOMEOSTASIS |

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
git clone https://github.com/jang1563/GeneLab_benchmark
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
        "fold_MHU-2_test": {"sample_id_1": 0.92, "...": "..."},
        "fold_RR-6_test":  {"...": "..."},
        "fold_RR-7_test":  {"...": "..."}
    }
}

# Evaluate (requires cloning GitHub repo)
# python scripts/evaluate_submission.py --submission my_submission.json --task A5
```

**GO criteria**: Mean AUROC > 0.70 AND 95% CI lower > 0.50 AND perm_p < 0.05.

### Model Tracks

| Track | Examples | Input Format |
|-------|---------|-------------|
| **Tier 1 — Classical ML** | LR, RF, XGBoost, PCA-LR | Tabular gene × sample |
| **Tier 2 — Foundation Models** | Geneformer, scGPT | Gene rank order (tokenized) |
| **Tier 3 — Text LLMs** | GPT-4o, Claude, Llama 3 | Natural language gene list |

---

## Source Data

All data derived from publicly available NASA OSDR datasets:

| OSD ID | Tissue | Mission | n (binary) |
|--------|--------|---------|-----------|
| OSD-48 | Liver | RR-1 | 18 |
| OSD-137 | Liver | RR-3 | 20 |
| OSD-245 | Liver | RR-6 | 48 |
| OSD-379 | Liver | RR-8 | 40 |
| OSD-242 | Liver | RR-9 | 39 |
| OSD-686 | Liver | MHU-2 | 28 |
| OSD-101 | Gastrocnemius | RR-1 | 12 |
| OSD-401 | Gastrocnemius | RR-5 | 12 |
| OSD-326 | Gastrocnemius | RR-9 | 8 |
| OSD-102 | Kidney | RR-1 | 47 |
| OSD-163 | Kidney | RR-3 | 32 |
| OSD-253 | Kidney | RR-7 | 39 |
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
  title   = {GeneLab Benchmark: A Multi-Tissue Spaceflight Transcriptomics Benchmark for AI/ML Models},
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

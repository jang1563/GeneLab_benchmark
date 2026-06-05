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
  - mouse
  - rna-seq
  - foundation-model
  - pathway-analysis
  - single-cell
  - spatial-transcriptomics
size_categories:
  - 1GB<n<10GB
language:
  - en
pretty_name: "GeneLab Spaceflight Transcriptomics Benchmark"
---

# GeneLab Spaceflight Transcriptomics Benchmark

**A comprehensive benchmark for evaluating ML models and foundation models on NASA spaceflight transcriptomics data.**

Version: v7.0 with v7.1 documentation consistency patch | Dataset freeze: 2026-03-01 | Code: [GitHub](https://github.com/jang1563/GeneLab_benchmark)

---

## Overview

GeneLab Benchmark provides standardized train/test splits for evaluating how well machine learning models generalize **spaceflight transcriptomic signatures across ISS missions**.

**Core challenge**: Given RNA-seq from one spaceflight mission, can a classifier detect spaceflight vs. ground control in samples from a different mission it has never seen?

**Data source**: NASA Open Science Data Repository ([OSDR](https://osdr.nasa.gov/bio/repo/)) -- mouse multi-tissue bulk RNA-seq from ISS rodent research missions (C57BL/6J strain).

### Scope

| Dimension | Coverage |
|---|---|
| Tissues | 8 (Liver, Gastrocnemius, Kidney, Thymus, Skin, Eye, Lung, Colon) |
| Core ISS mission labels | 9 named labels in the public task package (RR-1, RR-3, RR-5, RR-6, RR-7, RR-8, RR-9, MHU-1, MHU-2) |
| OSD Studies | 24+ source accessions across release layers |
| Samples | 600+ binary/control samples across the processed release layers |
| Classifiers | 8 (PCA-LR, ElasticNet-LR, RF, XGBoost, SVM-Linear, SVM-RBF, TabNet, LightGBM) |
| Feature types | 4 (Gene, Hallmark pathways, KEGG pathways, Combined) |
| Foundation / Language Models | 4 gene-expression FMs + 3 text LLMs |

---

The GitHub repository contains the full v1-v7 benchmark surface. This Hugging Face dataset card exposes the public feature-matrix package and reviewer-facing result summary; counts below separate full-release scope from specific analysis subsets.

## Transparency and Release Boundary

This Hugging Face card is the concise public card for the v1-v7 feature-matrix
package. Project-level transparency documentation is maintained in GitHub so the
system scope, evaluation protocol, release-readiness gates, and claim boundaries
can be reviewed together:

- Transparency card pack:
  https://github.com/jang1563/GeneLab_benchmark/blob/main/docs/SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md
- System card:
  https://github.com/jang1563/GeneLab_benchmark/blob/main/docs/SPACEBIOBENCH_SYSTEM_CARD.md
- Evaluation card:
  https://github.com/jang1563/GeneLab_benchmark/blob/main/docs/SPACEBIOBENCH_EVALUATION_CARD.md
- Release readiness card:
  https://github.com/jang1563/GeneLab_benchmark/blob/main/docs/SPACEBIOBENCH_RELEASE_READINESS_CARD.md
- Claim register:
  https://github.com/jang1563/GeneLab_benchmark/blob/main/docs/SPACEBIOBENCH_CLAIM_REGISTER.md

For DOI-oriented citation, use the final GitHub release tag and archive DOI
recorded in the release archive manifest once minted.

Current release wording: v1-v7 are the canonical historical benchmark and
reviewer-facing feature-matrix surface; v8 and v9 are documented as extension or
metadata-alpha surfaces and should not be described as frozen public benchmark
payload releases.

## Dataset Structure

This HuggingFace repository contains **feature matrices** (train_X.csv, test_X.csv) for 4 benchmark tasks that passed significance thresholds. Labels and metadata are in the [GitHub repository](https://github.com/jang1563/GeneLab_benchmark).

```
genelab-benchmark/
├── A2_gastrocnemius_lomo/        <- 3 missions, 32 samples
│   ├── fold_RR-1_test/
│   │   ├── train_X.csv           <- Training features (samples x genes)
│   │   └── test_X.csv            <- Test features
│   ├── fold_RR-5_test/
│   └── fold_RR-9_test/
│
├── A4_thymus_lomo/               <- 4 missions, 67 samples
│   ├── fold_MHU-1_test/
│   ├── fold_MHU-2_test/
│   ├── fold_RR-6_test/
│   └── fold_RR-9_test/
│
├── A5_skin_lomo/                 <- 3 missions, 102 samples
│   ├── fold_MHU-2_test/
│   ├── fold_RR-6_test/
│   └── fold_RR-7_test/
│
├── A6_eye_lomo/                  <- 3 missions, 37 samples
│   ├── fold_RR-1_test/
│   ├── fold_RR-3_test/
│   └── fold_OSD-397_test/        <- stable public label for OSD-397
│
├── v4/evaluation/                <- Multi-method evaluation results (JSON)
├── v5/evaluation/                <- Systems biology analysis results
└── v6/evaluation/                <- Human translation analysis results
```

Each fold holds out one mission as the test set and trains on the remaining missions. This **Leave-One-Mission-Out (LOMO)** cross-validation evaluates true cross-mission generalization.

---

## File Format

### Feature matrix (train_X.csv, test_X.csv)

- **Rows**: Sample IDs
- **Columns**: Ensembl mouse gene IDs (e.g., `ENSMUSG00000021969`)
- **Values**: Log2(DESeq2 size-factor normalized counts + 1)
- **Gene selection**: Top 75th percentile variance, computed on **training missions only** (no test leakage)
- **Typical shape**: ~20,000 genes per sample

### Labels (in GitHub repo)

| Value | Meaning |
|-------|---------|
| `1` | Flight (spaceflight / microgravity) |
| `0` | Ground Control (vivarium / ground control) |

---

## Results

### Multi-Method Evaluation (256 evaluations: 8 tissues x 8 methods x 4 feature types)

Best AUROC per tissue across all method-feature combinations, evaluated via LOMO cross-validation (6 tissues) or 5-fold stratified CV (Lung, Colon -- single-mission datasets):

| Tissue | Best AUROC | Method | Feature | perm_p | Significant |
|--------|-----------|--------|---------|--------|-------------|
| **Thymus** | **0.948** | PCA-LR | KEGG | 0.031 | Yes |
| **Colon** | **0.921** | PCA-LR | KEGG | 0.033 | Yes |
| **Lung** | **0.901** | PCA-LR | Gene | <0.05 | Yes |
| **Kidney** | **0.829** | ElasticNet-LR | Hallmark | <0.01 | Yes |
| **Eye** | **0.823** | PCA-LR | Hallmark | -- | - |
| **Skin** | **0.819** | PCA-LR | Gene | -- | - |
| **Gastrocnemius** | **0.776** | PCA-LR | Gene | -- | - |
| **Liver** | **0.670** | PCA-LR | Gene | -- | - |

Significance: best-row permutation p where available. Across all 256 configurations, 40 individual evaluations are significant at p<0.05 and 6/8 tissues have at least one significant configuration.

PCA-LR on gene features provides a strong baseline (8-tissue mean AUROC = 0.776).

### Foundation / Language Model Snapshot

| Model | Best Single-Tissue AUROC | vs PCA-LR Baseline (0.776) |
|-------|-----------|-------------------|
| scGPT | 0.666 (6-tissue mean, v1) | Below baseline |
| scFoundation | 0.635 (liver, p<0.01) | Below baseline |
| UCE | 0.632 (thymus, p=0.031) | Below baseline |
| Mouse-Geneformer | 0.476 (6-tissue mean, v1) | Below baseline |
| Text LLMs (GPT-4o, Claude, Llama 3) | 0.47-0.51 | Chance level |

scGPT and Mouse-Geneformer report 6-tissue v1 means; scFoundation and UCE rows show best single-tissue v3 results. All of these model families still underperform the classical PCA-LR benchmark surface.

### Negative Controls

- Permutation control: AUROC = 0.50 +/- 0.03 (label shuffling)
- Housekeeping gene control: AUROC = 0.49-0.55 (non-informative features)
- Held-out validation: Thymus RR-23 (0.905), Skin RR-7 (0.885)

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

---

## Evaluation

Evaluate predictions using the included script ([GitHub](https://github.com/jang1563/GeneLab_benchmark)):

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

# Run evaluation (requires cloning GitHub repo for labels)
# python scripts/evaluate_submission.py --submission my_submission.json --task A5
```

### Model Tracks

| Track | Examples | Input Format |
|-------|---------|-------------|
| **Classical ML** | PCA-LR, ElasticNet-LR, RF, XGBoost, SVM-Linear, SVM-RBF, TabNet, LightGBM | Tabular (gene x sample) |
| **Foundation Models** | Geneformer, scGPT, UCE, scFoundation | Gene rank order / embeddings |
| **Text LLMs** | GPT-4o, Claude, Llama 3 | Natural language gene list |

---

## Source Data

All data derived from publicly available NASA OSDR datasets:

| OSD ID | Tissue | Mission | n (Flight + Ground) |
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
| OSD-289 | Thymus | MHU-2 | 12 |
| OSD-244 | Thymus | RR-6 | 35 |
| OSD-421 | Thymus | RR-9 | 20 |
| OSD-238 | Skin (dorsal) | MHU-2 | 18 |
| OSD-239 | Skin (femoral) | MHU-2 | 17 |
| OSD-243 | Skin | RR-6 | 37 |
| OSD-254 | Skin | RR-7 | 30 |
| OSD-100 | Eye | RR-1 | 12 |
| OSD-194 | Eye | RR-3 | 9 |
| OSD-397 | Eye | OSD-397 | 16 |
| OSD-248 | Lung | RR-6 | 39 |
| OSD-247 | Colon | RR-6 | 36 |

Lung and Colon additionally include Basal Control samples treated as ground control (+19 and +18 respectively).

---

## Preprocessing

1. DESeq2 size-factor normalization (per mission)
2. Log2(counts + 1) transformation
3. Low-expression filter (>=20% samples with count > 1)
4. Top 75th percentile variance gene selection **per fold, training missions only** (prevents test leakage)
5. Pathway scores: gseapy ssGSEA (MSigDB Hallmark, KEGG gene sets)

---

## Citation

*(Manuscript in preparation)*

```bibtex
@dataset{kim2026genelab,
  title   = {GeneLab Benchmark: A Multi-Tissue Spaceflight Transcriptomics
             Benchmark for AI/ML Models},
  author  = {Kim, Jihoon},
  year    = {2026},
  url     = {https://huggingface.co/datasets/jang1563/genelab-benchmark},
  note    = {v7.0 with v7.1.1 documentation, public-card, and metadata consistency patch; data freeze 2026-03-01}
}
```

Data source: NASA Open Science Data Repository (OSDR) -- https://osdr.nasa.gov/bio/repo/

---

## License

- Dataset: CC-BY-4.0
- Code: MIT ([GitHub repository](https://github.com/jang1563/GeneLab_benchmark))
- Source data: NASA OSDR public data; see individual OSDR dataset licenses.

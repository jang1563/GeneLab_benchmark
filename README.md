# GeneLab Benchmark

**A public benchmark for evaluating AI/ML and Foundation Models on NASA OSDR spaceflight transcriptomics data.**

Version: v1.0-alpha (Dataset freeze: 2026-03-01)
Status: Phase 1 complete — **4 tissues GO** (A2 Gastrocnemius, A4 Thymus, A5 Skin, A6 Eye pathway-level)

---

## What Is This?

GeneLab Benchmark provides standardized tasks for evaluating how well machine learning models — from classical baselines to gene expression foundation models (Geneformer, scGPT) and text-based LLMs (GPT-4o, Claude) — generalize **spaceflight transcriptomic signatures across missions**.

**Core challenge**: Train a model on one spaceflight mission's RNA-seq data. Can it classify samples from a different mission it has never seen?

**Data source**: [NASA Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/) — mouse multi-tissue bulk RNA-seq from ISS and ground control missions.

---

## Key Features

- **Leave-One-Mission-Out (LOMO)** cross-validation — mission = independence unit, preventing cross-mission data leakage
- **Category A**: Spaceflight detection per tissue (binary: Flight vs. Ground)
- **Category B**: Cross-mission transfer matrix (train on mission i, test on mission j) for all 6 tissues
- **3-tier model evaluation**: Classical ML → Gene Expression Foundation Models → Text LLMs
- **Standardized submission format** with automatic AUROC/CI/p-value evaluation

---

## Phase 1 Results Summary

**Category A — Gene-level (primary)**:

| Task | Tissue | Missions | Method | Mean AUROC | 95% CI lower | perm_p | Decision |
|------|--------|---------|--------|-----------|-------------|--------|----------|
| A4 | Thymus | 4† | PCA-LR | **0.923** | 0.878 | 0.037 | ✓ **GO** |
| A2 | Gastrocnemius | 3 | LR | **0.907** | 0.717 | 0.026 | ✓ **GO** |
| **A5** | **Skin** | **3**§ | **LR** | **0.821** | **0.637** | **0.0023** | ✓ **GO** |
| A6 | Eye | 3 | LR | 0.811 | 0.470 | 0.063 | ✗ NO-GO‡ |
| A1 | Liver | 6 | LR | 0.653 | 0.457 | 0.091 | ✗ NO-GO |
| A3 | Kidney | 3 | LR | 0.593 | 0.431 | 0.281 | ✗ NO-GO |

**Category A — Pathway-level (GSVA Hallmark, secondary)**:

| Task | Tissue | Method | Mean AUROC | 95% CI lower | perm_p | Decision |
|------|--------|--------|-----------|-------------|--------|----------|
| **A6** | **Eye** | **PCA-LR** | **0.915** | **0.745** | **0.014** | ✓ **GO**⊕ |
| A3 | Kidney | LR | 0.755 | 0.481 | 0.071 | ✗ NO-GO |

†A4 includes MHU-1 (Track 2b, GC/FLT strain mismatch — see PHASE1_RESULTS.md)
§A5: MHU-2 = dorsal (OSD-238) + femoral (OSD-239) merged; RR-7 = OSD-254 C57BL/6J non-BSL subset (n=30)
‡A6 gene-level: AUROC passes but CI lower fails (n=9–16 per fold)
⊕A6 pathway-level: GSVA Hallmark 50-pathway scores rescue CI lower (0.470→0.745). Oxidative phosphorylation dominant.

**Category B — Cross-Mission Transfer (PCA-LR)**:

| Task | Tissue | N pairs | Mean AUROC | 95% CI | AUROC≥0.70 |
|------|--------|---------|-----------|--------|-----------|
| B4 | Thymus | 12 | **0.860** | [0.763, 0.953] | 9/12 |
| B5 | Skin | 6 | **0.772** | [0.691, 0.834] | 5/6 |
| B6 | Eye | 6 | 0.754 | [0.688, 0.838] | 5/6 |
| B2 | Gastrocnemius | 6 | 0.801 | [0.653, 0.944] | 4/6 |
| B1 | Liver | 30 | 0.577 | [0.492, 0.666] | 13/30 |
| B3 | Kidney | 6 | 0.555 | [0.397, 0.681] | 2/6 |

See [PHASE1_RESULTS.md](docs/development_history/PHASE1_RESULTS.md) for full results including per-fold tables, SHAP analysis, within-LOO sanity checks, and pathway analysis.

---

## Repository Structure

```
GeneLab_benchmark/
├── README.md                   ← This file
├── PLAN.md                     ← Benchmark design specification (v0.7)
├── DESIGN_DECISIONS.md         ← Architecture decisions log (DD-01 to DD-17)
├── docs/development_history/PHASE1_RESULTS.md ← Full Phase 1 analysis results
│
├── tasks/                      ← Public task inputs (features + labels)
│   ├── A2_gastrocnemius_lomo/
│   │   ├── task_info.json
│   │   ├── fold_RR-1_test/
│   │   │   ├── train_X.csv     ← Training features (samples × genes)
│   │   │   ├── train_y.csv     ← Training labels (1=Flight, 0=Ground)
│   │   │   ├── test_X.csv      ← Test features (PUBLIC)
│   │   │   └── test_y.csv      ← Test labels (PUBLIC for LOMO)
│   │   └── ...
│   ├── A4_thymus_lomo/
│   │   └── ...
│   └── A5_skin_lomo/           ← NEW
│       ├── fold_MHU-2_test/
│       ├── fold_RR-6_test/
│       └── fold_RR-7_test/
│
├── scripts/
│   ├── run_baselines.py        ← Classical ML baseline runner
│   ├── generate_tasks.py       ← Task split generator
│   ├── cross_mission_transfer.py ← Category B matrix generator
│   ├── shap_analysis.py        ← SHAP feature importance
│   └── evaluate_submission.py  ← Submission evaluator
│
├── docs/
│   └── submission_format.md    ← Submission format specification
│
├── evaluation/                 ← Baseline evaluation results
│   ├── A4_baseline_results.json
│   ├── A2_baseline_results.json
│   ├── A5_baseline_results.json    ← NEW
│   ├── A5_shap_rf.json             ← NEW
│   └── B_cross_mission_summary.json
│
└── processed/                  ← Processed data (intermediate)
    ├── A_detection/
    │   ├── gastrocnemius/
    │   ├── thymus/
    │   └── skin/               ← NEW
    └── B_cross_mission/
        ├── liver/ gastrocnemius/ kidney/ thymus/ eye/
        └── skin/               ← NEW
```

---

## Getting Started

Feature matrices (train_X.csv, test_X.csv) are hosted on HuggingFace due to size (~2 GB).
Labels, metadata, and fold structure are in this repository.

### Option A — Load from HuggingFace (recommended)

```bash
pip install -r requirements.txt huggingface_hub
```

```python
from huggingface_hub import hf_hub_download
import pandas as pd

train_X = pd.read_csv(
    hf_hub_download(
        repo_id="jang1563/genelab-benchmark",
        filename="A5_skin_lomo/fold_RR-7_test/train_X.csv",
        repo_type="dataset",
    ),
    index_col=0,
)
train_y = pd.read_csv("tasks/A5_skin_lomo/fold_RR-7_test/train_y.csv", index_col=0)
print(f"Train: {train_X.shape}")  # (72, 20110)
```

Or download a full task at once:

```bash
python scripts/upload_to_hf.py --task A5 --dry-run   # preview
# After cloning HF data locally, run baselines:
python scripts/run_baselines.py --task A5 --model lr
```

### Option B — Reproduce from OSDR raw data

Requires R 4.2+ with Bioconductor. See [docs/r_dependencies.md](docs/r_dependencies.md).

```bash
# 1. Download raw data from NASA OSDR
python scripts/fetch_osdr.py --osd OSD-238 OSD-239 OSD-243 OSD-254

# 2. Normalize (DESeq2)
Rscript scripts/normalize_rr7_skin.R   # example

# 3. Quality filter + build all_missions
python scripts/quality_filter.py --tissue skin

# 4. Generate LOMO folds
python scripts/generate_tasks.py --task A5
```

---

## Quick Start

### 1. Explore a task

```python
import pandas as pd

# Load A5 Skin — fold RR-7 test
train_X = pd.read_csv("tasks/A5_skin_lomo/fold_RR-7_test/train_X.csv", index_col=0)
train_y = pd.read_csv("tasks/A5_skin_lomo/fold_RR-7_test/train_y.csv", index_col=0)
test_X  = pd.read_csv("tasks/A5_skin_lomo/fold_RR-7_test/test_X.csv", index_col=0)
test_y  = pd.read_csv("tasks/A5_skin_lomo/fold_RR-7_test/test_y.csv", index_col=0)

print(f"Train: {train_X.shape}, Test: {test_X.shape}")
print(f"Train labels: {train_y.iloc[:,0].value_counts().to_dict()}")
# Features: Ensembl mouse gene IDs (e.g., ENSMUSG00000021969)
# Labels: 1.0 = Flight, 0.0 = Ground/Vivarium Control
```

### 2. Run a baseline model

```bash
python scripts/run_baselines.py --task A5 --model lr
python scripts/run_baselines.py --task A4 --model pca_lr
```

### 3. Submit your model's predictions

Prepare a JSON file (see [docs/submission_format.md](docs/submission_format.md)):

```json
{
  "task_id": "A5",
  "model_name": "MyModel_v1",
  "predictions": {
    "fold_MHU-2_test": {"sample_id_1": 0.92, "sample_id_2": 0.07, ...},
    "fold_RR-6_test":  {...},
    "fold_RR-7_test":  {...}
  }
}
```

Evaluate:

```bash
python scripts/evaluate_submission.py \
    --submission my_submission.json \
    --task A5
```

---

## Tasks (v1.0)

### Category A — Spaceflight Detection (LOMO)

**Goal**: Binary classification (Flight vs. Ground) using Leave-One-Mission-Out CV.

| Task | Tissue | Missions | Samples (binary) | Folds | Status |
|------|--------|---------|-----------------|-------|--------|
| A2 | Gastrocnemius | RR-1, RR-5, RR-9 | 32 | 3 | ✓ GO |
| A4 | Thymus | MHU-1†, MHU-2, RR-6, RR-9 | 67 | 4 | ✓ GO |
| A5 | Skin | MHU-2§, RR-6, RR-7 | 102 | 3 | ✓ GO |
| A6 | Eye | RR-1, RR-3, TBD | 37 | 3 | ✓ GO (pathway) |
| A1 | Liver | MHU-2, RR-1, RR-3, RR-6, RR-8, RR-9 | 193 | 6 | ✗ NO-GO |
| A3 | Kidney | RR-1, RR-3, RR-7 | 118 | 3 | ✗ NO-GO |

†MHU-1 = Track 2b (GC strain = C57BL/6CR, FLT = C57BL/6J mismatch — see [PHASE1_RESULTS.md](docs/development_history/PHASE1_RESULTS.md))
§MHU-2 = OSD-238 (dorsal) + OSD-239 (femoral) merged as single mission; RR-7 = OSD-254 C57BL/6J non-BSL subset

**Input**: Log2-normalized expression values for ~20,000 mouse genes (Ensembl IDs, e.g., ENSMUSG00000021969).
**Label**: `1.0` = Flight, `0.0` = Ground/Vivarium Control. Basal Control (BC) samples excluded.

### Category B — Cross-Mission Transfer

**Goal**: Train on one mission, evaluate generalization to another (all N×(N-1) ordered pairs).

See `processed/B_cross_mission/{tissue}/` for per-tissue AUROC matrices and `evaluation/B_cross_mission_summary.json` for aggregated results.

---

## Baseline Submissions

Pre-computed baseline predictions are available in `evaluation/` for reference and reproducibility.

**Category A (LOMO)**

| File | Task | Model | Mean AUROC | Go/No-Go |
|------|------|-------|------------|----------|
| `submission_LR_baseline_A5.json` | A5 Skin | LR (ElasticNet) | 0.821 | ✓ GO |
| `submission_PCALR_baseline_A4.json` | A4 Thymus | PCA-LR (L2, lbfgs) | 0.923 | ✓ GO |
| `submission_LR_baseline_A2.json` | A2 Gastrocnemius | LR-ElasticNet (SAGA) | 0.917 | ✓ GO |

**Category B (Cross-Mission Transfer)**

| Task | Tissue | N pairs | PCA-LR Mean AUROC | LFC Mean AUROC |
|------|--------|---------|------------------|---------------|
| B4 | Thymus | 12 | 0.860 | 0.868 |
| B5 | Skin | 6 | 0.772 | 0.750 |
| B6 | Eye | 6 | 0.754 | 0.696 |
| B2 | Gastrocnemius | 6 | 0.801 | 0.655 |
| B1 | Liver | 30 | 0.577 | 0.534 |
| B3 | Kidney | 6 | 0.555 | 0.465 |

> Category B does not report a single GO/NO-GO — see [DD-17](DESIGN_DECISIONS.md) for evaluation criteria.

Evaluate a baseline submission:

```bash
# Category A
python scripts/evaluate_submission.py \
    --submission evaluation/submission_LR_baseline_A5.json \
    --task A5

# Category B (summary across all tissues)
python scripts/cross_mission_transfer.py --tissue skin
```

> **Reproducibility note**: The official `A2_baseline_results.json` was computed with `max_iter=2000` (SAGA not fully converged for 15k genes). The baseline submission above uses `max_iter=10000` (converged); A2 mean AUROC improves from 0.907 → 0.917. GO/No-Go conclusion unchanged. See `PHASE1_RESULTS.md §B3` for details.

---

## Evaluation Protocol

All submissions are evaluated with:

| Metric | Description | Go threshold |
|--------|-------------|-------------|
| Mean AUROC | Average AUROC across folds | > 0.700 |
| 95% CI lower | Bootstrap CI (N=2000) lower bound | > 0.500 |
| perm_p | Permutation p-value (N=1000, pseudocount) | < 0.050 |

All three conditions must pass for a GO decision.

---

## Model Tracks

| Track | Examples | Input Format |
|-------|---------|-------------|
| **Tier 1 — Classical ML** | LR, RF, XGBoost, PCA-LR | Tabular gene × sample |
| **Tier 2 — Foundation Models** | Geneformer, scGPT | Gene rank order (tokenized) |
| **Tier 3 — Text LLMs** | GPT-4o, Claude, Llama 3 | Natural language gene list (see DD-16) |

For Tier 3 (Text LLM) input format specification, see [DESIGN_DECISIONS.md](DESIGN_DECISIONS.md) (DD-16).

---

## Data

All data is derived from publicly available NASA OSDR datasets.

| Tissue | OSD Accession | Mission | n samples | Note |
|--------|--------------|---------|-----------|------|
| Thymus | OSD-289 | MHU-1 | 6 | Track 2b (GC = C57BL/6CR) |
| Thymus | OSD-289 | MHU-2 | 6 | Track 2a |
| Thymus | OSD-244 | RR-6 | 35 | Track 2a |
| Thymus | OSD-421 | RR-9 | 20 | Track 2a |
| Gastrocnemius | OSD-101 | RR-1 | 12 | Track 2a |
| Gastrocnemius | OSD-401 | RR-5 | 12 | Track 2a |
| Gastrocnemius | OSD-326 | RR-9 | 8 | Track 2a |
| Skin | OSD-238 | MHU-2 (dorsal) | 18 | merged as "MHU-2" (6F+6GC+6VC; AG excluded) |
| Skin | OSD-239 | MHU-2 (femoral) | 17 | merged as "MHU-2" (5F+12GC; AG excluded) |
| Skin | OSD-243 | RR-6 | 37 | Track 2a |
| Skin | OSD-254 | RR-7 | 30 | C57BL/6J non-BSL subset only |

MHU-1 and MHU-2 are both sub-experiments within OSD-289 (GLDS-289), separated by mission label during preprocessing.
OSD-254 (RR-7 Skin) is a mixed-strain study (C57BL/6J + C3H/HeJ); only the C57BL/6J non-BSL samples (n=30) are included in Track 2a (A5).

Preprocessing: DESeq2 normalization (per-mission), log2(counts + 1), global low-expression filter (≥20% samples with count>1), top 75th percentile variance gene selection per fold (train missions only — DD-03).

---

## Design Decisions

Key methodological choices are documented in [DESIGN_DECISIONS.md](DESIGN_DECISIONS.md):

- **DD-03**: LOMO-aware variance filter (train missions only — no test leakage)
- **DD-04**: Mission = independence unit for LOMO (not sample)
- **DD-06**: Track 2a = C57BL/6J only; Track 2b = all strains
- **DD-08**: Evaluation metrics (AUROC + bootstrap CI + permutation p)
- **DD-11**: Go/No-Go decision criteria
- **DD-13**: Baseline model set (LR, RF, XGBoost, PCA-LR)
- **DD-16**: Text LLM evaluation track specification
- **DD-17**: Category B evaluation criteria (Transfer Pattern Summary, perm_p floor)

---

## Changelog

| Version | Date | Changes |
|---------|------|---------|
| v1.0-alpha | 2026-03-01 | Phase 1 complete. A2+A4+A5 GO (gene-level). A6 Eye pathway-level GO (GSVA Hallmark). Category B all 6 tissues (B1–B6). Submission format + evaluator. Dataset freeze. |

---

## Citation

*(Manuscript in preparation)*

Data source: NASA Open Science Data Repository (OSDR) — [osdr.nasa.gov](https://osdr.nasa.gov/bio/repo/)

---

## License

Code: MIT License
Data: NASA OSDR public data (see individual dataset licenses at OSDR)

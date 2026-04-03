# GeneLab Benchmark

**A comprehensive benchmark for evaluating AI/ML and Foundation Models on NASA OSDR spaceflight transcriptomics data.**

Version: v6.0 (2026-03-30) | Dataset freeze: 2026-03-01
Status: **v1–v6 Complete** | v7 Graph Neural Networks In Progress

[![Dataset on HuggingFace](https://img.shields.io/badge/HuggingFace-Dataset-yellow)](https://huggingface.co/datasets/jang1563/genelab-benchmark)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## What Is This?

GeneLab Benchmark provides standardized tasks for evaluating how well machine learning models — from classical baselines to gene expression foundation models (Geneformer, scGPT, UCE, scFoundation) and text-based LLMs (GPT-4o, Claude, Llama) — generalize **spaceflight transcriptomic signatures across missions**.

**Core challenge**: Train a model on one spaceflight mission's RNA-seq data. Can it classify samples from a different mission it has never seen?

**Data source**: [NASA Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/) — mouse multi-tissue bulk RNA-seq from ISS and ground control missions.

### Benchmark Scope

- **8 tissues**: Liver, Gastrocnemius, Kidney, Thymus, Skin, Eye, Lung, Colon
- **17+ ISS missions**: RR-1 through RR-9, MHU-1, MHU-2, and more
- **24+ verified OSD studies**, ~600+ samples (binary Flight/Ground, with BC/VC controls in v4)
- **25+ evaluation tasks** across 7 categories (A-D, J, NC, Validation)
- **v4**: 256 evaluations (8 tissues x 8 classifiers x 4 feature types)
- **Multi-species**: Mouse bulk, Drosophila KEGG, mouse scRNA-seq (RRRM-1/RRRM-2), spatial Visium
- **5 foundation models evaluated**: Geneformer, scGPT, UCE, scFoundation, + 3 text LLMs

---

## Key Features

- **Leave-One-Mission-Out (LOMO)** cross-validation — mission = independence unit, preventing cross-mission data leakage
- **Category A**: Spaceflight detection per tissue (binary: Flight vs. Ground)
- **Category B**: Cross-mission transfer matrix (train on mission i, test on mission j) for all tissues
- **Category C**: Cross-tissue transfer (3 methods: gene, DEG, pathway)
- **Category D**: Condition/confounder prediction (mission, strain, hardware, gravity)
- **3-tier model evaluation**: Classical ML -> Gene Expression Foundation Models -> Text LLMs
- **Standardized submission format** with automatic AUROC/CI/p-value evaluation
- **Biological validation**: NES pathway conservation, Cell 2020 concordance, negative controls
- **v4 multi-method**: 8 classifiers (PCA-LR, ElasticNet-LR, RF, XGBoost, SVM-Linear, SVM-RBF, TabNet, LightGBM)
- **v4 multi-feature**: 4 feature types (gene, Hallmark, KEGG, pathway-combined)

---

## Results Overview

### v4 Phase 1 — Multi-Method Evaluation (256 evaluations)

**Best overall classifier**: PCA-LR (gene mean AUROC 0.776), ElasticNet-LR 2nd (0.762)

| Tissue | Best AUROC | Best Method | Best Feature | Significant (p<0.05) |
|--------|-----------|-------------|-------------|---------------------|
| **Thymus** | **0.948*** | PCA-LR | KEGG | Yes |
| **Colon** | **0.921*** | PCA-LR | KEGG | Yes |
| **Lung** | **0.901*** | PCA-LR | Gene | Yes |
| **Kidney** | **0.829**** | ElasticNet-LR | Hallmark | Yes |
| **Skin** | 0.819 | PCA-LR | Gene | - |
| **Eye** | 0.823 | PCA-LR | Hallmark | - |
| **Liver** | 0.670 | PCA-LR | Gene | - |
| **Gastrocnemius** | 0.776 | PCA-LR | Gene | - |

*p<0.05, **p<0.01. 40/256 evaluations significant at p<0.05; 6/8 tissues have at least 1 significant result.

**Key v4 findings**:
- SVM-RBF worst on high-dimensional gene data (0.510 mean)
- TabNet disappointing (0.527 mean) — deep learning doesn't help small-n transcriptomics
- Pathway features improve kidney (0.584->0.829), thymus (0.908->0.948), eye (0.697->0.823)
- Gene features better for skin (0.819) and lung (0.901)
- v4 expanded controls: BC/VC included (e.g., liver 261 vs v1's 193 samples)

### v1 Results — Spaceflight Detection (LOMO, 6 tissues)

**Gene-level (primary)**:

| Task | Tissue | Missions | Method | Mean AUROC | 95% CI lower | perm_p | Decision |
|------|--------|---------|--------|-----------|-------------|--------|----------|
| A4 | Thymus | 4 | PCA-LR | **0.923** | 0.878 | 0.037 | GO |
| A2 | Gastrocnemius | 3 | LR | **0.907** | 0.717 | 0.026 | GO |
| A5 | Skin | 3 | LR | **0.821** | 0.637 | 0.0023 | GO |
| A6 | Eye | 3 | LR | 0.811 | 0.470 | 0.063 | NO-GO |
| A1 | Liver | 6 | LR | 0.653 | 0.457 | 0.091 | NO-GO |
| A3 | Kidney | 3 | LR | 0.593 | 0.431 | 0.281 | NO-GO |

**Pathway-level rescue**: Eye gene-level NO-GO -> GSVA Hallmark AUROC **0.915** (GO). Kidney 0.43->0.74.

### Cross-Mission Transfer (Category B)

| Tissue | Mean AUROC | 95% CI | AUROC>=0.70 | Tier |
|--------|-----------|--------|-----------|------|
| Thymus | **0.860** | [0.763, 0.953] | 9/12 | 1 |
| Gastrocnemius | **0.801** | [0.653, 0.944] | 4/6 | 1 |
| Skin | **0.772** | [0.691, 0.834] | 5/6 | 2 |
| Eye | 0.754 | [0.688, 0.838] | 5/6 | 2 |
| Liver | 0.577 | [0.492, 0.666] | 13/30 | 3 |
| Kidney | 0.555 | [0.397, 0.681] | 2/6 | 3 |

### v5 — Biological Interpretation

| Analysis | Key Finding |
|----------|------------|
| **Immune deconvolution** (mMCP-counter, 8 tissues) | Skin strongest: 6/14 cell types FDR<0.05. Kidney, thymus 2/14 each. Most tissues 0 significant |
| **Cross-organ signaling** (OmniPath, 111 L–R pairs) | 1 SHAP-active L–R pair. TF activity: thymus 240 sig, skin 241, kidney 177, liver 105 |
| **Metabolic flux** (iMM1865 E-Flux + pFBA) | Flight vs ground objectives differ in all 6 tissues (largest: thymus FLT 15,695 vs GC 14,696) |
| **Drug targets** (DGIdb + ChEMBL) | 271/834 spaceflight genes druggable; 1,284 FDA-approved drug–gene interactions |
| **Consensus biomarker panel** (20 genes) | Top genes: MUP22, Thrsp, Apoa1, NPAS2, PER2. AUROC: gastro 0.806, liver 0.754, eye 0.728 |

### Foundation Model Comparison (5 FMs vs PCA-LR)

| Model | Architecture | Mean AUROC (6 tissues) | vs PCA-LR |
|-------|------------|----------------------|-----------|
| **PCA-LR baseline** | Classical ML | **0.758** | — |
| scGPT | 12L Transformer, 33M human cells | 0.667 | -0.092 |
| scFoundation | 100M params, 50M cells | 0.635† (liver best, p<0.01) | Below baseline |
| UCE | 33L, 36M cells | 0.632† (thymus best, p=0.031) | Below baseline |
| Mouse-Geneformer | 6L BERT, 30M mouse cells | 0.476 | -0.283 |
| Text LLMs (3x) | GPT-4o, Claude, Llama 3 | 0.47-0.51 | Chance level |

**FM verdict**: All foundation models underperform classical PCA-LR (mean 0.758). Pre-trained cell atlas knowledge does not improve spaceflight detection. scFoundation > UCE overall. †Best single-tissue AUROC shown; mean across all 7 tissues is lower.

### Independent Held-Out Validation

| Tissue | Mission | AUROC | 95% CI | p-value | Duration | n_test |
|--------|---------|-------|--------|---------|---------|--------|
| Thymus | RR-23 | **0.905** | [0.672, 1.000] | 0.005 | 30 days | 18 |
| Skin | RR-7 | **0.885** | [0.745, 0.986] | <0.001 | 75 days | 30 |

---

## Key Scientific Findings

### Pre-registered Hypotheses

| Hypothesis | Statement | Verdict | Key Evidence |
|---|---|---|---|
| **H1** | Liver has the most consistent cross-mission transcriptome | **REFUTED** | Thymus (0.860) >> Liver (0.577). Thymus and Gastrocnemius = Tier 1. |
| **H2** | Transfer failure from biological diversity, not batch effects | **SUPPORTED** | NES conservation r=0.9 (5 tissues). D3 pathway F1=0.06 (batch-invariant). |
| **H3** | Pathway-level preserves spaceflight response better than gene-level | **CONDITIONALLY SUPPORTED** | Kidney rescue (0.43->0.74), Eye (0.79->0.92). But tissue-pair dependent. |

### Novel Findings

- **Tissue hierarchy**: Thymus > Gastrocnemius > Skin > Eye >> Liver > Kidney for spaceflight detection
- **Pathway rescue**: Pathway features rescue NO-GO tissues (kidney +0.31, eye +0.13)
- **FM failure mode**: scRNA-seq-pretrained FMs don't transfer to small-n bulk transcriptomics
- **Spaceflight amplifies aging**: OLD AUROC=0.945 vs YNG=0.679 (Delta=+0.266, RR-8 liver)
- **PBMC NK cells**: Strongest single-cell spaceflight signal (AUROC=0.845, RRRM-2)
- **Brain spatial**: Genuine negative — no detectable spaceflight signal (Visium AUROC=0.139)
- **Cross-species**: Drosophila-mouse NES correlations significant but negative (r=-0.19 to -0.59)

---

## Version Structure

| Version | Scope | Status | Location |
|---------|-------|--------|----------|
| **v1.0** | Mouse bulk RNA-seq, 6 tissues, 25+ tasks, 3 model tiers, 2 held-out validations | **Complete** | Project root |
| **v2.0** | Temporal dynamics (T1-T3), cross-species (E1-E3), single-cell (F1, F2 RRRM-1) | **Complete** | `v2/` |
| **v3.0** | Multi-species (E4), Spatial Visium (F3), RRRM-2 scRNA-seq (F5), UCE/scFoundation FM, radiation analogs, cross-tissue transfer 7×7 | **Complete** | `v3/` |
| **v4.0** | Multi-method benchmark: 8 tissues × 8 classifiers × 4 features (256 evals), ablation, SHAP, WGCNA, module preservation, STRING PPI | **Complete** | `v4/` |
| **v5.0** | Biological interpretation: immune deconvolution (mMCP-counter), metabolic flux (iMM1865 E-Flux), drug targets (DGIdb/ChEMBL), consensus 20-gene biomarker panel | **Complete** | `v5/` |
| **v6.0** | Cross-species validation: gene/pathway conservation, cross-species transfer, TF conservation, biomarker/drug target validation across species | **Complete** | `v6/` |

---

## Repository Structure

```
GeneLab_benchmark/
├── README.md                       <- This file
├── DATA_CATALOG.md                 <- OSDR inventory (24+ studies)
├── CITATION.cff                    <- Citation metadata
│
├── tasks/                          <- Public task inputs (17 directories)
│   ├── A1_liver_lomo/              <- 6 folds + 3 variants
│   ├── A2_gastrocnemius_lomo/      <- 3 folds
│   ├── A3_kidney_lomo/             <- 3 folds
│   ├── A4_thymus_lomo/             <- 4 folds + holdout
│   ├── A5_skin_lomo/               <- 3 folds
│   ├── A6_eye_lomo/                <- 3 folds
│   └── B1-B6_*_cross_mission/     <- N x (N-1) mission pairs per tissue
│
├── scripts/                        <- v1 pipeline scripts (35 Python/R/shell)
│   ├── run_baselines.py            <- Classical ML baseline runner
│   ├── evaluate_submission.py      <- Submission evaluator (AUROC, CI, perm_p)
│   ├── generate_tasks.py           <- LOMO split generator
│   └── ...                         <- See scripts/ for full list
│
├── docs/
│   ├── hf_dataset_card.md          <- HuggingFace dataset documentation
│   ├── submission_format.md        <- JSON submission specification
│   ├── text_llm_format.md          <- Text LLM evaluation format
│   └── BIOLOGICAL_GROUND_TRUTH.md  <- Validation reference (Cell 2020, SOMA 2024)
│
├── evaluation/                     <- Result JSON files (v1 baseline)
│   ├── RESULTS_SUMMARY.md          <- Comprehensive results table (v1–v5)
│   └── *.json                      <- Per-task results
│
├── v2/                             <- Temporal dynamics, cross-species, scRNA-seq
│   ├── README.md
│   ├── scripts/                    <- 19 Python scripts
│   ├── evaluation/                 <- v2 results
│   └── figures/                    <- 3 D3.js interactive HTML figures
│
├── v3/                             <- Multi-species, spatial Visium, FM evaluation
│   ├── README.md
│   ├── scripts/                    <- 30 scripts (Python/Bash/R)
│   ├── evaluation/                 <- 19 result JSONs
│   └── figures/                    <- 5 D3.js interactive HTML figures
│
├── v4/                             <- Multi-method benchmark + network biology
│   ├── scripts/                    <- 18 scripts (classifiers, SHAP, WGCNA, PPI)
│   ├── evaluation/                 <- 256+ result JSONs + SHAP/WGCNA outputs
│   ├── wgcna_outputs/              <- Per-tissue WGCNA module data (6 tissues)
│   └── figures/html/               <- 11 D3.js interactive HTML figures
│
├── v5/                             <- Biological interpretation layer
│   ├── scripts/                    <- Immune deconv, metabolic flux, drug targets
│   ├── evaluation/                 <- 25 result JSONs (immune, TF, metabolic, drugs)
│   └── figures/html/               <- 5 D3.js interactive HTML figures
│
├── v6/                             <- Cross-species validation
│   ├── scripts/                    <- 7 Python scripts (phases A–F)
│   └── evaluation/                 <- 5 result JSONs
│
└── processed/                      <- Intermediate analysis outputs
    ├── A_detection/                <- Per-tissue LOMO data
    ├── B_cross_mission/            <- Transfer matrices
    ├── fgsea/                      <- 60+ fGSEA results
    └── pathway_scores/             <- 48 ssGSEA pathway score files
```

---

## Getting Started

Feature matrices (train_X.csv, test_X.csv) are hosted on HuggingFace due to size (~2 GB).
Labels, metadata, and fold structure are in this repository.

### Option A -- Load from HuggingFace (recommended)

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

### Option B -- Reproduce from OSDR raw data

Requires R 4.2+ with Bioconductor. See [docs/r_dependencies.md](docs/r_dependencies.md).

```bash
# 1. Download raw data from NASA OSDR
python scripts/fetch_osdr.py --osd OSD-238 OSD-239 OSD-243 OSD-254

# 2. Normalize (DESeq2)
Rscript scripts/normalize_rr7_skin.R

# 3. Quality filter + build all_missions
python scripts/quality_filter.py --tissue skin

# 4. Generate LOMO folds
python scripts/generate_tasks.py --task A5
```

---

## Evaluation Protocol

All submissions are evaluated with:

| Metric | Description | Go threshold |
|--------|-------------|-------------|
| Mean AUROC | Average AUROC across folds | > 0.700 |
| 95% CI lower | Bootstrap CI (N=2000) lower bound | > 0.500 |
| perm_p | Permutation p-value (N=1000, pseudocount) | < 0.050 |

All three conditions must pass for a GO decision.

### Submission Format

Prepare a JSON file (see [docs/submission_format.md](docs/submission_format.md)):

```json
{
  "task_id": "A5",
  "model_name": "MyModel_v1",
  "predictions": {
    "fold_MHU-2_test": {"sample_id_1": 0.92, "sample_id_2": 0.07},
    "fold_RR-6_test":  {"...": "..."},
    "fold_RR-7_test":  {"...": "..."}
  }
}
```

```bash
python scripts/evaluate_submission.py \
    --submission my_submission.json \
    --task A5
```

---

## Model Tracks

| Track | Examples | Input Format |
|-------|---------|-------------|
| **Tier 1 -- Classical ML** | LR, RF, XGBoost, PCA-LR, SVM, LightGBM, TabNet | Tabular gene x sample |
| **Tier 2 -- Foundation Models** | Geneformer, scGPT, UCE, scFoundation | Gene rank order / embeddings |
| **Tier 3 -- Text LLMs** | GPT-4o, Claude, Llama 3 | Natural language gene list |

---

## Data

All data is derived from publicly available NASA OSDR datasets.

### v1 Core (6 tissues, 24 studies)

| Tissue | OSD Accession | Mission | n samples |
|--------|--------------|---------|-----------|
| Liver | OSD-48, 137, 245, 379, 242, 686 | RR-1, RR-3, RR-6, RR-8, RR-9, MHU-2 | 193 (v1) / 261 (v4 w/ BC/VC) |
| Gastrocnemius | OSD-101, 401, 326 | RR-1, RR-5, RR-9 | 32 |
| Kidney | OSD-102, 163, 253 | RR-1, RR-3, RR-7 | 118 |
| Thymus | OSD-289, 244, 421 | MHU-1, MHU-2, RR-6, RR-9 | 67 |
| Skin | OSD-238, 239, 243, 254 | MHU-2, RR-6, RR-7 | 102 |
| Eye | OSD-100, 194, 397 | RR-1, RR-3, TBD | 37 |

### v4 Additions (2 tissues)

| Tissue | OSD Accession | Mission | n samples |
|--------|--------------|---------|-----------|
| Lung | Multiple OSD studies | Multiple missions | See v4 results |
| Colon | Multiple OSD studies | Multiple missions | See v4 results |

Preprocessing: DESeq2 normalization (per-mission), log2(counts + 1), global low-expression filter (>=20% samples with count>1), top 75th percentile variance gene selection per fold (train missions only).

---

## Design Decisions

Key methodological choices underpinning this benchmark:

- **Feature encoding**: log2(DESeq2 normalized counts); LFC forbidden in Category A (label leakage)
- **LOMO splits**: Mission = independence unit; variance filter applied on training missions only (no test leakage)
- **Strain**: Track 2a = C57BL/6J only; Track 2b = all strains
- **Evaluation**: AUROC + bootstrap CI (N=2000) + permutation p (N=1000); GO requires all 3 AND conditions
- **Negative controls**: NC1 label permutation, NC2 housekeeping genes, NC3 cross-species gene set
- **v4 classifiers**: PCA-LR, ElasticNet-LR, RF, XGBoost, SVM-Linear, SVM-RBF, TabNet, LightGBM
- **v4 features**: Gene (log2-norm), Hallmark (ssGSEA), KEGG (ssGSEA), pathway-combined
- **Pathway analysis**: fGSEA group-level + gseapy ssGSEA sample-level
- **Text LLM track**: Gene list → natural language prompt → binary classification
- **Category B**: Transfer Pattern Summary (no single GO/NO-GO; use pairs ≥0.70 and perm_p<0.05 counts)

---

## Changelog

| Version | Date | Changes |
|---------|------|---------|
| v6.0 | 2026-03-30 | Cross-species validation complete: gene conservation, pathway conservation, cross-species transfer, TF conservation, biomarker validation, drug target validation (6 result JSONs). |
| v5.0 | 2026-03-29 | Biological interpretation complete: immune deconvolution (mMCP-counter, 8 tissues), cross-organ signaling (OmniPath, 111 L–R pairs), metabolic flux (iMM1865 E-Flux + pFBA), drug targets (DGIdb + ChEMBL, 1,284 FDA interactions), consensus 20-gene biomarker panel (gastro AUROC 0.806). 5 integration figures. |
| v4.0 | 2026-03-28 | v4 complete: 256 evaluations (8 tissues × 8 methods × 4 features). PCA-LR best (AUROC 0.776). Friedman p=0.015. Ablation (569 evals), SHAP multi-method, Python WGCNA (6 tissues), module preservation, STRING PPI. 11 publication figures. |
| v3.0 | 2026-03-20 | Multi-species (E4 Drosophila KEGG), spatial Visium brain (NEGATIVE), RRRM-2 scRNA-seq (PBMC NK 0.845), UCE + scFoundation FM eval, 7×7 cross-tissue transfer, radiation analogs. 5 publication figures. |
| v2.0 | 2026-03-18 | RRRM-1 scRNA-seq 4 tissues (38K cells). E1 cross-species r=0.352. T1–T3 temporal dynamics. 3 publication figures. |
| v1.3 | 2026-03-13 | Tier 3 LLM zero-shot. scGPT mean 0.667. Held-out: thymus RR-23 (0.905), skin RR-7 (0.885). 4 main + 4 supplementary figures. |
| v1.0 | 2026-03-07 | Initial release: 6 tissues, Categories A–D, Geneformer, Cell 2020 validation, fGSEA, GSVA. |

---

## Citation

*(Manuscript in preparation)*

```bibtex
@dataset{kang2026genelab,
  title   = {GeneLab Benchmark: A Multi-Tissue Spaceflight Transcriptomics Benchmark for AI/ML Models},
  author  = {Kang, Jaeyoung},
  year    = {2026},
  url     = {https://huggingface.co/datasets/jang1563/genelab-benchmark},
  note    = {v6.0}
}
```

Data source: NASA Open Science Data Repository (OSDR) -- [osdr.nasa.gov](https://osdr.nasa.gov/bio/repo/)

---

## License

Code: MIT License
Data: NASA OSDR public data (see individual dataset licenses at OSDR)

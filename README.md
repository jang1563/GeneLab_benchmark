# GeneLab Benchmark

**A public benchmark for evaluating AI/ML and Foundation Models on NASA OSDR spaceflight transcriptomics data.**

Version: v1.3 (Dataset freeze: 2026-03-01)
Status: **v1 Analysis complete** — 6 tissues, 3 model tiers, 2 independent held-out validations | **v2 in progress** — temporal dynamics, cross-species, RRRM-1 scRNA-seq

[![Dataset on HuggingFace](https://img.shields.io/badge/HuggingFace-Dataset-yellow)](https://huggingface.co/datasets/jang1563/genelab-benchmark)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

---

## What Is This?

GeneLab Benchmark provides standardized tasks for evaluating how well machine learning models — from classical baselines to gene expression foundation models (Geneformer, scGPT) and text-based LLMs (GPT-4o, Claude) — generalize **spaceflight transcriptomic signatures across missions**.

**Core challenge**: Train a model on one spaceflight mission's RNA-seq data. Can it classify samples from a different mission it has never seen?

**Data source**: [NASA Open Science Data Repository (OSDR)](https://osdr.nasa.gov/bio/repo/) — mouse multi-tissue bulk RNA-seq from ISS and ground control missions.

### Benchmark Scope

- **6 tissues**: Liver, Gastrocnemius, Kidney, Thymus, Skin, Eye
- **17 ISS missions**: RR-1 through RR-9, MHU-1, MHU-2, and more
- **24 verified OSD studies**, ~450 samples (binary Flight/Ground)
- **25+ evaluation tasks** across 7 categories (A–D, J, NC, Validation)

---

## Key Features

- **Leave-One-Mission-Out (LOMO)** cross-validation — mission = independence unit, preventing cross-mission data leakage
- **Category A**: Spaceflight detection per tissue (binary: Flight vs. Ground)
- **Category B**: Cross-mission transfer matrix (train on mission i, test on mission j) for all 6 tissues
- **Category C**: Cross-tissue transfer (3 methods: gene, DEG, pathway)
- **Category D**: Condition/confounder prediction (mission, strain, hardware, gravity)
- **3-tier model evaluation**: Classical ML → Gene Expression Foundation Models → Text LLMs
- **Standardized submission format** with automatic AUROC/CI/p-value evaluation
- **Biological validation**: NES pathway conservation, Cell 2020 concordance, negative controls

---

## Results Summary (v1.3)

### Category A — Spaceflight Detection (LOMO)

**Gene-level (primary)**:

| Task | Tissue | Missions | Method | Mean AUROC | 95% CI lower | perm_p | Decision |
|------|--------|---------|--------|-----------|-------------|--------|----------|
| A4 | Thymus | 4† | PCA-LR | **0.923** | 0.878 | 0.037 | ✓ **GO** |
| A2 | Gastrocnemius | 3 | LR | **0.907** | 0.717 | 0.026 | ✓ **GO** |
| A5 | Skin | 3§ | LR | **0.821** | 0.637 | 0.0023 | ✓ **GO** |
| A6 | Eye | 3 | LR | 0.811 | 0.470 | 0.063 | ✗ NO-GO‡ |
| A1 | Liver | 6 | LR | 0.653 | 0.457 | 0.091 | ✗ NO-GO |
| A3 | Kidney | 3 | LR | 0.593 | 0.431 | 0.281 | ✗ NO-GO |

**Pathway-level (GSVA Hallmark, secondary)**:

| Task | Tissue | Method | Mean AUROC | 95% CI lower | perm_p | Decision |
|------|--------|--------|-----------|-------------|--------|----------|
| **A6** | **Eye** | **PCA-LR** | **0.915** | **0.745** | **0.014** | ✓ **GO**⊕ |
| A3 | Kidney | LR | 0.755 | 0.481 | 0.071 | ✗ NO-GO |

<details>
<summary>Footnotes</summary>

- †A4 includes MHU-1 (Track 2b, GC/FLT strain mismatch — see PHASE1_RESULTS.md)
- §A5: MHU-2 = dorsal (OSD-238) + femoral (OSD-239) merged; RR-7 = OSD-254 C57BL/6J non-BSL subset (n=30)
- ‡A6 gene-level: AUROC passes but CI lower fails (n=9–16 per fold)
- ⊕A6 pathway-level: GSVA Hallmark 50-pathway scores rescue CI lower (0.470→0.745). Oxidative phosphorylation dominant.

</details>

### Category B — Cross-Mission Transfer (PCA-LR)

| Task | Tissue | N pairs | Mean AUROC | 95% CI | AUROC≥0.70 | Tier |
|------|--------|---------|-----------|--------|-----------|------|
| B4 | Thymus | 12 | **0.860** | [0.763, 0.953] | 9/12 | 1 |
| B2 | Gastrocnemius | 6 | **0.801** | [0.653, 0.944] | 4/6 | 1 |
| B5 | Skin | 6 | **0.772** | [0.691, 0.834] | 5/6 | 2 |
| B6 | Eye | 6 | 0.754 | [0.688, 0.838] | 5/6 | 2 |
| B1 | Liver | 30 | 0.577 | [0.492, 0.666] | 13/30 | 3 |
| B3 | Kidney | 6 | 0.555 | [0.397, 0.681] | 2/6 | 3 |

### Category C — Cross-Tissue Transfer (3 Methods)

| Pair | Method A (Gene) | Method B (DEG) | Method C (Pathway) | Best |
|------|--------|--------|--------|--------|
| C1: liver→kidney | **0.730** | 0.441 NS | 0.483 NS | A |
| C2: liver→gastro | 0.563 NS | 0.676 | **0.867** | C |
| C3: liver→thymus | 0.350 NS | **0.621** | 0.184 (anti) | B |
| C4: thymus→kidney | 0.585 NS | 0.539 NS | **0.690** | C |

### Category D — Condition/Confounder Prediction (macro-F1)

| Task | Tissue | N | Gene F1 | Pathway F1 | Gene p | Interpretation |
|------|--------|---|---------|-----------|--------|---------|
| D3 Mission ID (6-class) | Liver | 264 | **1.000** | 0.056 NS | <0.001 | Perfect batch separation; pathways batch-invariant |
| D4 Strain (2-class) | Thymus | 34 | **0.892** | 0.817 | 0.004 | Strain detectable from GC-only. EXPLORATORY (n_minority=3) |
| D5 Hardware (RR vs MHU) | Liver | 264 | **1.000** | 0.386 NS | <0.001 | Perfect gene separation; collinear with D3 |
| D5 Hardware (RR vs MHU) | Thymus | 92 | **1.000** | 0.352 NS | <0.001 | Perfect gene separation; collinear with D3 |
| D6 Gravity (3-class) | Liver | 9 | **0.886** | 0.413 NS | 0.002 | Microgravity separable from artificial gravity |
| D6 Gravity (3-class) | Thymus | 9 | **0.657** | 0.641 | 0.037 | Gene ≈ Pathway for gravity detection |

**Confounder hierarchy**: D3 (mission F1=1.0) ≥ D5 (hardware F1=1.0, collinear) ≥ D4 (strain F1=0.89, exploratory). All pathway F1 ≈ 0.05–0.41 → pathways resist confounder detection.

### J5 — Gene-level vs Pathway-level (15 comparisons)

| Category | N | Gene wins | Pathway wins | Mean diff |
|---|---|---|---|---|
| A (Detection) | 5 | 3 | 2 | +0.032 |
| C (Cross-tissue) | 4 | 2 | 2 | -0.001 |
| D (Condition, D3–D6) | 6 | 6 | 0 | -0.462 |
| **Total** | **15** | **11** | **4** | **-0.174** |

Notable finding — **"Kidney Rescue"**: gene-level AUROC=0.43 (fail) → pathway-level AUROC=0.74 (success, +0.31). Eye shows similar rescue (0.79→0.92, +0.13).

See [PHASE1_RESULTS.md](docs/development_history/PHASE1_RESULTS.md) for full results including per-fold tables, SHAP analysis, and pathway analysis.

---

## Key Scientific Findings

### Pre-registered Hypotheses

| Hypothesis | Statement | Verdict | Key Evidence |
|---|---|---|---|
| **H1** | Liver has the most consistent cross-mission transcriptome | **REFUTED** | Thymus (0.860) >> Liver (0.577). Thymus and Gastrocnemius = Tier 1. |
| **H2** | Transfer failure from biological diversity, not batch effects | **SUPPORTED** | NES conservation r=0.9 (5 tissues). D3 pathway F1=0.06 (batch-invariant). limma_rbe mean delta=0.01. |
| **H3** | Pathway-level preserves spaceflight response better than gene-level | **CONDITIONALLY SUPPORTED** | Kidney rescue (0.43→0.74), Eye (0.79→0.92). But tissue-pair dependent. |

### NES Pathway Conservation vs Transfer Success

Normalized Enrichment Score (NES) correlation between mission pairs predicts cross-mission transfer performance:

| Tissue | NES Mean r | Transfer AUROC | Spearman |
|---|---|---|---|
| Thymus | 0.619 | 0.860 | |
| Eye | 0.335 | 0.754 | |
| Skin | 0.147 | 0.772 | |
| Liver | 0.059 | 0.577 | |
| Kidney | 0.048 | 0.555 | |

5-tissue Spearman r = 0.9 (excluding gastrocnemius, which has incomplete fGSEA data). Original 4-tissue r = 1.0.

### External Validation (Cell 2020)

Validated against Beheshti et al. (Cell 2020, PMID 33242417) multi-omics consensus:
- **Pathway direction concordance**: 71.7% across 5 tissues (STRONG agreement)
- **Gene SHAP top-50 overlap**: 10.7% (47× above random chance)
- Tissue-specific: Thymus/Gastrocnemius 100%, Liver/Eye 67%, Kidney 25%

### Negative Controls (all PASS)

| Control | Method | Expected | Result |
|---|---|---|---|
| NC1 | Permutation test (28 entries) | AUROC ≈ 0.50 | 0.50 ± 0.03 |
| NC2 | Housekeeping genes only (50 genes) | AUROC ≈ 0.50 | 0.49–0.55 |

### Independent Held-Out Validation

Two missions were withheld entirely from all model development:

| Tissue | Mission | AUROC | 95% CI | p-value | Duration | n_test |
|--------|---------|-------|--------|---------|---------|--------|
| Thymus | RR-23 | **0.905** | [0.672, 1.000] | 0.005 | 30 days | 18 (10 FLT, 8 GC) |
| Skin | RR-7 | **0.885** | [0.745, 0.986] | <0.001 | 75 days | 30 (10 FLT, 20 GC) |

Both held-out evaluations pass the GO criteria (AUROC>0.700, CI_lower>0.500, perm_p<0.050), confirming that PCA-LR classifiers trained on 3–12 ISS missions generalize to completely unseen missions including the longest mission evaluated (75 days).

### Biological Validation (fGSEA Hallmark)

| Tissue | Top Enriched Pathways | Consistency |
|---|---|---|
| Liver | OXIDATIVE_PHOSPHORYLATION, FATTY_ACID_METABOLISM | Literature-concordant |
| Thymus | E2F_TARGETS, G2M_CHECKPOINT, IFN-gamma | Thymocyte proliferation |
| Gastrocnemius | OXIDATIVE_PHOSPHORYLATION, MYOGENESIS | Muscle metabolism |
| Kidney | MTORC1_SIGNALING, CHOLESTEROL_HOMEOSTASIS | Renal metabolism |
| Eye | OXIDATIVE_PHOSPHORYLATION (dominant 3/3 missions) | Retina metabolic demand |
| Skin | E2F_TARGETS, G2M_CHECKPOINT, EPITHELIAL_MESENCHYMAL_TRANSITION | Cell proliferation + ECM remodeling |

---

## Repository Structure

```
GeneLab_benchmark/
├── README.md                       ← This file
├── PLAN.md                         ← Benchmark design specification (v1.5)
├── DESIGN_DECISIONS.md             ← Architecture decisions log (DD-01 to DD-21)
├── DATA_CATALOG.md                 ← Auto-generated OSDR inventory (24 studies)
├── CITATION.cff                    ← Citation metadata
│
├── tasks/                          ← Public task inputs (17 directories)
│   ├── A1_liver_lomo/              ← 6 folds + 3 variants (standard, ComBat, ISS-only)
│   ├── A2_gastrocnemius_lomo/      ← 3 folds
│   ├── A3_kidney_lomo/             ← 3 folds
│   ├── A4_thymus_lomo/             ← 4 folds + holdout
│   ├── A5_skin_lomo/               ← 3 folds
│   ├── A6_eye_lomo/                ← 3 folds
│   └── B1–B6_*_cross_mission/     ← N×(N-1) mission pairs per tissue
│
├── scripts/                        ← v1 pipeline scripts (35 Python/R/shell)
│   ├── run_baselines.py            ← Classical ML baseline runner (LR, RF, XGBoost, PCA-LR)
│   ├── evaluate_submission.py      ← Submission evaluator (AUROC, CI, perm_p)
│   ├── generate_tasks.py           ← LOMO split generator
│   ├── cross_mission_transfer.py   ← Category B matrix generator
│   ├── cross_tissue_transfer.py    ← Category C: 3 methods
│   ├── condition_prediction.py     ← Category D: mission/strain/hardware/gravity
│   ├── gene_vs_pathway_comparison.py ← J5: feature representation
│   ├── dge_pipeline_comparison.py  ← J2: DESeq2 vs edgeR vs limma-voom
│   ├── shap_analysis.py            ← SHAP feature importance
│   ├── run_fgsea.R                 ← Group-level fGSEA enrichment
│   ├── compute_pathway_scores.R    ← Sample-level GSVA scores
│   ├── batch_correction_eval.py    ← J3: ComBat-seq, limma, RUVseq
│   ├── housekeeping_control.py     ← NC2: housekeeping gene baseline
│   ├── cell2020_validation.py      ← External validation vs Cell 2020
│   ├── compute_nes_conservation.py ← NES pathway conservation
│   ├── geneformer_tokenize.py      ← Tier 2: gene rank tokenization
│   ├── geneformer_finetune.py      ← Tier 2: Mouse-GF BERT fine-tuning
│   ├── scgpt_tokenize.py           ← Tier 2: scGPT ortholog mapping + tokenization
│   ├── scgpt_finetune.py           ← Tier 2: scGPT Transformer fine-tuning
│   ├── run_holdout.py              ← Independent held-out evaluation
│   └── utils.py                    ← Shared utilities
│
├── docs/
│   ├── BIOLOGICAL_GROUND_TRUTH.md  ← Validation reference (Cell 2020, SOMA 2024)
│   ├── submission_format.md        ← JSON submission specification
│   ├── text_llm_format.md          ← Text LLM evaluation format (DD-16)
│   ├── hf_dataset_card.md          ← HuggingFace dataset documentation
│   └── development_history/
│       └── PHASE1_RESULTS.md       ← Full Phase 1 analysis
│
├── evaluation/                     ← Result JSON files (~50)
│   ├── A*_baseline_results.json    ← Per-tissue classical ML results
│   ├── A*_shap_rf.json             ← SHAP feature rankings
│   ├── B_cross_mission_summary.json
│   ├── C_cross_tissue_summary.json
│   ├── D_condition_summary.json
│   ├── J2_dge_pipeline_comparison.json   ← DESeq2 vs edgeR vs limma
│   ├── J3_batch_correction_comparison.json
│   ├── J5_gene_vs_pathway.json
│   ├── multi_db_comparison.json    ← 4 pathway DBs × 6 tissues LOMO
│   ├── geneformer_mouse_gf_all_tissues_summary.json
│   ├── scgpt_whole_human_all_tissues_summary.json
│   ├── A4_holdout_results.json     ← Thymus RR-23 held-out
│   ├── A5_holdout_results.json     ← Skin RR-7 held-out
│   ├── NC1_permutation_summary.json
│   ├── NC2_housekeeping_summary.json
│   ├── cell2020_validation.json
│   ├── NES_conservation_vs_transfer.json
│   ├── RESULTS_SUMMARY.md          ← Comprehensive results table
│   └── submission_*.json           ← Baseline + FM submission files
│
├── v2/                             ← v2.0 extension (temporal, cross-species, scRNA-seq)
│   ├── README.md                   ← v2 overview and status
│   ├── scripts/                    ← v2-specific scripts (19 Python)
│   │   ├── temporal_analysis.py    ← T1/T2/T3 temporal dynamics
│   │   ├── cross_species_nes_comparison.py ← E1/E2 cross-species NES
│   │   ├── e3_cfrna_celltype_origin.py     ← E3 cfRNA origin
│   │   ├── generate_main_figures.py        ← v2 integrated figures (D3.js)
│   │   ├── llm_*.py                ← Tier 3 LLM pipeline (3 scripts)
│   │   ├── rrrm1_*.py              ← RRRM-1 scRNA-seq pipeline (7 scripts)
│   │   └── i4_*.py                 ← I4 PBMC analysis (2 scripts)
│   ├── evaluation/                 ← v2 results
│   │   └── V2_RESULTS_SUMMARY.md   ← v2 comprehensive results
│   ├── docs/                       ← v2 design and planning documents
│   │   ├── V2_PAPER_CONTENT.md     ← v2 manuscript draft
│   │   ├── DATA_CATALOG_V2.md      ← v2 data sources
│   │   ├── RRRM1_SC_BENCHMARK_PLAN_V3.md ← scRNA-seq benchmark plan
│   │   └── RRRM1_*.md             ← RRRM-1 annotation/marker docs
│   └── processed/                  ← v2 intermediate outputs
│       ├── T_temporal/             ← T1/T2/T3 results
│       ├── E_crossspecies/         ← E1/E2/E3 results
│       └── F1_scrna/               ← I4 PBMC cell-type fGSEA
│
└── processed/                      ← v1 intermediate analysis outputs
    ├── A_detection/                ← Per-tissue LOMO data
    ├── B_cross_mission/            ← Transfer matrices + CI
    ├── C_cross_tissue/             ← 4 pairs × 3 methods
    ├── D_condition/                ← Condition prediction
    ├── fgsea/                      ← 60 fGSEA results (6 tissues × missions × 3 DBs)
    ├── pathway_scores/             ← 54 GSVA files (5 tissues × missions × 3 DBs)
    └── qc_reports/
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
# A1 has multiple variants; select one explicitly
python scripts/run_baselines.py --task A1 --task-dir A1_liver_lomo --model lr
```

### 3. Submit your model's predictions

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

Evaluate:

```bash
python scripts/evaluate_submission.py \
    --submission my_submission.json \
    --task A5

# A1 example (variant must be explicit)
python scripts/evaluate_submission.py \
    --submission my_submission.json \
    --task A1 \
    --task-dir A1_liver_lomo
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

### Category C — Cross-Tissue Transfer

**Goal**: Train on tissue X, predict spaceflight status on tissue Y. Evaluates whether spaceflight signatures are shared across tissues.

Three transfer methods:
- **Method A (Gene)**: Direct gene intersection transfer
- **Method B (DEG)**: Differentially expressed gene overlap
- **Method C (Pathway)**: GSVA Hallmark pathway score transfer

### Category D — Condition/Confounder Prediction

**Goal**: Predict confounding variables (mission identity, strain, hardware, gravity level) to quantify batch effects and biological confounders.

- **D3**: Mission ID (liver, 6-class) — batch effect quantification
- **D4**: Strain (thymus GC, C57BL/6J vs C57BL/6CR) — exploratory (n=3)
- **D5**: Hardware (RR vs MHU, liver + thymus) — collinear with D3
- **D6**: Gravity (MHU-2, uG/AG/GC, liver + thymus) — biological signal

Key finding: D3 gene F1=1.0 (perfect mission separation) vs pathway F1=0.06 (batch-invariant) confirms pathways absorb batch effects.

---

## Baseline Submissions

Pre-computed baseline predictions are available in `evaluation/` for reference and reproducibility.

**Category A (LOMO)**

| File | Task | Model | Mean AUROC | Go/No-Go |
|------|------|-------|------------|----------|
| `submission_PCA-LR_baseline_v1_A4_eval.json` | A4 Thymus | PCA-LR (L2, lbfgs) | 0.923 | ✓ GO |
| `submission_LR-ElasticNet_baseline_v1_A2_eval.json` | A2 Gastrocnemius | LR-ElasticNet (SAGA) | 0.917 | ✓ GO |
| `submission_PCA-LR_baseline_v1_A5_eval.json` | A5 Skin | LR (ElasticNet) | 0.821 | ✓ GO |
| `submission_PCA-LR_baseline_v1_A6_eval.json` | A6 Eye | PCA-LR (pathway) | 0.915 | ✓ GO |

**Category B (Cross-Mission Transfer)**

| Task | Tissue | N pairs | PCA-LR Mean AUROC | LFC Mean AUROC |
|------|--------|---------|------------------|---------------|
| B4 | Thymus | 12 | 0.860 | 0.868 |
| B2 | Gastrocnemius | 6 | 0.801 | 0.655 |
| B5 | Skin | 6 | 0.772 | 0.750 |
| B6 | Eye | 6 | 0.754 | 0.696 |
| B1 | Liver | 30 | 0.577 | 0.534 |
| B3 | Kidney | 6 | 0.555 | 0.465 |

> Category B does not report a single GO/NO-GO — see [DD-17](DESIGN_DECISIONS.md) for evaluation criteria.

Evaluate a baseline submission:

```bash
# Category A
python scripts/evaluate_submission.py \
    --submission evaluation/submission_PCA-LR_baseline_v1_A5_eval.json \
    --task A5

# Category B (summary across all tissues)
python scripts/cross_mission_transfer.py --tissue skin
```

> **Reproducibility note**: `A2_baseline_results.json` was computed with `max_iter=2000` (SAGA not fully converged for 15k genes). The baseline submission (`submission_LR-ElasticNet_baseline_v1_A2_eval.json`) uses `max_iter=10000` (converged); A2 mean AUROC improves from 0.907 → 0.917. GO/No-Go conclusion unchanged.

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
| **Tier 2 — Foundation Models** | Geneformer (Mouse-GF) | Gene rank order (tokenized) |
| **Tier 3 — Text LLMs** | GPT-4o, Claude, Llama 3 | Natural language gene list (see DD-16) |

### Tier 2 Results: Foundation Models vs Classical ML (LOMO AUROC)

| Tissue | Baseline | scGPT | Δ scGPT | Mouse-GF | Δ GF | Winner |
|--------|---------|-------|---------|---------|------|--------|
| Liver | 0.588 | 0.628 | +0.040 | 0.486 | -0.102 | scGPT |
| Gastrocnemius | 0.907 | 0.685 | -0.222 | 0.382 | -0.525 | Baseline |
| Kidney | 0.521 | 0.556 | +0.035 | 0.452 | -0.069 | scGPT |
| Thymus | 0.923 | 0.782 | -0.141 | 0.495 | -0.428 | Baseline |
| Skin | 0.821 | 0.691 | -0.130 | 0.557 | -0.265 | Baseline |
| Eye | 0.789 | 0.650 | -0.139 | 0.484 | -0.305 | Baseline |
| **Mean** | **0.758** | **0.667** | **-0.092** | **0.476** | **-0.283** | **Baseline** |

**scGPT** (12-layer Transformer, 33M parameter pretrained on 33M human cells, whole-human model) outperforms Mouse-Geneformer by +0.191 despite requiring mouse→human ortholog mapping (18,398 shared genes). Scale > species alignment. Both FMs underperform classical ML on mean AUROC.

**Mouse-Geneformer** (6-layer BERT, 56K gene vocab, pretrained on 30M mouse scRNA-seq cells): lowest mean AUROC (0.476). Classical ML wins 6/6 tissues, scGPT wins 4/6.

Finding: **scRNA-seq-pretrained FMs do not transfer to small-n (n=30–100) bulk transcriptomics**, regardless of scale or species alignment.

### Tier 3 Results: Text LLM Zero-Shot (LOMO AUROC)

| Model | Mean AUROC (6 tasks) |
|-------|---------------------|
| DeepSeek-V3 | 0.505 |
| Gemini-2.5-Flash | 0.501 |
| Llama-3.3-70B | 0.471 |
| PCA-LR baseline | **0.757** |

All three LLMs perform near chance (0.47–0.51) on binary flight/ground classification from gene name lists. Text-based zero-shot classification provides no signal beyond random. Input format: sorted differential expression gene list (see DD-16).

For Tier 3 input format specification, see [DESIGN_DECISIONS.md](DESIGN_DECISIONS.md) (DD-16).

---

## Data

All data is derived from publicly available NASA OSDR datasets (24 studies, 6 tissues).

| Tissue | OSD Accession | Mission | n samples | Note |
|--------|--------------|---------|-----------|------|
| Liver | OSD-48 | RR-1 | 18 | Track 2a |
| Liver | OSD-137 | RR-3 | 20 | Track 2a |
| Liver | OSD-245 | RR-6 | 48 | Track 2a |
| Liver | OSD-379 | RR-8 | 40 | Track 2a |
| Liver | OSD-242 | RR-9 | 39 | Track 2a |
| Liver | OSD-686 | MHU-2 | 28 | Track 2a (uG/GC/AG 3-group) |
| Gastrocnemius | OSD-101 | RR-1 | 12 | Track 2a |
| Gastrocnemius | OSD-401 | RR-5 | 12 | Track 2a |
| Gastrocnemius | OSD-326 | RR-9 | 8 | Track 2a |
| Kidney | OSD-102 | RR-1 | 47 | Track 2a |
| Kidney | OSD-163 | RR-3 | 32 | Track 2a |
| Kidney | OSD-253 | RR-7 | 39 | Track 2a |
| Thymus | OSD-289 | MHU-1 | 6 | Track 2b (GC = C57BL/6CR) |
| Thymus | OSD-289 | MHU-2 | 6 | Track 2a |
| Thymus | OSD-244 | RR-6 | 35 | Track 2a |
| Thymus | OSD-421 | RR-9 | 20 | Track 2a |
| Skin | OSD-238 | MHU-2 (dorsal) | 18 | merged as "MHU-2" (6F+6GC+6VC; AG excluded) |
| Skin | OSD-239 | MHU-2 (femoral) | 17 | merged as "MHU-2" (5F+12GC; AG excluded) |
| Skin | OSD-243 | RR-6 | 37 | Track 2a |
| Skin | OSD-254 | RR-7 | 30 | C57BL/6J non-BSL subset only |
| Eye | OSD-100 | RR-1 | 12 | Track 2a |
| Eye | OSD-194 | RR-3 | 9 | Track 2a |
| Eye | OSD-397 | TBD | 16 | Track 2a |

Preprocessing: DESeq2 normalization (per-mission), log2(counts + 1), global low-expression filter (≥20% samples with count>1), top 75th percentile variance gene selection per fold (train missions only — DD-03).

---

## Execution Safety Defaults (2026-03)

- `run_baselines.py` and `shap_analysis.py` exclude `fold_*_holdout` by default.
- `shap_analysis.py` includes holdout only with `--include-holdout`.
- `evaluate_submission.py` accepts holdout predictions if provided, but does not require them.
- If one task ID matches multiple directories (for example A1), scripts now raise an ambiguity error unless `--task-dir` is provided.
- Geneformer `mouse_gf` path is configurable:
  - Tokenize: `--mouse-gf-base` or env `MOUSE_GF_BASE`
  - Finetune: `--mouse-gf-model-dir` or env `MOUSE_GF_MODEL_DIR`

---

## Design Decisions

Key methodological choices are documented in [DESIGN_DECISIONS.md](DESIGN_DECISIONS.md):

- **DD-01**: Feature = log2(DESeq2 normalized counts) — LFC forbidden in Category A (label leakage)
- **DD-03**: LOMO-aware variance filter (train missions only — no test leakage)
- **DD-04**: Mission = independence unit for LOMO (not sample)
- **DD-06**: Track 2a = C57BL/6J only; Track 2b = all strains
- **DD-08**: Evaluation metrics (AUROC + bootstrap CI + permutation p)
- **DD-11**: Go/No-Go decision criteria (3 AND conditions)
- **DD-12**: Negative controls (NC1 permutation, NC2 housekeeping, NC3 cross-species)
- **DD-13**: Baseline model set (LR, RF, XGBoost, PCA-LR)
- **DD-15**: Pathway analysis (fGSEA group-level + GSVA sample-level)
- **DD-16**: Text LLM evaluation track specification
- **DD-17**: Category B evaluation criteria (Transfer Pattern Summary, perm_p floor)
- **DD-21**: scGPT fine-tuning strategy (freeze 10/12 layers, ortholog mapping, 10 epochs)

---

## Changelog

| Version | Date | Changes |
|---------|------|---------|
| v1.3 | 2026-03-13 | v2 temporal/cross-species phases complete (T1-T3, E1-E3, F1). 3 integrated v2 figures (D3.js). RRRM-1 scRNA-seq pipeline complete: STARsolo → QC → annotation → doublet removal for 4 tissues (blood/eye/muscle/skin, 38K cells). F2 benchmark tasks defined (composition/pseudo-bulk/classifier/cross-species). |
| v1.2 | 2026-03-10 | v2 Phase 1-5 complete: T1/T2/T3 temporal dynamics, E1 cross-species NES (r=0.352), E2 duration effect (Δr=+0.446), E3 cfRNA origin, F1 PBMC cell-type fGSEA. 3 integrated main figures. |
| v1.1 | 2026-03-09 | Analysis complete. scGPT (mean 0.667) + 3-way FM comparison. Held-out: skin RR-7 (0.885). J2 DGE pipeline comparison (ρ=0.926, Jaccard=0.600). Publication figures: fig1–4 + figS1–5 (D3.js v7). V1_PAPER_CONTENT.md v1.1. |
| v1.0 | 2026-03-07 | Tier 3 LLM zero-shot complete (3 providers × 6 tasks, mean 0.47–0.51). Held-out thymus RR-23 (AUROC=0.905). scGPT Tier 2 complete (6 tissues, mean AUROC=0.667). Multi-DB LOMO (4 DBs × 6 tissues). Mouse-GF Tier 2 (22 LOMO folds, Cayuga A40, mean AUROC=0.476). Category A–D + J3 + J5 + NC1/NC2. Cell 2020 external validation (71.7% concordance). fGSEA 80 files, GSVA 88 files. Dataset freeze: 2026-03-01. |

---

## Version Structure

| Version | Scope | Status | Location |
|---------|-------|--------|----------|
| **v1.0** | Mouse bulk RNA-seq, 6 tissues, 25+ tasks, 3 model tiers, 2 held-out validations | **Complete** | Project root (`scripts/`, `evaluation/`, `tasks/`) |
| **v2.0** | Temporal dynamics (T1-T3), cross-species (E1-E3), single-cell (F1, F2) | **In progress** | `v2/` directory |

v1.0 code is in project root. v2.0 extends with temporal/cross-species analysis (complete) and RRRM-1 scRNA-seq benchmarks (in progress). See [`v2/README.md`](v2/README.md) for full v2 scope and current status.

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

Data source: NASA Open Science Data Repository (OSDR) — [osdr.nasa.gov](https://osdr.nasa.gov/bio/repo/)

---

## License

Code: MIT License
Data: NASA OSDR public data (see individual dataset licenses at OSDR)

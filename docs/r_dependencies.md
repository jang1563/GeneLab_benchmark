# R Dependencies

GeneLab benchmark preprocessing requires **R 4.2+** with Bioconductor packages.
Python analysis (`run_baselines.py`, `evaluate_submission.py`) does not require R.

---

## Installation

```r
# Install Bioconductor package manager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Core preprocessing packages
BiocManager::install(c(
    "DESeq2",    # Size-factor normalization (per-mission)
    "limma",     # Linear model for differential expression
    "edgeR",     # Count data normalization utilities
    "GSVA",      # Gene Set Variation Analysis (pathway scores)
    "fgsea"      # Fast gene set enrichment
))

# CRAN packages
install.packages(c(
    "sva",       # ComBat-seq batch correction
    "ggplot2",   # Visualization
    "dplyr",     # Data manipulation
    "readr"      # CSV I/O
))
```

---

## R Scripts and Their Dependencies

| Script | R Packages Used | Purpose |
|--------|----------------|---------|
| `scripts/normalize_*.R` | DESeq2, edgeR | Per-mission size-factor normalization |
| `scripts/batch_correct.R` | sva (ComBat-seq) | Batch correction (Category J analysis) |
| `scripts/compute_pathway_scores.R` | GSVA, limma | GSVA Hallmark 50-pathway scores |
| `scripts/run_fgsea.R` | fgsea, limma | Fast gene set enrichment analysis |
| `scripts/_install_deps.R` | — | Interactive dependency installer |

---

## Notes

- **DESeq2 normalization** is applied per-mission (not joint across missions). See `DESIGN_DECISIONS.md` DD-10.
- **Python preprocessing** (`quality_filter.py`, `generate_tasks.py`) is applied after R normalization.
- R scripts produce per-mission `*_log2_norm.csv` files in `processed/A_detection/{tissue}/`.
- The final `{tissue}_all_missions_log2_norm.csv` is assembled in Python.
- R ≥ 4.2 required; Bioconductor 3.16+ recommended.

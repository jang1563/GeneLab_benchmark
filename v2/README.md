# GeneLab Benchmark v2.0

**Status**: Planning (v1.0 publication pending)

## Scope

v2.0 extends the benchmark beyond mouse bulk RNA-seq:

| Category | Description | Data |
|----------|-------------|------|
| **E** | Cross-Species Conservation | Mouse, Human, C. elegans, Arabidopsis |
| **F** | Single-Cell & Spatial | snRNA-seq, scRNA-seq, Visium, GeoMx |
| **G** | Microbiome | ISS environmental, mouse gut, human timeseries |

## Prerequisites

- v1.0 paper accepted/published
- C. elegans OSDR data quality verified
- Human data SpaceOmicsBench integration plan finalized

## Directory Structure

```
v2/
├── scripts/        # v2-specific analysis scripts
├── evaluation/     # v2 results and metrics
├── data/           # v2 additional data (cross-species, single-cell)
└── docs/           # v2 design documents
```

## Relationship to v1

- v1 code lives in the project root (`scripts/`, `evaluation/`, etc.)
- v1 is frozen at git tag `v1.0`
- Shared utilities from `scripts/utils.py` can be imported by v2 scripts
- v2 scripts should NOT modify v1 files

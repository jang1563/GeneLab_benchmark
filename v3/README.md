# GeneLabBench v3

**Start date**: 2026-03-18
**Status**: Execution in progress

## Directory Structure

```
v3/
  scripts/       # Analysis scripts (Python/R/Bash)
  evaluation/    # Output JSONs with metrics
  figures/       # D3.js HTML interactive figures
  processed/     # Intermediate data files
  docs/          # Plans, paper content, design decisions
```

## Version History

| Version | Focus | Date | Status |
|---------|-------|------|--------|
| v1.0 | Bulk RNA-seq benchmark (6 tissues, LOMO, 25+ tasks) | 2025-12 ~ 2026-02 | Complete |
| v1.3 | +Temporal (T1-T3) +Foundation Models (scGPT, Geneformer) +LLM | 2026-02 ~ 2026-03-09 | Complete |
| v2.0 | +Cross-species (E1-E3) +I4 PBMC (F1) +RRRM-1 scRNA-seq (F2) | 2026-03-07 ~ 2026-03-18 | Complete |
| **v3.0** | TBD | 2026-03-18 ~ | **Execution in progress** |

## v3 Scope

Defined in `v3/docs/PLAN_V3.md`

## Current Execution Notes

- Main plan: `v3/docs/PLAN_V3.md`
- HPC/FM runbook: `v3/docs/HPC_FM_RUNBOOK.md`
- Current status is mixed: some v3 analyses already have outputs, `B_ext` is complete, UCE is verified on Cayuga, and scFoundation is blocked only by external model-weight access.

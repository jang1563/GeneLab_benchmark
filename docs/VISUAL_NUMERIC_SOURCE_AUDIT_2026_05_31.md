# Visual Numeric Source Audit

Date: 2026-05-31

Purpose: verify that regenerated premium figure/table values match their cited source files.

Total checks: 117
Passing checks: 117
Failing checks: 0

Generated audit files:

- `output/premium_audit/premium_visual_numeric_audit.csv`
- `output/premium_audit/premium_visual_numeric_audit.json`

## Figure Summary

| Figure | Checks | Fails | Status |
|---|---:|---:|---|
| fig1_tissue_transfer | 24 | 0 | PASS |
| fig2_pathway | 32 | 0 | PASS |
| fig3_model | 18 | 0 | PASS |
| fig4_v9_platform | 9 | 0 | PASS |
| fig5_public_bulk_boundary | 6 | 0 | PASS |
| fig6_organoid | 28 | 0 | PASS |

## Failing Checks

No failing numeric checks.

## Interpretation

- A PASS means the scripted figure value matches the cited source within the audit tolerance.
- This audit does not validate whether a source file is itself final or publication-frozen.
- Source freeze and caption wording remain separate L4 requirements.

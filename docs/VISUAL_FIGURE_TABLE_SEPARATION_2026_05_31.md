# Visual Figure/Table Separation Pass

Date: 2026-05-31

Purpose: separate true figures from table-like content after review that
readiness matrices, decision tables, and release-gate lists feel awkward when
presented as "figures."

## Decision

Use this convention going forward:

- `Figure`: plots, schematics, workflows, spatial/temporal structure, model
  behavior, and visual comparisons.
- `Table`: readiness matrices, release gates, source inventories, decision
  matrices, claim-boundary language, and detailed summary rows.
- `Status slide`: acceptable for deck communication, but should not be treated
  as a manuscript main figure unless it is converted into a schematic.

This does not mean tables are low quality. It means tables should be formatted
and cited as tables, not disguised as figure panels.

## P0 Separation Applied

### Figure 1

Preferred figure:

- `output/premium_figures/premium_fig1_tissue_transfer_hierarchy.png`
- `output/premium_figures/premium_fig1_tissue_transfer_hierarchy.pdf`

Separated table:

- `output/premium_tables/table1_tissue_transfer_summary.csv`
- `output/premium_tables/table1_tissue_transfer_summary.json`

Change:

- Removed the right-side `Main readout` table from the figure.
- Kept the dot/interval plot as the figure.
- Moved summary statements such as highest transfer, near-chance transfer, and
  thymus-minus-liver difference into a table.

### Figure 4

Preferred figure:

- `output/premium_figures/premium_fig4_v9_platform_architecture.png`
- `output/premium_figures/premium_fig4_v9_platform_architecture.pdf`

Legacy filename regenerated for continuity:

- `output/premium_figures/premium_fig4_v9_platform_status.png`
- `output/premium_figures/premium_fig4_v9_platform_status.pdf`

Separated tables:

- `output/premium_tables/table4_v9_lane_readiness.csv`
- `output/premium_tables/table4_v9_lane_readiness.json`
- `output/premium_tables/table4_v9_evidence_footprint.csv`
- `output/premium_tables/table4_v9_evidence_footprint.json`

Change:

- Removed the lane readiness matrix and evidence-footprint list from the figure
  body.
- Kept a source-to-release architecture schematic and a lane-readiness
  schematic.
- Moved exact readiness rows and footprint counts to tables.

### Figure 5

Preferred figure:

- `output/premium_figures/premium_fig5_public_bulk_release_boundary_schematic.png`
- `output/premium_figures/premium_fig5_public_bulk_release_boundary_schematic.pdf`

Legacy filename regenerated for continuity:

- `output/premium_figures/premium_fig5_public_bulk_alpha_boundary.png`
- `output/premium_figures/premium_fig5_public_bulk_alpha_boundary.pdf`

Separated tables:

- `output/premium_tables/table5_public_bulk_alpha_gates.csv`
- `output/premium_tables/table5_public_bulk_alpha_gates.json`
- `output/premium_tables/table5_public_bulk_release_options.csv`
- `output/premium_tables/table5_public_bulk_release_options.json`
- `output/premium_tables/table5_public_bulk_claim_language.csv`
- `output/premium_tables/table5_public_bulk_claim_language.json`

Change:

- Replaced the gate ladder, release-option table, and allowed/prohibited
  language table with a sparse release-boundary schematic.
- Moved the gate details and claim-language rows to separate tables.

### Figure 6

Preferred figure:

- `output/premium_figures/premium_fig6_organoid_diagnostic_surface.png`
- `output/premium_figures/premium_fig6_organoid_diagnostic_surface.pdf`

Separated table:

- `output/premium_tables/table6_organoid_biology_check_decisions.csv`
- `output/premium_tables/table6_organoid_biology_check_decisions.json`

Change:

- Removed the decision table panel from the figure.
- Kept source footprint, separated metric mini-axes, and top-ranked model-gene
  enrichment.
- Moved default/secondary/negative biology-check decisions to a table.

## Manuscript Implication

Likely main figures after this pass:

1. Tissue-transfer hierarchy.
2. Pathway unwanted-signal reduction and selected task rescue.
3. Model-family comparison under small training-set constraints.
4. v9 platform architecture schematic.
5. Human organoid biology-check extension.

Likely tables:

1. Tissue-transfer summary table.
2. v9 lane readiness table.
3. v9 evidence-footprint table.
4. Public bulk release-gate table.
5. Public bulk claim-language table.
6. Organoid biology-check decision table.

Figure 5 should be treated as a release/status schematic, not a main
scientific result figure.


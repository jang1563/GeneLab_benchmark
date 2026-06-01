# Visual Manuscript Variant Review

Date: 2026-06-01

Purpose: record the first manuscript-specific rebuild and enlarged visual QA
for the Fig2/Fig3/Fig6 candidates. These variants are not downsampled slide exports;
they are newly laid out for a two-column manuscript figure width.

## New Manuscript Variants

Generated files:

- `output/premium_figures/manuscript_variants/premium_fig2_pathway_rescue_manuscript.png`
- `output/premium_figures/manuscript_variants/premium_fig2_pathway_rescue_manuscript.pdf`
- `output/premium_figures/manuscript_variants/premium_fig3_model_tier_comparison_manuscript.png`
- `output/premium_figures/manuscript_variants/premium_fig3_model_tier_comparison_manuscript.pdf`
- `output/premium_figures/manuscript_variants/premium_fig6_organoid_biology_check_manuscript.png`
- `output/premium_figures/manuscript_variants/premium_fig6_organoid_biology_check_manuscript.pdf`

New supporting table:

- `output/premium_tables/table3_model_surface_coverage.csv`
- `output/premium_tables/table3_model_surface_coverage.json`

## What Changed

### Fig2 Pathway Rescue

The manuscript variant replaces the deck-oriented bar-heavy layout with:

- a full-width paired-dot/lollipop panel for coupled-label prediction;
- a compact pathway-minus-gene AUROC delta panel;
- an exploratory pathway activity agreement scatter panel without a fitted trend.

Text changes:

- removed the legacy "artifact" framing from visible figure text;
- replaced it with "unwanted label signals";
- removed provenance and operating-language footers from the figure body;
- moved the Gene/Pathway key out of the data region.

Visual QA findings after enlarged inspection:

- title, subtitle, panel labels, and footer fit at manuscript width;
- Panel A no longer has a legend/data-label collision;
- Panel B negative delta labels are readable after moving labels away from the y labels;
- Panel C labels are readable, with the skin/gastrocnemius labels separated enough for draft use.

Remaining caption note:

- state that mission, hardware, and gravity label prediction are diagnostic
  checks on coupled experimental labels, not independent causal estimates.

### Fig3 Model Comparison

The manuscript variant removes the text-heavy coverage panel from the figure and
keeps the coverage warning in a compact footer plus a separate table.

The manuscript variant now uses:

- a full-width dot/lollipop model score panel with a zero-origin axis;
- a full-width tissue delta panel for scGPT and Mouse-Geneformer;
- a separate `table3_model_surface_coverage` output for task-scope definitions.

Visual QA findings after enlarged inspection:

- title, subtitle, legends, panel labels, and footer fit at manuscript width;
- chance-line value labels no longer collide with the dashed line;
- the previous table-like Panel C has been removed from the figure;
- task-scope differences are visible without turning the figure into a table.

Remaining caption note:

- explicitly state that shared 6-task means, text-only zero-shot rows, and
  single-tissue extension rows are separate benchmark task sets.

### Fig6 Organoid Biology Check

The manuscript variant replaces the dense deck-oriented diagnostic figure with:

- a compact count strip for the public organoid dataset footprint;
- a dot-strip panel separating primary prediction metrics from secondary
  flight-response and model-gene checks;
- a top-k enrichment panel for flight-response genes among ranked model genes.

Text changes:

- replaced the visible "surface" label with "dataset";
- expanded abbreviated labels such as "dir." and "rank r" into "direction"
  and "rank corr.";
- shortened the subtitle so it fits at manuscript width without right-edge
  clipping.

Visual QA findings after enlarged inspection:

- title and subtitle fit without clipping;
- Panel A labels are readable and no longer use internal benchmark shorthand;
- Panel B y-axis labels fit after increasing the left margin;
- Panel C keeps the enrichment readout as a figure rather than a table;
- grayscale and color-vision transforms preserve the main interpretation
  because marker shape, row position, bar height, and direct numeric labels
  carry the signal.

Remaining caption note:

- state that this is a small public human neural-organoid extension dataset
  and should be read as a biology check, not a replacement for the primary
  mouse held-out mission benchmark.

## Numeric Verification

Command:

```bash
python3 scripts/audit_premium_visual_sources.py
```

Result:

- 117 checks run.
- 117 checks passed.
- 0 failures.

Interpretation:

- The visual layer remains numerically faithful to current source files.
- This does not replace source freeze or manuscript caption review.

## Current Readiness Verdict

| Figure | Manuscript variant | Enlarged visual QA | Current status |
|---|---|---|---|
| Fig2 pathway rescue | present | PASS with caption note | main-figure candidate |
| Fig3 model comparison | present | PASS with comparison table | main-figure candidate |
| Fig6 organoid biology check | present | PASS after manuscript relayout | extension figure candidate |

Next useful figure work:

1. Draft captions for Fig1/Fig2/Fig3/Fig6.
2. Build a grayscale/colorblind QA contact sheet.
3. Keep Fig4 as a resource overview and Fig5 as release/status material.
4. Decide whether the paper needs a separate table for Fig6 dataset metadata.

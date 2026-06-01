# v1-v9 Figure QA Pass 1

Date: 2026-05-31

Purpose: first pass at checking whether the highest-priority existing figure
assets can be reused for the v1-v9 slide deck.

Companion files:

- `docs/V1_V9_PRESENTATION_AND_MANUSCRIPT_MASTER_OUTLINE_2026_05_31.md`
- `docs/V1_V9_SLIDE_ASSET_MANIFEST_2026_05_31.md`

## Tools Tried

Playwright/browser route:

- `npx` was not available in the shell, so the bundled Playwright CLI wrapper
  could not be used in this pass.
- Node REPL dynamic import of `playwright` failed because the `playwright`
  package was not available.

QuickLook route:

- `qlmanage -t -s 1400` successfully produced thumbnails for the eight
  selected assets.
- Output directory:
  `output/slide_asset_qa/quicklook/`

Interpretation:

- QuickLook is enough for file-existence and basic PNG/PDF thumbnail checks.
- QuickLook is not enough for D3 HTML figures, because the selected HTML files
  depend on external D3 from `https://d3js.org/d3.v7.min.js`.

## Eight-Asset QA Results

| Asset | QuickLook thumbnail | Visual status | Deck action |
|---|---|---|---|
| `figures/fig1_tissue_hierarchy.html` | produced | QuickLook shows blank SVG container plus button, not the plot | Needs real browser render or figure regeneration |
| `figures/fig2_pathway_mechanism.html` | produced | D3 HTML; QuickLook not trusted | Needs real browser render or figure regeneration |
| `figures/fig3_model_comparison.html` | produced | D3 HTML; QuickLook not trusted | Needs real browser render or figure regeneration |
| `v4/figures/html/Fig1_benchmark.html` | produced | D3 HTML; QuickLook not trusted | Needs real browser render or figure regeneration |
| `v2/figures/Fig1_temporal.html` | produced | D3 HTML; QuickLook not trusted | Needs real browser render or figure regeneration |
| `v3/figures/v3_Fig2_spatial_overview.html` | produced | D3 HTML; QuickLook not trusted | Needs real browser render or figure regeneration |
| `v8/figures/Figure1_Species_Transfer.png` | produced | Visually usable; 2681 x 731 PNG, three-panel layout | Candidate for main deck |
| `v8/figures/Figure2_Stressor_Decomposition.png` | produced | Visually readable, but original figure has an empty right panel | Revise or crop before main deck |

## QuickLook Outputs

- `output/slide_asset_qa/quicklook/fig1_tissue_hierarchy.html.png`
- `output/slide_asset_qa/quicklook/fig2_pathway_mechanism.html.png`
- `output/slide_asset_qa/quicklook/fig3_model_comparison.html.png`
- `output/slide_asset_qa/quicklook/Fig1_benchmark.html.png`
- `output/slide_asset_qa/quicklook/Fig1_temporal.html.png`
- `output/slide_asset_qa/quicklook/v3_Fig2_spatial_overview.html.png`
- `output/slide_asset_qa/quicklook/Figure1_Species_Transfer.png.png`
- `output/slide_asset_qa/quicklook/Figure2_Stressor_Decomposition.png.png`

## HTML Figure Dependency Finding

The selected HTML figures use remote D3:

- `figures/fig1_tissue_hierarchy.html`
- `figures/fig2_pathway_mechanism.html`
- `figures/fig3_model_comparison.html`
- `v4/figures/html/Fig1_benchmark.html`
- `v2/figures/Fig1_temporal.html`
- `v3/figures/v3_Fig2_spatial_overview.html`

This matters because slide production should not depend on live CDN access. For
the final deck, use one of these routes:

1. Browser-render each HTML figure in an environment with D3 access and export
   static PNG/SVG.
2. Vendor D3 locally and patch figure HTML for local/offline rendering.
3. Regenerate the most important figures in a unified static style from the
   underlying CSV/JSON result sources.

Route 3 is probably best for a polished paper/deck because it also resolves
style consistency and canonical-value checks.

## v8 PNG Finding

`v8/figures/Figure1_Species_Transfer.png` is a strong candidate for reuse:

- It is already a high-resolution PNG.
- It has a clear three-panel structure.
- It supports the v8 BRIDGE story.
- It still needs a visible "hypothesis-generation only" boundary on the slide.

`v8/figures/Figure2_Stressor_Decomposition.png` needs revision before main-deck
use:

- It is high-resolution and readable.
- The rightmost panel is blank in the original PNG, not just in QuickLook.
- Use either a cropped two-panel version or regenerate the figure with the
  intended third panel populated/removed.

## Immediate Recommendation

Do not build the final deck by screenshotting the current HTML files through
QuickLook.

Instead, make the next block one of:

1. Static regeneration of Figures 1-4 from canonical data sources.
2. Browser-render/export setup with local D3 or Playwright.
3. A hybrid: reuse v8 PNG Figure1, revise/crop v8 Figure2, and regenerate the
   core benchmark figures in a unified style.

Preferred route:

> Regenerate a unified static figure set for the main story, using the existing
> HTML figures as design references and the canonical result files as numeric
> sources.

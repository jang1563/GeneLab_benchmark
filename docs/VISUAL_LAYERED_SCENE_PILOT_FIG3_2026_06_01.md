# Visual Layered Scene Pilot: Fig3 Model Comparison

Date: 2026-06-01

Purpose: test the layered-slide rule on the model-tier comparison without
turning the result into a leaderboard.

## Generated Outputs

- `output/premium_slide_scenes/fig3_model_tier_scene_plate.png`
- `output/premium_slide_scenes/fig3_model_tier_layered_scene.png`
- `output/premium_slide_scenes/fig3_model_tier_layered_scene.pdf`
- `output/premium_slide_scenes/fig3_model_tier_editable_overlay_spec.json`
- `output/premium_slide_scenes/fig3_model_tier_layered_scene_manifest.json`

Generator:

- `scripts/build_layered_slide_scenes.py --scene fig3`

Source proof object:

- `output/premium_figures/manuscript_variants/premium_fig3_model_tier_comparison_manuscript.png`

## Layer Implementation

| Layer | Implementation | Verdict |
|---|---|---|
| Z0 canvas / atmosphere | procedural dark benchmark texture | PASS; aligned with Fig1/Fig2 family |
| Z1 measurement layer | faint model-family lane: classical, single-cell, text LLM | PASS; restrained and non-dashboard-like |
| Z2 evidence surface | source-derived Panel A and Panel B crops with backplates | PASS after crop repair; model names and panel titles remain readable |
| Z3 interpretation layer | headline, short callout, two claim blocks, corner focus windows | PASS for pilot; frames comparability rather than model spectacle |
| Z4 trust/status layer | concise shared-row comparability label | PASS; the key caveat is visible but not dominant |
| Z5 motion/focus layer | short transition arrow from aggregate comparison to tissue deltas | PASS after replacing a long diagonal arrow |

## Render Iterations

First render issues:

- Panel A crop cut off model labels and the right side of the panel title;
- Panel B crop clipped the panel title;
- a long diagonal arrow entered the proof panels and felt too presentation-like;
- the Panel B focus window included too much y-label and legend territory.

Fixes applied:

- widened both proof crops to preserve model names, panel titles, and value
  labels;
- moved the evidence plates slightly left and kept the left interpretation
  column clean;
- replaced the diagonal arrow with a short transition between the two proof
  plates;
- narrowed the Panel B focus window to the plot region carrying negative
  tissue deltas.

## Current Verdict

Fig3 passes as a draft layered-scene candidate.

The slide claim should remain:

> Scale alone does not transfer.

The scientific boundary should remain visible:

> Shared 6-task rows are direct comparisons.

This prevents the slide from becoming a generic model leaderboard. The central
message is that direct shared-row comparison favors the classical baseline, and
that foundation-model gains are local rather than universal.

## Reusable Decisions

Carry forward to Fig6/v9:

- do not let motion arrows cross proof data unless the arrow is itself the
  measured relationship;
- when the source figure already contains a strong title, crop repair matters
  more than adding external labels;
- trust labels can be short if the proof panel contains the technical caveat;
- model and platform slides need especially restrained visual language.

Watch-outs:

- if this becomes a final deck slide, preserve enough width for long model
  names such as `Mouse-Geneformer`;
- do not recolor model families in a way that implies new grouping or ranking
  beyond the audited proof figure;
- if a talk audience is less technical, replace `PCA-LR` in the overlay
  callout with `classical baseline` and leave the model name in the proof panel.

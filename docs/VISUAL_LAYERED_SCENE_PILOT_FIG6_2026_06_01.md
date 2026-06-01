# Visual Layered Scene Pilot: Fig6 Human Organoid Biology Checks

Date: 2026-06-01

Purpose: test the layered-slide rule on the human organoid extension while
preventing overclaim.

## Generated Outputs

- `output/premium_slide_scenes/fig6_organoid_scene_plate.png`
- `output/premium_slide_scenes/fig6_organoid_layered_scene.png`
- `output/premium_slide_scenes/fig6_organoid_layered_scene.pdf`
- `output/premium_slide_scenes/fig6_organoid_editable_overlay_spec.json`
- `output/premium_slide_scenes/fig6_organoid_layered_scene_manifest.json`

Generator:

- `scripts/build_layered_slide_scenes.py --scene fig6`

Source proof object:

- `output/premium_figures/manuscript_variants/premium_fig6_organoid_biology_check_manuscript.png`

## Layer Implementation

| Layer | Implementation | Verdict |
|---|---|---|
| Z0 canvas / atmosphere | procedural dark cellular texture | PASS; consistent with Fig1-Fig3 and not decorative |
| Z1 measurement layer | faint source/sample/gene-check lane | PASS; supports the cautionary scope |
| Z2 evidence surface | source-derived footprint crop plus prediction/biology-check crop | PASS after tightening Panel A crop |
| Z3 interpretation layer | headline, short callout, two claim blocks, corner focus windows | PASS; keeps the slide diagnostic rather than promotional |
| Z4 trust/status layer | draft-extension caution label | PASS; visible enough to prevent validation overclaim |
| Z5 motion/focus layer | short transition from dataset footprint to biology checks | PASS; restrained and non-causal |

## Render Iterations

First render issues:

- the dataset-footprint crop included too much white space below Panel A;
- the combined Panel B/C crop risked cutting the long Panel C title and bottom
  labels;
- the slide needed a stronger caution layer than Fig1-Fig3 because organoid
  results can be misread as validation.

Fixes applied:

- tightened the Panel A crop around the source/sample footprint;
- widened and deepened the Panel B/C crop to preserve the long title and x-axis
  labels;
- used the headline `Organoids add biology checks`, not `Organoids validate`;
- added the visible scope line
  `Draft extension: source factors remain coupled`;
- focused the top crop on `2 sources, 42 samples` and the lower crop on
  gene-response overlap checks.

## Current Verdict

Fig6 passes as a draft layered-scene candidate.

The slide should be used as an extension/diagnostic result, not as a primary
validation claim. The final deck version should keep the caution layer
slide-native and editable.

## Reusable Decisions

Carry forward to v9 platform slides:

- extension slides need an explicit caution layer, even when the proof figure
  looks strong;
- compact source-footprint crops are better than full-figure pastes;
- one short transition from dataset footprint to diagnostic check is enough;
- avoid visual language that makes draft extension datasets look complete.

Watch-outs:

- if this slide is used in a manuscript talk, do not omit the source-factor
  coupling line;
- keep `biology checks` rather than `validation` in visible slide text;
- if space allows in final deck production, a small note can specify that the
  42 samples come from two public neural-organoid RNA-seq sources.

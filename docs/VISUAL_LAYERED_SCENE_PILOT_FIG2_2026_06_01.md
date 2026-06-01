# Visual Layered Scene Pilot: Fig2 Pathway Summaries

Date: 2026-06-01

Purpose: test the layered-slide rule on the pathway-summary result while
keeping the scientific boundary clear:

> Pathway summaries reduce unwanted coupled-label signal, with selected
> detection-task gains.

## Generated Outputs

- `output/premium_slide_scenes/fig2_pathway_scene_plate.png`
- `output/premium_slide_scenes/fig2_pathway_layered_scene.png`
- `output/premium_slide_scenes/fig2_pathway_layered_scene.pdf`
- `output/premium_slide_scenes/fig2_pathway_editable_overlay_spec.json`
- `output/premium_slide_scenes/fig2_pathway_layered_scene_manifest.json`

Generator:

- `scripts/build_layered_slide_scenes.py --scene fig2`

Source proof object:

- `output/premium_figures/manuscript_variants/premium_fig2_pathway_rescue_manuscript.png`

## Layer Implementation

| Layer | Implementation | Verdict |
|---|---|---|
| Z0 canvas / atmosphere | procedural dark transcriptomic texture | PASS; consistent with Fig1 without becoming decorative |
| Z1 measurement layer | faint gene-input to pathway-summary assay lane | PASS after removing top text labels that competed with the headline |
| Z2 evidence surface | source-derived Panel A and Panel B crops with backplates | PASS; proof crops remain readable on the scene plate alone |
| Z3 interpretation layer | headline, short callout, two claim blocks, corner focus windows | PASS for pilot; less jargon and less visually heavy than large rings |
| Z4 trust/status layer | concise diagnostic-scope label | PASS; no raw paths, internal decisions, or workflow language visible |
| Z5 motion/focus layer | one gene-to-pathway transition arrow inside Panel A | PASS; one dominant movement only |

## Render Iterations

First render issues:

- the headline and callout competed with the Panel A proof crop;
- top measurement labels added extra jargon near the headline;
- large oval focus rings looked coarse and over-annotated;
- the Panel B ring extended awkwardly and did not feel slide-native;
- visible wording used `label leakage`, which is accurate but less intuitive
  for a mixed scientific audience.

Fixes applied:

- moved and resized the proof crops so the left interpretation column has
  clean breathing room;
- replaced the top measurement arrow/text with a quieter assay-lane cue;
- changed the main headline to `Pathways suppress unwanted labels`;
- replaced large ellipses with corner focus windows and a single subtle
  transition arrow;
- changed the scope line to
  `Diagnostic check: mission, hardware, gravity labels`;
- narrowed the Panel B focus window to the selected positive-gain region.

## Current Verdict

Fig2 now passes as a draft layered-scene candidate.

It should be treated like Fig1: the current Matplotlib render is a faithful
prototype, but final deck production should keep Z0-Z2 as a raster scene plate
and rebuild Z3-Z5 as editable slide-native text, focus windows, and arrows.

## Reusable Decisions

Carry forward to Fig3/Fig6:

- avoid large oval annotations unless the underlying geometry genuinely calls
  for an ellipse;
- prefer crop windows, bracket corners, and one motion path for quantitative
  proof objects;
- replace ML jargon with reader-facing claims when the proof panel already
  carries the technical wording;
- keep status text about scientific scope, not internal build state.

Watch-outs:

- pathway summaries should not be visually framed as universally better;
- the selected gain region must remain tissue/task-specific;
- if slide-native rebuilding changes typography, recheck that the headline and
  Panel A proof crop do not collide.

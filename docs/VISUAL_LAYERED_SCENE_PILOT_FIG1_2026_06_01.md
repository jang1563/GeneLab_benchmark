# Visual Layered Scene Pilot: Fig1 Tissue Transfer

Date: 2026-06-01

Purpose: record the first practical test of the layered-slide rule:

> Layered PNG scene + editable scientific interpretation.

## Generated Outputs

- `output/premium_slide_scenes/fig1_tissue_transfer_scene_plate.png`
- `output/premium_slide_scenes/fig1_tissue_transfer_layered_scene.png`
- `output/premium_slide_scenes/fig1_tissue_transfer_layered_scene.pdf`
- `output/premium_slide_scenes/fig1_tissue_transfer_editable_overlay_spec.json`
- `output/premium_slide_scenes/fig1_tissue_transfer_layered_scene_manifest.json`

Generator:

- `scripts/build_layered_slide_scenes.py`

Source proof object:

- `output/premium_figures/premium_fig1_tissue_transfer_hierarchy.png`

## Layer Implementation

| Layer | Implementation | Verdict |
|---|---|---|
| Z0 canvas / atmosphere | procedural dark transcriptomic texture | PASS; contextual without becoming a space wallpaper |
| Z1 measurement layer | faint held-out mission arc and train/held-out ticks | PASS; readable but intentionally quiet |
| Z2 evidence surface | cropped, source-derived Fig1 proof object with backplate/shadow | PASS after crop; data region is larger than full-figure paste |
| Z3 interpretation layer | headline, short callout, high/near-chance focus rings | PASS for pilot; should become slide-native editable objects later |
| Z4 trust/status layer | concise scope label | PASS; no raw file paths or internal process wording visible |
| Z5 motion/focus layer | one transfer-break arrow | PASS; one dominant movement only |

## Render Iterations

First render issues:

- proof plate overlapped the left headline;
- focus ring targeted the mid-range region too strongly;
- full proof figure included too much white space, making the data region small.

Fixes applied:

- moved and resized the proof plate;
- replaced the full proof figure with a source-derived plot crop;
- expanded the crop leftward to preserve tissue labels;
- shortened the headline/callout;
- tuned focus rings for high-transfer and near-chance tissue groups.

## Current Verdict

Pilot passes as a proof of design grammar.

It is not yet the final presentation slide because the overlay is still drawn
by Matplotlib. For final deck production, keep the scene plate as raster
Z0-Z2 and rebuild Z3-Z5 as editable slide-native text, rings, and arrows.

## Reusable Decisions

Carry forward to Fig2/Fig3:

- use a source-derived proof crop rather than a full pasted manuscript figure;
- preserve a separate scene plate and editable overlay specification;
- keep only one focus/motion device per slide;
- keep trust text concise and reader-facing;
- avoid dashboard cards and decorative backgrounds.

Watch-outs:

- make sure proof-object crops do not cut axis labels;
- check that overlays do not obscure plotted values;
- do not let dark-scene atmosphere become the story.

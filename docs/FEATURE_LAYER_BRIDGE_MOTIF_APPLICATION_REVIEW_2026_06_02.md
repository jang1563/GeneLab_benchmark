# Feature-Layer Bridge Motif Application Review

Date: 2026-06-02

Purpose:

- apply the v0.2 biological motif pack to an actual methods bridge slide;
- explain gene/RNA versus pathway feature layers for a first-time viewer;
- test whether the motif pack can support a premium, layered slide scene rather
  than only a contact sheet.

## Output

Scene root:

- `output/premium_feature_layer_bridge/`

Generated assets:

- `output/premium_feature_layer_bridge/feature_layer_bridge_scene.png`
- `output/premium_feature_layer_bridge/feature_layer_bridge_scene.pdf`
- `output/premium_feature_layer_bridge/feature_layer_bridge_scene_plate.png`
- `output/premium_feature_layer_bridge/feature_layer_bridge_overlay_spec.json`
- `output/premium_feature_layer_bridge/feature_layer_bridge_manifest.json`
- `output/premium_feature_layer_bridge/feature_layer_bridge_qa.json`

Generator:

- `scripts/build_feature_layer_bridge_scene.py`

## Audience Question

The slide answers:

> What does gene versus pathway mean in this benchmark?

The visual answer is:

> Samples are measured as genes, gene activity is summarized into pathway
> biology, and both feature layers are scored under the same hidden-mission
> split.

## Visual Structure

The scene uses a layered PNG plate plus editable text shell.

Z-stack:

- Z0 canvas: warm scientific paper field with faint grid and grain;
- Z1 measurement layer: horizontal assay rail, mission split markers, and
  sample-by-gene orientation;
- Z2 evidence surface: v0.2 motif plates for sample matrix, RNA activity,
  pathway program, protein-complex context, and score sheet;
- Z3 interpretation layer: gene cloud, compression wedge, pathway bands, and
  directional transfer arrows;
- Z4 trust/status layer: same-split marker, training/hidden mission statement,
  and caveat that pathway summaries do not improve every task.

## Text Simplification

Earlier bridge language was too internal.

Removed or softened:

- "Feature layers turn expression into biology";
- "same fold rule";
- "summary layer";
- "not a new dataset".

Current visible text:

- headline: "Gene signals become pathway summaries";
- support: "Both layers use the same hidden-mission split.";
- primary labels: "Samples x genes", "Gene activity", "Pathway summary",
  "Hidden mission score";
- status: "fit on training missions; score on hidden missions";
- caveat: "Pathways summarize biology; not every task improves."

Visible word count:

- 44 words;
- budget: 48 words.

## Manual Visual QA

Full-size inspection:

- pass: 3840 x 2160 preview inspected;
- pass: no material text overlap;
- pass: headline, section labels, and status line are readable;
- pass: central compression wedge makes gene-to-pathway transformation visible.

Thumbnail inspection:

- pass: main left-to-right flow remains readable at fit-to-screen scale;
- conditional pass: M1-M4 split labels are intentionally secondary and may not
  be legible in small thumbnails.

Biological specificity:

- pass: sample matrix, RNA trace, gene cloud, pathway field, protein-complex
  context, and hidden-score marker make the scene biological rather than a
  generic workflow rail.

Premium quality gate:

- conditional pass: suitable as a draft premium bridge scene for deck assembly;
- limitation: this is still a symbolic explanation scene, not a source-derived
  proof figure.

## Claim Boundary

The slide should not imply that pathway summaries always improve performance.

Allowed claim:

- pathway features are train-only summaries of gene activity;
- both gene-level and pathway-level representations are evaluated with the same
  hidden-mission split;
- pathway summaries add biological interpretability, with task-dependent
  performance effects.

Avoid:

- claiming universal pathway superiority;
- implying a new dataset was introduced;
- treating this bridge as evidence for a result without nearby empirical
  panels.

## Next Polish Options

Keep this version as the deck assembly candidate.

Possible later upgrades:

- replace the symbolic score sheet with a small source-derived performance
  proof tile when final figures are locked;
- build a matching grayscale/colorblind QA variant;
- create a manuscript-safe version with fewer atmospheric layers;
- convert labels into native PowerPoint text when assembling the deck.

## Verdict

Draft premium bridge candidate.

This version is materially stronger than the earlier bridge direction because
it explains the method visually, uses biology-specific motif language, keeps
visible jargon low, and preserves a clear claim boundary.

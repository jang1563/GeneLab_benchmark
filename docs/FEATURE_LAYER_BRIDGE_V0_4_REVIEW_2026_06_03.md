# Feature-layer bridge v0.4 review

Date: 2026-06-03

## Anchor

This continues the v1-v9 deck-spine pass. The active slide is slide 6:

> Feature layers: genes versus pathway summaries before results.

The purpose is to help first-time viewers understand what changes between a
gene-level view and a pathway-summary view before the deck moves into v1-v7
benchmark result slides.

## Why v0.4 Was Needed

Existing drafts were useful but not ready for the 2026 deck spine:

- `output/premium_feature_layer_bridge/feature_layer_bridge_scene.png`
  explained the flow, but it was too light/flat relative to the current dark
  premium deck identity and used internal labels such as `M1`-style mission
  markers.
- `output/premium_feature_layer_bridge_v0_3_test/feature_layer_bridge_v0_3_motif_test.png`
  improved biological/data texture, but still carried a meta-test title and did
  not have enough contrast for a final slide candidate.

v0.4 combines the better parts:

- dark mission-grade canvas,
- larger matrix and pathway proof-like objects,
- simpler first-time-viewer text,
- scene plate separated from editable overlay,
- explicit caveat that pathways help some tasks, not all.

## Generated Outputs

Builder:

- `scripts/build_feature_layer_bridge_v4.py`

Output root:

- `output/premium_feature_layer_bridge_v0_4/`

Generated slide assets:

- `output/premium_feature_layer_bridge_v0_4/feature_layer_bridge_v0_4_scene_plate.png`
- `output/premium_feature_layer_bridge_v0_4/feature_layer_bridge_v0_4_rendered_preview.png`

Generated QA and specs:

- `output/premium_feature_layer_bridge_v0_4/qa/feature_layer_bridge_v0_4_rendered_preview_grayscale.png`
- `output/premium_feature_layer_bridge_v0_4/feature_layer_bridge_v0_4_overlay_spec.json`
- `output/premium_feature_layer_bridge_v0_4/feature_layer_bridge_v0_4_manifest.json`
- `output/premium_feature_layer_bridge_v0_4/feature_layer_bridge_v0_4_qa.json`

## Visible Text Policy

Audience-facing text is intentionally simple:

- "What the model sees can change"
- "Gene table"
- "Gene-level view"
- "Pathway summaries"
- "Score"
- "held-out mission"

Removed from the visible layer:

- `M1`, `M2`, `M3`, `M4`
- "same split"
- "replacement test"
- "generated schematic texture"
- "not source evidence"

Current visible text count is 39 words against a 45-word method-slide budget.

## Visual QA

Verdict: draft premium slide-6 candidate.

Passes:

- rendered preview inspected at full size;
- scene plate inspected separately and is suitable for editable overlay use;
- grayscale QA preserves the object sequence and arrow path;
- matrix and pathway objects are large enough after crop/scale adjustment;
- no visible text overlap;
- no result claim is introduced.

Remaining caution:

- This is a schematic methods bridge, not a quantitative result figure.
- The final deck should rebuild the headline, labels, status, and caveat as
  editable slide objects instead of relying only on the rendered PNG.
- The footer caveat is intentionally subordinate; it should remain visible but
  not compete with the central source-to-score flow.

## Deck Placement

Use as slide 6 in the 24-slide spine:

1. B1/B2 compressed methods bridge.
2. B3/B4 compressed held-out/train-only bridge.
3. v0.4 feature-layer bridge.
4. v1-v7 result core begins with tissue hierarchy and pathway rescue.

This placement prevents the pathway result slide from appearing before the
audience understands what a gene-level versus pathway-summary view means.

## Claim Boundary

This slide supports:

- methods explanation,
- feature-view intuition,
- train/held-out split continuity,
- transition into pathway result figures.

It does not support:

- a new benchmark result,
- a universal pathway-superiority claim,
- a source-specific proof claim,
- final manuscript quantitative evidence.

## Next Work

With slide 6 now represented, the next core-deck gap is slides 10-14:

- v4 hardening,
- temporal and RRRM lessons,
- negative results,
- biological interpretation,
- human translation.

Those slides need static export or premium rebuild planning so the v1-v7 result
core remains visually stronger than the v9 extension materials.

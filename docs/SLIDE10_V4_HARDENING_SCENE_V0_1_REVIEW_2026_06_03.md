# Slide 10 v4 hardening scene v0.1 review

Date: 2026-06-03

## Anchor

This continues the v1-v9 deck-spine visual pass after the slide 6 feature-layer
bridge.

The active slide is slide 10:

> v4 hardening: the benchmark result survives a wider method grid.

The goal is to make the first core result slide feel like high-profile
scientific evidence rather than a pasted HTML figure or a table in disguise.

## Why This Rebuild Was Needed

The old v4 HTML figures remain useful source material, but they are too dense
and too visually inconsistent for the premium deck spine.

This rebuild makes slide 10 do one job:

- show the canonical tissue signal surface;
- show that v4 expanded the benchmark across tissues, classifiers, and feature
  views;
- keep the significance boundary visible;
- avoid turning the slide into a new leaderboard.

## Generated Outputs

Builder:

- `scripts/build_slide10_v4_hardening_scene.py`

Output root:

- `output/premium_core_result_slides/slide10_v4_hardening_v0_1/`

Generated slide assets:

- `output/premium_core_result_slides/slide10_v4_hardening_v0_1/slide10_v4_hardening_scene_plate.png`
- `output/premium_core_result_slides/slide10_v4_hardening_v0_1/slide10_v4_hardening_rendered_preview.png`

Generated QA and metadata:

- `output/premium_core_result_slides/slide10_v4_hardening_v0_1/qa/slide10_v4_hardening_rendered_preview_grayscale.png`
- `output/premium_core_result_slides/slide10_v4_hardening_v0_1/slide10_v4_hardening_manifest.json`
- `output/premium_core_result_slides/slide10_v4_hardening_v0_1/slide10_v4_hardening_qa.json`

## Source And Claim Policy

Visible values use the canonical public-facing surface:

- `docs/CANONICAL_RESULTS_V7_1.md`

The manifest also records the related v4 evaluation source:

- `v4/evaluation/M1_summary.json`

Visible slide claims:

- v4 spans 8 tissues, 8 classifiers, and 4 feature views.
- This yields 256 evaluations.
- 40 evaluations are significant.
- 6/8 tissues have signal.
- The displayed best-row table is canonical v7.1 and not every best row is
  significant.

Claims this slide should not make:

- every displayed best tissue row is statistically significant;
- this is a new leaderboard;
- v4 proves a universal tissue or pathway rule;
- raw/canonical discrepancies are resolved on the visual layer.

## Visual QA

Verdict: draft premium slide-10 candidate.

Passes:

- rendered preview inspected at full size;
- scene plate inspected separately for future editable overlay use;
- grayscale QA preserves the tissue bars, hardening grid, and metric tiles;
- primary read order is clear: title, tissue surface, hardening grid, metric
  tiles, caveat;
- visible text count is 44 words against a 45-word result-slide budget;
- no visible text overlap detected;
- dark deck identity matches the slide 6 v0.4 bridge direction.

Remaining caution:

- small method and feature labels are supporting texture, not the main
  explanation layer;
- final PPTX implementation should rebuild title, labels, caveat, and source as
  editable objects over the scene plate;
- if slide 10 becomes too dense in the final deck, keep the 256/40/6-of-8 tiles
  and simplify the per-row method labels first.

## Deck Placement

Use after the feature-layer bridge and before temporal/RRRM lessons:

1. Feature-layer bridge explains what changes between gene and pathway views.
2. Slide 10 shows the result surface survives broader v4 hardening.
3. Slide 11 can then explain temporal, preservation, and RRRM lessons.

This keeps the deck from jumping straight from methods into biology without
showing robustness evidence.

## Next Work

The next P0 visual target is slide 14:

- human translation ladder;
- pathway-level partial conservation;
- weak TF/gene-level transfer boundary;
- target-evidence tiers without clinical overclaiming.

Slides 11-13 remain P1 rebuild targets after slide 14 is stabilized.

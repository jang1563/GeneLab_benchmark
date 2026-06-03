# Slide 12 negative boundary v0.1 review

Date: 2026-06-03

## Anchor

This follows slide 11 temporal/RRRM guardrails in the v1-v9 core result rebuild.

The active slide is slide 12:

> Negative and failure-mode results define benchmark limits.

The purpose is to make null/negative findings feel intentional and useful,
without implying that the underlying biology is absent or that all future
models must fail.

## Why This Rebuild Was Needed

The old v3 HTML figure and README contain valuable negative evidence, but the
deck needs a clearer story:

- spatial brain data do not show a detectable spaceflight signal in the analyzed
  OSD-352 setting;
- companion bulk and spatial results both warn against over-reading the small
  brain cohort;
- UCE and scFoundation do not automatically improve the benchmark surface;
- these negatives define current task boundaries and motivate better follow-up
  tests.

## Generated Outputs

Builder:

- `scripts/build_slide12_negative_boundary_scene.py`

Output root:

- `output/premium_core_result_slides/slide12_negative_boundary_v0_1/`

Generated slide assets:

- `output/premium_core_result_slides/slide12_negative_boundary_v0_1/slide12_negative_boundary_scene_plate.png`
- `output/premium_core_result_slides/slide12_negative_boundary_v0_1/slide12_negative_boundary_rendered_preview.png`

Generated QA and metadata:

- `output/premium_core_result_slides/slide12_negative_boundary_v0_1/qa/slide12_negative_boundary_rendered_preview_grayscale.png`
- `output/premium_core_result_slides/slide12_negative_boundary_v0_1/slide12_negative_boundary_source_summary.json`
- `output/premium_core_result_slides/slide12_negative_boundary_v0_1/slide12_negative_boundary_manifest.json`
- `output/premium_core_result_slides/slide12_negative_boundary_v0_1/slide12_negative_boundary_qa.json`

## Source And Claim Policy

Source files:

- `v3/evaluation/F3a_visium_classification.json`
- `v3/evaluation/F3b_visium_svg.json`
- `v3/evaluation/F3d_visium_cross_resolution.json`
- `v3/evaluation/FM_uce.json`
- `v3/evaluation/FM_scfoundation.json`
- `v3/README.md`

Visible slide claims:

- OSD-352 RR-3 brain spatial section AUROC is 0.139.
- Animal-level spatial AUROC is 0.444.
- Companion bulk AUROC is 0.000 in the same small n=6 setting.
- UCE multi-mission mean AUROC is 0.542 with one significant tissue.
- scFoundation multi-mission mean AUROC is 0.584 with two significant tissues.
- The PCA-LR reference mean from v3 documentation is 0.758.

Claims this slide should not make:

- brain has no spaceflight biology in general;
- all spatial transcriptomics tasks will fail;
- foundation models universally fail;
- PCA-LR is the final best possible model;
- the negative result is a weak or disposable finding.

## Visual QA

Verdict: draft premium slide-12 candidate.

Passes:

- rendered preview inspected at full size;
- grayscale QA preserves the spatial null, boundary corridor, and FM comparison;
- visible text count is 43 words against a 50-word budget;
- first render was revised to shorten the subtitle and reduce an overly heavy
  center corridor;
- left and right evidence surfaces are balanced and do not read as a dense
  table;
- footer caveat explicitly prevents universal absence or impossibility claims.

Remaining caution:

- the central boundary corridor is interpretive, not source evidence;
- final PPTX should keep the scene plate and rebuild headline, source, and
  caveat as editable text;
- if the deck feels too negative in sequence, add a presenter bridge into slide
  13: "boundaries make the biology interpretation more credible."

## Deck Placement

Use after slide 11 and before slide 13:

1. Slide 11: timepoint and RRRM labels need guardrails.
2. Slide 12: negative results define limits.
3. Slide 13: biological interpretation can proceed with clearer boundaries.

This placement prevents the biology slide from looking like a universal
interpretation engine.

## Next Work

The next rebuild target is slide 13:

- biological interpretation triad;
- immune/signaling, metabolism/target, and consensus marker evidence;
- no causal proof or treatment recommendation.

# Slide 11 temporal/RRRM guardrails v0.1 review

Date: 2026-06-03

## Anchor

This follows the slide 10 and slide 14 P0 rebuilds in the v1-v9 deck-spine
visual pass.

The active slide is slide 11:

> Temporal and RRRM lessons: timepoint labels need guardrails.

The purpose is to help first-time viewers understand why v2 is not just another
result panel. It protects the benchmark story from over-reading timepoint,
preservation, and single-cell pilot labels as direct biology.

## Why This Rebuild Was Needed

The v2 results contain many useful details, but a slide deck cannot present the
full temporal, recovery, age, human PBMC, and RRRM analysis surface at once.

This rebuild makes slide 11 do one job:

- show that preservation can dominate ISS-T versus LAR labels;
- show that recovery patterns are descriptive projections;
- show that RRRM single-cell data are benchmark-ready but underpowered for
  mature composition claims;
- keep v2 positioned as a guardrail layer before biological interpretation.

## Generated Outputs

Builder:

- `scripts/build_slide11_temporal_rrrm_guardrails_scene.py`

Output root:

- `output/premium_core_result_slides/slide11_temporal_rrrm_guardrails_v0_1/`

Generated slide assets:

- `output/premium_core_result_slides/slide11_temporal_rrrm_guardrails_v0_1/slide11_temporal_rrrm_guardrails_scene_plate.png`
- `output/premium_core_result_slides/slide11_temporal_rrrm_guardrails_v0_1/slide11_temporal_rrrm_guardrails_rendered_preview.png`

Generated QA and metadata:

- `output/premium_core_result_slides/slide11_temporal_rrrm_guardrails_v0_1/qa/slide11_temporal_rrrm_guardrails_rendered_preview_grayscale.png`
- `output/premium_core_result_slides/slide11_temporal_rrrm_guardrails_v0_1/slide11_temporal_rrrm_guardrails_source_summary.json`
- `output/premium_core_result_slides/slide11_temporal_rrrm_guardrails_v0_1/slide11_temporal_rrrm_guardrails_manifest.json`
- `output/premium_core_result_slides/slide11_temporal_rrrm_guardrails_v0_1/slide11_temporal_rrrm_guardrails_qa.json`

## Source And Claim Policy

Source files:

- `v2/evaluation/T_temporal_summary.json`
- `v2/evaluation/V2_RESULTS_SUMMARY.md`
- `v2/evaluation/F2A_composition.json`
- `v2/docs/RRRM1_BENCHMARK_READY_MANIFEST_2026-03-12.csv`
- `v2/docs/RRRM1_BROAD_ANNOTATION_SUMMARY_2026-03-12.md`

Visible slide claims:

- RR-8 ISS-T/LAR separation is stronger in ground controls than flight:
  GC AUROC 0.973 versus FLT AUROC 0.930.
- RR-8 baseline separation remains high: BSL AUROC 0.983.
- Recovery is descriptive: RR-8 PCA recovery ratio is 0.652 and flight
  probability shifts from 1.000 to 0.404.
- RR-6 also projects closer to baseline: ratio 0.842 and probability
  1.000 to 0.185.
- RRRM has four benchmark-ready tissues and 38,081 cells in the summarized
  post-hardening inventory.
- RRRM composition testing is underpowered: 0 FDR<0.05 composition hits.

Claims this slide should not make:

- ISS-T versus LAR is a clean biological timing label;
- LAR recovery is held-out validation evidence;
- RRRM composition shifts prove mature cell-type mechanism;
- first-pass RRRM annotation is final manuscript-level cell labeling.

## Visual QA

Verdict: draft premium slide-11 candidate.

Passes:

- rendered preview inspected at full size;
- grayscale QA preserves the three major lanes;
- visible text count is 49 words against a 52-word budget;
- subtitle shortened after first render to avoid an over-explained slide;
- the slide reads as confound -> projection -> RRRM readiness, not as a table;
- caveat remains visible but subordinate.

Remaining caution:

- left and right evidence surfaces are still framed plates; final PPTX should
  keep them as background proof shelves and rebuild editable text on top;
- do not add the age-amplification result to this slide unless slide 11 is split
  into two slides;
- do not mix the human PBMC temporal result here unless slide 13 or a
  supplementary bridge needs it.

## Deck Placement

Use after slide 10 and before negative/biology interpretation:

1. Slide 10 establishes v4 result hardening.
2. Slide 11 explains why temporal and single-cell labels need guardrails.
3. Slide 12 can then make negative/failure-mode results feel intentional.

This prevents the audience from interpreting later biological slides as if all
timepoint and cell-type labels were clean primary biology.

## Next Work

The next rebuild target should be slide 12:

- spatial brain negative result;
- foundation-model negative/control result;
- one clear boundary message: negative results define benchmark limits.

Slide 13 biological interpretation can follow once the failure-mode boundary is
visually stable.

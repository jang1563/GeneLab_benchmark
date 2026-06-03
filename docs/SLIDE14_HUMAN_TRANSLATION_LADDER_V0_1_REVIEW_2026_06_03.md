# Slide 14 human translation ladder v0.1 review

Date: 2026-06-03

## Anchor

This follows the slide 10 v4 hardening scene in the v1-v9 deck-spine rebuild.

The active slide is slide 14:

> Human translation: partial pathway and target-tier alignment, not direct
> gene-level transfer.

The purpose is to keep the translational story useful without overstating what
mouse spaceflight signatures can prove in human cfRNA.

## Why This Rebuild Was Needed

The existing v6 Fig10 is a useful six-panel analysis figure, but it is too dense
for the premium deck spine. A first-time viewer needs the translation boundary
before they need the full per-phase detail.

This rebuild makes slide 14 do one job:

- show that pathway-level evidence partially carries across species/data type;
- show that direct gene-level and TF-level transfer remain weak;
- show target evidence tiers as triage, not treatment guidance;
- keep human validation and clinical actionability clearly out of scope.

## Generated Outputs

Builder:

- `scripts/build_slide14_human_translation_ladder_scene.py`

Output root:

- `output/premium_core_result_slides/slide14_human_translation_ladder_v0_1/`

Generated slide assets:

- `output/premium_core_result_slides/slide14_human_translation_ladder_v0_1/slide14_human_translation_ladder_scene_plate.png`
- `output/premium_core_result_slides/slide14_human_translation_ladder_v0_1/slide14_human_translation_ladder_rendered_preview.png`

Generated QA and metadata:

- `output/premium_core_result_slides/slide14_human_translation_ladder_v0_1/qa/slide14_human_translation_ladder_rendered_preview_grayscale.png`
- `output/premium_core_result_slides/slide14_human_translation_ladder_v0_1/slide14_human_translation_ladder_source_summary.json`
- `output/premium_core_result_slides/slide14_human_translation_ladder_v0_1/slide14_human_translation_ladder_manifest.json`
- `output/premium_core_result_slides/slide14_human_translation_ladder_v0_1/slide14_human_translation_ladder_qa.json`

## Source And Claim Policy

Source files:

- `v6/evaluation/V6_A_gene_conservation.json`
- `v6/evaluation/V6_B_pathway_conservation.json`
- `v6/evaluation/V6_C_cross_species_transfer.json`
- `v6/evaluation/V6_D_biomarker_validation.json`
- `v6/evaluation/V6_E_tf_conservation.json`
- `v6/evaluation/V6_F_drug_target_validation.json`
- `docs/PROJECT_RESULTS_LOCATION_INVENTORY_2026_05_31.md`

Visible slide claims:

- translation is partial, not direct;
- pathway alignment is partial: mean rho 0.285 across five analyzed tissues;
- direct gene-transfer tests are weak: 0/8 significant transfer tests;
- TF conservation is weak overall: mean rho 0.0298;
- biomarker panel detection is incomplete evidence: 15/20 detected, 0 FDR hits,
  0 panel genes in the human response list;
- target tiers are evidence strata: A=3, B=7, C=186, D=11.

Claims this slide should not make:

- mouse signatures replace human validation;
- gene-level transfer is clean or reliable;
- TF conservation is broadly strong;
- Tier A/B target links are clinical validation;
- any treatment, countermeasure, or drug recommendation follows from v6.

## Visual QA

Verdict: draft premium slide-14 candidate after corridor revision.

Passes:

- rendered preview inspected at full size;
- grayscale QA preserves the main reading order and target-tier structure;
- visible text count is 45 words against a 50-word translation-slide budget;
- title and subtitle carry the plain-language interpretation before technical
  labels appear;
- the central outer container was removed after v0.1 first render because it
  looked too much like nested card boxes;
- biological symbols are used as small signposts, not as hero proof objects;
- source and caveat text remain visible but subordinate.

Remaining caution:

- the slide still uses evidence shelves; final PPTX should keep those as
  background surfaces and rebuild title/callouts as editable objects;
- `rho`, `AUROC`, and tier labels are acceptable evidence labels, but they
  should not be the first thing a presenter says aloud;
- if the final deck feels too busy, simplify the left mouse evidence checklist
  before removing the central translation ladder.

## Deck Placement

Use after the v4 hardening and intermediate biology slides:

1. Slide 10: v4 hardening protects against cherry-pick criticism.
2. Slides 11-13: temporal/RRRM, negative results, and biological interpretation.
3. Slide 14: human translation boundary closes the v1-v7 result core.

This keeps translation framed as a boundary-aware extension of the benchmark,
not as a clinical claim.

## Next Work

With P0 slide 10 and slide 14 now represented, the next visual targets are P1:

- slide 11 temporal and RRRM lessons;
- slide 12 negative and failure-mode results;
- slide 13 biological interpretation triad.

The next rebuild should start with slide 11 because it provides the bridge
between benchmark robustness and later biological interpretation.

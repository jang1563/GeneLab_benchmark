# Slides 4-5 Dark Methods Variant Review

Date: 2026-06-03

## Scope

This review records the dark methods variants for slides 4-5. The goal is to remove the visual discontinuity between the dark opening block, the feature-layer bridge, and the dark result-family slides while preserving the existing light methods scenes as backup/source-proof assets.

## Why This Was Needed

After slides 1-3 were upgraded into dark premium opening scenes, slides 4-5 stood out as light proof-stage surfaces. That contrast was defensible as a methods-proof rhythm, but it weakened the continuity of the first 14-slide production spine.

The new decision is to use dark methods variants as the main-deck candidates and retain the light versions for speaker notes, backup, or source-proof detail.

## Outputs

- Contact sheet: `output/premium_methods_dark_variants_slides_4_5_v0_1/slides_4_5_dark_methods_contact_sheet.png`
- Grayscale contact sheet: `output/premium_methods_dark_variants_slides_4_5_v0_1/qa/slides_4_5_dark_methods_contact_sheet_grayscale.png`
- Caption pack: `output/premium_methods_dark_variants_slides_4_5_v0_1/slides_4_5_dark_methods_caption_pack.md`
- Overlay spec: `output/premium_methods_dark_variants_slides_4_5_v0_1/slides_4_5_dark_methods_overlay_spec.json`
- Manifest: `output/premium_methods_dark_variants_slides_4_5_v0_1/slides_4_5_dark_methods_manifest.json`
- QA: `output/premium_methods_dark_variants_slides_4_5_v0_1/slides_4_5_dark_methods_qa.json`
- Updated assembly board: `output/slides_1_14_deck_assembly_bridge_v0_1/slides_1_14_deck_assembly_bridge.png`

## Slide Decisions

Slide 4 shows one dark evidence surface: public studies, attached source context, benchmark task, and score record. It answers a first-time viewer question: what exactly is being tested?

Slide 5 shows the training side, split boundary, hidden test side, and train-only guard. It answers the next question: does the model transfer to a mission it did not learn from?

The slide 5 first render used large green/red zones that looked too crude for the deck. Those were revised into subtle dark zones with a clearer boundary line, mission nodes, and a compact train-only guard strip.

## Backup Asset Policy

The original light scenes are not discarded:

- Slide 4 backup: `output/premium_bridge_scenes/b1_evaluation_layer/rendered_preview.png`
- Slide 4 notes/backup: `output/premium_bridge_rebuild_scenes/b2_study_to_task_premium/rendered_preview.png`
- Slide 5 backup: `output/premium_bridge_rebuild_scenes/b3_mission_held_out_premium/rendered_preview.png`
- Slide 5 notes/backup: `output/premium_bridge_rebuild_scenes/b4_train_only_guard_premium/rendered_preview.png`

The main deck should use the dark variants. The light versions remain useful where a paper-like proof-stage contrast or more explicit methods definition is needed.

## Z-Stack Rule

- Z0: dark field and molecular texture.
- Z1: grid, orbit arcs, task lanes, and split boundary.
- Z2: source/task/score plates and mission nodes.
- Z3: editable labels, arrows, and focus text.
- Z4: claim boundary and source/caveat line.

The rendered previews should not be treated as final all-in-one PPTX slides. The final deck should rebuild all text as editable objects.

## QA Findings

Automatic QA passes:

- 2 rendered previews exist.
- 2 scene plates exist.
- 2 grayscale QA outputs exist.
- Contact sheet and grayscale contact sheet exist.
- All previews are 3840 x 2160.
- Nonblank proxy passes for both slides.
- Forbidden visible terms are absent.

Manual visual review was performed on the color contact sheet, grayscale contact sheet, individual slide 4 and slide 5 previews, and the updated 1-14 assembly board.

Current status: pass as draft main-deck methods candidates.

## Remaining Watchpoints

Slide 4 should stay a first-time-viewer bridge, not an implementation diagram. B2 definitions belong in notes or backup.

Slide 5 should never be described as random-sample validation. The independence unit is the mission.

The final PPTX should keep the dark variants as scene plates and rebuild the title, subtitle, bridge, caveat, source, and small labels as editable text.

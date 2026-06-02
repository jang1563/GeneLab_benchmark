# BioVis source-proof module pack v0.4 review

Date: 2026-06-02

## Purpose

v0.4 moves the biology visual system from decorative biological motifs toward
evidence-ready figure infrastructure. The goal is not to add more icons. The
goal is to reserve a rigorous visual home for source objects before final
figures are available.

Design rule:

> Layered PNG scene + editable scientific interpretation, with proof object,
> status, source line, and caveat reserved before claims are written.

This prevents a common failure mode in scientific slide decks: a polished claim
is composed first, and the evidence screenshot, dataset accession, or result
plot is retrofitted later.

## Generated Assets

Asset root:

- `assets/biovis_source_proof_module_pack_v0_4/`

Generator:

- `scripts/build_biovis_source_proof_module_pack_v4.py`

The pack contains 6 base modules, each in light and dark variants:

- `source_dataset_record_plate`
- `expression_matrix_proof_plate`
- `single_cell_embedding_proof_plate`
- `held_out_task_proof_plate`
- `organoid_evidence_plate`
- `result_claim_plate`

Generated editable and rendered assets:

- `assets/biovis_source_proof_module_pack_v0_4/modules/svg/`
- `assets/biovis_source_proof_module_pack_v0_4/modules/png/`
- `assets/biovis_source_proof_module_pack_v0_4/modules_dark/svg/`
- `assets/biovis_source_proof_module_pack_v0_4/modules_dark/png/`

Preview and QA assets:

- `assets/biovis_source_proof_module_pack_v0_4/preview/biovis_source_proof_modules_v0_4_light.png`
- `assets/biovis_source_proof_module_pack_v0_4/preview/biovis_source_proof_modules_v0_4_dark.png`
- `assets/biovis_source_proof_module_pack_v0_4/qa/biovis_source_proof_modules_v0_4_light_grayscale.png`
- `assets/biovis_source_proof_module_pack_v0_4/qa/biovis_source_proof_modules_v0_4_dark_grayscale.png`
- `assets/biovis_source_proof_module_pack_v0_4/manifest.json`
- `assets/biovis_source_proof_module_pack_v0_4/manifest.csv`
- `assets/biovis_source_proof_module_pack_v0_4/qa/biovis_source_proof_modules_v0_4_qa.json`

## Bridge Application Stress Test

Bridge application generator:

- `scripts/build_biovis_source_proof_bridge_application_v4.py`

Output root:

- `output/biovis_source_proof_bridge_application_v0_4/`

Generated slide-scale stress tests:

- `output/biovis_source_proof_bridge_application_v0_4/01_light_source_proof_bridge.png`
- `output/biovis_source_proof_bridge_application_v0_4/02_dark_source_proof_bridge.png`
- `output/biovis_source_proof_bridge_application_v0_4/03_source_proof_bridge_contact_sheet.png`
- `output/biovis_source_proof_bridge_application_v0_4/biovis_source_proof_bridge_application_v0_4_manifest.json`
- `output/biovis_source_proof_bridge_application_v0_4/biovis_source_proof_bridge_application_v0_4_qa.json`

The stress test combines v0.4 source-proof modules with the existing v0.3
biological symbol/method modules. It tests whether the source-proof objects can
live inside a layered bridge layout without becoming a card-box collage.

## Visual Review

Verdict: pass with conditions.

What works:

- The proof object slot is explicit in every module.
- Status labels stay attached to the evidence surface.
- Source, processed, held-out, caveat, and validation states remain visible.
- The dark bridge variant has the strongest premium direction and best depth.
- The modules can accept real screenshots or result plots without redesign.

Conditions:

- Main slides should use one dominant proof module plus one method or
  claim-boundary module.
- Multiple large proof modules on one slide should be reserved for appendix,
  review-board, or QA material.
- Light variants are clear but read more like review boards; dark variants are
  stronger for high-profile presentation slides.

## Claim Boundary

The placeholders are not evidence. They are frames for evidence.

Before any empirical claim is made, each placeholder must be replaced with a
real proof object:

- OSDR, GeneLab, GEO, or repository record screenshot
- actual expression matrix or QC heatmap
- actual single-cell embedding, spatial map, or cell-state distribution
- actual train/test or mission-held-out split manifest
- cited organoid image, assay output, or source-derived panel
- final result plot with reviewed claim status

The source line and caveat/status label must remain visible after replacement.

## Production Use Rule

For a main-deck scientific slide:

1. Insert one proof module as the dominant evidence surface.
2. Replace the proof slot with the real source/result object.
3. Attach one v0.3 method module or claim-boundary module.
4. Add only the minimal editable interpretation text needed to explain the
   transition.
5. Keep source/date/status visible.

For appendix or reviewer backup:

1. Use multiple proof modules on one review board.
2. Keep accession, source date, processing state, and caveat visible.
3. Treat the board as evidence inventory, not as a narrative main figure.

## Next P0

Source-object replacement sprint:

1. Capture real OSDR/GeneLab/GEO dataset record proof objects.
2. Capture or generate actual matrix/QC, embedding, split-manifest, and result
   proof objects from the project outputs.
3. Build one production-grade dark slide using one source-proof module plus one
   method/claim module.
4. Run grayscale/colorblind and small-screen readability QA.
5. Only then decide whether the slide is main-figure, bridge, or appendix.

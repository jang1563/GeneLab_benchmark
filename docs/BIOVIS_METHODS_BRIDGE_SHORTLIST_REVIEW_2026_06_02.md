# BioVis methods-bridge shortlist review

Date: 2026-06-02

## Anchor

This implements the drift-control decision from
`docs/VISUAL_METHODS_BRIDGE_PURPOSE_DRIFT_CHECK_2026_06_02.md`.

The goal is not to create more source-proof variants. The goal is to classify
the existing bridge assets before moving back into v1-v9 deck assembly.

## Generated Outputs

Builder:

- `scripts/build_biovis_methods_bridge_shortlist.py`

Output root:

- `output/biovis_methods_bridge_shortlist_v0_1/`

Generated contact sheet:

- `output/biovis_methods_bridge_shortlist_v0_1/contact_sheet/methods_bridge_shortlist_contact_sheet.png`

Generated QA:

- `output/biovis_methods_bridge_shortlist_v0_1/qa/methods_bridge_shortlist_contact_sheet_grayscale.png`
- `output/biovis_methods_bridge_shortlist_v0_1/methods_bridge_shortlist_manifest.json`
- `output/biovis_methods_bridge_shortlist_v0_1/methods_bridge_slide_classification.json`
- `output/biovis_methods_bridge_shortlist_v0_1/methods_bridge_slide_classification.csv`

## Classification Decision

Main-deck candidates:

- Human organoid source-to-matrix:
  `output/biovis_organoid_audience_matrix_proof_v0_2/panels/01_dark_organoid_clean_source_to_matrix.png`
- OSD-120 source-to-task:
  `output/biovis_osd120_audience_split_proof_v0_2/panels/01_dark_osd120_clean_source_to_task.png`

Appendix / reviewer backup:

- OSDR source-record contact sheet:
  `output/biovis_osdr_source_record_proof_panel_v0_1/panels/01_osdr_source_record_contact_sheet.png`
- Extension source lanes:
  `output/biovis_osdr_source_record_proof_panel_v0_1/panels/04_dark_extension_source_lanes.png`

Superseded for audience-facing use:

- Earlier organoid source proof:
  `output/biovis_osdr_source_record_proof_panel_v0_1/panels/02_dark_organoid_source_record_proof.png`
- Earlier OSD-120 source proof:
  `output/biovis_osdr_source_record_proof_panel_v0_1/panels/03_dark_osd120_source_record_proof.png`

## Slide Placement Rule

Use at most one organoid bridge and at most one OSD-120 / multispecies bridge.

The v1-v7 benchmark result core should remain the spine of the talk or
manuscript narrative. Organoid and OSD-120 visuals should explain provenance and
extension readiness; they should not become the main result story.

Recommended placement:

- Short talk: combine as one extension-methods slide or move both to backup.
- Full talk: use the organoid bridge as one methods-extension slide and the
  OSD-120 bridge as one constrained-task slide.
- Manuscript: keep both as supplementary methods/provenance figures unless the
  v9 extension lane becomes a formal result section.

## Visual QA

Verdict: ready as a decision/contact sheet.

Checks:

- The contact sheet is clearly labeled as a review artifact, not a final slide.
- Main-deck, appendix, and superseded roles are visually distinct.
- Real source-proof PNGs are used as thumbnails, so the sheet remains tied to
  existing evidence rather than becoming a speculative diagram.
- Grayscale QA preserves the main hierarchy and role grouping.
- The visible text avoids the internal terms that were removed from the
  audience-facing proof panels.

Limit:

- This is intentionally text-heavy because it is a selection board. It should
  not be used as a high-profile main figure. The final deck should rebuild only
  the selected Z3-Z5 callouts as editable slide objects on top of the selected
  Z2 evidence scene.

## Drift Control

Stop open-ended proof generation here.

Allowed next work:

- deck-spine placement,
- caption pack for selected bridge assets,
- slide-level QA after insertion,
- targeted visual fixes only if insertion QA reveals overlap, low contrast, or
  confusing wording.

Not allowed without a new decision:

- more OSD-120 proof variants,
- more organoid proof variants,
- new extension tracks that push v1-v7 results out of the first narrative half.

## Next Step

Return to the v1-v9 deck spine and decide where the selected methods-bridge
slides sit relative to:

1. v1-v7 benchmark result core,
2. v8 hypothesis-only incubator,
3. v9 platform/provenance layer,
4. explicitly marked extension lanes.

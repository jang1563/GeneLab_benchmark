# Slides 1-14 PPTX Skeleton Build Review

Date: 2026-06-03

## Scope

This pass converts the audited slides 1-14 visual system into an actual editable PPTX skeleton. It is a production draft, not the final polished slide deck.

The build uses audited scene plates as visual proof surfaces and rebuilds the presentation shell as editable artifact-tool objects: eyebrow, title, subtitle, bridge line, claim boundary, source line, and speaker notes.

## Outputs

- PPTX: `output/slides_1_14_pptx_skeleton_v0_1/spacebiobench_slides_1_14_pptx_skeleton_v0_1.pptx`
- Contact sheet: `output/slides_1_14_pptx_skeleton_v0_1/slides_1_14_pptx_skeleton_contact_sheet.png`
- Manifest: `output/slides_1_14_pptx_skeleton_v0_1/slides_1_14_pptx_skeleton_manifest.json`
- Artifact build manifest: `output/slides_1_14_pptx_skeleton_v0_1/slides_1_14_pptx_skeleton_artifact_manifest.json`
- Speaker notes export: `output/slides_1_14_pptx_skeleton_v0_1/slides_1_14_pptx_skeleton_speaker_notes.md`

## Build Method

Script: `scripts/build_slides_1_14_pptx_skeleton.py`

The script reads `output/slides_1_14_pptx_readiness_audit_v0_1/slides_1_14_pptx_readiness_audit.json`, writes one artifact-tool slide module per slide in the Presentations workspace, then exports the deck using the bundled Codex Node runtime and artifact-tool builder.

The first attempt exposed two production lessons:

- Do not show internal viewer-language watch labels on the slide face. Those now live only in speaker notes.
- Opening slides need shorter, brand-forward titles. Slide 1 now uses `SpaceBio-Bench`; slide 3 now uses `From benchmark results to platform`.

## QA Result

PPTX export succeeded.

- Slide count: 14
- PPTX size: about 35 MB
- Contact sheet rendered and inspected
- Individual previews rendered through artifact-tool
- JSON manifests validate
- Script compiles with `PYTHONPYCACHEPREFIX=/private/tmp`
- PPTX package includes `ppt/notesSlides/notesSlide*.xml`, confirming speaker notes are embedded

Visible slide XML was checked for internal or drift terms:

- `Bridge terms`: absent
- `/Users`: absent from slide text
- `production brief`, `placeholder`, `wireframe`, `micro-plan`, `internal decision`, `draft`: absent
- `organoid`, `OSD-120`: absent from slides 1-14 visible text

Artifact-tool layout QA was run in warn-only mode. It reports expected text-image overlap because each slide intentionally uses a full-slide scene plate image behind editable overlays. After spacing revisions, filtered layout signals for `text-text`, `tight-text`, and `box-bottom-pad` no longer appear.

## Visual Review

The contact sheet now reads as one coherent dark scientific deck family. Slides 1-3 establish the opening thesis and positioning, slides 4-6 carry the methods bridge, slides 7-9 open the result family, and slides 10-14 hold the hardened core result spine.

The current skeleton is intentionally restrained. It avoids decorative card grids, avoids visible internal decision text, and keeps table-like evidence inside proof objects or source plates rather than turning the slide face into a boxed report.

## Remaining Work

This is ready for the next polishing pass, not final delivery. The next pass should:

1. Inspect each PPTX slide in PowerPoint or a full PPTX renderer.
2. Tighten the footer source and caveat hierarchy where it feels too small.
3. Decide whether to replace any PNG proof objects with native editable charts.
4. Add final chapter divider or transition rhythm if the full deck extends beyond slide 14.
5. Rerender and repeat color, grayscale, and overlap QA after any template or branding layer is added.

## Next Anchor

Run a PPTX polish pass on slides 1-14: slide-by-slide readable-size review, footer/source refinement, and deciding which proof objects should remain PNG proof plates versus be rebuilt as native editable charts.

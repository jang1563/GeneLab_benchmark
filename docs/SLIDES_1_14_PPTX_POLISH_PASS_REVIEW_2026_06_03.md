# Slides 1-14 PPTX Polish Pass Review

Date: 2026-06-03

## Scope

This pass continues from the slides 1-14 PPTX skeleton and tightens the actual
presentation draft. It does not change the scientific story or introduce v8/v9
extension material into slides 1-14. The goal is readability, footer/source
discipline, and proof-object decisions before a final template or full-talk
extension layer is added.

## Outputs

- PPTX: `output/slides_1_14_pptx_polish_v0_2/spacebiobench_slides_1_14_polished_v0_2.pptx`
- PDF smoke export: `output/slides_1_14_pptx_polish_v0_2/spacebiobench_slides_1_14_polished_v0_2.pdf`
- Contact sheet: `output/slides_1_14_pptx_polish_v0_2/slides_1_14_polished_contact_sheet.png`
- Manifest: `output/slides_1_14_pptx_polish_v0_2/slides_1_14_polished_manifest.json`
- Artifact build manifest: `output/slides_1_14_pptx_polish_v0_2/slides_1_14_polished_artifact_manifest.json`
- Speaker notes: `output/slides_1_14_pptx_polish_v0_2/slides_1_14_polished_speaker_notes.md`

## Polish Changes

- Expanded the bottom matte from 82 px to 96 px so footer content has breathing
  room.
- Increased caveat text from 10.5 pt to 11.5 pt and source text from 8.5 pt to
  9.2 pt.
- Corrected slide 7-9 source lines so they describe evidence provenance rather
  than internal scene labels.
- Replaced generic core-slide source lines with slide-specific v4-v6 evidence
  references for slides 10-14.
- Shortened core-slide caveats while preserving the claim boundaries:
  hardening is coverage evidence, temporal/RRRM is guardrail evidence, negative
  results are task-limit evidence, biological interpretation is hypothesis and
  target triage, and human translation is partial alignment only.

## QA Result

- Artifact-tool PPTX export succeeded.
- Slide count: 14.
- Rendered previews and contact sheet were generated for all slides.
- Layout JSON was generated for all slides.
- PPTX package contains 14 `ppt/notesSlides/notesSlide*.xml` files.
- Visible slide XML contains no `/Users`, `production brief`, `placeholder`,
  `wireframe`, `micro-plan`, `internal decision`, `draft`, `organoid`, or
  `OSD-120` hits.
- LibreOffice headless PDF conversion succeeded.

## Proof-Object Decision

Keep the current audited scene plates as PNG proof surfaces for this v0.2 draft.
They are the safest presentation surface because the slides combine proof
objects, claim-boundary labels, and visual hierarchy that were already audited
as a family. The editable PPTX layer remains the title, subtitle, bridge,
caveat, source line, and speaker notes.

Native editable chart rebuild is deferred. The best future candidates are slide
7 tissue hierarchy, slide 10 hardening surface, and slide 14 target-tier bars,
but rebuilding them now would add drift risk without changing the talk spine.

## Visual Decision

The v0.2 contact sheet keeps the v0.1 macro-rhythm:

- Slides 1-3: opening thesis and positioning.
- Slides 4-6: methods bridge.
- Slides 7-9: early result family.
- Slides 10-14: hardened core result spine.

The footer hierarchy is now readable at full-slide view while still subordinate
to the proof objects. No new chapter divider is needed inside slides 1-14.

## Next Anchor

Continue with one of two branches:

1. Extend the deck after slide 14 with the planned v7/v8/v9 boundary and
   platform/extension sections.
2. Apply a final institutional template or branding layer to the v0.2 polished
   slides, then rerender and repeat contact-sheet, XML, notes, and PDF smoke QA.

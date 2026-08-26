# SpaceBio-Bench Conference Deck v0.5 Review

Date: 2026-06-03

## Scope

This pass derives a 24-slide conference cut-down from the 30-slide v0.4 depth
master. The goal is to preserve the explanatory depth added in v0.4 while
compressing the deck back to a shorter talk format.

This is still not the final institutional-template deck.

## Outputs

- PPTX: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_v0_5.pptx`
- PDF smoke export: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_v0_5.pdf`
- Contact sheet: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_contact_sheet.png`
- Manifest: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_manifest.json`
- Artifact build manifest: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_artifact_manifest.json`
- Speaker notes: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_speaker_notes.md`

## Build Script

- Script: `scripts/build_spacebiobench_conference_deck.py`
- Base lineage: imports the v0.4 depth-master builder and v0.3 full-talk
  builder.
- Output role: 24-slide conference draft with native compressed teaching
  scenes, retained result scene plates, and visible data-card strips.

## Compression Decisions

- v0.4 slides 4-5 were compressed into slide 4:
  `Public data becomes a task contract`.
- v0.4 slide 6 mission-held-out protocol remains as slide 5.
- v0.4 slides 7-9 were compressed into slide 6:
  `Feature views feed one held-out score`.
- v0.4 result spine remains as slides 7-14 with compact visible data-card
  strips.
- v7/v8 boundary material remains as slides 15-17.
- v0.4 v9 glossary and payload-readiness ladder were compressed into slide 18:
  `Metadata is not payload readiness`.
- v9 platform/public bulk remain as slides 19-20.
- Organoid, OSD-120, and single-cell extension tracks remain as slides 21-23.
- Claim matrix and roadmap were compressed into slide 24:
  `Separate claims, then freeze`.

## Visual QA

- Contact sheet shows a 24-slide flow with the depth scaffolding preserved.
- Full-slide checks were performed for slides 4, 6, 18, and 24.
- Slide 24 title was shortened from `Keep claims separated, then freeze the
  package` to `Separate claims, then freeze` to remove title/subtitle overlap.

## Package QA

- Artifact-tool PPTX export succeeded.
- Slide count: 24.
- Notes count: 24 `ppt/notesSlides/notesSlide*.xml` files.
- Visible slide XML contains no `/Users`, `production brief`, `placeholder`,
  `wireframe`, `micro-plan`, or `internal decision` hits.
- Slides 1-20 contain no `organoid` or `OSD-120` hits, preserving extension
  introduction until slide 21.
- Visible XML contains the compressed teaching anchors:
  `Public data becomes a task contract`,
  `Held-out means a hidden mission`,
  `Feature views feed one held-out score`,
  `Metadata is not payload readiness`, and
  `Separate claims, then freeze`.
- LibreOffice headless PDF conversion succeeded.
- Output sizes: PPTX 32 MB, PDF 6.6 MB, contact sheet 3.5 MB.

## Residual Risks

- This version is faster and more conference-appropriate than v0.4, but the
  compressed teaching slides carry less explanatory redundancy. For a mixed
  seminar audience, v0.4 remains the safer master deck.
- The result proof surfaces remain raster scene plates, with editable title,
  subtitle, bridge, caveat, source, speaker notes, and data-card strips.
- Institutional branding/template application is still deferred.
- Timing has not been rehearsed; slide count is 24, but several result slides
  still need a disciplined talk track.

## Next Anchor

Continue from:

`SpaceBio-Bench v0.5 conference cut-down -> rehearsal timing / institutional template pass`

Recommended next branch: rehearse the 24-slide flow into a timed 15-minute or
20-minute talk script, then apply institutional template/branding only after
the talk-track timing is stable.

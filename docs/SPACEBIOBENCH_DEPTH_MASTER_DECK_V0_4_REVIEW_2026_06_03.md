# SpaceBio-Bench Depth Master Deck v0.4 Review

Date: 2026-06-03

## Scope

This pass follows the depth audit and creates a 30-slide seminar-friendly
SpaceBio-Bench deck. The goal is to make the deck teach the data unit,
evaluation protocol, feature views, score-record grammar, result interpretation,
and platform readiness boundaries directly on the slides.

This is a depth master, not the final 24-slide conference cut-down and not the
final institutional-template deck.

## Outputs

- PPTX: `output/spacebiobench_depth_master_deck_v0_4/spacebiobench_depth_master_deck_v0_4.pptx`
- PDF smoke export: `output/spacebiobench_depth_master_deck_v0_4/spacebiobench_depth_master_deck_v0_4.pdf`
- Contact sheet: `output/spacebiobench_depth_master_deck_v0_4/spacebiobench_depth_master_deck_contact_sheet.png`
- Manifest: `output/spacebiobench_depth_master_deck_v0_4/spacebiobench_depth_master_deck_manifest.json`
- Artifact build manifest: `output/spacebiobench_depth_master_deck_v0_4/spacebiobench_depth_master_deck_artifact_manifest.json`
- Speaker notes: `output/spacebiobench_depth_master_deck_v0_4/spacebiobench_depth_master_deck_speaker_notes.md`

## Build Script

- Script: `scripts/build_spacebiobench_depth_master_deck.py`
- Base lineage: imports the v0.3 builder and preserves the existing visual
  system, source parsing, and extension slide specs.
- Output role: native editable teaching scenes plus audited scene-plate proof
  surfaces with editable text overlays.

## Structure

- Slides 1-3: opening thesis and positioning.
- Slides 4-9: new teaching-first methods sequence.
- Slides 10-17: v1-v7 result spine with visible data-card strips.
- Slides 18-21: v7/v8 boundary and translation transition.
- Slides 22-25: v9 platform glossary, platform surface, public bulk alpha, and
  payload-readiness ladder.
- Slides 26-28: organoid, OSD-120, and single-cell extension tracks.
- Slides 29-30: closing claim matrix and roadmap.

## Depth Additions

- Slide 4 defines the data chain:
  `Repository -> Study -> Mission -> Tissue -> Sample -> Matrix -> Task`.
- Slide 5 shows one benchmark task as a contract with source, biology, samples,
  split, features, metric, and provenance.
- Slide 6 visualizes mission-held-out evaluation with train-side missions,
  freeze boundary, hidden mission, and score-once rule.
- Slide 7 shows gene matrix, pathway scores, and model input as parallel
  feature views feeding the same evaluator.
- Slide 8 defines the score-record grammar: task ID, split, model, feature
  view, metric, and caveat.
- Slide 9 teaches AUROC plot grammar before the first result.
- Slides 10-17 add compact data-card strips to help the audience read each
  result contract.
- Slide 19 explains why translation appears only after completed benchmark
  evidence.
- Slide 22 defines platform terms before v9 platform claims.
- Slide 25 explains that metadata readiness is not payload readiness.

## Visual QA

- Contact sheet shows a clearer teaching sequence before the first result.
- Full-slide checks were performed for slides 4, 5, 6, 7, 8, 9, 10, 19, 22,
  25, 28, and 30.
- Slide 6 bottom method claim was shortened and given more height.
- Slide 25 title was shortened from `Metadata readiness is not payload
  readiness` to `Metadata is not payload readiness` to remove title/subtitle
  overlap.

## Package QA

- Artifact-tool PPTX export succeeded.
- Slide count: 30.
- Notes count: 30 `ppt/notesSlides/notesSlide*.xml` files.
- Visible slide XML contains no `/Users`, `production brief`, `placeholder`,
  `wireframe`, `micro-plan`, or `internal decision` hits.
- Slides 1-25 contain no `organoid` or `OSD-120` hits, preserving the extension
  boundary until slide 26.
- Visible XML contains the new teaching anchors:
  `Start with the data unit`, `One task is a contract`,
  `Held-out means a hidden mission`,
  `A score is only readable with its contract`,
  `Platform terms before platform claims`, and
  `Metadata is not payload readiness`.
- LibreOffice headless PDF conversion succeeded.
- Output sizes: PPTX 32 MB, PDF 6.6 MB, contact sheet 3.9 MB.

## Residual Risks

- Result proof objects remain raster scene plates; text shell, data-card strips,
  and native teaching slides are editable.
- Some data-card strips are deliberately compact; they work as seminar cues, but
  a printed handout or manuscript supplement should carry full task tables.
- The deck is now better for a mixed seminar audience but may be too long for a
  short conference slot.
- Institutional branding/template application is still deferred.

## Next Anchor

Continue from:

`SpaceBio-Bench v0.4 depth master -> 24-slide conference cut-down`

Recommended next branch: derive a 24-slide conference version from this depth
master by preserving slides 4-9 concepts in fewer slides, keeping the AUROC
primer/data-card scaffolding, and compressing the v9 glossary/readiness ladder
into one platform-readiness slide.

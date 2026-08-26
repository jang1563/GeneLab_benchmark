# SpaceBio-Bench Conference Deck v0.5 - Extension Slide 04 Visual Review

Date: 2026-06-09

## Scope

- Deck: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_v0_5.pptx`
- Slide reviewed: slide 23, `Same-study, not mission-held-out`
- Review anchor: make OSD-120 read as a same-study diagnostic, not as mission-held-out transfer evidence.

## Change Made

Kept the original OSD-120 source-to-task proof panel and added three large read-order badges:

- `1 SOURCE STUDY` - `public OSDR record`
- `2 SAME-STUDY SPLIT` - `train/test inside study`
- `3 CLAIM LIMIT` - `not mission-held-out`

Follow-up crop fix:

- The proof panel was shifted down by 26 px on slide 23 only.
- The badge row was lowered and widened into a continuous row so the original tiny lower proof text does not show through between badges.
- Speaker notes were updated so the presenter reads the badges first.

## Review

The original proof panel has enough raster resolution, but the deck-level top matte clipped the upper OSDR source-card area when the asset was used as a full-cover scene. The source asset itself was not cropped; the issue came from the conference deck overlay.

The slide-specific offset fixes the crop while preserving the proof surface. The new badges make the claim path explicit:

- OSD-120 begins as one public OSDR source study,
- the train/test split stays inside the same study,
- the slide does not support mission-held-out or cross-species generalization.

## Zoom / Audience QA

Additional zoom crops were inspected from the rendered slide:

- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide23_zoom/01_header_claim.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide23_zoom/02_osdr_record.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide23_zoom/03_task_split.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide23_zoom/04_badges_and_strip.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide23_zoom/05_caveat_source.png`

Audience-readability judgment:

- Title, subtitle, top-right bridge, badges, bottom data strip, and footer caveat are readable.
- The OSD-120 card top is no longer clipped after the slide-specific asset offset.
- The OSDR description paragraph remains evidence texture; it is not needed for audience comprehension.
- The right-side task proof keeps the primary train/test counts and secondary check visible under zoom.
- Badge placement does not collide with the bottom data strip.

Resolution judgment:

- Rendered QA preview: 1280 x 720.
- Embedded proof panel asset: 3840 x 2160.
- Projection risk is moderate for OSDR microtext, low for the intended message because the badges and split diagram carry the required read path.

## QA

- PPTX build completed.
- Rendered slide 23 inspected visually.
- Zoom crops inspected.
- OSDR screenshot top clipping checked and corrected.
- PPTX XML: 25 slides.
- Speaker notes XML: 25 notes.
- Slide 23 badge text present: yes.
- Slide 23 note starts from badge-reading talk track: yes.
- Forbidden visible phrase `single-cell leaderboard`: none.
- Extension terms remain confined to later extension slides:
  - slide 22: `organoid`
  - slide 23: `OSD-120`
- PDF export completed.
- PDF page count: 25.
- Output sizes:
  - PPTX: 28,908,222 bytes
  - PDF: 6,273,973 bytes
  - contact sheet: 3,669,109 bytes
  - speaker notes: 15,049 bytes

## Decision

Keep the slide 23 edit with the 26 px asset offset. The OSDR screenshot is no longer clipped, and the badges make the same-study diagnostic boundary immediately understandable.

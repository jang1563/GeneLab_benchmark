# SpaceBio-Bench Conference Deck v0.5 - Extension Slide 03 Visual Review

Date: 2026-06-09

## Scope

- Deck: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_v0_5.pptx`
- Slide reviewed: slide 22, `Source records become local matrices`
- Review anchor: make the organoid source-to-matrix bridge legible for a mixed audience without promoting it as core benchmark evidence.

## Change Made

Kept the original organoid proof panel and added three large read-order badges:

- `1 SOURCE RECORD` - `public OSDR page`
- `2 LOCAL MATRIX` - `sample-by-gene table`
- `3 DRAFT BOUNDARY` - `outside core benchmark`

Speaker notes were updated so the presenter reads the badges first.

Follow-up crop fix:

- The proof panel was shifted down by 26 px on slide 22 only.
- The badge row was lowered and widened into a continuous row so the original tiny lower proof text does not show through between badges.

## Review

The original proof panel has enough raster resolution, but some source-record and explanatory labels are too small for a conference room. That is acceptable only if those details function as evidence texture, not as required reading.

Second-pass review found that the upper OSDR screenshot was too tight against the conference deck's top matte. The source asset itself was not cropped, but the deck-level matte covered the upper edge of the first source-record card. The slide-specific `assetYOffset` correction moves the proof panel down enough that both OSDR cards show their titles and top content.

The new badges make the intended transformation explicit:

- public source records are the starting objects,
- local matrices are the analysis objects,
- the organoid lane remains a draft extension outside the public bulk core.

This preserves the scientific proof surface while making the slide easier to read from the back of the room.

## Zoom / Audience QA

Additional zoom crops were inspected from the rendered slide:

- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide22_zoom/01_header_claim.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide22_zoom/02_source_records.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide22_zoom/03_local_matrices.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide22_zoom/04_badges_and_strip.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide22_zoom/05_caveat_source.png`

Audience-readability judgment:

- Title, subtitle, top-right bridge, badges, bottom data strip, and footer caveat are readable.
- The OSDR source-record crop is readable as a public record proof object, but accession-level details remain evidence texture.
- OSD-863 and OSD-871 card tops are no longer clipped after the slide-specific asset offset.
- The matrix proof crop keeps the two matrix objects and sample/gene counts visible enough under zoom.
- Badge placement does not collide with the bottom data strip.
- The original tiny `Viewer path` / `What this proves` boxes are visually superseded by the larger continuous badge row; this is a net improvement for live presentation.

Resolution judgment:

- Rendered QA preview: 1280 x 720.
- Embedded proof panel asset: 3840 x 2160.
- Projection risk is moderate for source-record microtext, low for the intended message because the badges carry the required read path.

## QA

- PPTX build completed.
- Rendered slide 22 inspected visually.
- Zoom crops inspected.
- OSDR screenshot top clipping checked and corrected.
- PPTX XML: 25 slides.
- Speaker notes XML: 25 notes.
- Slide 22 badge text present: yes.
- Slide 22 note starts from badge-reading talk track: yes.
- Forbidden visible phrase `single-cell leaderboard`: none.
- Extension terms remain confined to later extension slides:
  - slide 22: `organoid`
  - slide 23: `OSD-120`
- PDF export completed.
- PDF page count: 25.
- Output sizes:
  - PPTX: 28,907,522 bytes
  - PDF: 6,271,568 bytes
  - contact sheet: 3,663,351 bytes
  - speaker notes: 14,952 bytes

## Decision

Keep the slide 22 edit with the 26 px asset offset. The OSDR screenshots are no longer clipped, and the badges make the organoid source-to-matrix bridge immediately understandable while preserving the draft-extension boundary.

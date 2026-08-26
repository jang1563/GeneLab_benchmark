# SpaceBio-Bench Conference Deck v0.5 - Extension Slide 05 Visual Review

Date: 2026-06-09

## Scope

- Deck: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_v0_5.pptx`
- Slide reviewed: slide 24, `Metric spec exists; payload blocker remains`
- Review anchor: make the single-cell blocker legible to a mixed audience without promoting any RRRM score.

## Change Made

Kept the native status-card design and added a large bottom read rule:

- `RAW EXISTS` - `source files`
- `PAYLOAD MISSING` - `no h5ad / STARsolo`
- `SCORE BLOCKED` - `no public metric`

The speaker notes were updated so the presenter starts from this read rule.

## Review

The original slide was clean and readable, but it assumed the audience already knew why raw FASTQ files do not equal an analysis-ready single-cell payload. For a mixed ML/biology audience, the needed interpretation is:

- raw source files exist,
- processed h5ad or STARsolo payload is missing,
- public scoring is therefore blocked.

The new bottom flow makes that chain visible without adding a screenshot or overloading the card grid.

## Zoom / Audience QA

Additional zoom crops were inspected from the rendered slide:

- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide24_zoom/01_header_claim.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide24_zoom/02_status_cards.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide24_zoom/03_file_list.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide24_zoom/04_read_rule.png`
- `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-conference-deck/qa/slide24_zoom/05_footer.png`

Audience-readability judgment:

- Title, subtitle, top-right bridge, status cards, file-list counts, read rule, and footer caveat are readable.
- `0 processed payloads` is large enough to act as the numerical blocker.
- `RAW EXISTS -> PAYLOAD MISSING -> SCORE BLOCKED` is legible and does not collide with the footer.
- Technical terms `h5ad` and `STARsolo` remain present, but the flow explains their role as missing processed payloads.

Resolution judgment:

- Rendered QA preview: 1280 x 720.
- Slide 24 is native vector/text, with no embedded raster proof image.
- Projection risk is low because the key blocker is carried by large editable text.

## QA

- PPTX build completed.
- Rendered slide 24 inspected visually.
- Zoom crops inspected.
- PPTX XML: 25 slides.
- Speaker notes XML: 25 notes.
- Slide 24 read-rule text present: yes.
- Slide 24 note starts from bottom read-rule talk track: yes.
- Forbidden visible phrase `single-cell leaderboard`: none.
- Extension terms remain confined to later extension slides:
  - slide 22: `organoid`
  - slide 23: `OSD-120`
- PDF export completed.
- PDF page count: 25.
- Output sizes:
  - PPTX: 28,909,208 bytes
  - PDF: 6,277,615 bytes
  - contact sheet: 3,689,336 bytes
  - speaker notes: 15,096 bytes

## Decision

Keep the slide 24 edit. The bottom read rule makes the single-cell blocker immediately understandable while preserving the no-score claim boundary.

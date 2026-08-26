# SpaceBio-Bench Conference Deck v0.5 Speaker Notes Pass

Date: 2026-06-09

## Scope

This pass enriches the v0.5 conference deck speaker notes before any
institutional-template work. The deck remains the 24-slide conference cut-down;
visible slide content was not intentionally changed in this pass.

Primary timing source:

`docs/SPACEBIOBENCH_CONFERENCE_DECK_V0_5_REHEARSAL_TIMING_2026_06_03.md`

Primary output deck:

`output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_v0_5.pptx`

## Implementation

Script updated:

`scripts/build_spacebiobench_conference_deck.py`

The script now defines a per-slide `CONFERENCE_TALK_TRACK` map and applies it
through `enrich_speaker_notes(deck_spec)`. Each slide note now contains:

- purpose
- 20-minute timing
- 15-minute cut timing
- 20-minute talk track
- 15-minute cut instruction
- claim guardrail
- presenter cue

The same enriched notes are written into PPTX speaker notes via
`slide.speakerNotes.setText(spec.notes || "")` and into:

`output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_speaker_notes.md`

## Rebuild

Commands run:

- `python3 scripts/build_spacebiobench_conference_deck.py`
- `soffice --headless --convert-to pdf --outdir output/spacebiobench_conference_deck_v0_5 output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_v0_5.pptx`

## Notes QA

PPTX XML extraction after rebuild:

- Visible slide count: 24
- Notes slide count: 24
- Notes with `20-minute talk track`: 24
- Notes with `15-minute cut`: 24
- Notes forbidden phrase hits for `single-cell leaderboard`: none

Visible-slide XML gate remained clean:

- Visible forbidden phrase hits for `single-cell leaderboard`: none
- Extension term hits:
  - slide 21: `organoid`
  - slide 22: `OSD-120`
- Slides 1-20 extension hits: none

## PDF And Package QA

Updated PDF smoke export:

- Pages: 24
- Page size: 960.009 x 540 pt
- PDF version: 1.7
- Tagged: yes
- Encrypted: no
- JavaScript: no
- File size: 6,905,913 bytes

Output sizes after this pass:

- PPTX: 33,536,325 bytes
- PDF: 6,905,913 bytes
- Contact sheet: 3,705,338 bytes
- Speaker notes markdown: 13,216 bytes

## Visual QA

The regenerated contact sheet was reviewed after the notes rebuild. The visible
24-slide flow remained stable:

- slides 1-3: opening and positioning
- slides 4-6: compressed methods teaching
- slides 7-14: result spine
- slides 15-17: v7/v8 boundary
- slides 18-20: v9 platform/readiness
- slides 21-23: extension tracks
- slide 24: separated claim-status close

## Residual Risks

- The enriched notes are still concise rehearsal notes, not a fully scripted
  word-for-word manuscript.
- The deck still needs a final 15-minute versus 20-minute slot decision.
- Institutional branding/template application is still deferred.
- Slide 20 remains the most likely visual simplification candidate if the talk
  is shortened.

## Next Anchor

Continue from:

`SpaceBio-Bench v0.5 conference deck -> institutional template pass or final timing decision`

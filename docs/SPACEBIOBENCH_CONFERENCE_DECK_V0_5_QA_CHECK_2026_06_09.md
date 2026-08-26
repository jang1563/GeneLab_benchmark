# SpaceBio-Bench Conference Deck v0.5 QA Check

Date: 2026-06-09

## Scope

This check resumes the SpaceBio-Bench v0.5 conference deck work after the
session break. It focuses on visual QA, visible-slide XML text gates, package
integrity, and the specific close-slide wording fix needed before rehearsing or
template application.

The active deck remains the 24-slide conference cut-down from the v0.4 depth
master, not the final institutional-template deck.

## Checked Outputs

- PPTX: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_v0_5.pptx`
- PDF: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_v0_5.pdf`
- Contact sheet: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_contact_sheet.png`
- Manifest: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_manifest.json`
- Artifact manifest: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_artifact_manifest.json`
- Speaker notes: `output/spacebiobench_conference_deck_v0_5/spacebiobench_conference_deck_speaker_notes.md`

## Fix Applied

Script: `scripts/build_spacebiobench_conference_deck.py`

Slide 24 previously contained a negated caveat using the exact phrase
`single-cell leaderboard`. That phrase was removed from visible slide text to
avoid creating an accidental claim surface.

Updated slide 24 wording:

- Block card body: `No RRRM score promotion before processed payload audit passes.`
- Footer caveat: `No frozen v9 release, no v8 countermeasure, no RRRM score claim.`

The deck was rebuilt with `python3 scripts/build_spacebiobench_conference_deck.py`.

## Visual QA

- Contact sheet review: 24-slide flow is coherent.
- Slides 1-20 preserve the main benchmark-to-platform narrative.
- Slides 21-23 remain the explicit extension section.
- Slide 24 now closes with separated claim status and next-step blocks.
- Full-slide spot checks were performed for slides 4, 6, 18, 20, 23, and 24.

No major text overlap or incoherent layout break was observed in the checked
slides. Slide 20 remains visually acceptable for conference delivery, though
its proof object is small relative to the title/data-strip claim surface.

## XML Text QA

Visible PPTX XML extraction reported:

- Slide count: 24
- Notes count: 24
- Forbidden phrase hits: none for `single-cell leaderboard`
- Extension term hits:
  - slide 21: `organoid`
  - slide 22: `OSD-120`
- Slides 1-20 extension hits: none

This preserves the intended narrative boundary: extension topics are introduced
only after the core benchmark/platform story.

## PDF And Package QA

`pdfinfo` on the rebuilt PDF reported:

- Pages: 24
- Page size: 960.009 x 540 pt
- PDF version: 1.7
- Tagged: yes
- Encrypted: no
- JavaScript: no
- File size: 6,905,935 bytes

Output sizes:

- PPTX: 33,533,804 bytes
- PDF: 6,905,935 bytes
- Contact sheet: 3,705,338 bytes

## Residual Risks

- Institutional branding/template application is still deferred.
- Speaker notes have now been enriched in
  `docs/SPACEBIOBENCH_CONFERENCE_DECK_V0_5_SPEAKER_NOTES_PASS_2026_06_09.md`,
  but the deck still needs a final live timing choice before delivery.
- The deck needs a final 15-minute versus 20-minute timing decision.
- Slide 20 can stay for this cut, but it is the most likely candidate for
  enlargement or simplification if the talk is shortened.
- Existing V9/single-cell code work remains separate from this deck QA pass and
  should not be treated as part of this slide-deck checkpoint.

## Next Anchor

Continue from:

`SpaceBio-Bench v0.5 conference deck -> speaker notes enrichment or institutional template pass`

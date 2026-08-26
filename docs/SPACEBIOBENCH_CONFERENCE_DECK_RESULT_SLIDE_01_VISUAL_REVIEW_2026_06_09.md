# SpaceBio-Bench Conference Deck Result Slide 01 Visual Review

Date: 2026-06-09

## Scope

This pass applies the review criterion that slides should be visually
self-explanatory for a mixed audience, rather than depending mainly on presenter
speech.

Reviewed and updated slide:

- Slide 8: `Some tissues transfer`

No new slide was added. The goal was to make the first result slide work as a
visual worked example after the method bridge slides.

## Local Flow

- Slide 6: `What the model actually sees`
  - The model receives numerical sample views.
- Slide 7: `Feature views feed one held-out score`
  - Numerical views become one held-out score and AUROC grammar.
- Slide 8: `Some tissues transfer`
  - The first real plot now teaches how to read a result in place.

## Changes Made

Added a plot-reading guide panel to slide 8:

- `PLOT GUIDE`
- `Read this first`
- `0.5`: chance baseline
- `dot`: tissue score
- `line`: uncertainty
- `Right of 0.5 suggests transfer; claim stays tissue-specific.`

The panel points toward the result plot without covering the plot or data-card
strip.

Updated slide 8 speaker notes:

- The 20-minute track now tells the presenter to use slide 8 as the worked
  example.
- The 15-minute cut now names the minimum visual reading sequence:
  chance line, dot, uncertainty, tissue-specific claim.

## Review

Audience clarity:

- Pass. A viewer can now infer the plot-reading order from the slide itself:
  baseline, dot, uncertainty, claim boundary.
- The slide no longer relies entirely on the presenter to explain what the
  plot encodes.

Narrative flow:

- Pass. The sequence now reads as:
  `model receives numerical views -> evaluator produces one score -> first result is read on-slide`.
- This supports the idea that visual explanation can carry more of the teaching
  load than spoken explanation.

Visual fit:

- Pass. The guide panel uses the existing dark panel grammar and amber result
  accent.
- It occupies unused left-side space and does not obscure the proof plot.
- Contact sheet rhythm remains stable.

Claim safety:

- Pass. The guide explicitly keeps the result tissue-specific.
- It does not introduce a universal transfer claim, leaderboard claim, or
  mechanistic biology claim.

Timing:

- No slide count change.
- Slide 8 remains budgeted as the one result slide that should be read more
  slowly and deliberately.

## QA

PPTX XML extraction:

- Slide count: 25
- Notes count: 25
- Slide 8 contains `PLOT GUIDE` and `Read this first`.
- Slide 8 contains `Right of 0.5 suggests transfer`.
- Slide 8 notes contain `worked example`.
- Visible forbidden phrase hits for `single-cell leaderboard`: none
- Extension term hits:
  - slide 22: `organoid`
  - slide 23: `OSD-120`
- Slides 1-21 extension hits: none

PDF smoke export:

- Pages: 25
- Page size: 960.009 x 540 pt
- PDF version: 1.7
- Tagged: yes
- Encrypted: no
- JavaScript: no
- File size: 6,920,232 bytes

Output sizes:

- PPTX: 33,542,666 bytes
- PDF: 6,920,232 bytes
- Contact sheet: 3,847,533 bytes
- Speaker notes markdown: 13,969 bytes

## Decision

Keep the slide 8 visual guide.

This is the right direction for the deck: make the first result visually teach
how subsequent result slides should be read, then use presenter speech to pace
and emphasize rather than to rescue unclear slides.

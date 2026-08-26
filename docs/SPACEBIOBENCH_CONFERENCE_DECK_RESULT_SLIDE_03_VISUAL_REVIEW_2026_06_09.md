# SpaceBio-Bench Conference Deck Result Slide 03 Visual Review

Date: 2026-06-09

## Scope

This pass continues the visual-first review of the early result section.

Reviewed and updated slide:

- Slide 10: `Scale alone does not transfer`

No new slide was added.

## Local Flow

- Slide 8: plot-reading guide for tissue transfer.
- Slide 9: rescue guide for pathway summaries.
- Slide 10: model guide for tested-setting model comparison.
- Slide 11: hardening grid, next candidate for review.

## Changes Made

Added a left-side model guide panel to slide 10:

- `MODEL GUIDE`
- `Compare locally`
- `scope`: shared 6-task panel
- `baseline`: PCA-LR remains strong
- `gains`: local, often negative
- `No scale-only win; not a universal ranking.`

Updated slide 10 speaker notes:

- The 20-minute track now starts with the visual model guide.
- The 15-minute cut now says:
  `Read the guide only: shared tasks, strong baseline, local gains, no scale-only win.`

## Review

Audience clarity:

- Pass. The slide now explains what comparison scope is being used before the
  viewer tries to interpret the two model-comparison plots.
- `Scale` is reframed visually as a tested-setting comparison, not as a general
  claim about all larger or foundation models.

Narrative flow:

- Pass. Slides 8-10 now form a coherent guided-result sequence:
  `PLOT GUIDE -> RESCUE GUIDE -> MODEL GUIDE`.
- This makes the result section less dependent on presenter explanation.

Visual fit:

- Pass. The model guide uses the same left-side guide grammar as slides 8-9.
- It uses empty space and does not cover either plot or the task strip.
- At contact-sheet size, slide 10 now clearly reads as a bounded model
  comparison rather than a generic ranking slide.

Claim safety:

- Pass. The guide explicitly says `not a universal ranking`.
- It does not universalize foundation-model failure.
- It does not introduce a leaderboard, mechanism, or biological interpretation
  claim.

Timing:

- No slide count change.
- Slide 10 remains budgeted at 0:45 for the 20-minute track and 0:30 for the
  15-minute cut.

## QA

PPTX XML extraction:

- Slide count: 25
- Notes count: 25
- Slide 10 contains `MODEL GUIDE` and `Compare locally`.
- Slide 10 contains `No scale-only win; not a universal ranking`.
- Slide 10 notes contain `model guide`.
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
- File size: 6,926,668 bytes

Output sizes:

- PPTX: 33,544,412 bytes
- PDF: 6,926,668 bytes
- Contact sheet: 3,828,204 bytes
- Speaker notes markdown: 14,175 bytes

## Decision

Keep the slide 10 model guide.

The first three result slides now have enough on-slide reading structure to
support mixed-audience comprehension without relying on long spoken
explanations.

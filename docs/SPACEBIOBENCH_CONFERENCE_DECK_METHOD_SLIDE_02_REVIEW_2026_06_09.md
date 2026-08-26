# SpaceBio-Bench Conference Deck Method Slide 02 Review

Date: 2026-06-09

## Scope

This pass reviews and lightly polishes the slide immediately after the new
model-view bridge.

Reviewed slide:

- Slide 7: `Feature views feed one held-out score`

No new slide was added in this pass. The goal was to decide whether slide 7 can
carry the score-reading role after slide 6, or whether another explanatory slide
is needed before the first result.

## Local Flow

- Slide 6: `What the model actually sees`
  - Clarifies that the model receives numerical views of samples.
- Slide 7: `Feature views feed one held-out score`
  - Converts those numerical views into evaluator and AUROC score grammar.
- Slide 8: `Some tissues transfer`
  - First real result; should now function as a worked example.

## Changes Made

Slide 7 terminology was aligned with slide 6:

- Changed `Model input` to `Compressed input`.
- Changed detail text to `numeric model input`.
- Changed subtitle from `transformed model inputs` to `compressed model inputs`.
- Updated speaker notes to use `compressed model inputs`.

The AUROC primer was made more readable:

- Enlarged the primer panel slightly.
- Changed `0.5 chance` to `0.5 = chance`.
- Spaced the legend rows to reduce crowding.

The bottom reading rule was shortened:

`Read every result as: task contract + feature view + held-out score + caveat.`

## Review

Audience clarity:

- Pass. Slide 7 now inherits slide 6's `Compressed input` language cleanly.
- The score grammar is understandable without requiring a prior ML background.
- The slide still asks the presenter to explain AUROC verbally; it does not try
  to teach all of AUROC on the slide.

Narrative flow:

- Pass. Slide 6 now explains what the model receives, slide 7 explains how the
  received view becomes a score, and slide 8 can be read as the first example.
- This avoids adding a separate `How to read one benchmark result` slide for now.

Visual fit:

- Pass. The slide remains visually lighter than the result slides, which is
  appropriate for a teaching bridge.
- AUROC primer is still compact, but the legend is more legible than before.
- Contact sheet rhythm remains stable with slides 4-7 acting as the methods
  teaching block.

Claim safety:

- Pass. The slide is explicitly a methods and plot-reading primer.
- It does not introduce a real result, model superiority claim, leaderboard
  claim, or biological interpretation claim.

Timing:

- The slide remains budgeted at 1:00 for the 20-minute track and 0:45 for the
  15-minute cut.
- This is probably the maximum time it should receive; extra explanation should
  happen verbally on slide 8 as the worked example.

## QA

PPTX XML extraction:

- Slide count: 25
- Notes count: 25
- Slide 7 contains `Compressed input`.
- Slide 7 speaker notes use `compressed model inputs`.
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
- File size: 6,917,020 bytes

Output sizes:

- PPTX: 33,541,751 bytes
- PDF: 6,917,020 bytes
- Contact sheet: 3,856,492 bytes
- Speaker notes markdown: 13,860 bytes

## Decision

Keep the polished slide 7.

Do not add `How to read one benchmark result` yet. The current 6 -> 7 -> 8
sequence is strong enough to test in rehearsal:

`model receives numerical views -> evaluator produces one score -> first result is read as a worked example`

# SpaceBio-Bench Conference Deck Method Slide 01 Review

Date: 2026-06-09

## Scope

This pass adds and reviews one methodology-support slide for a mixed audience
that may include people who are not familiar with ML or LLM evaluation.

Only one slide was added:

- New slide 6: `What the model actually sees`

The deck is now a 25-slide methodology-insert conference draft. The prior
24-slide v0.5 QA notes remain useful as historical checks, but current package
counts are 25 slides and 25 notes.

## Placement

Inserted after:

- slide 5: `Held-out means a hidden mission`

Before:

- slide 7: `Feature views feed one held-out score`
- slide 8: `Some tissues transfer`

New local rhythm:

- slide 4: data becomes task contract
- slide 5: hidden mission protocol
- slide 6: model receives numerical views
- slide 7: feature views feed the held-out score
- slide 8: first result as a worked example

## Slide Intent

The slide translates ML vocabulary before the score appears.

Audience problem addressed:

Some viewers may hear "AI model" or "LLM" and imagine the model directly
understanding space biology. This slide reframes the benchmark as a numerical
evaluation pipeline:

`sample context -> numerical views -> evaluator -> hidden-mission score`

## Final Slide Text

Title:

`What the model actually sees`

Subtitle:

`AI models do not see space biology directly; the benchmark gives them numerical views of samples.`

Main labels:

- `Sample context`
- `Gene matrix`
- `Pathway scores`
- `Compressed input`
- `Train, then score`

Reading rule:

`A score is not a general AI/biology claim; it is one numerical view scored on one hidden task with a caveat.`

Caveat:

`Conceptual method bridge only; no performance or biology claim.`

## Review

Audience clarity:

- Pass. The slide uses plain nouns and a left-to-right pipeline.
- `Compressed input` is clearer than the first draft's `Embedding / input`.
- The slide should help both non-ML viewers and LLM-familiar viewers who may
  otherwise over-read the model as a general reasoning system.

Narrative flow:

- Pass. It fits naturally between hidden-mission protocol and feature/score
  primer.
- It turns the methods block from three compressed slides into four readable
  beats.
- First result now begins on slide 8, which is acceptable for a 20-minute mixed
  audience talk.

Visual fit:

- Pass. The slide matches the existing dark technical visual system.
- The three-panel flow is legible at full size and recognizable on the contact
  sheet.
- The bottom reading rule is dense but readable; it is acceptable as a
  presenter anchor rather than a read-aloud paragraph.

Claim safety:

- Pass. The slide explicitly avoids a general AI/biology claim.
- It does not introduce a performance claim, biological mechanism claim,
  leaderboard claim, or release claim.
- Visible XML still contains no `single-cell leaderboard` hit.

Timing:

- Adds 0:50 to the 20-minute track and 0:35 to the 15-minute cut.
- For a 20-minute mixed audience talk, this is worth the cost.
- For a 15-minute slot, recover time by compressing slides 9-10 or slides
  12-13 verbally rather than removing this slide immediately.

## QA

PPTX XML extraction:

- Slide count: 25
- Notes count: 25
- Slide 6 note purpose: `Plain-language ML bridge`
- Notes with `20-minute talk track`: 25
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
- File size: 6,917,036 bytes

Output sizes:

- PPTX: 33,541,743 bytes
- PDF: 6,917,036 bytes
- Contact sheet: 3,859,484 bytes
- Speaker notes markdown: 13,862 bytes

## Decision

Keep this slide as the current first methodology insert candidate.

Do not add the next explanatory slide yet. First review whether slide 7
(`Feature views feed one held-out score`) can carry the score-reading role
after this new model-view bridge. If slide 7 still feels too compressed, the
next one-slide candidate should be:

`How to read one benchmark result`

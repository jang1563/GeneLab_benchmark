# SpaceBio-Bench Conference Deck Result Slide 02 Visual Review

Date: 2026-06-09

## Scope

This pass continues the visual-first review of the result section. The goal is
to make each result understandable from the slide itself, with presenter speech
used for pacing rather than rescue.

Reviewed and updated slide:

- Slide 9: `Pathways suppress unwanted labels`

No new slide was added.

## Local Flow

- Slide 8: `Some tissues transfer`
  - First worked example with a plot-reading guide.
- Slide 9: `Pathways suppress unwanted labels`
  - Second result, now with a rescue guide that explains what is being reduced
    and how to bound the claim.
- Slide 10: `Scale alone does not transfer`
  - Next result candidate for visual-first review.

## Changes Made

Added a left-side rescue guide panel to slide 9:

- `RESCUE GUIDE`
- `Why this matters`
- `labels`: mission / hardware / gravity
- `view`: pathway summaries
- `result`: weaker selected labels
- `Selected rescue, not universal superiority.`

Updated slide 9 speaker notes:

- The 20-minute track now asks the presenter to start with the visual guide.
- The 15-minute cut now says: `Read the guide only: pathway summaries weaken selected unwanted labels.`

## Review

Audience clarity:

- Pass. Viewers can now see what the unwanted labels are, what view is being
  used, and why the result is bounded.
- The slide no longer requires the presenter to verbally define `selected
  rescue` before the audience can interpret the plots.

Narrative flow:

- Pass. Slide 8 and slide 9 now share a visual grammar:
  `guide panel -> proof plot -> compact task strip -> claim boundary`.
- This creates a stronger result-reading rhythm after the method section.

Visual fit:

- Pass. The guide uses unused left-side space and does not cover either plot.
- The amber guide border matches the early-result family accent.
- At contact-sheet size, the slide now reads as a guided result rather than a
  plot-only proof surface.

Claim safety:

- Pass. The panel explicitly says `Selected rescue, not universal superiority`.
- It does not imply pathway features always outperform gene-level models.
- It does not introduce a leaderboard, mechanism, treatment, or universal model
  claim.

Timing:

- No slide count change.
- Slide 9 remains budgeted at 0:45 for the 20-minute track and 0:30 for the
  15-minute cut.

## QA

PPTX XML extraction:

- Slide count: 25
- Notes count: 25
- Slide 9 contains `RESCUE GUIDE` and `Why this matters`.
- Slide 9 contains `Selected rescue, not universal superiority`.
- Slide 9 notes contain `visual guide`.
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
- File size: 6,923,450 bytes

Output sizes:

- PPTX: 33,543,526 bytes
- PDF: 6,923,450 bytes
- Contact sheet: 3,837,769 bytes
- Speaker notes markdown: 14,078 bytes

## Decision

Keep the slide 9 rescue guide.

The result section is moving in the right direction: each early result should
visually state the reading rule and the claim boundary, so the presenter does
not need to compensate for under-explained plots.

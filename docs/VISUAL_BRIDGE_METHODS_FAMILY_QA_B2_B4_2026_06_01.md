# Visual Bridge Methods Family QA: B2-B4

Date: 2026-06-01

Purpose: review whether B2, B3, and B4 work as a coherent methods bridge before
the result slides.

## Contact Sheet

- `output/premium_bridge_scenes/bridge_methods_b2_b4_contact_sheet.png`
- `output/premium_bridge_scenes/bridge_methods_b2_b4_contact_sheet_manifest.json`

Renderer:

- `scripts/render_bridge_methods_contact_sheet.py`

## Slides Reviewed

| Slide | Claim | Role |
|---|---|---|
| B2 | A study becomes a source-traceable task record | Explain the benchmark unit |
| B3 | The test set is a mission, not a random sample | Explain the independence unit |
| B4 | Feature processing stays on the training side | Explain leakage prevention |

## Family Verdict

Pass status:

- Draft methods bridge family candidate.

The three slides share:

- warm paper/provenance canvas;
- thin measurement rules;
- sparse editable interpretation labels;
- lower-left trust/caveat layer;
- one dominant movement per slide;
- no dashboard cards;
- no visible internal project-status jargon.

## Readability Review

B2:

- works as the object-definition bridge;
- makes the task record feel source-traceable without showing a table;
- should appear before any mission-held-out language.

B3:

- works as the split-definition bridge;
- makes the hidden mission visually obvious;
- should follow B2 immediately.

B4:

- works as the process-control bridge;
- explains why feature processing does not leak hidden-mission information;
- should follow B3 before result interpretation.

Together:

- B2 answers `what is a task?`;
- B3 answers `what is held out?`;
- B4 answers `how do we avoid learning from it?`;
- this sequence should precede Fig1/Fig2/Fig3 result slides.

## Design Risks

Remaining risks:

- B2 has the most conceptual burden. Speaker notes should name study, mission,
  sample, label, tissue, assay, and task in one clean sentence.
- B3 has the strongest visual boundary. Keep red restrained so it remains a
  trust-control cue rather than an alarm cue.
- B4 is visually quieter by design. It should not be separated from B3 in the
  final deck.
- Detailed task rows should remain outside figures. If the deck needs exact
  examples, use a separate table or appendix slide.

## Next Production Step

Proceed to B1 only after preserving this B2-B4 sequence:

1. B1: public space omics needs an evaluation layer;
2. B2: studies become source-traceable task records;
3. B3: one mission is hidden for evaluation;
4. B4: processing stays on the training side.

B1 should adapt to the same light methods/provenance grammar and should avoid
becoming a generic project overview.

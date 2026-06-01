# Visual Bridge Methods Family QA: B3-B4

Date: 2026-06-01

Purpose: review whether B3 and B4 work as a coherent methods-trust bridge
before the result slides.

## Contact Sheet

- `output/premium_bridge_scenes/bridge_methods_b3_b4_contact_sheet.png`
- `output/premium_bridge_scenes/bridge_methods_b3_b4_contact_sheet_manifest.json`

Renderer:

- `scripts/render_bridge_methods_contact_sheet.py`

## Slides Reviewed

| Slide | Claim | Role |
|---|---|---|
| B3 | The test set is a mission, not a random sample | Explain the independence unit |
| B4 | Feature processing stays on the training side | Explain leakage prevention |

## Family Verdict

Pass status:

- Draft methods bridge family candidate.

The two slides share:

- light paper/provenance canvas;
- thin measurement rules;
- restrained red trust boundary;
- large assertion-evidence headline;
- sparse editable interpretation text;
- lower-left caveat/status layer;
- no dashboard cards;
- no internal jargon in visible text.

## Readability Review

B3:

- works as the conceptual bridge;
- makes the hidden mission visually obvious;
- uses the stronger visual move and should appear first.

B4:

- works as the process-control bridge;
- explains why feature construction does not leak hidden-mission information;
- is quieter and should follow B3 immediately.

Together:

- B3 answers `what is held out?`;
- B4 answers `how do we avoid learning from it?`;
- this sequence should precede Fig1/Fig2 result interpretation.

## Design Risks

Remaining risks:

- B4 is visually drier than B3. This is acceptable for a trust-control slide,
  but it should not become the opening methods visual.
- Both slides use abstract mission labels (`M1-M4`). This is intentional for
  first-time comprehension; concrete accessions can move to speaker notes or
  appendix.
- The red boundary should remain thin. Thickening it would make the slides feel
  alarmist rather than premium.

## Next Production Step

Proceed to B2 or B1 only after deciding whether the deck needs the bridge order:

1. B1 public space omics needs an evaluation layer;
2. B2 studies become tasks;
3. B3 hidden mission;
4. B4 train-only guard.

For design development, B3/B4 should remain the grammar anchor. B1/B2 should
adapt to this light methods/provenance grammar rather than becoming generic
overview diagrams.

# v1-v9 deck spine contact sheet review

Date: 2026-06-02

## Anchor

This continues the deck-spine placement pass after the methods-bridge shortlist.
The goal is to make the full v1-v9 presentation order visible before building
final slides.

## Generated Outputs

Builder:

- `scripts/build_v1_v9_deck_spine_contact_sheet.py`

Output root:

- `output/v1_v9_deck_spine_contact_sheet_v0_1/`

Generated contact sheet:

- `output/v1_v9_deck_spine_contact_sheet_v0_1/v1_v9_deck_spine_contact_sheet.png`

Generated QA and manifests:

- `output/v1_v9_deck_spine_contact_sheet_v0_1/qa/v1_v9_deck_spine_contact_sheet_grayscale.png`
- `output/v1_v9_deck_spine_contact_sheet_v0_1/v1_v9_deck_spine_contact_sheet_manifest.json`
- `output/v1_v9_deck_spine_contact_sheet_v0_1/v1_v9_deck_spine_order.csv`

## Spine Decision

Use a 24-slide full-talk spine for the next deck assembly pass.

B1-B4 are compressed into two early methods slides:

- slide 4: evaluation layer / source-to-task bridge;
- slide 5: mission-held-out / train-only guard bridge.

This preserves room for:

- v1-v7 benchmark result core,
- v8 hypothesis-only incubator,
- v9 platform/release boundary,
- organoid, OSD-120, and single-cell extension lanes,
- final claimed/not-claimed and roadmap slides.

If the presentation becomes methods-focused, expand B1-B4 into four separate
slides and allow the full deck to become 26 slides. Do not squeeze the v1-v7
result core to make room for extension material.

## Visual QA

Verdict: ready as a deck-order planning board.

What works:

- The full 24-slide order is visible in one 4x6 board.
- v1-v7 result slides remain the middle spine.
- v8 remains visibly separated as a hypothesis-only incubator section.
- v9 platform slides appear before extension slides.
- Organoid and OSD-120 proof scenes appear only at slides 20 and 21.
- Asset gaps are visible without pretending that every slide already has a
  premium-ready figure.
- Grayscale QA keeps the sequence and grouping readable.

Limits:

- This is not a final slide and should not be pasted into the deck.
- Placeholder-like tiles are intentional asset-gap markers.
- Slides 10-14 still need static export or regenerated premium visuals.
- Slides 6, 15, 22, 23, and 24 need purpose-built bridge/boundary/roadmap
  visuals.

## Drift Check

The contact sheet reduces drift because it makes the extension limit visible:

- organoid is slide 20, not an early result slide;
- OSD-120 is slide 21 and explicitly same-study;
- single-cell remains a status/blocker slide;
- claim-boundary and roadmap still close the deck.

The remaining drift risk is now not organoid/OSD-120 overproduction. It is core
slide neglect: if slides 6-15 are not exported or rebuilt well, the deck could
become prettier around v9 than around the main benchmark results.

## Next Work

High-priority next block:

1. feature-layer bridge for genes versus pathway summaries;
2. model/metric caption bridge for direct shared-row comparisons;
3. static export or premium rebuild plan for slides 10-14;
4. v9 public bulk boundary slide QA;
5. slide-level insertion QA for slides 20-21.

The next visual sprint should favor slides 6-15 and the missing methods bridge,
not more extension proof variants.

# Core Result Slides 10-14 Integrated QA

Date: 2026-06-03

## Scope

This review covers the premium core-result slide family:

- Slide 10: v4 hardening
- Slide 11: temporal and RRRM guardrails
- Slide 12: negative boundary
- Slide 13: biological interpretation
- Slide 14: human translation boundary

The active goal is not to add new science. The goal is to verify that the five
slides work as one high-profile deck spine without visual drift, jargon drift,
or claim drift.

## Generated Outputs

Builder:

- `scripts/build_core_result_slide_set_qa.py`

Output root:

- `output/premium_core_result_slides/core_result_set_qa_v0_1/`

Generated review assets:

- `output/premium_core_result_slides/core_result_set_qa_v0_1/core_result_slides_10_14_contact_sheet.png`
- `output/premium_core_result_slides/core_result_set_qa_v0_1/core_result_slides_10_14_color_qa_sheet.png`
- `output/premium_core_result_slides/core_result_set_qa_v0_1/core_result_slides_10_14_caption_pack.md`
- `output/premium_core_result_slides/core_result_set_qa_v0_1/core_result_slides_10_14_caption_pack.json`
- `output/premium_core_result_slides/core_result_set_qa_v0_1/core_result_slides_10_14_color_variant_manifest.json`
- `output/premium_core_result_slides/core_result_set_qa_v0_1/core_result_slides_10_14_manifest.json`
- `output/premium_core_result_slides/core_result_set_qa_v0_1/core_result_slides_10_14_integrated_qa.json`

## Integrated Read Order

The five-slide sequence now reads coherently:

1. Slide 10 establishes robustness and coverage.
2. Slide 11 adds guardrails around temporal, preservation, and RRRM labels.
3. Slide 12 makes negative/failure-mode evidence intentional.
4. Slide 13 converts bounded benchmark signal into biological hypotheses.
5. Slide 14 closes the result core with a human-translation boundary.

This order avoids the original drift risk: jumping from model performance into
biology or translation before the viewer has seen robustness, limitations, and
interpretive guardrails.

## QA Findings

Automatic checks passed:

- All five source previews exist.
- All five source previews are 3840 x 2160 RGB PNGs.
- 25 color-vision QA variants were generated as 960 x 540 thumbnails.
- All visible text counts are within budget:
  - Slide 10: 44/45
  - Slide 11: 49/52
  - Slide 12: 40/50
  - Slide 13: 50/58
  - Slide 14: 45/50
- JSON manifest, caption pack, and integrated QA parse cleanly.
- Builder compiles.

Manual visual QA:

- Contact sheet inspected: family read order is coherent.
- Color QA sheet inspected: original, grayscale, deuteranopia, protanopia, and
  tritanopia variants preserve panel structure and key metric surfaces.
- First color-QA render had a footer overlap on the tritanopia row; layout was
  revised and regenerated.

Verdict:

- Pass as draft core-result slide family.
- Ready for caption-pack use and downstream PPTX assembly planning.

## Caption Pack Decision

The caption pack separates three levels:

- Presenter caption: short, plain-language bridge for a first-time viewer.
- Manuscript-style caption: longer but still claim-bounded.
- Do-not-say line: explicit overclaim guardrail.

This separation is important because several slides contain technical labels
that should not be the first spoken explanation. The spoken spine should stay
plain:

- robustness
- guardrails
- limits
- hypotheses
- translation boundary

## Claim-Drift Review

No new claim drift was introduced.

Locked boundaries:

- Slide 10: hardening and coverage evidence only, not a new leaderboard.
- Slide 11: recovery is descriptive; RRRM composition remains underpowered.
- Slide 12: negative result defines current task limits, not universal absence
  of biology or universal model failure.
- Slide 13: biological interpretation and target triage only, not mechanism,
  treatment, or countermeasure recommendation.
- Slide 14: partial human alignment only, not direct gene transfer or clinical
  proof.

Source-consistency caution retained:

- Slide 13 follows current JSON values over legacy HTML narrative text.
- `Fig7_immune_signaling.html` says 62 active L-R pairs, but the current JSON
  supports 1 active edge and 1 active pair.
- `Fig9_consensus_panel.html` says every panel gene has known drug
  interactions, but the current JSON supports 10/20 target-linked panel genes.

## Remaining Production Cautions

- Slide 10 is tight at 44/45 visible words; keep final PPTX text editable and
  avoid adding labels.
- Final deck assembly should use each scene plate as a background proof object
  and rebuild titles, source labels, caveats, and captions as editable text.
- Do not turn Slide 13 target links into a drug table inside the figure; detailed
  drug/target information belongs in backup tables or manuscript tables.
- Do not add age-amplification or PBMC temporal results into Slide 11 unless the
  slide is split.

## Next Recommended Step

Build the deck-assembly bridge for slides 1-14:

- connect methods bridge slides 4-6 into slide 10;
- use the new caption pack for slides 10-14;
- prepare PPTX overlay rules: source-proof PNG plate plus editable scientific
  interpretation text.

# Premium Bridge Color QA Review

Date: 2026-06-01

Scope:

- B1 evaluation-layer slide;
- B1-B4 premium bridge-family contact sheet;
- SpaceBio-Bench semantic color tokens.

Generated QA assets:

- B1 color modes:
  `output/premium_bridge_color_qa/b1_evaluation_layer_color_modes.png`
- B1-B4 family color modes:
  `output/premium_bridge_color_qa/b1_b4_family_color_modes.png`
- Semantic palette sheet:
  `output/premium_bridge_color_qa/semantic_palette_color_qa.png`
- Metrics:
  `output/premium_bridge_color_qa/semantic_palette_color_metrics.json`
- QA:
  `output/premium_bridge_color_qa/premium_bridge_color_qa.json`

Generator:

- `scripts/build_premium_bridge_color_qa.py`

## Executive Verdict

Conditional pass.

The current B1/B1-B4 color system is premium enough for a deck assembly mockup.
It is restrained, scientific, and more high-profile than the earlier bridge
wireframes. The palette avoids a one-note theme and does not look like generic
dashboard UI.

It is not yet final-color locked.

Two issues should be carried into the final polish pass:

1. `label_amber` is borderline for small text on the paper/canvas surface.
2. `model_purple` is visually clear but slightly louder than the rest of the
   palette; if repeated across many slides, it can drift toward SaaS/consumer
   tech rather than high-profile scientific editorial.

## Original Color Review

B1:

- pass: the four semantic zones are clear in the original render:
  `Public studies`, `Missions, tissues, samples, labels`, `Benchmark tasks`,
  and `Audited scores`;
- pass: the palette is not monochrome and does not depend on a decorative
  space background;
- pass: the headline and support copy feel calm and premium;
- conditional pass: the purple `Benchmark tasks` callout is the loudest visible
  color. It works as a single focal transition, but should not become the
  dominant deck accent.

B1-B4 family:

- pass: the family reads as one light methods/provenance system;
- pass: red is reserved for held-out/trust-boundary semantics in B3/B4;
- conditional pass: B1 is denser and more colorful than B2-B4. This is
  acceptable if B1 remains the opening concept bridge.

## Grayscale Review

The grayscale transform shows the main accessibility constraint:

- B1's overall structure survives because the slide uses position, object
  type, rail direction, and text labels;
- blue, green, teal, and amber become similar gray values;
- source dots and task-record dots should not be interpreted by color alone;
- the motion path still reads, but color contrast is not the reason it reads.

Verdict:

- conditional pass. The design follows the right rule: color supports meaning,
  but does not carry meaning alone.

## Color-Vision Simulation Review

Protanopia/deuteranopia:

- pass: B1 remains readable;
- pass: the task surface and score ledger remain separated by object type and
  placement;
- conditional pass: green and amber move toward olive tones, but this does not
  break the story.

Tritanopia:

- conditional pass: blue, green, and teal collapse more strongly;
- pass: the slide still works because each semantic zone is labeled and placed
  in a different region;
- watch item: do not add unlabeled multi-color dot legends to main slides.

## Quantitative Checks

Token contrast against paper:

- `source_blue`: 5.19
- `bio_green`: 4.64
- `assay_teal`: 3.34
- `label_amber`: 2.97
- `test_red`: 6.19
- `model_purple`: 5.98
- `muted`: 5.35
- `rule`: 1.92

Interpretation:

- headline, muted text, blue, green, red, and purple are comfortably readable;
- teal is acceptable for callouts but should stay bold or larger than tiny text;
- amber is below the preferred 3.0 threshold on warm/cool canvas and should not
  be used for small critical text;
- rule is intentionally low contrast and should remain infrastructure only.

Pairwise colorblind warnings:

- blue/green and green/teal become close under tritanopia;
- muted can become close to semantic colors under some transforms;
- purple is clear in normal and red-green simulations, but less distinct from
  muted under tritanopia.

These are not blockers because the rendered slides use spatial separation and
labels. They are design constraints for later slides.

## Premiumness Review

What feels premium:

- warm/cool off-white canvas with faint measurement texture;
- semantic colors used in small doses;
- no card-box palette or dashboard tiles;
- restrained shadows and evidence surfaces;
- red reserved for trust/test semantics;
- enough color variety to avoid a beige or monotone deck.

What weakens premiumness:

- overly bright purple can feel like product-marketing UI if overused;
- amber is a little too light for small text;
- B1's source field is richer than B2-B4, so final family polish should decide
  whether this is intentional pacing or an inconsistency.

## Design Decision

Do not rerender the slide set immediately.

Reason:

- B1-B4 pass current deck-assembly quality;
- the biggest color risks are polish-level, not concept-breaking;
- premature global color-token changes would require rerendering B1-B4 and
  checking downstream figure consistency.

Carry-forward rule:

- color may classify, but structure must explain;
- if a color appears as text, keep contrast above 3.0 minimum and preferably
  above 4.5 for small labels;
- amber should be darkened before final deck export if it remains a text color;
- purple should stay a secondary accent, never the dominant slide identity.

Recommended final-polish token candidates:

- `label_amber`: darken from `#C4861B` toward `#A36F13`;
- `model_purple`: calm from `#6D3EDB` toward `#5D4BA8` if it appears on more
  than a few slides;
- keep `source_blue`, `bio_green`, and `assay_teal` but avoid placing them as
  unlabeled adjacent categories in tiny legends.

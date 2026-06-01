# B1 Color Token A/B Review

Date: 2026-06-01

Purpose:

- test whether the final-polish color warnings from the bridge color QA should
  become token changes;
- compare darkened amber and calmer purple on the actual B1 render before
  changing global visual identity tokens.

Generator:

- `scripts/build_b1_color_token_ab.py`

Comparison sheets:

- original:
  `output/premium_bridge_color_ab/b1_color_token_ab_original.png`
- grayscale:
  `output/premium_bridge_color_ab/b1_color_token_ab_grayscale.png`
- tritanopia:
  `output/premium_bridge_color_ab/b1_color_token_ab_tritanopia.png`
- palette delta:
  `output/premium_bridge_color_ab/b1_color_token_ab_palette_delta.png`

## Candidates

Baseline:

- current visual identity v0.1 tokens;
- `label_amber`: `#C4861B`;
- `model_purple`: `#6D3EDB`.

Variant A:

- darken amber only;
- `label_amber`: `#A36F13`;
- `model_purple`: unchanged.

Variant B:

- darken amber and calm purple;
- `label_amber`: `#A36F13`;
- `model_purple`: `#5D4BA8`.

## Quantitative Findings

Amber contrast:

- baseline: 2.97 against paper;
- A/B: 4.15 against paper.

This fixes the main contrast warning from the previous color QA.

Purple contrast:

- baseline/A: 5.98 against paper;
- B: 6.62 against paper.

Purple contrast improves in B, but this is not the main decision criterion. The
more important issue is whether purple remains visually distinct from muted
infrastructure colors in color-vision simulations.

Purple vs muted under tritanopia:

- baseline/A: distance 31.6;
- B: distance 20.8.

B therefore makes purple calmer but less distinct from muted under tritanopia.

## Visual Findings

Original comparison:

- A is the best balance. It improves the amber cue while preserving the current
  B1 hierarchy.
- B feels more editorial, but the benchmark-task accent loses useful energy.
- Baseline remains acceptable but leaves amber too light for future small-text
  use.

Grayscale comparison:

- all variants remain structurally readable because B1 is built on position,
  rail direction, object type, and labels;
- A does not disrupt the grayscale hierarchy;
- B does not create a grayscale improvement large enough to justify the purple
  change.

Tritanopia comparison:

- blue, green, and teal already compress, so the slide must keep relying on
  spatial grouping and labels;
- A does not add a new accessibility penalty;
- B makes the task-purple/muted distinction weaker, so it should not become the
  global token yet.

## Decision

Recommended final-polish token change:

- adopt `label_amber = #A36F13`.

Do not change `model_purple` globally yet.

Reason:

- current purple is slightly loud, but it keeps task objects distinct;
- a calmer purple should be tested in a wider deck context before becoming the
  global model/task color;
- the immediate confirmed issue is amber contrast, and Variant A fixes it with
  minimal side effects.

## Implementation Recommendation

Initial recommendation before implementation:

- do not mutate the committed visual identity config until the A/B comparison
  was rendered and reviewed.

Implementation status:

- completed in the subsequent color-token application pass;
- `label_amber = #A36F13` was applied to the visual identity token file and the
  B2-B4 premium bridge renderer;
- `model_purple` remains unchanged.

The completed implementation pass:

1. updated the visual identity token file with `label_amber = #A36F13`;
2. updated the B2-B4 premium bridge renderer amber value to match;
3. rerendered B1-B4 together;
4. repeated color QA on the full family.

The A/B script now pins the pre-polish baseline amber value so the comparison
remains auditable after the global token change.

# Premium Bridge Color Token Application Review

Date: 2026-06-01

Purpose:

- apply the B1 color-token A/B decision to the actual bridge-rendering system;
- rerender B1-B4;
- rerun family and color accessibility QA.

## Applied Change

Adopted:

- `label_amber`: `#A36F13`

Not changed:

- `model_purple`: remains `#6D3EDB`

Files updated:

- `config/visual_identity/spacebio_bench_visual_identity_v0_1.json`
- `scripts/build_bridge_premium_rebuild_scenes.py`
- `scripts/build_b1_color_token_ab.py`

The A/B script now preserves the pre-polish amber baseline by explicitly
overriding the baseline candidate to `#C4861B`. This keeps the comparison
auditable after the global visual identity token changes.

## Regenerated Assets

B1:

- `output/premium_bridge_scenes/b1_evaluation_layer/rendered_preview.png`
- `output/premium_bridge_scenes/b1_evaluation_layer/qa.json`

B2-B4 premium rebuilds:

- `output/premium_bridge_rebuild_scenes/b2_study_to_task_premium/rendered_preview.png`
- `output/premium_bridge_rebuild_scenes/b3_mission_held_out_premium/rendered_preview.png`
- `output/premium_bridge_rebuild_scenes/b4_train_only_guard_premium/rendered_preview.png`

Family review:

- `output/premium_bridge_family_review/b1_b4_premium_family_contact_sheet.png`
- `output/premium_bridge_family_review/b1_b4_premium_family_contact_sheet_qa.json`

Color QA:

- `output/premium_bridge_color_qa/b1_evaluation_layer_color_modes.png`
- `output/premium_bridge_color_qa/b1_b4_family_color_modes.png`
- `output/premium_bridge_color_qa/semantic_palette_color_qa.png`
- `output/premium_bridge_color_qa/premium_bridge_color_qa.json`

## QA Results

Automatic color QA:

- pass: no remaining text-color contrast warning against the paper surface;
- pass: `label_amber` contrast improved from 2.97 to above the small-label
  threshold;
- expected warning: blue/green/teal remain close under grayscale/tritanopia;
- new non-blocking warning: amber and test red are closer under deuteranopia.

Manual visual review:

- pass: B1 remains premium and readable after rerender;
- pass: B2's label contrast looks more mature and less decorative;
- pass: B3/B4 remain semantically unchanged because red continues to carry
  held-out/trust-boundary meaning;
- conditional pass: amber and red should not be used as adjacent unlabeled
  categories in future tiny legends.

## Decision

Keep the amber-dark token.

Reason:

- it resolves the only direct contrast warning;
- it improves the visual maturity of label/task details;
- it does not disrupt B1-B4 family hierarchy;
- the new amber/red colorblind warning is mitigated by role separation, spatial
  separation, and shape separation in the current bridge slides.

Carry-forward rule:

- amber can be used for label/contrast semantics;
- red remains reserved for held-out/test/trust-boundary semantics;
- do not place amber and red as same-scale adjacent categorical dots without
  text labels or shape differences.

## Remaining Polish Work

Before final deck export:

- normalize status/source note placement across B1-B4;
- decide whether B1's richer source field should remain as opening-bridge
  pacing or be quieted to match B2-B4;
- keep purple as a secondary task/model accent until a wider deck-wide color
  pass proves it should be calmer.

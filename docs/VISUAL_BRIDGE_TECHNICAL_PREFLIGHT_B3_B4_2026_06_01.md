# Visual Bridge Technical Preflight: B3-B4

Date: 2026-06-01

Purpose: translate the B3/B4 content briefs into concrete technical production
contracts before any polished rendering. This is the checkpoint that prevents
method slides from becoming attractive but under-specified schematics.

Companion files:

- `docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md`
- `docs/VISUAL_METHODS_BRIDGE_AND_CONSULTING_BRIEF_2026_06_01.md`
- `docs/VISUAL_TECHNICAL_PRODUCTION_PROTOCOL_2026_06_01.md`
- `scripts/build_visual_bridge_scene_contracts.py`
- `scripts/audit_visual_scene_contracts.py`

## Why B3 And B4 First

B3 and B4 are the methodological trust hinge.

If a first-time viewer misunderstands these slides, every result figure can be
misread:

- Fig1 tissue transfer can look like ordinary sample-level CV.
- Fig2 pathway behavior can look like feature selection after seeing test data.
- Fig3 model comparisons can look like a generic leaderboard.
- v9 release/resource slides can look like platform polish without benchmark
  discipline.

Therefore B3/B4 must be technically precise before they are visually polished.

## B3 Contract

Slide ID:

- `b3_mission_held_out`

Output family:

- `output/premium_bridge_scenes/b3_mission_held_out/`

Decision headline:

> The test set is a mission, not a random sample

Audience question:

> Why is mission-held-out evaluation central to this benchmark?

Primary visual move:

One mission moves behind a clean boundary while the remaining missions feed the
training lane.

Required visible elements:

- `train on prior missions`
- `hide one mission`
- `score after training`
- `held-out mission evaluation`

Required caveat:

- `No hidden-mission samples are used for training.`

Evidence sources:

- `docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md`
- `docs/VISUAL_METHODS_EXPLANATION_GAP_MAP_2026_06_01.md`
- `docs/VISUAL_METHODS_STORYBOARD_2026_06_01.md`
- `docs/V1_V9_PRESENTATION_AND_MANUSCRIPT_MASTER_OUTLINE_2026_05_31.md`
- `docs/PROJECT_SLIDE_CONTENT_INVENTORY_V1_TO_V9_2026_05_31.md`

Forbidden visible terms:

- `LOMO`
- `random CV`
- `cross-validation`
- `payload`
- `RRRM`
- `alpha`
- `NES`
- `macro-F1`
- raw local paths

Technical risk:

- Red boundary can become alarmist. Use restrained red and thin-rule structure.

## B4 Contract

Slide ID:

- `b4_train_only_guard`

Output family:

- `output/premium_bridge_scenes/b4_train_only_guard/`

Decision headline:

> Feature processing stays on the training side

Audience question:

> How does the benchmark avoid learning from the mission it is supposed to test?

Primary visual move:

Two lanes move in parallel until the hidden mission joins only at scoring.

Required visible elements:

- `choose features`
- `scale`
- `fit model`
- `score hidden mission`
- `train-only processing`

Required caveat:

- `Feature choices and scaling are learned before the hidden mission is scored.`

Evidence sources:

- `docs/VISUAL_BRIDGE_CONTENT_BRIEFS_B1_B4_2026_06_01.md`
- `docs/VISUAL_METHODS_EXPLANATION_GAP_MAP_2026_06_01.md`
- `docs/VISUAL_METHODS_STORYBOARD_2026_06_01.md`
- `docs/VISUAL_METHODS_EXPLAINER_PILOT_2026_06_01.md`
- `output/premium_methods_scenes/methods_data_to_evaluation_overview_manifest.json`

Forbidden visible terms:

- package names;
- function names;
- `LOMO`;
- `payload`;
- `artifact`;
- `macro-F1`;
- raw local paths.

Technical risk:

- Process-control visuals can become software architecture. Use verbs, lanes,
  and a trust gate rather than package or code structure.

## Shared Pre-Render Checks

The generated contract must pass:

1. `manifest.json` exists.
2. `overlay_spec.json` exists.
3. `qa.json` exists.
4. slide IDs match across all three files.
5. canvas is `3840 x 2160`, 16:9.
6. content brief path exists.
7. all evidence source paths exist.
8. visible text is under the 45-word bridge-slide budget.
9. forbidden visible terms are absent.
10. normalized coordinates stay inside the safe area unless explicitly marked.
11. output image paths are declared even if not yet rendered.
12. `qa.json` records that this is `pre_render`, not a finished slide.

## What This Does Not Approve

Passing preflight does not mean the slide is visually acceptable.

It only means the slide is ready for a careful prototype render. After
rendering, the post-render gate still requires:

- full-size visual inspection;
- thumbnail/contact-sheet inspection;
- text overlap check;
- caveat visibility check;
- source/claim-boundary review;
- optional color/grayscale QA if color carries meaning.

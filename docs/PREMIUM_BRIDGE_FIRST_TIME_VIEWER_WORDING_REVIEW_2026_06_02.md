# Premium Bridge First-Time-Viewer Wording Review

Date: 2026-06-02

Purpose:

- remove avoidable jargon from B1-B4 visible text;
- keep methods-critical terms where they carry real scientific meaning;
- preserve the premium proof-plate style while making the workflow easier to
  understand on first view.

## Audit Result

P0 wording issues:

- B1 support line used `traceable tasks`, which repeated the traceability idea
  after the status/source layer had already been simplified;
- B2 headline used `traceable contract` and `loose matrix`, which read more
  like an internal benchmark-design metaphor than a first-view explanation;
- B2 endpoint label used `task contract`, reinforcing the same metaphor.

Kept for now:

- B3 `held-out mission` and `mission-held-out evaluation`, because the
  evaluation split is the method claim;
- B4 `choose features`, `scale`, and `fit model`, because those are the actual
  leakage-sensitive processing steps.

## Applied Changes

B1:

- changed support line from `traceable tasks` to `benchmark tasks`.

B2:

- changed headline to `A benchmark task preserves source context`;
- changed endpoint callout to `benchmark task`;
- changed metadata and contact-sheet captions from task-contract language to
  benchmark-task/source-context language.

## Regenerated Assets

- `output/premium_bridge_scenes/b1_evaluation_layer/rendered_preview.png`
- `output/premium_bridge_rebuild_scenes/b2_study_to_task_premium/rendered_preview.png`
- `output/premium_bridge_rebuild_scenes/b3_mission_held_out_premium/rendered_preview.png`
- `output/premium_bridge_rebuild_scenes/b4_train_only_guard_premium/rendered_preview.png`
- `output/premium_bridge_family_review/b1_b4_premium_family_contact_sheet.png`
- `output/premium_bridge_color_qa/b1_b4_family_color_modes.png`

## QA Results

Automatic checks:

- B1-B4 rendered previews remain 3840 x 2160;
- visible word counts remain below the 45-word bridge-slide limit:
  B1 39, B2 31, B3 27, B4 28;
- color QA produced no text contrast warnings.

Manual visual review:

- pass: B1 support line is cleaner and no longer repeats traceability jargon;
- pass: B2 reads more directly as a source-context-preserving task build;
- pass: B2 contact-sheet caption now matches the slide instead of retaining
  stale contract language;
- conditional pass: B3/B4 remain technical, but those terms belong in the
  method explanation and should be bridged in speaker notes.

## Verdict

Pass for deck assembly.

The P0 visible-wording issue is resolved. The next methods-explanation pass
should not rewrite B3/B4 into generic language; it should add concise captions
or speaker-note bridges that explain why mission-held-out evaluation and
train-only processing prevent leakage.

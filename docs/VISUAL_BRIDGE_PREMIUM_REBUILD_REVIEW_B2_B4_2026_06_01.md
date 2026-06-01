# Visual Bridge Premium Rebuild Review: B2-B4

Date: 2026-06-01

Scope:

- `b2_study_to_task_premium`
- `b3_mission_held_out_premium`
- `b4_train_only_guard_premium`

## Decision

The earlier B2-B4 bridge renders are demoted to logic wireframes. They should
not be used in a premium deck.

The new premium rebuild family is the current deck-assembly candidate, with one
important caveat: these are still slide candidates, not final locked art. They
must be reviewed again inside the actual deck rhythm.

## Outputs

Output root:

- `output/premium_bridge_rebuild_scenes/`

Renderer:

- `scripts/build_bridge_premium_rebuild_scenes.py`

Contact sheet:

- `output/premium_bridge_rebuild_scenes/bridge_methods_premium_rebuild_contact_sheet.png`
- `output/premium_bridge_rebuild_scenes/bridge_methods_premium_rebuild_contact_sheet_manifest.json`

Slide outputs:

- `output/premium_bridge_rebuild_scenes/b2_study_to_task_premium/rendered_preview.png`
- `output/premium_bridge_rebuild_scenes/b3_mission_held_out_premium/rendered_preview.png`
- `output/premium_bridge_rebuild_scenes/b4_train_only_guard_premium/rendered_preview.png`

## Validation

Commands run:

```bash
PYTHONPYCACHEPREFIX=output/.pycache python3 -m py_compile scripts/build_bridge_premium_rebuild_scenes.py scripts/audit_visual_scene_contracts.py
python3 scripts/build_bridge_premium_rebuild_scenes.py --slide all
python3 scripts/audit_visual_scene_contracts.py --root output/premium_bridge_rebuild_scenes --output-root output/premium_bridge_rebuild_scenes
```

Audit result:

- 277 checks;
- 0 error failures;
- 0 warnings.

Rendered PNG dimensions:

- each slide preview: 3840 x 2160;
- contact sheet: 4180 x 1254.

## What Changed From The Rejected Draft

Major improvements:

- larger evidence surfaces;
- stronger visual mass at contact-sheet scale;
- more visible depth through layered source objects, mission evidence, process
  stations, shadows, and translucent fields;
- reduced cardbox framing after removing the heavy outer panel border and
  shadow;
- clearer B2 source-to-contract flow;
- clearer B3 mission-level train/test boundary;
- clearer B4 hidden-mission bypass lane.

Specific fixes after first premium rebuild render:

- B2 task-contract strip moved right so the final flow cue no longer reads
  backward.
- B4 hidden mission moved to a lower bypass lane.
- B4 `hidden mission bypass` label moved away from the M4 node.
- Full-canvas evidence band replaced the large bordered panel to reduce the
  cheap card-box feeling.

## Slide-Level Review

### B2

Verdict:

- conditional pass as a premium rebuild candidate.

What works:

- the slide now reads as source evidence consolidating into a task contract;
- the task object is a document strip, not a table;
- the central flow remains visible at thumbnail scale;
- visible text avoids internal project vocabulary.

Remaining risk:

- the phrase `loose matrix` is plain, but may be too informal for a manuscript
  version. A deck version can keep it; a paper figure may prefer `expression
  matrix alone`.

### B3

Verdict:

- pass as the strongest current bridge candidate.

What works:

- mission-level split is visually obvious;
- red boundary is present but not alarmist;
- hidden mission and scoring relation are readable;
- this slide can serve as the grammar anchor for the bridge family.

Remaining risk:

- if later result slides use a much darker or more photographic grammar, B3 may
  need a stronger transition cue.

### B4

Verdict:

- conditional pass as a premium rebuild candidate.

What works:

- train-only processing now has a larger control surface;
- hidden mission bypass is separated from training choices;
- feature selection, scaling, model fitting, and scoring are distinct;
- no package/function names appear.

Remaining risk:

- it remains a procedural control slide and will feel drier than B3. Keep it
  immediately after B3 so the audience reads it as the leakage-control
  consequence of the hidden-mission split.

## Deck Use Guidance

Recommended bridge order:

1. B1: public space omics needs an evaluation layer.
2. B2: a study becomes a traceable task contract.
3. B3: the test set is an entire mission.
4. B4: training choices stop before the hidden mission.

Use the premium rebuild family, not the earlier B2-B4 draft family, for the
next deck assembly mockup.

Before final deck export:

- inspect these slides inside the deck template;
- check whether the headline scale matches surrounding figures;
- decide whether B2 keeps `loose matrix` or changes to `expression matrix
  alone`;
- run one grayscale/contact-sheet pass after B1 is rebuilt.

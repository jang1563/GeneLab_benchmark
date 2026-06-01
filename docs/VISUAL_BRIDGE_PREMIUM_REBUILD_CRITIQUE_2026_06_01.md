# Visual Bridge Premium Rebuild Critique

Date: 2026-06-01

Scope:

- `b2_study_to_task`
- `b3_mission_held_out`
- `b4_train_only_guard`

Verdict on the current draft bridge slides:

- The current B2-B4 renders are useful logic wireframes.
- They are not premium presentation slides.
- They should not be used as-is in a 2026 ISGP-level deck.

## Why The Current Draft Looks Weak

### 1. Insufficient Visual Mass

The central objects are too small and too thin relative to the 16:9 canvas. At
contact-sheet scale, each slide reads as scattered schematic marks on a pale
background rather than as a designed scientific scene.

Premium correction:

- make the evidence surface occupy the central field;
- use one dominant scene object per slide;
- give each slide a clear optical center of gravity.

### 2. Layering Exists In Theory, Not In Perception

The draft follows the Z-layer language in file structure, but the rendered
image does not visibly convey depth. Shadows are too faint, surfaces are too
flat, and the foreground interpretation layer is not anchored to an underlying
scientific proof object.

Premium correction:

- Z0: atmospheric canvas with subtle paper/evidence texture;
- Z1: measurement rail, assay lanes, mission ticks;
- Z2: large evidence surface with shadow/backplate;
- Z3: editable interpretation labels tied to specific objects;
- Z4: compact trust/caveat line;
- Z5: one dominant focus path or boundary.

### 3. The Slides Look Like Generic Workflow Diagrams

B2 looks like a simple left-to-right process. B3 and B4 look like abstract
method sketches. None of them yet conveys the specific scientific object:
public space-omics records becoming benchmark contracts under mission-held-out
evaluation.

Premium correction:

- B2 must look like source evidence being consolidated into a task contract,
  not a table row.
- B3 must look like mission-level evidence being separated by an evaluation
  boundary, not a thin split diagram.
- B4 must look like a train-only control surface, not software architecture.

### 4. Too Much Empty Space, Too Little Hierarchy

The headline is readable, but the body region has weak hierarchy. Important
objects are close in size to decorative ticks, making the slide feel tentative.

Premium correction:

- enlarge primary objects;
- reduce low-value micro marks;
- add stronger grouping through rails, lanes, and field shading;
- keep visible text sparse but make it purposeful.

### 5. The Current Contact Sheet Exposes The Problem

At thumbnail scale, the slides collapse into thin lines. This is a useful
failure signal: a premium slide must remain recognizable in a contact sheet
before it is worth polishing at full size.

Premium correction:

- require thumbnail recognition before declaring the slide usable;
- reject slides where the visual thesis disappears at 20-25% scale.

## Rebuild Standard

A rebuilt bridge slide must pass these criteria:

1. The visual thesis is recognizable in three seconds at contact-sheet scale.
2. The central proof object has visible depth, not just outline geometry.
3. The slide contains one dominant movement, not a row of small icons.
4. Text is sparse and plain-language.
5. No internal project-management terms appear in visible text.
6. Tables are not embedded inside figures. Exact rows stay in tables,
   appendix, or source files.
7. The slide works as an editable scientific interpretation shell over a
   source-like PNG scene.

## Production Decision

Do not iterate the existing draft in place.

Create a separate premium rebuild family:

- output root: `output/premium_bridge_rebuild_scenes/`
- renderer: `scripts/build_bridge_premium_rebuild_scenes.py`
- review: `docs/VISUAL_BRIDGE_PREMIUM_REBUILD_REVIEW_B2_B4_2026_06_01.md`

The old B2-B4 files remain as audit trail and logic references. The rebuilt
family becomes the candidate for deck assembly only if it passes both full-size
and contact-sheet review.

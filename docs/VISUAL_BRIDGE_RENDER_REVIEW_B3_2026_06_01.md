# Visual Bridge Render Review: B3

Date: 2026-06-01

Slide:

- `b3_mission_held_out`

Claim:

> The test set is a mission, not a random sample.

## Outputs

- `output/premium_bridge_scenes/b3_mission_held_out/scene_plate.png`
- `output/premium_bridge_scenes/b3_mission_held_out/rendered_preview.png`
- `output/premium_bridge_scenes/b3_mission_held_out/rendered_preview.pdf`
- `output/premium_bridge_scenes/b3_mission_held_out/manifest.json`
- `output/premium_bridge_scenes/b3_mission_held_out/overlay_spec.json`
- `output/premium_bridge_scenes/b3_mission_held_out/qa.json`

Renderer:

- `scripts/build_bridge_method_scenes.py`

## Technical Validation

Commands run:

```bash
PYTHONPYCACHEPREFIX=output/.pycache python3 -m py_compile scripts/build_visual_bridge_scene_contracts.py scripts/build_bridge_method_scenes.py scripts/audit_visual_scene_contracts.py
python3 scripts/build_visual_bridge_scene_contracts.py
python3 scripts/audit_visual_scene_contracts.py
python3 scripts/build_bridge_method_scenes.py --slide b3_mission_held_out
python3 scripts/audit_visual_scene_contracts.py --root output/premium_bridge_scenes/b3_mission_held_out --output-root output/premium_bridge_scenes/b3_mission_held_out
```

Results:

- pre-render contract audit: 169 checks, 0 errors, 0 warnings;
- B3 post-render contract audit: 85 checks, 0 errors, 0 warnings;
- rendered PNG size: 3840 x 2160;
- rendered PDF: one-page PDF.

## Manual Visual Review

Initial render issue:

- The renderer treated overlay `y` coordinates as Matplotlib bottom-origin
  coordinates, while the contract used slide-style top-origin coordinates.
- Result: the headline appeared at the bottom of the slide.

Fix:

- Added top-origin to Matplotlib y-coordinate conversion in
  `scripts/build_bridge_method_scenes.py`.

Second render issue:

- `train on prior missions` overlapped M2/rail.
- `hide one mission` overlapped M4.
- `score after training` sat too close to the score glyph.

Fix:

- Updated B3 overlay coordinates in
  `scripts/build_visual_bridge_scene_contracts.py`.

Current visual verdict:

- Draft bridge-slide candidate.
- The headline is readable in under five seconds.
- The train/test split is visually obvious.
- The red boundary is restrained enough to read as a trust boundary rather than
  an alarm graphic.
- The caveat is visible at lower-left without dominating the slide.
- No internal terms such as `LOMO`, `payload`, `RRRM`, `alpha`, `NES`, or
  `macro-F1` appear in visible text.

## Remaining Before Final Deck Use

B3 should still be reviewed again after B4 is rendered because the two methods
trust slides must feel like one family.

Potential refinements for a later pass:

- decide whether the headline should be slightly smaller for deck-wide
  consistency;
- consider adding a very subtle chapter cue only if the deck needs section
  navigation;
- test B3 and B4 together in a contact sheet before locking either slide.

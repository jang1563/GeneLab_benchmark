# Visual Identity Style-Frame Review

Date: 2026-06-01

Purpose: compare three visual identity directions before rebuilding B1-B4.

## Outputs

Renderer:

- `scripts/build_visual_identity_style_frames.py`

Output root:

- `output/visual_identity_style_frames/`

Contact sheet:

- `output/visual_identity_style_frames/visual_identity_style_frame_contact_sheet.png`
- `output/visual_identity_style_frames/visual_identity_style_frame_contact_sheet_manifest.json`

Frames:

- `output/visual_identity_style_frames/editorial/rendered_preview.png`
- `output/visual_identity_style_frames/mission_system/rendered_preview.png`
- `output/visual_identity_style_frames/hybrid/rendered_preview.png`

## Shared Test Message

All three style frames use the same conceptual message:

> public studies -> auditable task surface -> held-out mission boundary

This keeps the comparison about visual identity rather than content.

## Technical Validation

Commands run:

```bash
PYTHONPYCACHEPREFIX=output/.pycache python3 -m py_compile scripts/build_visual_identity_style_frames.py
python3 scripts/build_visual_identity_style_frames.py --style all
python3 -m json.tool output/visual_identity_style_frames/visual_identity_style_frame_contact_sheet_manifest.json
python3 -m json.tool output/visual_identity_style_frames/hybrid/qa.json
```

Rendered outputs:

- each style frame: 3840 x 2160 PNG plus PDF;
- contact sheet: 4180 x 1276 PNG.

## Scored Directions

| Direction | Trust | Premium | Thumbnail | Manuscript | Reproducibility | Verdict |
|---|---:|---:|---:|---:|---:|---|
| Scientific editorial | 5 | 3 | 3 | 5 | 5 | Not primary |
| Mission system | 3 | 5 | 5 | 2 | 4 | Accent only |
| Mission-grade editorial hybrid | 5 | 5 | 4 | 4 | 5 | Primary direction |

## Direction Review

### Scientific Editorial

Strengths:

- most manuscript-portable;
- lowest risk of over-stylization;
- strongest compatibility with journal figure conventions.

Weaknesses:

- too quiet for opening bridge slides;
- weak depth and weak premium first impression;
- risks looking like a plain methods figure rather than a branded deck system.

Decision:

- keep as the manuscript/export fallback, not the primary deck identity.

### Mission System

Strengths:

- strongest immediate identity;
- excellent thumbnail recognition;
- clear mission/control-room atmosphere.

Weaknesses:

- too dark for the whole deck;
- risks drifting into sci-fi or operational UI;
- weaker manuscript portability.

Decision:

- use selectively for mission/control accent moments, not as the default deck
  grammar.

### Mission-Grade Editorial Hybrid

Strengths:

- best balance of scientific trust and premium deck presence;
- preserves source/evidence discipline while adding z-depth;
- avoids NASA imitation and avoids decorative space wallpaper;
- suitable for B1-B4 bridge system after refinement.

Weaknesses:

- held-out/test object should sit farther inside the right safe area in final
  B1-B4 renders;
- depth must remain semantic, not decorative;
- if overused, the warm canvas can become too beige unless anchored by
  mission navy/source blue/bio green.

Decision:

- select hybrid as the primary deck visual identity direction.

## Z-Stack Implications

The selected hybrid direction should use this layer hierarchy:

| Layer | Use In B1-B4 |
|---|---|
| Z0 | warm/cool mission-grade canvas, subtle grain only |
| Z1 | mission rails, assay lanes, source ticks |
| Z2 | full-width evidence band, not a bordered card |
| Z3 | source stack, task document, mission nodes |
| Z4 | headline and sparse callouts |
| Z5 | held-out/test boundary or train-only guard |
| Z6 | caveat/source-traceability status |
| Z7 | one source-to-task or train-to-score movement |

Depth should be visible but not theatrical. The proof object should carry the
main elevation; text and background rules should remain flat.

## Implementation Decision

Next renderer work should:

1. load `config/visual_identity/spacebio_bench_visual_identity_v0_1.json`;
2. set `style_direction = hybrid`;
3. declare `reference_calibration_role` in each manifest, using
   `docs/VISUAL_PRIOR_PREMIUM_DECKS_AND_AGENTIC_WORKFLOW_REVIEW_2026_06_01.md`;
4. output z-layer/depth-token metadata in each manifest;
5. reject card-box panel framing;
6. create B1 first as the actual brand stress test;
7. then rebuild B2-B4 to match B1, instead of treating previous bridge renders
   as final.

## Current Asset Status

- `output/premium_bridge_scenes/`: logic wireframes only.
- `output/premium_bridge_rebuild_scenes/`: useful bridge candidates, but not
  the locked visual identity.
- `output/visual_identity_style_frames/`: selected identity direction for the
  next real B1-B4 rebuild.

# Visual V9 Provenance-Document Grammar

Date: 2026-06-01

Purpose: decide how v9 platform/resource slides should relate to the dark
result-slide grammar used for Fig1/Fig2/Fig3/Fig6.

## Decision

Use a separate light provenance-document grammar for v9 platform/resource
slides.

Rationale:

- Fig1/Fig2/Fig3/Fig6 are result slides. They benefit from a dark scientific
  scene with a source-derived proof object and an interpretation layer.
- Fig4/Fig5 are release-state and provenance-boundary slides. They should feel
  like an evidence paper, not another biological result.
- v9 has higher overclaim risk because public-bulk metadata, data-file mirrors,
  organoid diagnostics, multispecies pilots, and single-cell blockers are at
  different readiness levels.
- A lighter document grammar makes readiness boundaries and pending gates more
  legible than an atmospheric dark scene.

## Grammar

| Layer | Role | V9 implementation |
|---|---|---|
| Z0 document canvas | Evidence-paper surface | off-white paper texture, faint row/column rules |
| Z1 provenance rail | Release measurement layer | source/task/check or metadata/mirror/package rail |
| Z2 evidence surface | Audited proof object | source-derived Fig4/Fig5 proof crops on paper proof plates |
| Z3 interpretation | Reader-facing release claim | concise title, callout, and status blocks |
| Z4 trust/status | Claim boundary | metadata preview, not payload-frozen; platform status, not final result |
| Z5 motion/focus | One release path | one arrow or focus window around the active gate |

## Generated Outputs

Generator:

- `scripts/build_v9_provenance_document_scenes.py`

Platform/resource scene:

- `output/premium_v9_document_scenes/v9_platform_provenance_document_scene_plate.png`
- `output/premium_v9_document_scenes/v9_platform_provenance_document_scene.png`
- `output/premium_v9_document_scenes/v9_platform_provenance_document_scene.pdf`
- `output/premium_v9_document_scenes/v9_platform_provenance_document_overlay_spec.json`
- `output/premium_v9_document_scenes/v9_platform_provenance_document_scene_manifest.json`

Public-bulk boundary scene:

- `output/premium_v9_document_scenes/v9_public_bulk_boundary_document_scene_plate.png`
- `output/premium_v9_document_scenes/v9_public_bulk_boundary_document_scene.png`
- `output/premium_v9_document_scenes/v9_public_bulk_boundary_document_scene.pdf`
- `output/premium_v9_document_scenes/v9_public_bulk_boundary_document_overlay_spec.json`
- `output/premium_v9_document_scenes/v9_public_bulk_boundary_document_scene_manifest.json`

Review sheet:

- `output/premium_v9_document_scenes/v9_provenance_document_contact_sheet.png`
- `output/premium_v9_document_scenes/v9_provenance_document_contact_sheet_manifest.json`

## Scene Verdicts

| Scene | Visible claim | Boundary text | Verdict |
|---|---|---|---|
| v9 platform/resource | `V9 is a staged evidence resource` | `Platform status; not a final benchmark result` | PASS as resource/status slide |
| Public-bulk boundary | `Public bulk is metadata-ready` | `Metadata preview only; not payload-frozen` | PASS as release-boundary slide |

## What Changed From Result Slides

- Background moved from dark scientific atmosphere to light evidence paper.
- Emphasis moved from biological interpretation to release boundary.
- Focus windows highlight readiness/pending gates, not performance regions.
- Visible text avoids implying that v9 extension tracks are released or frozen.

## Watch-Outs

- Do not use these as main biological result figures.
- Do not remove the metadata/payload boundary language.
- Do not make the v9 platform slide look like a SaaS dashboard.
- If final deck production uses slide-native objects, preserve the light
  document grammar and keep Z4 status labels visible.

## Recommendation

Use the dark grammar for result slides:

- Fig1 tissue transfer;
- Fig2 pathway summaries;
- Fig3 model comparison;
- Fig6 organoid biology checks.

Use the light provenance-document grammar for v9 resource slides:

- Fig4 platform/resource overview;
- Fig5 public-bulk release boundary.

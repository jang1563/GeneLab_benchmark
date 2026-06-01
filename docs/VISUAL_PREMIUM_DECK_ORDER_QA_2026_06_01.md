# Visual Premium Deck-Order QA

Date: 2026-06-01

Purpose: review the current premium slide-scene candidates in deck order,
including the transition from biological result slides to v9 resource/status
slides.

Contact sheet:

- `output/premium_deck_review/premium_slide_deck_order_contact_sheet.png`
- `output/premium_deck_review/premium_slide_deck_order_contact_sheet_manifest.json`

Generator:

- `scripts/render_premium_slide_deck_contact_sheet.py`

## Reviewed Order

| Order | Slide | Grammar | Visible claim | Verdict |
|---|---|---|---|---|
| 1 | Fig1 tissue transfer | dark result scene | `Some tissues transfer` | PASS |
| 2 | Fig2 pathway summaries | dark result scene | `Pathways suppress unwanted labels` | PASS |
| 3 | Fig3 model comparison | dark result scene | `Scale alone does not transfer` | PASS |
| 4 | Fig6 organoid checks | dark result scene with caution layer | `Organoids add biology checks` | PASS |
| 5 | v9 platform/resource | light provenance-document scene | `V9 is a staged evidence resource` | PASS |
| 6 | public-bulk boundary | light provenance-document scene | `Public bulk is metadata-ready` | PASS |

## Verdict

The grammar transition works.

The first four slides read as scientific result/extension slides. The final two
slides read as resource and release-boundary slides. This is the right
separation for the v1-v9 manuscript/deck narrative because v9 should not look
like another finalized biological benchmark result.

## What Is Working

- Dark result slides share a coherent visual family.
- Fig6 stays in the result-slide family but carries the strongest caution
  layer, which is appropriate for organoid diagnostics.
- v9 switches to a light document canvas, making provenance and release status
  visually primary.
- The public-bulk boundary slide makes `metadata-ready` and `data-file mirror
  pending` separable at a glance.
- No visible slide text exposes internal workflow language or raw file paths.

## Remaining Watch-Outs

- v9 platform proof text is dense; keep this slide as a resource/status slide,
  not a stand-alone scientific result.
- Final deck production should rebuild all Z3-Z5 text, focus windows, and
  arrows as editable slide-native objects.
- If slides are presented in a large auditorium, the v9 platform proof crop may
  need one additional zoomed variant focused only on public-bulk readiness.
- Fig1 remains the only result slide with elliptical focus rings; acceptable
  for now, but compare a corner-window variant before final export.

## Production Decision

Proceed with two grammars:

1. Dark result grammar for Fig1/Fig2/Fig3/Fig6.
2. Light provenance-document grammar for v9 platform and public-bulk boundary.

This keeps the deck visually premium while preserving scientific and release
boundaries.

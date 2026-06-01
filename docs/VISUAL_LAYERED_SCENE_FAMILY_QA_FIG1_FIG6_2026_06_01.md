# Visual Layered Scene Family QA: Fig1-Fig6

Date: 2026-06-01

Purpose: review whether the first four layered slide-scene pilots form a
coherent high-profile deck family.

Contact sheet:

- `output/premium_slide_scenes/layered_scene_family_contact_sheet.png`
- `output/premium_slide_scenes/layered_scene_family_contact_sheet_manifest.json`

Generator:

- `scripts/render_layered_scene_contact_sheet.py`

## Reviewed Scenes

| Scene | Claim | Verdict |
|---|---|---|
| Fig1 tissue transfer | `Some tissues transfer` | PASS |
| Fig2 pathway summaries | `Pathways suppress unwanted labels` | PASS |
| Fig3 model comparison | `Scale alone does not transfer` | PASS |
| Fig6 organoid biology checks | `Organoids add biology checks` | PASS |

## Family-Level Verdict

Fig1/Fig2/Fig3/Fig6 now read as one draft slide family.

The shared grammar is clear:

- dark scientific canvas;
- quiet measurement lane;
- source-derived proof object with depth;
- left-side reader-facing interpretation;
- small trust/status line;
- one motion or focus device.

The family avoids the two main failure modes identified earlier:

- it does not look like a set of pasted white plots on blank slides;
- it does not look like dashboard cards or decorative card boxes.

## Remaining Tensions

| Issue | Severity | Decision |
|---|---|---|
| Fig1 uses ellipses while Fig2/Fig3/Fig6 use corner windows | Low | Accept for now; tissue clusters justify ellipses, but compare one Fig1 corner-window variant later |
| Fig6 is denser than Fig1-Fig3 | Medium-low | Accept; it carries source footprint plus diagnostics and has the strongest caution layer |
| Bottom Z1 labels are thumbnail-small | Low | Accept; they are measurement atmosphere, not primary content |
| Right-side evidence layout repeats across result slides | Medium-low | Accept for draft; vary layout later for section transitions or platform slides |

## Text And Claim Discipline

Visible claims are now reader-facing:

- Fig1: tissue transfer is heterogeneous.
- Fig2: pathway summaries reduce unwanted labels.
- Fig3: scale alone does not solve transfer.
- Fig6: organoids add biology checks.

The strongest caution layer is on Fig6:

- `Draft extension: source factors remain coupled`.

This is correct. Organoid slides have higher overclaim risk than the mouse
bulk benchmark slides.

## Production Recommendation

Move forward with these four as draft main deck candidates.

Before final deck export:

1. rebuild Z3-Z5 as editable slide-native objects;
2. recheck at full screen and presenter-view thumbnail scale;
3. run grayscale/color-vision QA on the final slide exports;
4. decide whether Fig1 should keep the ellipse treatment or switch to a
   corner-window variant for family consistency.

Next design block:

- v9 platform/resource slides. These may need a lighter provenance-document
  grammar rather than the same dark result-slide grammar.

# Visual External Text and Color QA

Date: 2026-06-01

Purpose: confirm that current external-facing figure text avoids internal
operating language, and record the first grayscale/color-vision QA pass for the
main figure candidates.

## Text Audit

Checked files:

- `scripts/build_premium_visuals.py`
- `docs/VISUAL_MANUSCRIPT_VARIANT_REVIEW_2026_06_01.md`
- `docs/VISUAL_EXTERNAL_CAPTION_PACK_2026_06_01.md`

High-risk phrases checked:

- `Source:`, `Sources:`
- `Claim boundary`, `Boundary:`
- `Coverage:`, `Coverage matters`
- `Do not read`
- `What to show`
- `Selected release path`
- `Allowed`, `Prohibited`
- `metadata-only`
- `not a leaderboard`
- `Release boundary`
- `draft main`
- `caption boundary`

Result:

- No matches in the current figure generator and external caption/review docs.
- Source-file paths and provenance labels remain appropriate in manifests and
  audit files, but they should not be used as manuscript figure text.

## Figure Text Changes

Removed from visible figures:

- bottom source/provenance footers;
- claim-boundary phrasing;
- release-path and package-decision language;
- allowed/prohibited wording;
- text-heavy coverage warnings.

Replaced with reader-facing scientific wording:

- "evidence depth differs by data layer";
- "metadata indexed";
- "pilot evidence";
- "data pending";
- "evaluation scope" and "task set";
- "descriptive association" and "no fitted trend" where needed to read a plot.

Additional visible-text cleanup after enlarged inspection:

- deck-style Fig3 now uses "Evaluation scope" and "task set" language and
  removes a lower-right text collision;
- deck-style Fig4/Fig6 now use "extension datasets" and "biology-check
  dataset" instead of compact "surface" wording;
- Fig6 manuscript expands shorthand labels into "direction" and "rank corr.".

## Caption Pack

New external caption pack:

- `docs/VISUAL_EXTERNAL_CAPTION_PACK_2026_06_01.md`

Captions were drafted for:

- Fig1 tissue-specific cross-mission generalization;
- Fig2 pathway summaries and coupled-label prediction;
- Fig3 model comparison across evaluation task sets;
- Fig6 human neural organoid biology-check dataset.

## Grayscale And Color-Vision QA

New script:

- `scripts/render_premium_color_qa.py`

Generated outputs:

- `output/premium_figures/color_qa/premium_color_qa_contact_sheet.png`
- `output/premium_figures/color_qa/premium_color_qa_manifest.csv`
- `output/premium_figures/color_qa/premium_color_qa_manifest.json`

QA modes:

- original;
- grayscale;
- approximate deuteranopia;
- approximate protanopia.

Visual inspection summary:

- Fig2 manuscript: passes; group differences are still readable by position,
  direct labels, and paired-line structure.
- Fig3 manuscript: passes; model ranking and tissue deltas are readable without
  relying only on color.
- Fig6 manuscript: passes after the manuscript relayout; the visible text now
  uses "dataset", "direction", and "rank corr." instead of compact internal
  shorthand, and bar heights, marker shape, row labels, and direct values carry
  the signal under grayscale and color-vision transforms.
- Fig1: passes after adding marker-shape groups and direct high/mid/near-chance
  labels; color emphasis is no longer the only cue for tissue ranking.

## Current External-Facing Set

Recommended manuscript-facing figure set after this pass:

- `premium_fig1_tissue_transfer_hierarchy`
- `premium_fig2_pathway_rescue_manuscript`
- `premium_fig3_model_tier_comparison_manuscript`
- `premium_fig6_organoid_biology_check_manuscript`

Use Fig4/Fig5 as resource/status visuals only if the manuscript needs a data
package or platform section. They are now worded more externally, but they are
not primary biological-result figures.

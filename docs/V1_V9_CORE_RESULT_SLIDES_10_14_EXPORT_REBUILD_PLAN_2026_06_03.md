# v1-v9 core result slides 10-14 export/rebuild plan

Date: 2026-06-03

## Anchor

This follows the slide 6 feature-layer bridge. The remaining visual risk in the
24-slide deck spine is that v9 extension materials look more polished than the
v1-v7 result core.

This note fixes the next work block for slides 10-14:

- slide 10: v4 hardening,
- slide 11: temporal and RRRM lessons,
- slide 12: negative results,
- slide 13: biological interpretation,
- slide 14: human translation.

## Current Verdict

Do not paste the old HTML figures directly into the final deck.

All five slides should be treated as premium rebuild or static-export targets.
The existing HTML figures and JSON summaries are source material, not final
deck assets.

## Slide 10: v4 Hardening

Purpose:

- show that the benchmark result is not a narrow model or tissue cherry-pick.

Source assets:

- `v4/figures/html/Fig1_benchmark.html`
- `v4/figures/html/Fig2_ablation.html`
- `v4/figures/html/Fig5_generalize.html`
- `v4/evaluation/M1_summary.json`
- `docs/CANONICAL_RESULTS_V7_1.md`

Source claim:

- v4 covers 8 tissues x 8 classifiers x 4 feature types = 256 evaluations.
- 40 evaluations are significant at p<0.05.
- 6/8 tissues have at least one significant configuration.
- Best-row table should use the canonical v7.1 wording and subset notes.

Recommended visual:

- rebuild, not screenshot;
- compact 8-tissue x feature/method surface plus 2-3 canonical callouts;
- use this as a "hardening surface", not another leaderboard.

Claim boundary:

- do not imply every best-row tissue is significant;
- keep canonical/raw discrepancy out of visible slide text unless the slide is
  explicitly about reconciliation.

## Slide 11: Temporal And RRRM Lessons

Purpose:

- show that benchmark value includes temporal confounds, recovery patterns, and
  single-cell lessons.

Source assets:

- `v2/figures/Fig1_temporal.html`
- `v2/figures/Fig4_scrna_summary.html`
- `v2/evaluation/V2_RESULTS_SUMMARY.md`

Source claim:

- ISS-T versus LAR differences are strongly confounded by preservation method.
- RR-6 and RR-8 LAR samples project closer to baseline than ISS-T in
  descriptive recovery analyses.
- RRRM single-cell summaries should be presented as lessons, not as the main
  benchmark score surface.

Recommended visual:

- two-lane slide:
  - temporal/preservation lesson,
  - RRRM single-cell lesson;
- avoid a dense table; use one confound marker and one cell-type signal object.

Claim boundary:

- preservation artifact is a confound lesson, not a failure;
- recovery projection is descriptive, not held-out validation evidence.

## Slide 12: Negative Results

Purpose:

- make negative/failure-mode results feel like benchmark value, not weak work.

Source assets:

- `v3/figures/v3_Fig2_spatial_overview.html`
- `v3/figures/v3_Fig5_fm_comparison.html`
- `v3/README.md`

Source claim:

- spatial brain result is negative: section-level AUROC=0.139 and
  animal-level AUROC=0.444.
- brain bulk RNA-seq overfits at n=6 and should be described carefully.
- UCE/scFoundation do not automatically improve spaceflight detection.

Recommended visual:

- rebuild as a "failure modes are informative" slide;
- left: spatial/batch/brain negative object;
- right: FM negative/control object;
- center caption: negative results define benchmark boundaries.

Claim boundary:

- do not say brain has no biology; say no detectable spaceflight signal in this
  analyzed setting;
- do not make a universal foundation-model failure claim.

## Slide 13: Biological Interpretation

Purpose:

- show that the benchmark can connect score surfaces to biology without
  overstating causality.

Source assets:

- `v5/figures/html/Fig7_immune_signaling.html`
- `v5/figures/html/Fig8_metabolism_drugs.html`
- `v5/figures/html/Fig9_consensus_panel.html`
- `v5/evaluation/consensus_biomarker_panel.json`
- `v5/evaluation/cross_organ_signaling.json`
- `v5/evaluation/drug_targets.json`

Source claim:

- consensus biomarker panel has 20 genes scored from 1,919 candidates.
- cross-organ signaling uses 111 ligand-receptor pairs and 56 directed tissue
  pairs.
- druggability/target links are interpretation and triage evidence, not causal
  validation.

Recommended visual:

- rebuild as a biology interpretation triad:
  - immune/signaling,
  - metabolism/drug-target,
  - consensus marker panel;
- show one dominant biological map, not three old HTML screenshots.

Claim boundary:

- interpretation layer only;
- no causal proof and no treatment recommendation.

## Slide 14: Human Translation

Purpose:

- show translation potential while preserving the "not clean gene transfer"
  boundary.

Source assets:

- `v6/figures/html/v6_Fig10_human_translation.html`
- `v6/figures/html/v6_FigS8_translation_detail.html`
- `v6/evaluation/V6_B_pathway_conservation.json`
- `v6/evaluation/V6_C_cross_species_transfer.json`
- `v6/evaluation/V6_E_tf_conservation.json`
- `v6/evaluation/V6_F_drug_target_validation.json`

Source claim:

- pathway conservation has mean rho=0.285 across 5 analyzed tissues.
- TF conservation mean rho is low at 0.0298, so this should not be framed as a
  broad regulatory transfer success.
- drug-target validation separates validated/promising/detected/untranslatable
  tiers.

Recommended visual:

- rebuild as a translation ladder:
  - pathway-level partial conservation,
  - weak/partial gene or TF transfer,
  - target evidence tiers;
- avoid detailed drug lists in the main slide.

Claim boundary:

- pathway/target evidence, not clean gene-level transfer;
- no clinical or countermeasure recommendation.

## Priority Order

P0:

- slide 10 v4 hardening rebuild,
- slide 14 human translation rebuild.

P1:

- slide 11 temporal/RRRM lessons,
- slide 12 negative results,
- slide 13 biological interpretation.

Reason:

- slide 10 protects the benchmark from cherry-pick criticism;
- slide 14 protects the translation claim boundary;
- slides 11-13 add nuance and biological depth after the core score story is
  stable.

## Next Implementation Step

Create a small slides 10-14 rebuild matrix/contact sheet using existing source
assets as thumbnails and explicit "rebuild/export" status labels. Then rebuild
slide 10 first.

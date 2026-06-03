# Slide 13 Biological Interpretation Scene Review

Date: 2026-06-03

## Role

Slide 13 is the v5 biological interpretation layer for the premium core-result
deck. Its job is to show how benchmark hits become follow-up biological
hypotheses across immune context, metabolic context, and biomarker/target
triage.

This is not a treatment, countermeasure, causal mechanism, or validated
translation claim.

## Source Files Read

- `v5/evaluation/immune_deconv_*.json`
- `v5/evaluation/cross_organ_signaling.json`
- `v5/evaluation/metabolic_flux_*.json`
- `v5/evaluation/consensus_biomarker_panel.json`
- `v5/evaluation/drug_targets.json`
- `v5/figures/html/Fig7_immune_signaling.html`
- `v5/figures/html/Fig8_metabolism_drugs.html`
- `v5/figures/html/Fig9_consensus_panel.html`
- `README.md`

## Evidence Extracted

Immune context:
- 8 tissues and 105 immune cell-type tests were read from the deconvolution
  JSON files.
- 10 FDR<0.05 cell-type shifts were detected.
- The strongest FDR-supported examples were thymus CD8 T cells, skin vessels,
  skin basophils, and kidney B-derived cells.
- `cross_organ_signaling.json` currently reports 1 active tissue edge and 1
  active ligand-receptor pair: thymus to colon, F12 to Gp9.

Metabolic context:
- 6 tissues had successful E-Flux/iMM1865 outputs.
- Mean gene coverage was 89.8%, with 84.6% minimum coverage.
- Each tissue reported 103 pathway/subsystem contexts.
- The visible slide uses simplified subsystem labels to avoid dense pathway
  jargon in the main figure.

Biomarker/target triage:
- 1,919 genes were scored into a 20-gene consensus panel.
- Mean panel AUROC across 8 tissues was 0.682; best tissue was gastrocnemius at
  0.806.
- 10 of 20 panel genes had drug links in the current JSON.
- `drug_targets.json` reports 271/834 mapped human genes with DGIdb
  interactions, 1,284 tier-1 approved links, 200 tier-3 preclinical links, and
  26 ChEMBL targets.

## Source-Consistency Flags

Two legacy HTML narrative sentences were not reused in the slide:

- `Fig7_immune_signaling.html` says 62 active L-R pairs, but the current JSON
  reports 1 active edge and 1 active pair.
- `Fig9_consensus_panel.html` says every panel gene has known drug interactions,
  but the current JSON supports 10/20 drug-linked panel genes.

Resolution: the slide follows current JSON values and records these flags in
the manifest.

## Visual Decision

The slide avoids table-like layouts inside the figure. The main figure is a
layered PNG scene plus editable-style interpretation shell:

- Z0 canvas: dark field, grid, orbital arcs.
- Z1 measurement layer: tissue nodes, flux bars, AUROC strip.
- Z2 evidence surfaces: three source-derived panels.
- Z3 interpretation layer: a restrained bottom trajectory from measurement to
  triage.
- Z4 trust layer: footer caveat and source label.

Visible wording was reduced from jargon-heavy labels:
- `FDR cell-type shifts` -> `significant shifts`
- `active L-R pair` -> `signaling link`
- `drug-linked` -> `target-linked`
- long subsystem names -> simplified pathway labels

## Outputs

- Scene plate:
  `output/premium_core_result_slides/slide13_biological_interpretation_v0_1/slide13_biological_interpretation_scene_plate.png`
- Rendered preview:
  `output/premium_core_result_slides/slide13_biological_interpretation_v0_1/slide13_biological_interpretation_rendered_preview.png`
- Grayscale QA:
  `output/premium_core_result_slides/slide13_biological_interpretation_v0_1/qa/slide13_biological_interpretation_rendered_preview_grayscale.png`
- Source summary:
  `output/premium_core_result_slides/slide13_biological_interpretation_v0_1/slide13_biological_interpretation_source_summary.json`
- Manifest:
  `output/premium_core_result_slides/slide13_biological_interpretation_v0_1/slide13_biological_interpretation_manifest.json`
- QA:
  `output/premium_core_result_slides/slide13_biological_interpretation_v0_1/slide13_biological_interpretation_qa.json`

## QA Result

Automatic checks passed:
- 3840 x 2160 RGB PNG outputs.
- JSON QA/manifest parse cleanly.
- Python builder compiles.
- Visible text count: 51 words, under the 58-word budget.

Manual visual review:
- First render had an immune bubble overflow, a too-heavy bottom trajectory, and
  residual `FDR/L-R/drug-linked` jargon.
- Second render fixed panel containment, softened the trajectory, simplified
  labels, and passed rendered/grayscale inspection.

## Claim Boundary

Recommended spoken/written caption:

> v5 converts benchmark signal into biological hypotheses: immune shifts,
> pathway-level metabolic context, and target-linked biomarker triage. These are
> ranked follow-up signals, not causal mechanisms or countermeasure
> recommendations.

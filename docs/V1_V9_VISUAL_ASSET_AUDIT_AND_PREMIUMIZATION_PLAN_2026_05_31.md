# v1-v9 Visual Asset Audit And Premiumization Plan

Date: 2026-05-31

Purpose: careful first-pass audit of visual assets from v1 through current v9,
with a premiumization plan for future sessions.

Companion standard:

- `docs/VISUAL_PREMIUM_QUALITY_STANDARD_2026_05_31.md`
- `docs/VISUAL_LAYERED_SCENE_BLUEPRINT_2026_06_01.md`

Related prior docs:

- `docs/V1_V9_PRESENTATION_AND_MANUSCRIPT_MASTER_OUTLINE_2026_05_31.md`
- `docs/V1_V9_SLIDE_ASSET_MANIFEST_2026_05_31.md`
- `docs/V1_V9_FIGURE_QA_PASS1_2026_05_31.md`

## Executive Verdict

The current asset collection is scientifically rich but not yet visually
premium.

The repository contains many useful figure prototypes, but the visual system is
fragmented across versions:

- v1-v7 are mostly D3 HTML figures that require live/remote D3 and are not
  static-publication assets yet.
- v8 has static PNG/PDF figure pairs, but several need design cleanup,
  cropping, or claim-boundary labels.
- v9 has very strong evidence tables/reports but almost no final visual
  figures; v9 should be visualized as architecture, status matrices, and
  extension-track dashboards.

Best strategy:

> Use existing figures as scientific/design references, then regenerate the
> main deck/manuscript figures in a unified static premium style from canonical
> source files.

Slide-deck refinement:

> For the premium presentation, turn source-verified proof figures into
> layered scenes: Z0-Z2 carry canvas, measurement context, and audited evidence;
> Z3-Z5 carry editable interpretation, trust/status language, and one focus
> movement.

JK decision:

> If a current visual looks weak, we should not spend energy polishing it in
> place. Rebuild the important visuals from scratch when that produces a
> stronger deck/paper.

## Asset Count

Focused figure directories scanned:

- `figures/`
- `v2/figures/`
- `v3/figures/`
- `v4/figures/`
- `v5/figures/`
- `v6/figures/`
- `v7/figures/`
- `v8/figures/`

Counts:

| Asset type | Count | Meaning |
|---|---:|---|
| HTML | 47 | D3 figure outputs, mostly interactive/static hybrid |
| PNG | 5 | v8 static raster figures |
| PDF | 5 | v8 static vector/page figures |
| Python figure scripts | 19 | v4-v8 figure generation sources |
| Total focused figure assets | 76 | Excluding v9 evidence tables |

v9 evidence/visualization candidates:

| v9 artifact type | Count in reports/scaffold areas | Meaning |
|---|---:|---|
| CSV | 303 | status tables, baselines, audits, diagnostics |
| JSON | 507 | metrics, manifests, run records, audit data |
| Markdown | 61 | review notes and decision records |
| Static figure files | 0 found in v9 report/scaffold areas | needs new design work |

## Cross-Cutting Technical Finding

All selected v1-v7 HTML figures depend on remote D3:

```html
<script src="https://d3js.org/d3.v7.min.js"></script>
```

Impact:

- They are not offline-stable.
- QuickLook thumbnails do not render them correctly.
- They are not immediately manuscript-ready.
- They should not be screenshot through QuickLook.

Resolution options:

1. Vendor D3 locally and browser-export each HTML figure.
2. Use a browser environment with internet access and export SVG/PNG.
3. Regenerate the main figures from canonical data sources in a static plotting
   system.

Preferred:

- Regenerate the main figures.
- Keep HTML figures as prototypes and backup references.

## Version-Level Asset Review

### Root / v1 Core Figures

Files:

- `figures/fig1_tissue_hierarchy.html`
- `figures/fig2_pathway_mechanism.html`
- `figures/fig3_model_comparison.html`
- `figures/fig4_validation.html`
- `figures/figS1_data_controls.html`
- `figures/figS2_nes_multidb.html`
- `figures/figS3_transfer_detail.html`
- `figures/figS4_temporal_bio.html`
- `figures/figS5_pipeline_comparison.html`

Strengths:

- Strong story structure already exists: tissue hierarchy, pathway mechanism,
  model comparison, validation, controls.
- The root figures map well to the main paper/deck narrative.
- Several files include built-in SVG download logic.

Premium gaps:

- Remote D3 dependency.
- HTML/UI button not appropriate in final export.
- Style is version-local, not unified with v2-v9.
- Needs canonical-value check against `docs/CANONICAL_RESULTS_V7_1.md` and
  `evaluation/RESULTS_SUMMARY.md`.

Decision:

- Regenerate after numeric verification for main story use. Browser export is
  acceptable only for backup/reference use.
- Highest priority for premium rebuild.

Premium priority:

- P0: `fig1_tissue_hierarchy.html`
- P0: `fig2_pathway_mechanism.html`
- P0: `fig3_model_comparison.html`
- P1: `fig4_validation.html`
- P1/P2: supplementary figures

### v2 Figures

Files:

- `v2/figures/Fig1_temporal.html`
- `v2/figures/Fig2_crossspecies.html`
- `v2/figures/Fig3_pbmc_celltype.html`
- `v2/figures/Fig4_scrna_summary.html`
- `v2/figures/F1_temporal_heatmap.html`
- `v2/figures/F1_celltype_pathway_heatmap.html`
- `v2/figures/F2A_composition.html`
- `v2/figures/F2B_celltype_nes_heatmap.html`
- `v2/figures/F2C_loao_auroc.html`
- `v2/figures/F2D_crossspecies_scatter.html`
- `v2/figures/E1_crossspecies_scatter.html`
- `v2/figures/E2_duration_conservation.html`
- `v2/figures/E3_cfrna_origin.html`

Strengths:

- Rich temporal, cross-species, PBMC, and RRRM-1 single-cell content.
- Strong deck material for "extensions reveal confounds and cell-type signal."
- Good source coverage from `v2/evaluation/V2_RESULTS_SUMMARY.md`.

Premium gaps:

- Many figures are too detailed for main deck.
- Long embedded lines indicate self-contained D3/data-heavy HTML, but still
  remote-D3 dependent.
- v2 story can sprawl; only a small number should enter the main deck.

Decision:

- Select 2-3 representative panels for premium rebuild.
- Keep the rest as backup/supplement.

Premium priority:

- P1: temporal/preservation/recovery summary.
- P1: PBMC cell-type pathway response.
- P1: RRRM-1 scRNA-seq summary.
- P2: detailed heatmaps and cross-species scatter.

### v3 Figures

Files:

- `v3/figures/v3_Fig1_multispecies_concordance.html`
- `v3/figures/v3_Fig2_spatial_overview.html`
- `v3/figures/v3_Fig3_rrrm2_scrna.html`
- `v3/figures/v3_Fig4_extensions.html`
- `v3/figures/v3_Fig5_fm_comparison.html`

Strengths:

- v3 has compact figure count and clear extension topics.
- Spatial negative result and RRRM-2 PBMC/NK positive signal are good premium
  storytelling material.
- FM comparison reinforces the benchmark message.

Premium gaps:

- Same remote D3 issue.
- Needs claim-boundary wording for negative results.
- Cross-species directionality should not be overgeneralized.

Decision:

- Rebuild one "failure modes and extension signals" figure from v3.
- Use spatial negative + RRRM-2 positive + FM extension as a balanced trio.

Premium priority:

- P1: `v3_Fig2_spatial_overview.html`
- P1: `v3_Fig3_rrrm2_scrna.html`
- P2: `v3_Fig5_fm_comparison.html`

### v4 Figures

Files:

- `v4/figures/html/Fig1_benchmark.html`
- `v4/figures/html/Fig2_ablation.html`
- `v4/figures/html/Fig3_consensus.html`
- `v4/figures/html/Fig4_network.html`
- `v4/figures/html/Fig5_generalize.html`
- `v4/figures/html/Fig6_multiplatform.html`
- `v4/figures/html/FigS1_pertissue.html`
- `v4/figures/html/FigS2_delong.html`
- `v4/figures/html/FigS3_wgcna_overview.html`
- `v4/figures/html/FigS4_shap_cross.html`
- `v4/figures/html/FigS5_controls.html`
- matching Python scripts in `v4/figures/`

Strengths:

- v4 is the strongest "not cherry-picked" evidence layer.
- Matching Python scripts exist, which is valuable for regeneration.
- Figures cover benchmark grid, ablations, networks, generalizability, controls.

Premium gaps:

- Remote D3 outputs.
- Current v4 canonical public table differs from raw
  `v4/evaluation/M1_summary.json` best-row extraction.
- Some figures are likely too dense for deck main slides.

Decision:

- Do not finalize any v4 best-method visual until canonical/raw table decision
  is resolved.
- Rebuild a clean v4 hardening figure from one selected source-of-truth table.

Premium priority:

- P0: v4 hardening grid/table after numeric reconciliation.
- P1: ablation/generalizability.
- P2: WGCNA/network/control details.

### v5 Figures

Files:

- `v5/figures/html/Fig7_immune_signaling.html`
- `v5/figures/html/Fig8_metabolism_drugs.html`
- `v5/figures/html/Fig9_consensus_panel.html`
- `v5/figures/html/FigS6_immune_detail.html`
- `v5/figures/html/FigS7_metabolic_detail.html`
- matching Python scripts in `v5/figures/`

Strengths:

- Provides biological interpretation: immune, TF/signaling, metabolism,
  drug-target, biomarker panel.
- Useful for "prediction to mechanism" narrative.

Premium gaps:

- Dense mechanism panels can distract from benchmark story.
- These should not be presented as causal proof.
- Same D3/export issue.

Decision:

- Use one compact premium figure/table for v5/v6 combined "biology and
  translation" rather than many mechanism slides.

Premium priority:

- P1: `Fig7_immune_signaling.html`
- P1: `Fig8_metabolism_drugs.html`
- P2: `Fig9_consensus_panel.html`

### v6 Figures

Files:

- `v6/figures/html/v6_Fig10_human_translation.html`
- `v6/figures/html/v6_FigS8_translation_detail.html`
- matching Python scripts in `v6/figures/`

Strengths:

- Strong compact source for partial translation: gene weak, pathway/target
  partial, drug target tiers.
- Clear candidate for one manuscript figure or extended-data figure.

Premium gaps:

- Remote D3.
- Needs careful "partial translation" framing.
- Should not imply clean human validation of all mouse findings.

Decision:

- Rebuild one clean translation summary panel.

Premium priority:

- P1: `v6_Fig10_human_translation.html`
- P2: `v6_FigS8_translation_detail.html`

### v7 Figures

Files:

- `v7/figures/html/v7_methods_comparison.html`
- `v7/figures/html/v7_unified_signal_hierarchy.html`

Strengths:

- Captures advanced-method comparison and broader ecosystem linkage.
- Useful for showing v7/v7.1 public consolidation and method hierarchy.

Premium gaps:

- Remote D3.
- `v7_unified_signal_hierarchy` includes SpaceOmicsBench; since JK clarified
  SpaceOmicsBench is self/sibling, the visual framing must not imply external
  competition.
- Need subset labels and claim-boundary text.

Decision:

- Use cautiously.
- Rebuild if included in main deck; otherwise keep as backup.

Premium priority:

- P1/P2 depending on final deck story.

### v8 Figures

Files:

- `v8/figures/Figure1_Species_Transfer.png`
- `v8/figures/Figure1_Species_Transfer.pdf`
- `v8/figures/Figure2_Stressor_Decomposition.png`
- `v8/figures/Figure2_Stressor_Decomposition.pdf`
- `v8/figures/Figure3_Mars_Extrapolation.png`
- `v8/figures/Figure3_Mars_Extrapolation.pdf`
- `v8/figures/Figure4_Countermeasure_Pareto.png`
- `v8/figures/Figure4_Countermeasure_Pareto.pdf`
- `v8/figures/Figure5_Causal_DAG.png`
- `v8/figures/Figure5_Causal_DAG.pdf`
- `v8/figures/generate_main_figures.py`

Strengths:

- Static PNG/PDF pairs already exist.
- `Figure1_Species_Transfer.png` is visually readable and story-relevant.
- The figure script provides a path to regenerate.

Visual findings:

- Figure 1: usable, high-resolution, good three-panel layout, but needs
  hypothesis-only boundary text on the slide.
- Figure 2: readable, but original rightmost panel is blank; revise/crop before
  main use.
- Figure 3: visually readable, but high overclaim risk because Mars
  extrapolation can look predictive; likely backup unless reframed.
- Figure 4: has strong hypothesis content, but panel C is empty and label
  clipping/large whitespace reduce polish; backup or regenerate.
- Figure 5: causal DAG is readable but too sparse and visually raw for premium
  main deck; good conceptual reference, not final art.

Decision:

- Reuse Figure 1 as a near-term slide candidate.
- Regenerate/crop Figures 2 and 4 if used.
- Treat Figures 3 and 5 as backup/reference unless redesigned.

Premium priority:

- P1: Figure 1.
- P2: Figure 2 after crop/regeneration.
- P2/P3: Figures 3-5 with strong claim-boundary redesign.

### v9 Visual State

No static v9 figure set was found.

Relevant evidence tables/reports include:

- `v9/reports/bulk_lomo_baseline_summary.csv`
- `v9/source_inventory.csv`
- `v9/source_checksum_audit.csv`
- `v9/reports/public_bulk_alpha_gap_matrix/public_bulk_alpha_gap_matrix.csv`
- `v9/reports/public_bulk_alpha_snapshot_decision/snapshot_decision_summary.csv`
- `v9/human_organoid/task_manifest_index.draft.csv`
- `v9/human_organoid/expression_matrix_audit.draft.csv`
- `v9/human_organoid/signature_reference_audit.draft.csv`
- `v9/multispecies/task_manifest_index.draft.csv`
- `v9/multispecies/interaction_task_manifest_index.draft.csv`
- `v9/sc_spaceflight/asset_inventory_summary.csv`
- `v9/sc_spaceflight/anndata_manifest_draft_summary.csv`
- `v9/sc_spaceflight/metric_spec_summary.csv`
- `v9/sc_spaceflight/obs_var_audit_summary.csv`

Strengths:

- Excellent evidence substrate for status dashboards.
- Strong provenance and release-boundary story.
- Many CSV/JSON reports can feed compact visuals.

Premium gaps:

- No final static v9 visual language yet.
- Needs architecture diagram, status matrix, extension-track table, and blocker
  diagram.

Decision:

- Build v9 visuals from scratch.
- Do not try to repurpose raw CSV/JSON as slides.

Premium priority:

- P0: platform architecture.
- P0: public bulk alpha status matrix.
- P1: organoid/multispecies/single-cell extension matrix.
- P1: single-cell blocker/status diagram.

## Premium Main-Figure Candidate Set

Recommended main premium visuals for the full deck/paper:

| Priority | Figure | Source family | Reuse/regenerate decision |
|---|---|---|---|
| P0 | Mission-held-out benchmark design | root/v1 + v9 platform | New schematic |
| P0 | Tissue hierarchy / thymus vs liver | root `fig1` | Regenerate from canonical values |
| P0 | Pathway rescue and batch resistance | root `fig2` | Regenerate from canonical values |
| P0 | Model tier comparison | root `fig3`, v7 | Regenerate with subset labels |
| P0 | v4 hardening grid | v4 | Regenerate after canonical/raw reconciliation |
| P1 | Temporal/single-cell lessons | v2/v3 | Regenerate compressed composite |
| P1 | Biology/translation | v5/v6 | Regenerate compact composite |
| P1 | v8 translational incubator | v8 Figure 1 plus revised Figure 2 | Hybrid: reuse/regenerate |
| P0 | v9 platform/status dashboard | v9 CSV/JSON/MD evidence | New schematic/table figure |

## Retire Or Backup Candidates

Likely backup/supplement only:

- Dense v2 heatmaps.
- v4 SHAP/WGCNA details.
- v5 metabolic subsystem detail.
- v8 Mars extrapolation unless reframed carefully.
- v8 countermeasure Pareto unless panel C is fixed and hypothesis-only language
  is prominent.
- v7 unified signal hierarchy unless self/sibling ecosystem framing is
  redesigned.

## Immediate Next Sessions

Recommended order:

1. Build one Fig1 layered-scene pilot from the already audited tissue-transfer
   proof object.
2. Reuse the same z-stack grammar for Fig2 and Fig3.
3. Resolve v4 canonical/raw discrepancy and build premium v4 hardening figure.
4. Build v9 platform architecture and status matrix as layered evidence-trail
   scenes, not dashboards.
5. Return to v8 and decide reuse/crop/regenerate for translational figures with
   mandatory hypothesis-only trust labels.
6. Reconcile any remaining numeric source-of-truth issues before promoting new
   proof objects.

## Session Guardrails

- Do not improve aesthetics before verifying the source number.
- Do not use HTML QuickLook screenshots as final assets.
- Do not mix v8 intervention or Mars visuals into benchmark result claims.
- Do not make v9 extension tracks look released/frozen.
- Do not let the deck inherit old per-version styling; the final output should
  look like one coherent product.
- Do not paste a manuscript figure onto a blank slide and call it premium; each
  main slide needs a deliberate layered-scene brief.

# v1-v9 Slide Asset Manifest

Date: 2026-05-31

Purpose: map the proposed v1-v9 presentation to existing figures, source data,
new figure needs, and claim-boundary notes. This is the execution bridge from
content inventory to an actual slide deck.

Companion outline:

- `docs/V1_V9_PRESENTATION_AND_MANUSCRIPT_MASTER_OUTLINE_2026_05_31.md`

Companion QA:

- `docs/V1_V9_FIGURE_QA_PASS1_2026_05_31.md`

## Asset Inventory Summary

Existing visual assets:

- Root benchmark figures: 9 HTML figures in `figures/`.
- v2 figures: 13 HTML figures in `v2/figures/`.
- v3 figures: 5 HTML figures in `v3/figures/`.
- v4 figures: 11 HTML renderings in `v4/figures/html/` plus figure scripts.
- v5 figures: 5 HTML renderings in `v5/figures/html/` plus figure scripts.
- v6 figures: 2 HTML renderings in `v6/figures/html/` plus figure scripts.
- v7 figures: 2 HTML renderings in `v7/figures/html/`.
- v8 figures: 5 PNG/PDF figure pairs in `v8/figures/`.
- v9: mostly CSV/JSON/MD evidence; platform and extension-track slides need
  new schematic/table figures.

Primary caution:

- Existing figures were inventoried but not visually rendered/reviewed in the
  current pass. Any deck build should render or screenshot each selected HTML
  figure before final slide placement.
- QA pass 1 found that the selected HTML figures rely on remote D3 and do not
  render correctly through QuickLook; core benchmark figures should be browser
  exported with local D3 or regenerated as static figures.

## Recommended Long-Talk Slide Manifest

| Slide | Proposed title | Existing asset | Source data / text source | New work needed | Claim boundary |
|---:|---|---|---|---|---|
| 1 | SpaceBio-Bench | New title visual | master outline | Create title slide with mission-held-out/platform thesis | No "first AI benchmark" claim |
| 2 | Why Mission-Held-Out? | New schematic | `docs/V1_PAPER_CONTENT.md`, `evaluation/RESULTS_SUMMARY.md` | Diagram random split vs LOMO split by mission | Mission is independence unit |
| 3 | External Landscape | New positioning matrix | `docs/SIMILAR_PROJECTS_AND_POSITIONING_SCAN_2026_05_31.md`, `docs/LITERATURE_QUERY_LOG_SPACEBIO_BENCH_2026_05_31.md` | 2x2 or table: benchmark vs spaceflight specificity | SOMA/SpaceOmicsBench are self/sibling ecosystem |
| 4 | Project Evolution | New timeline | `docs/PROJECT_SLIDE_CONTENT_INVENTORY_V1_TO_V9_2026_05_31.md` | v1-v9 horizontal timeline | v8/v9 status tiers must be visible |
| 5 | Data And Task Design | `figures/figS1_data_controls.html` or new schematic | `docs/CANONICAL_RESULTS_V7_1.md`, root `tasks/` | Prefer new simplified design figure | Separate full v1-v7 scope from v1/FM subsets |
| 6 | Thymus Beats Liver | `figures/fig1_tissue_hierarchy.html` | `evaluation/RESULTS_SUMMARY.md` | Screenshot/render; maybe redraw static bar/CI plot | Scope as original mouse bulk LOMO/transfer hierarchy |
| 7 | Pathway Rescue | `figures/fig2_pathway_mechanism.html` | `evaluation/RESULTS_SUMMARY.md` | Verify exact kidney/eye values before final caption | "Can rescue selected tasks", not all tissues |
| 8 | Confounder Resistance | `figures/fig2_pathway_mechanism.html` or `figures/figS1_data_controls.html` | `evaluation/RESULTS_SUMMARY.md` | Possibly split from slide 7 into clean gene-vs-pathway panel | Gene features useful for QC; pathways not universally superior |
| 9 | Model Comparison | `figures/fig3_model_comparison.html` | `docs/CANONICAL_RESULTS_V7_1.md`, `evaluation/RESULTS_SUMMARY.md` | Add subset labels for 6-tissue vs best single-tissue rows | No universal FM failure claim |
| 10 | v4 Hardening | `v4/figures/html/Fig1_benchmark.html`, `v4/figures/html/Fig2_ablation.html`, `v4/figures/html/Fig5_generalize.html` | `docs/CANONICAL_RESULTS_V7_1.md`, `v4/evaluation/M1_summary.json` | Decide canonical vs raw best-row table before final slide | Important v4 canonical/raw difference unresolved |
| 11 | Temporal And Single-Cell Lessons | `v2/figures/Fig1_temporal.html`, `v2/figures/Fig4_scrna_summary.html` | `v2/evaluation/V2_RESULTS_SUMMARY.md` | Combine temporal/recovery and RRRM-1 into one summary slide | Preservation artifact is a confound lesson |
| 12 | Negative Results Matter | `v3/figures/v3_Fig2_spatial_overview.html`, `v3/figures/v3_Fig5_fm_comparison.html` | `v3/README.md` | Emphasize negative controls/results as benchmark value | Negative result, not failed project |
| 13 | Biological Interpretation | `v5/figures/html/Fig7_immune_signaling.html`, `v5/figures/html/Fig8_metabolism_drugs.html`, `v5/figures/html/Fig9_consensus_panel.html` | `v5/evaluation/` | Select one clean panel; avoid overcrowding | Interpretation layer, not causal proof |
| 14 | Human Translation | `v6/figures/html/v6_Fig10_human_translation.html`, `v6/figures/html/v6_FigS8_translation_detail.html` | `v6/evaluation/` | Use pathway/target conservation framing | Gene-level transfer is weak/partial |
| 15 | v7.1 Public Boundary | `v7/figures/html/v7_unified_signal_hierarchy.html`, `v7/figures/html/v7_methods_comparison.html` | `docs/CANONICAL_RESULTS_V7_1.md` | Add release-boundary callout box | Keep v8 out of v7.1 benchmark claims |
| 16 | v8 SpaceMed | `v8/figures/Figure1_Species_Transfer.png`, `v8/figures/Figure2_Stressor_Decomposition.png` | `v8/RESULTS_SUMMARY.md`, `v8/bridge/evaluation/`, `v8/decompose/evaluation/` | Possibly use two-panel summary | Hypothesis-generation only |
| 17 | v8 Intervention/Causal Boundaries | `v8/figures/Figure4_Countermeasure_Pareto.png`, `v8/figures/Figure5_Causal_DAG.png` | `v8/intervene/evaluation/`, `v8/causal/evaluation/` | Add red/gray "not clinical recommendation" boundary | No countermeasure recommendation |
| 18 | v9 Platform Architecture | New schematic | `v9/README.md`, `spacebio_bench/`, `v9/task_manifests/` | Draw package/manifest/evaluator/run-manifest flow | Platform scaffold, not frozen benchmark release |
| 19 | Public Bulk Alpha | New status table | `v9/reports/public_bulk_alpha_gap_matrix/`, `v9/reports/public_bulk_alpha_snapshot_decision/`, `v9/datapackage.draft.json` | Pass/blocker matrix; metadata-only alpha box | Not payload-hash-verified |
| 20 | Organoid Track | New compact table/flow | `v9/human_organoid/`, organoid reports | Table: OSD-863/871, 42 samples, DE reference, diagnostics | Draft diagnostic, not leaderboard |
| 21 | Multispecies Track | New compact table/flow | `v9/multispecies/` reports | Table: OSD-207, OSD-37, OSD-120 status | OSD-120 diagnostic-only |
| 22 | Single-Cell Track | New blocker/status visual | `v9/sc_spaceflight/` reports | RRRM asset inventory + metric contract + missing h5ad blocker | Not runnable yet |
| 23 | Claimed / Not Claimed | New boundary matrix | master outline | Convert claim matrix into clean slide | Reviewer-safety slide |
| 24 | Roadmap | New execution roadmap | master outline, `docs/V9_PURPOSE_DRIFT_AUDIT_2026_05_26.md` | Figure QA, v4 reconciliation, payload staging, formal search export | Future work, not current result |

## Short-Talk Slide Asset Plan

| Slide | Keep from long version | Asset choice |
|---:|---|---|
| 1 | SpaceBio-Bench | New title visual |
| 2 | Why mission-held-out? | New LOMO schematic |
| 3 | External landscape | New positioning matrix |
| 4 | Data/task design | New design schematic or `figures/figS1_data_controls.html` |
| 5 | Thymus beats liver | `figures/fig1_tissue_hierarchy.html` |
| 6 | Pathway abstraction | `figures/fig2_pathway_mechanism.html` |
| 7 | Model comparison | `figures/fig3_model_comparison.html` |
| 8 | Benchmark hardening | `v4/figures/html/Fig1_benchmark.html` plus canonical table |
| 9 | v8 incubator | `v8/figures/Figure1_Species_Transfer.png` plus boundary callout |
| 10 | v9 platform | New architecture schematic |
| 11 | Extension tracks | New organoid/multispecies/sc status table |
| 12 | Current status | New roadmap / blockers slide |

## Existing Figure Candidates By Section

### Root / v1-v4 Core

- `figures/fig1_tissue_hierarchy.html`
- `figures/fig2_pathway_mechanism.html`
- `figures/fig3_model_comparison.html`
- `figures/fig4_validation.html`
- `figures/figS1_data_controls.html`
- `figures/figS2_nes_multidb.html`
- `figures/figS3_transfer_detail.html`
- `figures/figS4_temporal_bio.html`
- `figures/figS5_pipeline_comparison.html`

Use:

- Good candidates for slides 5-9 and manuscript Figures 2-4.

Needs:

- Render QA.
- Static export into unified deck style.
- Check whether the displayed values match `docs/CANONICAL_RESULTS_V7_1.md`.

### v2 Temporal / Human / RRRM-1

- `v2/figures/Fig1_temporal.html`
- `v2/figures/F1_temporal_heatmap.html`
- `v2/figures/Fig2_crossspecies.html`
- `v2/figures/E1_crossspecies_scatter.html`
- `v2/figures/E2_duration_conservation.html`
- `v2/figures/E3_cfrna_origin.html`
- `v2/figures/Fig3_pbmc_celltype.html`
- `v2/figures/F1_celltype_pathway_heatmap.html`
- `v2/figures/Fig4_scrna_summary.html`
- `v2/figures/F2A_composition.html`
- `v2/figures/F2B_celltype_nes_heatmap.html`
- `v2/figures/F2C_loao_auroc.html`
- `v2/figures/F2D_crossspecies_scatter.html`

Use:

- Best deck candidates are `Fig1_temporal.html`, `Fig3_pbmc_celltype.html`, and
  `Fig4_scrna_summary.html`.

Needs:

- Select one or two panels only; v2 can easily sprawl.

### v3 Extensions

- `v3/figures/v3_Fig1_multispecies_concordance.html`
- `v3/figures/v3_Fig2_spatial_overview.html`
- `v3/figures/v3_Fig3_rrrm2_scrna.html`
- `v3/figures/v3_Fig4_extensions.html`
- `v3/figures/v3_Fig5_fm_comparison.html`

Use:

- Good "negative results and failure modes" slide.
- Spatial brain negative and RRRM-2 PBMC/NK positive make a balanced pair.

Needs:

- Keep cross-species Drosophila-mouse result carefully worded as directionality
  / nontrivial transfer, not global species impossibility.

### v4-v7 Hardening / Mechanism / Translation

- `v4/figures/html/Fig1_benchmark.html`
- `v4/figures/html/Fig2_ablation.html`
- `v4/figures/html/Fig3_consensus.html`
- `v4/figures/html/Fig4_network.html`
- `v4/figures/html/Fig5_generalize.html`
- `v4/figures/html/Fig6_multiplatform.html`
- `v5/figures/html/Fig7_immune_signaling.html`
- `v5/figures/html/Fig8_metabolism_drugs.html`
- `v5/figures/html/Fig9_consensus_panel.html`
- `v6/figures/html/v6_Fig10_human_translation.html`
- `v7/figures/html/v7_methods_comparison.html`
- `v7/figures/html/v7_unified_signal_hierarchy.html`

Use:

- v4 should support "not cherry-picked."
- v5/v6 should support "prediction to biology and translation."
- v7 should support "public-safe consolidation."

Needs:

- Resolve v4 canonical/raw discrepancy before any "best per tissue" slide.
- Avoid too many detailed mechanism panels in the main deck.

### v8 SpaceMed

- `v8/figures/Figure1_Species_Transfer.png`
- `v8/figures/Figure2_Stressor_Decomposition.png`
- `v8/figures/Figure3_Mars_Extrapolation.png`
- `v8/figures/Figure4_Countermeasure_Pareto.png`
- `v8/figures/Figure5_Causal_DAG.png`

Use:

- These are presentation-ready file types.
- `Figure1` and `Figure2` are safest for main narrative.
- `Figure4` and `Figure5` are useful only with explicit hypothesis-only and no
  recommendation labels.
- `Figure3` is visually useful but carries high overclaim risk.

Needs:

- Decide whether Mars extrapolation appears in main deck or backup only.

### v9 Platform / Extension Tracks

No primary figure directory equivalent was found for v9. Recommended new visual
assets:

1. v9 platform architecture:
   - `task_manifest -> registry -> submission validation -> evaluator -> run_manifest`.
2. public bulk alpha status matrix:
   - passes, blockers, needs-update.
3. extension track status table:
   - public bulk, organoid, multispecies, single-cell.
4. single-cell blocker diagram:
   - asset inventory exists, manifest exists, metric spec exists, canonical
     h5ad payload missing.

Primary data sources:

- `v9/README.md`
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

## Figure QA Checklist

For each selected figure:

1. Open/render the HTML, PNG, or PDF.
2. Confirm it is not blank and text is legible.
3. Confirm numeric labels match the canonical source selected for the slide.
4. Export to static PNG/SVG/PDF in a deck-ready resolution.
5. Record the exact source file and export file.
6. Add a one-line claim-boundary caption.

## Open Decisions

1. Should the long deck be 24 slides, or should it be compressed to 18 slides?
2. Should v8 Mars extrapolation be main-deck or backup?
3. Should v5/v6 get separate slides or one combined "biology/translation" slide?
4. Should v9 organoid and multispecies be separate slides, or one extension
   matrix?
5. Should the deck use existing HTML figure screenshots or regenerated figures
   in one consistent visual system?

## Recommended Next Work Block

Run visual QA and static export for the highest-priority existing figures:

1. `figures/fig1_tissue_hierarchy.html`
2. `figures/fig2_pathway_mechanism.html`
3. `figures/fig3_model_comparison.html`
4. `v4/figures/html/Fig1_benchmark.html`
5. `v2/figures/Fig1_temporal.html`
6. `v3/figures/v3_Fig2_spatial_overview.html`
7. `v8/figures/Figure1_Species_Transfer.png`
8. `v8/figures/Figure2_Stressor_Decomposition.png`

This eight-asset pass is enough to decide whether the deck should reuse
existing figures or regenerate a unified figure set.

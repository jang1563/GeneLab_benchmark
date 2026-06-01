# Visual Premium Rebuild Pass 1

Date: 2026-05-31

Purpose: record the first rebuild pass after the decision that weak existing
figures can be replaced rather than polished.

## Decision Applied

Existing visual assets are now treated as prototypes/references. Main
presentation and manuscript visuals may be rebuilt from scratch when that gives
better clarity, consistency, and claim control.

JK style correction:

- Do not use dashboard/card-box visual language for scientific figures.
- High-profile journal style should prefer thin-rule tables, dot/lollipop
  summaries, restrained schematics, axis-based plots, and quiet whitespace.
- Status/release figures must look like paper figures, not product UI.

## New Script

Created:

- `scripts/build_premium_visuals.py`

Current capability:

- Builds a first premium static core figure from canonical values in
  `evaluation/RESULTS_SUMMARY.md`.
- Exports PNG and PDF.
- Writes CSV and JSON source manifests.
- Uses `Agg` headless Matplotlib backend and local `MPLCONFIGDIR`.
- Builds six first-pass premium static figures:
  - core tissue / pathway result;
  - pathway artifact / rescue result;
  - model-tier comparison with explicit benchmark-surface boundaries;
  - v9 platform-status figure separating public-bulk, organoid, multispecies,
    and single-cell readiness lanes;
  - public-bulk metadata-alpha release-boundary figure;
  - human organoid diagnostic-surface figure.

P0 revision applied:

- Figure 1 was rebuilt as a tissue-transfer hierarchy:
  - `output/premium_figures/premium_fig1_tissue_transfer_hierarchy.png`
  - `output/premium_figures/premium_fig1_tissue_transfer_hierarchy.pdf`
  - `output/premium_figures/manifests/premium_fig1_tissue_transfer_hierarchy_manifest.csv`
  - `output/premium_figures/manifests/premium_fig1_tissue_transfer_hierarchy_manifest.json`
- The legacy Figure 1 path is also regenerated for continuity, but the current
  preferred asset is `premium_fig1_tissue_transfer_hierarchy`.
- Figure 2 Panel C no longer uses a fitted line or an outlier-exclusion claim.
- Figure 6 Panel B now separates prediction and biology-check metrics into
  separate mini-axes.
- Figure-visible jargon was audited and reduced; full decisions are recorded in
  `docs/VISUAL_TEXT_JARGON_AUDIT_2026_05_31.md`.

Figure/table separation applied:

- Full decisions are recorded in
  `docs/VISUAL_FIGURE_TABLE_SEPARATION_2026_05_31.md`.
- Table-like panels were moved to `output/premium_tables/`.
- Figure 1 is now a dot/interval plot only; its main-readout summary is a
  table.
- Figure 4 is now a platform architecture schematic; lane readiness and
  footprint counts are tables.
- Figure 5 is now a sparse release-boundary schematic; release gates, release
  options, and claim language are tables.
- Figure 6 no longer contains the decision table panel; organoid biology-check
  decisions are a table.

L4 finalization QA applied:

- Numeric source audit script added:
  `scripts/audit_premium_visual_sources.py`.
- Manuscript-width preview script added:
  `scripts/render_premium_manuscript_previews.py`.
- Numeric audit result: 117/117 checks passed.
- Manuscript preview and L4-readiness decisions are recorded in
  `docs/VISUAL_L4_FINALIZATION_QA_2026_05_31.md`.
- One rounding fix was applied: Mouse-Geneformer skin delta is now `-0.264`,
  matching `0.557 - 0.821` from `evaluation/RESULTS_SUMMARY.md`.

## Generated Prototype

Generated outputs:

- `output/premium_figures/premium_fig1_core_tissue_pathway.png`
- `output/premium_figures/premium_fig1_core_tissue_pathway.pdf`
- `output/premium_figures/manifests/premium_fig1_core_tissue_pathway_manifest.csv`
- `output/premium_figures/manifests/premium_fig1_core_tissue_pathway_manifest.json`

Figure content:

- Panel A: Category B cross-mission transfer AUROC by tissue with confidence
  intervals.
- Panel B: Category A LOMO detection AUROC by tissue with FDR marker.
- Panel C: gene-vs-pathway detection comparison, highlighting task-specific
  pathway rescue.

Source:

- `evaluation/RESULTS_SUMMARY.md`

Claim boundary:

- Mouse bulk GeneLab/OSDR LOMO benchmark surface.
- Not a universal tissue, model, or foundation-model claim.
- Pathway rescue is selected/task-specific.

## Visual QA

First render:

- Successfully generated 3520 x 1980 PNG and PDF.
- Visual inspection found two issues:
  - title/subtitle overlapped with panel A title;
  - panel C legend overlapped lower-right labels.

Fix applied:

- Lowered plot grid top margin.
- Moved panel C legend above the panel.
- Re-rendered successfully.

Pending:

- Second render visual inspection completed in the next session.
- Title/subtitle and panel C legend overlap were resolved.
- Current status: usable L3 prototype, not final L4.

## Quality Assessment

Current level:

- L2/L3 prototype.

What is already better than legacy HTML:

- Static PNG/PDF outputs.
- No remote D3 dependency.
- Unified layout and palette.
- Source manifest emitted with the figure.
- Claim boundary embedded as footer.

What is still needed for L4 premium:

- Recheck second render visually.
- Verify every number against the final source-of-truth table.
- Consider simplifying panel B or moving it to a supplement if the figure feels
  too dense.
- Decide whether Panel C belongs in this figure or should be separated into a
  dedicated pathway-rescue figure.
- Add formal QA score against the premium rubric.

## Next Recommended Block

1. Decide whether Figure 1 should keep panel C or move pathway rescue entirely
   into Figure 2.
2. Run a numeric audit from final source-of-truth tables for Figures 1-3.
3. Start the v9 architecture/status visual family.

Do not proceed to a full deck until the first six premium figures share a
consistent design system and each has a source manifest.

## Second Generated Prototype

Generated outputs:

- `output/premium_figures/premium_fig2_pathway_artifact_rescue.png`
- `output/premium_figures/premium_fig2_pathway_artifact_rescue.pdf`
- `output/premium_figures/manifests/premium_fig2_pathway_artifact_rescue_manifest.csv`
- `output/premium_figures/manifests/premium_fig2_pathway_artifact_rescue_manifest.json`

Figure content:

- Panel A: gene versus pathway feature performance for condition/confounder
  prediction.
- Panel B: pathway-minus-gene AUROC delta for Category A detection tasks.
- Panel C: pathway NES concordance versus cross-mission transfer AUROC.

Visual QA:

- Static PNG/PDF generation succeeded.
- Layout is substantially stronger than the legacy HTML/QuickLook route.
- Panel A is clear and high-value for explaining artifact suppression.
- Panel C is useful but should be checked for label placement in manuscript
  size.

Current level:

- L3 prototype.

Needed for L4:

- Numeric source verification against a final claim table.
- Maybe shorten the title for slide use.
- Decide whether gastrocnemius outlier label should be visually de-emphasized
  more strongly.

## Third Generated Prototype

Generated outputs:

- `output/premium_figures/premium_fig3_model_tier_comparison.png`
- `output/premium_figures/premium_fig3_model_tier_comparison.pdf`
- `output/premium_figures/manifests/premium_fig3_model_tier_comparison_manifest.csv`
- `output/premium_figures/manifests/premium_fig3_model_tier_comparison_manifest.json`

Figure content:

- Panel A: classical ML, single-cell FM, and zero-shot text LLM mean scores.
- Panel B: scGPT and Mouse-Geneformer tissue-level AUROC deltas versus matched
  classical baselines.
- Panel C: benchmark-surface boundary note separating 6-tissue means,
  zero-shot text tasks, and v3 best-row extension results.

Source:

- `docs/CANONICAL_RESULTS_V7_1.md`
- `evaluation/RESULTS_SUMMARY.md`

Visual QA:

- Static PNG/PDF generation succeeded.
- Output size: 3520 x 1980 PNG.
- Visual inspection found no blocking overlaps.
- Panel A is strong for the talk/manuscript story because the classical
  baseline remains visibly above current FM/LLM rows.
- Panel B is useful because it prevents the figure from becoming only a mean
  leaderboard; it shows local scGPT gains but mostly negative tissue deltas.
- Panel C is essential because the comparison mixes 6-tissue v1 means,
  zero-shot text tasks, and v3 extension best rows.

Current level:

- L3 prototype.

Needed for L4:

- Keep the mixed-surface warning in the caption and/or callout.
- Consider whether the Panel A x-axis should start at 0.0 for a conservative
  manuscript version, while keeping the current 0.4 zoom for slides.
- Make final caption explicit that this is not a universal 8-tissue FM
  leaderboard.
- Run a final numeric audit after freezing the canonical results table.

## Fourth Generated Prototype

Generated outputs:

- `output/premium_figures/premium_fig4_v9_platform_status.png`
- `output/premium_figures/premium_fig4_v9_platform_status.pdf`
- `output/premium_figures/manifests/premium_fig4_v9_platform_status_manifest.csv`
- `output/premium_figures/manifests/premium_fig4_v9_platform_status_manifest.json`

Figure content:

- Panel A: v9 lane readiness matrix for public mouse bulk, human organoid,
  multispecies, and single-cell RRRM.
- Panel B: evidence-footprint tiles showing key counts across the four lanes.
- Panel C: platform architecture from public sources through manifests,
  metrics, and release gates.

Source:

- `v9/README.md`
- `docs/V9_LONG_HORIZON_EXECUTION_PLAN.md`
- `v9/reports/public_bulk_alpha_snapshot_decision/snapshot_decision_summary.csv`
- `v9/reports/public_bulk_alpha_gap_matrix/public_bulk_alpha_gap_summary.csv`
- `v9/human_organoid/*` draft summaries and reports
- `v9/multispecies/*` draft summaries and reports
- `v9/sc_spaceflight/*` asset and obs/var audit summaries

Visual QA:

- Static PNG/PDF generation succeeded.
- Output size: 3520 x 1980 PNG.
- First render had a small Panel C title/body spacing issue in the final
  architecture box.
- Fix applied by lowering body text and setting explicit top alignment.
- Second render has no blocking overlaps.
- After JK style review, removed dashboard card-box treatment:
  - status matrix is now a thin-rule table;
  - footprint panel is now a compact evidence list;
  - architecture panel is now a numbered schematic spine.

Current level:

- L3 prototype.

Needed for L4:

- Decide whether this is a deck-only status slide or also a manuscript overview
  figure.
- If used in manuscript, replace some tiles with a formally sourced table
  panel or move the footprint into supplementary material.
- Keep the claim boundary: v9 platform status, not final benchmark results.

## Fifth Generated Prototype

Generated outputs:

- `output/premium_figures/premium_fig5_public_bulk_alpha_boundary.png`
- `output/premium_figures/premium_fig5_public_bulk_alpha_boundary.pdf`
- `output/premium_figures/manifests/premium_fig5_public_bulk_alpha_boundary_manifest.csv`
- `output/premium_figures/manifests/premium_fig5_public_bulk_alpha_boundary_manifest.json`

Figure content:

- Panel A: metadata-alpha gate ladder showing pass, needs-update, and blocked
  gates.
- Panel B: selected release path: metadata-only alpha selected, payload mirror
  deferred, no-alpha-until-payload-frozen rejected.
- Panel C: allowed versus prohibited language for slides, dataset card, and
  manuscript.

Source:

- `v9/reports/public_bulk_alpha_gap_matrix/public_bulk_alpha_gap_matrix.csv`
- `v9/reports/public_bulk_alpha_snapshot_decision/snapshot_option_matrix.csv`
- `v9/reports/public_bulk_alpha_snapshot_decision/snapshot_claim_boundary.csv`

Visual QA:

- Static PNG/PDF generation succeeded.
- Output size: 3520 x 1980 PNG.
- First render had Panel C text overflow from allowed-language cards into the
  prohibited-language column.
- Fix applied by shortening the allowed-language phrases while preserving the
  underlying claim boundary in the manifest.
- Second render has no blocking overlaps.
- After JK style review, removed card/box treatment:
  - gate ladder is now a thin-rule table;
  - release path is now a vertical decision spine;
  - language boundary is now a two-column rule table.

Current level:

- L3 prototype.

Needed for L4:

- If used in a manuscript, back this with a formal release-boundary table.
- If used in a deck, keep it as a clear governance/status slide rather than a
  scientific result slide.
- Re-render after the dataset card and Data Package alpha-boundary updates are
  complete, because two gates are explicitly marked `needs update`.

## Sixth Generated Prototype

Generated outputs:

- `output/premium_figures/premium_fig6_organoid_diagnostic_surface.png`
- `output/premium_figures/premium_fig6_organoid_diagnostic_surface.pdf`
- `output/premium_figures/manifests/premium_fig6_organoid_diagnostic_surface_manifest.csv`
- `output/premium_figures/manifests/premium_fig6_organoid_diagnostic_surface_manifest.json`

Figure content:

- Panel A: human organoid provenance footprint for OSD-863/OSD-871.
- Panel B: primary classification metrics versus diagnostic-only biological
  metrics.
- Panel C: enrichment among top-ranked model genes.
- Panel D: release-boundary decision table for default, secondary, and negative
  diagnostic artifacts.

Source:

- `v9/human_organoid/source_inventory.draft.csv`
- `v9/human_organoid/sample_factors.draft.csv`
- `v9/human_organoid/de_references/human_organoid_de_reference_manifest.draft.json`
- `v9/human_organoid/reports/nearest_centroid/human_organoid_baseline_summary.csv`
- `v9/human_organoid/reports/source_transfer_signature/metrics.json`
- `v9/human_organoid/reports/logistic_feature_effect/metrics.json`
- `v9/human_organoid/reports/ORGANOID_FEATURE_EFFECT_NULL_CALIBRATION_REVIEW.md`
- `v9/human_organoid/reports/ORGANOID_DIAGNOSTIC_CONSOLIDATION_AND_RELEASE_BOUNDARY.md`

Visual QA:

- Static PNG/PDF generation succeeded.
- Output size: 3520 x 1980 PNG.
- First render had three visual issues:
  - Panel A tile numbers were too close to the status strip;
  - Panel B warning text overlapped the first bar;
  - Panel C and D titles visually competed.
- Fixes applied:
  - lowered tile value text;
  - removed the redundant Panel B warning text;
  - shortened Panel C and D titles and increased inter-panel spacing.
- Final render has no blocking overlaps.
- After JK style review, removed card/tile treatment:
  - Panel A is now an evidence spine instead of tiles;
  - Panel D is now a thin-rule decision table instead of row boxes.

Current level:

- L3 prototype.

Needed for L4:

- Keep the caption explicit that organoid metrics are diagnostic-only.
- Do not compare response-signature rank correlation, feature-effect rank
  correlation, and classification AUROC as one leaderboard.
- For manuscript use, consider converting Panel D into a concise supplementary
  table and keeping Panels A-C as the main figure.

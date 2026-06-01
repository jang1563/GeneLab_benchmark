# Project Results Location Inventory

Date: 2026-05-31

Purpose: slide-deck preparation source map. This document records where the
project results live, what each result surface contains, and which files should
be read first before turning v1-v9 into a deck.

Companion audit: `docs/PROJECT_RESULTS_DEEP_READ_AUDIT_2026_05_31.md` records
what was manually read, what was programmatically parsed, and what remains
unreviewed.

## Executive Read Order

Read in this order when building the slide narrative:

1. `docs/PROJECT_SLIDE_CONTENT_INVENTORY_V1_TO_V9_2026_05_31.md`
   - Current high-level v1-v9 story, headline results, and slide angles.
2. `docs/CANONICAL_RESULTS_V7_1.md`
   - Public-safe v1-v7 benchmark scope, canonical counts, and release-boundary
     language.
3. `evaluation/RESULTS_SUMMARY.md`
   - Root/v1-v4 style benchmark results, foundation-model/LLM comparison, and
     detailed core benchmark numbers.
4. `v2/evaluation/V2_RESULTS_SUMMARY.md`
   - Temporal, preservation/recovery, age, human cfRNA/PBMC, and RRRM-1
     single-cell expansion.
5. `v3/README.md`
   - Multi-species, spatial, RRRM-2, radiation, and UCE/scFoundation result
     summary.
6. `v8/RESULTS_SUMMARY.md`
   - SpaceMed v8 bridge, stressor, intervention, Mars-regime, and exploratory
     hypothesis-generation result surface.
7. `v9/README.md` and `v9/OPERATING_BACKLOG.md`
   - SpaceBio-Bench v9 platform state, public bulk alpha, organoid,
     multispecies, and single-cell scaffold.
8. Figure directories:
   - `figures/`
   - `v2/figures/`
   - `v3/figures/`
   - `v4/figures/html/`
   - `v5/figures/html/`
   - `v6/figures/html/`
   - `v7/figures/html/`
   - `v8/figures/`

## Result Surface Counts

These counts are file counts from the current working tree, excluding
`submissions/`.

| Area | Files | What it is |
|---|---:|---|
| `evaluation/` | 121 | Root/v1-v4 style benchmark JSON summaries, FM/LLM outputs, controls |
| `figures/` | 9 | Root manuscript-style HTML figures and supplements |
| `processed/` | 418 | Root processed matrices, pathway scores, QC, intermediate task data |
| `docs/` | 35 | Canonical summaries, release notes, v8/v9 research/planning/review docs |
| `v2/evaluation/` | 46 | Temporal, cross-species, PBMC, RRRM-1, LLM parse-aware outputs |
| `v2/figures/` | 13 | v2 temporal/cross-species/single-cell HTML figures |
| `v2/processed/` | 112 | v2 processed temporal, cfRNA, PBMC, RRRM-1, LLM prompt artifacts |
| `v3/evaluation/` | 21 | Multi-species, spatial, RRRM-2, radiation, FM JSON outputs |
| `v3/figures/` | 5 | v3 interactive HTML figures |
| `v4/evaluation/` | 302 | 8-tissue x method x feature benchmark grid and controls |
| `v4/figures/` | 22 | v4 figure scripts and rendered HTML figures |
| `v4/wgcna_outputs/` | 42 | Tissue-level WGCNA module outputs |
| `v5/evaluation/` | 25 | Immune, TF, metabolism, drug target, biomarker, cross-organ biology |
| `v5/figures/` | 10 | v5 biological-mechanism figure scripts and HTML outputs |
| `v6/evaluation/` | 6 | Human translation, pathway/gene/TF/drug target validation summaries |
| `v6/figures/` | 4 | v6 translation figure script and HTML outputs |
| `v7/evaluation/` | 24 | GNN and scPRINT-2 comparator outputs |
| `v7/figures/` | 2 | v7 method comparison and unified hierarchy figures |
| `v8/bridge/evaluation/` | 10 | Mouse-human bridge and species-transfer result files |
| `v8/decompose/evaluation/` | 18 | Stressor decomposition, Mars extrapolation, saturation sensitivity |
| `v8/intervene/evaluation/` | 34 | LINCS, offline reversal, CRISPR orthogonal, Pareto, safety triage |
| `v8/causal/evaluation/` | 4 | DAG and ICP/stressor pathway score outputs |
| `v8/multiomics/evaluation/` | 2 | RNA-to-clinical / propagation-style pilot outputs |
| `v8/provenance/` | 21 | v8 run manifests, freeze/audit manifests, schemas |
| `v8/release/` | 1 | v8 beta artifact manifest |
| `v8/figures/` | 11 | v8 PDF/PNG main figures plus generator |
| `v9/reports/` | 106 | Public bulk baselines, alpha gap/snapshot/recenter reports |
| `v9/human_organoid/` | 159 | Human organoid inventories, matrices, DE references, diagnostics |
| `v9/multispecies/` | 626 | Multi-species inventories, task manifests, OSD-120 diagnostic branch |
| `v9/sc_spaceflight/` | 35 | Single-cell RRRM inventory, manifest draft, metric spec, payload audit |
| `v9/task_manifests/` | 8 | v9 public bulk LOMO task manifest exports |

## Root / v1 Result Surface

Primary files:

- `README.md`
  - Project-level overview, public benchmark framing, v7.1/v8 boundary.
- `docs/V1_PAPER_CONTENT.md`
  - Original paper-ready v1 methods/results narrative.
- `evaluation/RESULTS_SUMMARY.md`
  - Main root benchmark result source, including A/B/C/D/J controls, v4
    headline table, foundation-model and LLM comparison.
- `docs/CANONICAL_RESULTS_V7_1.md`
  - Safest public-facing result source for v1-v7 scope and claim language.
- `tasks/`
  - Task definitions and fold data for legacy LOMO benchmark tasks.
- `processed/`
  - Processed task matrices, pathway scores, QC reports, and intermediate
    analysis outputs.
- `evaluation/*.json`
  - Machine-readable benchmark outputs.
- `figures/*.html`
  - Root HTML figures:
    - `fig1_tissue_hierarchy.html`
    - `fig2_pathway_mechanism.html`
    - `fig3_model_comparison.html`
    - `fig4_validation.html`
    - `figS1`-`figS5`

What is here:

- Original mouse bulk RNA-seq LOMO benchmark.
- Core categories: spaceflight detection, cross-mission transfer,
  cross-tissue transfer, confounder prediction, gene-vs-pathway comparison,
  negative controls, and held-out validation.
- Classical ML, gene-expression foundation model, and text LLM result surfaces.

Slide-use claims:

- The benchmark is mission-held-out, not random-split-only.
- Thymus and other immune-rich tissues generalize better than liver.
- Pathway abstraction can reduce mission/batch artifacts.
- Current FM/LLM rows do not beat tuned classical baselines on this small-n
  domain-shift benchmark.

Claim caution:

- v1 FM rows are mostly 6-tissue surfaces; do not present them as a uniform
  8-tissue v4 leaderboard.

## v2 Result Surface

Primary files:

- `v2/README.md`
- `v2/evaluation/V2_RESULTS_SUMMARY.md`
- `v2/evaluation/*.json`
- `v2/processed/`
- `v2/figures/*.html`

Main result families:

- T1: ISS-Terminal versus Live Animal Return timing and preservation artifact.
- T2: LAR recovery and overshoot signatures.
- T3: age by spaceflight interaction in RR-8 liver.
- J1: GLDS pipeline-version comparison.
- LLM zero-shot parse-aware benchmark.
- E1-E3: mouse-human pathway/cfRNA conservation and duration effects.
- F1: I4 PBMC cell-type pathway analysis.
- F2: RRRM-1 single-cell hardening and benchmark task definitions.

Key slide-ready results:

- T1 shows preservation-method artifact can dominate timing effects.
- T2 shows return samples project toward baseline, with overshoot in many
  pathways.
- T3 supports an age-amplified spaceflight signature: OLD liver AUROC 0.945 vs
  YNG 0.679.
- Mouse liver pathway NES partially conserves with human JAXA cfRNA
  (Spearman r=0.352).
- RRRM-1 single-cell outputs cover 38,081 cells after hardening.

Claim caution:

- T1 should be framed as a confound/artifact lesson, not a clean biological
  timing effect.
- LLM parse failures are part of end-to-end benchmark behavior and should be
  reported separately from parsed-only behavior.

## v3 Result Surface

Primary files:

- `v3/README.md`
- `v3/evaluation/*.json`
- `v3/figures/*.html`
- `v3/scripts/`

Main result families:

- E4/E5: mouse, Drosophila, and cross-species NES/phylogenetic analysis.
- F3: spatial Visium brain benchmark.
- F5: RRRM-2 single-cell benchmark.
- A7/A7b/A8: extended lung, colon, and skin bulk tasks.
- R1/R2/R3: radiation analog analyses.
- B_ext: extended cross-tissue transfer matrix.
- FM: UCE and scFoundation comparisons.

Key slide-ready results:

- Spatial brain is a strong negative result: section AUROC 0.139 and
  animal-level AUROC 0.444, with PC1 dominated by slide batch.
- RRRM-2 PBMC has strong cell-type signal, especially NK AUROC 0.845.
- Bone marrow is near chance across cell types.
- UCE/scFoundation remain below the classical PCA-LR benchmark surface.

Claim caution:

- Spatial Visium negative result should be treated as evidence of tissue/task
  specificity, not a failed pipeline.
- Cross-species Drosophila-mouse directionality is not a direct mammalian
  translational claim.

## v4 Result Surface

Primary files:

- `v4/evaluation/`
- `v4/figures/html/`
- `v4/wgcna_outputs/`
- `evaluation/RESULTS_SUMMARY.md`
- `docs/CANONICAL_RESULTS_V7_1.md`

Main result families:

- 8 tissues x 8 classifiers x 4 feature types = 256 evaluation grid.
- Ablation and negative controls.
- Consensus summaries.
- WGCNA module outputs by tissue.

Key slide-ready results:

- Best rows:
  - Thymus 0.948, PCA-LR + KEGG.
  - Colon 0.921, PCA-LR + KEGG.
  - Lung 0.901, PCA-LR + gene.
  - Kidney 0.829, ElasticNet-LR + Hallmark.
  - Eye 0.823, PCA-LR + Hallmark.
  - Skin 0.819, PCA-LR + gene.
  - Gastrocnemius 0.776, PCA-LR + gene.
  - Liver 0.670, PCA-LR + gene.
- PCA-LR is best overall by 8-tissue gene-level mean AUROC: 0.776.
- ElasticNet-LR is second: 0.762.
- 40/256 configurations are significant at p<0.05.
- 6/8 tissues have at least one significant configuration.

Claim caution:

- Not every row in the best-method table has a significant best-row
  permutation p-value. Use the canonical note from
  `docs/CANONICAL_RESULTS_V7_1.md`.

## v5 Result Surface

Primary files:

- `v5/evaluation/`
- `v5/figures/html/`
- `v5/scripts/`

Main result families:

- Immune deconvolution by tissue.
- TF activity by tissue.
- Metabolic flux.
- Drug target mapping.
- Consensus biomarker panel.
- Cross-organ signaling.

Key slide-ready results:

- Skin has the strongest immune-deconvolution signal among v5 tissues.
- TF activity is especially rich in thymus, skin, kidney, and liver.
- Drug-target mapping links mouse candidate genes to human targetability:
  1,919 mouse target genes, 834 mapped to human, 271 with DGIdb interactions,
  and 3,154 total DGIdb interactions.
- Consensus panel contains 20 genes, with validation AUROC varying by tissue
  from stronger gastrocnemius/liver/eye/kidney to weaker colon/lung.

Claim caution:

- v5 is interpretation and hypothesis prioritization. It should not be framed
  as drug efficacy or clinical validation.

## v6 Result Surface

Primary files:

- `v6/evaluation/V6_A_gene_conservation.json`
- `v6/evaluation/V6_B_pathway_conservation.json`
- `v6/evaluation/V6_C_cross_species_transfer.json`
- `v6/evaluation/V6_D_biomarker_validation.json`
- `v6/evaluation/V6_E_tf_conservation.json`
- `v6/evaluation/V6_F_drug_target_validation.json`
- `v6/figures/html/`

Main result families:

- Mouse-human gene conservation.
- Pathway conservation.
- Cross-species transfer.
- Biomarker validation in human cfRNA.
- TF conservation.
- Drug target validation.

Key slide-ready results:

- Pathway-level conservation is partial: mean rho 0.285 across five human
  pathway surfaces.
- Human cfRNA detects 15/20 consensus panel genes, but no DE FDR<0.05 panel
  genes in the tested setting.
- TF conservation is weak overall: mean rho about 0.030.
- Drug target validation yields 3 Tier A candidates and 7 Tier B promising
  candidates among the evaluated target set.

Claim caution:

- v6 supports partial translational conservation, not direct replacement of
  human validation with mouse data.

## v7 / v7.1 Result Surface

Primary files:

- `docs/CANONICAL_RESULTS_V7_1.md`
- `docs/V7_V8_CLOSURE_STATUS_2026_05_10.md`
- `v7/evaluation/`
- `v7/figures/html/`

Main result families:

- GNN variants using no edges, random edges, or WGCNA edges.
- scPRINT-2 comparator rows.
- Unified signal hierarchy.
- Public-release consistency and closure status.

Key slide-ready results:

- WGCNA/GNN variants do not overturn the classical-baseline conclusion.
- scPRINT-2 rows are below PCA-LR comparisons across the available tissues.
- v7.1 is the public-facing benchmark cleanup and release-boundary layer.
- Clean-checkout HPC validation passed at the documented v7/v8 closure commit:
  47/47 tests passed.

Claim caution:

- Keep v7.1 benchmark paper claims separate from v8 intervention, Mars, or
  countermeasure language.

## v8 Result Surface

Primary files:

- `v8/README.md`
- `v8/RESULTS_SUMMARY.md`
- `v8/RESULTS_SUMMARY_TABLE.csv`
- `v8/bridge/evaluation/`
- `v8/decompose/evaluation/`
- `v8/intervene/evaluation/`
- `v8/causal/evaluation/`
- `v8/multiomics/evaluation/`
- `v8/provenance/`
- `v8/release/v8_beta_artifact_manifest.json`
- `v8/figures/*.png` and `v8/figures/*.pdf`

Main result families:

- BRIDGE: rodent-human pathway/NES bridge into human I4/Twins style outputs.
- DECOMPOSE: stressor decomposition and Mars-regime extrapolation.
- INTERVENE: LINCS perturbation reversal, CRISPR orthogonal support, safety
  triage, Pareto prioritization.
- CAUSAL: DAG/ICP result surface.
- MULTIOMICS: early RNA-to-clinical or propagation-style result surface.
- PROVENANCE: run manifests, input-freeze/audit files, beta artifact manifest.

Key slide-ready results:

- Mouse NES improves human I4 post-flight pathway prediction: RF supervised
  AUROC 0.888 vs baseline 0.712, improvement +0.175 [0.134, 0.219].
- Top ICP-stability stressor is T, with ICP=0.5401.
- Interaction terms dominate 44-61% variance in top-responsive genes.
- Low-LET gamma and high-LET HZE show opposite correlation directions in one
  stressor result surface.
- Top multi-tissue signature-reversal hypothesis is CGP-60474.
- Linear Mars extrapolation breaks above 5x dose amplification and should be
  interpreted as exploratory regime flagging.

Claim caution:

- v8 is an incubator/beta-facing translational extension. Intervention and
  Mars outputs are hypothesis-generation and validation-prioritization
  artifacts, not clinical, operational, or countermeasure recommendations.
- v8 beta is not complete unless clean-checkout recomputation, external source
  freeze, manifest validation, and artifact split are finished.

## v9 Result Surface

Primary files:

- `v9/README.md`
- `v9/OPERATING_BACKLOG.md`
- `docs/V9_LONG_HORIZON_EXECUTION_PLAN.md`
- `docs/V9_LONG_RUN_OPERATING_PROTOCOL.md`
- `docs/V9_PURPOSE_DRIFT_AUDIT_2026_05_26.md`
- `docs/V9_ORGANOID_AND_SPECIES_EXTENSION_REVIEW_2026_05_21.md`
- `docs/V9_PUBLIC_BULK_ALPHA_METADATA_SNAPSHOT_DECISION.md`
- `docs/V9_SC_METRIC_SPEC.md`
- `docs/V9_SC_PAYLOAD_STAGING_PLAN.md`
- `docs/V9_SC_OBS_VAR_AUDIT.md`

Core platform files:

- `spacebio_bench/`
  - v9 package spine: manifests, profiles, metrics, evaluator, reports,
    public-bulk loader, organoid/multispecies/single-cell helpers.
- `v9/task_manifests/*.json`
  - Eight generated public bulk LOMO manifests.
- `v9/task_manifest_index.csv` and `.json`
  - Compact manifest registry.
- `v9/task_data_index.csv` and `.json`
  - Fold-level data path and count summaries.
- `v9/source_inventory.csv` and `.json`
  - Source-level inventory for 22 public bulk source rows.
- `v9/source_checksum_audit.csv` and `.json`
  - OSDR API and checksum-manifest evidence.
- `v9/datapackage.draft.json`
  - Draft Frictionless Data Package descriptor.

Public bulk result files:

- `v9/reports/nearest_centroid/`
  - Baseline predictions, metrics, and run manifests.
- `v9/reports/sklearn_baselines/`
  - PCA-LR and L2 logistic baseline predictions, metrics, and run manifests.
- `v9/reports/bulk_lomo_baseline_summary.csv` and `.json`
  - 24 baseline rows across nearest-centroid, PCA-LR, and L2 logistic.
- `v9/reports/public_bulk_alpha_gap_matrix/`
  - Public bulk freeze-gap review artifacts.
- `v9/reports/public_bulk_alpha_snapshot_decision/`
  - Metadata-only alpha snapshot decision artifacts.
- `v9/reports/recenter_decision/`
  - Decision to refocus from OSD-120 branch toward public bulk alpha and
    single-cell flagship.

Human organoid result files:

- `v9/human_organoid/source_inventory.draft.csv`
- `v9/human_organoid/task_manifests/*.json`
- `v9/human_organoid/sample_factors.draft.csv`
- `v9/human_organoid/geo_sample_metadata.draft.csv`
- `v9/human_organoid/expression_matrix_audit.draft.csv`
- `v9/human_organoid/signature_reference_audit.draft.csv`
- `v9/human_organoid/de_references/human_organoid_de_reference.draft.csv.gz`
- `v9/human_organoid/reports/nearest_centroid/`
- `v9/human_organoid/reports/sensitivity/`
- `v9/human_organoid/reports/source_transfer_signature/`
- `v9/human_organoid/reports/microglia_source_transfer_signature/`
- `v9/human_organoid/reports/logistic_feature_effect/`
- `v9/human_organoid/reports/pca_lr_feature_effect/`
- `v9/human_organoid/reports/ORGANOID_DIAGNOSTIC_CONSOLIDATION_AND_RELEASE_BOUNDARY.md`

Human organoid slide-ready status:

- Public human neural organoid extension is present as a draft/non-primary v9
  track.
- OSD-863 and OSD-871 are included.
- GEO GSE259421 metadata recovers 42/42 public sample-level Subject/iPSC-line
  identifiers.
- Derived DE reference contains 242,708 gene/contrast rows and 2,368
  FDR<=0.05 rows.
- Diagnostic response-signature and feature-effect metrics exist, but they
  remain secondary or diagnostic surfaces.

Multispecies result files:

- `v9/multispecies/source_inventory.draft.csv`
- `v9/multispecies/sample_factors.draft.csv`
- `v9/multispecies/expression_matrix_audit.draft.csv`
- `v9/multispecies/source_checksum_audit.draft.csv`
- `v9/multispecies/task_manifests/`
- `v9/multispecies/reports/nearest_centroid/`
- `v9/multispecies/reports/sensitivity/`
- `v9/multispecies/interaction_task_manifests/`
- `v9/multispecies/reports/interaction_nearest_centroid/`
- `v9/multispecies/reports/interaction_sensitivity/`
- `v9/multispecies/reports/interaction_logistic_l2/`
- `v9/multispecies/reports/interaction_logistic_sparse_l1/`
- `v9/multispecies/reports/interaction_public_alpha_card/`
- `v9/multispecies/reports/interaction_diagnostic_metadata_release_note/`

Multispecies slide-ready status:

- OSD-37 Arabidopsis and OSD-207 Drosophila are species-native draft task
  candidates.
- OSD-120 is handled separately as a light/genotype interaction diagnostic
  case, not merged into the simple plant task.
- OSD-120 sparse L1 logistic `tvg2000_log1p_zscore_l1_c1` is the leading
  draft diagnostic candidate in the OSD-120 branch.
- The OSD-120 branch is now closed as diagnostic metadata unless owner metadata
  appears.

Single-cell result files:

- `v9/sc_spaceflight/asset_inventory_summary.csv`
- `v9/sc_spaceflight/asset_inventory.csv`
- `v9/sc_spaceflight/local_payload_scan.csv`
- `v9/sc_spaceflight/task_manifests/draft_rrrm1_blood_single_cell_spaceflight.json`
- `v9/sc_spaceflight/anndata_manifest_draft_summary.csv`
- `v9/sc_spaceflight/metric_spec_summary.csv`
- `v9/sc_spaceflight/metric_spec_metrics.csv`
- `v9/sc_spaceflight/metric_spec_required_inputs.csv`
- `v9/sc_spaceflight/metric_spec_skip_policy.csv`
- `v9/sc_spaceflight/payload_staging_plan_summary.csv`
- `v9/sc_spaceflight/payload_staging_candidates.csv`
- `v9/sc_spaceflight/obs_var_audit_requirements.csv`
- `v9/sc_spaceflight/obs_var_audit_summary.csv`
- `v9/sc_spaceflight/obs_var_audit_results.csv`
- `v9/sc_spaceflight/payload_manifest.csv`

Single-cell slide-ready status:

- V9-SC-001 records 54 legacy RRRM/single-cell asset paths.
- V9-SC-002 selects OSD-918 RRRM-1 blood as the first non-runnable AnnData
  task-manifest contract.
- Draft contract records 8 samples, 4 Flight, 4 Ground, 4,395 QC cells, and
  19,064 genes.
- V9-SC-003 fixes six `genelab_sc` metric formulas and skip policies.
- V9-SC-004 defines the canonical future h5ad target and obs/var audit gate.
- V9-SC-005 currently emits a skip-only audit because the canonical h5ad payload
  is missing.

v9 claim caution:

- v9 is platform/alpha scaffold, not a frozen benchmark-result release.
- Public bulk alpha is metadata-only because local payload hash verification is
  not complete.
- Human organoid results are draft/diagnostic and should not be promoted as
  primary benchmark conclusions yet.
- Multispecies OSD-120 is diagnostic metadata, not a leaderboard or release
  archive.
- Single-cell is manifest/metric/audit scaffolding until the canonical payload
  is staged or regenerated.

## Slide Deck Source Buckets

Use these buckets when converting the inventory into slides.

| Deck bucket | Best source files | What to extract |
|---|---|---|
| Problem and benchmark design | `README.md`, `docs/V1_PAPER_CONTENT.md`, `docs/CANONICAL_RESULTS_V7_1.md` | Mission-held-out benchmark concept, task families, public scope |
| Core result | `evaluation/RESULTS_SUMMARY.md`, `v4/evaluation/`, root `figures/` | Tissue hierarchy, pathway rescue, FM/LLM underperformance |
| Robustness and extensions | `v2/evaluation/V2_RESULTS_SUMMARY.md`, `v3/README.md` | Temporal artifacts, recovery, age, cfRNA, PBMC, RRRM, spatial negative |
| Method hardening | `v4/evaluation/`, `v7/evaluation/`, `v7/figures/html/` | 256-grid result, GNN/scPRINT-2 comparisons |
| Biological interpretation | `v5/evaluation/`, `v6/evaluation/` | Immune, TF, metabolism, drug targets, human translation |
| Translational incubator | `v8/RESULTS_SUMMARY.md`, `v8/figures/`, `v8/provenance/` | Bridge, stressor, intervention hypothesis, Mars regime caution |
| Platform future | `v9/README.md`, `v9/OPERATING_BACKLOG.md`, `spacebio_bench/` | Manifests, evaluators, public bulk alpha, organoid/multispecies/single-cell |

## Current Boundary Notes For Slides

- Use v1-v7 as the main benchmark story.
- Use v8 as a clearly labeled translational incubator.
- Use v9 as the platformization and future-release story.
- Keep operational/clinical countermeasure language out of the main deck unless
  explicitly framed as hypothesis generation.
- Do not claim a frozen v9 release; present it as metadata alpha plus active
  extension scaffolds.
- Do not touch or summarize `submissions/` unless it becomes an explicit task.

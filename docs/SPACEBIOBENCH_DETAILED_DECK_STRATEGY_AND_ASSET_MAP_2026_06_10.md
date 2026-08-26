# SpaceBio-Bench Detailed Deck Strategy And Asset Map

Date: 2026-06-10

Purpose: planning artifact for a premium, consulting-grade detailed deck. This
deck is separate from the concise conference deck. The concise deck argues the
thesis. The detailed deck shows the reasoning chain, evidence base, robustness
checks, claim boundaries, and platform-readiness logic.

## North Star

The detailed deck should not be a longer chronology of v1-v9. It should be a
guided technical walkthrough:

> Can spaceflight transcriptomic signatures generalize to unseen missions, and
> can biological AI models learn that signal without overclaiming beyond the
> evidence?

Each slide should answer one of three questions:

1. What was tested?
2. What happened?
3. What can and cannot be claimed?

## Target Format

- Primary use: technical seminar, lab meeting, reviewer walkthrough, or long
  conference workshop.
- Target length: 60-65 slides.
- Talk length: 40-50 minutes plus Q&A.
- Audience: mixed computational biology, space biology, and ML/LLM familiarity.
- Design target: premium scientific consulting deck; high visual hierarchy,
  strong proof objects, sparse interpretation overlays, explicit claim
  boundaries.

## Narrative Architecture

### Act 1: Problem To Benchmark Method

Goal: make the audience understand why the benchmark is a mission-held-out
evaluation problem, not a generic classification exercise.

Recommended slides:

1. Title / thesis.
2. External gap: OSDR has many studies, but generalization is not guaranteed.
3. Project map: v1-v9 as claim layers, not equal-strength discoveries.
4. What counts as a benchmark task.
5. Data-to-task contract: source record, study, mission, sample, matrix, task.
6. What the model sees: sample-by-feature matrices.
7. Mission-held-out split: the hidden unit is an entire mission.
8. Leakage guard: feature choices stop before the hidden mission.
9. Metric primer: AUROC, uncertainty, and chance.
10. Model-family primer: classical ML, expression foundation models, text LLMs.
11. Claim-boundary primer: benchmark evidence is not mechanism or clinical
    evidence.

Primary existing visual assets:

- `output/premium_methods_scenes/methods_data_to_evaluation_overview.png`
- `output/premium_methods_dark_variants_slides_4_5_v0_1/slide04_evaluation_layer_dark_rendered_preview.png`
- `output/premium_methods_dark_variants_slides_4_5_v0_1/slide05_mission_heldout_dark_rendered_preview.png`
- `output/premium_bridge_rebuild_scenes/b2_study_to_task_premium/rendered_preview.png`
- `output/premium_bridge_rebuild_scenes/b3_mission_held_out_premium/rendered_preview.png`
- `output/premium_bridge_rebuild_scenes/b4_train_only_guard_premium/rendered_preview.png`
- `output/premium_feature_layer_bridge/feature_layer_bridge_scene.png`

### Act 2: Core Bulk Benchmark Evidence

Goal: establish the main scientific result, then explain why it happens.

Recommended slides:

12. Worked example: reading one tissue score.
13. Tissue hierarchy: some tissues generalize better than others.
14. Why liver is not the winner despite more data.
15. Transfer matrix view: mission-pair evidence behind the tissue ranking.
16. Pathway rescue: pathways suppress selected unwanted labels.
17. Gene vs pathway interpretation: complementary, not universally superior.
18. New slide required: pathway conservation predicts transfer.
19. New slide required: NES conservation scatter or rank bridge.
20. Practical implication: fGSEA can screen transfer feasibility before model
    training.
21. Held-out validation: RR-23 thymus and RR-7 skin.
22. Negative controls: chance behavior where signal should be absent.

Primary existing visual assets:

- `output/premium_figures/premium_fig1_tissue_transfer_hierarchy.png`
- `output/premium_slide_scenes/fig1_tissue_transfer_layered_scene.png`
- `output/premium_figures/premium_fig2_pathway_artifact_rescue.png`
- `output/premium_slide_scenes/fig2_pathway_layered_scene.png`

New premium assets required:

- `NES conservation vs transfer AUROC` scatter/rank bridge.
- `Held-out validation` slide for RR-23 thymus and RR-7 skin.
- `Transfer matrix explainer` if existing figures are too manuscript-like.

Progress:

- Built `NES conservation vs transfer AUROC` premium proof asset:
  `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-detailed-deck/assets/nes_conservation/nes_conservation_predicts_transfer_premium.png`
- Grayscale QA render:
  `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-detailed-deck/assets/nes_conservation/nes_conservation_predicts_transfer_premium_grayscale.png`
- Manifest:
  `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-detailed-deck/assets/nes_conservation/nes_conservation_predicts_transfer_premium_manifest.json`
- Key readout: all six tissues shown; primary 5-tissue Spearman `r_s=0.90`
  excludes gastrocnemius because only two fGSEA missions have NES data; strict
  4-tissue check gives `r_s=1.00`.

### Act 3: Models, Robustness, And External Validation

Goal: keep the model comparison from becoming a leaderboard slide. It should
read as a controlled benchmark result under fixed task conditions.

Recommended slides:

23. Why baseline strength matters.
24. Classical ML result surface.
25. Foundation model setup: pretrained expression models adapted to small-n bulk.
26. Text LLM setup: prompt-only diagnostic, not direct matrix learning.
27. Model comparison: scale alone does not transfer.
28. Why bulk RNA-seq is a hard domain for single-cell-pretrained models.
29. v4 hardening: 8 tissues x 8 classifiers x 4 feature types.
30. v7/scPRINT/GNN checks: newer model ideas preserve the benchmark lesson.
31. New slide required: DGE pipeline robustness.
32. New slide required: external biological validation.
33. Robustness summary: leakage, split, metric, negative controls, external
    concordance.

Primary existing visual assets:

- `output/premium_figures/premium_fig3_model_tier_comparison.png`
- `output/premium_slide_scenes/fig3_model_tier_layered_scene.png`
- `v4/figures/html/Fig1_benchmark.html`
- `v7/figures/html/v7_methods_comparison.html`

New premium assets required:

- `DGE robustness` proof slide: DESeq2 / edgeR / limma-voom, Log2FC rank rho
  0.926, DEG Jaccard 0.600.
- `External validation` proof slide: Cell 2020 pathway concordance 71.7% and
  Gene SHAP overlap 47x above chance.

### Act 4: Biology And Translation Layers

Goal: show that the benchmark is biologically interpretable, while keeping
mechanism and translation claims bounded.

Recommended slides:

34. Temporal and RRRM guardrails: preservation, recovery, underpowered pilots.
35. Spatial and negative biological boundaries.
36. v5 systems-biology interpretation overview.
37. Immune deconvolution and TF activity.
38. Metabolic / drug target / biomarker layers as triage, not treatment.
39. v6 mouse-to-human translation overview.
40. Gene-level weak transfer vs pathway partial conservation.
41. Drug target validation tiers as prioritization evidence.
42. Biology summary: prediction -> interpretation -> bounded hypothesis.

Primary existing visual/source assets:

- `v2/figures/Fig1_temporal.html`
- `v2/figures/Fig3_pbmc_celltype.html`
- `v3/figures/v3_Fig2_spatial_overview.html`
- `v5/figures/html/Fig7_immune_signaling.html`
- `v5/figures/html/Fig8_metabolism_drugs.html`
- `v5/figures/html/Fig9_consensus_panel.html`
- `v6/figures/html/v6_Fig10_human_translation.html`

Likely remake requirement:

- Most v2-v6 HTML figures are better used as source material than dropped
  directly into a premium deck. Rebuild the most important ones as native
  visual summaries with one claim per slide.

### Act 5: v8 Translational Incubator Boundary

Goal: include the translational work without implying interventions or crew
health recommendations.

Recommended slides:

43. Why translation appears after benchmark evidence.
44. v8 pillar map: BRIDGE, DECOMPOSE, INTERVENE, CAUSAL.
45. BRIDGE: mouse NES improves human pathway prediction.
46. DECOMPOSE: stressor interaction terms and radiation-context fragility.
47. INTERVENE: perturbation hits are hypothesis-only.
48. Mars / countermeasure boundary: no operational recommendation.
49. v8 summary: useful incubator, not validated intervention evidence.

Primary existing visual assets:

- `v8/figures/Figure1_Species_Transfer.png`
- `v8/figures/Figure2_Stressor_Decomposition.png`
- `v8/figures/Figure4_Countermeasure_Pareto.png`
- `v8/figures/Figure5_Causal_DAG.png`

Risk note:

- Countermeasure and Mars visuals should be backup or heavily bounded. In the
  main detailed deck, they should read as prioritization/hypothesis surfaces.

### Act 6: v9 Platform, Provenance, And Extension Readiness

Goal: show platform maturity, provenance discipline, and release blockers.

Recommended slides:

50. Platform turn: from analysis branch to reproducible benchmark system.
51. Object chain: manifest, evaluator, run record.
52. Task registry and metric profiles.
53. Public bulk metadata alpha.
54. Payload/readiness ladder: parsed source vs mirrored payload vs hash-verified
    release.
55. Source checksum audit and Data Package boundary.
56. Organoid extension: useful biology-check dataset, non-primary.
57. OSD-120 multispecies: same-study diagnostic, not mission-held-out.
58. Single-cell RRRM: metric spec exists, processed payload missing.
59. Roadmap: freeze payloads, QA package, release and paper alignment.
60. Final claim boundary slide.

Primary existing visual assets:

- `output/premium_figures/premium_fig4_v9_platform_architecture.png`
- `output/premium_figures/premium_fig4_v9_platform_status.png`
- `output/premium_figures/premium_fig5_public_bulk_alpha_boundary.png`
- `output/premium_figures/premium_fig5_public_bulk_release_boundary_schematic.png`
- `output/premium_v9_document_scenes/v9_platform_provenance_document_scene.png`
- `output/premium_v9_document_scenes/v9_public_bulk_boundary_document_scene.png`
- `output/premium_figures/premium_fig6_organoid_diagnostic_surface.png`
- `output/premium_slide_scenes/fig6_organoid_layered_scene.png`
- `output/biovis_organoid_audience_matrix_proof_v0_2/panels/01_dark_organoid_clean_source_to_matrix.png`
- `output/biovis_source_object_proof_crops_v0_1/production_mocks/02_dark_osd120_split_proof_replacement.png`

## Visual Asset Register

Generated asset register:

- CSV: `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-detailed-deck/assets/spacebiobench_detailed_deck_asset_register.csv`
- JSON: `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-detailed-deck/assets/spacebiobench_detailed_deck_asset_register.json`
- Contact sheet: `outputs/019e8bd1-3b42-7821-9e1c-beebfa8e2ece/presentations/spacebiobench-detailed-deck/assets/spacebiobench_detailed_deck_asset_contact_sheet.png`

Register summary:

- 45 total asset records.
- 30 image assets inspected in contact-sheet form.
- Existing premium assets are strongest for methods, core tissue/pathway/model
  results, v9 platform, organoid, and source-proof visuals.
- New asset work is most important for NES conservation, held-out validation,
  DGE pipeline robustness, and external biological validation.

## Consulting-Grade Design Rules

1. One slide, one claim, one proof object.
2. Keep dense manuscript figures as source material unless they are already
   readable at thumbnail scale.
3. Use concept interludes before hard results: task, split, metric, feature,
   model family, claim boundary.
4. Prefer native diagrams for explanations and high-resolution raster charts
   for proof objects.
5. Carry a consistent data strip: source, split, metric, claim status.
6. Use explicit guardrail badges where claim drift is likely.
7. Separate core benchmark evidence from extension and release-readiness lanes.
8. Treat v8 intervention/Mars surfaces as hypothesis-only.
9. Treat v9 as metadata-alpha/platform readiness unless payload verification
   passes.
10. QA every final slide at full size, contact-sheet size, and grayscale.

## Immediate Build Sequence

Recommended next work block:

1. Build the detailed deck spine as a slide table with slide number, claim,
   proof object, asset path, speaker purpose, and claim boundary.
2. Produce the four missing premium proof assets:
   - NES conservation predicts transfer.
   - RR-23/RR-7 held-out validation.
   - DGE pipeline robustness.
   - Cell 2020 / SHAP external validation.
3. Rebuild Act 1 methods primer first, because all later result slides depend
   on the audience understanding task, split, metric, and model families.
4. Only then assemble Act 2 core benchmark evidence.

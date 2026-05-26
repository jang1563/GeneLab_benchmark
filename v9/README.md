# SpaceBio-Bench v9 Planning Artifacts

Status: planning/alpha scaffold. This directory does not claim v9 benchmark
results yet.

## Current contents

- `task_manifests/*.json`: v9 task manifests exported from the existing
  `tasks/*/task_info.json` bulk LOMO tasks.
- `task_manifest_index.csv` and `task_manifest_index.json`: compact registry
  summaries generated from the task manifests.
- `task_data_index.csv` and `task_data_index.json`: fold-level data path/count
  summaries for the public bulk LOMO scaffold.
- `reports/nearest_centroid/`: first runnable v9-alpha baseline report across
  all eight generated bulk LOMO manifests.
- `reports/sklearn_baselines/`: PCA-LR and L2 logistic-regression reports
  across all eight generated bulk LOMO manifests.
- `reports/bulk_lomo_baseline_summary.csv` and `.json`: normalized
  cross-baseline summary for nearest-centroid plus sklearn baselines.
- `source_inventory.csv` and `source_inventory.json`: source-level inventory
  deduplicated from all generated bulk LOMO task manifests.
- `source_checksum_audit.csv` and `source_checksum_audit.json`: OSDR API
  file-list and checksum-manifest evidence for the 22 public bulk source rows.
- `datapackage.draft.json`: draft Frictionless Data Package descriptor for
  v9 public bulk metadata, provenance evidence, and baseline outputs.
- `human_organoid/source_inventory.draft.csv` and `.json`: draft source
  inventory for public human neural organoid RNA-seq extension candidates
  OSD-863 and OSD-871.
- `human_organoid/task_manifests/*.json`: draft human organoid task manifest.
- `human_organoid/task_manifest_index.draft.csv` and `.json`: draft index for
  human organoid extension task manifests.
- `human_organoid/source_checksum_audit.draft.csv` and `.json`: OSDR API
  file-list and checksum-manifest evidence for OSD-863 and OSD-871.
- `human_organoid/sample_table_audit.draft.csv` and `.json`: OSDR sample-table
  discovery and row-count evidence for OSD-863 and OSD-871.
- `human_organoid/sample_factors.draft.csv` and `.json`: parsed sample-level
  disease, Ground/LEO, microglia, organoid-type, GEO Subject, and iPSC-line
  factors for the 42 public GSE259421 samples.
- `human_organoid/geo_sample_metadata.draft.csv` and `.json`: GEO GSE259421
  series-matrix-derived sample titles, Subject identifiers, iPSC-line ids,
  BioSample/SRA accessions, cell type, and treatment metadata.
- `human_organoid/expression_matrix_audit.draft.csv` and `.json`: OSDR
  normalized-count matrix discovery, download, SHA-256, row-count, and
  sample-column alignment evidence for OSD-863 and OSD-871.
- `human_organoid/signature_reference_audit.draft.csv` and `.json`: OSDR
  DE/signature reference audit showing public differential-expression tables,
  contrast-definition files, direct/reversed spaceflight contrast counts, and
  the current metric policy.
- `human_organoid/de_references/human_organoid_de_reference.draft.csv.gz` and
  `_manifest.draft.json`: derived all-gene DE reference table for the eight
  direct OSD-863/OSD-871 Ground Control versus Space Flight contrasts, with
  log2FC normalized to `LEO_or_ISS - Ground`.
- `human_organoid/matrices/`: downloaded human organoid normalized-count
  matrices referenced by the draft expression-matrix audit.
- `human_organoid/reports/nearest_centroid/`: first draft human organoid pilot
  baseline outputs, including predictions, metrics, run manifest, and summary.
- `human_organoid/reports/sensitivity/`: draft sensitivity grid for the human
  organoid pilot baseline across transform, scaling, and variable-gene settings.
- `human_organoid/reports/ORGANOID_BASELINE_ROBUSTNESS_REVIEW.md`: conservative
  robustness and confounding review for the draft organoid pilot baseline.
- `human_organoid/sample_scale_diagnostics.draft.csv` and `.json`: per-sample
  expression scale diagnostics over common human gene features.
- `human_organoid/group_scale_diagnostics.draft.csv` and `.json`: grouped
  sample-total and sparsity summaries by source, label, organoid type,
  microglia condition, disease context, donor/Subject id, and iPSC-line id.
- `human_organoid/donor_metadata_audit.draft.csv` and `.json`: donor/iPSC-line
  metadata availability audit showing OSDR SampleTables lack donor fields but
  GEO series metadata recovers Subject/iPSC-line ids for all 42 samples.
- `human_organoid/reports/ORGANOID_DONOR_LIBRARY_SIZE_AUDIT.md`: interpretation
  note for donor metadata gaps and library-size diagnostics.
- `human_organoid/reports/donor_diagnostics/`: diagnostic-only GEO
  donor/Subject holdout baseline outputs.
- `human_organoid/reports/ORGANOID_DONOR_AWARE_SPLIT_DECISION.md`: decision
  note keeping donor holdouts diagnostic-only rather than promoting them to the
  default benchmark split.
- `human_organoid/reports/ORGANOID_SIGNATURE_METRIC_REFERENCE_DECISION.md`:
  decision note keeping DE/signature metrics reference-backed but non-primary
  after the V9-ORG-014 input contract and V9-ORG-015 diagnostic scorer.
- `human_organoid/reports/ORGANOID_DE_REFERENCE_CONTRACT.md`: V9-ORG-014
  contract for the derived DE reference, response-signature artifact schema,
  and skip-aware metric policy.
- `human_organoid/reports/ORGANOID_RESPONSE_SIGNATURE_SCORER.md`: V9-ORG-015
  diagnostic scorer contract for supplied `response_signature.csv` artifacts.
- `human_organoid/reports/response_signature_smoke/`: V9-ORG-016 smoke-test
  report exercising the response-signature scorer against the real derived DE
  reference with a clearly marked non-model fixture.
- `human_organoid/reports/ORGANOID_RESPONSE_SIGNATURE_ADAPTER_DESIGN.md`:
  V9-ORG-017 design for the first non-leaky model-produced
  `response_signature.csv` adapter.
- `human_organoid/reports/source_transfer_signature/`: V9-ORG-018
  source-transfer empirical response-signature baseline report, using
  train-fold-only signatures and post hoc DE-reference scoring.
- `human_organoid/reports/ORGANOID_SOURCE_TRANSFER_SIGNATURE_REVIEW.md`:
  V9-ORG-019 diagnostic review of the source-transfer signature baseline,
  including per-contrast behavior, sign-imbalance checks, leakage boundaries,
  and the next adapter decision.
- `human_organoid/reports/ORGANOID_PER_CONDITION_SIGNATURE_ADAPTER_DESIGN.md`:
  V9-ORG-020 design for a microglia-matched source-transfer response-signature
  adapter.
- `human_organoid/reports/microglia_source_transfer_signature/`: V9-ORG-021
  microglia-matched source-transfer empirical response-signature report and
  comparison against the global source-transfer diagnostic.
- `human_organoid/reports/ORGANOID_SOURCE_TRANSFER_ADAPTER_COMPARISON_REVIEW.md`:
  V9-ORG-022 comparison review for global versus microglia-matched
  source-transfer diagnostics.
- `human_organoid/reports/ORGANOID_SHARED_CONTROL_SIGNATURE_ADAPTER_DESIGN.md`:
  V9-ORG-023 design for a partial shared-control disease+microglia
  source-transfer diagnostic.
- `human_organoid/reports/shared_control_source_transfer_signature/`:
  V9-ORG-024 partial shared-control disease+microglia source-transfer
  response-signature report.
- `human_organoid/reports/ORGANOID_CLASSIFIER_FEATURE_EFFECT_CONTRACT.md`:
  V9-ORG-025 contract separating classifier-derived feature effects from
  log2FC response signatures.
- `human_organoid/reports/logistic_feature_effect/`: V9-ORG-026 L2 logistic
  gene-space feature-effect pilot report.
- `human_organoid/reports/ORGANOID_FEATURE_EFFECT_DIAGNOSTIC_REVIEW.md`:
  V9-ORG-027 review of L2 logistic feature-effect diagnostics versus
  response-signature diagnostics.
- `human_organoid/reports/ORGANOID_FEATURE_EFFECT_NULL_CALIBRATION_REVIEW.md`:
  V9-ORG-028 review of calibrated top-k null diagnostics for L2 logistic
  feature effects.
- `human_organoid/reports/ORGANOID_PCA_LR_FEATURE_EFFECT_DESIGN.md`:
  V9-ORG-029 design for PCA-LR reconstructed gene-space feature effects.
- `human_organoid/reports/pca_lr_feature_effect/`: V9-ORG-030 PCA-LR
  reconstructed gene-space feature-effect pilot report.
- `human_organoid/reports/ORGANOID_PCA_LR_FEATURE_EFFECT_REVIEW.md`:
  V9-ORG-030 comparison review for PCA-LR reconstructed feature effects versus
  the L2 logistic feature-effect report.
- `human_organoid/reports/ORGANOID_DIAGNOSTIC_CONSOLIDATION_AND_RELEASE_BOUNDARY.md`:
  V9-ORG-031 consolidation note separating draft-alpha organoid diagnostics
  from exploratory and negative-comparison reports.
- `multispecies/source_inventory.draft.csv` and `.json`: draft source
  inventory for non-mouse pilot sources OSD-207, OSD-37, and OSD-120.
- `multispecies/sample_factors.draft.csv` and `.json`: parsed sample-level
  Ground/LEO label, genotype/ecotype, condition-stratum, and light-treatment
  factors for the 124 local OSD-207, OSD-37, and OSD-120 samples.
- `multispecies/expression_matrix_audit.draft.csv` and `.json`: local
  normalized-count matrix row-count and sample-column alignment audit for
  OSD-207, OSD-37, and OSD-120.
- `multispecies/source_checksum_audit.draft.csv` and `.json`: OSDR API
  file-list and checksum-manifest evidence for OSD-207, OSD-37, and OSD-120.
- `multispecies/task_manifests/*.json`: draft species-native multispecies task
  manifests for OSD-37 Arabidopsis and OSD-207 Drosophila.
- `multispecies/task_manifest_index.draft.csv` and `.json`: draft index for
  species-native multispecies task manifests.
- `multispecies/interaction_task_manifests/*.json`: draft OSD-120 Arabidopsis
  root genotype/ecotype by light-treatment interaction task manifest.
- `multispecies/interaction_task_manifest_index.draft.csv` and `.json`: draft
  index for interaction-task manifests.
- `multispecies/reports/MULTISPECIES_CANDIDATE_SOURCE_DEEP_AUDIT.md`:
  V9-MULTI-003 deep audit of OSD-207, OSD-37, and OSD-120 task readiness.
- `multispecies/reports/MULTISPECIES_CHECKSUM_AND_LOCAL_PAYLOAD_AUDIT.md`:
  V9-MULTI-005 checksum and local payload spot-check review for the draft
  multispecies inputs.
- `multispecies/reports/MULTISPECIES_SPECIES_NATIVE_TASK_MANIFEST_DESIGN.md`:
  V9-MULTI-006 design note for the OSD-37 and OSD-207 species-native task
  manifests and OSD-120 interaction-task deferral.
- `multispecies/reports/nearest_centroid/`: V9-MULTI-007 draft-only
  nearest-centroid baseline feasibility outputs for the OSD-37 and OSD-207
  species-native task manifests.
- `multispecies/reports/MULTISPECIES_BASELINE_FEASIBILITY_REVIEW.md`:
  V9-MULTI-007 review of multispecies species-native baseline feasibility and
  condition-stratum robustness.
- `multispecies/reports/sensitivity/`: V9-MULTI-008 preprocessing sensitivity
  grid for the OSD-37 and OSD-207 species-native nearest-centroid baselines.
- `multispecies/reports/MULTISPECIES_BASELINE_SENSITIVITY_REVIEW.md`:
  V9-MULTI-008 review selecting the conservative default baseline setting and
  recording OSD-207's recurring weak `w1118_KCNQ370` stratum.
- `multispecies/reports/OSD120_INTERACTION_TASK_DESIGN.md`: V9-MULTI-009
  design note for treating OSD-120 as a separate Arabidopsis root
  genotype/ecotype by light-treatment interaction task.
- `multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json`:
  V9-MULTI-010 task manifest carrying genotype/ecotype, light-treatment, and
  condition-stratum holdout metadata for OSD-120.
- `multispecies/reports/interaction_nearest_centroid/`: V9-MULTI-011
  draft-only OSD-120 interaction nearest-centroid baseline outputs, separated
  by genotype/ecotype primary, light-treatment secondary, and
  condition-stratum diagnostic fold families.
- `multispecies/reports/OSD120_INTERACTION_BASELINE_FEASIBILITY_REVIEW.md`:
  V9-MULTI-011 review of OSD-120 interaction baseline feasibility and fold
  heterogeneity.
- `multispecies/reports/interaction_sensitivity/`: V9-MULTI-012 preprocessing
  sensitivity grid for the OSD-120 interaction nearest-centroid baseline across
  20 variants and three fold families, plus V9-MULTI-013
  `fold_detail_summary.csv` and `.json` aggregation tables for fold-level
  fragile-stratum inspection.
- `multispecies/reports/OSD120_INTERACTION_SENSITIVITY_REVIEW.md`:
  V9-MULTI-012/013 review confirming repeated OSD-120 fragile strata, keeping
  the conservative default baseline setting, and documenting the fold-detail
  summary table.
- `multispecies/reports/OSD120_INTERACTION_LOGISTIC_BASELINE_DESIGN.md`:
  V9-MULTI-014 design note for the first OSD-120 interaction L2 logistic
  diagnostic baseline.
- `multispecies/reports/interaction_logistic_l2/`: V9-MULTI-014 draft-only
  OSD-120 interaction L2 logistic baseline outputs, separated by the same three
  fold families, plus V9-MULTI-015 logistic `fold_detail_summary.csv` and
  side-by-side `fold_detail_comparison_vs_nearest_centroid.csv`.
- `multispecies/reports/OSD120_INTERACTION_LOGISTIC_BASELINE_REVIEW.md`:
  V9-MULTI-014/015 comparison review for L2 logistic versus nearest-centroid
  OSD-120 interaction diagnostics.
- `multispecies/reports/interaction_logistic_l2_sensitivity/`: V9-MULTI-016
  compact top-variable-gene and L2 `C` sensitivity grid for the OSD-120
  interaction logistic baseline, including fold-detail and nearest-centroid
  comparison tables.
- `multispecies/reports/OSD120_INTERACTION_LOGISTIC_SENSITIVITY_REVIEW.md`:
  V9-MULTI-016 review showing that the `Col.0.PhyD|Dark.Treatment` failure is
  feature-count sensitive rather than regularization-sensitive.
- `multispecies/reports/interaction_logistic_feature_audit/`: V9-MULTI-017
  selected-feature and fold-level coefficient audit for the top-500 versus
  top-2,000 OSD-120 logistic tradeoff.
- `multispecies/reports/OSD120_INTERACTION_LOGISTIC_FEATURE_SET_AUDIT_REVIEW.md`:
  V9-MULTI-017 review showing that the top-2,000 extra feature set strongly
  changes coefficient rankings across the fragile OSD-120 folds.
- `multispecies/reports/OSD120_INTERACTION_SPARSE_BRANCH_DESIGN.md`:
  V9-MULTI-018 design note for the next draft-only sparse L1 logistic branch
  and its fragile-fold gates.
- `multispecies/reports/interaction_logistic_sparse_l1/`: V9-MULTI-019
  draft-only sparse L1 logistic pilot outputs, including fold-detail,
  nearest-centroid comparison, and sparse coefficient audit tables.
- `multispecies/reports/OSD120_INTERACTION_SPARSE_L1_REVIEW.md`:
  V9-MULTI-019 review selecting `tvg2000_log1p_zscore_l1_c1` as the current
  best transparent OSD-120 diagnostic candidate.
- `multispecies/reports/interaction_logistic_sparse_l1_stability/`:
  V9-MULTI-020 deterministic train-fold subsampling stability audit for sparse
  L1 `C=0.3` and `C=1.0`.
- `multispecies/reports/OSD120_INTERACTION_SPARSE_L1_STABILITY_REVIEW.md`:
  V9-MULTI-020 review keeping sparse L1 `C=1.0` as the performance-leading
  transparent diagnostic and `C=0.3` as the compact stability comparator.
- `multispecies/reports/interaction_baseline_ladder/`: V9-MULTI-021
  consolidated OSD-120 interaction ladder tables for nearest centroid, dense
  L2, top-500 L2, sparse L1 `C=0.3`, and sparse L1 `C=1.0`.
- `multispecies/reports/OSD120_INTERACTION_BASELINE_LADDER_REVIEW.md`:
  V9-MULTI-021 decision note advancing sparse L1 `C=1.0` as the primary draft
  transparent diagnostic candidate.
- `multispecies/reports/interaction_diagnostic_candidate_package/`:
  V9-MULTI-022 figure-ready OSD-120 sparse L1 `C=1.0` diagnostic candidate
  package with summary, focus evidence, stable-feature evidence, and claim map.
- `multispecies/reports/OSD120_INTERACTION_DIAGNOSTIC_CANDIDATE_PACKAGE_REVIEW.md`:
  V9-MULTI-022 review binding the candidate package to local evidence and
  external OSD-120 light/genotype context while preserving draft-only claim
  boundaries.
- `multispecies/reports/interaction_figure_table_package/`: V9-MULTI-023
  human-facing OSD-120 sparse L1 `C=1.0` figure/table draft package with main
  focus table, stable-feature appendix, caption, and claim-boundary box.
- `multispecies/reports/OSD120_INTERACTION_FIGURE_TABLE_DRAFT_REVIEW.md`:
  V9-MULTI-023 review selecting the main focus table as the primary draft table
  surface and the stable-feature table as a conservative appendix.
- `multispecies/reports/interaction_paired_comparator_table/`: V9-MULTI-024
  optional paired comparator table for sparse L1 `C=1.0` versus compact sparse
  L1 `C=0.3`.
- `multispecies/reports/OSD120_INTERACTION_PAIRED_COMPARATOR_REVIEW.md`:
  V9-MULTI-024 decision note keeping `C=0.3` as appendix/supplement comparator
  rather than a second primary figure panel.
- `multispecies/reports/interaction_diagnostic_artifact_manifest/`:
  V9-MULTI-025 machine-readable OSD-120 diagnostic artifact manifest and
  claim-to-artifact map.
- `multispecies/reports/OSD120_INTERACTION_DIAGNOSTIC_ARTIFACT_MANIFEST_REVIEW.md`:
  V9-MULTI-025 review confirming the OSD-120 diagnostic evidence set is
  internally traceable while remaining draft-only.
- `multispecies/reports/interaction_release_readiness_gap_audit/`:
  V9-MULTI-026 OSD-120 release-readiness summary, gap table, and external
  release-reference table.
- `multispecies/reports/OSD120_INTERACTION_RELEASE_READINESS_GAP_AUDIT.md`:
  V9-MULTI-026 review concluding that the OSD-120 diagnostic package is
  internally traceable but not public-alpha ready until payload freeze,
  source-release target, and OSD-120-specific card/citation blockers are fixed.
- `multispecies/reports/interaction_payload_freeze_manifest/`:
  V9-MULTI-027 OSD-120 diagnostic-required payload freeze summary and
  533-row processed checksum manifest classification.
- `multispecies/reports/OSD120_INTERACTION_PAYLOAD_FREEZE_MANIFEST_REVIEW.md`:
  V9-MULTI-027 review confirming that the two diagnostic-required OSDR
  processed payloads match local MD5 hashes while the broader OSDR processed
  payload set remains outside the current freeze scope.
- `multispecies/reports/interaction_public_alpha_card/`:
  V9-MULTI-028 OSD-120 diagnostic public-alpha card draft and card summary.
- `multispecies/reports/OSD120_INTERACTION_PUBLIC_ALPHA_CARD_REVIEW.md`:
  V9-MULTI-028 review keeping the card source-specific, diagnostic-only,
  payload-boundary-aware, and explicitly not a frozen benchmark release.
- `multispecies/reports/interaction_rebuild_gate/`:
  V9-MULTI-029 OSD-120 packaging rebuild preflight summary, step manifest, and
  environment lock.
- `multispecies/reports/OSD120_INTERACTION_REBUILD_GATE_REVIEW.md`:
  V9-MULTI-029 review confirming that the packaging layer from V9-MULTI-021
  through V9-MULTI-028 is script-backed, present, and content-hashed without
  rerunning sparse L1 model grids.
- `multispecies/reports/interaction_public_metadata_package/`:
  V9-MULTI-030 OSD-120 public metadata summary, source-release target decision,
  metadata field table, external-reference table, and JSON skeleton.
- `multispecies/reports/OSD120_INTERACTION_PUBLIC_METADATA_PACKAGE_REVIEW.md`:
  V9-MULTI-030 review separating the public-now diagnostic metadata draft from
  blocked source-alpha, full-mirror, and frozen-benchmark release targets.
- `multispecies/reports/interaction_ro_crate_citation_scaffold/`:
  V9-MULTI-031 draft RO-Crate metadata, Data Package descriptor, validation
  table, citation-freeze checklist, and export summary.
- `multispecies/reports/OSD120_INTERACTION_RO_CRATE_CITATION_SCAFFOLD_REVIEW.md`:
  V9-MULTI-031 review confirming that the export scaffold is inspectable but
  archive-blocked by identifier, creator/contributor, and license placeholders.
- `multispecies/reports/interaction_archive_decision_gate/`: V9-MULTI-032
  archive identifier option matrix, license/rights decision table,
  creator/contributor decision table, external-reference table, and summary.
- `multispecies/reports/OSD120_INTERACTION_ARCHIVE_IDENTIFIER_LICENSE_DECISION_REVIEW.md`:
  V9-MULTI-032 review selecting no local archive identifier for the current
  diagnostic draft while blocking release DOI, creator, version, and license
  fields until the release owner supplies them.
- `multispecies/reports/interaction_citation_metadata_fill/`: V9-MULTI-033
  release-owner citation metadata intake template, fill-status table, summary,
  and descriptor patch preview.
- `multispecies/reports/OSD120_INTERACTION_CITATION_METADATA_FILL_REVIEW.md`:
  V9-MULTI-033 review confirming that no owner-supplied metadata was available,
  no RO-Crate/Data Package descriptor was mutated, and release blockers remain
  explicit.
- `multispecies/reports/interaction_archive_release_deferral_guard/`:
  V9-MULTI-034 archive-release deferral summary, owner metadata application
  guard table, deferral action table, and descriptor mutation guard JSON.
- `multispecies/reports/OSD120_INTERACTION_ARCHIVE_RELEASE_DEFERRAL_GUARD_REVIEW.md`:
  V9-MULTI-034 review deferring archive release and explicitly blocking
  descriptor mutation while owner metadata, archive route, version,
  creator/publisher, license, and exact OSDR citation fields are missing.
- `multispecies/reports/interaction_diagnostic_metadata_release_note/`:
  V9-MULTI-035 diagnostic metadata release-note summary, included sections,
  public claim boundary table, owner metadata retry checklist, and Markdown
  note.
- `multispecies/reports/OSD120_INTERACTION_DIAGNOSTIC_METADATA_RELEASE_NOTE_REVIEW.md`:
  V9-MULTI-035 closeout review confirming that the current OSD-120 public
  surface is diagnostic metadata only, not a DOI-backed archive release or
  leaderboard promotion, and that the next block should recenter v9.
- `reports/recenter_decision/`: V9-REFOCUS-001 public bulk alpha versus
  single-cell flagship candidate matrix, next-block action table, summary, and
  single-cell asset path scan.
- `docs/V9_REFOCUS_001_POST_OSD120_RECENTER_DECISION.md`: V9-REFOCUS-001
  decision note selecting `V9-BULK-ALPHA-001: public bulk alpha freeze-gap
  matrix` as the next active implementation block.
- `reports/public_bulk_alpha_gap_matrix/`: V9-BULK-ALPHA-001 public bulk
  alpha readiness summary, pass/blocker gap matrix, payload hash boundary
  table, claim-boundary table, and package update plan.
- `docs/V9_PUBLIC_BULK_ALPHA_FREEZE_GAP_MATRIX_REVIEW.md`:
  V9-BULK-ALPHA-001 review showing that public bulk remains the active lane but
  frozen release language is blocked by 0/22 local payload-hash-verified
  sources and an unresolved metadata-only alpha snapshot decision.
- `OPERATING_BACKLOG.md`: living long-run backlog for continuing v9 work across
  sessions.

## Planning docs

- `docs/V9_LONG_HORIZON_EXECUTION_PLAN.md`: long-term execution plan,
  workstreams, phases, decision gates, and risks.
- `docs/V9_PUBLIC_BULK_PACKAGE_DESIGN.md`: public bulk package boundary and
  draft Data Package design.
- `docs/v9_hf_dataset_card.md`: draft Hugging Face-style dataset card for the
  public bulk scaffold.
- `docs/V9_DESIGN_OPTIONS.md`: v9 option matrix and initial 30/60/90-day
  design path.
- `docs/V9_EXTERNAL_DEEP_RESEARCH_2026_05_21.md`: external ecosystem research.
- `docs/V9_SOURCE_AND_COMPETITOR_MATRIX.md`: source and competitor matrix.
- `docs/V9_ORGANOID_AND_SPECIES_EXTENSION_REVIEW_2026_05_21.md`: evidence
  review for human organoid and non-mouse species extension candidates.
- `docs/V9_MULTISPECIES_FEATURE_STRATEGY.md`: feature-namespace strategy for
  species-native and cross-species bridge tasks.
- `docs/V9_PURPOSE_DRIFT_AUDIT_2026_05_26.md`: strategic alignment audit after
  the OSD-120 metadata/deferral chain, concluding that v9 remains on mission.
  V9-MULTI-035 now closes the current OSD-120 metadata branch unless owner
  release metadata appears, and the active next path is a recenter decision.

## Regeneration

```bash
python scripts/export_v9_task_manifests.py
python scripts/build_v9_task_index.py
python scripts/build_v9_task_data_index.py
python scripts/run_v9_nearest_centroid.py \
  --output-dir v9/reports/nearest_centroid
python scripts/run_v9_sklearn_baselines.py \
  --output-dir v9/reports/sklearn_baselines
python scripts/build_v9_baseline_summary.py
python scripts/build_v9_source_inventory.py
python scripts/build_v9_extension_source_inventories.py
python scripts/audit_v9_multispecies_inputs.py
python scripts/audit_v9_source_checksums.py \
  --source-inventory v9/multispecies/source_inventory.draft.csv \
  --csv v9/multispecies/source_checksum_audit.draft.csv \
  --json v9/multispecies/source_checksum_audit.draft.json
python scripts/build_v9_multispecies_task_manifests.py
python scripts/build_v9_osd120_interaction_task_manifest.py
python scripts/run_v9_multispecies_baseline.py
python scripts/run_v9_osd120_interaction_baseline.py
python scripts/run_v9_osd120_interaction_sensitivity.py
python scripts/build_v9_osd120_interaction_fold_details.py
python scripts/run_v9_osd120_interaction_logistic.py
python scripts/build_v9_osd120_logistic_fold_comparison.py
python scripts/run_v9_osd120_interaction_logistic_sensitivity.py
python scripts/audit_v9_osd120_logistic_feature_sets.py
python scripts/run_v9_osd120_interaction_sparse_l1.py
python scripts/audit_v9_osd120_sparse_l1_stability.py
python scripts/build_v9_osd120_baseline_ladder.py
python scripts/build_v9_osd120_diagnostic_candidate_package.py
python scripts/build_v9_osd120_figure_table_package.py
python scripts/build_v9_osd120_paired_comparator_table.py
python scripts/build_v9_osd120_diagnostic_artifact_manifest.py
python scripts/build_v9_osd120_release_readiness_gap_audit.py
python scripts/build_v9_osd120_payload_freeze_manifest.py
python scripts/build_v9_osd120_public_alpha_card.py
python scripts/rebuild_v9_osd120_diagnostic_package.py
python scripts/build_v9_osd120_public_metadata_package.py
python scripts/build_v9_osd120_ro_crate_citation_scaffold.py
python scripts/build_v9_osd120_archive_decision_gate.py
python scripts/build_v9_osd120_citation_metadata_fill.py
python scripts/build_v9_osd120_archive_release_deferral_guard.py
python scripts/build_v9_osd120_diagnostic_metadata_release_note.py
python scripts/build_v9_recenter_decision.py
python scripts/build_v9_public_bulk_alpha_gap_matrix.py
python scripts/run_v9_multispecies_sensitivity.py
python scripts/audit_v9_source_checksums.py \
  --source-inventory v9/human_organoid/source_inventory.draft.csv \
  --csv v9/human_organoid/source_checksum_audit.draft.csv \
  --json v9/human_organoid/source_checksum_audit.draft.json
python scripts/audit_v9_sample_tables.py \
  --source-inventory v9/human_organoid/source_inventory.draft.csv \
  --csv v9/human_organoid/sample_table_audit.draft.csv \
  --json v9/human_organoid/sample_table_audit.draft.json
python scripts/build_v9_human_organoid_sample_factors.py \
  --source-inventory v9/human_organoid/source_inventory.draft.csv \
  --csv v9/human_organoid/sample_factors.draft.csv \
  --json v9/human_organoid/sample_factors.draft.json
python scripts/build_v9_human_organoid_geo_metadata.py \
  --series-matrix-gz path/to/GSE259421_series_matrix.txt.gz
python scripts/audit_v9_expression_matrices.py \
  --source-inventory v9/human_organoid/source_inventory.draft.csv \
  --sample-factors v9/human_organoid/sample_factors.draft.csv \
  --matrix-dir v9/human_organoid/matrices \
  --csv v9/human_organoid/expression_matrix_audit.draft.csv \
  --json v9/human_organoid/expression_matrix_audit.draft.json
python scripts/audit_v9_human_organoid_signature_references.py
python scripts/build_v9_human_organoid_de_reference.py
python scripts/build_v9_human_organoid_task_manifest.py
python scripts/run_v9_human_organoid_baseline.py
python scripts/run_v9_human_organoid_sensitivity.py
python scripts/audit_v9_human_organoid_diagnostics.py
python scripts/run_v9_human_organoid_donor_diagnostics.py
python scripts/run_v9_human_organoid_response_signature_smoke.py
python scripts/run_v9_human_organoid_source_transfer_signature.py
python scripts/run_v9_human_organoid_microglia_source_transfer_signature.py
python scripts/run_v9_human_organoid_shared_control_source_transfer_signature.py
python scripts/run_v9_human_organoid_logistic_feature_effect.py
python scripts/run_v9_human_organoid_pca_lr_feature_effect.py
python scripts/audit_v9_source_checksums.py
python scripts/build_v9_datapackage_draft.py
python scripts/evaluate_v9_submission.py \
  --task-manifest v9/task_manifests/A2_gastrocnemius_bulk_lomo.json \
  --submission path/to/predictions.csv \
  --output-dir path/to/report
python scripts/evaluate_v9_submission.py \
  --task-manifest v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json \
  --submission path/to/predictions.csv \
  --response-signature path/to/response_signature.csv \
  --output-dir path/to/report
```

The exporter reads legacy task metadata from `tasks/`, attaches public OSDR
source records from `scripts/fetch_osdr.py`, and validates each generated JSON
against the minimal `spacebio_bench` manifest contract.

The index builder reads the generated manifests and writes compact CSV/JSON
summaries for quick inspection and downstream loaders.

The data-index builder validates legacy fold files for each bulk LOMO manifest
and writes fold-level path/count summaries.

The evaluator validates a prediction CSV, computes available `genelab_minimal`
metrics, and writes `metrics.json` plus `run_manifest.json`. For organoid
tasks, it can also score diagnostic DE/signature metrics when
`--response-signature path/to/response_signature.csv` is supplied.

The nearest-centroid runner reads the legacy processed fold matrices through
the v9 loader, writes one `predictions.csv`, `metrics.json`, and
`run_manifest.json` per task, then writes a cross-task summary at
`v9/reports/nearest_centroid/bulk_lomo_summary.csv`.

The sklearn runner currently evaluates PCA-LR and L2 logistic regression using
train-fold-fitted preprocessing only. The cross-baseline summary builder
normalizes nearest-centroid and sklearn outputs into
`v9/reports/bulk_lomo_baseline_summary.csv`.

The source-inventory builder deduplicates `source_records` across all task
manifests and writes one row per OSDR source accession. The current inventory
covers 22 public bulk source rows and now includes organism, material,
model-system, assay-modality, and feature-namespace fields. It keeps
`checksum_status` at
`legacy_task_source_unfrozen` because the generated manifests have not yet been
promoted to a frozen v9 release state.

The extension-inventory builder writes draft inventories for public human
organoid and multispecies candidates. These files are planning artifacts, not
frozen benchmark tasks. They keep organoids, plants, and fly data outside the
mouse bulk LOMO leaderboard until split rules, feature namespaces, checksum
audits, and baselines are separately defined.

The human-organoid task-manifest builder converts the draft OSD-863/OSD-871
source inventory into a draft `human_organoid_spaceflight` manifest and a
separate extension task index. It records a blocked LEO/ISS-versus-ground pilot
design, donor/organoid/microglia blocking requirements, the single-mission
limitation, sample-factor-backed draft folds, and normalized expression matrix
alignment status. It also records that public OSDR DE tables and contrast
definitions are available for future DE/signature metrics, and points to the
V9-ORG-014 derived DE reference plus `response_signature.csv` contract. It
still treats payload checksum verification as pending and keeps organoid
baseline language draft-only.

The human-organoid checksum audit uses the same OSDR file-list/checksum parser
as the public mouse bulk audit. It currently records API-ok file listings and
parsed raw/processed MD5 manifests for OSD-863 and OSD-871, while keeping
`freeze_ready=false` because payload files have not been downloaded and hashed.

The human-organoid sample-table audit finds and parses one OSDR SampleTable per
source. It currently records 19 sample-table rows for OSD-863 and 23 for
OSD-871, matching the 42 public samples reported for GSE259421. Manual condition
factor mapping is still required before fold generation.

The human-organoid sample-factor builder parses the compact OSDR `condition`
field into `disease_context`, `spaceflight_condition`, `true_label`, and
`microglia_condition`. The GEO metadata builder then augments those rows from
GSE259421 series-matrix metadata with `donor_or_line_id`, `ipsc_line_id`, sample
title, cell type, and treatment. The draft task manifest carries four default
sample-count-backed draft folds: hold out cortical organoids, hold out
dopaminergic organoids, hold out with-microglia samples, and hold out
without-microglia samples. GEO-derived donor/Subject holdouts are recorded as
diagnostic folds only because donor, source, organoid fate, and disease context
are not independently crossed.

The human-organoid expression-matrix audit discovers and downloads the canonical
OSDR `*_rna_seq_Normalized_Counts_GLbulkRNAseq.csv` matrices, excluding rRNArm
and unnormalized files. The current audit records one aligned matrix per source:
OSD-863 has 30,408 feature rows and 19/19 matched sample columns, while OSD-871
has 30,269 feature rows and 23/23 matched sample columns.

The human-organoid nearest-centroid pilot baseline loads both aligned matrices,
uses their 27,986 common human gene features, applies train-fold-only log1p,
variable-gene selection, and z-scoring, then evaluates the four
sample-count-backed draft folds. Current summary output is
`v9/human_organoid/reports/nearest_centroid/human_organoid_baseline_summary.csv`
with `release_status=draft_not_frozen` and
`claim_boundary=pilot_baseline_only_not_leaderboard`.

The human-organoid sensitivity runner evaluates 20 preprocessing variants across
log1p/no-transform, z-score/no-scaling, and five train-fold variable-gene
settings. The robustness review keeps the log1p/z-score/2,000-gene baseline as
the conservative default, despite higher raw-count sensitivity scores, because
raw-scale variants show worse calibration and need library-size/confounding
audits before interpretation.

The human-organoid diagnostics audit records that current OSDR SampleTable files
expose only `condition`, while GEO GSE259421 series metadata recovers
Subject/iPSC-line ids for all 42 public samples. Donor-aware folds remain
diagnostic-only because the design is small and coupled. The audit also shows
label-associated sample-total and sparsity differences, which keeps raw-count
sensitivity variants out of the default baseline.

The donor diagnostic runner evaluates the GEO-derived Subject/iPSC-line holdouts
without changing the default benchmark fold family. Current donor diagnostic
outputs stay under `v9/human_organoid/reports/donor_diagnostics/` with
`claim_boundary=donor_diagnostic_only_not_leaderboard`.

The human-organoid signature-reference audit records that OSDR lists public
differential-expression tables and contrast-definition CSVs for both organoid
sources. Each source has 56 parsed contrast pairs, including four direct
matched Ground Control versus Space Flight contrasts and four reversed matches.
The draft task manifest therefore marks `de_direction_match` and
`signature_rank_correlation` as reference-backed future metrics, but keeps
classification metrics primary. V9-ORG-014 adds the derived DE reference and
`response_signature.csv` contract. V9-ORG-015 adds diagnostic scoring when that
artifact is supplied, while keeping these metrics non-primary.

The response-signature smoke runner mirrors a small number of rows from the
real derived DE reference into an example `response_signature.csv` and evaluates
that fixture under `v9/human_organoid/reports/response_signature_smoke/`. The
result is intentionally a scorer plumbing check, not a model-performance or
biological claim.

The response-signature adapter design selects a conservative source-transfer
empirical baseline as the first real model-produced signature path. OSD-863
target contrasts should be predicted from OSD-871 training samples, and OSD-871
target contrasts should be predicted from OSD-863 training samples. This avoids
using target-source expression or DE references to generate the scored
signatures.

The source-transfer signature runner implements that design under
`v9/human_organoid/reports/source_transfer_signature/`. It writes compressed
`response_signature.csv.gz`, records
`reference_not_used_for_signature_generation`, and keeps the report
diagnostic-only.

The source-transfer diagnostic review keeps this baseline as the first
conservative response-signature diagnostic. Its direction score is stronger
than trivial sign baselines, but weak rank correlation keeps the claim narrow.
The next organoid signature block should design a per-condition source-transfer
adapter before promoting any more complex model signature.

The per-condition adapter design selects a microglia-matched source-transfer
signature as the next comparative diagnostic. Microglia status is crossed in
both organoid sources with Ground and LEO/ISS samples, while disease context is
only partially shared, so disease-matched signatures are deferred as partial
coverage diagnostics.

The microglia-matched source-transfer runner implements that comparison under
`v9/human_organoid/reports/microglia_source_transfer_signature/`. It keeps the
same non-leaky source-transfer boundary, emits 223,888 response-signature rows,
and reports `de_direction_match=0.7889125799573561` and
`signature_rank_correlation=0.1500722829461316`. Compared with the global
source-transfer diagnostic, direction agreement improves by 0.018239 but rank
correlation decreases by 0.025936, so the microglia-matched adapter remains a
secondary diagnostic rather than a replacement.

The source-transfer adapter comparison review keeps both adapters: global as the
first conservative diagnostic and microglia-matched as a secondary
condition-sensitivity diagnostic. It selects a partial shared-control
disease+microglia design as the next organoid signature block before any
classifier-coefficient signature work.

The shared-control adapter design limits disease+microglia matching to the four
`no_known_diseases` target contrasts where the opposite-source training stratum
exists. PPMS and sporadic Parkinson disease contrasts are explicitly skipped
rather than silently falling back to a broader signature.

The shared-control source-transfer runner implements that partial diagnostic
under `v9/human_organoid/reports/shared_control_source_transfer_signature/`. It
emits 111,944 response-signature rows for four shared-control contrasts, skips
four disease-specific contrasts in metadata, and reports partial-coverage
`de_direction_match=0.5899419729206963` plus
`signature_rank_correlation=0.03566701395309356`. The result is weaker than the
global and microglia-matched diagnostics on the same emitted contrasts, so
disease+microglia stratification is not promoted further.

The classifier feature-effect contract keeps model coefficients out of
`response_signature.csv`. It defines a separate optional `feature_effect.csv`
artifact and recommends L2 logistic gene-space coefficients as the first
implementation path, with rank/sign/top-k diagnostics kept non-primary.

The L2 logistic feature-effect pilot implements that contract under
`v9/human_organoid/reports/logistic_feature_effect/`. It emits 16,000
feature-effect rows from source-transfer training samples, records
`reference_not_used_for_effect_generation`, keeps `de_direction_match` skipped
because no `response_signature.csv` is supplied, and reports diagnostic
`feature_effect_direction_match=0.6078431372549019` plus
`feature_effect_rank_correlation=0.08672800238082004`. Top-k DE overlap is low:
1/50, 5/100, 10/250, and 14/500.

The feature-effect diagnostic review keeps `feature_effect.csv` as an optional
secondary diagnostic, but does not promote it. It recommends top-k/null
calibration before PCA-LR reconstructed gene-weight work.

The feature-effect top-k null calibration adds expected overlap, enrichment,
and exact hypergeometric upper-tail p-values to
`feature_effect_top_k_de_overlap`. After regenerating the L2 logistic report,
K=100, K=250, and K=500 show aggregate enrichment above a random top-k null,
while K=50 does not. The result strengthens the secondary diagnostic claim but
keeps feature effects out of primary benchmark metrics and out of
`response_signature.csv`.

The PCA-LR reconstructed feature-effect design defines the auditable gene-space
mapping `pca.components_.T @ logistic.coef_[0]` in train-standardized
transformed gene space. It approves one constrained pilot implementation as an
optional secondary diagnostic, with calibrated comparison against the L2
logistic feature-effect report required before any further promotion.

The PCA-LR reconstructed feature-effect pilot implements that design under
`v9/human_organoid/reports/pca_lr_feature_effect/`. It emits 16,000
feature-effect rows and reports
`feature_effect_direction_match=0.5980392156862745` plus
`feature_effect_rank_correlation=0.08668664748698189`. Aggregate calibrated
top-k overlap is identical to the L2 logistic report, while direction agreement
is slightly weaker, so the PCA-LR report is kept only as an optional negative
comparison.

The organoid diagnostic consolidation note defines the draft-alpha boundary:
source/task/provenance artifacts, the derived DE reference, the global
source-transfer response-signature diagnostic, and the calibrated L2 logistic
feature-effect diagnostic are the main organoid surface. Microglia-matched,
shared-control, PCA-LR reconstructed, donor-holdout, and smoke-fixture reports
remain secondary, exploratory, or plumbing-only. The active scientific lane now
returns to multispecies expansion.

The multispecies candidate source deep audit reviews local OSD-207, OSD-37, and
OSD-120 processed matrices and sample tables. OSD-37 is the first plant
go-after-source-audit candidate, OSD-207 is the first fly go-after-source-audit
candidate, and OSD-120 is deferred to a later Arabidopsis light/genotype
interaction task. Follow-up multispecies blocks parsed sample factors, audited
matrix-column alignment, generated checksum evidence, drafted OSD-37/OSD-207
species-native manifests, and ran a draft-only nearest-centroid feasibility
baseline. OSD-37 is currently the cleaner first plant feasibility example,
while OSD-207 remains valid but more condition-stratum heterogeneous. The
multispecies sensitivity grid keeps the conservative default at `log1p`,
train-fold `zscore`, and top 2,000 train-fold variable genes, rather than
chasing the highest raw-count score.

The OSD-120 interaction ladder now consolidates nearest centroid, dense L2
logistic, top-500 L2 sensitivity, sparse L1, and sparse L1 stability evidence
under `multispecies/reports/interaction_baseline_ladder/`. The ladder keeps
`sparse_l1_c1` as the leading transparent draft diagnostic candidate
(`primary BA=0.9167`, `secondary BA=0.8333`, `diagnostic BA=0.8889`, 9/11 folds
improved versus nearest centroid, 2/11 tied, 0/11 worse) and keeps
`sparse_l1_c0p3` as the compact stability comparator.

The OSD-120 diagnostic candidate package turns `sparse_l1_c1` into a compact
evidence bundle for manuscript or poster drafting. It emits one summary row,
three fragile-focus rows, 19 stable-feature rows selected in at least half of
balanced train-fold subsamples, and five claim-map rows that separate local
benchmark evidence from external OSD-120 CARA light/genotype context.

The OSD-120 figure/table package translates that evidence bundle into a
human-facing table surface. The main table has three rows: `Col.0.PhyD|Dark`
improves from 0.500 to 0.667 BA, `Light.Treatment` improves from 0.556 to
0.833 BA, and `Col.0.PhyD` improves from 0.500 to 0.917 BA. The appendix keeps
19 stable sparse-feature rows clearly labeled as model evidence, not validated
biomarkers.

The paired comparator table confirms that sparse L1 `C=0.3` reaches the same
three focus-fold balanced accuracies as `C=1.0` with fewer nonzero
coefficients, but `C=1.0` remains the primary figure candidate because it has
stronger mean/minimum family BA, no nearest-centroid-worse fold, and stronger
stable-feature evidence. `C=0.3` stays appendix/supplement only.

The OSD-120 diagnostic artifact manifest indexes the complete diagnostic
evidence set: 26 artifacts across task design, sparse L1 outputs, stability
audits, ladder, candidate package, figure/table package, paired comparator,
review notes, and validation tests. Its claim map links seven local claims to
artifact ids, artifact paths, tests, external context URLs where relevant, and
limitations.

The checksum-audit script queries the OSDR biological data API file-list
endpoint for each source, records API response hashes, discovers checksum
manifest-like files, downloads the manifest text, parses MD5 entries in the
OSDR formats observed so far, and records match counts against the OSDR file
listing. The current audit covers all 22 public bulk source rows, finds 39
checksum manifests, parses 8,439 checksum entries, and matches 8,275 entries
back to listed OSDR payload names.

The Data Package draft builder writes `v9/datapackage.draft.json`. It describes
11 current resources: metadata indexes, source/provenance tables, baseline
summary tables, eight task manifests, 24 prediction CSVs, 24 metrics JSON files,
and 24 run-manifest JSON files. It deliberately excludes large fold payload
hashing until a dedicated public payload manifest block.

## Interpretation

- `checksum_status=legacy_task_source_unfrozen` means these manifests preserve
  the current legacy source identity but are not yet a v9 frozen release.
- `freeze_ready=false` in `source_checksum_audit.csv` means checksum-manifest
  evidence exists, but payload files have not yet been downloaded and hashed in
  this audit.
- The public benchmark core does not depend on gated human datasets.
- Intervention or Mars-regime claims are intentionally out of scope for these
  bulk LOMO task manifests.

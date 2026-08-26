# Project Results Deep-Read Audit

Date: 2026-05-31

Purpose: record what has actually been read, parsed, and verified while
preparing the v1-v9 slide-deck source inventory. This is a companion to
`docs/PROJECT_RESULTS_LOCATION_INVENTORY_2026_05_31.md`.

## Honest Scope Statement

The first inventory pass was a strong location-and-summary pass, not a literal
line-by-line reading of every result file. After review, this audit upgrades the
evidence base by separating three levels:

1. Manually read narrative/core result documents.
2. Programmatically parsed machine-readable JSON/CSV outputs.
3. Still-not-done visual or full provenance checks.

Current status:

- Core narrative/result documents were read and keyword-checked across 7,162
  lines.
- Machine-readable result surfaces were parsed across 1,531 JSON/CSV/Markdown
  files.
- All 1,087 JSON files parsed successfully.
- All 377 CSV files parsed successfully.
- 67 Markdown report files were counted/inventoried inside result directories.
- Figures were inventoried by file path, but not visually rendered/reviewed in
  this pass.
- Analyses were not rerun; this is a source-reading and result-surface audit.

## Manually Read Core Documents

These files were used as the primary narrative/result source layer:

| File | Lines | Read role |
|---|---:|---|
| `README.md` | 419 | project overview, v1-v7 scope, public-facing table |
| `docs/V1_PAPER_CONTENT.md` | 628 | original v1 paper narrative |
| `evaluation/RESULTS_SUMMARY.md` | 503 | root/v1-v5 result summary and canonical v4 table |
| `v2/evaluation/V2_RESULTS_SUMMARY.md` | 406 | temporal, human, PBMC, RRRM-1, LLM summary |
| `v3/README.md` | 105 | multi-species, spatial, RRRM-2, FM summary |
| `docs/CANONICAL_RESULTS_V7_1.md` | 118 | public-safe v7.1 scope and claim boundaries |
| `docs/V7_V8_CLOSURE_STATUS_2026_05_10.md` | 111 | v7/v8 validation and closure status |
| `v8/RESULTS_SUMMARY.md` | 37 | v8 bridge/decompose/intervene summary |
| `v9/README.md` | 837 | v9 platform and extension inventory |
| `v9/OPERATING_BACKLOG.md` | 3,826 | v9 workstream status, claims, and next block |
| `docs/V9_PURPOSE_DRIFT_AUDIT_2026_05_26.md` | 172 | v9 mission-drift and claim-boundary audit |

Total manually read/checkpointed narrative lines: 7,162.

## Machine-Readable Result Parse

Result directories scanned:

| Directory | JSON/CSV/MD files parsed or inventoried |
|---|---:|
| `evaluation/` | 121 |
| `v2/evaluation/` | 46 |
| `v3/evaluation/` | 19 |
| `v4/evaluation/` | 302 |
| `v5/evaluation/` | 25 |
| `v6/evaluation/` | 6 |
| `v7/evaluation/` | 24 |
| `v8/bridge/evaluation/` | 10 |
| `v8/decompose/evaluation/` | 18 |
| `v8/intervene/evaluation/` | 34 |
| `v8/causal/evaluation/` | 4 |
| `v8/multiomics/evaluation/` | 2 |
| `v9/reports/` | 106 |
| `v9/human_organoid/` | 153 |
| `v9/multispecies/` | 626 |
| `v9/sc_spaceflight/` | 35 |

Parse summary:

| Type | Count | Parse result |
|---|---:|---|
| JSON | 1,087 | 1,087 OK, 0 failed |
| CSV | 377 | 377 OK, 0 failed |
| Markdown reports in result dirs | 67 | inventoried |
| Total | 1,531 | scanned/inventoried |

Largest parsed CSVs:

| Rows | Columns | File |
|---:|---:|---|
| 103,022 | 8 | `v8/causal/evaluation/icp_stressor_pathway_scores.csv` |
| 63,140 | 20 | `v9/human_organoid/matrices/GLDS-716_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv` |
| 63,140 | 20 | `v9/human_organoid/matrices/GLDS-716_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv` |
| 63,140 | 24 | `v9/human_organoid/matrices/GLDS-720_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv` |
| 63,140 | 24 | `v9/human_organoid/matrices/GLDS-720_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv` |
| 30,408 | 20 | `v9/human_organoid/matrices/GLDS-716_rna_seq_Normalized_Counts_GLbulkRNAseq.csv` |
| 30,269 | 24 | `v9/human_organoid/matrices/GLDS-720_rna_seq_Normalized_Counts_GLbulkRNAseq.csv` |
| 30,000 | 16 | `v9/multispecies/reports/interaction_logistic_sparse_l1/feature_coefficient_audit.csv` |

## Important Finding: v4 Canonical vs Raw Summary Difference

The public-facing v4 table in `evaluation/RESULTS_SUMMARY.md`,
`README.md`, and `docs/CANONICAL_RESULTS_V7_1.md` is not identical to the raw
best-row extraction from `v4/evaluation/M1_summary.json`.

Public canonical table:

| Tissue | AUROC | Method | Feature |
|---|---:|---|---|
| Thymus | 0.948 | PCA-LR | KEGG |
| Colon | 0.921 | PCA-LR | KEGG |
| Lung | 0.901 | PCA-LR | Gene |
| Kidney | 0.829 | ElasticNet-LR | Hallmark |
| Eye | 0.823 | PCA-LR | Hallmark |
| Skin | 0.819 | PCA-LR | Gene |
| Gastrocnemius | 0.776 | PCA-LR | Gene |
| Liver | 0.670 | PCA-LR | Gene |

Raw best-row extraction from `v4/evaluation/M1_summary.json`:

| Tissue | AUROC | Method | Feature | perm_p | CV |
|---|---:|---|---|---:|---|
| Liver | 0.7659 | ElasticNet-LR | KEGG | 0.0931 | LOMO |
| Gastrocnemius | 0.8981 | ElasticNet-LR | Gene | 0.0576 | LOMO |
| Kidney | 0.8291 | PCA-LR | Hallmark | 0.0097 | LOMO |
| Thymus | 0.9483 | PCA-LR | KEGG | 0.0314 | LOMO |
| Eye | 0.8227 | PCA-LR | Hallmark | 0.0416 | LOMO |
| Skin | 0.8186 | ElasticNet-LR | Gene | 0.0037 | LOMO |
| Lung | 0.9009 | ElasticNet-LR | Gene | 0.0284 | 5-fold stratified |
| Colon | 0.9214 | PCA-LR | KEGG | 0.0328 | 5-fold stratified |

Interpretation:

- For deck writing, the current public-safe source remains
  `docs/CANONICAL_RESULTS_V7_1.md`.
- Before finalizing a slide that says "best method per tissue", resolve whether
  the deck should use the canonical table, the raw `M1_summary.json` best row,
  or a carefully captioned hybrid.
- This is the most important inconsistency found by the deeper pass.

## Version-Level Verified Notes

### Root / v1

Verified source layer:

- `README.md`
- `docs/V1_PAPER_CONTENT.md`
- `evaluation/RESULTS_SUMMARY.md`
- root `evaluation/*.json`

Verified claims:

- v1-v7 public benchmark scope is 8 tissues, 24+ OSD accessions, and 600+
  binary/control samples in the full release layer.
- v1 LOMO/FM core is a 6-tissue subset.
- Category B cross-mission transfer ranks thymus above liver: thymus 0.860 vs
  liver 0.577 in `evaluation/RESULTS_SUMMARY.md`.
- Category A gene-level rows include thymus 0.923, skin 0.821, gastrocnemius
  0.824, eye 0.789, liver 0.670, kidney 0.432 in
  `evaluation/RESULTS_SUMMARY.md`.
- Foundation-model/LLM summary remains: PCA-LR reference 0.758, scGPT 0.666,
  Mouse-Geneformer 0.476, UCE/scFoundation best single-tissue v3 rows below the
  classical benchmark surface, text LLMs near chance.

### v2

Verified source layer:

- `v2/evaluation/V2_RESULTS_SUMMARY.md`
- `v2/evaluation/*.json`
- `v2/processed/`
- `v2/figures/`

Verified claims:

- T1 preservation artifact is explicit: GC and BSL timing separations are often
  as strong as or stronger than FLT timing separation.
- T2 recovery ratios are RR-6 0.842 and RR-8 0.652.
- T3 age-stratified RR-8 liver spaceflight detection is OLD 0.945 vs YNG 0.679.
- Zero-shot text LLM means remain near chance: Llama 0.485, DeepSeek 0.507,
  Gemini 0.493 end-to-end.
- E1 mouse liver to human JAXA cfRNA Hallmark NES Spearman r is 0.352.
- RRRM-1 scRNA-seq hardening reports 38,081 cells.
- F2-D reports strong T-cell cross-species pathway conservation in selected
  mouse RRRM-1 to human I4 PBMC comparisons.

### v3

Verified source layer:

- `v3/README.md`
- `v3/evaluation/*.json`

Verified claims:

- E4 includes mouse and Drosophila NES analysis; Drosophila-mouse correlations
  are negative in the reported README summary.
- F3 spatial Visium brain is a negative result: section AUROC 0.139 and animal
  AUROC 0.444.
- F5 RRRM-2 PBMC cell-type signal includes NK AUROC 0.845.
- Bone marrow is reported near chance.
- UCE and scFoundation do not beat the classical baseline surface.

### v4

Verified source layer:

- `evaluation/RESULTS_SUMMARY.md`
- `docs/CANONICAL_RESULTS_V7_1.md`
- `v4/evaluation/M1_summary.json`
- `v4/evaluation/` parse sweep

Verified claims:

- v4 evaluation grid is 8 tissues x 8 classifiers x 4 feature types = 256
  evaluations.
- Canonical public table and raw `M1_summary.json` need explicit separation, as
  described above.
- Public-safe v4 deck text should also include the canonical note that 40/256
  configurations are significant at p<0.05, and 6/8 tissues have at least one
  significant configuration.

### v5

Verified source layer:

- `v5/evaluation/*.json`
- `evaluation/RESULTS_SUMMARY.md`

Parsed immune deconvolution counts:

| Tissue | n_samples | Cell types | FDR<0.05 |
|---|---:|---:|---:|
| Colon | 54 | 16 | 0 |
| Eye | 41 | 13 | 0 |
| Gastrocnemius | 32 | 11 | 0 |
| Kidney | 142 | 14 | 2 |
| Liver | 261 | 14 | 0 |
| Lung | 58 | 16 | 0 |
| Skin | 118 | 14 | 6 |
| Thymus | 86 | 16 | 2 |

Parsed TF activity counts:

| Tissue | n_samples | TFs tested | FDR<0.05 |
|---|---:|---:|---:|
| Colon | 54 | 758 | 0 |
| Eye | 41 | 746 | 0 |
| Gastrocnemius | 32 | 710 | 0 |
| Kidney | 142 | 734 | 177 |
| Liver | 261 | 708 | 105 |
| Lung | 58 | 758 | 0 |
| Skin | 118 | 740 | 241 |
| Thymus | 86 | 756 | 240 |

Other v5 parsed claims:

- `v5/evaluation/consensus_biomarker_panel.json`: panel size 20; 1,919
  candidates scored.
- Top panel genes include MUP22, Thrsp/THRSP, Apoa1, Top2a.
- Panel validation AUROCs are strongest in gastrocnemius, liver, eye, and
  kidney among the reported tissue rows.
- `v5/evaluation/drug_targets.json` and related v6 drug validation support
  hypothesis-prioritization only, not intervention efficacy.

### v6

Verified source layer:

- `v6/evaluation/*.json`

Parsed claims:

- `V6_A_gene_conservation.json`: universe size 14,633; DRR genes in universe
  389.
- `V6_B_pathway_conservation.json`: 50 human pathways, 5 tissues analyzed,
  mean rho 0.285.
- `V6_C_cross_species_transfer.json`: ortholog map size 17,103; cfRNA genes
  500.
- `V6_D_biomarker_validation.json`: panel size 20; 15/20 detected in cfRNA;
  detection rate 0.75; no DE FDR<0.05 panel genes in that setting.
- `V6_E_tf_conservation.json`: 728 human TFs, 60 significant human TFs,
  8 tissues analyzed, mean rho 0.0298.
- `V6_F_drug_target_validation.json`: 207 drug-target genes, 196 in cfRNA,
  tier counts A=3, B=7, C=186, D=11.

### v7 / v7.1

Verified source layer:

- `docs/CANONICAL_RESULTS_V7_1.md`
- `docs/V7_V8_CLOSURE_STATUS_2026_05_10.md`
- `v7/evaluation/*.json`

Parsed claims:

- GNN rows do not overturn the classical-baseline conclusion. A few liver GNN
  rows slightly exceed the local PCA-LR comparison, but the overall result does
  not change the headline.
- scPRINT-2 rows are below local PCA-LR comparisons:
  - eye 0.3095, delta -0.3877
  - gastrocnemius 0.5880, delta -0.2639
  - kidney 0.4842, delta -0.0996
  - liver 0.5051, delta -0.1776
  - skin 0.6283, delta -0.1678
  - thymus 0.6639, delta -0.2444
- v7/v8 closure note reports clean-checkout validation at commit
  `b486aee77ae956b108a59ac762ebeb7b302e7928` with 47/47 tests passing.

### v8

Verified source layer:

- `v8/RESULTS_SUMMARY.md`
- `v8/bridge/evaluation/*.json`
- `v8/decompose/evaluation/*.json`
- `v8/intervene/evaluation/*.json`
- `v8/causal/evaluation/*.csv/json`
- `v8/provenance/`

Parsed claims:

- `v8/bridge/evaluation/supervised_conservation.json`: RF base AUROC 0.7123,
  augmented AUROC 0.8873, delta 0.1749 with CI approximately 0.1342 to 0.2189.
- `v8/decompose/evaluation/factorial_flight_decomposition.json` covers
  hze_endocrine, spleen, skin_analog, and brain.
- `v8/decompose/evaluation/mars_saturation_summary.json` records linear and
  bounded Mars-regime sensitivity logic.
- `v8/intervene/evaluation/offline_reversal_summary.json`,
  `crispr_orthog_summary.json`, and `safety_triage_summary.json` support
  hypothesis-generation only.
- `v8/causal/evaluation/icp_stressor_pathway_scores.csv` is the largest v8 CSV
  surface parsed here: 103,022 rows.

v8 claim boundary:

- Bridge/stressor/intervention/Mars outputs must remain incubator or
  hypothesis-generation claims unless beta freeze/recompute/provenance gates
  are completed.

### v9

Verified source layer:

- `v9/README.md`
- `v9/OPERATING_BACKLOG.md`
- `docs/V9_PURPOSE_DRIFT_AUDIT_2026_05_26.md`
- `v9/reports/`
- `v9/human_organoid/`
- `v9/multispecies/`
- `v9/sc_spaceflight/`

Public bulk:

- `v9/task_manifests/` has 8 public bulk LOMO manifests.
- `v9/reports/bulk_lomo_baseline_summary.json` has 24 rows across
  nearest-centroid, PCA-LR, and L2 logistic baselines.
- `v9/source_inventory.csv/json` records 22 public bulk source rows.
- Current alpha boundary is metadata-only because 0/22 local payload hashes are
  verified for release-style payload freeze.

Human organoid:

- OSD-863 and OSD-871 are present as draft human neural organoid extension
  sources.
- GEO GSE259421 sample metadata recovers Subject/iPSC-line metadata for 42/42
  public samples according to the v9 README/backlog and audit artifacts.
- Normalized matrix surfaces include OSD-863 30,408 feature rows x 20 columns
  and OSD-871 30,269 feature rows x 24 columns in the parsed matrices.
- Derived DE reference contains 242,708 rows, with 2,368 FDR<=0.05 rows in the
  v9 backlog summary.
- `source_transfer_signature/metrics.json` reports 223,888 joined response
  signature rows and diagnostic `de_direction_match` above trivial baseline,
  but it is non-primary.
- `logistic_feature_effect/metrics.json` emits 16,000 feature-effect rows and
  remains a secondary diagnostic.

Multispecies:

- OSD-37 Arabidopsis and OSD-207 Drosophila are species-native draft task
  candidates.
- OSD-120 is a separate light/genotype interaction diagnostic branch.
- `v9/multispecies/reports/interaction_logistic_sparse_l1/tvg2000_log1p_zscore_l1_c1/multispecies_baseline_summary.json`
  has 3 fold-family rows:
  - primary genotype/ecotype holdout balanced accuracy 0.9167
  - secondary light-treatment holdout balanced accuracy 0.8333
  - diagnostic condition-stratum holdout balanced accuracy 0.8889
- OSD-120 remains draft diagnostic metadata, not a leaderboard or frozen public
  release.

Single-cell:

- `v9/sc_spaceflight/anndata_manifest_draft_summary.json` records OSD-918
  RRRM-1 blood, 8 samples, 4 Flight, 4 Ground, 4,395 QC cells, and 19,064
  genes.
- `v9/sc_spaceflight/obs_var_audit_summary.json` records a skip-only audit:
  27 skipped requirements, 17 blockers, missing canonical h5ad payload, no
  checksum, and no score claim.
- `v9/sc_spaceflight/osdr_processed_payload_discovery_summary.json` records
  19 official OSD-918 files, including 16 raw FASTQ files covering 8/8 expected
  blood SRX pairs, but no processed h5ad, STARsolo bundle, or processed
  checksum manifest.
- Active next block is V9-SC-006d owner scratch path intake or raw FASTQ
  regeneration feasibility.

## What Is Still Not Fully Reviewed

Not yet done in this audit:

- Visual rendering/QA of every HTML, PNG, and PDF figure.
- Manual line-by-line reading of all 1,531 machine-readable/result report files.
  Instead, JSON/CSV files were parser-validated and key fields were extracted.
- Re-running analyses from raw data.
- Verifying all external-source claims against live OSDR/GEO/NASA pages.
- Resolving the v4 canonical table versus raw `M1_summary.json` discrepancy.

## Deck Consequence

For the first slide deck draft:

- Use `docs/CANONICAL_RESULTS_V7_1.md` as the public-safe v1-v7 claim source.
- Add an internal speaker-note caveat that raw v4 `M1_summary.json` best rows
  differ from the canonical table and need a deliberate decision.
- Use v8 only as a translational incubator.
- Use v9 as platformization and alpha scaffold, not a frozen benchmark release.
- Keep organoid/multispecies/single-cell as future-facing or diagnostic tracks
  unless a later pass promotes them.

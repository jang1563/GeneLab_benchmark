# Project Slide Content Inventory: v1 To v9

Date: 2026-05-31

Purpose: slide-deck preparation source note. This file summarizes the project
story, version-by-version outputs, and headline results before deck formatting
or visual design work.

## One-Sentence Project Story

GeneLab Benchmark started as a mission-held-out mouse spaceflight
transcriptomics benchmark and has grown into SpaceBio-Bench: a provenance-first
platform for testing biological AI under spaceflight domain shift across
missions, tissues, species, modalities, organoids, single-cell tasks, and
stressor regimes.

## Narrative Arc For Slides

1. Benchmark problem: NASA OSDR contains many spaceflight omics studies, but the
   field lacked a standardized way to ask whether models generalize to unseen
   missions.
2. v1-v3 answered the scientific benchmark question: spaceflight signatures are
   tissue-specific; thymus generalizes better than liver; pathway biology
   often transfers better than raw genes; foundation models and text LLMs do
   not beat tuned classical baselines.
3. v4-v7 hardened the benchmark: broader method sweep, biological mechanism
   layers, translation to human data, graph/foundation-model stress tests, and
   public result-surface cleanup.
4. v8 moved from benchmark to translational hypothesis generation:
   rodent-human pathway bridge, stressor decomposition, intervention
   prioritization, causal evidence maps, and provenance-controlled beta
   packaging.
5. v9 reframed the work as a platform: typed manifests, metric profiles,
   registries, public metadata alpha, organoid and multispecies extensions, and
   the first single-cell flagship scaffold with explicit no-overclaim
   boundaries.

## Version Timeline

| Version | Main Goal | Main Outputs | Headline Result |
|---|---|---|---|
| v1 | Build the original mouse bulk RNA-seq cross-mission benchmark | root `tasks/`, `processed/`, `evaluation/`, `docs/V1_PAPER_CONTENT.md` | Thymus, not liver, is the most generalizable tissue; all FMs/LLMs underperform classical ML. |
| v2 | Add temporal, cross-species, and first single-cell resolution | `v2/processed`, `v2/evaluation`, `v2/figures`, RRRM-1 h5ad pipeline docs | Preservation artifacts, recovery overshoot, age-amplified response, human PBMC/cfRNA links, RRRM-1 scRNA-seq tasks. |
| v3 | Add multi-species, spatial Visium, RRRM-2, and more FMs | `v3/evaluation`, `v3/figures`, `v3/scripts` | Brain spatial is a genuine negative; RRRM-2 PBMC/NK has strong signal; UCE/scFoundation still underperform PCA-LR. |
| v4 | Run definitive multi-method benchmark and network biology | `v4/evaluation`, `v4/wgcna_outputs`, v4 figures | 8 tissues x 8 classifiers x 4 feature types = 256 evals; PCA-LR best, pathway features rescue selected tissues. |
| v5 | Add systems-biology interpretation | `v5/evaluation`, v5 figures | Skin has strongest immune deconvolution signal; TF activity and metabolic/drug-target layers make results biologically interpretable. |
| v6 | Test mouse-to-human translational conservation | `v6/evaluation` | Gene-level transfer is weak/negative; pathway-level conservation is partial; drug target validation yields Tier A candidates. |
| v7 / v7.1 | Unify benchmark and test newer FM/GNN ideas; clean public release surface | `v7/evaluation`, `docs/CANONICAL_RESULTS_V7_1.md` | scPRINT-2 and WGCNA-GNN variants do not overturn the classical-baseline conclusion; v7.1 locks public counts and claims. |
| v8 | SpaceMed translational incubator | `v8/bridge`, `v8/decompose`, `v8/intervene`, `v8/causal`, `v8/provenance` | Mouse NES improves human I4/Twins pathway prediction; stressor extrapolation is fragile; perturbation hits are hypothesis-only. |
| v9 | Platformize as SpaceBio-Bench | `spacebio_bench`, `v9`, `docs/V9_*` | Public bulk metadata alpha, human organoid and multispecies draft tracks, and single-cell RRRM scaffold now exist with manifest-backed claim boundaries. |

## v1: Original GeneLab Benchmark

Core question: can a model trained on some spaceflight missions classify
Flight vs Ground samples from a completely unseen mission?

Scope:

- 6 primary tissues: liver, gastrocnemius, kidney, thymus, skin, eye.
- 17 mission labels and roughly 450 binary samples in the original LOMO/FM
  core.
- Categories: A spaceflight detection, B cross-mission transfer, C
  cross-tissue transfer, D confounder prediction, J gene-vs-pathway comparison,
  negative controls, and held-out validation.
- Model tiers: classical ML, gene-expression foundation models, and text LLMs.

Key results:

- Category B transfer ranks thymus first: thymus mean AUROC 0.860 versus liver
  0.577; thymus vs liver p=0.001.
- Category A LOMO gene-level results: thymus 0.923, gastrocnemius 0.907, skin
  0.821, eye 0.811, liver 0.653, kidney 0.593.
- Pathway rescue: eye improves to 0.915 and kidney improves substantially
  from weak gene-level performance.
- Gene-level features perfectly encode mission identity in liver
  (macro-F1=1.000), while pathway features resist mission/batch detection
  (macro-F1 about 0.056).
- Foundation models and text LLMs fail to beat classical ML: PCA-LR/classical
  mean about 0.758, scGPT 0.666, Mouse-Geneformer 0.476, text LLMs near
  chance.
- Held-out validation supports generalization: RR-23 thymus AUROC 0.905 and
  RR-7 skin AUROC 0.885.

Slide angle:

- "The benchmark overturns the liver-centric assumption."
- "Mission-held-out evaluation is the key methodological contribution."
- "Pathway abstraction is a batch-resistant biological feature layer."

## v2: Temporal, Cross-Species, And RRRM-1 Single-Cell

Core goal: extend the benchmark beyond static mouse bulk classification into
temporal recovery, human comparison, and cell-type-resolved spaceflight biology.

Key results:

- T1: ISS-Terminal vs Live Animal Return differences are dominated by
  preservation method; GC samples separate as strongly as FLT samples.
- T2: RR-8 post-return liver profiles show recovery/overshoot. PCA recovery
  ratio is 0.652 and 25/27 Hallmark pathways recover or overshoot.
- T3: spaceflight amplifies aging. OLD RR-8 liver AUROC is 0.945 versus YNG
  0.679, delta +0.266.
- E1/E2: mouse liver pathway NES partially conserves with human JAXA cfRNA
  (Spearman r=0.352); short-duration I4 cfRNA does not show the same
  conservation.
- F1: I4 PBMC snRNA-seq shows cell-type-specific response; innate immune
  populations are most pathway-responsive and R+45 can be more divergent than
  R+1.
- F2: RRRM-1 scRNA-seq pipeline completed for blood, eye, muscle, and skin:
  38,081 cells after hardening, with broad cell-type annotation and benchmark
  task definitions.

Slide angle:

- "Spaceflight signal is entangled with handling and recovery biology."
- "Single-cell resolution reveals which immune/cell-type compartments carry
  the pathway signal."

## v3: Multi-Species, Spatial, RRRM-2, And Additional FMs

Core goal: stress-test the benchmark across species, spatial transcriptomics,
RRRM-2 single-cell data, radiation analogs, and newer foundation models.

Key results:

- E4 multi-species: mouse intra-species pathway correlations show tissue-pair
  structure; Drosophila-mouse comparisons are negative, showing nontrivial
  cross-species directionality.
- F3 spatial Visium brain is a strong negative result: section-level AUROC
  0.139, animal-level 0.444, PC1 dominated by slide batch rather than
  condition.
- F5 RRRM-2 scRNA-seq: PBMC/NK cell signal is strong (NK AUROC 0.845), while
  bone marrow cell types are near chance.
- Extended tissues: colon is strong (about 0.900), lung weaker, and skin
  extended analyses remain moderate.
- UCE and scFoundation do not beat PCA-LR. Best rows are UCE thymus 0.632 and
  scFoundation liver 0.635/gastrocnemius 0.691, still below the classical
  benchmark surface.

Slide angle:

- "Negative controls and negative biological results are part of the benchmark,
  not failures."
- "Cell atlas FMs still do not solve small-n bulk spaceflight generalization."

## v4: Multi-Method Benchmark Hardening

Core goal: turn the original benchmark into a broad method comparison and
biology-aware result surface.

Scope:

- 8 tissues: liver, gastrocnemius, kidney, thymus, skin, eye, lung, colon.
- 8 classifiers: PCA-LR, ElasticNet-LR, Random Forest, XGBoost, SVM-linear,
  SVM-RBF, TabNet, LightGBM.
- 4 feature types: gene, Hallmark, KEGG, pathway-combined.
- 256 evaluations.

Headline table:

- Thymus 0.948 with PCA-LR + KEGG.
- Colon 0.921 with PCA-LR + KEGG.
- Lung 0.901 with PCA-LR + gene.
- Kidney 0.829 with ElasticNet-LR + Hallmark.
- Eye 0.823 with PCA-LR + Hallmark.
- Skin 0.819 with PCA-LR + gene.
- Gastrocnemius 0.776 with PCA-LR + gene.
- Liver 0.670 with PCA-LR + gene.

Key result:

- PCA-LR is best overall, with 8-tissue gene-level mean AUROC 0.776.
- ElasticNet-LR is second at 0.762.
- 40/256 configurations are significant at p<0.05.
- Deep learning/kernel variants do not dominate; small-n transcriptomics favors
  simple regularized baselines.

Slide angle:

- "The conclusion is not cherry-picked: it survives an 8x8x4 method grid."

## v5: Systems-Biology Interpretation

Core goal: explain what the benchmark has learned biologically.

Key results:

- Immune deconvolution: skin has the strongest signal with 6/14 cell types FDR
  < 0.05; kidney and thymus each have 2/14; most other tissues have none.
- TF activity: thymus 240 significant TFs, skin 241, kidney 177, liver 105;
  colon/eye/gastrocnemius/lung have no significant TF rows in the current
  summary.
- Metabolic flux: E-Flux/pFBA outputs exist for six tissues and support
  tissue-specific metabolic interpretation.
- Drug target mapping: 1,919 mouse target genes, 834 mapped to human, 271 with
  DGIdb interactions, and 3,154 total DGIdb interactions.
- Consensus 20-gene biomarker panel: top panel genes include MUP22, Thrsp,
  Apoa1, and Top2a; panel AUROC is highest in gastrocnemius (0.806), then
  liver (0.754), eye (0.728), kidney (0.704), and lower in thymus/skin/lung/
  colon.

Slide angle:

- "v5 converts prediction into biological explanation: immune, TF, metabolic,
  druggability, and biomarker layers."

## v6: Mouse-To-Human Translational Validation

Core goal: test whether mouse benchmark signals translate to human spaceflight
data and drug-target evidence.

Key results:

- Gene conservation: 389 DRR genes in a universe of 14,633; gene-level transfer
  is generally weak.
- Pathway conservation: 50 human pathways, five tissues analyzed, mean rho
  0.285.
- Cross-species transfer uses 17,103 ortholog mappings and cfRNA gene features.
- Biomarker validation: 20-gene panel, 15 detected in cfRNA, but no FDR<0.05
  DE genes in the panel.
- TF conservation: 728 human TFs, 60 significant human TFs, mean rho 0.0298;
  limited direct conservation.
- Drug target validation: 207 drug-target genes, 196 detected in cfRNA, 3 Tier
  A validated genes and 7 Tier B promising genes.

Slide angle:

- "Translation is partial: pathways and target classes carry better signal than
  direct gene-level transfer."

## v7 And v7.1: Unified Benchmark/Future-FM Cleanup

Core goal: unify public v1-v7 results, test additional AI ideas, and clean up
the public benchmark release surface.

Key results:

- scPRINT-2 embeddings do not beat PCA-LR: best tissue rows include thymus
  0.664, skin 0.628, gastrocnemius 0.588, liver 0.505, kidney 0.484, eye 0.310;
  all remain below the paired PCA-LR reference.
- GNN/WGCNA graph variants are mixed and do not overturn classical baselines.
  Skin WGCNA GNN reaches 0.728, liver WGCNA 0.701, kidney WGCNA 0.551, but
  thymus WGCNA 0.644 is below the random graph 0.763.
- v7.1 canonical result file locks the public scope: full v1-v7 release has 8
  tissues, 24+ OSD accessions, 600+ samples, and 256 v4 multi-method
  evaluations.
- v7.1 also separates the v1-v7 benchmark paper from v8 translational claims.

Slide angle:

- "Even newer model families do not erase the core lesson: classical baselines
  and careful task design matter more than model fashion."

## v8: SpaceMed Translational Incubator

Core goal: move from detection benchmark to translational hypothesis generation
while keeping claim boundaries strict.

Pillars:

- BRIDGE: rodent pathway signals to human spaceflight data.
- DECOMPOSE: stressor factorial and Mars-regime extrapolation sensitivity.
- INTERVENE: perturbation reversal hypotheses from LINCS/CRISPR signatures.
- CAUSAL: invariant-causal-prediction evidence map.

Key results:

- BRIDGE: adding mouse tissue NES features improves I4/Twins supervised pathway
  conservation prediction. RF AUROC goes from 0.712 to 0.888, delta +0.175.
- Tissue bridge: gastrocnemius is uniquely positive against I4 PBMC in the
  refined full-MSigDB bridge; many solid-tissue programs invert in circulating
  compartments.
- DECOMPOSE: interaction terms dominate 44-61% variance in top-responsive genes;
  low-LET gamma and high-LET HZE can have opposite-sign concordance.
- Mars extrapolation: linear extrapolation breaks at high dose amplification;
  bounded sensitivity separates robust flags from saturation-sensitive flags.
- INTERVENE: LINCS/CRISPR triage nominates mechanism axes such as CDK, HSP90,
  AMPK/BMP, MEK/mTOR; CGP-60474 and quinacrine sit on the current Pareto front.
- All countermeasure/intervention language is hypothesis-only, not clinical or
  operational recommendation.

Release/provenance state:

- v8 alpha gate passed in clean HPC validation.
- v8 beta is not frozen: external source versions, API archives, one-command
  recomputation, provenance validation, and public artifact split remain
  release blockers.

Slide angle:

- "v8 is where the project becomes SpaceMed: useful translational hypotheses,
  but no overclaiming."

## v9: SpaceBio-Bench Platform

Core goal: turn the accumulated work into a public, living,
provenance-first benchmark platform.

Platform spine:

- `spacebio_bench` package skeleton.
- Task manifests, metric profiles, task registry, submission validation,
  evaluator, run manifests, baseline runners.
- Eight public bulk LOMO task manifests.
- 24 baseline rows across nearest-centroid, PCA-LR, and L2 logistic baselines.
- Source inventory: 22 public bulk source rows.
- Source checksum audit: 22/22 API-ok and checksum-manifest-parsed; 39 checksum
  manifests, 8,439 parsed MD5 entries, 8,275 payload-name matches.
- Data Package draft and Hugging Face-style metadata-alpha dataset card.

Public bulk alpha:

- Public bulk alpha freeze-gap matrix: 6 pass rows, 2 blockers, 2 needs-update.
- Main blocker: 0/22 public bulk sources have local payload-hash verification.
- Decision: metadata-only alpha snapshot, no frozen payload release.
- Data Package has 21 resources, including 10 alpha-boundary metadata tables.

Human organoid track:

- Public human neural organoid sources OSD-863 and OSD-871.
- 42 public GSE259421 samples; OSD-863 has 19 rows, OSD-871 has 23 rows.
- Normalized matrix audit: OSD-863 30,408 features with 19/19 samples matched;
  OSD-871 30,269 features with 23/23 samples matched.
- Derived DE reference: eight direct Ground vs Space Flight contrasts, 242,708
  gene/contrast rows, 2,368 FDR<=0.05 rows.
- Diagnostic response signatures and feature-effect reports are useful but
  remain non-primary and non-leaderboard.

Multispecies track:

- Source inventory for OSD-207, OSD-37, and OSD-120.
- OSD-37 Arabidopsis and OSD-207 Drosophila species-native task manifests.
- OSD-120 Arabidopsis root light/genotype interaction task developed as a
  diagnostic Phase 4A case study.
- OSD-120 model ladder advances sparse L1 C=1.0 as the leading transparent
  draft diagnostic candidate: primary BA 0.9167, secondary BA 0.8333,
  diagnostic BA 0.8889, 9 improved / 2 tied / 0 worse versus nearest centroid.
- OSD-120 release branch is closed at diagnostic metadata only; no DOI,
  archive identifier, license, creator metadata, or leaderboard promotion.

Single-cell flagship scaffold:

- V9-SC-001: 54 RRRM/single-cell legacy asset paths indexed, including 41
  RRRM-1 and 13 RRRM-2 paths; no local h5ad/loom/MTX payload found.
- V9-SC-002: draft non-runnable RRRM-1 blood `sc_spaceflight` manifest for
  OSD-918, 8 samples, 4 Flight, 4 Ground, 4,395 QC cells, 19,064 genes.
- V9-SC-003: `genelab_sc` metric specification with six metrics, seven input
  contracts, and six skip-policy rows.
- V9-SC-004: payload staging plan fixes canonical future h5ad target and 27
  obs/var/uns/matrix audit requirements.
- V9-SC-005: skip-aware AnnData auditor records the canonical h5ad as missing,
  with 27 skipped requirement rows and 17 blockers.
- Next block: V9-SC-006 canonical payload staging or RRRM-1 h5ad regeneration.

Slide angle:

- "v9 is not another analysis branch; it is the benchmark platform layer."
- "The project now has task manifests, metric contracts, provenance, and
  release boundaries."

## Cross-Version Scientific Claims

Safe headline claims:

- Spaceflight transcriptomic generalization is tissue-dependent.
- Thymus and gastrocnemius are more robust cross-mission benchmark tissues than
  liver in the original mouse bulk task family.
- Pathway abstractions can rescue some weak gene-level tasks and resist
  mission/batch confounders.
- Current gene-expression foundation models and text LLMs do not reliably beat
  tuned classical baselines on small-n bulk spaceflight transcriptomics.
- Single-cell and organoid extensions are promising but must be kept separate
  from bulk leaderboard claims until payloads, metrics, and provenance are
  frozen.
- v8 translational outputs are hypothesis-generation and prioritization tools,
  not crew-health countermeasure recommendations.
- v9 is correctly framed as an alpha scaffold / metadata-alpha platform, not a
  frozen benchmark release.

Unsafe or overclaiming language to avoid:

- "Foundation models failed universally in biology." The claim is narrower:
  they fail to beat classical baselines on these small-n spaceflight benchmark
  surfaces.
- "v8 identifies countermeasures." Safer: v8 prioritizes perturbation
  hypotheses for follow-up validation.
- "v9 is released/frozen." Safer: v9 has a metadata-only alpha scaffold with
  explicit payload blockers.
- "Organoid or OSD-120 is a leaderboard task." Safer: they are draft
  diagnostic or extension tracks.

## Likely Slide Sections

1. Title / motivation.
2. Problem: spaceflight omics needs mission-held-out benchmarks.
3. Data and benchmark design.
4. v1 key result: thymus > liver.
5. Pathway rescue and batch resistance.
6. Foundation model and LLM comparison.
7. Temporal and single-cell extensions.
8. Multi-species/spatial negative and positive controls.
9. v4 method hardening.
10. v5/v6 biology and translation.
11. v7/v7.1 public benchmark consolidation.
12. v8 SpaceMed translational incubator.
13. v9 SpaceBio-Bench platform architecture.
14. v9 extension tracks: organoids, multispecies, single-cell.
15. Current status, blockers, and next steps.

## Deck-Building Notes

- Use v1-v7 as the completed scientific benchmark story.
- Use v8 as a separate incubator / translational hypothesis layer.
- Use v9 as the current platformization and release-boundary story.
- Keep benchmark result slides separate from release/provenance slides.
- Put "not claimed" boundaries directly on slides for v8 and v9 so the story
  stays reviewer-safe.

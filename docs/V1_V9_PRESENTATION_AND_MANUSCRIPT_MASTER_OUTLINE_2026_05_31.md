# v1-v9 Presentation And Manuscript Master Outline

Date: 2026-05-31

Purpose: convert the full GeneLab Benchmark / SpaceBio-Bench v1-v9 project into
a coherent presentation and manuscript plan, using the internal result
inventories and the external positioning scan.

Companion documents:

- `docs/PROJECT_SLIDE_CONTENT_INVENTORY_V1_TO_V9_2026_05_31.md`
- `docs/PROJECT_RESULTS_LOCATION_INVENTORY_2026_05_31.md`
- `docs/PROJECT_RESULTS_DEEP_READ_AUDIT_2026_05_31.md`
- `docs/SIMILAR_PROJECTS_AND_POSITIONING_SCAN_2026_05_31.md`
- `docs/LITERATURE_QUERY_LOG_SPACEBIO_BENCH_2026_05_31.md`
- `docs/CANONICAL_RESULTS_V7_1.md`
- `docs/V9_PURPOSE_DRIFT_AUDIT_2026_05_26.md`
- `docs/VISUAL_LAYERED_SCENE_BLUEPRINT_2026_06_01.md`

## Bottom Line

Yes: v1-v9 can be turned into one strong presentation.

For one manuscript, the safest form is not "all results from v1 through v9 as
equal scientific discoveries." The safer and stronger form is:

> SpaceBio-Bench: a mission-held-out, provenance-first benchmark platform for
> evaluating biological AI under GeneLab/OSDR spaceflight omics domain shift.

In that framing:

- v1-v7 are the completed benchmark result core.
- v8 is a translational hypothesis-generation incubator, not a clinical or
  operational countermeasure result.
- v9 is the platform/release/provenance layer, with public bulk alpha plus
  organoid, multispecies, and single-cell draft extension tracks.

The paper can include all v1-v9, but the result hierarchy must be explicit.
The presentation can be broader and more narrative; the manuscript should be
more conservative and evidence-tiered.

## Recommended Product Split

### Product 1: Full v1-v9 Presentation

Audience:

- ASGSR / NASA / space biology audience.
- Computational biology / AI-for-biology audience.
- Lab/internal strategy audience.

Recommended title:

> SpaceBio-Bench: Testing Biological AI Under Spaceflight Domain Shift

Presentation thesis:

> Spaceflight omics needs benchmarks that test generalization across unseen
> missions. GeneLab Benchmark established that task, showed tissue-specific
> generalization and classical-baseline strength, then matured into
> SpaceBio-Bench: a provenance-first platform with translational, organoid,
> multispecies, and single-cell extension lanes.

This presentation can cover the whole arc from v1 to v9 because decks can show
program evolution, negative results, blockers, and future tracks more naturally
than a paper can.

### Product 2: One Manuscript, If Kept Unified

Recommended manuscript type:

- Benchmark/resource/platform paper.
- Better fit: Nature Methods / Genome Biology / Scientific Data + Methods-style
  resource framing.
- Riskier fit: pure discovery biology paper, because v8/v9 contain many
  diagnostic and alpha-stage outputs that should not be promoted as finalized
  biological claims.

Recommended title options:

1. `SpaceBio-Bench: a mission-held-out benchmark for biological AI in spaceflight omics`
2. `A provenance-first benchmark platform for cross-mission generalization in GeneLab spaceflight transcriptomics`
3. `Evaluating biological AI under spaceflight domain shift with mission-held-out GeneLab/OSDR benchmarks`

Manuscript thesis:

> We define a benchmark platform for spaceflight omics generalization, show
> that tissue identity and pathway abstraction strongly shape cross-mission
> transfer, benchmark classical ML against foundation-model and LLM baselines,
> and provide provenance-backed task manifests and extension tracks for future
> organoid, multispecies, and single-cell evaluation.

### Product 3: Cleaner Two-Manuscript Split

If time and venue strategy allow, the cleanest publication split is:

1. v1-v7 benchmark paper:
   - Multi-tissue mission-held-out mouse transcriptomics benchmark.
   - Tissue hierarchy, pathway rescue, confounder/batch resistance, FM/LLM
     comparison, held-out validation.
2. v8-v9 platform/resource paper:
   - SpaceMed translational hypothesis layer plus SpaceBio-Bench v9 platform.
   - Public bulk alpha, manifests, provenance, organoid, multispecies,
     single-cell task contracts.

This split reduces reviewer load. But one unified paper remains possible if
v8/v9 are framed as platform extensions rather than equal-strength discoveries.

## Core Narrative Spine

### Act 1: The Benchmark Problem

Field context:

- NASA OSDR/GeneLab contains rich public space biology omics data.
- Existing resources and studies are strong, but they do not directly solve
  the benchmark question: can a model trained on prior missions generalize to
  an unseen mission?
- External neighbors exist, especially NASA BPS RNA-seq benchmark datasets,
  GLARE, NASA AI4LS method studies, OpenProblems, and cell-eval.
- No one-to-one external duplicate was found for a multi-tissue,
  mission-held-out, provenance-first GeneLab/OSDR benchmark platform.

Safe positioning:

- Do not claim "first AI benchmark for space biology."
- Claim a distinct mission-held-out GeneLab/OSDR benchmark/platform niche.

### Act 2: v1-v7 Completed Benchmark Result Core

Main claim:

> Cross-mission spaceflight transcriptomic generalization is tissue-specific,
> pathway-dependent, and not automatically solved by current foundation models
> or text LLMs.

Evidence:

- Thymus generalizes better than liver in the original v1 transfer hierarchy.
- Category A/LOMO and Category B transfer results show strong tissue
  dependence.
- Pathway features can rescue weak gene-level tasks and reduce mission/batch
  detectability.
- v4 broad method grid shows the result is not a narrow model choice artifact.
- v5/v6 add biological interpretation and partial mouse-human translation.
- v7/v7.1 lock public-safe counts, claims, and model-comparison boundaries.

### Act 3: v8 Translational Incubator

Main claim:

> The benchmark can support translational hypothesis generation, but v8 outputs
> must remain hypothesis-only.

Evidence:

- Mouse NES features improve I4/Twins pathway prediction in v8 BRIDGE.
- Stressor decomposition shows interaction-dominated response structure and
  radiation-quality sensitivity.
- Intervention triage nominates perturbation hypotheses but not validated
  countermeasures.
- Mars-regime outputs are exploratory flags, not point predictions.

### Act 4: v9 Platformization

Main claim:

> The project has moved from a benchmark analysis into a reproducible benchmark
> platform with task manifests, metric profiles, registries, provenance, and
> explicit alpha-stage release boundaries.

Evidence:

- `spacebio_bench` package skeleton and task infrastructure exist.
- Eight public bulk LOMO task manifests exist.
- Public bulk baseline summaries exist across three simple baselines.
- Source inventory and checksum-manifest audit exist for 22 public bulk source
  rows.
- Human organoid, multispecies, and single-cell extension lanes exist, but are
  not primary leaderboard claims.
- Public bulk alpha is metadata-only because local payload-hash verification is
  still a blocker.

## Presentation Outline

Recommended length: 20-24 slides for a full talk; 12-15 slides for conference
short format.

### Full Talk Version

| Slide | Title | Core content | Evidence source |
|---:|---|---|---|
| 1 | SpaceBio-Bench | One-line thesis and ecosystem positioning | this outline |
| 2 | Why Mission-Held-Out? | Random splits are not enough; mission is the independence unit | `docs/V1_PAPER_CONTENT.md`, `evaluation/RESULTS_SUMMARY.md` |
| 3 | External Landscape | NASA BPS, GLARE, AI4LS, OpenProblems, cell-eval; no 1:1 duplicate | `docs/SIMILAR_PROJECTS_AND_POSITIONING_SCAN_2026_05_31.md` |
| 4 | Project Evolution | v1-v9 timeline from benchmark to platform | `docs/PROJECT_SLIDE_CONTENT_INVENTORY_V1_TO_V9_2026_05_31.md` |
| 5 | Data And Task Design | OSDR/GeneLab, tissues, missions, LOMO, feature layers | `docs/CANONICAL_RESULTS_V7_1.md` |
| 6 | Thymus Beats Liver | v1 tissue hierarchy and transfer result | `evaluation/RESULTS_SUMMARY.md` |
| 7 | Pathway Rescue | Eye/kidney rescue and pathway abstraction | `evaluation/RESULTS_SUMMARY.md` |
| 8 | Confounder Resistance | Mission identity/batch detectability: genes vs pathways | `evaluation/RESULTS_SUMMARY.md` |
| 9 | Model Comparison | Classical ML vs scGPT, Geneformer, UCE, scFoundation, text LLMs | `docs/CANONICAL_RESULTS_V7_1.md` |
| 10 | v4 Hardening | 8 tissues x 8 classifiers x 4 feature types | `docs/CANONICAL_RESULTS_V7_1.md` |
| 11 | Temporal And Single-Cell Lessons | preservation, recovery, age, RRRM-1 | `v2/evaluation/V2_RESULTS_SUMMARY.md` |
| 12 | Negative Results Matter | spatial brain negative, Drosophila-mouse directionality, FM negatives | `v3/README.md` |
| 13 | Biological Interpretation | immune, TF, metabolism, drug targets, biomarkers | `v5/evaluation/` and inventory |
| 14 | Human Translation | pathway conservation and target evidence, not clean gene transfer | `v6/evaluation/` and inventory |
| 15 | v7.1 Public Boundary | canonical counts and claim discipline | `docs/CANONICAL_RESULTS_V7_1.md` |
| 16 | v8 SpaceMed | bridge, decompose, intervene, causal pillars | `v8/RESULTS_SUMMARY.md` |
| 17 | v8 Claim Boundary | hypotheses only; no countermeasure recommendation | `docs/PROJECT_SLIDE_CONTENT_INVENTORY_V1_TO_V9_2026_05_31.md` |
| 18 | v9 Platform Architecture | package, manifests, metrics, evaluator, run manifests | `v9/README.md` |
| 19 | Public Bulk Alpha | metadata-only alpha, source inventory, checksum audit, blockers | `docs/V9_PURPOSE_DRIFT_AUDIT_2026_05_26.md`, `v9/README.md` |
| 20 | Organoid Track | OSD-863/871, 42 samples, DE references, diagnostic-only status | `v9/human_organoid/` |
| 21 | Multispecies Track | OSD-207, OSD-37, OSD-120 diagnostic branch | `v9/multispecies/` |
| 22 | Single-Cell Track | RRRM inventory, draft manifest, metric spec, missing payload blocker | `v9/sc_spaceflight/` |
| 23 | What Is Claimed / Not Claimed | safe and unsafe claim table | this outline |
| 24 | Roadmap | payload freeze, figure QA, PubMed/Scholar log, v9 single-cell staging | this outline |

### Short Talk Version

Use 12 slides:

1. Title and thesis.
2. Why mission-held-out benchmarks matter.
3. External landscape and gap.
4. Data/task design.
5. Thymus beats liver.
6. Pathway rescue and confounder resistance.
7. Classical baselines beat current FMs/LLMs.
8. v4-v7 hardening and validation.
9. v8 translational incubator.
10. v9 platform architecture.
11. Organoid/multispecies/single-cell extension lanes.
12. Current status, blockers, and next steps.

## Manuscript Outline

### Abstract

Structure:

1. Problem: spaceflight omics lacks standardized mission-held-out AI
   evaluation.
2. Approach: build GeneLab/OSDR multi-tissue benchmark with LOMO tasks,
   gene/pathway features, model tiers, and provenance-backed manifests.
3. Result 1: tissue-specific hierarchy, with thymus outperforming liver in
   original cross-mission transfer.
4. Result 2: pathway abstraction rescues selected tasks and reduces
   mission/batch confounding.
5. Result 3: current FM/LLM/GNN-style comparators do not overturn classical
   baseline strength on this small-n bulk domain-shift surface.
6. Platform: v9 provides task manifests, metric profiles, public bulk
   metadata-alpha, and organoid/multispecies/single-cell extension lanes with
   explicit claim boundaries.

### Introduction

Paragraph plan:

1. Spaceflight biology has a reproducibility/generalization problem, not just a
   data-availability problem.
2. OSDR/GeneLab, SOMA, SpaceOmicsBench, NASA BPS benchmark resources, GLARE,
   and AI4LS are the key ecosystem context.
3. The missing evaluation layer is mission-held-out generalization across
   tissues and task families.
4. Biological AI evaluation should include classical baselines, pathway
   abstraction, foundation-model comparators, negative controls, and
   provenance.
5. This study introduces SpaceBio-Bench and reports the v1-v9 development arc.

### Results

Recommended main Results sections:

1. `SpaceBio-Bench defines mission-held-out evaluation for GeneLab/OSDR omics`
   - Explain v1-v7 data, tasks, LOMO design, full-release scope, and v9
     manifest/platform layer.

2. `Spaceflight transcriptomic generalization is strongly tissue-dependent`
   - v1 tissue hierarchy.
   - Thymus versus liver.
   - Held-out validation.

3. `Pathway abstraction separates transferable biology from mission artifacts`
   - Gene versus pathway rescue.
   - Confounder/batch resistance.
   - NES conservation / transfer relationship if retained in final claim set.

4. `Simple classical baselines remain strong under small-n spaceflight domain shift`
   - Classical ML, scGPT, Mouse-Geneformer, UCE, scFoundation, text LLMs.
   - Caption all subset differences carefully.
   - Include v7/scPRINT/GNN conclusion only as extension/comparator summary.

5. `Benchmark hardening confirms the pattern across methods, tissues, and controls`
   - v4 256-evaluation grid.
   - v5 mechanism layers.
   - v6 human translation.
   - v7.1 public-safe scope.

6. `Temporal, single-cell, spatial, and multispecies extensions expose both signal and failure modes`
   - v2 preservation/recovery/age.
   - v2 RRRM-1 scRNA-seq.
   - v3 spatial brain negative.
   - v3 RRRM-2 PBMC/NK positive.
   - Keep this section concise to avoid diluting the main paper.

7. `SpaceMed v8 turns benchmark outputs into bounded translational hypotheses`
   - BRIDGE result.
   - Stressor decomposition.
   - Intervention and Mars sensitivity as exploratory/hypothesis-only.
   - Strong no-overclaim language.

8. `v9 turns the benchmark into a provenance-first public platform`
   - Task manifests, metric profiles, registry, evaluator, run manifests.
   - Public bulk metadata-only alpha.
   - Organoid, multispecies, and single-cell tracks as extension scaffolds.
   - Blockers and future release path.

### Discussion

Recommended discussion messages:

1. Mission-held-out evaluation changes the biological story: more data in a
   tissue does not guarantee more transferable signal.
2. Pathway abstractions are useful because they can preserve biology while
   reducing mission/batch artifacts.
3. Current foundation models are not automatically superior in small-n bulk
   spaceflight transcriptomics; benchmark design remains decisive.
4. Negative results are valuable: spatial brain, Drosophila-mouse direction,
   donor/organoid confounding, and missing single-cell payloads define honest
   boundaries.
5. v8/v9 show how benchmark infrastructure can support translation and platform
   release without conflating hypotheses with validated countermeasures.
6. The next phase should prioritize payload verification, visual figure QA,
   single-cell payload staging, and formal query-log export before submission.

## Main Figure Plan

| Figure | Role | Core panels | Primary sources | Status / risk |
|---|---|---|---|---|
| Fig 1 | Benchmark/platform overview | OSDR/GeneLab sources, LOMO design, v1-v9 timeline, manifest/provenance layer | inventories, `docs/CANONICAL_RESULTS_V7_1.md`, `v9/README.md` | Needs new schematic |
| Fig 2 | Tissue hierarchy | Category A/B tissue AUROCs, thymus vs liver, held-out validation | `evaluation/RESULTS_SUMMARY.md`, `docs/CANONICAL_RESULTS_V7_1.md` | Strong |
| Fig 3 | Pathway abstraction | gene vs pathway rescue, mission/batch detection, pathway conservation | `evaluation/RESULTS_SUMMARY.md` | Strong, but verify exact values |
| Fig 4 | Model tier comparison | classical ML vs FMs/LLMs/GNN-style comparators with subset labels | `docs/CANONICAL_RESULTS_V7_1.md`, v7 docs | Strong if captioned carefully |
| Fig 5 | Benchmark hardening and biological interpretation | v4 grid, v5 immune/TF/metabolic/drug layers, v6 translation | v4-v6 summaries | Medium; could become crowded |
| Fig 6 | Extension results and failure modes | v2 temporal, v3 spatial negative, RRRM scRNA positive/negative | v2/v3 summaries | Medium; select only best examples |
| Fig 7 | v8 translational incubator | BRIDGE AUROC improvement, stressor decomposition, hypothesis triage | `v8/RESULTS_SUMMARY.md` | High overclaim risk; label hypothesis-only |
| Fig 8 | v9 platform and extension tracks | task manifests, public bulk alpha, organoid/multispecies/sc tracks, blockers | `v9/README.md`, v9 reports | Strong as platform/status figure |

Minimum paper version:

- Main Figures 1-6.
- Put v8/v9 in Fig 6 or Fig 7 plus extended data.

Full platform paper version:

- Main Figures 1-8.

## Tables

Recommended main tables:

1. `Table 1: Benchmark scope and task families`
   - v1-v7 completed bulk benchmark.
   - v8 translational incubator.
   - v9 platform/extension tracks.

2. `Table 2: External ecosystem and positioning`
   - SOMA and SpaceOmicsBench as self/sibling ecosystem.
   - NASA BPS, GLARE, AI4LS, OSDR/GeneLab, OpenProblems, cell-eval as external
     neighbors.

3. `Table 3: Claim boundary matrix`
   - Completed benchmark claim.
   - Hypothesis-only claim.
   - Metadata-alpha claim.
   - Diagnostic-only extension.
   - Blocked/not claimed.

Supplementary tables:

- v4 full 256-evaluation result table.
- FM/LLM subset coverage table.
- v8 source/provenance table.
- v9 source inventory/checksum audit summary.
- Organoid/multispecies/single-cell extension status tables.

## Claim Boundary Matrix

| Claim | Status | Can use in title/abstract? | Notes |
|---|---|---|---|
| Mission-held-out GeneLab/OSDR benchmark | Core claim | Yes | Main novelty. |
| Thymus generalizes better than liver in original mouse transfer hierarchy | Core result | Yes | Keep scope tied to v1/v7.1 subsets. |
| Pathway features can rescue selected tasks and resist mission artifacts | Core result | Yes | Use exact values only after final value audit. |
| Classical baselines outperform tested FMs/LLMs on these benchmark surfaces | Core result | Yes | Do not claim universal FM failure. |
| v4 shows pattern survives broad method sweep | Core result | Yes | Resolve canonical vs raw best-row wording. |
| v8 suggests countermeasure hypotheses | Hypothesis-only | Maybe, but not as main claim | Never phrase as validated countermeasure. |
| Mars-regime extrapolation predicts risk | Not safe | No | Use exploratory flag language only. |
| v9 public bulk alpha is frozen payload release | Not safe | No | Current status is metadata-only alpha. |
| Organoid task is a leaderboard | Not safe | No | Draft/diagnostic extension only. |
| OSD-120 multispecies result is platform proof | Diagnostic only | No | Useful case study, not core alpha. |
| RRRM single-cell flagship is runnable | Not safe yet | No | Payload missing; manifest/metric/audit scaffold exists. |

## Evidence Map

Use these as source-of-truth layers when drafting:

| Topic | Primary source | Secondary source |
|---|---|---|
| v1-v7 public counts and release boundaries | `docs/CANONICAL_RESULTS_V7_1.md` | `README.md` |
| original v1 paper narrative | `docs/V1_PAPER_CONTENT.md` | `evaluation/RESULTS_SUMMARY.md` |
| root benchmark results | `evaluation/RESULTS_SUMMARY.md` | root `evaluation/*.json` |
| v2 temporal/human/single-cell | `v2/evaluation/V2_RESULTS_SUMMARY.md` | `v2/figures/`, `v2/evaluation/*.json` |
| v3 multispecies/spatial/FM extensions | `v3/README.md` | `v3/evaluation/*.json` |
| v4 broad grid | `docs/CANONICAL_RESULTS_V7_1.md` | `v4/evaluation/M1_summary.json` |
| v5 biological mechanism | `v5/evaluation/` | `v5/figures/` |
| v6 human translation | `v6/evaluation/` | `v6/figures/` |
| v7/scPRINT/GNN | `v7/evaluation/` | `docs/CANONICAL_RESULTS_V7_1.md` |
| v8 SpaceMed | `v8/RESULTS_SUMMARY.md` | `v8/*/evaluation/`, `v8/provenance/` |
| v9 platform | `v9/README.md` | `docs/V9_PURPOSE_DRIFT_AUDIT_2026_05_26.md` |
| external positioning | `docs/SIMILAR_PROJECTS_AND_POSITIONING_SCAN_2026_05_31.md` | `docs/LITERATURE_QUERY_LOG_SPACEBIO_BENCH_2026_05_31.md` |

## Immediate Pre-Deck Work

Before building slides:

1. Render and visually inspect existing HTML/PDF figures.
2. Decide whether to use existing figures or regenerate a unified visual style.
3. Resolve v4 canonical table versus raw `M1_summary.json` best-row differences.
4. Select one short-talk and one long-talk outline.
5. Build a slide asset manifest mapping each slide to a figure/table/source.

## Immediate Pre-Manuscript Work

Before drafting the full paper:

1. Freeze the manuscript type: benchmark paper, platform paper, or split
   two-paper strategy.
2. Repeat the formal literature search with PubMed, Google Scholar, and
   Semantic Scholar export.
3. Make a final numeric claim table with one source-of-truth file per number.
4. Render figures and perform visual QA.
5. Decide how the public Hugging Face/GitHub artifacts are cited under the
   target venue's review policy.
6. Keep v8/v9 status language in a reviewer-safe claim boundary table.

## Recommended Next Step

The next concrete work block should be:

> Build a slide asset manifest and figure audit: for each proposed slide, list
> the exact existing figure file or the new figure that must be generated,
> the source data file, the headline claim, and the claim-boundary note.

This will prevent the deck from becoming a beautiful but unverifiable story,
and it will simultaneously create the figure inventory needed for the paper.

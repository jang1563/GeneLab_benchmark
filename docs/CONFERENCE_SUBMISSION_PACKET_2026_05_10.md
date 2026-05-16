# GeneLab Benchmark v7.1 Conference Submission Packet

Status: working submission packet
Date: 2026-05-10
Primary target: MLCB 2026 8-page paper

This packet translates the v7.1 canonical benchmark surface into a compact
conference-submission plan. It is not a new analysis plan. It is a writing and
packaging plan for the already validated GeneLab Benchmark v1-v7 results.

## Recommendation

Lead with **GeneLab Benchmark v7.1 -> MLCB 2026 8-page paper**.

The strongest conference-facing hook is not "space biology dataset release" by
itself. The hook is:

> Under cross-mission bulk RNA-seq shift, tuned classical baselines remain
> stronger than current gene-expression foundation models, and pathway-level
> biology explains when generalization succeeds.

This is a better MLCB story than a broad journal-style atlas because it centers
the machine-learning claim:

- benchmark design with mission-held-out evaluation;
- small-n bulk domain shift where model choice matters;
- direct comparison of classical ML, gene-expression foundation models, and
  text LLMs;
- biological interpretation of why transfer works in some tissues but not
  others.

Keep v8 SpaceMed intervention, countermeasure, and Mars-extrapolation material
out of this submission. v8 can be mentioned only as future work if needed.

## Venue Fit

| Rank | Venue / Track | Current Public Dates Checked | Fit | Recommendation |
|---:|---|---|---|---|
| 1 | MLCB 2026, 8-page paper | MLCB home page: Nov 16-17, 2026 at New York Genome Center. Submit page: submissions due July 1 AoE, 8-page paper or 2-page abstract tracks. | Excellent | Primary target. Use 8-page paper, not 2-page abstract. |
| 2 | NeurIPS 2026 workshop contribution | NeurIPS workshop guidance: suggested contribution date Aug 29, mandatory notifications Sep 29, workshops Dec 11-13 depending site. | Strong backup | Wait for accepted workshop list, then submit to the best biology / AI-for-science / evaluation workshop. |
| 3 | AI4Science 2026 full paper | Full paper Aug 8, notification Sep 8, conference Oct 23-25. | Good but less proven | Backup only if proceedings and template are acceptable for this project. |
| 4 | Journal follow-up | Flexible | Strong | Use expanded Genome Biology / Nature Methods version after conference packet. |

MLCB note: as of 2026-05-10, the MLCB main page is updated for 2026, while the
submit page still contains some stale year text. Treat July 1 AoE as the active
working deadline, but verify the CMT page before final formatting.

Sources checked:

- https://www.mlcb.org/
- https://www.mlcb.org/submit
- https://neurips.cc/Conferences/2026/WorkshopsGuidance
- https://neurips.cc/Conferences/2026/CallForWorkshops
- https://ic-ai4science.org/

## Submission Positioning

Working title:

> GeneLab Benchmark: Cross-Mission Evaluation of Foundation Models and Classical
> Baselines for Spaceflight Transcriptomics

Alternate shorter title:

> Cross-Mission Benchmarks Reveal Foundation Model Limits in Spaceflight
> Transcriptomics

One-sentence contribution:

> We introduce a mission-held-out benchmark for NASA OSDR mouse spaceflight
> transcriptomics and show that classical baselines outperform current
> gene-expression foundation models under small-n bulk RNA-seq shift, with
> pathway conservation explaining tissue-specific transfer.

What reviewers should remember:

- **Task**: predict flight vs ground across held-out missions.
- **Main result**: PCA-LR is the best 8-tissue gene-level baseline; current
  gene-expression FMs and text LLMs do not beat it on this benchmark surface.
- **Biology**: thymus generalizes better than liver; pathway-level features
  resist mission/batch effects and rescue selected tissues.
- **Validation**: held-out thymus RR-23 and skin RR-7 support the benchmark
  signal beyond the original LOMO folds.

## Canonical Result Surface

Use `docs/CANONICAL_RESULTS_V7_1.md` as the source of truth. Do not copy values
from older paper drafts without checking that file.

Required headline table:

| Tissue | Best AUROC | Best Method | Best Feature | Significance note |
|---|---:|---|---|---|
| Thymus | 0.948 | PCA-LR | KEGG | p<0.05 |
| Colon | 0.921 | PCA-LR | KEGG | p<0.05 |
| Lung | 0.901 | PCA-LR | Gene | p<0.05 |
| Kidney | 0.829 | ElasticNet-LR | Hallmark | p<0.01 |
| Eye | 0.823 | PCA-LR | Hallmark | best-row p not significant |
| Skin | 0.819 | PCA-LR | Gene | best-row p not significant |
| Gastrocnemius | 0.776 | PCA-LR | Gene | best-row p not significant |
| Liver | 0.670 | PCA-LR | Gene | best-row p not significant |

Required caveats:

- v1-v7 full release spans 8 tissues and 24+ OSD accessions.
- Several original LOMO and FM rows use the 6-tissue v1 core.
- Do not present scGPT, Mouse-Geneformer, scFoundation, and UCE as one uniform
  8-tissue FM leaderboard.
- Do not mix v8 translational claims into the v7.1 benchmark paper.

## 8-Page MLCB Outline

### Abstract

Problem, benchmark, result, biological explanation, release.

Suggested abstract claim:

> Spaceflight transcriptomics is an extreme small-n domain-shift setting:
> models must generalize across missions, hardware, duration, and tissue
> contexts. GeneLab Benchmark standardizes mission-held-out evaluation across
> NASA OSDR mouse spaceflight RNA-seq and compares classical baselines,
> gene-expression foundation models, and text LLMs. Across the v7.1 public
> benchmark surface, PCA-LR and ElasticNet-LR remain the strongest baselines,
> while current foundation models fail to consistently outperform tuned
> classical methods. Pathway-level features reduce mission/batch separability
> and explain tissue-specific transfer, with thymus outperforming liver under
> cross-mission generalization. The benchmark provides reproducible tasks,
> feature matrices, and result artifacts for evaluating future biological
> foundation models under real distribution shift.

### 1. Introduction

- Spaceflight biology has many datasets but few cross-mission generalization
  benchmarks.
- Foundation models are often evaluated in same-distribution or cell-level
  settings; bulk RNA-seq small-n shift is different.
- Contributions:
  - mission-held-out benchmark;
  - multi-method comparison over 8 tissues x 8 classifiers x 4 feature types;
  - FM and text LLM comparison;
  - pathway explanation and held-out validation.

### 2. Benchmark Design

- Dataset scope: 8 tissues, 24+ OSD accessions, 600+ binary/control samples
  across release layers.
- Evaluation unit: mission, not random sample.
- Feature surfaces: gene, Hallmark, KEGG, combined pathway.
- Leakage control: feature selection inside fold; held-out validation for thymus
  and skin.

### 3. Results

Recommended result order:

1. Multi-method headline: PCA-LR gene mean AUROC 0.776, ElasticNet-LR second at
   0.762.
2. Thymus and colon lead the best-row tissue table; liver is not the strongest
   generalizer.
3. Pathway features rescue kidney, thymus, and eye contexts and reduce
   confounder separability.
4. Gene-expression FMs and text LLMs underperform the classical benchmark
   surface.
5. Held-out validation supports the strongest tissue claims.

### 4. Discussion

- Why a simple baseline wins here: small-n, bulk RNA-seq, mission shift, and
  domain mismatch with single-cell pretraining.
- Why this is useful: it is a stress test for future biological FMs.
- Limitations: mouse-only core, heterogeneous mission metadata, subset-specific
  FM coverage, small folds.
- Future work: v8 translational extension, but not as evidence for v7.1 claims.

## Figure Plan

For an 8-page paper, keep to 4 main figures:

| Figure | Content | Purpose |
|---|---|---|
| 1 | Benchmark schematic: OSDR -> tissue/task construction -> mission-held-out evaluation | Establish task and leakage controls |
| 2 | v4 multi-method result heatmap plus best-row tissue table | Show classical baselines and tissue hierarchy |
| 3 | Pathway-vs-gene explanation: batch/confounder F1 and rescue examples | Explain why transfer differs by feature layer |
| 4 | FM / text LLM comparison plus held-out validation inset | Land the ML benchmark message |

Move systems biology, v5/v6 translation, v8, and long mechanism detail to
supplement or future journal version.

## Immediate Checklist

1. Confirm MLCB CMT 2026 details and whether PMLR publication is desired.
2. Create a single-column 8-page LaTeX draft using 11 pt and 1 inch margins.
3. Build static figure exports from existing HTML/PNG/PDF artifacts.
4. Freeze all public claims against `docs/CANONICAL_RESULTS_V7_1.md`.
5. Add a short artifact availability paragraph:
   - GitHub for code;
   - Hugging Face for processed benchmark feature matrices;
   - NASA OSDR for source data;
   - future Zenodo DOI if available before submission.
6. Run final stale-value scan for old table values:
   - pre-v7.1 gastrocnemius and liver best-row values from older HF card text;
   - retired classifier examples from older HF card text;
   - unsupported language that frames all FM rows as one uniform 8-tissue
     leaderboard.

## Decision

Proceed with MLCB 2026 as the primary target unless the CMT page contradicts
the public MLCB deadline. The paper should be written as a benchmark and model
evaluation contribution, not as a broad space biology resource paper.

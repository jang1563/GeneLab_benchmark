# GeneLab Benchmark v7.1 Canonical Results and Scope

Status: canonical v7.1 result source
Dataset freeze: 2026-03-01
Patch date: 2026-05-09
Current public documentation patch: v7.1.2 (2026-06-05)

This file is the canonical public-facing source for v7.1 scope accounting,
headline result tables, and release-boundary notes. README, Hugging Face card,
paper outline, and citation metadata should be updated from this file in the
same edit batch when any public count or result value changes.

## Release Boundary

Use the following boundaries in public text:

- **v1-v7 GeneLab Benchmark**: cross-mission transcriptomics benchmark,
  multi-method evaluation, foundation-model comparisons, and result synthesis.
- **v7.1 documentation patch**: public metadata and result-surface consistency;
  no new benchmark result generation.
- **v7.1.2 public-card/metadata patch**: reviewer-facing card pack,
  release metadata, author metadata, and access documentation for the canonical
  v7.1 result surface; no new benchmark result generation.
- **v8 SpaceMed**: incubating translational extension. Do not mix v8
  intervention, countermeasure, or Mars-extrapolation claims into the v7.1
  benchmark paper.

## Scope Accounting

Keep full-release scope separate from analysis subsets:

| Surface | Scope |
|---|---|
| Full v1-v7 benchmark release | 8 tissues: liver, gastrocnemius, kidney, thymus, skin, eye, lung, colon |
| Public OSDR source catalog | 24+ OSD accessions, with excluded and track-specific rows documented in `DATA_CATALOG.md` |
| Processed sample scope | 600+ binary/control samples across the release layers |
| v4 multi-method evaluation | 8 tissues x 8 classifiers x 4 feature types = 256 evaluations |
| v1 LOMO/FM core | 6 tissues, used for several original LOMO and FM comparison rows |
| HF feature-matrix package | Reviewer-facing feature matrices for selected LOMO tasks plus versioned result artifacts |

Do not write "8-tissue benchmark" and then report a 6-tissue FM mean without a
subset note. The FM table below intentionally labels subset coverage.

## Canonical v4 Multi-Method Headline Table

Source of truth: `evaluation/RESULTS_SUMMARY.md`, mirrored in the root README.

| Tissue | Best AUROC | Best Method | Best Feature | Significance note |
|---|---:|---|---|---|
| Thymus | 0.948 | PCA-LR | KEGG | p<0.05 |
| Colon | 0.921 | PCA-LR | KEGG | p<0.05 |
| Lung | 0.901 | PCA-LR | Gene | p<0.05 |
| Kidney | 0.829 | ElasticNet-LR | Hallmark | p<0.01 |
| Eye | 0.823 | PCA-LR | Hallmark | best-row p not significant in README/evaluation summary |
| Skin | 0.819 | PCA-LR | Gene | best-row p not significant in README/evaluation summary |
| Gastrocnemius | 0.776 | PCA-LR | Gene | best-row p not significant in README/evaluation summary |
| Liver | 0.670 | PCA-LR | Gene | best-row p not significant in README/evaluation summary |

Additional v4 note:

- PCA-LR has the best 8-tissue gene-level mean AUROC: 0.776.
- ElasticNet-LR is second: 0.762.
- Across all 256 configurations, 40 evaluations are significant at p<0.05 and
  6/8 tissues have at least one significant configuration. This is not the same
  as saying every row in the best-method headline table has a significant
  best-row permutation p-value.

## Canonical FM / LLM Snapshot

Source of truth: root README plus `evaluation/RESULTS_SUMMARY.md`.

| Model | Reported surface | AUROC / score | Comparison |
|---|---|---:|---|
| PCA-LR baseline | 6-tissue v1 mean | 0.758 | classical reference |
| scGPT | 6-tissue v1 mean | 0.666 | below PCA-LR |
| scFoundation | best single-tissue v3 row | 0.635 | below PCA-LR surface |
| UCE | best single-tissue v3 row | 0.632 | below PCA-LR surface |
| Mouse-Geneformer | 6-tissue v1 mean | 0.476 | below PCA-LR |
| Text LLMs | zero-shot text reasoning | 0.47-0.51 | chance-level |

Required caption language:

> scGPT and Mouse-Geneformer report 6-tissue v1 means; scFoundation and UCE
> rows show best single-tissue v3 results, with full cross-tissue means below
> the classical baseline. The comparison is a benchmark-surface summary, not a
> single uniform 8-tissue FM leaderboard.

## Held-Out Validation

| Tissue | Mission | AUROC | 95% CI | p-value | n_test |
|---|---|---:|---|---:|---:|
| Thymus | RR-23 | 0.905 | [0.672, 1.000] | 0.005 | 18 |
| Skin | RR-7 | 0.885 | [0.745, 0.986] | <0.001 | 30 |

## License and Citation Notes

Use separate license statements:

- Code: MIT.
- Processed benchmark data / feature matrices: use the license declared on the
  public dataset card.
- Source data: NASA OSDR public data; cite OSDR and follow individual dataset
  licenses where applicable.

If a target venue is double blind, do not rely on named public GitHub/Hugging
Face URLs in the anonymized review PDF unless the venue policy permits it.

## Submission-Safe Claim Set

Safe headline:

> GeneLab Benchmark evaluates cross-mission generalization of mouse spaceflight
> transcriptomic signatures and shows that current gene-expression foundation
> models do not automatically outperform tuned classical baselines under
> small-n bulk RNA-seq shift.

Avoid in the v7.1 benchmark paper:

- v8 SpaceMed intervention or countermeasure claims;
- operational astronaut-health recommendations;
- Mars-risk extrapolation as benchmark evidence;
- unsupported claims that all FM rows are a single uniform 8-tissue comparison.

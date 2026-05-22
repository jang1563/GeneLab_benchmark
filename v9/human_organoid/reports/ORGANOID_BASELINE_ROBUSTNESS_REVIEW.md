# Human Organoid Baseline Robustness Review

Status: draft robustness checkpoint
Date: 2026-05-21
Task: `draft_human_organoid_spaceflight`

## Scope

This note reviews the first draft human organoid nearest-centroid pilot
baseline. It is not release language and is not a leaderboard claim.

The task uses OSD-863/GLDS-716 and OSD-871/GLDS-720 normalized count matrices
aligned to `v9/human_organoid/sample_factors.draft.csv`. The loader currently
sees 42 samples and 27,986 common human gene features across the two processed
matrices.

## External evidence checked

- [GSE259421 via OmicsDI](https://www.omicsdi.org/dataset/geo/GSE259421):
  confirms the study used human iPSC-derived neural organoids, four individuals,
  cortical and dopaminergic fates, with/without isogenic microglia, ISS versus
  ground culture, and transcriptome analysis.
- [GEO sample GSM8115910](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8115910):
  confirms RNA-seq processing details including STAR alignment, mouse-sequence
  contamination filtering via XenofilteR, and featureCounts-derived count
  matrices.
- [Marotta et al. 2024 PubMed record](https://pubmed.ncbi.nlm.nih.gov/39441987/):
  confirms the peer-reviewed organoid study context and publication metadata.

## Sensitivity grid

Command:

```bash
python scripts/run_v9_human_organoid_sensitivity.py
```

Output:

- `v9/human_organoid/reports/sensitivity/human_organoid_baseline_summary.csv`
- `v9/human_organoid/reports/sensitivity/human_organoid_baseline_summary.json`
- Per-variant `predictions.csv`, `metrics.json`, and `run_manifest.json` under
  `v9/human_organoid/reports/sensitivity/<variant_id>/`.

Grid:

- `transform`: `log1p`, `none`
- `scaling`: `zscore`, `none`
- `top_variable_genes`: 100, 500, 2000, 5000, 27986
- Total variants: 20

Metric ranges across the 20 draft variants:

| metric | min | median | max |
|---|---:|---:|---:|
| balanced_accuracy | 0.5159090909 | 0.6789772728 | 0.7693181818 |
| AUROC | 0.5676136364 | 0.6622159091 | 0.9295454545 |
| calibration_error | 0.0242551690 | 0.1521081357 | 0.2445559995 |
| macro_f1 | 0.5062365591 | 0.6610509184 | 0.7584818861 |

Top balanced-accuracy variants:

| variant_id | balanced_accuracy | AUROC | calibration_error |
|---|---:|---:|---:|
| `tvg2000_none_zscore` | 0.7693181818 | 0.9295454545 | 0.2445559995 |
| `tvg500_none_zscore` | 0.7465909091 | 0.9204545455 | 0.2184822121 |
| `tvg100_none_zscore` | 0.7227272727 | 0.8136363636 | 0.1806039092 |
| `tvg5000_none_zscore` | 0.6988636364 | 0.8806818182 | 0.1761045148 |
| `tvg27986_log1p_zscore` | 0.6886363636 | 0.8789772727 | 0.1661161537 |

## Default configuration decision

Keep the default as:

```text
transform=log1p
scaling=zscore
top_variable_genes=2000
```

Default metrics:

| metric | value |
|---|---:|
| n_predictions | 84 |
| balanced_accuracy | 0.5295454545 |
| AUROC | 0.6147727273 |
| calibration_error | 0.02537185971 |
| macro_f1 | 0.5194508009 |

Rationale:

- The raw normalized-count z-score variants score higher, but calibration error
  is much worse. That makes them useful sensitivity probes, not a safe default.
- The GEO processing note records mouse-contamination filtering and count-matrix
  generation. Until we audit library-size behavior, gene filtering, and donor
  metadata more deeply, raw-scale advantages should be treated as fragile.
- The log1p default is conservative and keeps the pilot closer to common
  expression-model preprocessing.

## Fold behavior

Default fold-level balanced accuracy:

| fold_id | n_test | balanced_accuracy |
|---|---:|---:|
| `holdout_organoid_type_cortical_neural_organoid` | 19 | 0.644 |
| `holdout_organoid_type_dopaminergic_neural_organoid` | 23 | 0.580 |
| `holdout_microglia_condition_with_microglia` | 20 | 0.449 |
| `holdout_microglia_condition_without_microglia` | 22 | 0.455 |

Best sensitivity variant, `tvg2000_none_zscore`, fold-level balanced accuracy:

| fold_id | n_test | balanced_accuracy |
|---|---:|---:|
| `holdout_organoid_type_cortical_neural_organoid` | 19 | 0.744 |
| `holdout_organoid_type_dopaminergic_neural_organoid` | 23 | 0.625 |
| `holdout_microglia_condition_with_microglia` | 20 | 0.798 |
| `holdout_microglia_condition_without_microglia` | 22 | 0.909 |

Interpretation:

- The default baseline struggles on microglia-condition holdouts.
- The high raw-count variant improves those folds, but the jump is large enough
  to require confounding and scale checks before interpretation.

## Confounding and leakage review

Observed factor structure in `sample_factors.draft.csv`:

- `source_id` and `organoid_type` are perfectly coupled:
  OSD-863 is cortical only; OSD-871 is dopaminergic only.
- Disease context is also coupled with organoid type:
  PPMS appears only in cortical organoids, and sporadic Parkinson disease
  appears only in dopaminergic organoids.
- Ground/LEO labels are reasonably balanced within organoid type, microglia
  condition, and disease context, but the task is still small.
- All samples are from one mission context, so mission-held-out generalization
  cannot be claimed.
- Donor/Subject and iPSC-line metadata are now parsed from GEO GSE259421 series
  metadata for all 42 samples, but donor is coupled to source, organoid fate,
  and disease context. Donor holdouts therefore remain diagnostic-only.

Leakage assessment:

- Train-fold feature selection and scaling are fit only on training samples.
- Test samples enter feature selection only through shared gene identifiers, not
  through variance estimates.
- The current task still risks source/organoid/disease confounding because the
  two source accessions encode distinct biological systems.

## Prohibited claims

Do not claim:

- a frozen benchmark release;
- a stable leaderboard;
- mission-held-out performance;
- donor-generalization;
- human brain or astronaut outcome prediction;
- disease-specific inference for PPMS or Parkinson disease;
- superiority of one preprocessing variant.

Allowed language:

- "draft pilot baseline";
- "sample-factor-backed, matrix-aligned organoid task";
- "sensitivity analysis shows preprocessing-dependent behavior";
- "additional donor, library-size, and DE-reference audits are required."

## Next action

Before manuscript-facing interpretation:

1. Decide whether GEO-derived donor holdouts should remain diagnostic-only or
   become a separate conservatively worded task family.
2. Add fold-level metric summaries to the sensitivity report.
3. Generate a frozen compact DE reference table and signature-output contract
   before computing `de_direction_match` or `signature_rank_correlation`.

## Follow-up diagnostic update

`v9/human_organoid/reports/ORGANOID_DONOR_LIBRARY_SIZE_AUDIT.md` now records the
local follow-up. OSDR SampleTable files expose only `condition`, but GEO
GSE259421 series metadata recovers Subject/iPSC-line ids for 42/42 samples. The
sample-scale diagnostics also show label-associated sample-total and
zero-fraction differences, plus donor-level scale variation coupled to source
and disease context, reinforcing the decision not to promote raw-count
sensitivity variants to the default.

`v9/human_organoid/reports/ORGANOID_DONOR_AWARE_SPLIT_DECISION.md` records the
follow-up split decision. GEO-derived donor holdouts have been evaluated under
`v9/human_organoid/reports/donor_diagnostics/`, but remain diagnostic-only with
`claim_boundary=donor_diagnostic_only_not_leaderboard`.

`v9/human_organoid/reports/ORGANOID_SIGNATURE_METRIC_REFERENCE_DECISION.md`
records the DE/signature metric reference decision. OSDR lists public
differential-expression tables and contrast-definition CSVs for both OSD-863 and
OSD-871. Each source has 56 parsed contrast pairs, including four direct
matched Ground Control versus Space Flight contrasts and four reversed matches.
The draft metric profile can therefore keep `de_direction_match` and
`signature_rank_correlation` as reference-backed future metrics, but current
baseline outputs remain classifier-primary until a frozen contrast subset,
log2 fold-change orientation, and gene-level response-signature submission
format are defined.

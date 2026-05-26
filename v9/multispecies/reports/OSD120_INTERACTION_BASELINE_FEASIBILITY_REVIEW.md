# OSD-120 Interaction Baseline Feasibility Review

Status: draft-only V9-MULTI-011 checkpoint  
Date: 2026-05-23

## Scope

This review covers the first transparent nearest-centroid feasibility baseline
for the OSD-120 Arabidopsis root interaction task. It keeps OSD-120 separate
from the OSD-37 and OSD-207 species-native reports because OSD-120 is a
single-source genotype/ecotype by light-treatment interaction task, not a
simple species-native spaceflight-vs-ground benchmark.

External context supports this split. NASA OSDR describes GeneLab/OSDR as a
standardized repository for space-relevant multi-omics data across organisms,
including plants. The 2024 CARA Arabidopsis transcriptomics report treats light
condition and genotype/ecotype as central axes of spaceflight transcriptional
response, which matches the task design here:

- NASA OSDR overview: https://science.nasa.gov/biological-physical/data/osdr/
- Zhou, Ferl, and Paul 2024, npj Microgravity:
  https://www.nature.com/articles/s41526-024-00417-0

## Inputs

- Task manifest:
  `v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json`
- Sample factors:
  `v9/multispecies/sample_factors.draft.csv`
- Expression matrix:
  `data/multispecies/arabidopsis/GLDS-120_rna_seq_Normalized_Counts_GLbulkRNAseq.csv`
- Baseline summary:
  `v9/multispecies/reports/interaction_nearest_centroid/multispecies_baseline_summary.csv`

The loader aligns 36 parsed OSD-120 samples to 24,740 species-local feature rows.
All folds are sample-count backed and label-balanced:

| Fold family | Folds | Train/test shape | Role |
|---|---:|---|---|
| `primary_genotype_or_ecotype_holdout` | 3 | 24 train / 12 test | Primary candidate |
| `secondary_light_treatment_holdout` | 2 | 18 train / 18 test | Secondary candidate |
| `diagnostic_condition_stratum_holdout` | 6 | 30 train / 6 test | Diagnostic |

## Baseline

The baseline is the same transparent nearest-centroid scaffold used for the
species-native multispecies drafts:

- transform: `log1p`
- scaling: train-fold `zscore`
- feature selection: top 2,000 train-fold variable genes
- classifier: class centroids in transformed/scaled gene space
- output: `predictions.csv`, `metrics.json`, and `run_manifest.json` per fold
  family
- claim boundary: `draft_interaction_feasibility_only_not_leaderboard`

## Summary

| Fold family | Balanced accuracy | AUROC | Calibration error | Delta metric |
|---|---:|---:|---:|---:|
| Primary genotype/ecotype | 0.6667 | 0.7346 | 0.1408 | genotype delta 0.2500 |
| Secondary light treatment | 0.6667 | 0.7840 | 0.1364 | light delta 0.2222 |
| Diagnostic condition stratum | 0.6667 | 0.7654 | 0.1415 | condition delta 0.3333 |

The global balanced accuracy is similar across fold families, but the fold-level
details show meaningful heterogeneity:

- Primary genotype/ecotype holdout:
  - `Col.0`: balanced accuracy 0.75
  - `Col.0.PhyD`: balanced accuracy 0.50
  - `Wassilewskija.ecotype`: balanced accuracy 0.75
- Secondary light-treatment holdout:
  - `Dark.Treatment`: balanced accuracy 0.7778
  - `Light.Treatment`: balanced accuracy 0.5556
- Diagnostic condition-stratum holdout:
  - weakest strata: `Col.0.PhyD|Dark.Treatment` and
    `Wassilewskija.ecotype|Dark.Treatment`, each balanced accuracy 0.50
  - strongest strata: `Col.0|Dark.Treatment` and
    `Wassilewskija.ecotype|Light.Treatment`, each balanced accuracy 0.8333

## Interpretation

The result is useful as a feasibility check, not as a benchmark claim. The
nearest-centroid baseline can run end to end from the manifest, and the three
fold families produce valid prediction/evaluation reports. However, OSD-120 is
small, single-source, and explicitly structured by genotype/ecotype and light.
The fold deltas show that the apparent aggregate score hides condition-specific
fragility.

The primary split should remain genotype/ecotype holdout because it asks the
cleanest transfer question for this dataset: can a model distinguish
spaceflight from ground when one Arabidopsis genotype/ecotype context is unseen?
The light-treatment split is scientifically important but harsher and more
entangled with the CARA design. The condition-stratum split is valuable as a
diagnostic stress test, but each test fold has only six samples, so it should
not be promoted as a leaderboard split.

## Decision

Keep OSD-120 in the `multispecies_light_interaction_spaceflight` family with
three separated evaluation surfaces:

1. Primary: genotype/ecotype holdout.
2. Secondary: light-treatment holdout.
3. Diagnostic: genotype/ecotype by light-treatment condition-stratum holdout.

Do not merge OSD-120 into the OSD-37/OSD-207 species-native summary, and do not
use the current baseline report as release-grade evidence. The next useful work
is a small OSD-120 sensitivity grid to check whether the observed fragile folds
are stable under transform/scaling/feature-count choices.

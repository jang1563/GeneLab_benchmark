# OSD-120 Interaction Task Design

Status: V9-MULTI-009 design note  
Date: 2026-05-23

## Decision

OSD-120 should become a separate draft interaction task, not a third
species-native task in the existing OSD-37/OSD-207 manifest set.

Rationale:

- OSD-120 has a balanced genotype/ecotype by light-treatment design.
- Its scientific question is not just Arabidopsis LEO/Ground classification.
- It should test whether a model handles spaceflight response across both
  genotype/ecotype context and light treatment.
- Mixing it into the simpler OSD-37 species-native task would blur the
  benchmark claim.

Recommended next implementation: create a separate task family and manifest
block, tentatively named `multispecies_light_interaction_spaceflight`, instead
of adding OSD-120 to `v9/multispecies/task_manifests/` immediately.

## Evidence Base

Inputs already available:

- `v9/multispecies/source_inventory.draft.csv`
- `v9/multispecies/sample_factors.draft.csv`
- `v9/multispecies/expression_matrix_audit.draft.csv`
- `v9/multispecies/source_checksum_audit.draft.csv`
- `data/multispecies/arabidopsis/GLDS-120_rna_seq_Normalized_Counts_GLbulkRNAseq.csv`
- `data/multispecies/arabidopsis/GLDS-120_rna_seq_SampleTable_GLbulkRNAseq.csv`

OSD-120 local audit status:

- 36 parsed samples.
- 24,740 feature rows.
- 36/36 normalized-count matrix columns aligned to sample factors.
- OSDR API OK.
- 1 processed checksum manifest.
- 533 parsed MD5 entries.
- Local SampleTable and normalized-count matrix MD5s match OSDR processed
  manifest entries.

## Factor Structure

OSD-120 has a complete 3 x 2 x 2 x 3 structure:

- genotype/ecotype: `Col.0`, `Col.0.PhyD`, `Wassilewskija.ecotype`
- light treatment: `Dark.Treatment`, `Light.Treatment`
- label: `Ground`, `LEO_or_ISS`
- replicates: 3 per genotype/ecotype x light x label cell

Condition-stratum counts:

| condition_stratum | Ground | LEO_or_ISS | total |
|---|---:|---:|---:|
| Col.0\|Dark.Treatment | 3 | 3 | 6 |
| Col.0\|Light.Treatment | 3 | 3 | 6 |
| Col.0.PhyD\|Dark.Treatment | 3 | 3 | 6 |
| Col.0.PhyD\|Light.Treatment | 3 | 3 | 6 |
| Wassilewskija.ecotype\|Dark.Treatment | 3 | 3 | 6 |
| Wassilewskija.ecotype\|Light.Treatment | 3 | 3 | 6 |

## Candidate Splits

### Genotype/Ecotype Holdout

Three folds:

- hold out `Col.0`: train 24, test 12
- hold out `Col.0.PhyD`: train 24, test 12
- hold out `Wassilewskija.ecotype`: train 24, test 12

Each fold has balanced train/test labels. This split asks whether the model can
transfer spaceflight discrimination across genotype/ecotype context while
seeing both light treatments in training.

### Light-Treatment Holdout

Two folds:

- hold out `Dark.Treatment`: train 18, test 18
- hold out `Light.Treatment`: train 18, test 18

Each fold has balanced train/test labels. This split is stricter in one sense
because the model must transfer across lighting context, but it has only two
folds.

### Genotype/Ecotype x Light Condition-Stratum Holdout

Six folds:

- hold out one genotype/ecotype-light condition at a time.
- train 30, test 6 per fold.
- each test fold has 3 Ground and 3 LEO/ISS samples.

This split is the closest analog to the OSD-37/OSD-207 condition-stratum
robustness folds, but each test fold is small. It is useful as a secondary
diagnostic, not as the sole primary split.

## Recommended Manifest Contract

Create one OSD-120 interaction manifest in a separate directory, for example:

- `v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json`

Recommended manifest fields:

- `task_family`: `multispecies_light_interaction_spaceflight`
- `task_id`: `draft_osd120_arabidopsis_root_light_interaction_spaceflight`
- `organism`: `Arabidopsis thaliana`
- `biospecimen_type`: `root`
- `feature_namespace`: `species_gene_ids_pending`
- `release_status`: `draft_not_frozen`

Recommended split structure:

- primary candidate folds: genotype/ecotype holdout
- secondary candidate folds: light-treatment holdout
- tertiary diagnostic folds: genotype/ecotype x light condition-stratum holdout

Recommended metrics:

- `balanced_accuracy`
- `auroc`
- `calibration_error`
- `genotype_holdout_delta`
- `light_treatment_holdout_delta`
- `condition_stratum_holdout_delta`

The interaction-specific delta metrics should be treated as robustness
diagnostics, not standalone leaderboard metrics.

## Boundary

Do not merge OSD-120 with OSD-37 yet, even though both are Arabidopsis. OSD-37
is a cleaner species-native plant feasibility task. OSD-120 is an interaction
design task with a different scientific claim.

Do not claim cross-species generalization from OSD-120. It remains
species-native Arabidopsis until a pathway/NES bridge task is explicitly
designed.

## Next Action

Implement the OSD-120 interaction manifest in a separate manifest directory,
then decide whether the existing multispecies loader/baseline runner should
support interaction tasks directly or whether interaction tasks should receive a
separate loader and baseline runner.

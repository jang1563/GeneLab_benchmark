# Multispecies Species-Native Task Manifest Design

Status: V9-MULTI-006 design note  
Date: 2026-05-23

## Decision

Create two draft species-native task manifests:

- `draft_osd37_arabidopsis_seedling_spaceflight`
- `draft_osd207_drosophila_whole_body_spaceflight`

Do not create an OSD-120 task manifest in this block. OSD-120 has adequate
source, sample, matrix, and local payload checksum evidence, but its scientific
shape is an Arabidopsis genotype/ecotype by light-treatment interaction task.
It should not be merged into the first OSD-37 plant task or treated as a simple
replicate of the same split design.

## Inputs Used

- `v9/multispecies/source_inventory.draft.csv`
- `v9/multispecies/sample_factors.draft.csv`
- `v9/multispecies/expression_matrix_audit.draft.csv`
- `v9/multispecies/source_checksum_audit.draft.csv`
- `v9/multispecies/reports/MULTISPECIES_CHECKSUM_AND_LOCAL_PAYLOAD_AUDIT.md`

These inputs establish that OSD-37 and OSD-207 have parsed condition factors,
aligned normalized-count matrices, OSDR checksum-manifest evidence, and local
SampleTable/matrix MD5 matches for the files used by the scaffold.

## Split Contract

The first multispecies manifests are within-source species-native pilots.

They are not:

- leave-one-mission-out tasks;
- raw-gene cross-species tasks;
- ortholog/pathway bridge tasks;
- frozen public leaderboard tasks.

The manifest split contract is therefore conservative:

- binary target: `LEO_or_ISS` versus `Ground`;
- unit: sample;
- feature namespace: species-native processed gene matrix;
- candidate robustness folds: hold out one `condition_stratum` at a time;
- interpretation: condition-stratum holdout tests robustness across
  genotype/ecotype context, not mission-held-out generalization.

## OSD-37 Shape

OSD-37 has 56 samples, 28 Ground and 28 LEO/ISS, across four Arabidopsis
ecotypes:

- Col.0: 16 samples
- Cvi.0: 12 samples
- Ler.0: 12 samples
- Ws.2: 16 samples

Each condition-stratum candidate fold preserves both labels in train and test.

## OSD-207 Shape

OSD-207 has 32 samples, 16 Ground and 16 LEO/ISS, across four Drosophila
background/variant strata:

- Canton.S_Sei.ts1: 8 samples
- Canton.S_Wild.Type: 8 samples
- w1118_KCNQ370: 8 samples
- w1118_Wild.Type: 8 samples

Each condition-stratum candidate fold has 24 train and 8 test samples with
balanced train/test labels.

## Generated Artifacts

- `v9/multispecies/task_manifests/draft_osd37_arabidopsis_seedling_spaceflight.json`
- `v9/multispecies/task_manifests/draft_osd207_drosophila_whole_body_spaceflight.json`
- `v9/multispecies/task_manifest_index.draft.csv`
- `v9/multispecies/task_manifest_index.draft.json`

## Next Design Block

The next scientific step should be a baseline feasibility design for these two
species-native tasks. Keep it simple:

- read the species-native normalized-count matrices already audited;
- align to the manifest sample ids;
- test conservative train-only preprocessing;
- report baseline outputs as draft diagnostics, not leaderboard claims.

OSD-120 should receive a separate interaction-task design note before any
manifest is generated for it.

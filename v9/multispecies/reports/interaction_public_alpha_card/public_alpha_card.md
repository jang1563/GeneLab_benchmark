---
pretty_name: OSD-120 Arabidopsis Root Light Interaction Diagnostic Draft
license: other
task_categories:
- tabular-classification
tags:
- biology
- space-biology
- genomics
- OSDR
- GeneLab
- Arabidopsis
- draft
---

# OSD-120 Arabidopsis Root Light Interaction Diagnostic Draft

Release status: draft diagnostic alpha card, not a frozen benchmark release.

This card describes the current OSD-120/GLDS-120 diagnostic evidence package in
SpaceBio-Bench v9. It is source-specific and within-source: it evaluates an
Arabidopsis thaliana root bulk RNA-seq light/genotype interaction task derived
from public NASA OSDR/GeneLab processed data.

## Source And Scope

- Source: OSD-120 / GLDS-120
- Source URL: https://osdr.nasa.gov/bio/repo/data/studies/OSD-120
- Organism: Arabidopsis thaliana
- Biospecimen: root
- Assay modality: bulk_rna_seq
- Task id: draft_osd120_arabidopsis_root_light_interaction_spaceflight
- Candidate: sparse_l1_c1
- Candidate variant: tvg2000_log1p_zscore_l1_c1

The external OSD-120 literature and source metadata support the task framing:
Arabidopsis root-tip spaceflight response with genotype/ecotype and light/dark
context. They do not support benchmark-performance claims by themselves.

## Payload Freeze Boundary

| Payload Scope | Count | Status |
| --- | --- | --- |
| Diagnostic-required processed payloads | 2 | True |
| Required payload MD5 matches | 2 | matched |
| Required payload missing | 0 | none |
| OSDR processed payloads outside current diagnostic scope | 531 | not downloaded for this card |
| Full OSDR processed payload freeze | False | not claimed |

The card's input freeze is intentionally narrow. The diagnostic package uses
the OSDR processed SampleTable and normalized-count matrix. Both local files
match the OSDR processed MD5 manifest. The broader OSDR processed payload set is
listed in the freeze manifest but is outside this diagnostic card's required
payload scope.

## Diagnostic Result Surface

| Focus Fold | Nearest BA | Sparse L1 BA | Delta | Stable >=0.5 |
| --- | --- | --- | --- | --- |
| Condition stratum: Col.0.PhyD \| Dark | 0.500 | 0.667 | +0.167 | 10 |
| Light treatment: Light | 0.556 | 0.833 | +0.278 | 2 |
| Genotype/ecotype: Col.0.PhyD | 0.500 | 0.917 | +0.417 | 7 |

Summary metrics:

- Primary genotype/ecotype balanced accuracy:
  0.9166666667
- Secondary light-treatment balanced accuracy:
  0.8333333333
- Diagnostic condition-stratum balanced accuracy:
  0.8888888889
- Nearest-centroid fold comparison:
  9 improve /
  2 tie /
  0 worse
- Stable sparse-model feature rows at selection frequency >=0.5:
  19

## Allowed Claims

- sparse_l1_c1 is a draft transparent diagnostic candidate, not a frozen leaderboard baseline.
- The candidate improves 9 of 11 nearest-centroid fold rows, ties 2, and worsens 0.
- The candidate recovers all three predefined fragile focus folds relative to nearest centroid.
- Nineteen candidate features are selected in at least 50% of balanced train-fold subsamples across the three focus folds.
- OSD-120 is an Arabidopsis CARA root-tip RNA-seq study where genotype/ecotype and lighting context are core design factors.
- The figure/table package is a draft within-source diagnostic surface with explicit disallowed claims.
- sparse_l1_c0p3 remains an appendix/supplement comparator rather than a second primary figure panel.

## Disallowed Claims

- Do not call this a frozen v9 benchmark baseline.
- Do not claim leave-one-mission-out or cross-study generalization.
- Do not describe selected genes as validated biomarkers.
- Do not make operational plant-growth recommendations.
- Do not claim a complete local OSD-120 OSDR payload mirror.

## Files To Inspect

- `interaction_diagnostic_artifact_manifest/diagnostic_artifact_manifest.csv`
- `interaction_diagnostic_artifact_manifest/diagnostic_claim_artifact_map.csv`
- `interaction_payload_freeze_manifest/payload_freeze_manifest.csv`
- `interaction_release_readiness_gap_audit/release_readiness_gap_table.csv`
- `interaction_figure_table_package/figure_main_focus_table.csv`
- `interaction_figure_table_package/figure_stable_feature_appendix.csv`
- `interaction_rebuild_gate/rebuild_gate_summary.csv`
- `interaction_rebuild_gate/rebuild_gate_steps.csv`
- `interaction_rebuild_gate/rebuild_gate_environment.csv`
- `interaction_public_metadata_package/public_metadata_summary.csv`
- `interaction_public_metadata_package/source_release_target_decision.csv`
- `interaction_public_metadata_package/public_metadata_skeleton.json`
- `interaction_ro_crate_citation_scaffold/ro_crate_export_summary.csv`
- `interaction_ro_crate_citation_scaffold/citation_freeze_checklist.csv`
- `interaction_ro_crate_citation_scaffold/ro-crate-metadata.draft.json`

## External Context Links

- https://doi.org/10.1371/journal.pone.0180186
- https://osdr.nasa.gov/bio/repo/data/studies/OSD-120
- https://www.nature.com/articles/s41526-024-00417-0

## Remaining Release Work

- Source release-target promotion is still pending.
- Full OSDR processed payload mirror is not claimed.
- A broader machine-readable metadata skeleton and draft RO-Crate/Data Package
  scaffold are available; archive identifier, creator, and license decisions
  are still pending.
- Packaging rebuild preflight is available through
  `python3 scripts/rebuild_v9_osd120_diagnostic_package.py --repo-root . --mode preflight`.
  It hashes packaging outputs and records environment context, but does not
  rerun model grids or freeze the benchmark release.

## Citation And Credit

Credit NASA OSDR/GeneLab and cite the upstream OSD-120 source when using this
diagnostic package. This local card is a draft SpaceBio-Bench diagnostic
surface and should not replace the upstream OSDR study metadata or citation.

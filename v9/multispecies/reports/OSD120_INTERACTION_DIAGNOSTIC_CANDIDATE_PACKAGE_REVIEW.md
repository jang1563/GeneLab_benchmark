# OSD-120 Diagnostic Candidate Package Review

Status: V9-MULTI-022 complete as of 2026-05-26.

This package turns the V9-MULTI-021 ladder decision into a compact,
figure-ready, claim-boundary-ready evidence bundle for the draft OSD-120
interaction diagnostic. It does not add a new model family and does not promote
the result as a frozen benchmark or leaderboard baseline.

Primary outputs:

- `interaction_diagnostic_candidate_package/candidate_package_summary.csv`
- `interaction_diagnostic_candidate_package/candidate_focus_evidence.csv`
- `interaction_diagnostic_candidate_package/candidate_stable_feature_evidence.csv`
- `interaction_diagnostic_candidate_package/candidate_claim_map.csv`

## External Context

OSD-120/GLDS-120 is a public NASA OSDR/GeneLab Arabidopsis thaliana CARA
root-tip RNA-seq study. The source study compares Col-0, Wassilewskija, and
Col-0 phyD root-tip responses to spaceflight, with light and dark treatment as
core experimental context.

Primary context sources:

- NASA OSDR OSD-120: https://osdr.nasa.gov/bio/repo/data/studies/OSD-120
- Paul et al. 2017 PLOS ONE: https://doi.org/10.1371/journal.pone.0180186
- 2024 npj Microgravity CARA light-context paper:
  https://www.nature.com/articles/s41526-024-00417-0

The literature supports the package scope: genotype/ecotype, phyD status, root
tips, and light/dark context are central to the experiment. It does not support
any benchmark-performance claim by itself; performance claims stay local to the
draft v9 artifacts.

## Package Summary

| field | value |
|---|---:|
| candidate | sparse_l1_c1 |
| variant | tvg2000_log1p_zscore_l1_c1 |
| primary genotype/ecotype BA | 0.9167 |
| secondary light-treatment BA | 0.8333 |
| diagnostic condition-stratum BA | 0.8889 |
| mean family BA | 0.8796 |
| minimum family BA | 0.8333 |
| nearest-centroid comparison | 9 improve / 2 tie / 0 worse |
| focus-fold comparison | 3 improve / 0 tie / 0 worse |
| stable features >=0.5 frequency | 19 |
| stable features >=0.8 frequency | 10 |

## Focus Evidence

| focus fold | nearest BA | candidate BA | delta | nonzero coefficients | stable >=0.5 | stable >=0.8 |
|---|---:|---:|---:|---:|---:|---:|
| Col.0.PhyD\|Dark.Treatment | 0.5000 | 0.6667 | 0.1667 | 10 | 10 | 7 |
| Light.Treatment | 0.5556 | 0.8333 | 0.2778 | 10 | 2 | 1 |
| Col.0.PhyD | 0.5000 | 0.9167 | 0.4167 | 12 | 7 | 2 |

The focus table is the right compact figure/table surface: it shows that the
candidate recovers all predefined fragile sentinels while keeping the evidence
traceable to the fold-detail, sparse coefficient, and stability tables.

## Stable Feature Evidence

`candidate_stable_feature_evidence.csv` emits 19 feature rows selected in at
least half of deterministic balanced train-fold subsamples. Ten are
`stable_ge_0_8`; nine are `stable_ge_0_5`. The strongest repeated signals are
fold-specific rather than universal, which is expected for an interaction task
whose source biology is explicitly genotype/ecotype and light-context dependent.

The package should not describe these genes as validated biomarkers. The safe
claim is narrower: they are auditable sparse-model evidence supporting the
current diagnostic candidate under this draft split.

## Claim Boundary

Supported:

- `sparse_l1_c1` is the primary draft transparent OSD-120 diagnostic candidate.
- It improves every predefined fragile focus fold versus nearest centroid.
- It has no worse fold than nearest centroid across the 11-fold ladder
  comparison.
- Its selected-feature evidence is sparse and auditable, with 19 stable feature
  rows at selection frequency >=0.5 across the three focus folds.
- OSD-120 is an appropriate local task for genotype/ecotype and light-treatment
  interaction diagnostics.

Not supported:

- No frozen v9 release claim.
- No leave-one-mission-out generalization claim.
- No cross-species or cross-study transfer claim.
- No biological validation or causal biomarker claim for selected genes.
- No operational plant-growth recommendation.

## Decision

Keep `sparse_l1_c1` as the packaged primary draft OSD-120 diagnostic candidate.
Use `candidate_package_summary.csv` for the one-row result table,
`candidate_focus_evidence.csv` for the main figure/table, and
`candidate_claim_map.csv` as the manuscript/poster guardrail.

Next block: build an OSD-120 figure/table draft or move to a paired control
package for the compact `sparse_l1_c0p3` comparator if the presentation needs a
stability-versus-performance contrast.

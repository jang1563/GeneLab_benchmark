# OSD-120 Interaction Baseline Ladder Review

Status: V9-MULTI-021 complete as of 2026-05-26.

This review consolidates the OSD-120 interaction nearest-centroid, dense L2
logistic, top-500 L2 sensitivity, sparse L1, and sparse L1 stability outputs
into one draft-only decision surface. It does not introduce a new model family
and does not promote a leaderboard claim.

Primary machine-readable outputs:

- `interaction_baseline_ladder/baseline_ladder_summary.csv`
- `interaction_baseline_ladder/baseline_ladder_summary.json`
- `interaction_baseline_ladder/baseline_ladder_focus_folds.csv`
- `interaction_baseline_ladder/baseline_ladder_focus_folds.json`

## Ladder Summary

| candidate | role | primary BA | secondary BA | diagnostic BA | vs nearest | focus dark | focus light | focus genotype |
|---|---|---:|---:|---:|---|---:|---:|---:|
| nearest_centroid_default | reference floor | 0.6667 | 0.6667 | 0.6667 | 0 improve / 11 tie / 0 worse | 0.5000 | 0.5556 | 0.5000 |
| l2_logistic_default | dense control | 0.7778 | 0.7500 | 0.8611 | 8 improve / 2 tie / 1 worse | 0.3333 | 0.7222 | 0.5833 |
| l2_logistic_top500_c1 | top-500 dense control | 0.8333 | 0.6667 | 0.8611 | 8 improve / 2 tie / 1 worse | 0.6667 | 0.5000 | 0.6667 |
| sparse_l1_c0p3 | compact sparse comparator | 0.8611 | 0.7778 | 0.9167 | 8 improve / 2 tie / 1 worse | 0.6667 | 0.8333 | 0.9167 |
| sparse_l1_c1 | leading transparent candidate | 0.9167 | 0.8333 | 0.8889 | 9 improve / 2 tie / 0 worse | 0.6667 | 0.8333 | 0.9167 |

## Interpretation

The nearest-centroid default remains the reference floor. It is balanced across
the three fold families at 0.6667 BA, but its three focus folds expose the
fragility that motivated the logistic branch: `Col.0.PhyD|Dark.Treatment` and
`Col.0.PhyD` sit at 0.5000 BA, while `Light.Treatment` sits at 0.5556 BA.

The dense L2 default improves aggregate performance, but it worsens the
`Col.0.PhyD|Dark.Treatment` diagnostic focus fold from 0.5000 to 0.3333. The
top-500 L2 sensitivity restores that dark focus fold to 0.6667, but it worsens
the `Light.Treatment` focus fold to 0.5000. Both dense variants remain useful
controls, not the main OSD-120 diagnostic candidate.

Sparse L1 `C=0.3` is the compact comparator. It reaches the strongest diagnostic
fold-family BA in this ladder at 0.9167 and gives strong focus-fold behavior,
but it still has one fold worse than nearest centroid in the full 11-fold
comparison.

Sparse L1 `C=1.0` is the leading transparent candidate. It has the best primary
and secondary family scores, no worse fold versus nearest centroid, and strong
focus-fold recovery. Stability audit evidence is also stronger than `C=0.3` for
`Col.0.PhyD|Dark.Treatment` and `Col.0.PhyD`: stable-feature counts at
frequency >=0.5 are 10 and 7 respectively, compared with 5 and 2 for `C=0.3`.

## Decision

Advance `sparse_l1_c1` as the primary draft OSD-120 sparse diagnostic candidate.
Retain `sparse_l1_c0p3` as the compact stability comparator. Retain nearest
centroid, dense L2 default, and top-500 L2 as controls for explaining the model
ladder and fragile-fold tradeoffs.

Next block: consolidate the OSD-120 diagnostic surface into a figure-ready
candidate package with explicit claim boundaries, rather than adding a new
model family immediately.

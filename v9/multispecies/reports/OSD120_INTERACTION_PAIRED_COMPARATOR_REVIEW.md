# OSD-120 Paired Comparator Review

Status: V9-MULTI-024 complete as of 2026-05-26.

This review decides whether the compact sparse L1 `C=0.3` comparator should be
shown beside the primary sparse L1 `C=1.0` candidate in the OSD-120 figure/table
surface.

Primary outputs:

- `interaction_paired_comparator_table/paired_comparator_summary.csv`
- `interaction_paired_comparator_table/paired_comparator_summary.json`
- `interaction_paired_comparator_table/paired_focus_comparator_table.csv`
- `interaction_paired_comparator_table/paired_focus_comparator_table.json`
- `interaction_paired_comparator_table/paired_comparator_decision.md`

## Summary

| comparison | sparse L1 C=1.0 | sparse L1 C=0.3 |
|---|---:|---:|
| mean family BA | 0.880 | 0.852 |
| minimum family BA | 0.833 | 0.778 |
| diagnostic family BA | 0.889 | 0.917 |
| nearest-centroid fold comparison | 9 improve / 2 tie / 0 worse | 8 improve / 2 tie / 1 worse |
| focus-fold BA comparison | 3 tied / 0 C=1 better / 0 C=0.3 better | 3 tied / 0 C=1 better / 0 C=0.3 better |
| focus nonzero coefficients | 32 | 19 |
| focus stable features >=0.5 | 19 | 8 |
| focus stable features >=0.8 | 10 | 6 |

## Focus Comparison

| focus fold | C=1.0 BA | C=0.3 BA | C=1 - C=0.3 | C=1 nonzero | C=0.3 nonzero |
|---|---:|---:|---:|---:|---:|
| Col.0.PhyD\|Dark.Treatment | 0.667 | 0.667 | +0.000 | 10 | 9 |
| Light.Treatment | 0.833 | 0.833 | +0.000 | 10 | 3 |
| Col.0.PhyD | 0.917 | 0.917 | +0.000 | 12 | 7 |

The paired table shows that `C=0.3` is genuinely useful as a compact comparator:
it reaches the same three focus-fold balanced accuracies with fewer nonzero
coefficients, especially for `Light.Treatment`.

However, `C=1.0` remains the better primary candidate. It has stronger mean and
minimum fold-family BA, stronger full-ladder behavior versus nearest centroid,
no nearest-centroid-worse fold, and more stable selected-feature evidence in the
dark and genotype focus folds.

## Decision

Do not add `sparse_l1_c0p3` as a second primary panel in the main OSD-120
figure. Keep the primary figure/table centered on `sparse_l1_c1`.

Use the paired comparator table only as an appendix or supplement when the
presentation needs to explain the stability-versus-compactness tradeoff.

Next block: build a final OSD-120 diagnostic artifact manifest that indexes the
ladder, candidate package, figure/table package, paired comparator, source
manifest, and tests as one traceable evidence set.

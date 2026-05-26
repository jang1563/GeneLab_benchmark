# OSD-120 Sparse L1 Paired Comparator Decision

`sparse_l1_c0p3` should remain an appendix or supplement comparator,
not a second primary panel in the main OSD-120 figure.

Reason:

- The three focus-fold balanced accuracies are tied between `c1` and `c0p3`.
- `c1` has stronger full-ladder behavior: 9 improve, 2 tie, 0 worse versus nearest centroid.
- `c0p3` has one nearest-centroid-worse fold in the full 11-fold comparison.
- `c0p3` is more compact in the focus folds, especially `Light.Treatment`, so it remains a useful control.

Recommended display:

- Main figure: keep the simpler `sparse_l1_c1` focus table.
- Appendix or supplement: include the paired comparator table when explaining the stability-versus-compactness tradeoff.

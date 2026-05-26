# OSD-120 Figure/Table Draft Review

Status: V9-MULTI-023 complete as of 2026-05-26.

This block turns the packaged `sparse_l1_c1` OSD-120 diagnostic candidate into
human-facing figure/table draft artifacts. It does not change the model, rerun
training, or expand the claim boundary.

Primary outputs:

- `interaction_figure_table_package/figure_main_focus_table.csv`
- `interaction_figure_table_package/figure_main_focus_table.json`
- `interaction_figure_table_package/figure_stable_feature_appendix.csv`
- `interaction_figure_table_package/figure_stable_feature_appendix.json`
- `interaction_figure_table_package/figure_caption.md`
- `interaction_figure_table_package/figure_claim_boundary_box.md`

## Main Table

| focus fold | nearest BA | sparse L1 BA | delta | nonzero coefficients | stable >=0.5 | stable >=0.8 |
|---|---:|---:|---:|---:|---:|---:|
| Condition stratum: Col.0.PhyD \| Dark | 0.500 | 0.667 | +0.167 | 10 | 10 | 7 |
| Light treatment: Light | 0.556 | 0.833 | +0.278 | 10 | 2 | 1 |
| Genotype/ecotype: Col.0.PhyD | 0.500 | 0.917 | +0.417 | 12 | 7 | 2 |

The main table is suitable for a poster or manuscript supplement because it
puts the three fragile sentinels on one surface and keeps each row tied to the
nearest-centroid baseline, sparse-model performance, non-zero coefficient
count, and stability-selection evidence.

## Appendix Table

`figure_stable_feature_appendix.csv` contains 19 stable sparse-model feature
rows. Ten are selected in at least 80% of deterministic balanced train-fold
subsamples and nine are selected in at least 50%. The appendix intentionally
labels these as sparse-model evidence, not validated biology.

Coefficient directions are split across the candidate feature evidence: 12 rows
are positive for the LEO/ISS class and 7 are negative for the LEO/ISS class.
This is an interpretability aid only; it is not a causal gene claim.

## Caption And Claim Box

`figure_caption.md` is written for direct reuse in a draft figure or supplement.
It states the task, candidate, fold sentinels, full ladder comparison, family
balanced accuracies, feature-stability caveat, and external context sources.

`figure_claim_boundary_box.md` is a shorter guardrail block. It keeps allowed
language to "packaged primary draft transparent diagnostic candidate" and
explicitly blocks frozen benchmark, leave-one-mission-out, cross-study,
validated biomarker, and operational recommendation claims.

## Decision

Use `figure_main_focus_table.csv` as the primary OSD-120 table surface and
`figure_stable_feature_appendix.csv` as the conservative appendix. The current
draft is ready for a poster/manuscript-style table.

Next block: build an optional paired comparator table for `sparse_l1_c0p3` if a
stability-versus-performance contrast is needed beside the primary `sparse_l1_c1`
candidate.

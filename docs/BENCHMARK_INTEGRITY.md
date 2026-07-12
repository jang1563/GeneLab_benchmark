# Benchmark Integrity And Label Visibility

## Public-label status

All labels in the released task package are public. This includes
`test_y.csv` in the two directories whose historical names end in `_holdout`:

- `tasks/A4_thymus_lomo/fold_RR-23_holdout/test_y.csv`
- `tasks/A5_skin_lomo/fold_RR-7_holdout/test_y.csv`

The `_holdout` suffix means that the mission was excluded from model training
for the recorded analysis. It does **not** mean that labels are hidden, private,
or available only to an evaluation server.

## Supported interpretation

The recorded RR-23 and RR-7 results are retrospective open-validation results.
They are useful for reproducing the released analysis and for checking that an
implementation follows the same mission-separated protocol. Because both
labels and historical results are public, they are not blind estimates for new
model development and should not be presented as private-holdout performance.

The current result files remain part of the scientific release:

- `evaluation/A4_holdout_results.json`
- `evaluation/A5_holdout_results.json`

Their filenames are retained for compatibility and provenance. Public prose
should call them **retrospective open validation**.

## New-model evaluation

Researchers may report scores on these folds if they explicitly say that the
folds and labels were public during model development. Such scores should be
treated as open benchmark comparisons, with the usual risk of adaptation to a
known test set.

A future blind evaluation requires a rotated evaluation set: a new eligible
mission or a versioned split whose labels are withheld from participants and
managed by an independent evaluator. No current repository file should be
described as providing that blind holdout.

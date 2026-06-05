# SpaceBio-Bench v9 Public Bulk Metadata Alpha

Status: metadata-only alpha subset, not a frozen payload release.

This directory contains the small public evidence subset referenced by the
SpaceBio-Bench system, evaluation, release-readiness, and claim-boundary cards.
It is intended to make the public bulk metadata-alpha claims inspectable from
the default GitHub branch without publishing payload matrices or draft extension
lanes.

## Included

- `task_manifests/*.json`: eight public bulk LOMO task manifests.
- `task_manifest_index.csv` and `.json`: task registry summary.
- `task_data_index.csv` and `.json`: fold-level row-count and path registry.
- `source_inventory.csv` and `.json`: 22 public OSDR source rows.
- `source_checksum_audit.csv` and `.json`: OSDR API and checksum-manifest
  evidence for the public bulk source rows.
- `datapackage.draft.json`: draft Frictionless Data Package descriptor for the
  metadata-only alpha subset.
- `reports/bulk_lomo_baseline_summary.csv` and `.json`: normalized summary of
  scaffold baselines.
- `reports/nearest_centroid/bulk_lomo_summary.csv` and `.json`: nearest-centroid
  scaffold baseline summary.
- `reports/sklearn_baselines/bulk_lomo_summary.csv` and `.json`: PCA-LR and
  L2 logistic-regression scaffold baseline summary.
- Per-task scaffold baseline `predictions.csv`, `metrics.json`, and
  `run_manifest.json` files referenced by `datapackage.draft.json`.
- `reports/public_bulk_alpha_gap_matrix/`: alpha-boundary gap and payload-hash
  blocker tables.
- `reports/public_bulk_alpha_snapshot_decision/`: metadata-alpha decision,
  allowed language, blocked language, and next-action tables.

## Excluded

- Fold-level payload matrices such as `train_X.csv` and `test_X.csv`.
- Local payload mirrors and payload-level SHA-256 manifests.
- Single-cell, organoid, multispecies, and other draft extension-lane outputs.
- Any DOI/archive-ready payload bundle.
- Any state-of-the-art leaderboard or clinical, crew-health, intervention, or
  countermeasure claim.

## Interpretation

The files here support claims about task/source/provenance metadata, scaffold
baseline plumbing, and release-readiness blockers. They do not support frozen
payload, biological mechanism, translational, or operational-readiness claims.

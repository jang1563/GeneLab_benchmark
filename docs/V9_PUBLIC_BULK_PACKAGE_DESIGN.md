# v9 Public Bulk Package Design

Status: draft design
Date: 2026-05-21

## Purpose

This document defines the public bulk packaging boundary for SpaceBio-Bench v9
before dataset-card, RO-Crate, Zenodo, or Hugging Face release language is
expanded.

The current v9 scaffold has enough structure to describe a package, but it is
not a frozen release. The correct next move is to separate metadata, data
payloads, provenance evidence, and benchmark outputs so future freeze work can
be explicit rather than implicit.

## External Rules Applied

Frictionless Data Package:

- A package descriptor is the central JSON object that lists data resources.
- Resources can point to local files or remote URLs.
- Resource paths must be POSIX-style relative paths, not absolute local paths.
- Tabular resources should carry machine-readable field schemas.

BagIt:

- A payload is not integrity-checked only because an upstream checksum manifest
  exists.
- A release-quality package needs payload-level hashes for the files actually
  distributed in the package.

OSDR:

- OSDR source identity should stay tied to OSD accessions and public OSDR API
  evidence.
- Public and controlled human data must remain separate.

## Current Artifact Classes

### Metadata Spine

These files are small, Git-friendly, and should remain in the repository:

- `v9/task_manifests/*.json`
- `v9/task_manifest_index.csv`
- `v9/task_manifest_index.json`
- `v9/task_data_index.csv`
- `v9/task_data_index.json`
- `v9/source_inventory.csv`
- `v9/source_inventory.json`
- `v9/source_checksum_audit.csv`
- `v9/source_checksum_audit.json`
- `v9/datapackage.draft.json`

Rationale:

- They define task identity, source identity, fold references, and provenance
  evidence.
- They are the minimum package descriptor layer an external user needs before
  downloading larger data payloads.

### Public Bulk Payload Bundle

These files are public-bulk benchmark payloads, but should be treated as a
downloadable bundle rather than the core Git metadata spine:

- legacy processed fold matrices under `tasks/*/fold_*/*`
- `train_X.csv`, `test_X.csv`
- `train_y.csv`, `test_y.csv`
- `train_meta.csv`, `test_meta.csv`
- `selected_genes.txt`
- `fold_info.json`

Current status:

- These payloads are indexed by `v9/task_data_index.csv`.
- They are not yet represented by payload-level package hashes.
- They should not be called frozen until the package builder computes hashes for
  every distributed payload file.

### Benchmark Output Bundle

These files are generated benchmark outputs:

- `v9/reports/bulk_lomo_baseline_summary.csv`
- `v9/reports/nearest_centroid/bulk_lomo_summary.csv`
- `v9/reports/sklearn_baselines/bulk_lomo_summary.csv`
- per-task `predictions.csv`
- per-task `metrics.json`
- per-task `run_manifest.json`

Current decision:

- Include summary tables and per-run output files in the draft descriptor.
- Keep them separate from public bulk payload data because they are derived
  outputs, not benchmark inputs.

### Deferred Or Excluded

Not part of the v9 public bulk alpha package:

- raw OSDR sequencing payload downloads
- gated or controlled human data
- local caches and virtual environments
- `submissions/`
- v8 intervention outputs
- Mars-regime or countermeasure artifacts

Raw OSDR payloads remain upstream source evidence. The v9 public bulk alpha
package should distribute the curated public benchmark inputs needed to run the
current bulk LOMO tasks, not mirror every OSDR raw file.

## Draft Data Package Descriptor

Generated file:

- `v9/datapackage.draft.json`

Generator:

```bash
python scripts/build_v9_datapackage_draft.py
```

Descriptor scope:

- package-level draft metadata
- sources pointing to OSDR API, source inventory, and checksum audit
- CSV resources with inferred Table Schema fields
- resource byte counts and SHA-256 hashes for small Git-tracked v9 artifacts
- file collections for task manifests and baseline outputs
- custom `spacebio_bench:*` fields for draft status and bundle partitioning

Current resources:

| Resource | Bundle part | Count |
|---|---:|---:|
| `task_manifest_index` | metadata spine | 1 |
| `task_data_index` | data index | 1 |
| `source_inventory` | metadata spine | 1 |
| `source_checksum_audit` | provenance evidence | 1 |
| `bulk_lomo_baseline_summary` | benchmark outputs | 1 |
| `nearest_centroid_summary` | benchmark outputs | 1 |
| `sklearn_baseline_summary` | benchmark outputs | 1 |
| `task_manifests` | metadata spine | 8 |
| `baseline_predictions` | benchmark outputs | 24 |
| `baseline_metrics` | benchmark outputs | 24 |
| `baseline_run_manifests` | benchmark outputs | 24 |

The descriptor deliberately says:

```json
"spacebio_bench:release_status": "draft_not_frozen"
```

and:

```json
"spacebio_bench:payload_verification_status": "checksum_manifests_parsed_payloads_not_hashed"
```

This language prevents checksum-manifest evidence from being mistaken for a
payload freeze.

## Future Package Layout

Recommended eventual release directory:

```text
spacebio-bench-v9-public-bulk/
  datapackage.json
  README.md
  task_manifests/
  indexes/
  provenance/
  data/
    bulk_lomo/
      tasks/
  reports/
    baselines/
  schemas/
  scripts/
```

Mapping from current repo:

| Current path | Future package path | Role |
|---|---|---|
| `v9/task_manifests/*.json` | `task_manifests/*.json` | task definitions |
| `v9/task_manifest_index.csv` | `indexes/task_manifest_index.csv` | task registry |
| `v9/task_data_index.csv` | `indexes/task_data_index.csv` | fold data index |
| `v9/source_inventory.csv` | `provenance/source_inventory.csv` | source table |
| `v9/source_checksum_audit.csv` | `provenance/source_checksum_audit.csv` | upstream checksum evidence |
| `tasks/*/fold_*/*` | `data/bulk_lomo/tasks/*/fold_*/*` | benchmark payload |
| `v9/reports/*` | `reports/baselines/*` | derived baseline outputs |
| `spacebio_bench/schemas/*.json` | `schemas/*.json` | validation schemas |

## Freeze Gates

### Gate P0: Draft Descriptor

Status: complete.

Evidence:

- `v9/datapackage.draft.json`
- resource paths are relative
- small v9 metadata and output artifacts have SHA-256 hashes
- CSV resources have schema fields

### Gate P1: Payload Manifest

Status: not started.

Required:

- enumerate every distributed public bulk payload file
- compute SHA-256 for each payload file
- record byte counts
- record source task and fold
- exclude local absolute paths

Expected output:

- `v9/public_bulk_payload_manifest.csv`
- `v9/public_bulk_payload_manifest.json`

### Gate P2: Payload Verification

Status: not started.

Required:

- verify local distributed payload files against the payload manifest
- optionally compare source-derived processed files against OSDR checksum
  manifest evidence where filenames align
- record mismatches as data, not crashes

Expected output:

- `v9/public_bulk_payload_verification.csv`
- `v9/public_bulk_payload_verification.json`

### Gate P3: Release Descriptor

Status: blocked by P1/P2.

Required:

- rename `datapackage.draft.json` to release `datapackage.json`
- finalize license and citation language
- replace draft status fields
- attach package-level checksum or BagIt manifest

## Next Work Block

Recommended next block:

`V9-THEN-004: Dataset card draft`

Reason:

- The public package boundary is now explicit.
- The dataset card can say what the draft package contains without claiming a
  payload freeze.
- Payload hashing can be a later dedicated block because it will touch large
  files and should be run deliberately.

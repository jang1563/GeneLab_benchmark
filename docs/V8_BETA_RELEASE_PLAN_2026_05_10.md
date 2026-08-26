# v8.0-beta Release-Candidate Plan

Date: 2026-05-10
Branch: v3
Current state: beta reproducibility hardening in progress.

## What Changed for Beta

The v8 alpha analysis outputs remain the scientific baseline. The beta work
adds release machinery around them:

- manifest validation for every promoted v8 result,
- an external input-freeze record,
- a public artifact split manifest,
- one-command HPC rebuild orchestration,
- release-gate integration so provenance failures block future releases.

## New Beta Entry Points

```bash
python scripts/validate_v8_provenance.py
bash scripts/hpc_v8_beta_rebuild.sh
```

For a final frozen beta, use:

```bash
python scripts/validate_v8_provenance.py --require-frozen
bash scripts/hpc_v8_beta_rebuild.sh --require-frozen
```

## Current Release-Candidate Status

Passed/ready:

- All existing v8 run manifests are structurally validated by the new checker.
- Exact tracked-output SHA-256 checks are verified for files with single-file
  checksums.
- `input_freeze.json` records the current OSDR, SpaceOmicsBench, processed
  fGSEA, L1000CDS2, and Enrichr/CRISPR input surfaces.
- `v8/provenance/external_source_audit.json` records the current
  SpaceOmicsBench processed-input bundle checksum
  (`fce14de513a6c10587c8b03d21e0f510ea8081afb99573f0067ed5983abaa774`) and
  file-level SHA-256s. The local source bundle has no git metadata, so this
  checksum is the release-candidate identity until an upstream tag or archived
  source bundle supersedes it.
- `v8/provenance/api_raw_archive_audit.json` records the current INTERVENE raw
  API-response archive checksum
  (`31327597a42184cb6bafb8058dd611725f888876ac5f010fad7ee11f3118fefd`) for
  30 L1000CDS2/Enrichr JSON response files. Re-parsing the archive produced
  byte-identical INTERVENE parsed outputs against the tracked CSV/JSON files.
- `v8_beta_artifact_manifest.json` records the GitHub, Hugging Face, Zenodo,
  and HPC/object-storage release split.
- `hpc_release_validate.sh` now runs provenance validation before summary
  regeneration.

Still open before declaring v8.0-beta frozen:

1. Copy the release-candidate raw API archive and other raw external inputs to
   durable Hugging Face, Zenodo, or object-storage targets and update the
   recorded archive root hints.
2. Run `scripts/hpc_v8_beta_rebuild.sh` from a fresh clean checkout after the
   input freeze is promoted.
3. Upload/record the Hugging Face dataset bundle revision.
4. Mint or reserve the Zenodo DOI and record it in the artifact manifest.

## Artifact Split

GitHub:

- code and HPC scripts,
- compact v8 evaluation JSON/CSV files,
- `v8/intervene/signatures/`,
- figures,
- provenance manifests and beta metadata.

Hugging Face dataset:

- compact public result tables,
- v8 evaluation summaries,
- intervention signatures,
- regenerated pathway summaries that are too data-like for the code-only
  release narrative.

Zenodo:

- immutable snapshot of the beta release package,
- DOI metadata linking code, data bundle, and manuscript draft.

HPC/object storage:

- raw OSDR caches,
- SpaceOmicsBench/SOMA source inputs,
- raw API response dumps if they are too large or license-sensitive for Git.

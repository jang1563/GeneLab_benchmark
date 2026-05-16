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
- `v8_beta_artifact_manifest.json` records the GitHub, Hugging Face, Zenodo,
  and HPC/object-storage release split.
- `hpc_release_validate.sh` now runs provenance validation before summary
  regeneration.

Still open before declaring v8.0-beta frozen:

1. Replace the SpaceOmicsBench pending marker with an exact upstream commit,
   release tag, or archived source checksum.
2. Pin a concrete L1000CDS2 db-version if the API exposes one, or archive raw
   response payloads externally with checksums.
3. Run `scripts/hpc_v8_beta_rebuild.sh` from a fresh clean checkout after the
   input freeze is promoted.
4. Upload/record the Hugging Face dataset bundle revision.
5. Mint or reserve the Zenodo DOI and record it in the artifact manifest.

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

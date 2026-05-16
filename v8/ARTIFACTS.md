# v8 Artifact Policy

v8 is an incubating analysis layer. Keep the Git repository focused on code,
lightweight result summaries, provenance, and publication-ready figures. Large
external downloads and regenerable caches belong outside Git.

## Track In Git

- `v8/**/*.py`: analysis drivers and figure-generation scripts.
- `v8/README.md`, `v8/RESULTS_SUMMARY.md`, `v8/ARTIFACTS.md`: human-readable
  orientation and current results.
- `v8/provenance/`: small JSON manifests that map each promoted result to its
  input data, command, output file, and validation status.
- Small evaluation summaries under `v8/**/evaluation/`: JSON/CSV outputs that
  are needed to reproduce manuscript tables and figures, provided each file is
  below the repository size gate.
- `v8/intervene/signatures/`: compact exported gene signatures used for LINCS
  and Perturb-seq queries.
- `v8/figures/`: publication figures if they are final or required for review.

## Do Not Track

- `v8/bridge/geo_cache/`: GEO annotation and series-matrix downloads. These can
  exceed GitHub and Git LFS file limits and are regenerable.
- `v8/**/__pycache__/` and `*.pyc`: Python bytecode.
- Local Codex settings such as `.claude/`.
- Raw SpaceOmicsBench or SOMA data. Configure those via `SPACEOMICS_ROOT`.
- Any single file larger than 50 MB unless it has been explicitly approved for a
  data-release repository instead of this code repository.

## External Artifact Targets

- Hugging Face dataset repo: compact public benchmark inputs, result tables, and
  small regenerated pathway summaries.
- Zenodo release archive: versioned release bundle and DOI snapshot.
- HPC scratch or object storage: raw external downloads, GEO caches, and
  intermediate matrices that are cheap to regenerate but expensive to version.

## Reproducibility Contract

Each v8 result used in the manuscript should have:

- source input path or external dataset accession,
- script path and command line,
- output file path,
- figure or table where it is consumed,
- checksum or release tag when the result becomes part of a frozen release.

Use `v8/provenance/manifest.schema.json` as the lightweight contract for these
records. Full RO-Crate export can wait until v8 beta, but missing manifests
should block promotion of exploratory results into manuscript claims.

## v8 Beta Release-Candidate Metadata

- `scripts/validate_v8_provenance.py` validates all promoted run manifests,
  exact tracked-output checksums, `v8/provenance/input_freeze.json`, and
  `v8/release/v8_beta_artifact_manifest.json`.
- `scripts/hpc_v8_beta_rebuild.sh` is the one-command HPC orchestrator for a
  beta rebuild: pillar recomputation, figures, summary regeneration, provenance
  validation, and final release hygiene.
- `v8/provenance/input_freeze.json` is currently `release_candidate`, not
  `frozen`, because SpaceOmicsBench and L1000CDS2 still need exact external
  version pins or archived raw payloads.
- `v8/release/v8_beta_artifact_manifest.json` defines the release split across
  GitHub, Hugging Face, Zenodo, and HPC/object storage.

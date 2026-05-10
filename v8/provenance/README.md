# v8 Provenance

This directory stores small, tracked provenance manifests for v8 result files.
The manifests are designed to be written on HPC after a run finishes and then
committed with the corresponding lightweight summaries.

Do not place raw OSDR downloads, GEO caches, SpaceOmicsBench inputs, API response
dumps, or intermediate matrices here. Store those in HPC scratch, object storage,
or the eventual public dataset artifact.

## Required Pattern

Each promoted v8 result should have a JSON manifest under `runs/` that records:

- claim or table/figure supported,
- pillar and script,
- exact command line,
- input accessions or file paths,
- output file paths,
- checksum or release tag when frozen,
- HPC job metadata when available,
- validation status.

The schema is intentionally lightweight and RO-Crate-compatible enough to export
later when v8 reaches beta.

## Beta Gate

Run the standard validator before promoting v8 artifacts:

```bash
python scripts/validate_v8_provenance.py
```

For final beta freezing, the stricter gate is:

```bash
python scripts/validate_v8_provenance.py --require-frozen
```

The stricter mode requires `v8/provenance/input_freeze.json` and
`v8/release/v8_beta_artifact_manifest.json` to be marked `frozen` with no open
release blockers. Until then, v8 remains a beta release candidate rather than a
frozen beta snapshot.

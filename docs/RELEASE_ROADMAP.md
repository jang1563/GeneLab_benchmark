# Release Roadmap

This document records the near-term release plan after the v7 public metadata
cleanup and the start of the v8 SpaceMed incubator.

Status update, 2026-05-10: the combined v7/v8 clean-checkout gate passed on
HPC at commit `b486aee77ae956b108a59ac762ebeb7b302e7928` with
`bash scripts/hpc_release_validate.sh --v8-summary`. Remaining release work is
artifact packaging, optional v7.1/v8 alpha tagging, and the separate v8 beta
reproducibility program.

## v7.1 Patch Release

Goal: make the current public benchmark surface consistent, portable, and safe
to publish.

Scope:

- Promote the stable public eye mission label `OSD-397`.
- Keep `TBD` only as a legacy storage alias where old processed filenames still
  require it.
- Refresh README, Hugging Face card, citation metadata, and submission docs.
- Keep v7 scPRINT2/GNN results portable by avoiding machine-specific checkpoint
  paths.
- Add release hygiene guards for local settings and v8 cache artifacts.

Recommended commit groups:

1. Public metadata and documentation consistency.
2. `A6_eye_lomo` / v2 LLM fold rename from `fold_TBD_test` to
   `fold_OSD-397_test`.
3. v3/v4/v7 script portability and result-surface refresh.
4. Regression tests and release hygiene checks.

Validation gate:

- Run review-fix tests on HPC.
- Run whitespace/diff checks before commit.
- Inspect staged diff with rename detection enabled:
  `git diff --cached -M --summary`.
- Confirm no staged `.claude/`, `v8/bridge/geo_cache/`, `__pycache__/`, or
  single file larger than 50 MB.

Current gate status: passed on 2026-05-10 from a clean clone of `origin/v3`.

## v8.0-alpha

Goal: keep v8 as a clearly labeled incubating layer while preserving enough
code and lightweight outputs for scientific review.

Implementation plan: see `docs/V8_IMPLEMENTATION_RESEARCH.md`.

Track:

- v8 analysis scripts.
- compact evaluation outputs required for tables and figures.
- exported intervention signatures.
- final or review-needed figures.
- `v8/README.md`, `v8/ARTIFACTS.md`, and result summaries.
- `v8/provenance/`: small run manifests that map claims to scripts, inputs,
  outputs, and validation status.

Do not track:

- GEO cache downloads.
- raw SpaceOmicsBench/SOMA inputs.
- Python bytecode or local Codex settings.
- large intermediate matrices that can be regenerated on HPC.

Validation gate:

- `SPACEOMICS_ROOT` must point to a local SpaceOmicsBench checkout.
- Each pillar should have a documented command line and output manifest.
- Claims in the v8 manuscript draft should map to a result file and script.
- Consolidated summary regeneration passed on 2026-05-10 from a clean clone of
  `origin/v3`.

## v8.0-beta

Goal: turn the incubator into a reproducible release candidate.

Required before beta:

- Add a machine-readable provenance manifest for each pillar.
- Freeze exact input dataset versions or accessions.
- Add one-command HPC entrypoints for BRIDGE, DECOMPOSE, INTERVENE, causal
  integration, and summary regeneration.
- Separate public release artifacts into GitHub code, Hugging Face dataset
  files, and Zenodo DOI archive.

Current beta progress, 2026-05-10:

- `scripts/validate_v8_provenance.py` validates all promoted run manifests plus
  beta input/artifact metadata and is part of `scripts/hpc_release_validate.sh`.
- `scripts/hpc_v8_beta_rebuild.sh` orchestrates BRIDGE, DECOMPOSE, INTERVENE,
  causal, figures, summary, provenance validation, and the release gate from a
  clean HPC checkout.
- `v8/provenance/input_freeze.json` records the current release-candidate
  external input freeze and the remaining blockers before it can be promoted to
  `frozen`.
- `v8/release/v8_beta_artifact_manifest.json` records the GitHub, Hugging Face,
  Zenodo, and HPC/object-storage artifact split for beta packaging.

Remaining blockers: exact SpaceOmicsBench upstream commit/release tag, concrete
L1000CDS2 db-version or archived raw API responses, one fresh
`hpc_v8_beta_rebuild.sh --require-frozen` pass after freeze promotion, and final
Hugging Face/Zenodo upload records.

## Long-Term Direction

- Keep GeneLab Benchmark v1-v7 focused on cross-mission transcriptomics model
  evaluation.
- Let v8 expand into translational astronaut-health modeling, but preserve a
  hard boundary between benchmark tasks, external human data dependencies, and
  countermeasure discovery.
- Add a cell-eval-aligned single-cell subtrack only where the metric semantics
  fit: DE overlap, direction match, and mission discrimination rather than
  replacing bulk LOMO AUROC/F1 metrics.

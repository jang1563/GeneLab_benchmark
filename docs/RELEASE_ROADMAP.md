# Release Roadmap

This document records the near-term release plan after the v7 public metadata
cleanup. The v8 SpaceMed incubator has been split out and lives on the
separate `v3` branch; it is intentionally excluded from this public release
surface and from the v7.1 patch.

Status update, 2026-05-10: the combined v7/v8 clean-checkout gate passed on
HPC at commit `b486aee77ae956b108a59ac762ebeb7b302e7928` (then on the
combined branch). Remaining public-release work is artifact packaging and
optional v7.1 alpha tagging. v8 beta reproducibility work continues on the
`v3` branch.

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

## v8.0-alpha and v8.0-beta (on `v3` branch)

The v8 SpaceMed incubator (BRIDGE / DECOMPOSE / INTERVENE / causal pillars)
is being developed on the `v3` branch and is not part of this public
release. Implementation plan, analysis scripts, evaluation outputs,
provenance manifests, and the beta-rebuild orchestration live there. v8
will be merged back into this line only after the beta freeze gate passes
and the GitHub / Hugging Face / Zenodo artifact split is finalized.

Until then, treat v8 claims (intervention, countermeasure, Mars
extrapolation) as out of scope for the v1-v7 benchmark paper and HF
dataset card.

## Long-Term Direction

- Keep GeneLab Benchmark v1-v7 focused on cross-mission transcriptomics model
  evaluation.
- Let v8 expand into translational astronaut-health modeling, but preserve a
  hard boundary between benchmark tasks, external human data dependencies, and
  countermeasure discovery.
- Add a cell-eval-aligned single-cell subtrack only where the metric semantics
  fit: DE overlap, direction match, and mission discrimination rather than
  replacing bulk LOMO AUROC/F1 metrics.

# v7 / v8 Closure Status

Date: 2026-05-10
Branch: v3
Scope: GeneLab Benchmark v7 public-release closure and v8 SpaceMed incubator
readiness. MLCB submission packaging is intentionally out of scope for this
status note.

## Executive Status

| Layer | Current status | Meaning |
|---|---|---|
| v7.1 public benchmark cleanup | Release-candidate gate passed | The public-facing values, OSD-397 naming, README/HF card, and canonical result table are aligned. Clean-checkout HPC release validation passed; public artifact packaging and optional tagging remain. |
| v8.0-alpha incubator | Alpha gate passed | Pillar scripts, lightweight outputs, figures, and provenance manifests are present and marked `hpc_validated`; v8 summary regeneration passed from a clean HPC clone. |
| v8.0-beta | Not complete | Beta still requires clean-checkout recomputation, frozen external dataset versions, manifest validation, and public artifact split. |

## Validation Record

Final v7/v8 gate completed on 2026-05-10 from a clean HPC clone of
`origin/v3` at commit `b486aee77ae956b108a59ac762ebeb7b302e7928`.

Command:

```bash
PYTHON_BIN=.hpc-venv/bin/python bash scripts/hpc_release_validate.sh --v8-summary
```

Result:

- Python regression tests: 47/47 passing.
- Whitespace and conflict-marker diff checks: passing.
- Release hygiene check: passing.
- v8 consolidated summary regeneration: passing.
- Clean clone status after regeneration: no tracked changes.

## v7.1 Closure

Completed:

- Canonical v7.1 public result source added:
  `docs/CANONICAL_RESULTS_V7_1.md`.
- README, Hugging Face card, paper outline, and data catalog synchronized to
  the canonical v7.1 result surface.
- Stable eye mission label `OSD-397` is visible in README, submission docs,
  tests, task directory names, and task metadata.
- `tasks/A6_eye_lomo/fold_OSD-397_test` exists; `fold_TBD_test` is not present
  under `tasks/A6_eye_lomo/`.
- Legacy `TBD` remains only where code maps old processed file/directory names
  to the public `OSD-397` label.
- v8 cache policy is documented so large local downloads remain outside Git.

Remaining before calling v7.1 fully released:

1. Confirm public artifact destinations:
   - GitHub for code and small result artifacts;
   - Hugging Face for public feature matrices and compact result bundles;
   - Zenodo DOI if a frozen v7.1 archive is desired.
2. Decide whether v7.1 gets a git tag.

## v8.0-alpha Readiness

Completed:

- v8 pillar layout exists for BRIDGE, DECOMPOSE, INTERVENE, causal, figures,
  evaluation outputs, and provenance.
- HPC entrypoints exist:
  - `scripts/hpc_v8_bridge.sh`
  - `scripts/hpc_v8_decompose.sh`
  - `scripts/hpc_v8_intervene.sh`
  - `scripts/hpc_v8_causal.sh`
  - `scripts/hpc_v8_summary.sh`
- Provenance manifests under `v8/provenance/runs/` are marked
  `hpc_validated`, include `generated_at`, and record
  `local_tests_run=false` plus `hpc_tests_run=true`.
- `v8/bridge/geo_cache/` may exist locally, but no files from that cache are
  tracked in Git.
- `v8/MANUSCRIPT_DRAFT.md` now frames Mars projections and intervention hits as
  hypothesis-generation and validation-prioritization results, not clinical or
  operational recommendations.

Remaining before calling v8.0-alpha fully closed:

1. Confirm `v8/RESULTS_SUMMARY.md` and `v8/MANUSCRIPT_DRAFT.md` only promote
   claims with manifest-backed outputs.
2. Do a final manuscript-language scan for:
   - operational crew-health recommendations;
   - Mars point-prediction language;
   - intervention or countermeasure claims not labeled exploratory.
3. If alpha is declared, tag or otherwise record the exact commit used for the
   alpha review package.

## v8.0-beta Gap

v8 beta is intentionally not done. Required before beta:

- Clean-checkout HPC recomputation of all promoted pillars.
- Frozen external source versions or accessions for OSDR, SOMA, and
  SpaceOmicsBench inputs.
- One-command pillar regeneration from the documented HPC entrypoints.
- Machine-checked manifest validation against
  `v8/provenance/manifest.schema.json`.
- Final figures regenerated from manifest-linked result files.
- Hugging Face dataset bundle and Zenodo DOI release plan.

## Immediate Next Actions

1. Keep MLCB submission materials separate from v7/v8 closure.
2. Decide whether to tag v7.1 from the validated commit.
3. Decide whether v8 should be tagged as `v8.0-alpha` or remain an
   incubator until beta requirements are scheduled.
4. Schedule the v8 beta work separately from the v7.1 patch release.

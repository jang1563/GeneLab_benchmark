#!/usr/bin/env bash
# One-command v8 beta rebuild gate. This orchestrates the pillar entrypoints,
# figure generation, summary regeneration, provenance validation, and the
# standard release hygiene gate from a clean HPC checkout.
set -euo pipefail

RUN_BRIDGE=1
RUN_DECOMPOSE=1
RUN_INTERVENE=1
RUN_CAUSAL=1
RUN_FIGURES=1
RUN_SUMMARY=1
RUN_RELEASE_GATE=1
RUN_INTERVENE_API=0
RUN_SIGNATURE_REFRESH=0
REQUIRE_FROZEN=0
DECOMPOSE_EXTRA_ARGS=()

for arg in "$@"; do
  case "$arg" in
    --skip-bridge)
      RUN_BRIDGE=0
      ;;
    --skip-decompose)
      RUN_DECOMPOSE=0
      ;;
    --skip-intervene)
      RUN_INTERVENE=0
      ;;
    --skip-causal)
      RUN_CAUSAL=0
      ;;
    --skip-figures)
      RUN_FIGURES=0
      ;;
    --skip-summary)
      RUN_SUMMARY=0
      ;;
    --skip-release-gate)
      RUN_RELEASE_GATE=0
      ;;
    --with-api)
      RUN_INTERVENE_API=1
      ;;
    --refresh-signatures)
      RUN_SIGNATURE_REFRESH=1
      ;;
    --skip-bootstrap)
      DECOMPOSE_EXTRA_ARGS+=("--skip-bootstrap")
      ;;
    --allow-incomplete-cache)
      DECOMPOSE_EXTRA_ARGS+=("--allow-incomplete-cache")
      ;;
    --require-frozen)
      REQUIRE_FROZEN=1
      ;;
    -h|--help)
      cat <<'USAGE'
Usage: bash scripts/hpc_v8_beta_rebuild.sh [options]

Options:
  --skip-bridge             Skip BRIDGE recomputation
  --skip-decompose          Skip DECOMPOSE recomputation
  --skip-intervene          Skip INTERVENE refresh
  --skip-causal             Skip causal DAG rebuild
  --skip-figures            Skip v8 figure regeneration
  --skip-summary            Skip v8 summary regeneration
  --skip-release-gate       Skip scripts/hpc_release_validate.sh --v8-summary
  --with-api                Re-call L1000CDS2 and Enrichr APIs
  --refresh-signatures      Refresh INTERVENE tissue signatures before API work
  --skip-bootstrap          Pass through to hpc_v8_decompose.sh
  --allow-incomplete-cache  Pass through to hpc_v8_decompose.sh
  --require-frozen          Require frozen beta manifests with no open blockers
USAGE
      exit 0
      ;;
    *)
      echo "Unknown argument: $arg" >&2
      exit 2
      ;;
  esac
done

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

PYTHON_BIN="${PYTHON_BIN:-python3}"
if [[ "$PYTHON_BIN" == "python3" && -x "$REPO_ROOT/.hpc-venv/bin/python" ]]; then
  PYTHON_BIN="$REPO_ROOT/.hpc-venv/bin/python"
fi

if [[ -n "$(git status --porcelain --untracked-files=no)" ]]; then
  echo "Tracked worktree changes are present; use a clean checkout for beta rebuilds." >&2
  git status --short --untracked-files=no >&2
  exit 3
fi

echo "[v8 BETA 0/8] Provenance preflight"
if [[ "$REQUIRE_FROZEN" -eq 1 ]]; then
  "$PYTHON_BIN" scripts/validate_v8_provenance.py --require-frozen
else
  "$PYTHON_BIN" scripts/validate_v8_provenance.py
fi

if [[ "$RUN_BRIDGE" -eq 1 ]]; then
  echo "[v8 BETA 1/8] BRIDGE"
  bash scripts/hpc_v8_bridge.sh
else
  echo "[v8 BETA 1/8] Skipping BRIDGE"
fi

if [[ "$RUN_DECOMPOSE" -eq 1 ]]; then
  echo "[v8 BETA 2/8] DECOMPOSE"
  bash scripts/hpc_v8_decompose.sh "${DECOMPOSE_EXTRA_ARGS[@]}"
else
  echo "[v8 BETA 2/8] Skipping DECOMPOSE"
fi

if [[ "$RUN_INTERVENE" -eq 1 ]]; then
  echo "[v8 BETA 3/8] INTERVENE"
  INTERVENE_ARGS=()
  if [[ "$RUN_SIGNATURE_REFRESH" -eq 1 ]]; then
    INTERVENE_ARGS+=("--refresh-signatures")
  fi
  if [[ "$RUN_INTERVENE_API" -eq 1 ]]; then
    INTERVENE_ARGS+=("--with-api")
  fi
  bash scripts/hpc_v8_intervene.sh "${INTERVENE_ARGS[@]}"
else
  echo "[v8 BETA 3/8] Skipping INTERVENE"
fi

if [[ "$RUN_CAUSAL" -eq 1 ]]; then
  echo "[v8 BETA 4/8] CAUSAL"
  bash scripts/hpc_v8_causal.sh
else
  echo "[v8 BETA 4/8] Skipping CAUSAL"
fi

if [[ "$RUN_FIGURES" -eq 1 ]]; then
  echo "[v8 BETA 5/8] Figures"
  "$PYTHON_BIN" v8/figures/generate_main_figures.py
else
  echo "[v8 BETA 5/8] Skipping figures"
fi

if [[ "$RUN_SUMMARY" -eq 1 ]]; then
  echo "[v8 BETA 6/8] Summary"
  bash scripts/hpc_v8_summary.sh
else
  echo "[v8 BETA 6/8] Skipping summary"
fi

echo "[v8 BETA 7/8] Provenance postflight"
if [[ "$REQUIRE_FROZEN" -eq 1 ]]; then
  "$PYTHON_BIN" scripts/validate_v8_provenance.py --require-frozen
else
  "$PYTHON_BIN" scripts/validate_v8_provenance.py
fi

if [[ "$RUN_RELEASE_GATE" -eq 1 ]]; then
  echo "[v8 BETA 8/8] Release validation"
  PYTHON_BIN="$PYTHON_BIN" bash scripts/hpc_release_validate.sh --v8-summary
else
  echo "[v8 BETA 8/8] Skipping release validation"
fi

echo "v8 beta rebuild gate complete. Review git diff, update manifests with run metadata, then tag only after --require-frozen passes."

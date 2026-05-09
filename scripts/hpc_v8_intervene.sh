#!/usr/bin/env bash
# HPC entrypoint for v8 INTERVENE. Signature export is local-data only. External
# API calls are opt-in because they are time/version/rate-limit sensitive.
set -euo pipefail

RUN_API=0
TOP_N=150

for arg in "$@"; do
  case "$arg" in
    --with-api)
      RUN_API=1
      ;;
    --top-n=*)
      TOP_N="${arg#*=}"
      ;;
    -h|--help)
      echo "Usage: bash scripts/hpc_v8_intervene.sh [--with-api] [--top-n=150]"
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

echo "[v8 INTERVENE 1/4] Export tissue signatures"
"$PYTHON_BIN" v8/intervene/export_signatures.py --top-n "$TOP_N"

if [[ "$RUN_API" -eq 1 ]]; then
  echo "[v8 INTERVENE 2/4] L1000CDS2 chemical reversal query"
  "$PYTHON_BIN" v8/intervene/lincs_query.py

  echo "[v8 INTERVENE 3/4] Multi-tissue Pareto scoring"
  "$PYTHON_BIN" v8/intervene/pareto_multi_tissue.py

  echo "[v8 INTERVENE 4/4] Enrichr CRISPR KO orthogonal query"
  "$PYTHON_BIN" v8/intervene/perturb_seq_orthog.py
else
  echo "[v8 INTERVENE 2/4] Skipping L1000CDS2 query; pass --with-api to run"
  echo "[v8 INTERVENE 3/4] Skipping Pareto refresh; requires fresh L1000CDS2 outputs"
  echo "[v8 INTERVENE 4/4] Skipping Enrichr CRISPR query; pass --with-api to run"
fi

echo "v8 INTERVENE complete. Update v8/provenance/runs/intervene_*.json with HPC metadata and checksums."

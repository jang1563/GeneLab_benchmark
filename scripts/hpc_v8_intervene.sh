#!/usr/bin/env bash
# HPC entrypoint for v8 INTERVENE. Signature refresh and external API calls are
# opt-in because raw-data availability, API versions, and rate limits can drift.
set -euo pipefail

RUN_SIGNATURE_EXPORT=0
RUN_API=0
TOP_N=150
API_RAW_DIR=""

for arg in "$@"; do
  case "$arg" in
    --refresh-signatures)
      RUN_SIGNATURE_EXPORT=1
      ;;
    --with-api)
      RUN_API=1
      ;;
    --top-n=*)
      TOP_N="${arg#*=}"
      ;;
    --api-raw-dir=*)
      API_RAW_DIR="${arg#*=}"
      ;;
    -h|--help)
      echo "Usage: bash scripts/hpc_v8_intervene.sh [--refresh-signatures] [--with-api] [--top-n=150] [--api-raw-dir=PATH]"
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

if [[ "$RUN_SIGNATURE_EXPORT" -eq 1 ]]; then
  echo "[v8 INTERVENE 1/6] Export tissue signatures"
  "$PYTHON_BIN" v8/intervene/export_signatures.py --top-n "$TOP_N"
else
  echo "[v8 INTERVENE 1/6] Skipping signature refresh; pass --refresh-signatures to run"
fi

if [[ "$RUN_API" -eq 1 ]]; then
  API_RAW_ARGS=()
  if [[ -n "$API_RAW_DIR" ]]; then
    API_RAW_ARGS+=("--raw-dir" "$API_RAW_DIR")
  fi

  echo "[v8 INTERVENE 2/6] L1000CDS2 chemical reversal query"
  "$PYTHON_BIN" v8/intervene/lincs_query.py "${API_RAW_ARGS[@]}"

  echo "[v8 INTERVENE 3/6] Multi-tissue Pareto scoring"
  "$PYTHON_BIN" v8/intervene/pareto_multi_tissue.py

  echo "[v8 INTERVENE 4/6] Enrichr CRISPR KO orthogonal query"
  "$PYTHON_BIN" v8/intervene/perturb_seq_orthog.py "${API_RAW_ARGS[@]}"
else
  echo "[v8 INTERVENE 2/6] Skipping L1000CDS2 query; pass --with-api to run"
  echo "[v8 INTERVENE 3/6] Skipping Pareto refresh; requires fresh L1000CDS2 outputs"
  echo "[v8 INTERVENE 4/6] Skipping Enrichr CRISPR query; pass --with-api to run"
fi

echo "[v8 INTERVENE 5/6] API snapshot manifest for tracked outputs"
"$PYTHON_BIN" v8/intervene/api_snapshot_manifest.py

echo "[v8 INTERVENE 6/6] Safety triage for tracked candidates"
"$PYTHON_BIN" v8/intervene/safety_triage.py

echo "v8 INTERVENE complete. Update v8/provenance/runs/intervene_*.json with HPC metadata and checksums."

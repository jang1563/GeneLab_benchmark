#!/bin/bash
set -euo pipefail

CONDA_ENV_BASE="/home/fs01/jak4013/miniconda3/miniconda3/envs/seurat.v5.R.4.3.3.python.3.11.nfcore"
PYTHON_BIN="${CONDA_ENV_BASE}/bin/python"
SCRIPT="/athena/masonlab/scratch/users/jak4013/rrrm1_scrna/rrrm1_singlecell_hardening.py"

if [[ ! -x "${PYTHON_BIN}" ]]; then
    echo "ERROR: Python env not found at ${PYTHON_BIN}"
    exit 1
fi

if [[ ! -f "${SCRIPT}" ]]; then
    echo "ERROR: Hardening script not found at ${SCRIPT}"
    exit 1
fi

export LD_LIBRARY_PATH="${CONDA_ENV_BASE}/lib:${LD_LIBRARY_PATH:-}"
export MPLBACKEND=Agg

exec "${PYTHON_BIN}" "${SCRIPT}" "$@"

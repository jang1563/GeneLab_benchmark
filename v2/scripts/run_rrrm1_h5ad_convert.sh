#!/bin/bash
set -euo pipefail

# Wrapper for RRRM-1 STARsolo -> h5ad conversion on Cayuga.

CONDA_ENV_BASE="/home/fs01/jak4013/miniconda3/miniconda3/envs/seurat.v5.R.4.3.3.python.3.11.nfcore"
PYTHON_BIN="${CONDA_ENV_BASE}/bin/python"
CONVERTER="/athena/masonlab/scratch/users/jak4013/rrrm1_scrna/rrrm1_h5ad_convert.py"

if [[ ! -x "${PYTHON_BIN}" ]]; then
    echo "ERROR: Python env not found at ${PYTHON_BIN}"
    exit 1
fi

if [[ ! -f "${CONVERTER}" ]]; then
    echo "ERROR: Converter script not found at ${CONVERTER}"
    exit 1
fi

export LD_LIBRARY_PATH="${CONDA_ENV_BASE}/lib:${LD_LIBRARY_PATH:-}"

if [[ $# -eq 0 ]]; then
    echo "Running RRRM-1 h5ad conversion for all completed OSDs"
    exec "${PYTHON_BIN}" "${CONVERTER}"
fi

for osd in "$@"; do
    echo "Running RRRM-1 h5ad conversion for OSD-${osd}"
    "${PYTHON_BIN}" "${CONVERTER}" --osd "${osd}"
done

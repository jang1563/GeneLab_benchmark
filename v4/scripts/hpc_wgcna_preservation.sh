#!/bin/bash
#SBATCH --job-name=wgcna_pres
#SBATCH --partition=scu-gpu
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --output=v4/logs/wgcna_preservation_%j.out
#SBATCH --error=v4/logs/wgcna_preservation_%j.err

# wgcna_preservation_py.py — module preservation across 15 tissue pairs (Python)
# Run AFTER all 6 wgcna_analysis_py.py jobs complete.
#
# Usage:
#   sbatch v4/scripts/hpc_wgcna_preservation.sh

set -eo pipefail

PROJECT_DIR="/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark"
LOG_DIR="$PROJECT_DIR/v4/logs"
mkdir -p "$LOG_DIR"

echo "========================================"
echo "WGCNA Preservation (C(6,2)=15 pairs)"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "========================================"

CONDA_BASE="/home/fs01/jak4013/miniconda3/miniconda3"
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate perturb_seq_new
set -u

# Validate all 6 tissue outputs present
LOMO_TISSUES=(liver gastrocnemius kidney thymus eye skin)
for TISSUE in "${LOMO_TISSUES[@]}"; do
  MOD_FILE="$PROJECT_DIR/v4/wgcna_outputs/$TISSUE/module_assignments.csv"
  if [ ! -f "$MOD_FILE" ]; then
    echo "WARNING: $MOD_FILE not found — tissue $TISSUE will be skipped"
  fi
done

SCRIPT="$PROJECT_DIR/v4/scripts/wgcna_preservation_py.py"
echo "Running wgcna_preservation_py.py ..."
python -u "$SCRIPT"

echo "========================================"
echo "Done at $(date)"
echo "Output: $PROJECT_DIR/v4/evaluation/WGCNA_preservation.json"
echo "========================================"

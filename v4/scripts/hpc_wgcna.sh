#!/bin/bash
#SBATCH --job-name=wgcna_py
#SBATCH --array=0-5
#SBATCH --partition=scu-gpu
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --output=v4/logs/wgcna_%A_%a.out
#SBATCH --error=v4/logs/wgcna_%A_%a.err

# wgcna_analysis_py.py SLURM array job (PyWGCNA — Python implementation)
# 6 LOMO tissues: liver(0), gastrocnemius(1), kidney(2), thymus(3), eye(4), skin(5)
#
# Prerequisites:
#   python v4/scripts/wgcna_prep.py   (run on login node first)
#   pip install PyWGCNA dynamicTreeCut (already done in perturb_seq_new)
#
# Usage:
#   sbatch v4/scripts/hpc_wgcna.sh              # all 6 tissues in parallel
#   sbatch --array=0 v4/scripts/hpc_wgcna.sh    # liver only

set -eo pipefail

TISSUES=(liver gastrocnemius kidney thymus eye skin)
TISSUE=${TISSUES[$SLURM_ARRAY_TASK_ID]}

PROJECT_DIR="/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark"
SCRIPT="$PROJECT_DIR/v4/scripts/wgcna_analysis_py.py"
LOG_DIR="$PROJECT_DIR/v4/logs"
mkdir -p "$LOG_DIR" "$PROJECT_DIR/v4/wgcna_outputs/$TISSUE"

echo "========================================"
echo "WGCNA (PyWGCNA): $TISSUE  (task $SLURM_ARRAY_TASK_ID)"
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "========================================"

# Activate conda (disable -u to allow conda's unbound vars)
CONDA_BASE="/home/fs01/jak4013/miniconda3/miniconda3"
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate perturb_seq_new
set -u

# Validate packages
python -c "
from PyWGCNA import WGCNA
import dynamicTreeCut
print('PyWGCNA OK')
print('dynamicTreeCut OK')
"

# Check input file
EXPR_FILE="$PROJECT_DIR/v4/wgcna_inputs/${TISSUE}_expr.csv"
if [ ! -f "$EXPR_FILE" ]; then
  echo "ERROR: $EXPR_FILE not found — run wgcna_prep.py first"
  exit 1
fi
echo "Input: $EXPR_FILE"

# Run WGCNA (unbuffered output for SLURM log)
echo "Running wgcna_analysis_py.py --tissue $TISSUE ..."
python -u "$SCRIPT" --tissue "$TISSUE"

echo "========================================"
echo "Done: $TISSUE at $(date)"
OUTPUT_DIR="$PROJECT_DIR/v4/wgcna_outputs/$TISSUE"
ls -lh "$OUTPUT_DIR"/ 2>/dev/null || echo "No output directory found"
echo "========================================"

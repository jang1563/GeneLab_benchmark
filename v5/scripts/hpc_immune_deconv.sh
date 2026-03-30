#!/bin/bash
#SBATCH --job-name=v5_immune
#SBATCH --array=0-7
#SBATCH --partition=scu-cpu
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=v5/logs/immune_%a.out
#SBATCH --error=v5/logs/immune_%a.err

set -eo pipefail

TISSUES=(liver gastrocnemius kidney thymus eye skin lung colon)
TISSUE=${TISSUES[$SLURM_ARRAY_TASK_ID]}

echo "=== v5 Phase 1: Immune Deconvolution ==="
echo "Tissue: $TISSUE"
echo "Node: $(hostname)"
echo "Date: $(date)"

# Activate conda
source ~/.bashrc
conda activate perturb_seq_new

python -u v5/scripts/immune_deconv.py --tissue "$TISSUE"

echo "=== Done: $TISSUE ==="

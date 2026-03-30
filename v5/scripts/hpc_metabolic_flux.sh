#!/bin/bash
#SBATCH --job-name=v5_flux
#SBATCH --array=0-5
#SBATCH --partition=scu-cpu
#SBATCH --time=04:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --output=v5/logs/flux_%a.out
#SBATCH --error=v5/logs/flux_%a.err

set -eo pipefail

# 6 LOMO tissues only (lung/colon excluded — single-mission)
TISSUES=(liver gastrocnemius kidney thymus eye skin)
TISSUE=${TISSUES[$SLURM_ARRAY_TASK_ID]}

echo "=== v5 Phase 3: Metabolic Flux Modeling ==="
echo "Tissue: $TISSUE"
echo "Node: $(hostname)"
echo "Date: $(date)"

# Activate conda
source ~/.bashrc
conda activate perturb_seq_new

python -u v5/scripts/metabolic_flux.py --tissue "$TISSUE"

echo "=== Done: $TISSUE ==="

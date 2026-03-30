#!/bin/bash
#SBATCH --job-name=v4_shap_int
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=v4/logs/shap_int_%A_%a.out
#SBATCH --error=v4/logs/shap_int_%A_%a.err
#SBATCH --array=0-3

# GeneLabBench v4 Phase 4: SHAP interactions (tree methods only)
# 4 combinations: kidney×{xgb,rf}, eye×xgb, colon×xgb
#
# Usage: sbatch v4/scripts/hpc_shap_interactions.sh

set -eo pipefail
source ~/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate perturb_seq_new

cd /athena/cayuga_0003/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark

TISSUES=("kidney" "kidney" "eye"  "colon")
METHODS=("xgb"    "rf"     "xgb"  "xgb")

TISSUE=${TISSUES[$SLURM_ARRAY_TASK_ID]}
METHOD=${METHODS[$SLURM_ARRAY_TASK_ID]}

echo "=== Task $SLURM_ARRAY_TASK_ID: ${TISSUE} × ${METHOD} interactions ==="
echo "Start: $(date)"

python -u v4/scripts/shap_interactions.py --tissue "$TISSUE" --method "$METHOD"

echo "End: $(date)"

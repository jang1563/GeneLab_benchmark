#!/bin/bash
#SBATCH --job-name=v4_collect_pred
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=v4/logs/collect_pred_%A_%a.out
#SBATCH --error=v4/logs/collect_pred_%A_%a.err
#SBATCH --array=1-7  # 7 tissues (liver already done)

# GeneLabBench v4: Collect predictions per tissue (for DeLong test)
# Usage: sbatch v4/scripts/hpc_collect_predictions.sh

set -eo pipefail
source ~/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate perturb_seq_new

cd ${GENELAB_ROOT:?Set GENELAB_ROOT to your project directory}

TISSUES=("liver" "gastrocnemius" "kidney" "thymus" "eye" "skin" "lung" "colon")
TISSUE=${TISSUES[$SLURM_ARRAY_TASK_ID]}

echo "=== Task $SLURM_ARRAY_TASK_ID: $TISSUE ==="
echo "Start: $(date)"

python -u v4/scripts/collect_predictions.py --tissue "$TISSUE" --output "v4/evaluation/META_predictions_${TISSUE}.json" --no-verify

echo "End: $(date)"

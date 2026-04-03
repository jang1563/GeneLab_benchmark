#!/bin/bash
#SBATCH --job-name=v4_shap
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=v4/logs/shap_%A_%a.out
#SBATCH --error=v4/logs/shap_%A_%a.err
#SBATCH --array=0-15

# GeneLabBench v4 Phase 4: SHAP per tissue × method (16 combinations)
# Fast methods (tree/linear): ~5-30 min each
# Kernel methods (tasks 0,1,6): ~2-6 hours each (top-2000 genes)
#
# Usage: sbatch v4/scripts/hpc_shap.sh

set -eo pipefail
source ~/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate perturb_seq_new

cd ${GENELAB_ROOT:?Set GENELAB_ROOT to your project directory}

# 16 tissue × method combinations (ordered: kernel first for priority)
TISSUES=(    "liver"    "liver"    "gastrocnemius" "gastrocnemius" "kidney" "kidney" "thymus" "thymus" "eye"    "eye" "skin"          "skin"   "lung"          "lung"   "colon"  "colon")
METHODS=(    "svm_rbf"  "mlp"     "elasticnet_lr" "pca_lr"        "xgb"    "rf"     "knn"    "pca_lr" "pca_lr" "xgb" "elasticnet_lr" "pca_lr" "elasticnet_lr" "pca_lr" "pca_lr" "xgb")

TISSUE=${TISSUES[$SLURM_ARRAY_TASK_ID]}
METHOD=${METHODS[$SLURM_ARRAY_TASK_ID]}

echo "=== Task $SLURM_ARRAY_TASK_ID: ${TISSUE} × ${METHOD} ==="
echo "Start: $(date)"

python -u v4/scripts/shap_multi_method.py --tissue "$TISSUE" --method "$METHOD"

echo "End: $(date)"

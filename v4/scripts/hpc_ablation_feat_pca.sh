#!/bin/bash
#SBATCH --job-name=v4_abl_fp
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=v4/logs/ABL_FP_%A_%a.out
#SBATCH --error=v4/logs/ABL_FP_%A_%a.err
#SBATCH --array=0-159  # feat(112) + pca(48) = 160 jobs

# GeneLabBench v4: Feature Count + PCA Ablation
# Feature count: 7 K-values × 2 methods × 8 tissues = 112 (indices 0-111)
# PCA components: 6 nc-values × 8 tissues = 48 (indices 112-159)

set -eo pipefail
source ~/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate perturb_seq_new

cd ${GENELAB_ROOT:?Set GENELAB_ROOT to your project directory}

TISSUES=(liver gastrocnemius kidney thymus eye skin lung colon)
FEAT_METHODS=(pca_lr elasticnet_lr)
TOP_K_VALUES=(100 500 1000 2000 5000 10000 all)
NC_VALUES=(5 10 20 50 100 200)

N_TISSUES=${#TISSUES[@]}
N_FEAT_METHODS=${#FEAT_METHODS[@]}
N_TOP_K=${#TOP_K_VALUES[@]}
N_NC=${#NC_VALUES[@]}

FEAT_SIZE=$((N_TOP_K * N_FEAT_METHODS * N_TISSUES))  # 112

TASK_ID=${SLURM_ARRAY_TASK_ID}

if [ $TASK_ID -lt $FEAT_SIZE ]; then
    # Feature count ablation
    TISSUE_IDX=$((TASK_ID / (N_FEAT_METHODS * N_TOP_K)))
    REMAINDER=$((TASK_ID % (N_FEAT_METHODS * N_TOP_K)))
    METHOD_IDX=$((REMAINDER / N_TOP_K))
    K_IDX=$((REMAINDER % N_TOP_K))

    TISSUE=${TISSUES[$TISSUE_IDX]}
    METHOD=${FEAT_METHODS[$METHOD_IDX]}
    TOP_K=${TOP_K_VALUES[$K_IDX]}

    echo "=========================================="
    echo "Feature Ablation: ${TISSUE} / ${METHOD} / K=${TOP_K}"
    echo "Job: ${TASK_ID}, Time: $(date), Node: $(hostname)"
    echo "=========================================="

    python -u v4/scripts/ablation_feature_count.py \
        --tissue "$TISSUE" \
        --method "$METHOD" \
        --top-k "$TOP_K" \
        --n-bootstrap 2000 \
        --n-perm 1000
else
    # PCA component ablation
    IDX=$((TASK_ID - FEAT_SIZE))
    TISSUE_IDX=$((IDX / N_NC))
    NC_IDX=$((IDX % N_NC))

    TISSUE=${TISSUES[$TISSUE_IDX]}
    NC=${NC_VALUES[$NC_IDX]}

    echo "=========================================="
    echo "PCA Ablation: ${TISSUE} / nc=${NC}"
    echo "Job: ${TASK_ID}, Time: $(date), Node: $(hostname)"
    echo "=========================================="

    python -u v4/scripts/ablation_pca_components.py \
        --tissue "$TISSUE" \
        --n-components "$NC" \
        --n-bootstrap 2000 \
        --n-perm 1000
fi

echo "Done: $(date)"

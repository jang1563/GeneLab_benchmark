#!/bin/bash
#SBATCH --job-name=annot_404
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/logs/annotate_404.log
#SBATCH --error=/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/v3/logs/annotate_404.err

echo "=== GLDS-404 PBMC Cell Type Annotation ==="
echo "Node: $(hostname)"
echo "Date: $(date)"

WORK_DIR="/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark"

export LD_LIBRARY_PATH=""
source /home/fs01/jak4013/miniconda3/miniconda3/etc/profile.d/conda.sh
set -eo pipefail
conda activate scgpt_env_new

python3 -u "${WORK_DIR}/v3/scripts/annotate_glds404_pbmc.py" \
    --data-dir "${WORK_DIR}/v3/data/rrrm2" \
    --method per_cluster \
    --cluster-col seurat_clusters

echo ""
echo "=== Annotation Complete ==="
echo "Date: $(date)"

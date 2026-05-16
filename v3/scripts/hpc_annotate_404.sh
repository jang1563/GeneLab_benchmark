#!/bin/bash
#SBATCH --job-name=annot_404
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

echo "=== GLDS-404 PBMC Cell Type Annotation ==="
echo "Node: $(hostname)"
echo "Date: $(date)"

WORK_DIR="${SCRATCH_DIR:?Set SCRATCH_DIR}/huggingface/benchmark/GeneLab_benchmark"

export LD_LIBRARY_PATH=""
source ${CONDA_PREFIX:-$HOME/miniconda3}/etc/profile.d/conda.sh
set -eo pipefail
conda activate scgpt_env_new

python3 -u "${WORK_DIR}/v3/scripts/annotate_glds404_pbmc.py" \
    --data-dir "${WORK_DIR}/v3/data/rrrm2" \
    --method per_cluster \
    --cluster-col seurat_clusters

echo ""
echo "=== Annotation Complete ==="
echo "Date: $(date)"

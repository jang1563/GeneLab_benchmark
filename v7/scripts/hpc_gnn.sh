#!/bin/bash
#SBATCH --job-name=gnn_wgcna
#SBATCH --array=0-17              # 6 tissues × 3 graph types = 18 jobs
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=16:00:00
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err

# GNN on WGCNA topology — v7 GeneLabBench
# Array job: 6 tissues × 3 graph types = 18 combinations
#
# Graph types: wgcna (biological TOM), random (same edge count, random),
#              no_edges (global pooling = MLP baseline)
#
# Results in: GeneLab_benchmark/v7/evaluation/GNN_{tissue}_{graph_type}.json

set -eo pipefail
mkdir -p ${HOME}/genelab_logs

TISSUES=(liver gastrocnemius kidney thymus eye skin)
GRAPH_TYPES=(wgcna random no_edges)

# Map SLURM array index to tissue × graph_type
TISSUE_IDX=$(( SLURM_ARRAY_TASK_ID / 3 ))
GRAPH_IDX=$(( SLURM_ARRAY_TASK_ID % 3 ))

TISSUE=${TISSUES[$TISSUE_IDX]}
GRAPH=${GRAPH_TYPES[$GRAPH_IDX]}

echo "Job ID: $SLURM_JOB_ID"
echo "Array task: $SLURM_ARRAY_TASK_ID → tissue=$TISSUE, graph=$GRAPH"
echo "Node: $SLURMD_NODENAME"
echo "GPU: $CUDA_VISIBLE_DEVICES"
date

source ${CONDA_PREFIX:-$HOME/miniconda3}/etc/profile.d/conda.sh
conda activate perturb_seq_new

SCRIPT_DIR="${GENELAB_ROOT}/v7/unified"
cd "${GENELAB_ROOT}"

python -u "$SCRIPT_DIR/gnn_wgcna.py" \
    --tissue "$TISSUE" \
    --graph-type "$GRAPH" \
    --topology-scope train_fold \
    --n-edges-per-gene 10 \
    --n-bootstrap 1000 \
    --n-perm 100 \
    --seed 42

echo "Done: $TISSUE / $GRAPH"
date

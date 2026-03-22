#!/bin/bash
#SBATCH --job-name=v4_gsva
#SBATCH --partition=scu-cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=v4/logs/gsva_%j.out
#SBATCH --error=v4/logs/gsva_%j.err

# Generate GSVA/ssGSEA pathway scores for all 8 tissues

set -eo pipefail
source ~/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate perturb_seq_new

cd /athena/cayuga_0003/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark

echo "Starting GSVA generation: $(date)"
python -u v4/scripts/generate_pathway_scores.py --all-tissues --db hallmark,kegg
echo "Done: $(date)"

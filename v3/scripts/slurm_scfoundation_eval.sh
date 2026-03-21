#!/bin/bash
#SBATCH --job-name=scf_eval
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=6:00:00
#SBATCH --output=v3/logs/scf_eval_%j.out
#SBATCH --error=v3/logs/scf_eval_%j.err

set -euo pipefail

PROJECT_DIR="/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark"
CONDA_BASE="/home/fs01/jak4013/miniconda3/miniconda3"
SCF_ENV="$CONDA_BASE/envs/scfoundation_env"
SCF_DIR="/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/scFoundation"
SCF_CKPT="$SCF_DIR/model/models/models.ckpt"
DEFAULT_SCF_ARGS="--tissue liver"

echo "=== scFoundation Evaluation ==="
echo "Date: $(date)"
echo "Node: $(hostname)"
echo "GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo 'N/A')"
echo "Project dir: $PROJECT_DIR"
echo "Conda base: $CONDA_BASE"
echo "scFoundation env: $SCF_ENV"
echo "scFoundation dir: $SCF_DIR"

source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$SCF_ENV"

which python
python -V
python - <<'PY'
import sys
for name in ["torch", "scanpy", "pandas", "numpy"]:
    __import__(name)
print("python_executable", sys.executable)
PY

# Clone scFoundation if not exists
if [ ! -d "$SCF_DIR" ]; then
    echo "Cloning scFoundation..."
    git clone https://github.com/biomap-research/scFoundation.git "$SCF_DIR"
fi

if [ ! -f "$SCF_CKPT" ]; then
    echo "[ERROR] Missing scFoundation weights: $SCF_CKPT"
    if [ -f "$SCF_DIR/model/models/download.txt" ]; then
        echo "Download instructions:"
        cat "$SCF_DIR/model/models/download.txt"
    fi
    exit 1
fi

cd "$PROJECT_DIR"
mkdir -p v3/evaluation v3/logs

read -r -a SCF_ARGV <<< "${SCF_ARGS:-$DEFAULT_SCF_ARGS}"
echo "scFoundation args: ${SCF_ARGV[*]}"

python -u v3/scripts/scfoundation_eval.py \
    --scf-dir "$SCF_DIR" \
    "${SCF_ARGV[@]}" \
    2>&1

echo "=== scFoundation evaluation complete: $(date) ==="

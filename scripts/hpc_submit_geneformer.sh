#!/bin/bash
# ── GeneLab Benchmark: Geneformer LOMO Fine-Tuning (Cornell Cayuga HPC) ───────
#
# Cluster: Cayuga (Cornell CAC) — Slurm v25.05.0
# VPN required; contact scu@med.cornell.edu for access.
#
# GPU nodes (scu-gpu partition covers g0001-g0003):
#   g0001:       128 CPUs, 1024GB RAM, 4× NVIDIA A100 80GB PCIe
#   g0002-g0003: 128 CPUs, 1024GB RAM, 4× NVIDIA A40  48GB PCIe each
#   g0004:       128 CPUs, 1032GB RAM, 4× NVIDIA H100 (restricted — NOT in scu-gpu)
#
# Max time: 7 days on scu-gpu partition
# GRES syntax: --gres=gpu:a40:<n>  or  --gres=gpu:a100:<n>
#
# Storage:
#   Home dir:  /home/fs01/<cwid>     (NFS, NOT backed up — small files only)
#   Athena:    /athena/cayuga_XXXX/scratch/$USER/  ← use this for all job I/O
#   Lab link:  /athena/<labname>     → /athena/cayuga_XXXX
#
# ONE-TIME SETUP (run on login node, NOT in a job):
#   1. Miniconda already installed at:
#      /athena/masonlab/scratch/users/jak4013/miniconda3
#
#   2. mouse_gf environment (Python 3.8.10, CUDA 12.1):
#      conda create -n mouse_gf python=3.8.10 -y
#      conda activate mouse_gf
#      pip install torch==2.4.1+cu121 --index-url https://download.pytorch.org/whl/cu121
#      pip install transformers datasets accelerate scikit-learn safetensors tqdm huggingface_hub
#
#   3. Mouse-Geneformer source (no setup.py — use PYTHONPATH):
#      cd /home/fs01/jak4013
#      git clone https://github.com/zou-group/Mouse-Geneformer.git
#      # Model weights in: ${PROJECT_DIR}/models/mouse_gf_base/
#      # Dict files in:    ${PROJECT_DIR}/models/
#
#   4. Pre-tokenize on login node (CPU only, ~5 min):
#      cd /athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark
#      export PYTHONPATH="/home/fs01/jak4013/Mouse-Geneformer:$PYTHONPATH"
#      conda activate mouse_gf
#      python scripts/geneformer_tokenize.py --task A1 --model-version mouse_gf
#      # Tokenized data saved to tasks/A1_liver_lomo/fold_*/geneformer_tokens/mouse_gf/
#
# SUBMIT ARRAY JOB (6 LOMO folds simultaneously):
#   sbatch scripts/hpc_submit_geneformer.sh
#
# CHECK JOB STATUS:
#   squeue -u $USER
#   sacct -j <JOBID> --format=JobID,State,ExitCode,Elapsed,MaxRSS
#
# INTERACTIVE DEBUG SESSION (A40, 1h):
#   srun -p scu-gpu --gres=gpu:a40:1 --mem=32G --time=01:00:00 --pty bash
# ─────────────────────────────────────────────────────────────────────────────

#SBATCH --job-name=mouse_gf_lomo
#SBATCH --partition=scu-gpu
#SBATCH --output=/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/logs/geneformer_%A_%a.log
#SBATCH --error=/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/logs/geneformer_%A_%a.err
#SBATCH --time=04:00:00                   # 4h per fold; increase to 08:00:00 if needed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G                         # A40/A100 nodes: 1024GB total — 32G is safe
#SBATCH --gres=gpu:a40:1                 # 1x A40 48GB. Swap to gpu:a100:1 for A100.
#SBATCH --array=0-5                      # 6 LOMO folds (A1 liver: RR-1,3,6,8,9,MHU-2)

# ── Configuration ─────────────────────────────────────────────────────────────
# A1 liver LOMO: 6 missions, all used as held-out in turn.
FOLDS=("RR-1" "RR-3" "RR-6" "RR-8" "RR-9" "MHU-2")
FOLD=${FOLDS[$SLURM_ARRAY_TASK_ID]}
TASK="A1"
TASK_DIR="A1_liver_lomo"
MODEL_VERSION="mouse_gf"
EPOCHS=10
BATCH_SIZE=16        # A40 48GB: 16 comfortable; Mouse-Geneformer 6L is smaller than human
LR="2e-5"
FREEZE_LAYERS=4      # Freeze bottom 4/6 BERT layers (top-2 + head trainable)
                     # Recommended for small-n (n≈30-60): limits overfitting
                     # Use 0 for full fine-tuning if n>200; use 6 for head-only
SEED=42

# ── Cayuga / MasonLab account details ─────────────────────────────────────────
PROJECT_DIR="/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark"
CONDA_ENV="mouse_gf"
# Mouse-Geneformer has no setup.py → use PYTHONPATH (installed in home dir via git clone)
MOUSE_GF_SRC="/home/fs01/jak4013/Mouse-Geneformer"
# ─────────────────────────────────────────────────────────────────────────────

echo "======================================"
echo "Mouse-Geneformer LOMO Fine-Tuning — Cayuga"
echo "Job ID: ${SLURM_JOB_ID} | Array: ${SLURM_ARRAY_TASK_ID}"
echo "Fold: ${FOLD} | Task: ${TASK} | Model: Geneformer-${MODEL_VERSION}"
echo "Epochs: ${EPOCHS} | Batch: ${BATCH_SIZE} | LR: ${LR} | freeze_layers: ${FREEZE_LAYERS}"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader 2>/dev/null \
    || echo "nvidia-smi unavailable"
echo "======================================"

# Activate conda via ~/.bashrc (conda init writes to bashrc on this cluster)
source ~/.bashrc
conda activate "${CONDA_ENV}" \
    || { echo "ERROR: conda env '${CONDA_ENV}' not found"; exit 1; }

# Mouse-Geneformer: no setup.py, add to PYTHONPATH
export PYTHONPATH="${MOUSE_GF_SRC}:${PYTHONPATH}"

cd "${PROJECT_DIR}" || { echo "ERROR: project dir not found: ${PROJECT_DIR}"; exit 1; }
mkdir -p "${PROJECT_DIR}/logs"

echo "Working dir: $(pwd)"
echo "Python: $(which python) — $(python --version)"
echo "PyTorch: $(python -c 'import torch; print(torch.__version__, "CUDA:", torch.cuda.is_available())')"

# Step 1: Tokenize fold if not already done
# (Safe to re-run; will skip if tokenized data already exists)
echo ""
echo "[Step 1] Tokenizing fold ${FOLD} (will skip if already done)..."
python scripts/geneformer_tokenize.py \
    --task "${TASK}" \
    --task-dir "${TASK_DIR}" \
    --fold "${FOLD}" \
    --model-version "${MODEL_VERSION}" \
    2>&1 | tee "logs/tokenize_${FOLD}.log"

# Step 2: Fine-tune on GPU
echo ""
echo "[Step 2] Fine-tuning fold ${FOLD} on GPU..."
python scripts/geneformer_finetune.py \
    --task "${TASK}" \
    --task-dir "${TASK_DIR}" \
    --fold "${FOLD}" \
    --model-version "${MODEL_VERSION}" \
    --device cuda \
    --epochs "${EPOCHS}" \
    --batch-size "${BATCH_SIZE}" \
    --lr "${LR}" \
    --freeze-layers "${FREEZE_LAYERS}" \
    --seed "${SEED}" \
    2>&1 | tee "logs/finetune_${FOLD}.log"

EXIT_CODE=$?
echo ""
echo "======================================"
echo "Done: Fold ${FOLD} — exit code: ${EXIT_CODE}"
echo "======================================"
exit ${EXIT_CODE}

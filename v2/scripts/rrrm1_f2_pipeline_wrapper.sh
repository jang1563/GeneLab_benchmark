#!/bin/bash
#
# rrrm1_f2_pipeline_wrapper.sh
# Submit the full F2 downstream pipeline as a SLURM dependency chain.
#
# Prerequisites:
#   - Per-SRX STARsolo job 2703854 is running (or already complete)
#   - rrrm1_merge_per_srx.py, rrrm1_initial_scanpy.py, rrrm1_broad_annotate.py,
#     rrrm1_singlecell_hardening.py, rrrm1_f2a_composition.py,
#     rrrm1_f2b_pseudobulk_fgsea.py, rrrm1_f2c_loao_classifier.py
#     are all present in SCRATCH_DIR
#
# Usage:
#   bash rrrm1_f2_pipeline_wrapper.sh [STARSOLO_JOB_ID]
#   # Default STARSOLO_JOB_ID = 2703854 (the running per-SRX job)
#
# Submit via:
#   /opt/ohpc/pub/software/slurm/24.05.2/bin/sbatch rrrm1_f2_pipeline_wrapper.sh

set -euo pipefail

SBATCH=/opt/ohpc/pub/software/slurm/24.05.2/bin/sbatch
SQUEUE=/opt/ohpc/pub/software/slurm/24.05.2/bin/squeue

STARSOLO_JOB="${1:-2703854}"
SCRATCH="/athena/masonlab/scratch/users/jak4013/rrrm1_scrna"
HOME_DIR="/home/fs01/jak4013/rrrm1_scrna"
CONDA_ACTIVATE="source /home/fs01/jak4013/miniconda3/miniconda3/etc/profile.d/conda.sh && conda activate scgpt_env_new"

echo "=== F2 Pipeline Wrapper ==="
echo "Waiting for per-SRX STARsolo job: ${STARSOLO_JOB}"
echo "Date: $(date)"

# ── Step 1: Merge per-SRX STARsolo outputs → labeled h5ads ─────────────────
MERGE_SCRIPT="${SCRATCH}/rrrm1_merge_per_srx_slurm.sh"
cat > "${MERGE_SCRIPT}" << 'INNEREOF'
#!/bin/bash
#SBATCH --job-name=rrrm1_merge
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/athena/masonlab/scratch/users/jak4013/rrrm1_scrna/logs/merge_%j.out
#SBATCH --error=/athena/masonlab/scratch/users/jak4013/rrrm1_scrna/logs/merge_%j.err
set -euo pipefail
source /home/fs01/jak4013/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate scgpt_env_new
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
echo "=== Merging per-SRX outputs ==="
date
python3 /athena/masonlab/scratch/users/jak4013/rrrm1_scrna/rrrm1_merge_per_srx.py \
    --all \
    --map_csv /home/fs01/jak4013/rrrm1_scrna/RRRM1_SRX_CONDITION_MAP.csv
echo "=== Merge done ===" && date
INNEREOF

JID_MERGE=$(${SBATCH} --dependency=afterok:${STARSOLO_JOB} "${MERGE_SCRIPT}" | awk '{print $NF}')
echo "Step 1 (merge) submitted: job ${JID_MERGE}"


# ── Step 2: Initial scanpy processing (array 0-3, each tissue sequential) ──
PROCESS_SCRIPT="${SCRATCH}/rrrm1_process_labeled_slurm.sh"
cat > "${PROCESS_SCRIPT}" << 'INNEREOF'
#!/bin/bash
#SBATCH --job-name=rrrm1_process
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=3:00:00
#SBATCH --output=/athena/masonlab/scratch/users/jak4013/rrrm1_scrna/logs/process_%A_%a.out
#SBATCH --error=/athena/masonlab/scratch/users/jak4013/rrrm1_scrna/logs/process_%A_%a.err
#SBATCH --array=0-3
set -euo pipefail
source /home/fs01/jak4013/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate scgpt_env_new
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

TISSUES=(blood eye muscle skin)
TISSUE="${TISSUES[$SLURM_ARRAY_TASK_ID]}"
echo "=== Processing ${TISSUE} ==="
date

# Step 2a: initial scanpy (QC + HVG + PCA + UMAP, now with .raw before normalization)
python3 /athena/masonlab/scratch/users/jak4013/rrrm1_scrna/rrrm1_initial_scanpy.py \
    --tissue "${TISSUE}"

# Step 2b: broad annotation
python3 /athena/masonlab/scratch/users/jak4013/rrrm1_scrna/rrrm1_broad_annotate.py \
    --tissue "${TISSUE}"

# Step 2c: hardening (scrublet doublet detection)
python3 /athena/masonlab/scratch/users/jak4013/rrrm1_scrna/rrrm1_singlecell_hardening.py \
    --tissue "${TISSUE}"

echo "=== ${TISSUE} fully processed ===" && date
INNEREOF

JID_PROCESS=$(${SBATCH} --dependency=afterok:${JID_MERGE} "${PROCESS_SCRIPT}" | awk '{print $NF}')
echo "Step 2 (process) submitted: job ${JID_PROCESS}"


# ── Step 3: F2 analysis (A + B + C in parallel) ────────────────────────────
F2_SCRIPT="${SCRATCH}/rrrm1_f2_analysis_slurm.sh"
cat > "${F2_SCRIPT}" << 'INNEREOF'
#!/bin/bash
#SBATCH --job-name=rrrm1_f2
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:0
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=/athena/masonlab/scratch/users/jak4013/rrrm1_scrna/logs/f2_%A_%a.out
#SBATCH --error=/athena/masonlab/scratch/users/jak4013/rrrm1_scrna/logs/f2_%A_%a.err
#SBATCH --array=0-2
set -euo pipefail
source /home/fs01/jak4013/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate scgpt_env_new
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

TASKS=(f2a f2b f2c)
TASK="${TASKS[$SLURM_ARRAY_TASK_ID]}"

echo "=== Running ${TASK} ===" && date

case "${TASK}" in
  f2a)
    python3 /athena/masonlab/scratch/users/jak4013/rrrm1_scrna/rrrm1_f2a_composition.py --all
    ;;
  f2b)
    python3 /athena/masonlab/scratch/users/jak4013/rrrm1_scrna/rrrm1_f2b_pseudobulk_fgsea.py --all
    ;;
  f2c)
    python3 /athena/masonlab/scratch/users/jak4013/rrrm1_scrna/rrrm1_f2c_loao_classifier.py --all
    ;;
esac

echo "=== ${TASK} done ===" && date
INNEREOF

JID_F2=$(${SBATCH} --dependency=afterok:${JID_PROCESS} "${F2_SCRIPT}" | awk '{print $NF}')
echo "Step 3 (F2 analysis) submitted: job ${JID_F2}"


echo ""
echo "=== Pipeline chain submitted ==="
echo "  STARsolo (prerequisite): ${STARSOLO_JOB}"
echo "  Step 1 (merge):          ${JID_MERGE}  [afterok:${STARSOLO_JOB}]"
echo "  Step 2 (process ×4):     ${JID_PROCESS} [afterok:${JID_MERGE}]"
echo "  Step 3 (F2-A/B/C ×3):   ${JID_F2}     [afterok:${JID_PROCESS}]"
echo ""
echo "Monitor: ${SQUEUE} -u jak4013"
echo "Results: /home/fs01/jak4013/rrrm1_scrna/evaluation/"
echo "Figures: /home/fs01/jak4013/rrrm1_scrna/figures/"

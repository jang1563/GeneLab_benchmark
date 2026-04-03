#!/bin/bash
# Master submit script for v7 analyses
# Run from GeneLab_benchmark root directory on Cayuga login node.
# Usage: bash v7/scripts/hpc_submit_v7.sh [--with-setup|--setup-only] [--gnn-only] [--scprint-only]

set -eo pipefail

SETUP_ONLY=0
WITH_SETUP=0
GNN_ONLY=0
SCPRINT_ONLY=0
for arg in "$@"; do
    case $arg in
        --setup-only)  SETUP_ONLY=1 ;;
        --with-setup)  WITH_SETUP=1 ;;
        --gnn-only)    GNN_ONLY=1 ;;
        --scprint-only) SCPRINT_ONLY=1 ;;
    esac
done

BASE_DIR="${GENELAB_ROOT:?Set GENELAB_ROOT}"
SLURM=/opt/ohpc/pub/software/slurm/24.05.2/bin/sbatch

echo "v7 HPC Submission"
echo "Base: $BASE_DIR"
echo ""

# ── Optional setup (install PyG + scPRINT-2) ─────────────────────────────────
if [ $SETUP_ONLY -eq 1 ] || [ $WITH_SETUP -eq 1 ]; then
    echo ">>> Step 1: Running setup..."
    bash "$BASE_DIR/v7/scripts/hpc_setup_v7.sh"
    echo ""
fi

if [ $SETUP_ONLY -eq 1 ]; then
    echo "Setup complete. Exiting (--setup-only)."
    exit 0
fi

# ── Step 2: Submit GNN jobs ───────────────────────────────────────────────────
if [ $SCPRINT_ONLY -eq 0 ]; then
    echo ">>> Step 2: Submitting GNN array jobs (18 tasks)..."
    GNN_JOB=$($SLURM "$BASE_DIR/v7/scripts/hpc_gnn.sh")
    GNN_ID=$(echo $GNN_JOB | grep -oP '(?<=job )\d+')
    echo "  GNN job ID: $GNN_ID"
    echo "  Monitor: squeue -j $GNN_ID"
fi

# ── Step 3: Submit scPRINT jobs ────────────────────────────────────────────────
if [ $GNN_ONLY -eq 0 ]; then
    echo ">>> Step 3: Submitting scPRINT array jobs (6 tasks)..."
    SCPRINT_JOB=$($SLURM "$BASE_DIR/v7/scripts/hpc_scprint.sh")
    SCPRINT_ID=$(echo $SCPRINT_JOB | grep -oP '(?<=job )\d+')
    echo "  scPRINT job ID: $SCPRINT_ID"
    echo "  Monitor: squeue -j $SCPRINT_ID"
fi

echo ""
echo "All jobs submitted!"
echo ""
echo "Tip:"
echo "  If setup was already run, the normal fast path is:"
echo "  bash v7/scripts/hpc_submit_v7.sh"
echo ""
echo "Check progress:"
echo "  squeue -u \$USER"
echo ""
echo "Expected outputs:"
echo "  v7/evaluation/GNN_{tissue}_{graph_type}.json  (18 files)"
echo "  v7/evaluation/SCPRINT2_{tissue}.json           (6 files)"
echo ""
echo "After completion, sync back and run:"
echo "  python v7/unified/fig_v7_methods.py"

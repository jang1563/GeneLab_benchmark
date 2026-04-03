#!/bin/bash
# GeneLabBench v4: Submit all evaluations
# Usage: bash v4/scripts/hpc_submit_all.sh

set -eo pipefail
cd ${GENELAB_ROOT:?Set GENELAB_ROOT to your project directory}

SBATCH=/opt/ohpc/pub/software/slurm/24.05.2/bin/sbatch

# Create log directory
mkdir -p v4/logs

echo "=== GeneLabBench v4: Submitting evaluations ==="
echo ""

# Phase 1A: Gene features (CPU methods)
echo "Submitting: 8 tissues × 7 CPU methods × gene features = 56 jobs"
JOB_GENE_CPU=$($SBATCH --parsable v4/scripts/hpc_multi_method_cpu.sh gene)
echo "  CPU gene: Job array ${JOB_GENE_CPU}"

# Phase 1B: Gene features (TabNet/GPU)
echo "Submitting: 8 tissues × 1 GPU method × gene features = 8 jobs"
JOB_GENE_GPU=$($SBATCH --parsable v4/scripts/hpc_multi_method_gpu.sh gene)
echo "  GPU gene: Job array ${JOB_GENE_GPU}"

echo ""
echo "Total: 64 gene-feature evaluations submitted"
echo ""
echo "Monitor: squeue -u \$USER"
echo "Results: ls v4/evaluation/M1_*.json"
echo ""

# Pathway features (submit only if pathway_scores exist)
if [ -d "processed/pathway_scores" ]; then
    echo "Pathway scores found. Submitting pathway evaluations..."

    JOB_HM_CPU=$($SBATCH --parsable v4/scripts/hpc_multi_method_cpu.sh pathway_hallmark)
    echo "  CPU hallmark: Job array ${JOB_HM_CPU}"

    JOB_KEGG_CPU=$($SBATCH --parsable v4/scripts/hpc_multi_method_cpu.sh pathway_kegg)
    echo "  CPU KEGG: Job array ${JOB_KEGG_CPU}"

    JOB_COMB_CPU=$($SBATCH --parsable v4/scripts/hpc_multi_method_cpu.sh combined)
    echo "  CPU combined: Job array ${JOB_COMB_CPU}"

    JOB_HM_GPU=$($SBATCH --parsable v4/scripts/hpc_multi_method_gpu.sh pathway_hallmark)
    JOB_KEGG_GPU=$($SBATCH --parsable v4/scripts/hpc_multi_method_gpu.sh pathway_kegg)
    JOB_COMB_GPU=$($SBATCH --parsable v4/scripts/hpc_multi_method_gpu.sh combined)
    echo "  GPU pathway: 3 arrays submitted"

    echo ""
    echo "Total with pathways: 256 evaluations"
else
    echo "No pathway_scores/ found. Run generate_gsva.py first for pathway features."
    echo "Gene-only evaluations: 64 total"
fi

#!/bin/bash
# ── GeneLab Benchmark: fGSEA regeneration on Cayuga HPC ──────────────────────
#
# Runs run_fgsea.R for one tissue (or launches all six in parallel) on
# cayuga-phobos. Used for the 2026-08-14 flight-positive regeneration of
# processed/fgsea/ — see processed/fgsea/README.md.
#
# Usage (on phobos):
#   bash scripts/hpc_run_fgsea.sh <tissue>   # one tissue
#   bash scripts/hpc_run_fgsea.sh --all      # all six, in parallel (nohup)
#
# Staging a run from a local checkout:
#   BASE=/athena/cayuga_0003/scratch/users/jak4013/genelab_fgsea
#   ssh cayuga-login1 "mkdir -p $BASE/{scripts,logs,processed/fgsea/summary}"
#   rsync -a --prune-empty-dirs --include='*/' \
#         --include='*differential_expression*' --exclude='*' \
#         data/mouse/ cayuga-login1:$BASE/data/mouse/
#   rsync -a processed/ensembl_symbol_map.csv cayuga-login1:$BASE/processed/
#   rsync -a processed/gene_sets/ cayuga-login1:$BASE/processed/gene_sets/
#   rsync -a scripts/run_fgsea.R scripts/hpc_run_fgsea.sh cayuga-login1:$BASE/scripts/
#   # then pull results back (per-mission CSVs only; rebuild summaries locally):
#   rsync -a --exclude='summary/' cayuga-login1:$BASE/processed/fgsea/ processed/fgsea/
#
# NOTE: the HPC conda env has msigdbr 7.5.1, which cannot serve the current
# collections. Gene sets are therefore loaded from the pre-exported
# processed/gene_sets/*.gmt via --gmt-dir, keeping them identical to a local
# msigdbr 25.1.1 run. Regenerate those GMTs locally if MSigDB is updated.
#
# NOTE: each per-tissue run writes processed/fgsea/summary/all_fgsea_*.csv
# covering only its own tissue, so parallel runs overwrite each other there.
# Rebuild the combined summaries from the per-mission CSVs after pulling back.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

BASE="${GENELAB_FGSEA_BASE:-/athena/cayuga_0003/scratch/users/jak4013/genelab_fgsea}"
CONDA_ENV="${GENELAB_CONDA_ENV:-$HOME/miniconda3/miniconda3/envs/seurat.v5.R.4.3.3.JupyterLab.kmp.ml}"
DBS="${GENELAB_FGSEA_DBS:-hallmark,kegg,reactome,c2cgp,c5bp,mitocarta}"
TISSUES=(liver kidney thymus gastrocnemius eye skin)

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <tissue|--all>" >&2
    exit 1
fi

cd "$BASE"
mkdir -p logs processed/fgsea/summary

if [[ "$1" == "--all" ]]; then
    for t in "${TISSUES[@]}"; do
        nohup bash "$0" "$t" > "logs/nohup_${t}.out" 2>&1 &
        echo "launched $t (pid $!)"
    done
    echo "ALL LAUNCHED — poll with: grep -c 'Saved:' logs/fgsea_*.log"
    exit 0
fi

TISSUE="$1"
export LD_LIBRARY_PATH="${CONDA_ENV}/lib:${LD_LIBRARY_PATH:-}"
# Six tissues share 96 cores; cap workers per job.
export BIOCPARALLEL_WORKER_NUMBER="${GENELAB_FGSEA_WORKERS:-14}"

"${CONDA_ENV}/bin/Rscript" scripts/run_fgsea.R \
    --tissue "$TISSUE" --all \
    --db "$DBS" \
    --gmt-dir processed/gene_sets \
    > "logs/fgsea_${TISSUE}.log" 2>&1
echo "TISSUE ${TISSUE} DONE"

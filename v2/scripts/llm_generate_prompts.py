#!/usr/bin/env python3
"""
llm_generate_prompts.py — GeneLab_benchmark v2.0: Tier 3 Prompt Generation

Generates per-sample prompts for LLM-based spaceflight detection.
Each prompt includes:
  - Top-50 PROTEIN-CODING genes by RAW variance (from LOMO training data)
  - Per-sample z-scores (already computed in task data)
  - Gene symbols only (no ENSMUSG IDs — LLMs need recognizable names)
  - Genes ordered by |z-score| descending (most extreme first)

Protein-coding filter excludes:
  - Gm##### (predicted/pseudogenes)
  - LOC##### (unannotated loci)
  - ###Rik (RIKEN cDNA clones)
  - Multi-mapped genes (containing |)
  - Numeric-prefix genes (e.g., 6820431F20Rik)
  - ncRNAs (Rn7s, Rn7sk, Rn18s, etc.)
  - Unmapped ENSMUSG IDs (no symbol available)

Output:
  v2/processed/llm_prompts/{task_id}/
    fold_{mission}_prompts.json  — per-sample prompts for each fold

Usage:
  python v2/scripts/llm_generate_prompts.py --task A4
  python v2/scripts/llm_generate_prompts.py --all
"""

import json
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

V2_DIR = Path(__file__).resolve().parent.parent
PROJECT_DIR = V2_DIR.parent
sys.path.insert(0, str(PROJECT_DIR / "scripts"))

from utils import load_gene_features

# ── Config ─────────────────────────────────────────────────────────────────────
N_GENES = 50
TASKS_DIR = PROJECT_DIR / "tasks"
PROMPT_DIR = V2_DIR / "processed" / "llm_prompts"

# Task → tissue mapping
TASK_TISSUE = {
    "A1": "liver",
    "A2": "gastrocnemius",
    "A3": "kidney",
    "A4": "thymus",
    "A5": "skin",
    "A6": "eye",
}

TASK_DIR_MAP = {
    "A1": "A1_liver_lomo",
    "A2": "A2_gastrocnemius_lomo",
    "A3": "A3_kidney_lomo",
    "A4": "A4_thymus_lomo",
    "A5": "A5_skin_lomo",
    "A6": "A6_eye_lomo",
}

# ── Gene symbol mapping ──────────────────────────────────────────────────────

def load_symbol_map():
    """Load ENSEMBL → symbol mapping."""
    path = PROJECT_DIR / "processed" / "ensembl_symbol_map.csv"
    df = pd.read_csv(path)
    return dict(zip(df["ENSEMBL"], df["SYMBOL"]))


# ── Prompt templates ─────────────────────────────────────────────────────────

SYSTEM_PROMPT = (
    "You are a bioinformatics expert specializing in spaceflight transcriptomics "
    "and mouse gene expression analysis.\n\n"
    "Key biological effects of spaceflight/microgravity on mouse tissues include:\n"
    "- Oxidative stress and DNA damage response (upregulation of Nrf2 pathway, DNA repair genes)\n"
    "- Immune dysregulation (altered T-cell function, inflammatory cytokine changes)\n"
    "- Metabolic reprogramming (lipid metabolism changes, mitochondrial dysfunction)\n"
    "- Muscle and bone atrophy-related gene expression changes\n"
    "- Liver: xenobiotic metabolism changes, bile acid pathway alterations\n"
    "- Kidney: fluid homeostasis and electrolyte balance gene changes\n"
    "- Thymus: T-cell development and selection pathway changes\n\n"
    "Use your knowledge of these spaceflight-associated transcriptomic signatures "
    "to classify samples."
)


def make_user_prompt(tissue, gene_entries, n_genes=N_GENES):
    """Generate user prompt with per-sample z-scores.

    gene_entries: list of (symbol, zscore) tuples, already sorted by |zscore| desc.
    """
    gene_lines = []
    for symbol, zscore in gene_entries:
        direction = "UP" if zscore > 0 else "DOWN"
        gene_lines.append(f"  {symbol}: z = {zscore:+.2f} ({direction})")

    genes_formatted = "\n".join(gene_lines)

    return (
        f"Below are the {n_genes} most variable genes from a mouse {tissue} RNA-seq sample, "
        f"ordered by absolute z-score (most extreme first). "
        f"Z-scores are relative to a cross-study training cohort of ground control and flight mice.\n\n"
        f"{genes_formatted}\n\n"
        f"Based on the expression pattern above and your knowledge of spaceflight biology, "
        f"classify this sample:\n"
        f"(A) Flight — from mice aboard the International Space Station\n"
        f"(B) Ground — from ground control mice under normal gravity\n\n"
        f"Answer with the letter (A or B) followed by your confidence (0.0 to 1.0).\n"
        f"Format: \"A 0.82\" or \"B 0.65\""
    )


# ── Protein-coding gene filter ─────────────────────────────────────────────

import re

# ncRNA and pseudogene patterns
_EXCLUDE_PATTERNS = re.compile(
    r'^('
    # ncRNAs
    r'Rn7s[lk]?|Rn18s|Rn28s|Rn45s|Rn4\.5s|Rmrp|Rpph1|Terc|Malat1|Neat1|Xist|Tsix'
    r'|Snor[a-z]\d|Snhg\d|Mir\d|Scarna\d|Snord\d|Bc1$'
    # Predicted genes
    r'|Gm\d'
    # RNase P RNA-like
    r'|Rprl\d'
    r')',
    re.IGNORECASE
)

# Processed pseudogene suffix (e.g., Rps2-ps10, Ubb-ps, Gvin-ps4)
_PSEUDOGENE_SUFFIX = re.compile(r'-ps\d*$')


def is_llm_interpretable(symbol):
    """Check if a gene symbol is interpretable by an LLM.

    Returns True for protein-coding genes with recognizable symbols.
    Returns False for predicted genes, ncRNAs, pseudogenes, RIKEN clones, etc.
    """
    if not symbol or not isinstance(symbol, str):
        return False

    s = symbol.strip()
    if not s:
        return False

    # ENSMUSG IDs (unmapped)
    if s.startswith("ENSMUSG"):
        return False

    # LOC##### unannotated loci
    if s.startswith("LOC"):
        return False

    # RIKEN cDNA clones (e.g., 9030619P08Rik, A930018M24Rik)
    if s.endswith("Rik"):
        return False

    # Multi-mapped genes (e.g., Gm22634|Gm23804|Gm25679)
    if "|" in s:
        return False

    # Numeric-prefix genes (often RIKEN-like)
    if s[0].isdigit():
        return False

    # Processed pseudogenes (e.g., Rps2-ps10, Rpl15-ps6, Ubb-ps, Gvin-ps4)
    if _PSEUDOGENE_SUFFIX.search(s):
        return False

    # ncRNAs, predicted genes, and technical artifacts
    if _EXCLUDE_PATTERNS.match(s):
        return False

    return True


# ── Gene selection ──────────────────────────────────────────────────────────

def select_top_genes_by_raw_variance(tissue, train_sample_ids, selected_genes,
                                      symbol_map, n_genes=N_GENES):
    """Select top-N protein-coding genes by variance on RAW expression data.

    Filters out pseudogenes, ncRNAs, RIKEN clones, and unmapped genes
    BEFORE variance ranking, ensuring all selected genes are interpretable
    by an LLM.

    Args:
        tissue: tissue name for loading raw features
        train_sample_ids: sample IDs from the training fold
        selected_genes: genes that passed variance filter (from selected_genes.txt)
        symbol_map: ENSEMBL → symbol dict
        n_genes: number of top genes to select

    Returns:
        List of top-N gene IDs sorted by raw variance (descending),
        all with LLM-interpretable symbols
    """
    raw_feat = load_gene_features(tissue)

    # Match training samples
    common_samples = [s for s in train_sample_ids if s in raw_feat.index]
    if len(common_samples) < 5:
        # Fallback: try matching without mission prefix
        raw_idx_map = {idx.split(".")[-1] if "." in idx else idx: idx
                       for idx in raw_feat.index}
        common_samples = []
        for s in train_sample_ids:
            s_bare = s.split(".")[-1] if "." in s else s
            if s_bare in raw_idx_map:
                common_samples.append(raw_idx_map[s_bare])
            elif s in raw_feat.index:
                common_samples.append(s)

    if len(common_samples) < 5:
        print(f"      WARNING: Only {len(common_samples)} train samples matched in raw data")
        # Fallback: filter selected_genes for interpretable ones
        interpretable = [g for g in selected_genes
                         if is_llm_interpretable(symbol_map.get(g, g))]
        return interpretable[:n_genes]

    # Filter to selected genes present in raw data AND interpretable
    interpretable_genes = []
    n_filtered = 0
    for g in selected_genes:
        if g not in raw_feat.columns:
            continue
        sym = symbol_map.get(g, g)
        if is_llm_interpretable(sym):
            interpretable_genes.append(g)
        else:
            n_filtered += 1

    if len(interpretable_genes) < n_genes:
        print(f"      WARNING: Only {len(interpretable_genes)} interpretable genes "
              f"({n_filtered} filtered out)")
        return interpretable_genes

    # Compute variance on RAW training data (interpretable genes only)
    raw_train = raw_feat.loc[common_samples, interpretable_genes]
    gene_var = raw_train.var(axis=0)

    # Sort by variance descending, take top N
    top_genes = gene_var.sort_values(ascending=False).head(n_genes).index.tolist()
    return top_genes


# ── Core ─────────────────────────────────────────────────────────────────────

def generate_prompts_for_task(task_id):
    """Generate prompts for all folds of a task."""
    tissue = TASK_TISSUE[task_id]
    task_dir_name = TASK_DIR_MAP[task_id]
    task_dir = TASKS_DIR / task_dir_name

    if not task_dir.exists():
        print(f"  [SKIP] Task directory not found: {task_dir}")
        return None

    symbol_map = load_symbol_map()
    output_dir = PROMPT_DIR / task_id
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all fold directories
    folds = sorted([d for d in task_dir.iterdir()
                    if d.is_dir() and d.name.startswith("fold_") and "holdout" not in d.name])

    print(f"\n  Task {task_id} ({tissue}): {len(folds)} folds")

    total_prompts = 0
    task_meta = {
        "task_id": task_id,
        "tissue": tissue,
        "n_genes": N_GENES,
        "gene_selection": "top-50 protein-coding genes by raw variance",
        "gene_filter": "exclude Gm/LOC/Rik/ncRNA/multi-mapped/unmapped",
        "gene_ordering": "|z-score| descending per sample",
        "generated_at": datetime.now().isoformat(),
        "folds": {},
    }

    for fold_dir in folds:
        fold_name = fold_dir.name
        print(f"    {fold_name}...")

        # Load z-scored train and test data
        train_X = pd.read_csv(fold_dir / "train_X.csv", index_col=0)
        test_X = pd.read_csv(fold_dir / "test_X.csv", index_col=0)
        test_y = pd.read_csv(fold_dir / "test_y.csv", index_col=0)

        # Load all genes that passed variance filter
        selected_genes_path = fold_dir / "selected_genes.txt"
        selected_genes = selected_genes_path.read_text().strip().split("\n")

        # Select top-N protein-coding genes by RAW variance
        top_genes = select_top_genes_by_raw_variance(
            tissue, train_X.index.tolist(), selected_genes,
            symbol_map, n_genes=N_GENES
        )

        # Ensure genes are available in test data
        available = [g for g in top_genes if g in test_X.columns]
        if len(available) < N_GENES:
            # Fill with additional interpretable genes from raw-variance ranking
            all_raw_ranked = select_top_genes_by_raw_variance(
                tissue, train_X.index.tolist(), selected_genes,
                symbol_map, n_genes=N_GENES * 3
            )
            for g in all_raw_ranked:
                if g not in available and g in test_X.columns:
                    available.append(g)
                if len(available) >= N_GENES:
                    break

        # All genes should now have symbols (filter guarantees it)
        n_mapped = sum(1 for g in available
                       if is_llm_interpretable(symbol_map.get(g, g)))

        # Generate per-sample prompts
        fold_prompts = []
        for sample_id in test_X.index:
            # Use z-scores directly (data is already z-scored by generate_tasks.py)
            sample_zscores = test_X.loc[sample_id, available]

            gene_entries = []
            for gene_id in available:
                z = sample_zscores[gene_id]
                if not np.isfinite(z):
                    continue
                # Use symbol if available, otherwise ENSMUSG ID
                symbol = symbol_map.get(gene_id, gene_id)
                gene_entries.append((symbol, float(z)))

            # Sort by |z-score| descending (most extreme first)
            gene_entries.sort(key=lambda x: abs(x[1]), reverse=True)

            user_prompt = make_user_prompt(tissue, gene_entries, n_genes=len(gene_entries))

            # True label
            true_label = int(test_y.loc[sample_id].values[0]) if sample_id in test_y.index else None

            fold_prompts.append({
                "sample_id": sample_id,
                "true_label": true_label,
                "system_prompt": SYSTEM_PROMPT,
                "user_prompt": user_prompt,
            })

        total_prompts += len(fold_prompts)

        # Save fold prompts
        fold_output = {
            "fold": fold_name,
            "task_id": task_id,
            "tissue": tissue,
            "n_genes": len(available),
            "n_genes_with_symbol": n_mapped,
            "n_test_samples": len(fold_prompts),
            "gene_selection_method": "protein_coding_raw_variance_top50",
            "gene_ordering": "abs_zscore_descending",
            "prompts": fold_prompts,
        }

        fold_path = output_dir / f"{fold_name}_prompts.json"
        fold_path.write_text(json.dumps(fold_output, indent=2))

        task_meta["folds"][fold_name] = {
            "n_test_samples": len(fold_prompts),
            "n_genes": len(available),
            "n_genes_with_symbol": n_mapped,
        }

        print(f"      {len(fold_prompts)} prompts, {n_mapped}/{len(available)} genes mapped")

    # Save task metadata
    task_meta["total_prompts"] = total_prompts
    meta_path = output_dir / "task_meta.json"
    meta_path.write_text(json.dumps(task_meta, indent=2))

    print(f"  Total: {total_prompts} prompts → {output_dir}")
    return task_meta


def parse_args():
    parser = argparse.ArgumentParser(description="Generate LLM prompts for Tier 3")
    parser.add_argument("--task", type=str, choices=list(TASK_TISSUE.keys()),
                        help="Task ID (e.g., A4)")
    parser.add_argument("--all", action="store_true", help="Generate for all tasks")
    return parser.parse_args()


def main():
    args = parse_args()

    if not args.task and not args.all:
        print("Specify --task A1/A2/.../A6 or --all")
        return

    print("=" * 70)
    print("Tier 3: LLM Prompt Generation (v2 — protein-coding, raw variance)")
    print(f"  Top {N_GENES} protein-coding genes per prompt (by raw variance)")
    print(f"  Filter: exclude Gm/LOC/Rik/ncRNA/multi-mapped/unmapped")
    print(f"  Gene ordering: |z-score| descending per sample")
    print("=" * 70)

    tasks = list(TASK_TISSUE.keys()) if args.all else [args.task]
    all_meta = {}

    for task in tasks:
        meta = generate_prompts_for_task(task)
        if meta:
            all_meta[task] = meta

    # Summary
    total = sum(m["total_prompts"] for m in all_meta.values())
    print(f"\n  Grand total: {total} prompts across {len(all_meta)} tasks")

    # Save overall summary
    summary_path = PROMPT_DIR / "generation_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps({
        "generated_at": datetime.now().isoformat(),
        "n_genes_per_prompt": N_GENES,
        "gene_selection": "protein_coding_raw_variance_top50",
        "gene_filter": "exclude_Gm_LOC_Rik_ncRNA_multimapped_unmapped",
        "gene_ordering": "abs_zscore_descending",
        "total_prompts": total,
        "tasks": {k: v["total_prompts"] for k, v in all_meta.items()},
    }, indent=2))
    print(f"  Summary: {summary_path}")


if __name__ == "__main__":
    main()

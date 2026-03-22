#!/usr/bin/env python3
"""
generate_pathway_scores.py — GeneLabBench v4: Generate pathway scores via ssGSEA

Python fallback for R GSVA. Uses gseapy ssgsea to compute per-sample pathway
activity scores for all 8 tissues × {hallmark, kegg}.

Input:  processed/A_detection/{tissue}/{tissue}_*_log2_norm.csv
Output: processed/pathway_scores/{tissue}/{mission}_gsva_{db}.csv

Usage:
  python generate_pathway_scores.py --tissue liver --all
  python generate_pathway_scores.py --all-tissues
  python generate_pathway_scores.py --tissue lung --db hallmark
"""

import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

warnings.filterwarnings("ignore")

import sys
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from v4_utils import TISSUE_MISSIONS, BASE_DIR, load_metadata

PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
PATHWAY_DIR = BASE_DIR / "processed" / "pathway_scores"
SYMBOL_MAP_FILE = BASE_DIR / "processed" / "ensembl_symbol_map.csv"


def load_ensembl_to_symbol():
    """Load ENSMUSG → gene symbol mapping."""
    if SYMBOL_MAP_FILE.exists():
        df = pd.read_csv(SYMBOL_MAP_FILE)
        return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))

    # Fallback: extract from any log2_norm file columns
    # (gene names starting with ENSMUSG don't have symbols — need mapping)
    print("[WARN] No ensembl_symbol_map.csv found. Will try to use ENSMUSG IDs directly.")
    return {}


def get_gene_sets(db="hallmark", species="Mus musculus"):
    """Get gene sets from msigdbr-compatible gseapy sources."""
    import gseapy as gp

    if db == "hallmark":
        # MSigDB Hallmark via gseapy
        try:
            gs = gp.get_library("MSigDB_Hallmark_2020", organism="Mouse")
            print(f"  Loaded MSigDB_Hallmark_2020: {len(gs)} gene sets")
            return gs
        except Exception:
            gs = gp.get_library("MSigDB_Hallmark_2020", organism="Human")
            print(f"  Loaded MSigDB_Hallmark_2020 (Human): {len(gs)} gene sets")
            return gs

    elif db == "kegg":
        try:
            gs = gp.get_library("KEGG_2019_Mouse", organism="Mouse")
            print(f"  Loaded KEGG_2019_Mouse: {len(gs)} gene sets")
            return gs
        except Exception:
            gs = gp.get_library("KEGG_2021_Human", organism="Human")
            print(f"  Loaded KEGG_2021_Human: {len(gs)} gene sets")
            return gs

    else:
        raise ValueError(f"Unknown db: {db}. Use 'hallmark' or 'kegg'.")


def compute_ssgsea(expr_df, gene_sets, min_size=15, max_size=500):
    """Compute ssGSEA scores for each sample.

    Args:
        expr_df: genes × samples DataFrame (gene symbols as index)
        gene_sets: dict of {pathway_name: [gene_list]}

    Returns:
        samples × pathways DataFrame
    """
    import gseapy as gp

    # gseapy ssgsea expects genes × samples
    if expr_df.shape[0] < expr_df.shape[1]:
        expr_df = expr_df.T  # samples × genes → genes × samples

    ss = gp.ssgsea(
        data=expr_df,
        gene_sets=gene_sets,
        outdir=None,
        min_size=min_size,
        max_size=max_size,
        sample_norm_method="rank",
        permutation_num=0,
        no_plot=True,
        threads=4,
        verbose=False,
    )

    # Result: res2d has columns [Name, Term, ES, NES]
    # Name = sample name, Term = pathway name
    if ss.res2d.empty:
        raise ValueError("ssGSEA returned empty results — check gene name overlap")

    # Pivot to samples × pathways
    scores = ss.res2d.pivot(index="Name", columns="Term", values="NES")
    scores.index.name = None
    scores.columns.name = None

    return scores


def process_tissue_mission(tissue, mission, dbs, symbol_map):
    """Process one tissue-mission pair."""
    # Load expression data
    expr_file = PROCESSED_DIR / tissue / f"{tissue}_all_missions_log2_norm.csv"
    if not expr_file.exists():
        expr_file = PROCESSED_DIR / tissue / f"{tissue}_{mission}_log2_norm.csv"
    if not expr_file.exists():
        print(f"  [SKIP] No expression file for {tissue}/{mission}")
        return

    expr = pd.read_csv(expr_file, index_col=0)

    # Orient: ensure genes × samples (index = genes)
    if expr.shape[0] < expr.shape[1]:
        # More columns than rows → likely samples × genes, transpose
        expr = expr.T

    # Convert to numeric (CSV may load as object type)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr.dropna(how="all")

    # Map ENSMUSG to symbols if needed
    if expr.index[0].startswith("ENSMUSG"):
        if symbol_map:
            expr.index = [symbol_map.get(g, g) for g in expr.index]
            # Drop genes without symbol mapping
            expr = expr[~expr.index.str.startswith("ENSMUSG")]
            # Drop duplicates
            expr = expr[~expr.index.duplicated(keep="first")]
            print(f"  Mapped ENSMUSG → symbols: {len(expr)} genes")
        else:
            print(f"  [WARN] ENSMUSG IDs without symbol map — pathway scoring may fail")

    # Uppercase gene names for matching with MSigDB (human uppercase convention)
    expr.index = expr.index.str.upper()

    # Filter to mission samples if multi-mission file
    meta = load_metadata(tissue)
    mission_samples = meta[meta["mission"] == mission].index.tolist()

    # Try to find matching columns
    available_cols = set(expr.columns)
    matching = [s for s in mission_samples if s in available_cols]

    if not matching:
        # Try without mission prefix
        for s in mission_samples:
            parts = str(s).split(".", 1)
            stripped = parts[1] if len(parts) == 2 else s
            if stripped in available_cols:
                matching.append(stripped)

    if len(matching) < 3:
        print(f"  [SKIP] {tissue}/{mission}: only {len(matching)} matching samples")
        return

    expr_mission = expr[matching]
    print(f"  {tissue}/{mission}: {expr_mission.shape[0]} genes × {expr_mission.shape[1]} samples")

    for db in dbs:
        out_file = PATHWAY_DIR / tissue / f"{mission}_gsva_{db}.csv"
        if out_file.exists():
            print(f"  [EXISTS] {out_file.name}")
            continue

        try:
            gene_sets = get_gene_sets(db)
            scores = compute_ssgsea(expr_mission, gene_sets)
            out_file.parent.mkdir(parents=True, exist_ok=True)
            scores.to_csv(out_file)
            print(f"  [SAVED] {out_file.name}: {scores.shape[0]} samples × {scores.shape[1]} pathways")
        except Exception as e:
            print(f"  [ERROR] {tissue}/{mission}/{db}: {e}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tissue", default=None)
    parser.add_argument("--all-tissues", action="store_true")
    parser.add_argument("--db", default="hallmark,kegg")
    parser.add_argument("--all", action="store_true",
                        help="All missions for given tissue")
    args = parser.parse_args()

    dbs = [d.strip() for d in args.db.split(",")]
    symbol_map = load_ensembl_to_symbol()
    print(f"Symbol map: {len(symbol_map)} entries")

    if args.all_tissues:
        tissues = list(TISSUE_MISSIONS.keys())
    elif args.tissue:
        tissues = [args.tissue]
    else:
        parser.error("Specify --tissue or --all-tissues")

    for tissue in tissues:
        print(f"\n{'='*60}")
        print(f"Tissue: {tissue}")
        print(f"{'='*60}")
        missions = TISSUE_MISSIONS[tissue]
        for mission in missions:
            process_tissue_mission(tissue, mission, dbs, symbol_map)


if __name__ == "__main__":
    main()

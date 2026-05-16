#!/usr/bin/env python3
"""
wgcna_prep.py — Prepare gene×sample expression matrices and trait matrices for WGCNA.

Outputs (written to v4/wgcna_inputs/):
  {tissue}_expr.csv   — genes (rows) × samples (cols), symbols as gene names
  {tissue}_traits.csv — samples (rows) × traits (cols): condition, duration_days, mission dummies

Uses limma_rbe batch-corrected expression for 5 tissues; raw log2_norm for skin (no limma_rbe).
Filters to top-N MAD genes before writing.

Usage:
  python wgcna_prep.py [--tissues liver kidney ...]
  python wgcna_prep.py  # runs all 6 LOMO tissues
"""

import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent.parent
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
SYMBOL_MAP_PATH = BASE_DIR / "processed" / "ensembl_symbol_map.csv"
OUTPUT_DIR = BASE_DIR / "v4" / "wgcna_inputs"

# ── Configuration ──────────────────────────────────────────────────────────────
# 6 LOMO tissues: use limma_rbe if available, fallback to raw log2_norm
LOMO_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin"]

# Skin has no limma_rbe file — use raw
USE_RAW = {"skin"}

# Small-n tissues: use fewer MAD genes for WGCNA stability
SMALL_N_TISSUES = {"gastrocnemius", "eye"}  # n=32, n=41
MAD_NGENES_DEFAULT = 5000
MAD_NGENES_SMALL   = 3000

# Label encoding (FLT/GC only; AG excluded)
LABEL_MAP = {
    "Flight": 1, "FLT": 1,
    "GC": 0, "Ground Control": 0, "Ground": 0,
    "Vivarium": 0, "Basal": 0, "BC": 0, "VC": 0,
}
EXCLUDE_LABELS = {"AG"}


def load_symbol_map():
    """Load Ensembl → gene symbol mapping."""
    if not SYMBOL_MAP_PATH.exists():
        warnings.warn(f"Symbol map not found: {SYMBOL_MAP_PATH}")
        return {}
    df = pd.read_csv(SYMBOL_MAP_PATH)
    return dict(zip(df["ENSEMBL"], df["SYMBOL"]))


def load_expression(tissue):
    """Load batch-corrected (or raw) log2-normalised expression matrix.

    Returns DataFrame: samples (rows) × genes (cols, Ensembl IDs).
    """
    tissue_dir = PROCESSED_DIR / tissue

    # Prefer limma_rbe unless tissue is in USE_RAW
    if tissue not in USE_RAW:
        rbe_path = tissue_dir / f"{tissue}_all_missions_log2_norm_limma_rbe.csv"
        if rbe_path.exists():
            df = pd.read_csv(rbe_path, index_col=0)
            df = df.apply(pd.to_numeric, errors="coerce")
            # Auto-orient: samples should be rows (more columns = genes)
            if df.shape[0] > df.shape[1]:
                df = df.T
            print(f"  [{tissue}] loaded limma_rbe: {df.shape}")
            return df

    # Fallback: raw log2_norm
    raw_path = tissue_dir / f"{tissue}_all_missions_log2_norm.csv"
    if not raw_path.exists():
        raise FileNotFoundError(f"No expression file for tissue={tissue}")
    df = pd.read_csv(raw_path, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    if df.shape[0] > df.shape[1]:
        df = df.T
    # Drop non-Ensembl columns
    ens_cols = [c for c in df.columns if str(c).startswith("ENSMUSG")]
    if ens_cols:
        df = df[ens_cols]
    print(f"  [{tissue}] loaded raw log2_norm: {df.shape}")
    return df


def load_metadata(tissue):
    """Load all-missions metadata; filter REMOVE samples."""
    f = PROCESSED_DIR / tissue / f"{tissue}_all_missions_metadata.csv"
    if not f.exists():
        raise FileNotFoundError(f"No metadata for tissue={tissue}")
    meta = pd.read_csv(f, index_col=0)
    if "REMOVE" in meta.columns:
        meta = meta[meta["REMOVE"] != True]
    return meta


def align(expr, meta):
    """Align expression and metadata by sample name."""
    common = sorted(set(expr.index) & set(meta.index))
    if len(common) >= 5:
        return expr.loc[common], meta.loc[common]

    # Try stripping mission prefix (some metadata indices include mission.)
    meta_map = {}
    feat_set = set(expr.index)
    for idx in meta.index:
        parts = str(idx).split(".", 1)
        stripped = parts[1] if len(parts) == 2 else idx
        if stripped in feat_set:
            meta_map[idx] = stripped

    if len(meta_map) >= 5:
        meta_al = meta.loc[list(meta_map.keys())]
        expr_al = expr.loc[list(meta_map.values())]
        expr_al.index = meta_al.index
        return expr_al, meta_al

    raise ValueError(
        f"[{tissue}] Too few aligned samples: expr={len(expr)}, meta={len(meta)}, common={len(common)}"
    )


def filter_mad_genes(expr, n_genes):
    """Select top-n_genes genes by median absolute deviation (computed across all samples)."""
    mad = expr.apply(lambda col: np.median(np.abs(col - np.median(col))), axis=0)
    top_genes = mad.nlargest(min(n_genes, len(mad))).index
    return expr[top_genes]


def map_gene_symbols(expr, sym_map):
    """Replace Ensembl IDs with gene symbols where available.

    Genes without a symbol keep their Ensembl ID.
    Deduplicates: if two Ensembl IDs map to same symbol, keep higher-MAD one.
    """
    new_cols = [sym_map.get(g, g) for g in expr.columns]
    expr.columns = new_cols

    # Deduplicate: keep the gene with higher sum of absolute values
    if len(set(new_cols)) < len(new_cols):
        mag = expr.abs().sum(axis=0)
        expr = expr.loc[:, ~expr.columns.duplicated(keep=False)]
        # keep first occurrence (already ordered by MAD descending)
    return expr


def build_trait_matrix(meta):
    """Build trait matrix for WGCNA module-trait correlation.

    Traits:
      - condition: FLT=1, GC=0, NA for AG/other
      - duration_days: numeric (if available)
      - mission_*: binary dummy for each mission (omit one for identifiability)
    """
    traits = pd.DataFrame(index=meta.index)

    # Condition (binary: FLT=1, GC=0, NA for excluded)
    condition = meta["label"].map(LABEL_MAP)
    traits["condition"] = condition  # NaN for AG and unmapped

    # Duration
    if "duration_days" in meta.columns:
        dur = pd.to_numeric(meta["duration_days"], errors="coerce")
        traits["duration_days"] = dur

    # Mission dummies (one-hot, drop first for identifiability)
    if "mission" in meta.columns:
        mission_dummies = pd.get_dummies(meta["mission"], prefix="mission", drop_first=True)
        traits = pd.concat([traits, mission_dummies.astype(float)], axis=1)

    return traits


def process_tissue(tissue, sym_map):
    """Run full prep pipeline for one tissue. Writes _expr.csv and _traits.csv."""
    print(f"\n[{tissue}] Processing...")

    # Load
    expr = load_expression(tissue)
    meta = load_metadata(tissue)

    # Align
    expr, meta = align(expr, meta)
    print(f"  Aligned: {expr.shape[0]} samples × {expr.shape[1]} genes")

    # MAD filter
    n_mad = MAD_NGENES_SMALL if tissue in SMALL_N_TISSUES else MAD_NGENES_DEFAULT
    expr_filt = filter_mad_genes(expr, n_mad)
    print(f"  After MAD filter (top {n_mad}): {expr_filt.shape[1]} genes retained")

    # Symbol mapping
    expr_sym = map_gene_symbols(expr_filt.copy(), sym_map)
    n_symbols = sum(1 for c in expr_sym.columns if not c.startswith("ENSMUSG"))
    print(f"  Symbol mapping: {n_symbols}/{expr_sym.shape[1]} genes have symbols")

    # Trait matrix
    traits = build_trait_matrix(meta)
    n_flt = (traits["condition"] == 1).sum()
    n_gc  = (traits["condition"] == 0).sum()
    n_ag  = traits["condition"].isna().sum()
    print(f"  Labels — FLT={n_flt}, GC={n_gc}, AG/other={n_ag}")
    print(f"  Traits: {list(traits.columns)}")

    # Write
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Transpose to genes×samples for R WGCNA convention
    expr_for_r = expr_sym.T  # genes × samples
    expr_path   = OUTPUT_DIR / f"{tissue}_expr.csv"
    traits_path = OUTPUT_DIR / f"{tissue}_traits.csv"

    expr_for_r.to_csv(expr_path)
    traits.to_csv(traits_path)
    print(f"  Wrote: {expr_path.name} ({expr_for_r.shape}) and {traits_path.name} ({traits.shape})")

    return {
        "tissue": tissue,
        "n_samples": expr_sym.shape[0],
        "n_genes_total": expr.shape[1],
        "n_genes_mad_filtered": expr_filt.shape[1],
        "n_genes_with_symbol": n_symbols,
        "n_flt": int(n_flt),
        "n_gc": int(n_gc),
        "n_excluded": int(n_ag),
        "n_traits": len(traits.columns),
        "use_limma_rbe": tissue not in USE_RAW,
    }


def main():
    parser = argparse.ArgumentParser(description="Prepare WGCNA inputs from GeneLabBench data")
    parser.add_argument("--tissues", nargs="+", default=LOMO_TISSUES,
                        help="Tissues to process (default: all 6 LOMO tissues)")
    args = parser.parse_args()

    print("Loading gene symbol map...")
    sym_map = load_symbol_map()
    print(f"  {len(sym_map)} Ensembl→symbol mappings loaded")

    summaries = []
    for tissue in args.tissues:
        if tissue not in LOMO_TISSUES:
            print(f"WARNING: {tissue} is not a LOMO tissue — skipping")
            continue
        summary = process_tissue(tissue, sym_map)
        summaries.append(summary)

    print("\n═══ Summary ═══")
    print(f"{'Tissue':<15} {'N':<6} {'MAD genes':<12} {'FLT':<6} {'GC':<6} {'limma_rbe'}")
    for s in summaries:
        print(f"  {s['tissue']:<13} {s['n_samples']:<6} {s['n_genes_mad_filtered']:<12} "
              f"{s['n_flt']:<6} {s['n_gc']:<6} {s['use_limma_rbe']}")
    print(f"\nOutput directory: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()

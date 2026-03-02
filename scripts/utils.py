#!/usr/bin/env python3
"""
utils.py — GeneLab_benchmark: Shared utility functions

Common data loading and alignment functions used across benchmark scripts.
"""

import numpy as np
import pandas as pd
from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
PATHWAY_DIR = BASE_DIR / "processed" / "pathway_scores"

# ── Config ─────────────────────────────────────────────────────────────────────
TISSUE_MISSIONS = {
    "liver": ["RR-1", "RR-3", "RR-6", "RR-8", "RR-9", "MHU-2"],
    "gastrocnemius": ["RR-1", "RR-5", "RR-9"],
    "kidney": ["RR-1", "RR-3", "RR-7"],
    "thymus": ["RR-6", "MHU-1", "MHU-2", "RR-9"],
    "eye": ["RR-1", "RR-3", "TBD"],
}


# ── Data Loading ───────────────────────────────────────────────────────────────

def load_metadata(tissue):
    """Load all-missions metadata for a tissue."""
    f = PROCESSED_DIR / tissue / f"{tissue}_all_missions_metadata.csv"
    meta = pd.read_csv(f, index_col=0)
    if "REMOVE" in meta.columns:
        meta = meta[meta["REMOVE"] != True]
    return meta


def load_gene_features(tissue):
    """Load gene-level log2 normalized counts (samples x genes)."""
    f = PROCESSED_DIR / tissue / f"{tissue}_all_missions_log2_norm.csv"
    df = pd.read_csv(f, index_col=0)
    if df.shape[0] > df.shape[1]:
        df = df.T
    meta_cols = [c for c in df.columns if not str(c).startswith("ENSMUSG")]
    if meta_cols:
        df = df.drop(columns=meta_cols)
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def load_pathway_features(tissue, db="hallmark"):
    """Load GSVA pathway scores across all missions for a tissue."""
    all_scores = []
    for mission in TISSUE_MISSIONS.get(tissue, []):
        f = PATHWAY_DIR / tissue / f"{mission}_gsva_{db}.csv"
        if not f.exists():
            continue
        scores = pd.read_csv(f, index_col=0)
        all_scores.append(scores)
    if not all_scores:
        return None
    combined = pd.concat(all_scores)
    if "mission" in combined.columns:
        combined = combined.drop(columns=["mission"])
    return combined


def align_features_with_meta(features, meta):
    """Align feature matrix with metadata by sample name."""
    feat_set = set(features.index)
    meta_set = set(meta.index)

    common = sorted(feat_set & meta_set)
    if len(common) >= 5:
        return features.loc[common], meta.loc[common]

    # Try stripping mission prefix from metadata index
    meta_map = {}
    for idx in meta.index:
        parts = str(idx).split(".", 1)
        stripped = parts[1] if len(parts) == 2 else idx
        if stripped in feat_set:
            meta_map[idx] = stripped

    if len(meta_map) >= 5:
        meta_aligned = meta.loc[list(meta_map.keys())]
        feat_aligned = features.loc[list(meta_map.values())]
        feat_aligned.index = meta_aligned.index
        return feat_aligned, meta_aligned

    raise ValueError(
        f"Too few aligned samples: features={len(feat_set)}, "
        f"meta={len(meta_set)}, common={len(common)}"
    )

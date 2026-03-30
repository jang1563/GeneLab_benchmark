"""
v6_utils.py — Shared utilities for GeneLabBench v6: Mouse→Human Translation

Functions for loading ortholog maps, cfRNA data, SHAP results, TF activity,
drug targets, and biomarker panel from both GeneLab_benchmark and SpaceOmicsBench.
"""

import os
import json
import numpy as np
import pandas as pd
from pathlib import Path

# ── Base paths ──────────────────────────────────────────────────────────────
GENELAB_BASE = Path(os.path.expanduser(
    "~/Dropbox/Bioinformatics/Claude/GeneLab_benchmark"))
SPACEOMICS_BASE = Path(os.path.expanduser(
    "~/Dropbox/Bioinformatics/Claude/SpaceOmicsBench"))

V4_EVAL = GENELAB_BASE / "v4" / "evaluation"
V5_EVAL = GENELAB_BASE / "v5" / "evaluation"
V6_EVAL = GENELAB_BASE / "v6" / "evaluation"
PROCESSED = GENELAB_BASE / "processed"

# SpaceOmicsBench data paths (v2.1 is latest with metabolomics)
SOB_PROCESSED = SPACEOMICS_BASE / "v2_public" / "data" / "processed"
SOB_FEATURE_MATRIX = (SPACEOMICS_BASE /
    "2025_01_08_v2.1.metabolomics" / "SpaceOmicsBench_v2.1" / "data" /
    "cfrna_feature_matrix.csv")

# All 8 mouse tissues
ALL_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus",
               "eye", "skin", "lung", "colon"]

# ── Ortholog mapping ────────────────────────────────────────────────────────
def load_ortholog_map():
    """Load mouse_symbol → human_symbol mapping.
    Returns dict with ~17K pairs.
    """
    path = GENELAB_BASE / "v3" / "data" / "mouse_human_orthologs.tsv"
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df["mouse_symbol"], df["human_symbol"]))


def load_ensembl_to_symbol():
    """Load ENSMUSG → mouse_symbol mapping.
    Returns dict with ~27K entries.
    """
    path = PROCESSED / "ensembl_symbol_map.csv"
    df = pd.read_csv(path)
    return dict(zip(df["ENSEMBL"], df["SYMBOL"]))


def map_mouse_to_human(gene_list, orth_map=None, ens2sym=None):
    """Map a list of mixed ENSMUSG/symbol mouse genes to human symbols.

    Returns:
        mapped: dict {mouse_gene: human_symbol}
        unmapped: list of genes with no ortholog
    """
    if orth_map is None:
        orth_map = load_ortholog_map()
    if ens2sym is None:
        ens2sym = load_ensembl_to_symbol()

    mapped = {}
    unmapped = []

    for gene in gene_list:
        # If ENSMUSG, convert to symbol first
        symbol = ens2sym.get(gene, gene) if gene.startswith("ENSMUSG") else gene

        # Try direct mapping
        human = orth_map.get(symbol)
        if human is None:
            # Try case variations
            human = orth_map.get(symbol.capitalize())
        if human is None:
            human = orth_map.get(symbol.upper())

        if human:
            mapped[gene] = human
        else:
            unmapped.append(gene)

    return mapped, unmapped


# ── cfRNA data loading ──────────────────────────────────────────────────────
def load_cfrna_de():
    """Load full cfRNA differential expression (26,845 genes).
    Returns DataFrame with gene as index.
    """
    path = SOB_PROCESSED / "cfrna_3group_de.csv"
    df = pd.read_csv(path)
    df = df.set_index("gene")
    return df


def load_cfrna_drr():
    """Load 466 DRR (differentially regulated response) gene names.
    Returns set of gene names.
    """
    path = SOB_PROCESSED / "cfrna_466drr.csv"
    df = pd.read_csv(path)
    return set(df["gene"].values)


def load_cfrna_feature_matrix():
    """Load cfRNA feature matrix (24 samples × 500 genes).
    Returns: X (DataFrame, samples×genes), meta (DataFrame with phase/day).
    """
    df = pd.read_csv(SOB_FEATURE_MATRIX)
    meta_cols = ["sample_id", "mission", "crew_id", "phase", "day"]
    meta = df[meta_cols].copy()
    X = df.drop(columns=meta_cols)
    return X, meta


# ── SHAP data loading ──────────────────────────────────────────────────────
def load_shap_consensus():
    """Load SHAP consensus genes (5 genes across ≥3 tissues, ≥2 methods).
    Returns dict with consensus_genes, tissue_specific_genes, etc.
    """
    path = V4_EVAL / "SHAP_consensus.json"
    with open(path) as f:
        return json.load(f)


def load_shap_top_genes(tissue, top_n=100):
    """Load SHAP top genes for a tissue (union across available methods).
    Returns set of gene symbols.
    """
    genes = set()
    # Find all SHAP files for this tissue (excluding interactions)
    shap_files = sorted(V4_EVAL.glob(f"SHAP_{tissue}_*.json"))
    shap_files = [f for f in shap_files if "interaction" not in f.name
                  and "WGCNA" not in f.name and "consensus" not in f.name]

    for path in shap_files:
        with open(path) as f:
            data = json.load(f)
        # Primary format: top_100_genes as list of strings
        if "top_100_genes" in data:
            for g in data["top_100_genes"][:top_n]:
                genes.add(g)
        elif "top_genes" in data:
            for g in data["top_genes"][:top_n]:
                if isinstance(g, dict):
                    genes.add(g.get("gene", g.get("name", "")))
                else:
                    genes.add(g)
    return genes


def load_all_shap_top_genes(top_n=100):
    """Load SHAP top genes for all 8 tissues.
    Returns dict {tissue: set of gene symbols}.
    """
    result = {}
    for tissue in ALL_TISSUES:
        result[tissue] = load_shap_top_genes(tissue, top_n)
    return result


# ── WGCNA data loading ─────────────────────────────────────────────────────
def load_wgcna_hub_genes():
    """Load WGCNA hub genes (566 across 6 tissues).
    Returns dict structure from JSON.
    """
    path = V4_EVAL / "WGCNA_hub_genes.json"
    with open(path) as f:
        return json.load(f)


# ── TF activity loading ────────────────────────────────────────────────────
def load_tf_activity(tissue):
    """Load TF activity results for a tissue.
    Returns dict with tf_results containing per-TF statistics.
    """
    path = V5_EVAL / f"tf_activity_{tissue}.json"
    with open(path) as f:
        return json.load(f)


# ── Drug target loading ────────────────────────────────────────────────────
def load_drug_targets():
    """Load v5 drug target mapping results.
    Returns dict with tiers, per-tissue druggability, etc.
    """
    path = V5_EVAL / "drug_targets.json"
    with open(path) as f:
        return json.load(f)


# ── Biomarker panel loading ────────────────────────────────────────────────
def load_biomarker_panel():
    """Load v5 consensus biomarker panel (20 genes).
    Returns dict with panel list, each gene having mouse + human_symbol.
    """
    path = V5_EVAL / "consensus_biomarker_panel.json"
    with open(path) as f:
        return json.load(f)


# ── Pathway scores loading ─────────────────────────────────────────────────
def load_pathway_scores_hallmark(tissue):
    """Load ssGSEA Hallmark pathway scores for a tissue.
    Returns DataFrame (samples × pathways).
    """
    pw_dir = PROCESSED / "pathway_scores" / tissue
    files = list(pw_dir.glob("*hallmark*.csv")) + list(pw_dir.glob("*Hallmark*.csv"))
    if not files:
        return None

    dfs = []
    for f in files:
        df = pd.read_csv(f, index_col=0)
        # Ensure numeric
        df = df.apply(pd.to_numeric, errors="coerce")
        dfs.append(df)

    if len(dfs) == 1:
        return dfs[0]

    # Concat and intersect columns
    common_cols = set(dfs[0].columns)
    for df in dfs[1:]:
        common_cols &= set(df.columns)
    common_cols = sorted(common_cols)
    return pd.concat([df[common_cols] for df in dfs], axis=0)


# ── Metadata loading ───────────────────────────────────────────────────────
LABEL_MAP = {
    "Flight": 1, "FLT": 1,
    "Ground Control": 0, "GC": 0, "Ground": 0,
    "Basal": 0, "BC": 0, "Basal Control": 0,
    "Vivarium": 0, "VC": 0, "Vivarium Control": 0,
}
EXCLUDE_LABELS = {"AG", "Artificial Gravity"}


def load_tissue_metadata(tissue):
    """Load metadata for a mouse tissue with encoded labels.
    Returns DataFrame with 'label' column (0/1) and 'mission' column.
    """
    import sys
    sys.path.insert(0, str(GENELAB_BASE / "v4" / "scripts"))
    from v4_utils import load_metadata
    meta = load_metadata(tissue)
    return meta


# ── JSON output helper ─────────────────────────────────────────────────────
class SafeEncoder(json.JSONEncoder):
    """Handle numpy and other non-standard types for JSON serialization."""
    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, set):
            return sorted(obj)
        return super().default(obj)


def save_json(data, path):
    """Save dict to JSON with numpy-safe encoding."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(data, f, indent=2, cls=SafeEncoder)
    print(f"  Saved: {path}")

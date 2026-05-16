#!/usr/bin/env python3
"""E4/E5: Multi-species KEGG NES concordance.

Compares spaceflight pathway-level responses across species:
  Mouse (v1 6 tissues) ↔ Drosophila (GLDS-207) ↔ Arabidopsis (GLDS-37, GLDS-120)

Strategy:
  1. Per species: Welch t-test ranking (FLT vs GC) on normalized counts
  2. fGSEA with Enrichr KEGG libraries (species-specific)
  3. Cross-species matching via normalized pathway names
  4. Spearman r with bootstrap CI + permutation p

Notes:
  - Mouse + Drosophila: Both have Enrichr KEGG libraries (KEGG_2019_Mouse,
    KEGG_2019_Drosophila) with matching pathway names.
  - Arabidopsis: No Enrichr KEGG available. WikiPathways attempted as fallback.
    Limited cross-species matching expected (a finding, not a bug).
  - v1 fGSEA used KEGG_MEDICUS (incompatible names) → fresh Enrichr KEGG run.

Usage:
  python e4_multispecies_nes.py [--skip-fgsea]
"""
import argparse
import gzip
import json
import re
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Paths ──────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent  # v3/
V1_BASE = BASE_DIR.parent  # GeneLab_benchmark/
EVAL_DIR = BASE_DIR / "evaluation"
FIG_DIR = BASE_DIR / "figures"
EVAL_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

MULTISPECIES_DIR = V1_BASE / "data" / "multispecies"
MOUSE_DATA_DIR = V1_BASE / "data" / "mouse"
ENSEMBL_MAP = V1_BASE / "processed" / "ensembl_symbol_map.csv"

# Mouse tissue → representative mission (largest n, most-used in v1)
# Each tissue has its own GLDS dataset number
MOUSE_TISSUES = {
    "liver":         [("RR-1", "GLDS-48"), ("RR-3", "GLDS-137")],
    "thymus":        [("MHU-2", "GLDS-289"), ("RR-6", "GLDS-244")],
    "kidney":        [("RR-1", "GLDS-102"), ("RR-3", "GLDS-163")],
    "eye":           [("RR-1", "GLDS-100")],
    "skin":          [("RR-7", "GLDS-254"), ("RR-6", "GLDS-243")],
    "gastrocnemius": [("RR-1", "GLDS-101"), ("RR-9", "GLDS-326")],
}

RANDOM_SEED = 42


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def save_json(data: dict, path: Path):
    with open(path, "w") as f:
        json.dump(data, f, indent=2, cls=NumpyEncoder)
    print(f"Saved: {path}")


# ── Shared fGSEA column normalizer ────────────────────────────────────
def normalize_fgsea_columns(result_df: pd.DataFrame) -> pd.DataFrame:
    """Normalize gseapy res2d column names to: pathway, NES, pval, fdr."""
    col_map = {}
    for col in result_df.columns:
        lc = col.lower()
        if lc == "nes":
            col_map[col] = "NES"
        elif "pval" in lc or "nom p" in lc:
            col_map[col] = "pval"
        elif "fdr" in lc:
            col_map[col] = "fdr"
    result_df = result_df.rename(columns=col_map)

    # Identify pathway column
    # gseapy res2d: "Name" = library/run name (e.g. "prerank"), "Term" = actual pathway
    # Must prioritize "Term" over "Name" to get pathway names, not "prerank"
    pw_col = None
    for priority in ["term", "gene_set"]:
        for c in result_df.columns:
            if c.lower() == priority:
                pw_col = c
                break
        if pw_col is not None:
            break
    if pw_col is None:
        # Fallback: skip "Name" (it's "prerank"), use first non-stat column
        for c in result_df.columns:
            if c.lower() not in ("name", "nes", "pval", "fdr", "es",
                                  "nom p-val", "fdr q-val", "fwer p-val",
                                  "tag %", "gene %", "lead_genes"):
                pw_col = c
                break
    if pw_col is None:
        pw_col = result_df.columns[0]
    if pw_col != "pathway":
        result_df = result_df.rename(columns={pw_col: "pathway"})

    return result_df


# ── Gene Ranking ──────────────────────────────────────────────────────
def rank_genes_ttest(expr_df: pd.DataFrame, meta: pd.DataFrame,
                     label: str = "") -> pd.Series:
    """Welch t-test ranking (FLT vs GC). Returns t-statistics per gene."""
    common = sorted(set(expr_df.columns) & set(meta.index))
    if len(common) < 4:
        raise ValueError(f"Too few common samples: {len(common)}")

    expr = expr_df[common]
    meta_a = meta.loc[common]

    flt_samples = meta_a[meta_a["condition"] == "FLT"].index.tolist()
    gc_samples = meta_a[meta_a["condition"] == "GC"].index.tolist()
    print(f"  {label}Ranking: {len(flt_samples)} FLT vs {len(gc_samples)} GC")

    flt_arr = expr[flt_samples].values.astype(float).T  # samples × genes
    gc_arr = expr[gc_samples].values.astype(float).T

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        t_stats, _ = stats.ttest_ind(flt_arr, gc_arr, equal_var=False, axis=0)

    t_stats = np.where(np.isfinite(t_stats), t_stats, 0.0)
    ranking = pd.Series(t_stats, index=expr.index)
    ranking = ranking[ranking != 0]
    return ranking


def run_gseapy_prerank(ranking: pd.Series, gene_set_lib: str,
                       label: str = "") -> pd.DataFrame:
    """Run gseapy preranked fGSEA. Returns normalized DataFrame."""
    import gseapy as gp

    ranking_clean = ranking.copy()
    ranking_clean = ranking_clean[~ranking_clean.index.duplicated(keep="first")]
    ranking_clean = ranking_clean.dropna()

    if len(ranking_clean) < 100:
        print(f"  {label}Too few ranked genes ({len(ranking_clean)})")
        return pd.DataFrame()

    print(f"  {label}fGSEA: {gene_set_lib} ({len(ranking_clean)} genes)...")
    try:
        res = gp.prerank(
            rnk=ranking_clean,
            gene_sets=gene_set_lib,
            outdir=None,
            no_plot=True,
            min_size=10,
            max_size=500,
            permutation_num=1000,
            seed=RANDOM_SEED,
            verbose=False,
        )
        result_df = normalize_fgsea_columns(res.res2d)
        print(f"  {label}Got {len(result_df)} pathways")
        return result_df
    except Exception as e:
        print(f"  {label}fGSEA error: {e}")
        return pd.DataFrame()


# ── 1. Mouse ──────────────────────────────────────────────────────────
def load_ensembl_symbol_map() -> dict:
    """Load ENSMUSG → uppercase gene symbol mapping."""
    if not ENSEMBL_MAP.exists():
        print(f"  ENSMUSG→symbol map not found: {ENSEMBL_MAP}")
        return {}
    df = pd.read_csv(ENSEMBL_MAP)
    # Expected columns: ensembl_gene_id, mgi_symbol (or similar)
    id_col = None
    sym_col = None
    for c in df.columns:
        cl = c.lower()
        if "ensembl" in cl or "ensmusg" in cl:
            id_col = c
        elif "symbol" in cl or "mgi" in cl or "gene_name" in cl:
            sym_col = c
    if id_col is None or sym_col is None:
        # Fallback: assume first two columns
        id_col = df.columns[0]
        sym_col = df.columns[1]
    return dict(zip(df[id_col].astype(str), df[sym_col].astype(str)))


def run_mouse_kegg(tissue: str) -> dict:
    """Run fresh Enrichr KEGG fGSEA for one mouse tissue.

    Loads normalized counts from v1 data directory.
    Returns: {pathway: NES} or empty dict.
    """
    tissue_dir = MOUSE_DATA_DIR / tissue
    if not tissue_dir.exists():
        print(f"  Mouse {tissue}: data dir not found")
        return {}

    missions = MOUSE_TISSUES.get(tissue, [])
    if not missions:
        return {}

    # Use first available mission
    for mission, glds in missions:
        mission_dir = tissue_dir / mission
        if not mission_dir.exists():
            continue

        # Find normalized counts and sample table (with or without _GLbulkRNAseq suffix)
        count_files = list(mission_dir.glob(f"{glds}_rna_seq_Normalized_Counts*.csv"))
        # Prefer _GLbulkRNAseq version; exclude _rRNArm_ and _ERCC_ variants
        count_files = [f for f in count_files if "rRNArm" not in f.name and "ERCC" not in f.name]
        sample_files = list(mission_dir.glob(f"{glds}_rna_seq_SampleTable*.csv"))

        if not count_files or not sample_files:
            continue

        print(f"  Mouse {tissue}/{mission} ({glds})...")

        # Load expression (genes × samples, index_col=0)
        expr = pd.read_csv(count_files[0], index_col=0)
        print(f"    Expression: {expr.shape[0]} genes × {expr.shape[1]} samples")

        # Load sample table (index_col=0 for unnamed first column)
        st = pd.read_csv(sample_files[0], index_col=0)

        # Parse condition: find FLT vs GC
        # GeneLab sample tables use "condition" column
        cond_col = None
        for c in st.columns:
            if c.lower() in ("condition", "factor value[spaceflight]"):
                cond_col = c
                break

        if cond_col is None:
            # Try to infer from column values
            for c in st.columns:
                vals = st[c].astype(str).str.lower().unique()
                if any("flight" in v or "space" in v for v in vals):
                    cond_col = c
                    break

        if cond_col is None:
            print(f"    No condition column found in {sample_files[0].name}")
            continue

        # Map condition to FLT/GC
        meta = pd.DataFrame(index=st.index)
        cond_vals = st[cond_col].astype(str)
        meta["condition"] = "UNKNOWN"
        meta.loc[cond_vals.str.contains("Flight|Space|FLT", case=False, na=False), "condition"] = "FLT"
        meta.loc[cond_vals.str.contains("Ground|Control|GC|Vivarium|Basal", case=False, na=False), "condition"] = "GC"

        n_flt = (meta["condition"] == "FLT").sum()
        n_gc = (meta["condition"] == "GC").sum()
        if n_flt < 2 or n_gc < 2:
            print(f"    Insufficient samples: FLT={n_flt}, GC={n_gc}")
            continue

        # Filter out non-gene rows (metadata contamination)
        gene_mask = expr.index.astype(str).str.startswith("ENSMUSG")
        if gene_mask.sum() > 0:
            expr = expr[gene_mask]

        # Rank genes
        ranking = rank_genes_ttest(expr, meta, label=f"{tissue}: ")

        # Convert ENSMUSG → gene symbols (uppercase for Enrichr KEGG matching)
        if len(ranking) > 0 and str(ranking.index[0]).startswith("ENSMUSG"):
            sym_map = load_ensembl_symbol_map()
            if sym_map:
                ranking.index = ranking.index.map(
                    lambda x: sym_map.get(x, x))
                ranking = ranking[~ranking.index.astype(str).str.startswith("ENSMUSG")]
                ranking = ranking[~ranking.index.duplicated(keep="first")]
                print(f"    Mapped to {len(ranking)} gene symbols")

        # Uppercase for Enrichr matching
        ranking.index = ranking.index.astype(str).str.upper()
        ranking = ranking[~ranking.index.duplicated(keep="first")]

        # Run KEGG fGSEA
        fgsea_df = run_gseapy_prerank(ranking, "KEGG_2019_Mouse",
                                       label=f"{tissue}: ")

        if not fgsea_df.empty and "NES" in fgsea_df.columns:
            return dict(zip(fgsea_df["pathway"], fgsea_df["NES"].astype(float)))

    return {}


# ── 2. Drosophila ─────────────────────────────────────────────────────
def download_flybase_symbols() -> dict:
    """Download FBgn → gene symbol mapping from FlyBase.

    Uses fbgn_annotation_ID.tsv.gz which has columns:
      gene_symbol, organism_abbreviation, primary_FBgn, ...
    The gene symbol is in column 0, FBgn in column 2.
    """
    local_map = MULTISPECIES_DIR / "drosophila" / "fbgn_to_symbol.tsv"
    if local_map.exists():
        df = pd.read_csv(local_map, sep="\t", header=None, names=["fbgn", "symbol"])
        return dict(zip(df["fbgn"], df["symbol"]))

    print("  Downloading FlyBase gene symbol mapping...")
    url = "https://ftp.flybase.net/releases/current/precomputed_files/genes/fbgn_annotation_ID.tsv.gz"
    try:
        response = urllib.request.urlopen(url, timeout=60)
        data = gzip.decompress(response.read())
        lines = data.decode().strip().split("\n")

        fbgn_map = {}
        for line in lines:
            if line.startswith("#") or line.startswith("##"):
                continue
            parts = line.split("\t")
            # Format: gene_symbol \t organism_abbrev \t primary_FBgn# \t ...
            if len(parts) >= 3:
                symbol = parts[0]
                fbgn = parts[2]  # primary FBgn#
                if fbgn.startswith("FBgn") and symbol:
                    fbgn_map[fbgn] = symbol

        # Save locally
        local_map.parent.mkdir(parents=True, exist_ok=True)
        with open(local_map, "w") as f:
            for fbgn, sym in fbgn_map.items():
                f.write(f"{fbgn}\t{sym}\n")

        print(f"  Downloaded {len(fbgn_map)} FBgn→symbol mappings")
        return fbgn_map

    except Exception as e:
        print(f"  FlyBase download failed: {e}")
        return {}


def parse_glds207_metadata(sample_table_path: Path) -> pd.DataFrame:
    """Parse GLDS-207 (Drosophila) sample table."""
    df = pd.read_csv(sample_table_path, index_col=0)
    records = []
    for name in df.index:
        cond_str = str(df.loc[name, "condition"]) if "condition" in df.columns else str(name)

        if "Space" in str(name) or "Space" in cond_str:
            condition = "FLT"
        elif "Ground" in str(name) or "Ground" in cond_str:
            condition = "GC"
        else:
            condition = "UNKNOWN"

        records.append({"sample_name": str(name), "condition": condition})

    return pd.DataFrame(records).set_index("sample_name")


def _biomart_query(query_xml: str, timeout: int = 120) -> str:
    """Execute a BioMart query and return raw TSV text."""
    url = "http://www.ensembl.org/biomart/martservice?query=" + urllib.request.quote(query_xml)
    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode()


def download_drosophila_human_orthologs() -> dict:
    """Download Drosophila → human ortholog mapping from Ensembl BioMart.

    Uses two separate 2-attribute queries (more reliable than 3-attribute):
      1. FBgn → gene_symbol
      2. gene_symbol → human_symbol
    Combines them: FBgn → symbol → human

    Returns: {fbgn_or_symbol_upper: human_symbol_upper}
    """
    local_path = MULTISPECIES_DIR / "drosophila" / "dmel_human_orthologs_v3.tsv"
    if local_path.exists():
        df = pd.read_csv(local_path, sep="\t")
        mapping = {}
        for _, row in df.iterrows():
            fbgn = str(row.iloc[0]).strip().upper()
            human = str(row.iloc[1]).strip().upper()
            if fbgn and human and human != "" and human != "NAN":
                mapping[fbgn] = human
        print(f"  Loaded {len(mapping)} FBgn→human orthologs (cached)")
        return mapping

    print("  Downloading Drosophila orthologs from Ensembl BioMart (2-step)...")

    # Step 1: FBgn → gene symbol
    q1 = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1">
  <Dataset name="dmelanogaster_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
  </Dataset>
</Query>"""

    # Step 2: gene symbol → human ortholog
    q2 = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1">
  <Dataset name="dmelanogaster_gene_ensembl" interface="default">
    <Attribute name="external_gene_name"/>
    <Attribute name="hsapiens_homolog_associated_gene_name"/>
  </Dataset>
</Query>"""

    try:
        # Query 1: FBgn → symbol
        data1 = _biomart_query(q1)
        fbgn_to_sym = {}
        for line in data1.strip().split("\n")[1:]:
            parts = line.split("\t")
            if len(parts) >= 2 and parts[0].strip() and parts[1].strip():
                fbgn_to_sym[parts[0].strip()] = parts[1].strip()
        print(f"  Step 1: {len(fbgn_to_sym)} FBgn→symbol mappings")

        # Query 2: symbol → human
        data2 = _biomart_query(q2)
        sym_to_human = {}
        for line in data2.strip().split("\n")[1:]:
            parts = line.split("\t")
            if len(parts) >= 2 and parts[0].strip() and parts[1].strip():
                sym = parts[0].strip().upper()
                human = parts[1].strip().upper()
                # Keep first mapping (one-to-many: pick first)
                if sym not in sym_to_human:
                    sym_to_human[sym] = human
        print(f"  Step 2: {len(sym_to_human)} symbol→human mappings")

        # Combine: FBgn → symbol → human
        mapping = {}
        rows_saved = []
        for fbgn, sym in fbgn_to_sym.items():
            human = sym_to_human.get(sym.upper())
            if human:
                mapping[fbgn.upper()] = human
                mapping[sym.upper()] = human
                rows_saved.append(f"{fbgn}\t{human}")

        # Cache FBgn→human (deduped)
        local_path.parent.mkdir(parents=True, exist_ok=True)
        seen = set()
        with open(local_path, "w") as f:
            f.write("fbgn\thuman_symbol\n")
            for fbgn, sym in fbgn_to_sym.items():
                human = sym_to_human.get(sym.upper())
                if human and fbgn not in seen:
                    f.write(f"{fbgn}\t{human}\n")
                    seen.add(fbgn)

        n_fbgn = sum(1 for k in mapping if k.startswith("FBGN"))
        print(f"  Combined: {n_fbgn} FBgn→human, {len(mapping)} total entries")
        return mapping

    except Exception as e:
        print(f"  BioMart download failed: {e}")
        return {}


def run_drosophila_kegg() -> dict:
    """Run Enrichr KEGG fGSEA for Drosophila (GLDS-207).

    Strategy: FBgn → Drosophila symbol → human ortholog → KEGG_2019_Mouse
    Returns: {pathway: NES} or empty dict.
    """
    dros_dir = MULTISPECIES_DIR / "drosophila"
    counts_path = dros_dir / "GLDS-207_rna_seq_Normalized_Counts_GLbulkRNAseq.csv"
    samples_path = dros_dir / "GLDS-207_rna_seq_SampleTable_GLbulkRNAseq.csv"

    if not counts_path.exists() or not samples_path.exists():
        print("  Drosophila data not found")
        return {}

    meta = parse_glds207_metadata(samples_path)
    expr = pd.read_csv(counts_path, index_col=0)
    print(f"  Expression: {expr.shape[0]} genes × {expr.shape[1]} samples")

    ranking = rank_genes_ttest(expr, meta, label="Dros: ")
    print(f"  Ranked {len(ranking)} genes (FBgn IDs)")

    # Map FBgn IDs → human orthologs via Ensembl BioMart (includes both FBgn and symbol keys)
    ortho_map = download_drosophila_human_orthologs()
    if ortho_map:
        ranking_upper = ranking.copy()
        ranking_upper.index = ranking_upper.index.astype(str).str.upper()
        ranking_upper = ranking_upper[~ranking_upper.index.duplicated(keep="first")]

        mapped_idx = ranking_upper.index.map(lambda x: ortho_map.get(x, ""))
        has_map = mapped_idx != ""
        ranking_human = ranking_upper[has_map].copy()
        ranking_human.index = mapped_idx[has_map]
        ranking_human = ranking_human[~ranking_human.index.duplicated(keep="first")]
        print(f"  Mapped {len(ranking_human)} genes to human orthologs "
              f"({len(ranking_human)/len(ranking_upper)*100:.1f}% coverage)")
        ranking = ranking_human
    else:
        # Fallback: try FBgn→symbol then uppercase
        fbgn_map = download_flybase_symbols()
        if fbgn_map:
            ranking.index = ranking.index.map(lambda x: fbgn_map.get(x, x))
            ranking = ranking[~ranking.index.astype(str).str.startswith("FBgn")]
            ranking = ranking[~ranking.index.duplicated(keep="first")]
        ranking.index = ranking.index.astype(str).str.upper()
        ranking = ranking[~ranking.index.duplicated(keep="first")]
        print(f"  No ortholog map — using uppercase symbols ({len(ranking)} genes)")

    fgsea_df = run_gseapy_prerank(ranking, "KEGG_2019_Mouse", label="Dros: ")

    if not fgsea_df.empty and "NES" in fgsea_df.columns:
        csv_out = EVAL_DIR / "E4_drosophila_GLDS207_kegg_nes.csv"
        fgsea_df.to_csv(csv_out, index=False)
        return dict(zip(fgsea_df["pathway"], fgsea_df["NES"].astype(float)))

    return {}


# ── 3. Arabidopsis ────────────────────────────────────────────────────
def parse_arabidopsis_metadata(sample_table_path: Path) -> pd.DataFrame:
    """Parse Arabidopsis sample table (GLDS-37 or GLDS-120)."""
    df = pd.read_csv(sample_table_path, index_col=0)
    records = []
    for name in df.index:
        cond_str = str(df.loc[name, "condition"]) if "condition" in df.columns else str(name)

        if "Space" in cond_str or "FLT" in str(name):
            condition = "FLT"
        elif "Ground" in cond_str or "GC" in str(name):
            condition = "GC"
        else:
            condition = "UNKNOWN"

        records.append({"sample_name": str(name), "condition": condition})

    return pd.DataFrame(records).set_index("sample_name")


def download_arabidopsis_human_orthologs() -> dict:
    """Download Arabidopsis → human ortholog mapping from Ensembl Plants BioMart.

    Returns: {AT_gene_upper: human_symbol_upper}
    Note: Cross-kingdom orthologs are very limited. This is expected.
    """
    local_path = MULTISPECIES_DIR / "arabidopsis" / "ath_human_orthologs.tsv"
    if local_path.exists():
        df = pd.read_csv(local_path, sep="\t")
        mapping = {}
        for _, row in df.iterrows():
            ath = str(row.iloc[0]).strip().upper()
            human = str(row.iloc[1]).strip().upper()
            if ath and human and human != "" and human != "NAN":
                mapping[ath] = human
        print(f"  Loaded {len(mapping)} Arabidopsis→human orthologs (cached)")
        return mapping

    # Ensembl Plants BioMart for Arabidopsis orthologs
    print("  Downloading Arabidopsis→human orthologs from Ensembl Plants BioMart...")
    query = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="plants_mart" formatter="TSV" header="1" uniqueRows="1">
  <Dataset name="athaliana_eg_gene" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="hsapiens_homolog_associated_gene_name"/>
  </Dataset>
</Query>"""

    url = "https://plants.ensembl.org/biomart/martservice?query=" + urllib.request.quote(query)
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = resp.read().decode()

        lines = data.strip().split("\n")
        mapping = {}
        header = True
        rows_saved = []
        for line in lines:
            if header:
                header = False
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                ath_id = parts[0].strip()
                human_sym = parts[1].strip()
                if ath_id and human_sym:
                    mapping[ath_id.upper()] = human_sym.upper()
                    rows_saved.append(f"{ath_id}\t{human_sym}")

        local_path.parent.mkdir(parents=True, exist_ok=True)
        with open(local_path, "w") as f:
            f.write("ath_gene_id\thuman_symbol\n")
            f.write("\n".join(rows_saved) + "\n")

        print(f"  Downloaded {len(mapping)} Arabidopsis→human orthologs")
        return mapping

    except Exception as e:
        print(f"  Plants BioMart download failed: {e}")
        return {}


def run_arabidopsis_pathway(glds_id: str, label: str) -> dict:
    """Run pathway fGSEA for Arabidopsis.

    Strategy: AT gene IDs → human orthologs → KEGG_2019_Mouse
    Cross-kingdom coverage expected to be limited (a finding, not a bug).
    Returns: {pathway: NES} or empty dict.
    """
    arab_dir = MULTISPECIES_DIR / "arabidopsis"
    counts_path = arab_dir / f"GLDS-{glds_id}_rna_seq_Normalized_Counts_GLbulkRNAseq.csv"
    samples_path = arab_dir / f"GLDS-{glds_id}_rna_seq_SampleTable_GLbulkRNAseq.csv"

    if not counts_path.exists() or not samples_path.exists():
        print(f"  Arabidopsis GLDS-{glds_id} data not found")
        return {}

    meta = parse_arabidopsis_metadata(samples_path)
    expr = pd.read_csv(counts_path, index_col=0)
    print(f"  Expression: {expr.shape[0]} genes × {expr.shape[1]} samples")

    ranking = rank_genes_ttest(expr, meta, label=f"Arab-{glds_id}: ")
    print(f"  Ranked {len(ranking)} genes (AT IDs)")

    # Map Arabidopsis gene IDs → human orthologs for KEGG matching
    ortho_map = download_arabidopsis_human_orthologs()
    if ortho_map:
        ranking_upper = ranking.copy()
        ranking_upper.index = ranking_upper.index.astype(str).str.upper()
        ranking_upper = ranking_upper[~ranking_upper.index.duplicated(keep="first")]

        mapped_idx = ranking_upper.index.map(lambda x: ortho_map.get(x, ""))
        has_map = mapped_idx != ""
        ranking_human = ranking_upper[has_map].copy()
        ranking_human.index = mapped_idx[has_map]
        ranking_human = ranking_human[~ranking_human.index.duplicated(keep="first")]
        print(f"  Mapped {len(ranking_human)} genes to human orthologs "
              f"({len(ranking_human)/len(ranking_upper)*100:.1f}% coverage)")

        if len(ranking_human) >= 100:
            fgsea_df = run_gseapy_prerank(ranking_human, "KEGG_2019_Mouse",
                                           label=f"Arab-{glds_id}: ")
            if not fgsea_df.empty and "NES" in fgsea_df.columns:
                csv_out = EVAL_DIR / f"E4_{label}_kegg_nes.csv"
                fgsea_df.to_csv(csv_out, index=False)
                return dict(zip(fgsea_df["pathway"], fgsea_df["NES"].astype(float)))
        else:
            print(f"  Too few mapped genes ({len(ranking_human)}) for fGSEA")

    print(f"  No pathway DB worked for Arabidopsis GLDS-{glds_id} (cross-kingdom limitation)")
    return {}


# ── Pathway Name Normalization ────────────────────────────────────────
def normalize_pathway_name(pw: str) -> str:
    """Normalize KEGG/pathway names for cross-species matching.

    Examples:
        'Oxidative phosphorylation Mus musculus (mouse)' → 'oxidative phosphorylation'
        'Oxidative phosphorylation Drosophila melanogaster (fly)' → 'oxidative phosphorylation'
        'dme00190 Oxidative phosphorylation' → 'oxidative phosphorylation'
    """
    pw = str(pw)
    # Remove species suffixes like "Mus musculus (mouse)"
    pw = re.sub(r'\s+(Mus musculus|Homo sapiens|Drosophila melanogaster|'
                r'Arabidopsis thaliana)\s*(\(.*?\))?\s*$', '', pw, flags=re.IGNORECASE)
    # Remove KEGG pathway IDs (hsa00010, dme00190, mmu00010, ath00010)
    pw = re.sub(r'^[a-z]{2,4}\d{5}\s+', '', pw)
    # Remove "KEGG_" prefix
    pw = re.sub(r'^KEGG[_\s]+', '', pw, flags=re.IGNORECASE)
    # Normalize
    pw = pw.strip().lower()
    pw = re.sub(r'\s+', ' ', pw)
    pw = pw.rstrip('.')
    # Remove common trailing noise
    pw = re.sub(r'\s*-\s*$', '', pw)
    return pw


# ── Cross-Species Concordance ─────────────────────────────────────────
def compute_concordance(nes_a: dict, nes_b: dict,
                        name_a: str, name_b: str,
                        min_pathways: int = 5) -> dict:
    """Compute Spearman concordance between two NES dicts."""
    norm_a = {normalize_pathway_name(pw): nes for pw, nes in nes_a.items()}
    norm_b = {normalize_pathway_name(pw): nes for pw, nes in nes_b.items()}

    common = sorted(set(norm_a.keys()) & set(norm_b.keys()))

    if len(common) < min_pathways:
        return {
            "species_a": name_a,
            "species_b": name_b,
            "n_pathways_a": len(nes_a),
            "n_pathways_b": len(nes_b),
            "n_common_pathways": len(common),
            "common_pathways": common,
            "status": "insufficient_common_pathways",
        }

    x = np.array([norm_a[pw] for pw in common])
    y = np.array([norm_b[pw] for pw in common])

    r, p = stats.spearmanr(x, y)

    # Bootstrap CI
    rng = np.random.default_rng(RANDOM_SEED)
    boot_rs = []
    for _ in range(2000):
        idx = rng.choice(len(x), size=len(x), replace=True)
        if len(np.unique(idx)) < 3:
            continue
        br, _ = stats.spearmanr(x[idx], y[idx])
        if np.isfinite(br):
            boot_rs.append(br)
    ci_low = float(np.percentile(boot_rs, 2.5)) if boot_rs else float("nan")
    ci_high = float(np.percentile(boot_rs, 97.5)) if boot_rs else float("nan")

    # Permutation p-value (two-sided)
    perm_count = 0
    n_perm = 5000
    for _ in range(n_perm):
        perm_r, _ = stats.spearmanr(x, rng.permutation(y))
        if abs(perm_r) >= abs(r):
            perm_count += 1
    p_perm = (perm_count + 1) / (n_perm + 1)

    return {
        "species_a": name_a,
        "species_b": name_b,
        "n_pathways_a": len(nes_a),
        "n_pathways_b": len(nes_b),
        "n_common_pathways": len(common),
        "spearman_r": float(r),
        "spearman_p": float(p),
        "bootstrap_ci": [ci_low, ci_high],
        "perm_p": float(p_perm),
        "common_pathways": common,
    }


# ── Main ──────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="E4/E5 multi-species NES")
    parser.add_argument("--skip-fgsea", action="store_true",
                        help="Skip fGSEA, load cached NES CSVs")
    args = parser.parse_args()

    all_nes = {}  # {key: {pathway: NES}}

    if not args.skip_fgsea:
        # ── Mouse (6 tissues, fresh Enrichr KEGG) ─────────────────
        print("\n=== Mouse (fresh Enrichr KEGG_2019_Mouse) ===")
        for tissue in MOUSE_TISSUES:
            nes = run_mouse_kegg(tissue)
            if nes:
                all_nes[f"mouse_{tissue}"] = nes
                # Save CSV
                df = pd.DataFrame([{"pathway": k, "NES": v} for k, v in nes.items()])
                df.to_csv(EVAL_DIR / f"E4_mouse_{tissue}_kegg_nes.csv", index=False)

        # ── Drosophila ────────────────────────────────────────────
        print("\n=== Drosophila (GLDS-207) ===")
        nes = run_drosophila_kegg()
        if nes:
            all_nes["drosophila_GLDS207"] = nes

        # ── Arabidopsis ───────────────────────────────────────────
        for glds_id, label in [("37", "arabidopsis_GLDS37"),
                                ("120", "arabidopsis_GLDS120")]:
            print(f"\n=== Arabidopsis (GLDS-{glds_id}) ===")
            nes = run_arabidopsis_pathway(glds_id, label)
            if nes:
                all_nes[label] = nes

    else:
        # Load cached NES from CSVs
        print("\n=== Loading cached NES ===")
        for csv_path in sorted(EVAL_DIR.glob("E4_*_nes.csv")):
            df = pd.read_csv(csv_path)
            if "NES" in df.columns and "pathway" in df.columns:
                label = csv_path.stem.replace("E4_", "").replace("_kegg_nes", "").replace("_pathway_nes", "")
                all_nes[label] = dict(zip(df["pathway"], df["NES"].astype(float)))
                print(f"  {label}: {len(all_nes[label])} pathways")

    # ── Pairwise Concordance (E4) ─────────────────────────────────
    print("\n=== E4: Pairwise NES Concordance ===")
    keys = sorted(all_nes.keys())
    print(f"  Available: {keys}")

    concordance = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            print(f"\n  {a} vs {b}")
            result = compute_concordance(all_nes[a], all_nes[b], a, b)
            concordance.append(result)
            if "spearman_r" in result:
                print(f"    r={result['spearman_r']:.3f}, "
                      f"n={result['n_common_pathways']}, "
                      f"p={result['perm_p']:.4f}")
            else:
                print(f"    {result.get('status', 'no data')}")

    e4_output = {
        "task": "E4",
        "description": "Multi-species KEGG NES concordance",
        "species_data": {k: {"n_pathways": len(v)} for k, v in all_nes.items()},
        "pairwise_concordance": concordance,
    }
    save_json(e4_output, EVAL_DIR / "E4_multispecies_nes.json")

    # ── Phylogenetic Summary (E5) ─────────────────────────────────
    print("\n=== E5: Phylogenetic Distance Summary ===")
    phylo_distances = {
        ("mouse", "drosophila"): 800,
        ("mouse", "arabidopsis"): 1500,
        ("drosophila", "arabidopsis"): 1500,
    }

    e5_results = []
    for result in concordance:
        sp_a = result["species_a"].split("_")[0]
        sp_b = result["species_b"].split("_")[0]
        dist = phylo_distances.get((sp_a, sp_b)) or phylo_distances.get((sp_b, sp_a))

        entry = {
            "pair": f"{result['species_a']} vs {result['species_b']}",
            "species_pair": f"{sp_a}-{sp_b}",
            "phylo_distance_mya": dist,
            "n_common_pathways": result["n_common_pathways"],
        }
        if "spearman_r" in result:
            entry["spearman_r"] = result["spearman_r"]
            entry["perm_p"] = result["perm_p"]
        else:
            entry["status"] = result.get("status", "no_data")

        e5_results.append(entry)

    save_json({"task": "E5", "description": "Phylogenetic distance vs concordance",
               "results": e5_results}, EVAL_DIR / "E5_phylogenetic.json")

    print("\n=== E4/E5 Complete ===")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
rrrm1_f2b_pseudobulk_fgsea.py  —  F2-B: Cell-type pseudo-bulk fGSEA

For each tissue × cell_type:
  1. Aggregate raw counts per animal → pseudo-bulk matrix [8 × n_genes]
  2. PyDESeq2: FLT(4) vs GC(4) Wald test
  3. gseapy.prerank on Wald statistic ranking → Hallmark NES per pathway
  4. Compare cell-type NES vs v1.0 bulk NES (same tissue, RR-8 mission)

Requires: pydeseq2, gseapy (pip install pydeseq2 gseapy)

Usage:
    python3 rrrm1_f2b_pseudobulk_fgsea.py --all
    python3 rrrm1_f2b_pseudobulk_fgsea.py --tissue muscle

Inputs:
  *_labeled.h5ad    (raw counts + condition + srx + animal_id)
  *_hardened.h5ad   (broad_celltype labels — matched by obs_names)

Outputs:
  v2/evaluation/F2B_pseudobulk_fgsea.json
  v2/figures/F2B_celltype_nes_heatmap.html
"""

import argparse
import json
import os
import re
import sys
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────
SCRATCH = Path(os.environ.get("SCRATCH_DIR", "/path/to/scratch")) / "rrrm1_scrna"
HARDENED_DIR = SCRATCH / "downstream_initial" / "hardening" / "objects"
BASE_DIR = Path(os.environ.get("HOME")) / "rrrm1_scrna"

EVAL_DIR = BASE_DIR / "evaluation"
FIG_DIR = BASE_DIR / "figures"
PROCESSED_DIR = BASE_DIR / "processed" / "F2B"

BULK_FGSEA_DIR = Path(os.environ.get("HOME")) / "rrrm1_scrna" / "v1_bulk_fgsea"
GMT_DIR = Path(os.environ.get("HOME")) / "rrrm1_scrna" / "gene_sets"

TISSUE_OSD = {
    "blood":  "OSD-918",
    "eye":    "OSD-920",
    "muscle": "OSD-924",
    "skin":   "OSD-934",
}

# Minimum cells per condition to include a cell type
MIN_CELLS_PER_CONDITION = 20
# Minimum pseudo-bulk samples per group for DESeq2
MIN_SAMPLES_PER_GROUP = 3


# ── Gene set utilities ─────────────────────────────────────────────────────

def download_mouse_hallmark_gmt(out_dir: Path) -> Path:
    """Download Hallmark mouse gene set if not cached."""
    gmt_path = out_dir / "mh.all.v2024.1.Mm.symbols.gmt"
    if gmt_path.exists():
        return gmt_path
    out_dir.mkdir(parents=True, exist_ok=True)
    url = ("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
           "2024.1.Mm/mh.all.v2024.1.Mm.symbols.gmt")
    print(f"  Downloading mouse Hallmark GMT from MSigDB...")
    urllib.request.urlretrieve(url, gmt_path)
    print(f"  Saved: {gmt_path}")
    return gmt_path


def normalize_pathway_name(name: str) -> str:
    """Standardize pathway names for cross-study matching."""
    name = name.upper().strip().replace(" ", "_").strip("_")
    if not name.startswith("HALLMARK_"):
        name = "HALLMARK_" + name
    return name


def load_hallmark_gmt(gmt_path: Path) -> dict:
    """Parse GMT → {pathway_name: set(genes)}."""
    gene_sets = {}
    with open(gmt_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = normalize_pathway_name(parts[0])
            genes = set(parts[2:])
            gene_sets[name] = genes
    return gene_sets


# ── Data loading ───────────────────────────────────────────────────────────

def load_labeled(tissue: str) -> sc.AnnData:
    osd = TISSUE_OSD[tissue]
    path = SCRATCH / osd / f"{osd}_{tissue}_labeled.h5ad"
    if not path.exists():
        raise FileNotFoundError(f"Labeled h5ad not found: {path}")
    print(f"  Loading labeled: {path}")
    return sc.read_h5ad(path)


def load_hardened(tissue: str) -> sc.AnnData:
    path = HARDENED_DIR / f"RRRM1_{tissue}_hardened.h5ad"
    if not path.exists():
        raise FileNotFoundError(f"Hardened h5ad not found: {path}")
    print(f"  Loading hardened: {path}")
    return sc.read_h5ad(path)


def transfer_celltype_labels(labeled: sc.AnnData, hardened: sc.AnnData) -> sc.AnnData:
    """
    Transfer broad_celltype from hardened h5ad to labeled h5ad by obs_name matching.
    Barcodes in labeled are prefixed as {SRX}_{barcode}, while hardened uses the
    same prefix scheme from the original merged pipeline.
    Returns labeled adata with 'broad_celltype' column added.
    """
    # Build lookup: barcode → broad_celltype from hardened
    ct_series = hardened.obs["broad_celltype"]

    # Try direct match first
    common = labeled.obs_names.intersection(ct_series.index)
    if len(common) > 100:
        labeled.obs["broad_celltype"] = (
            ct_series.reindex(labeled.obs_names).astype(object).fillna("unknown")
        )
        n_assigned = (labeled.obs["broad_celltype"] != "unknown").sum()
        print(f"  Transferred broad_celltype to {n_assigned}/{len(labeled)} cells")
        return labeled

    # If no direct match, try barcode-only matching (strip SRX prefix)
    # labeled obs_names: SRX28491856_ATCGATCG-1
    # hardened obs_names might be from merged pipeline (different prefix)
    # Fall back: use per-SRX barcode without prefix for matching
    hardened_bc_to_ct = {}
    for obs_name, ct in ct_series.items():
        # Extract bare barcode: everything after last underscore (or last -)
        bc = obs_name.split("_")[-1] if "_" in obs_name else obs_name
        hardened_bc_to_ct[bc] = ct

    def get_ct(obs_name):
        bc = obs_name.split("_")[-1] if "_" in obs_name else obs_name
        return hardened_bc_to_ct.get(bc, "unknown")

    labeled.obs["broad_celltype"] = [get_ct(n) for n in labeled.obs_names]
    n_assigned = (labeled.obs["broad_celltype"] != "unknown").sum()
    print(f"  Transferred broad_celltype (barcode-only match): {n_assigned}/{len(labeled)} cells")
    return labeled


# ── Pseudo-bulk aggregation ────────────────────────────────────────────────

def aggregate_pseudobulk(adata: sc.AnnData) -> pd.DataFrame:
    """
    Sum raw counts per animal. Returns DataFrame [animals × genes].
    Uses .raw if available (raw counts before normalization), else .X.
    """
    import scipy.sparse as sp

    if adata.raw is not None:
        X = adata.raw.X
        gene_names = adata.raw.var_names
    else:
        X = adata.X
        gene_names = adata.var_names

    if sp.issparse(X):
        X = X.toarray()

    animals = adata.obs["animal_id"].values
    unique_animals = adata.obs["animal_id"].unique()

    rows = {}
    for animal in unique_animals:
        mask = animals == animal
        rows[animal] = X[mask, :].sum(axis=0)

    bulk = pd.DataFrame(rows, index=gene_names).T  # [animals × genes]
    return bulk


# ── DESeq2-style ranking via PyDESeq2 ─────────────────────────────────────

def run_pydeseq2(bulk_df: pd.DataFrame, condition_series: pd.Series) -> pd.Series:
    """
    Run PyDESeq2 on pseudo-bulk matrix [animals × genes].
    Returns Series: gene → Wald stat (for preranked fGSEA).
    Falls back to t-test if pydeseq2 not available.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        print("  pydeseq2 not found — falling back to Welch t-test ranking")
        return _ttest_ranking(bulk_df, condition_series)

    # Build sample metadata
    metadata = pd.DataFrame({"condition": condition_series})
    metadata = metadata.loc[bulk_df.index]

    # Count matrix: integer (DESeq2 requires counts)
    counts = bulk_df.astype(int)

    # Filter: remove all-zero genes
    nonzero = (counts > 0).any(axis=0)
    counts = counts.loc[:, nonzero]

    try:
        dds = DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design_factors="condition",
            ref_level=["condition", "GC"],
            quiet=True,
        )
        dds.deseq2()
        ds = DeseqStats(dds, contrast=["condition", "FLT", "GC"], quiet=True)
        ds.summary()
        res = ds.results_df.dropna(subset=["stat"])
        ranking = res["stat"].astype(float)
        ranking = ranking.replace([np.inf, -np.inf], np.nan).dropna()
        print(f"  PyDESeq2: {len(ranking)} genes ranked")
        return ranking
    except Exception as e:
        print(f"  PyDESeq2 failed ({e}) — falling back to t-test")
        return _ttest_ranking(bulk_df, condition_series)


def _ttest_ranking(bulk_df: pd.DataFrame, condition_series: pd.Series) -> pd.Series:
    """Fallback: Welch t-test, rank by -log10(p) × sign(mean_FLT - mean_GC)."""
    flt_animals = condition_series[condition_series == "FLT"].index
    gc_animals = condition_series[condition_series == "GC"].index

    flt = np.log1p(bulk_df.loc[bulk_df.index.isin(flt_animals)].values.astype(float))
    gc = np.log1p(bulk_df.loc[bulk_df.index.isin(gc_animals)].values.astype(float))

    genes = bulk_df.columns
    stats_arr = []
    for i in range(len(genes)):
        a, b = flt[:, i], gc[:, i]
        if np.std(a) == 0 and np.std(b) == 0:
            stats_arr.append(0.0)
            continue
        t, p = stats.ttest_ind(a, b, equal_var=False)
        sign = np.sign(np.mean(a) - np.mean(b))
        stats_arr.append(sign * (-np.log10(max(p, 1e-300))))

    ranking = pd.Series(stats_arr, index=genes)
    print(f"  t-test ranking: {len(ranking)} genes")
    return ranking


# ── fGSEA via gseapy ──────────────────────────────────────────────────────

def run_preranked_fgsea(ranking: pd.Series, gene_sets: dict,
                        min_size=10, max_size=500) -> pd.DataFrame:
    """
    Run gseapy prerank on gene ranking → NES per pathway.
    Returns DataFrame with columns: pathway, NES, pval, fdr.
    """
    try:
        import gseapy as gp
    except ImportError:
        print("  gseapy not found — cannot run fGSEA")
        return pd.DataFrame(columns=["pathway", "NES", "pval", "fdr"])

    rnk = ranking.sort_values(ascending=False).dropna()
    rnk = rnk[~rnk.index.duplicated(keep="first")]

    try:
        res = gp.prerank(
            rnk=rnk,
            gene_sets=gene_sets,
            min_size=min_size,
            max_size=max_size,
            permutation_num=1000,
            seed=42,
            verbose=False,
        )
        df = res.res2d.reset_index()
        # Column names vary by gseapy version
        nes_col = "NES" if "NES" in df.columns else "nes"
        p_col = "NOM p-val" if "NOM p-val" in df.columns else "pval"
        fdr_col = "FDR q-val" if "FDR q-val" in df.columns else "fdr"
        term_col = "Term" if "Term" in df.columns else "term"

        out = pd.DataFrame({
            "pathway": df[term_col].apply(normalize_pathway_name),
            "NES": pd.to_numeric(df[nes_col], errors="coerce"),
            "pval": pd.to_numeric(df[p_col], errors="coerce"),
            "fdr": pd.to_numeric(df[fdr_col], errors="coerce"),
        }).dropna(subset=["NES"])
        print(f"  fGSEA: {len(out)} pathways")
        return out
    except Exception as e:
        print(f"  gseapy failed: {e}")
        return pd.DataFrame(columns=["pathway", "NES", "pval", "fdr"])


# ── Bulk NES comparison ────────────────────────────────────────────────────

def load_bulk_nes(tissue: str, mission: str = "RR-8") -> pd.Series:
    """Load v1.0 bulk fGSEA NES for given tissue and mission."""
    tissue_map = {
        "blood": None,       # no v1.0 bulk blood
        "eye": "eye",
        "muscle": "gastrocnemius",
        "skin": "skin",
    }
    bulk_tissue = tissue_map.get(tissue)
    if bulk_tissue is None:
        return pd.Series(dtype=float)

    path = BULK_FGSEA_DIR / bulk_tissue / f"{mission}_fgsea_hallmark.csv"
    if not path.exists():
        # Try alternative missions
        for alt in ["RR-8", "RR-1", "RR-3", "RR-6", "RR-9"]:
            alt_path = BULK_FGSEA_DIR / bulk_tissue / f"{alt}_fgsea_hallmark.csv"
            if alt_path.exists():
                path = alt_path
                print(f"  Using bulk NES from {alt} (mission {mission} not found)")
                break
        else:
            print(f"  No bulk NES found for {tissue}")
            return pd.Series(dtype=float)

    df = pd.read_csv(path)
    df["pathway"] = df["pathway"].apply(normalize_pathway_name)
    return df.set_index("pathway")["NES"].astype(float)


def spearman_with_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 1000, seed: int = 42):
    if len(x) < 5:
        return {"r": float("nan"), "ci_low": float("nan"), "ci_high": float("nan"),
                "p": float("nan"), "n": len(x)}
    r, _ = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    n = len(x)
    boot_r = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            r_b, _ = stats.spearmanr(x[idx], y[idx])
            boot_r.append(r_b)
        except Exception:
            pass
    # Permutation p
    perm_r = []
    for _ in range(1000):
        y_perm = rng.permutation(y)
        r_p, _ = stats.spearmanr(x, y_perm)
        perm_r.append(r_p)
    p_perm = float(np.mean(np.abs(perm_r) >= np.abs(r)))
    p_perm = max(p_perm, 1 / 1000)
    return {
        "r": float(r),
        "ci_low": float(np.percentile(boot_r, 2.5)) if boot_r else float("nan"),
        "ci_high": float(np.percentile(boot_r, 97.5)) if boot_r else float("nan"),
        "p": p_perm,
        "n": int(len(x)),
    }


# ── HTML figure ────────────────────────────────────────────────────────────

def make_heatmap_html(all_nes: dict) -> str:
    """NES heatmap: x=cell_type, y=pathway, one panel per tissue."""
    data_json = json.dumps(all_nes)
    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>F2-B RRRM-1 Cell-Type Pseudo-bulk fGSEA NES</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; font-size: 11px; margin: 20px; }}
  h2 {{ font-size: 13px; font-weight: bold; }}
  .panel {{ display: inline-block; vertical-align: top; margin: 20px; }}
  .panel-title {{ font-weight: bold; font-size: 12px; margin-bottom: 5px; }}
</style>
</head>
<body>
<h2>F2-B: RRRM-1 scRNA-seq Pseudo-bulk Hallmark fGSEA NES (FLT vs GC)</h2>
<div id="container"></div>
<script>
const allNES = {data_json};

const colorScale = d3.scaleDiverging()
  .domain([-2, 0, 2])
  .interpolator(d3.interpolateRdBu).clamp(true);

const container = d3.select("#container");

for (const [tissue, ctData] of Object.entries(allNES)) {{
  const cellTypes = Object.keys(ctData);
  if (cellTypes.length === 0) continue;

  // Get all pathways
  const allPathways = new Set();
  cellTypes.forEach(ct => Object.keys(ctData[ct]).forEach(p => allPathways.add(p)));
  const pathways = Array.from(allPathways).sort();

  const cellW = 50, rowH = 12;
  const marginL = 200, marginT = 80;
  const W = marginL + cellTypes.length * cellW + 40;
  const H = marginT + pathways.length * rowH + 20;

  const div = container.append("div").attr("class","panel");
  div.append("div").attr("class","panel-title").text(tissue);
  const svg = div.append("svg").attr("width", W).attr("height", H);

  // Heatmap cells
  cellTypes.forEach((ct, ci) => {{
    pathways.forEach((pw, pi) => {{
      const nes = ctData[ct][pw] || 0;
      svg.append("rect")
        .attr("x", marginL + ci * cellW)
        .attr("y", marginT + pi * rowH)
        .attr("width", cellW - 1)
        .attr("height", rowH - 1)
        .attr("fill", colorScale(-nes));  // red=positive NES (FLT > GC)
    }});
  }});

  // Column labels (cell types)
  cellTypes.forEach((ct, ci) => {{
    svg.append("text")
      .attr("x", marginL + ci * cellW + cellW/2)
      .attr("y", marginT - 4)
      .attr("text-anchor","end")
      .attr("transform", `rotate(-45, ${{marginL + ci * cellW + cellW/2}}, ${{marginT - 4}})`)
      .style("font-size","9px")
      .text(ct.replace(/_/g," "));
  }});

  // Row labels (pathways)
  pathways.forEach((pw, pi) => {{
    svg.append("text")
      .attr("x", marginL - 4)
      .attr("y", marginT + pi * rowH + rowH * 0.75)
      .attr("text-anchor","end")
      .style("font-size","8px")
      .text(pw.replace("HALLMARK_","").replace(/_/g," ").toLowerCase());
  }});
}}
</script>
</body>
</html>
"""


# ── Main ───────────────────────────────────────────────────────────────────

def run_tissue(tissue: str, gmt_path: Path, gene_sets: dict) -> dict:
    print(f"\n[F2-B] {tissue}")

    # Load data
    labeled = load_labeled(tissue)
    hardened = load_hardened(tissue)

    # Transfer cell type labels
    labeled = transfer_celltype_labels(labeled, hardened)

    # Condition series
    condition_series = labeled.obs.set_index("animal_id")["condition"].drop_duplicates()

    tissue_results = {}
    all_nes_tissue = {}  # cell_type → {pathway: NES}

    cell_types = [ct for ct in labeled.obs["broad_celltype"].unique() if ct != "unknown"]

    for ct in sorted(cell_types):
        sub = labeled[labeled.obs["broad_celltype"] == ct]
        n_flt = (sub.obs["condition"] == "FLT").sum()
        n_gc = (sub.obs["condition"] == "GC").sum()

        if n_flt < MIN_CELLS_PER_CONDITION or n_gc < MIN_CELLS_PER_CONDITION:
            print(f"  SKIP {ct}: n_FLT={n_flt} n_GC={n_gc} (< {MIN_CELLS_PER_CONDITION})")
            continue

        # Check animal coverage per group
        n_animals_flt = sub.obs[sub.obs["condition"] == "FLT"]["animal_id"].nunique()
        n_animals_gc = sub.obs[sub.obs["condition"] == "GC"]["animal_id"].nunique()
        if n_animals_flt < MIN_SAMPLES_PER_GROUP or n_animals_gc < MIN_SAMPLES_PER_GROUP:
            print(f"  SKIP {ct}: n_animals FLT={n_animals_flt} GC={n_animals_gc}")
            continue

        print(f"  Processing {ct} (n={len(sub)}, {n_animals_flt}+{n_animals_gc} animals)")

        # Pseudo-bulk aggregation
        bulk_df = aggregate_pseudobulk(sub)

        # Map animal_id → condition
        animal_cond = sub.obs[["animal_id", "condition"]].drop_duplicates().set_index("animal_id")["condition"]
        animal_cond = animal_cond.loc[bulk_df.index]

        # DESeq2 ranking
        ranking = run_pydeseq2(bulk_df, animal_cond)

        # Convert ENSMUSG → gene symbols for gseapy/GMT compatibility
        source = sub.raw.to_adata() if sub.raw is not None else sub
        if (ranking.index[0].startswith("ENSMUSG")
                and "gene_symbols" in source.var.columns):
            sym_map = dict(zip(source.var_names.astype(str),
                               source.var["gene_symbols"].astype(str)))
            ranking.index = ranking.index.map(lambda x: sym_map.get(x, x))
            ranking = ranking[~ranking.index.duplicated(keep="first")]

        # fGSEA
        fgsea_df = run_preranked_fgsea(ranking, gene_sets)

        if fgsea_df.empty:
            continue

        # Save per-cell-type fGSEA CSV
        out_csv_dir = PROCESSED_DIR / tissue
        out_csv_dir.mkdir(parents=True, exist_ok=True)
        ct_safe = ct.replace("/", "_").replace(" ", "_")
        fgsea_df.to_csv(out_csv_dir / f"{ct_safe}_fgsea_hallmark.csv", index=False)

        # Build NES dict for heatmap
        nes_dict = fgsea_df.set_index("pathway")["NES"].to_dict()
        all_nes_tissue[ct] = nes_dict

        # Compare to v1.0 bulk NES
        bulk_nes = load_bulk_nes(tissue)
        if not bulk_nes.empty:
            common_pw = fgsea_df["pathway"].values
            common_pw = [p for p in common_pw if p in bulk_nes.index]
            if len(common_pw) >= 5:
                sc_nes = fgsea_df.set_index("pathway")["NES"]
                x = sc_nes.reindex(common_pw).values
                y = bulk_nes.reindex(common_pw).values
                valid = ~(np.isnan(x) | np.isnan(y))
                if valid.sum() >= 5:
                    corr = spearman_with_ci(x[valid], y[valid])
                    print(f"    vs v1.0 bulk: r={corr['r']:.3f} p={corr['p']:.3f} n={corr['n']}")
                else:
                    corr = {"r": float("nan"), "n": 0}
            else:
                corr = {"r": float("nan"), "n": 0}
        else:
            corr = {"r": float("nan"), "n": 0}

        tissue_results[ct] = {
            "n_cells": int(len(sub)),
            "n_animals_flt": int(n_animals_flt),
            "n_animals_gc": int(n_animals_gc),
            "n_pathways": int(len(fgsea_df)),
            "vs_bulk_r": corr.get("r"),
            "vs_bulk_p": corr.get("p"),
            "vs_bulk_n_pathways": corr.get("n"),
            "top_enriched": fgsea_df.nlargest(3, "NES")["pathway"].tolist(),
            "top_depleted": fgsea_df.nsmallest(3, "NES")["pathway"].tolist(),
        }

    return tissue_results, all_nes_tissue


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tissue", choices=list(TISSUE_OSD.keys()))
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if not args.tissue and not args.all:
        parser.error("Specify --tissue or --all")

    for d in [EVAL_DIR, FIG_DIR, PROCESSED_DIR]:
        d.mkdir(parents=True, exist_ok=True)

    # Download / load gene sets
    gmt_path = download_mouse_hallmark_gmt(GMT_DIR)
    gene_sets = load_hallmark_gmt(gmt_path)
    print(f"Loaded {len(gene_sets)} Hallmark pathways from {gmt_path.name}")

    tissues = list(TISSUE_OSD.keys()) if args.all else [args.tissue]

    all_results = {}
    all_nes = {}
    for tissue in tissues:
        try:
            res, nes = run_tissue(tissue, gmt_path, gene_sets)
            all_results[tissue] = res
            all_nes[tissue] = nes
        except FileNotFoundError as e:
            print(f"  SKIP {tissue}: {e}")

    if not all_results:
        print("No tissues processed.")
        sys.exit(1)

    # Save JSON
    out_json = EVAL_DIR / "F2B_pseudobulk_fgsea.json"
    with open(out_json, "w") as f:
        json.dump({"task": "F2-B", "results": all_results}, f, indent=2)
    print(f"\nSaved: {out_json}")

    # Save HTML heatmap
    out_html = FIG_DIR / "F2B_celltype_nes_heatmap.html"
    out_html.write_text(make_heatmap_html(all_nes))
    print(f"Saved: {out_html}")


if __name__ == "__main__":
    main()

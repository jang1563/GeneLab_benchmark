#!/usr/bin/env python3
"""
rrrm1_f2a_composition.py  —  F2-A: Cell-type composition shift analysis

For each tissue, compute per-animal cell-type proportions and test
FLT vs GC with Wilcoxon rank-sum (n=4 vs n=4, two-sided).
BH FDR correction across all tissue × cell_type pairs.

Usage:
    python3 rrrm1_f2a_composition.py --all
    python3 rrrm1_f2a_composition.py --tissue blood

Inputs:  *_hardened.h5ad  (obs: condition, animal_id, broad_celltype)
Outputs: v2/evaluation/F2A_composition.json
         v2/figures/F2A_composition.html
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from statsmodels.stats.multitest import multipletests

# ── Paths ──────────────────────────────────────────────────────────────────
SCRATCH = Path("/athena/masonlab/scratch/users/jak4013/rrrm1_scrna")
HARDENED_DIR = SCRATCH / "downstream_initial" / "hardening" / "objects"
BASE_DIR = Path("/home/fs01/jak4013/rrrm1_scrna")

EVAL_DIR = BASE_DIR / "evaluation"
FIG_DIR = BASE_DIR / "figures"

TISSUE_OSD = {
    "blood":  "OSD-918",
    "eye":    "OSD-920",
    "muscle": "OSD-924",
    "skin":   "OSD-934",
}

# Okabe-Ito palette (cell types)
OKABE_ITO = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999",
    "#000000", "#88CCEE",
]


def load_hardened(tissue: str) -> sc.AnnData:
    """Load hardened h5ad; try two possible paths."""
    osd = TISSUE_OSD[tissue]
    # Prefer per-SRX labeled + re-processed path
    candidates = [
        HARDENED_DIR / f"RRRM1_{tissue}_hardened.h5ad",
        SCRATCH / osd / f"{osd}_{tissue}_labeled.h5ad",  # fallback: labeled only
    ]
    for p in candidates:
        if p.exists():
            print(f"  Loading {p}")
            return sc.read_h5ad(p)
    raise FileNotFoundError(
        f"No hardened h5ad found for {tissue}. Tried:\n" +
        "\n".join(f"  {c}" for c in candidates)
    )


def compute_proportions(adata: sc.AnnData, tissue: str) -> pd.DataFrame:
    """Return DataFrame: rows=animals, cols=cell_types, values=proportions."""
    required = ["condition", "animal_id", "broad_celltype"]
    missing = [c for c in required if c not in adata.obs.columns]
    if missing:
        raise ValueError(f"[{tissue}] Missing obs columns: {missing}")

    records = []
    for animal in adata.obs["animal_id"].unique():
        sub = adata[adata.obs["animal_id"] == animal]
        condition = sub.obs["condition"].iloc[0]
        total = len(sub)
        counts = sub.obs["broad_celltype"].value_counts()
        row = {"animal_id": animal, "condition": condition, "n_cells": total}
        for ct, n in counts.items():
            row[ct] = n / total
        records.append(row)

    df = pd.DataFrame(records).fillna(0.0)
    return df


def wilcoxon_flt_vs_gc(prop_df: pd.DataFrame, cell_types: list) -> dict:
    """Wilcoxon rank-sum per cell type (FLT vs GC)."""
    flt = prop_df[prop_df["condition"] == "FLT"]
    gc = prop_df[prop_df["condition"] == "GC"]
    results = {}
    for ct in cell_types:
        a = flt[ct].values if ct in flt.columns else np.zeros(len(flt))
        b = gc[ct].values if ct in gc.columns else np.zeros(len(gc))
        stat, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        results[ct] = {
            "flt_mean": float(np.mean(a)),
            "flt_std": float(np.std(a)),
            "gc_mean": float(np.mean(b)),
            "gc_std": float(np.std(b)),
            "delta": float(np.mean(a) - np.mean(b)),
            "p_raw": float(p),
            "n_flt": int(len(a)),
            "n_gc": int(len(b)),
        }
    return results


def bh_correct(all_results: dict) -> dict:
    """Apply BH FDR across all tissue × cell_type pairs."""
    keys = []
    pvals = []
    for tissue, ct_dict in all_results.items():
        for ct, vals in ct_dict.items():
            keys.append((tissue, ct))
            pvals.append(vals["p_raw"])

    _, padj, _, _ = multipletests(pvals, method="fdr_bh")
    for (tissue, ct), q in zip(keys, padj):
        all_results[tissue][ct]["padj"] = float(q)
    return all_results


def make_composition_html(all_results: dict, prop_tables: dict) -> str:
    """Generate D3.js v7 stacked bar chart — one panel per tissue."""
    tissues = list(all_results.keys())

    # Build data for each tissue
    tissue_data = {}
    for tissue, prop_df in prop_tables.items():
        ct_cols = [c for c in prop_df.columns if c not in ("animal_id", "condition", "n_cells")]
        ct_cols = sorted(ct_cols, key=lambda c: -prop_df[ct_cols].mean()[c])

        rows = []
        for _, row in prop_df.iterrows():
            rows.append({
                "animal_id": row["animal_id"],
                "condition": row["condition"],
                "props": {ct: float(row[ct]) for ct in ct_cols},
            })
        tissue_data[tissue] = {"rows": rows, "cell_types": ct_cols}

    data_json = json.dumps(tissue_data)
    results_json = json.dumps(all_results)

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>F2-A RRRM-1 Cell-Type Composition</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; font-size: 12px; margin: 20px; background: #fff; }}
  h2 {{ font-size: 14px; font-weight: bold; }}
  .panel {{ display: inline-block; vertical-align: top; margin: 10px; }}
  .bar {{ stroke: none; }}
  .axis text {{ font-size: 9px; }}
  .legend-item {{ font-size: 9px; }}
  .sig {{ font-weight: bold; color: #D55E00; }}
  .stats-table {{ font-size: 9px; border-collapse: collapse; margin-top: 8px; }}
  .stats-table td, .stats-table th {{ border: 1px solid #ddd; padding: 3px 6px; }}
  .stats-table th {{ background: #f5f5f5; }}
</style>
</head>
<body>
<h2>F2-A: RRRM-1 scRNA-seq Cell-Type Composition (FLT vs GC)</h2>
<p style="font-size:10px;color:#666">Wilcoxon rank-sum, n=4 FLT vs n=4 GC per tissue. BH FDR correction across all comparisons.</p>
<div id="charts"></div>
<script>
const tissueData = {data_json};
const results = {results_json};

const palette = {json.dumps(OKABE_ITO)};
const W = 280, H = 300, margin = {{top: 30, right: 10, bottom: 80, left: 50}};
const innerW = W - margin.left - margin.right;
const innerH = H - margin.top - margin.bottom;

const container = d3.select("#charts");

for (const [tissue, data] of Object.entries(tissueData)) {{
  const cellTypes = data.cell_types;
  const colorScale = d3.scaleOrdinal().domain(cellTypes).range(palette);

  const panel = container.append("div").attr("class", "panel");
  panel.append("div").style("font-weight","bold").style("font-size","11px")
    .text(tissue.charAt(0).toUpperCase() + tissue.slice(1));

  const svg = panel.append("svg").attr("width", W).attr("height", H);
  const g = svg.append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const animals = data.rows.map(r => r.animal_id);
  const x = d3.scaleBand().domain(animals).range([0, innerW]).padding(0.2);
  const y = d3.scaleLinear().domain([0, 1]).range([innerH, 0]);

  // Stacked bars
  const stackData = d3.stack()
    .keys(cellTypes)
    .value((d, key) => d.props[key] || 0)(data.rows);

  g.selectAll(".layer")
    .data(stackData)
    .join("g")
    .attr("fill", d => colorScale(d.key))
    .selectAll("rect")
    .data(d => d)
    .join("rect")
    .attr("x", (d,i) => x(data.rows[i].animal_id))
    .attr("y", d => y(d[1]))
    .attr("height", d => y(d[0]) - y(d[1]))
    .attr("width", x.bandwidth());

  // Axes
  g.append("g").attr("transform", `translate(0,${{innerH}})`).call(
    d3.axisBottom(x).tickSize(2)
  ).selectAll("text").attr("transform","rotate(-40)").style("text-anchor","end").style("font-size","8px");

  g.append("g").call(d3.axisLeft(y).ticks(5).tickFormat(d3.format(".0%")).tickSize(2));

  // Y label
  svg.append("text").attr("transform","rotate(-90)")
    .attr("x", -(H/2)).attr("y", 12)
    .attr("text-anchor","middle").style("font-size","9px")
    .text("Proportion");

  // Color legend (right side)
  const legendX = innerW + margin.left + 2;
  cellTypes.forEach((ct, i) => {{
    svg.append("rect").attr("x", legendX).attr("y", margin.top + i*14)
      .attr("width",10).attr("height",10).attr("fill",colorScale(ct));
    svg.append("text").attr("x", legendX+13).attr("y", margin.top + i*14+9)
      .style("font-size","8px").text(ct.replace(/_/g," "));
  }});

  // Stats table
  const tissueRes = results[tissue];
  let tableHtml = '<table class="stats-table"><tr><th>Cell type</th><th>FLT mean</th><th>GC mean</th><th>Δ</th><th>p</th><th>padj</th></tr>';
  for (const [ct, v] of Object.entries(tissueRes)) {{
    const sig = v.padj < 0.05 ? ' class="sig"' : '';
    tableHtml += `<tr${{sig}}><td>${{ct.replace(/_/g," ")}}</td>`
      + `<td>${{(v.flt_mean*100).toFixed(1)}}%</td>`
      + `<td>${{(v.gc_mean*100).toFixed(1)}}%</td>`
      + `<td>${{v.delta >= 0 ? '+' : ''}}${{(v.delta*100).toFixed(1)}}%</td>`
      + `<td>${{v.p_raw < 0.001 ? '<0.001' : v.p_raw.toFixed(3)}}</td>`
      + `<td>${{v.padj < 0.001 ? '<0.001' : v.padj.toFixed(3)}}</td></tr>`;
  }}
  tableHtml += '</table>';
  panel.append("div").html(tableHtml);
}}
</script>
</body>
</html>
"""
    return html


def run_tissue(tissue: str) -> tuple:
    print(f"\n[F2-A] {tissue}")
    adata = load_hardened(tissue)

    prop_df = compute_proportions(adata, tissue)
    cell_types = [c for c in prop_df.columns if c not in ("animal_id", "condition", "n_cells")]
    print(f"  {len(adata)} cells, {len(cell_types)} cell types, "
          f"{prop_df['condition'].value_counts().to_dict()}")

    tissue_results = wilcoxon_flt_vs_gc(prop_df, cell_types)
    return tissue_results, prop_df


def main():
    parser = argparse.ArgumentParser(description="F2-A cell type composition analysis")
    parser.add_argument("--tissue", choices=list(TISSUE_OSD.keys()))
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()

    if not args.tissue and not args.all:
        parser.error("Specify --tissue or --all")

    EVAL_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    tissues = list(TISSUE_OSD.keys()) if args.all else [args.tissue]

    all_results = {}
    prop_tables = {}
    for tissue in tissues:
        try:
            res, prop_df = run_tissue(tissue)
            all_results[tissue] = res
            prop_tables[tissue] = prop_df
        except FileNotFoundError as e:
            print(f"  SKIP: {e}")

    if not all_results:
        print("No tissues processed.")
        sys.exit(1)

    # BH FDR correction
    all_results = bh_correct(all_results)

    # Print significant findings
    print("\n=== Significant (padj < 0.05) ===")
    any_sig = False
    for tissue, ct_dict in all_results.items():
        for ct, v in ct_dict.items():
            if v["padj"] < 0.05:
                print(f"  {tissue} {ct}: FLT={v['flt_mean']:.3f} GC={v['gc_mean']:.3f} "
                      f"Δ={v['delta']:+.3f} padj={v['padj']:.3f}")
                any_sig = True
    if not any_sig:
        print("  None (n=4+4 is small; report raw p-values)")

    # Save JSON
    out_json = EVAL_DIR / "F2A_composition.json"
    with open(out_json, "w") as f:
        json.dump({"task": "F2-A", "results": all_results}, f, indent=2)
    print(f"\nSaved: {out_json}")

    # Save HTML figure
    out_html = FIG_DIR / "F2A_composition.html"
    html = make_composition_html(all_results, prop_tables)
    out_html.write_text(html)
    print(f"Saved: {out_html}")


if __name__ == "__main__":
    main()

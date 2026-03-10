#!/usr/bin/env python3
"""
F1: Cell-type-specific pathway heatmap for I4 PBMC snRNA-seq
Input:  v2/processed/F1_scrna/i4_snrnaseq_celltype_fgsea.csv
Output: v2/figures/F1_celltype_pathway_heatmap.html
"""

import json
import math
import numpy as np
import pandas as pd
from pathlib import Path

# ── Load data ──────────────────────────────────────────────────────────────
df = pd.read_csv("v2/processed/F1_scrna/i4_snrnaseq_celltype_fgsea.csv")

# Cell type display order (innate → adaptive → pseudobulk)
CT_ORDER = [
    "CD14+ Monocyte", "CD16+ Monocyte", "Dendritic Cell",
    "Natural Killer Cell",
    "B Cell", "CD4+ T Cell", "CD8+ T Cell", "Other T Cell",
    "Other", "PBMC Pseudobulk"
]

# Pathway short labels (strip HALLMARK_ prefix)
def short_name(p):
    return p.replace("HALLMARK_", "").replace("_", " ").title()

# Pivot NES and padj
nes_pivot  = df.pivot(index="pathway", columns="Cell_Type", values="NES")
padj_pivot = df.pivot(index="pathway", columns="Cell_Type", values="padj")

# Reorder columns
nes_pivot  = nes_pivot.reindex(columns=[c for c in CT_ORDER if c in nes_pivot.columns])
padj_pivot = padj_pivot.reindex(columns=nes_pivot.columns)

# Order rows: pathways with any sig result first, then by mean |NES|
has_sig   = (padj_pivot < 0.05).any(axis=1)
mean_anes = nes_pivot.abs().mean(axis=1).fillna(0)
sort_idx  = np.lexsort([-mean_anes.values, -has_sig.astype(int).values])
row_order = nes_pivot.index[sort_idx]
nes_pivot  = nes_pivot.loc[row_order]
padj_pivot = padj_pivot.loc[row_order]

pathways   = [short_name(p) for p in nes_pivot.index]
cell_types = list(nes_pivot.columns)

# Build cell data for D3
cells = []
for ri, pw in enumerate(nes_pivot.index):
    for ci, ct in enumerate(nes_pivot.columns):
        nes  = nes_pivot.loc[pw, ct]
        padj = padj_pivot.loc[pw, ct]
        cells.append({
            "row":  ri,
            "col":  ci,
            "nes":  None if (isinstance(nes, float) and math.isnan(nes)) else round(float(nes), 3),
            "sig":  bool(isinstance(padj, float) and not math.isnan(padj) and padj < 0.05),
            "padj": None if (isinstance(padj, float) and math.isnan(padj)) else round(float(padj), 4),
        })

# Summary stats
n_sig_total = (padj_pivot < 0.05).sum().sum()
sig_by_ct   = (padj_pivot < 0.05).sum().to_dict()

data = {
    "pathways":   pathways,
    "cell_types": cell_types,
    "cells":      cells,
    "n_sig":      int(n_sig_total),
    "sig_by_ct":  {k: int(v) for k, v in sig_by_ct.items()},
}

# ── HTML template ──────────────────────────────────────────────────────────
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>F1: I4 PBMC Cell-Type Pathway Heatmap</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; font-size: 12px; margin: 20px; background: #fff; }}
  h2   {{ font-size: 14px; font-weight: bold; margin-bottom: 4px; }}
  .subtitle {{ font-size: 10px; color: #555; margin-bottom: 16px; }}
  .cell-rect  {{ stroke: #fff; stroke-width: 0.5px; }}
  .cell-rect:hover {{ stroke: #333; stroke-width: 1.5px; cursor: pointer; }}
  .sig-star {{ font-size: 10px; fill: #222; font-weight: bold; pointer-events: none; }}
  .axis text {{ font-size: 9px; }}
  .y-label {{ font-size: 9px; }}
  .legend-title {{ font-size: 9px; fill: #333; }}
  .legend-tick  {{ font-size: 8px; fill: #555; }}
  .bar-ct {{ fill: #0072B2; opacity: 0.85; }}
  .bar-ct:hover {{ opacity: 1; }}
  #tooltip {{
    position: absolute; background: rgba(0,0,0,0.78); color: #fff;
    padding: 6px 10px; border-radius: 4px; font-size: 10px;
    pointer-events: none; display: none; max-width: 220px;
  }}
  .panel-title {{ font-size: 11px; font-weight: bold; fill: #222; }}
  .na-text {{ font-size: 7px; fill: #bbb; }}
</style>
</head>
<body>
<h2>F1 · I4 PBMC snRNA-seq: Cell-Type Pathway Response (R+1 vs. Pre-flight)</h2>
<div class="subtitle">
  fGSEA · Hallmark gene sets · 10 cell types · {len(pathways)} pathways computed (minSize=10) · ★ padj&lt;0.05
</div>
<div id="tooltip"></div>
<svg id="main"></svg>

<script>
const DATA = {json.dumps(data)};

// ── Layout ────────────────────────────────────────────────────────────────
const cellW  = 52, cellH = 18;
const marginLeft = 200, marginTop = 80, marginBottom = 60;
const nRow = DATA.pathways.length, nCol = DATA.cell_types.length;
const heatW = nCol * cellW, heatH = nRow * cellH;
const barH  = 120;
const legendW = 180, legendH = 16;

const totalW = marginLeft + heatW + legendW + 40;
const totalH = marginTop + heatH + marginBottom + barH + 60;

const svg = d3.select("#main")
  .attr("width", totalW)
  .attr("height", totalH);

const tooltip = d3.select("#tooltip");

// ── Color scale (blue–white–red, NES range ±3) ───────────────────────────
const colorScale = d3.scaleDiverging()
  .domain([-2.5, 0, 2.5])
  .interpolator(d3.interpolateRdBu)
  .clamp(true);
// Note: RdBu: low=blue, high=red → reverse for spaceflight (down=blue, up=red is standard)

// ── Panel A: Heatmap ──────────────────────────────────────────────────────
const heatG = svg.append("g")
  .attr("transform", `translate(${{marginLeft}},${{marginTop}})`);

// Column headers (cell types)
heatG.selectAll(".colLabel")
  .data(DATA.cell_types)
  .join("text")
    .attr("class", "colLabel")
    .attr("x", (d, i) => i * cellW + cellW / 2)
    .attr("y", -6)
    .attr("text-anchor", "end")
    .attr("font-size", "9px")
    .attr("transform", (d, i) => `rotate(-45, ${{i * cellW + cellW/2}}, -6)`)
    .text(d => d);

// Row labels (pathways)
heatG.selectAll(".rowLabel")
  .data(DATA.pathways)
  .join("text")
    .attr("class", "y-label")
    .attr("x", -6)
    .attr("y", (d, i) => i * cellH + cellH / 2 + 3.5)
    .attr("text-anchor", "end")
    .text(d => d);

// Cells
const cellG = heatG.selectAll(".cellG")
  .data(DATA.cells)
  .join("g")
    .attr("class", "cellG")
    .attr("transform", d => `translate(${{d.col * cellW}}, ${{d.row * cellH}})`);

cellG.append("rect")
  .attr("class", "cell-rect")
  .attr("width", cellW)
  .attr("height", cellH)
  .attr("fill", d => d.nes === null ? "#eee" : colorScale(d.nes))
  .on("mouseover", function(event, d) {{
    const pw = DATA.pathways[d.row];
    const ct = DATA.cell_types[d.col];
    const nesStr = d.nes !== null ? d.nes.toFixed(3) : "N/A";
    const padjStr = d.padj !== null ? d.padj.toExponential(2) : "N/A";
    tooltip.style("display", "block")
      .html(`<b>${{pw}}</b><br/><b>${{ct}}</b><br/>NES: ${{nesStr}}<br/>padj: ${{padjStr}}`);
  }})
  .on("mousemove", function(event) {{
    tooltip.style("left", (event.pageX + 12) + "px")
           .style("top",  (event.pageY - 20) + "px");
  }})
  .on("mouseout", () => tooltip.style("display", "none"));

// NA label
cellG.filter(d => d.nes === null)
  .append("text")
    .attr("class", "na-text")
    .attr("x", cellW / 2).attr("y", cellH / 2 + 2.5)
    .attr("text-anchor", "middle")
    .text("—");

// Significance stars
cellG.filter(d => d.sig)
  .append("text")
    .attr("class", "sig-star")
    .attr("x", cellW / 2).attr("y", cellH / 2 + 4)
    .attr("text-anchor", "middle")
    .text("★");

// Panel label
svg.append("text").attr("class", "panel-title")
  .attr("x", 8).attr("y", marginTop - 12)
  .text("A");

// ── Color legend ─────────────────────────────────────────────────────────
const legG = svg.append("g")
  .attr("transform", `translate(${{marginLeft + heatW + 20}}, ${{marginTop + 10}})`);

const legScale = d3.scaleLinear().domain([0, legendW]).range([-2.5, 2.5]);
const legGrad = svg.append("defs").append("linearGradient")
  .attr("id", "legGrad").attr("x1", "0%").attr("x2", "100%");
d3.range(0, 1.01, 0.05).forEach(t => {{
  legGrad.append("stop")
    .attr("offset", `${{t * 100}}%`)
    .attr("stop-color", colorScale(t * 5 - 2.5));
}});

legG.append("rect")
  .attr("width", legendW).attr("height", legendH)
  .style("fill", "url(#legGrad)");

legG.append("text").attr("class", "legend-title")
  .attr("y", -4).text("NES");

[-2.5, -1.25, 0, 1.25, 2.5].forEach(v => {{
  const x = (v + 2.5) / 5 * legendW;
  legG.append("line")
    .attr("x1", x).attr("x2", x).attr("y1", legendH).attr("y2", legendH + 4)
    .attr("stroke", "#666");
  legG.append("text").attr("class", "legend-tick")
    .attr("x", x).attr("y", legendH + 12).attr("text-anchor", "middle")
    .text(v.toFixed(1));
}});

// ── Panel B: Bar chart — significant pathways per cell type ───────────────
const barTop = marginTop + heatH + marginBottom + 20;
const barG = svg.append("g")
  .attr("transform", `translate(${{marginLeft}}, ${{barTop}})`);

svg.append("text").attr("class", "panel-title")
  .attr("x", 8).attr("y", barTop - 8)
  .text("B");

const ctNames = DATA.cell_types;
const sigCounts = ctNames.map(ct => DATA.sig_by_ct[ct] || 0);
const maxSig = Math.max(...sigCounts) + 1;

const xBar = d3.scaleBand().domain(ctNames).range([0, heatW]).padding(0.25);
const yBar = d3.scaleLinear().domain([0, maxSig]).range([barH, 0]);

barG.selectAll(".bar")
  .data(ctNames)
  .join("rect")
    .attr("class", "bar-ct")
    .attr("x", d => xBar(d))
    .attr("y", d => yBar(DATA.sig_by_ct[d] || 0))
    .attr("width", xBar.bandwidth())
    .attr("height", d => barH - yBar(DATA.sig_by_ct[d] || 0));

barG.selectAll(".bar-label")
  .data(ctNames)
  .join("text")
    .attr("x", d => xBar(d) + xBar.bandwidth() / 2)
    .attr("y", d => yBar(DATA.sig_by_ct[d] || 0) - 3)
    .attr("text-anchor", "middle")
    .attr("font-size", "9px")
    .text(d => (DATA.sig_by_ct[d] || 0) > 0 ? (DATA.sig_by_ct[d] || 0) : "");

// x-axis
barG.append("g").attr("transform", `translate(0,${{barH}})`)
  .call(d3.axisBottom(xBar).tickSize(3))
  .selectAll("text")
    .attr("transform", "rotate(-35)")
    .attr("text-anchor", "end")
    .attr("font-size", "9px");

// y-axis
barG.append("g").call(d3.axisLeft(yBar).ticks(maxSig).tickFormat(d3.format("d")));

barG.append("text")
  .attr("x", -barH / 2).attr("y", -40)
  .attr("transform", "rotate(-90)")
  .attr("text-anchor", "middle")
  .attr("font-size", "9px")
  .text("Sig. pathways (padj<0.05)");

barG.append("text")
  .attr("x", heatW / 2).attr("y", barH + 52)
  .attr("text-anchor", "middle")
  .attr("font-size", "10px")
  .text("Cell Type");

</script>
</body>
</html>
"""

out_path = Path("v2/figures/F1_celltype_pathway_heatmap.html")
out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text(html, encoding="utf-8")
print(f"Saved: {out_path}")
print(f"Heatmap: {len(pathways)} pathways × {len(cell_types)} cell types")
print(f"Total sig (padj<0.05): {n_sig_total}")
print("Sig by cell type:", sig_by_ct)

#!/usr/bin/env python3
"""
F1-temporal: I4 PBMC temporal pathway heatmap
Shows sig pathway counts across 5 timepoint comparisons × 10 cell types
Output: v2/figures/F1_temporal_heatmap.html
"""
import json, math
import numpy as np
import pandas as pd
from pathlib import Path

df = pd.read_csv("v2/processed/F1_scrna/i4_snrnaseq_temporal_fgsea.csv")

SHEET_ORDER = ["I4-FP1", "I4-LP1", "I4-LP3", "I4-RP1", "I4-RP3"]
SHEET_LABEL = {
    "I4-FP1": "FP1\nR+1 vs pre",
    "I4-LP1": "LP1\nR+45 vs pre",
    "I4-LP3": "LP3\nR+45/82 vs pre",
    "I4-RP1": "RP1\nR+45 vs R+1",
    "I4-RP3": "RP3\nR+45/82 vs R+1",
}
CT_ORDER = ["CD14+ Monocyte","CD16+ Monocyte","Dendritic Cell","Natural Killer Cell",
            "B Cell","CD4+ T Cell","CD8+ T Cell","Other T Cell","Other","PBMC Pseudobulk"]

# Pivot: sig pathway count per cell type × sheet
sig_pivot = (df[df["padj"] < 0.05]
             .groupby(["sheet","Cell_Type"])["pathway"]
             .count()
             .unstack(fill_value=0)
             .reindex(index=SHEET_ORDER, columns=CT_ORDER, fill_value=0))

# Also get mean |NES| for heatmap color (second panel)
nes_mean = (df.groupby(["sheet","Cell_Type"])["NES"]
              .apply(lambda x: x.abs().mean())
              .unstack(fill_value=0)
              .reindex(index=SHEET_ORDER, columns=CT_ORDER, fill_value=0))

# Build D3 data
cells_sig = []
for ri, sheet in enumerate(SHEET_ORDER):
    for ci, ct in enumerate(CT_ORDER):
        cells_sig.append({"row": ri, "col": ci,
                          "n": int(sig_pivot.loc[sheet, ct]),
                          "nes_mean": round(float(nes_mean.loc[sheet, ct]), 3)})

fig_data = {
    "sheets": SHEET_ORDER,
    "sheet_labels": [SHEET_LABEL[s].replace("\n", " | ") for s in SHEET_ORDER],
    "cell_types": CT_ORDER,
    "cells": cells_sig,
    "max_sig": int(sig_pivot.values.max()),
}

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>F1 Temporal: I4 PBMC Pathway Response</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; font-size: 12px; margin: 20px; background:#fff; }}
  h2   {{ font-size:14px; font-weight:bold; margin-bottom:4px; }}
  .subtitle {{ font-size:10px; color:#555; margin-bottom:16px; }}
  .cell-rect {{ stroke:#fff; stroke-width:1px; cursor:pointer; }}
  .cell-rect:hover {{ stroke:#333; stroke-width:1.5px; }}
  .cell-label {{ font-size:9px; fill:#fff; pointer-events:none; font-weight:bold; }}
  .panel-title {{ font-size:11px; font-weight:bold; fill:#222; }}
  #tooltip {{ position:absolute; background:rgba(0,0,0,0.78); color:#fff;
    padding:5px 9px; border-radius:4px; font-size:10px; pointer-events:none; display:none; }}
  .phase-label {{ font-size:8.5px; fill:#555; }}
</style>
</head>
<body>
<h2>F1 Temporal · I4 PBMC Cell-Type Pathway Response Across Timepoints</h2>
<div class="subtitle">
  Hallmark fGSEA · padj&lt;0.05 · 5 comparisons × 10 cell types<br>
  <b>FP1–LP3:</b> vs. pre-flight &nbsp;|&nbsp; <b>RP1–RP3:</b> vs. R+1 (recovery)
</div>
<div id="tooltip"></div>
<svg id="main"></svg>
<script>
const D = {json.dumps(fig_data)};

const cellW = 58, cellH = 36;
const mL = 180, mT = 90, mB = 30, mR = 40;
const nRow = D.sheets.length, nCol = D.cell_types.length;
const W = mL + nCol*cellW + mR, H = mT + nRow*cellH + mB + 40;

const svg = d3.select("#main").attr("width", W).attr("height", H);
const tip = d3.select("#tooltip");

const maxN = D.max_sig;
// Color: 0=white, max=deep blue
const colorN = d3.scaleSequential().domain([0, maxN]).interpolator(d3.interpolateBlues);

const g = svg.append("g").attr("transform", `translate(${{mL}},${{mT}})`);

// Column headers
g.selectAll(".colLabel").data(D.cell_types).join("text")
  .attr("x", (d,i) => i*cellW + cellW/2)
  .attr("y", -8)
  .attr("text-anchor","end")
  .attr("font-size","9px")
  .attr("transform", (d,i) => `rotate(-40, ${{i*cellW+cellW/2}}, -8)`)
  .text(d => d);

// Row labels
const rowLabels = D.sheet_labels;
g.selectAll(".rowLabel").data(rowLabels).join("text")
  .attr("x", -8).attr("y", (d,i) => i*cellH + cellH/2 + 3.5)
  .attr("text-anchor","end").attr("font-size","9px")
  .text(d => d);

// Phase brackets
const fpY = 0, rpY = 3*cellH;
g.append("text").attr("class","phase-label")
  .attr("x", -mL+4).attr("y", fpY + 1.5*cellH)
  .attr("dominant-baseline","middle").attr("font-size","9px").attr("fill","#0072B2")
  .text("← vs. pre-flight");
g.append("text").attr("class","phase-label")
  .attr("x", -mL+4).attr("y", rpY + cellH)
  .attr("dominant-baseline","middle").attr("font-size","9px").attr("fill","#E69F00")
  .text("← recovery");

// Divider line between FP/LP and RP
g.append("line")
  .attr("x1",-4).attr("x2", nCol*cellW+4)
  .attr("y1", rpY-2).attr("y2", rpY-2)
  .attr("stroke","#aaa").attr("stroke-dasharray","4,3").attr("stroke-width",1);

// Cells
const cellG = g.selectAll(".cellG").data(D.cells).join("g")
  .attr("transform", d => `translate(${{d.col*cellW}}, ${{d.row*cellH}})`);

cellG.append("rect").attr("class","cell-rect")
  .attr("width", cellW).attr("height", cellH)
  .attr("fill", d => d.n === 0 ? "#f5f5f5" : colorN(d.n))
  .on("mouseover", (ev,d) => {{
    const sheet = D.sheets[d.row], ct = D.cell_types[d.col];
    tip.style("display","block")
       .html(`<b>${{ct}}</b><br><b>${{D.sheet_labels[d.row]}}</b><br>Sig pathways: ${{d.n}}<br>Mean |NES|: ${{d.nes_mean.toFixed(2)}}`);
  }})
  .on("mousemove", ev => tip.style("left",(ev.pageX+12)+"px").style("top",(ev.pageY-20)+"px"))
  .on("mouseout", () => tip.style("display","none"));

cellG.filter(d => d.n > 0).append("text").attr("class","cell-label")
  .attr("x", cellW/2).attr("y", cellH/2 + 4)
  .attr("text-anchor","middle")
  .style("fill", d => d.n >= maxN*0.6 ? "#fff" : "#333")
  .text(d => d.n);

// Legend
const legG = svg.append("g").attr("transform", `translate(${{mL}}, ${{mT + nRow*cellH + 15}})`);
const legW = 120, legH = 12;
const legGrad = svg.append("defs").append("linearGradient").attr("id","legG2")
  .attr("x1","0%").attr("x2","100%");
d3.range(0,1.01,0.1).forEach(t => {{
  legGrad.append("stop").attr("offset",`${{t*100}}%`).attr("stop-color", colorN(t*maxN));
}});
legG.append("rect").attr("width",legW).attr("height",legH).style("fill","url(#legG2)");
legG.append("text").attr("y",-3).attr("font-size","9px").text("Sig pathways (padj<0.05)");
[0, Math.round(maxN/2), maxN].forEach(v => {{
  const x = v/maxN*legW;
  legG.append("line").attr("x1",x).attr("x2",x).attr("y1",legH).attr("y2",legH+4).attr("stroke","#666");
  legG.append("text").attr("x",x).attr("y",legH+12).attr("text-anchor","middle").attr("font-size","8px").text(v);
}});
</script>
</body>
</html>
"""

out = Path("v2/figures/F1_temporal_heatmap.html")
out.write_text(html, encoding="utf-8")
print(f"Saved: {out}")
print("\nSig counts (rows=timepoint, cols=cell type):")
print(sig_pivot.to_string())

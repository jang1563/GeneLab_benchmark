#!/usr/bin/env python
"""Fig S6: Per-tissue immune deconvolution detail.

8 panels (one per tissue): violin/box plots showing mMCP-counter scores
for each cell type, split by FLT vs GC condition.
Since we only have summary statistics (mean, p-values), we show
grouped bar charts with significance annotations.
"""
import json
import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent.parent
V5_EVAL_DIR = BASE_DIR / "v5" / "evaluation"
OUT_DIR = BASE_DIR / "v5" / "figures" / "html"

TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin", "lung", "colon"]
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastrocnemius", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin", "lung": "Lung†", "colon": "Colon†",
}


def load_all_immune():
    """Load immune deconvolution results for all tissues."""
    all_tissues = {}
    for tissue in TISSUES:
        path = V5_EVAL_DIR / f"immune_deconv_{tissue}.json"
        if not path.exists():
            continue
        with open(path) as f:
            d = json.load(f)
        cell_data = []
        for ct, info in d.get("cell_types", {}).items():
            if "note" in info:
                continue
            cell_data.append({
                "cell_type": ct,
                "mean_flt": round(info.get("mean_flt", 0), 4),
                "mean_gc": round(info.get("mean_gc", 0), 4),
                "wilcoxon_p": info.get("wilcoxon_p", 1.0),
                "fdr_p": info.get("fdr_p", info.get("wilcoxon_p", 1.0)),
                "cliffs_delta": info.get("cliffs_delta", 0),
                "direction": info.get("direction", ""),
            })
        # Sort by absolute Cliff's delta descending
        cell_data.sort(key=lambda x: abs(x.get("cliffs_delta", 0)), reverse=True)
        all_tissues[tissue] = {
            "label": TISSUE_LABELS[tissue],
            "n_flt": d.get("n_flt", 0),
            "n_gc": d.get("n_gc", 0),
            "n_sig": d.get("n_significant_fdr05", 0),
            "cells": cell_data,
        }
    return all_tissues


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    all_tissues = load_all_immune()
    total_cells = sum(len(v["cells"]) for v in all_tissues.values())
    print(f"Fig S6: {len(all_tissues)} tissues, {total_cells} total cell type entries")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig S6: Per-Tissue Immune Deconvolution Detail</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
body {{ font-family: Arial, Helvetica, sans-serif; margin: 0; padding: 20px; background: #fff; }}
.figure-container {{ width: 1400px; margin: 0 auto; }}
.figure-title {{ font-size: 14px; font-weight: bold; text-align: center; margin-bottom: 4px; }}
.figure-subtitle {{ font-size: 10px; color: #666; text-align: center; margin-bottom: 16px; }}
.panels {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
.panel {{ border: 1px solid #e0e0e0; border-radius: 4px; padding: 10px; background: #fafafa; }}
.panel-label {{ font-size: 11px; font-weight: bold; margin-bottom: 2px; }}
.panel-desc {{ font-size: 8px; color: #666; margin-bottom: 6px; }}
.tooltip {{ position: absolute; background: rgba(0,0,0,0.85); color: #fff; padding: 6px 10px;
            border-radius: 3px; font-size: 9px; pointer-events: none; z-index: 100; }}
.legend-box {{ text-align: center; margin: 12px 0 6px 0; font-size: 9px; }}
.legend-box span {{ margin: 0 12px; }}
</style>
</head>
<body>
<div class="figure-container">
  <div class="figure-title">Fig S6: Per-Tissue Immune Cell Type Scores (mMCP-counter)</div>
  <div class="figure-subtitle">FLT vs GC mean deconvolution scores across 8 tissues. *FDR &lt; 0.05, **FDR &lt; 0.01</div>
  <div class="legend-box">
    <span style="color:#D55E00;">&#9632; FLT (spaceflight)</span>
    <span style="color:#0072B2;">&#9632; GC (ground control)</span>
  </div>
  <div class="panels" id="panels"></div>
</div>
<div class="tooltip" id="tooltip" style="display:none;"></div>

<script>
const OI = {{
  black: "#000000", orange: "#E69F00", skyblue: "#56B4E9", green: "#009E73",
  yellow: "#F0E442", blue: "#0072B2", vermilion: "#D55E00", pink: "#CC79A7"
}};
const tooltip = d3.select("#tooltip");

const allData = {json.dumps(all_tissues)};
const tissueOrder = {json.dumps(TISSUES)};

tissueOrder.forEach((tissue, idx) => {{
  const tData = allData[tissue];
  if (!tData) return;

  const panelDiv = d3.select("#panels").append("div").attr("class", "panel");
  panelDiv.append("div").attr("class", "panel-label")
    .text(`${{String.fromCharCode(65 + idx)}}  ${{tData.label}} (n=${{tData.n_flt}}F + ${{tData.n_gc}}GC, ${{tData.n_sig}} sig)`);

  const cells = tData.cells;
  if (cells.length === 0) {{
    panelDiv.append("div").style("font-size", "9px").style("color", "#999").text("No cell types scored");
    return;
  }}

  const margin = {{top: 8, right: 10, bottom: 60, left: 50}};
  const barW = 25;
  const width = cells.length * (barW * 2 + 10) + 20;
  const height = 180;

  const svg = panelDiv.append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const maxVal = d3.max(cells, d => Math.max(d.mean_flt, d.mean_gc)) || 1;
  const x = d3.scaleBand().domain(cells.map(d => d.cell_type)).range([0, width]).padding(0.2);
  const y = d3.scaleLinear().domain([0, maxVal * 1.15]).range([height, 0]);

  // Grouped bars
  cells.forEach(d => {{
    const cx = x(d.cell_type);
    const bw = x.bandwidth() / 2.2;

    // FLT bar
    svg.append("rect")
      .attr("x", cx).attr("y", y(d.mean_flt))
      .attr("width", bw).attr("height", height - y(d.mean_flt))
      .attr("fill", OI.vermilion).attr("opacity", 0.85)
      .on("mouseover", (e) => {{
        tooltip.style("display", "block")
          .html(`${{d.cell_type}}<br>FLT: ${{d.mean_flt.toFixed(3)}}<br>GC: ${{d.mean_gc.toFixed(3)}}<br>FDR: ${{d.fdr_p.toFixed(4)}}<br>δ: ${{d.cliffs_delta}}`)
          .style("left", (e.pageX+10)+"px").style("top", (e.pageY-10)+"px");
      }}).on("mouseout", () => tooltip.style("display", "none"));

    // GC bar
    svg.append("rect")
      .attr("x", cx + bw + 2).attr("y", y(d.mean_gc))
      .attr("width", bw).attr("height", height - y(d.mean_gc))
      .attr("fill", OI.blue).attr("opacity", 0.85)
      .on("mouseover", (e) => {{
        tooltip.style("display", "block")
          .html(`${{d.cell_type}}<br>FLT: ${{d.mean_flt.toFixed(3)}}<br>GC: ${{d.mean_gc.toFixed(3)}}<br>FDR: ${{d.fdr_p.toFixed(4)}}<br>δ: ${{d.cliffs_delta}}`)
          .style("left", (e.pageX+10)+"px").style("top", (e.pageY-10)+"px");
      }}).on("mouseout", () => tooltip.style("display", "none"));

    // Significance stars
    if (d.fdr_p < 0.01) {{
      svg.append("text").attr("x", cx + bw).attr("y", y(Math.max(d.mean_flt, d.mean_gc)) - 4)
        .attr("text-anchor", "middle").attr("font-size", "10px").attr("fill", "#c00").text("**");
    }} else if (d.fdr_p < 0.05) {{
      svg.append("text").attr("x", cx + bw).attr("y", y(Math.max(d.mean_flt, d.mean_gc)) - 4)
        .attr("text-anchor", "middle").attr("font-size", "10px").attr("fill", "#c00").text("*");
    }}
  }});

  svg.append("g").attr("transform", `translate(0,${{height}})`).call(d3.axisBottom(x))
    .selectAll("text").attr("transform", "rotate(-45)").style("text-anchor", "end").style("font-size", "7px");
  svg.append("g").call(d3.axisLeft(y).ticks(4)).selectAll("text").style("font-size", "7px");
  svg.append("text").attr("transform", "rotate(-90)").attr("y", -35).attr("x", -height/2)
    .attr("text-anchor", "middle").attr("font-size", "8px").text("Score");
}});
</script>
</body>
</html>"""

    out_path = OUT_DIR / "FigS6_immune_detail.html"
    with open(out_path, "w") as f:
        f.write(html)
    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()

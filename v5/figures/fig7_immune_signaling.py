#!/usr/bin/env python
"""Fig 7: Immune Deconvolution + Cross-Organ Signaling.

Panel A: Immune cell type heatmap (8 tissues × cell types, log2FC FLT/GC)
Panel B: Cross-organ L-R signaling chord overview
Panel C: TF activity heatmap (top-20 TFs × tissues with sig results)
Panel D: Tissue druggability bar chart (from Phase 4)
"""
import json
import os
from pathlib import Path

import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent.parent
V5_EVAL_DIR = BASE_DIR / "v5" / "evaluation"
OUT_DIR = BASE_DIR / "v5" / "figures" / "html"

TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin", "lung", "colon"]
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastrocnemius", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin", "lung": "Lung†", "colon": "Colon†",
}


def load_immune_data():
    """Load immune deconvolution results for all tissues."""
    all_data = []
    cell_types_union = set()
    for tissue in TISSUES:
        path = V5_EVAL_DIR / f"immune_deconv_{tissue}.json"
        if not path.exists():
            continue
        with open(path) as f:
            d = json.load(f)
        for ct, info in d.get("cell_types", {}).items():
            if "note" in info:
                continue
            cell_types_union.add(ct)
            mean_flt = info.get("mean_flt", 0)
            mean_gc = info.get("mean_gc", 0)
            log2fc = np.log2((mean_flt + 0.01) / (mean_gc + 0.01))
            all_data.append({
                "tissue": TISSUE_LABELS[tissue],
                "cell_type": ct,
                "log2fc": round(float(log2fc), 3),
                "fdr_p": info.get("fdr_p", info.get("wilcoxon_p", 1.0)),
                "cliffs_delta": info.get("cliffs_delta", 0),
            })
    return all_data, sorted(cell_types_union)


def load_signaling_data():
    """Load cross-organ L-R signaling data."""
    path = V5_EVAL_DIR / "cross_organ_signaling.json"
    if not path.exists():
        return []
    with open(path) as f:
        d = json.load(f)
    edges = []
    for key, info in d.get("tissue_pairs", {}).items():
        n = info.get("n_active_pairs", 0)
        if n > 0:
            edges.append({
                "source": TISSUE_LABELS.get(info["source"], info["source"]),
                "target": TISSUE_LABELS.get(info["target"], info["target"]),
                "n_pairs": n,
            })
    return edges


def load_tf_data():
    """Load TF activity data for tissues with significant results."""
    tf_scores = {}  # tf_name → {tissue: t_stat}
    for tissue in TISSUES:
        path = V5_EVAL_DIR / f"tf_activity_{tissue}.json"
        if not path.exists():
            continue
        with open(path) as f:
            d = json.load(f)
        for tf, info in d.get("tf_results", {}).items():
            if info.get("fdr_p", 1) < 0.05:
                if tf not in tf_scores:
                    tf_scores[tf] = {}
                tf_scores[tf][TISSUE_LABELS[tissue]] = {
                    "t_stat": info.get("t_stat_approx", 0),
                    "fdr_p": info.get("fdr_p", 1),
                    "direction": info.get("direction", ""),
                }

    # Top-20 TFs by number of significant tissues
    ranked = sorted(tf_scores.items(), key=lambda x: -len(x[1]))[:20]
    heatmap = []
    for tf, tissue_data in ranked:
        for tissue_label in [TISSUE_LABELS[t] for t in TISSUES]:
            if tissue_label in tissue_data:
                heatmap.append({
                    "tf": tf,
                    "tissue": tissue_label,
                    "t_stat": round(tissue_data[tissue_label]["t_stat"], 3),
                    "fdr_p": tissue_data[tissue_label]["fdr_p"],
                })
            else:
                heatmap.append({
                    "tf": tf, "tissue": tissue_label, "t_stat": 0, "fdr_p": 1,
                })
    return heatmap, [tf for tf, _ in ranked]


def load_druggability():
    """Load per-tissue druggability from Phase 4."""
    path = V5_EVAL_DIR / "drug_targets.json"
    if not path.exists():
        return []
    with open(path) as f:
        d = json.load(f)
    bars = []
    for tissue in TISSUES:
        td = d.get("tissue_druggability", {}).get(tissue, {})
        bars.append({
            "tissue": TISSUE_LABELS[tissue],
            "total": td.get("total", 0),
            "druggable": td.get("druggable", 0),
            "pct": td.get("pct_druggable", 0),
        })
    return bars


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    immune_data, cell_types = load_immune_data()
    signaling_edges = load_signaling_data()
    tf_heatmap, tf_names = load_tf_data()
    druggability = load_druggability()

    print(f"Panel A: {len(immune_data)} cells ({len(cell_types)} cell types × 8 tissues)")
    print(f"Panel B: {len(signaling_edges)} L-R edges")
    print(f"Panel C: {len(tf_heatmap)} TF×tissue cells ({len(tf_names)} TFs)")
    print(f"Panel D: {len(druggability)} tissue bars")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig 7: Immune Deconvolution & Cross-Organ Signaling</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
body {{ font-family: Arial, Helvetica, sans-serif; margin: 0; padding: 20px; background: #fff; }}
.figure-container {{ width: 1200px; margin: 0 auto; }}
.figure-title {{ font-size: 14px; font-weight: bold; text-align: center; margin-bottom: 4px; }}
.figure-subtitle {{ font-size: 10px; color: #666; text-align: center; margin-bottom: 16px; }}
.panels {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
.panel {{ border: 1px solid #e0e0e0; border-radius: 4px; padding: 12px; background: #fafafa; }}
.panel-label {{ font-size: 12px; font-weight: bold; margin-bottom: 2px; }}
.panel-desc {{ font-size: 9px; color: #666; margin-bottom: 8px; }}
.interpretation {{ margin-top: 16px; padding: 12px; background: #f0f7ff; border-radius: 4px;
                   font-size: 10px; line-height: 1.5; }}
.tooltip {{ position: absolute; background: rgba(0,0,0,0.85); color: #fff; padding: 6px 10px;
            border-radius: 3px; font-size: 9px; pointer-events: none; z-index: 100; }}
.dagger-note {{ font-size: 8px; color: #888; margin-top: 4px; }}
</style>
</head>
<body>
<div class="figure-container">
  <div class="figure-title">Fig 7: Immune Microenvironment & Cross-Organ Signaling in Spaceflight</div>
  <div class="figure-subtitle">mMCP-counter deconvolution, OmniPath L-R pairs, and decoupler-py TF activity across 8 tissues</div>
  <div class="panels">
    <div class="panel" id="panelA">
      <div class="panel-label">A</div>
      <div class="panel-desc">Immune cell type changes (log2 FC, FLT vs GC)</div>
      <div id="immune-heatmap"></div>
      <div class="dagger-note">*FDR &lt; 0.05 &nbsp; **FDR &lt; 0.01 &nbsp; †Single-mission (5-fold CV)</div>
    </div>
    <div class="panel" id="panelB">
      <div class="panel-label">B</div>
      <div class="panel-desc">Cross-organ L-R signaling (SHAP-filtered active pairs)</div>
      <div id="signaling-net"></div>
    </div>
    <div class="panel" id="panelC">
      <div class="panel-label">C</div>
      <div class="panel-desc">TF activity (FLT − GC, top-20 by cross-tissue significance)</div>
      <div id="tf-heatmap"></div>
    </div>
    <div class="panel" id="panelD">
      <div class="panel-label">D</div>
      <div class="panel-desc">Druggability of spaceflight-affected genes per tissue</div>
      <div id="druggability-bar"></div>
    </div>
  </div>
  <div class="interpretation">
    <strong>Key findings:</strong>
    Skin shows the strongest immune remodeling (6/14 cell types FDR&lt;0.05), followed by kidney and thymus (2 each).
    Cross-organ L-R analysis identifies 62 active ligand-receptor pairs linking tissues during spaceflight.
    TF activity reveals widespread transcriptional reprogramming: thymus (240 TFs), skin (241), kidney (177), and liver (105) show significant FLT vs GC differences.
    Thymus and kidney are the most druggable tissues (24.8% and 21.5% of spaceflight-affected genes targeted by known drugs).
  </div>
</div>
<div class="tooltip" id="tooltip" style="display:none;"></div>

<script>
const OI = {{
  black: "#000000", orange: "#E69F00", skyblue: "#56B4E9", green: "#009E73",
  yellow: "#F0E442", blue: "#0072B2", vermilion: "#D55E00", pink: "#CC79A7"
}};
const tissueColors = {{
  "Liver": OI.blue, "Gastrocnemius": OI.vermilion, "Kidney": OI.green,
  "Thymus": OI.orange, "Eye": OI.pink, "Skin": OI.skyblue,
  "Lung†": "#888888", "Colon†": "#555555"
}};
const tooltip = d3.select("#tooltip");

// Panel A: Immune heatmap
(function() {{
  const data = {json.dumps(immune_data)};
  const tissues = {json.dumps([TISSUE_LABELS[t] for t in TISSUES])};
  const cellTypes = {json.dumps(cell_types)};

  const margin = {{top: 10, right: 60, bottom: 80, left: 130}};
  const cellW = 50, cellH = 22;
  const width = tissues.length * cellW, height = cellTypes.length * cellH;

  const svg = d3.select("#immune-heatmap").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(tissues).range([0, width]).padding(0.05);
  const y = d3.scaleBand().domain(cellTypes).range([0, height]).padding(0.05);
  const color = d3.scaleDiverging(d3.interpolateRdBu).domain([1, 0, -1]);

  svg.selectAll("rect").data(data).join("rect")
    .attr("x", d => x(d.tissue)).attr("y", d => y(d.cell_type))
    .attr("width", x.bandwidth()).attr("height", y.bandwidth())
    .attr("fill", d => color(d.log2fc))
    .attr("stroke", d => d.fdr_p < 0.05 ? "#000" : "none")
    .attr("stroke-width", d => d.fdr_p < 0.01 ? 2 : 1)
    .on("mouseover", (e, d) => {{
      tooltip.style("display", "block")
        .html(`${{d.tissue}} — ${{d.cell_type}}<br>log2FC: ${{d.log2fc}}<br>FDR p: ${{d.fdr_p.toFixed(4)}}<br>Cliff's δ: ${{d.cliffs_delta}}`)
        .style("left", (e.pageX+10)+"px").style("top", (e.pageY-10)+"px");
    }}).on("mouseout", () => tooltip.style("display", "none"));

  svg.append("g").attr("transform", `translate(0,${{height}})`).call(d3.axisBottom(x))
    .selectAll("text").attr("transform", "rotate(-45)").style("text-anchor", "end").style("font-size", "8px");
  svg.append("g").call(d3.axisLeft(y)).selectAll("text").style("font-size", "8px");
}})();

// Panel B: Signaling network
(function() {{
  const edges = {json.dumps(signaling_edges)};
  const tissues = {json.dumps([TISSUE_LABELS[t] for t in TISSUES])};

  const width = 500, height = 350;
  const svg = d3.select("#signaling-net").append("svg")
    .attr("width", width).attr("height", height);

  const nodeRadius = 25;
  const cx = width / 2, cy = height / 2, r = 130;
  const nodes = tissues.map((t, i) => ({{
    name: t, x: cx + r * Math.cos(2 * Math.PI * i / tissues.length - Math.PI/2),
    y: cy + r * Math.sin(2 * Math.PI * i / tissues.length - Math.PI/2)
  }}));
  const nodeMap = Object.fromEntries(nodes.map(n => [n.name, n]));

  const maxPairs = d3.max(edges, d => d.n_pairs) || 1;

  svg.selectAll("line").data(edges).join("line")
    .attr("x1", d => nodeMap[d.source]?.x || 0).attr("y1", d => nodeMap[d.source]?.y || 0)
    .attr("x2", d => nodeMap[d.target]?.x || 0).attr("y2", d => nodeMap[d.target]?.y || 0)
    .attr("stroke", "#999").attr("stroke-opacity", 0.6)
    .attr("stroke-width", d => Math.max(1, d.n_pairs / maxPairs * 6));

  const g = svg.selectAll("g.node").data(nodes).join("g").attr("class", "node")
    .attr("transform", d => `translate(${{d.x}},${{d.y}})`);
  g.append("circle").attr("r", nodeRadius)
    .attr("fill", d => tissueColors[d.name] || "#ccc").attr("stroke", "#333").attr("stroke-width", 1);
  g.append("text").attr("dy", nodeRadius + 14).attr("text-anchor", "middle")
    .attr("font-size", "9px").text(d => d.name.replace("†", ""));

  svg.append("text").attr("x", width/2).attr("y", height - 5)
    .attr("text-anchor", "middle").attr("font-size", "8px").attr("fill", "#666")
    .text(`${{edges.length}} directed edges, edge width ∝ active L-R pairs`);
}})();

// Panel C: TF activity heatmap
(function() {{
  const data = {json.dumps(tf_heatmap)};
  const tissues = {json.dumps([TISSUE_LABELS[t] for t in TISSUES])};
  const tfs = {json.dumps(tf_names)};

  const margin = {{top: 10, right: 20, bottom: 80, left: 80}};
  const cellW = 50, cellH = 18;
  const width = tissues.length * cellW, height = tfs.length * cellH;

  const svg = d3.select("#tf-heatmap").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(tissues).range([0, width]).padding(0.05);
  const y = d3.scaleBand().domain(tfs).range([0, height]).padding(0.05);
  const maxAbs = d3.max(data, d => Math.abs(d.t_stat)) || 1;
  const color = d3.scaleDiverging(d3.interpolateRdBu).domain([maxAbs, 0, -maxAbs]);

  svg.selectAll("rect").data(data).join("rect")
    .attr("x", d => x(d.tissue)).attr("y", d => y(d.tf))
    .attr("width", x.bandwidth()).attr("height", y.bandwidth())
    .attr("fill", d => d.fdr_p < 0.05 ? color(d.t_stat) : "#f5f5f5")
    .attr("stroke", d => d.fdr_p < 0.05 ? "#333" : "#eee")
    .attr("stroke-width", 0.5)
    .on("mouseover", (e, d) => {{
      tooltip.style("display", "block")
        .html(`${{d.tf}} in ${{d.tissue}}<br>t-stat: ${{d.t_stat}}<br>FDR p: ${{d.fdr_p.toFixed(4)}}`)
        .style("left", (e.pageX+10)+"px").style("top", (e.pageY-10)+"px");
    }}).on("mouseout", () => tooltip.style("display", "none"));

  svg.append("g").attr("transform", `translate(0,${{height}})`).call(d3.axisBottom(x))
    .selectAll("text").attr("transform", "rotate(-45)").style("text-anchor", "end").style("font-size", "8px");
  svg.append("g").call(d3.axisLeft(y)).selectAll("text").style("font-size", "7px");
}})();

// Panel D: Druggability bars
(function() {{
  const data = {json.dumps(druggability)};

  const margin = {{top: 10, right: 40, bottom: 50, left: 120}};
  const width = 380, height = data.length * 35;

  const svg = d3.select("#druggability-bar").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleLinear().domain([0, 30]).range([0, width]);
  const y = d3.scaleBand().domain(data.map(d => d.tissue)).range([0, height]).padding(0.3);

  svg.selectAll("rect").data(data).join("rect")
    .attr("x", 0).attr("y", d => y(d.tissue))
    .attr("width", d => x(d.pct)).attr("height", y.bandwidth())
    .attr("fill", d => tissueColors[d.tissue] || "#ccc");

  svg.selectAll("text.pct").data(data).join("text").attr("class", "pct")
    .attr("x", d => x(d.pct) + 4).attr("y", d => y(d.tissue) + y.bandwidth()/2 + 4)
    .attr("font-size", "9px").text(d => `${{d.pct}}% (${{d.druggable}}/${{d.total}})`);

  svg.append("g").call(d3.axisLeft(y)).selectAll("text").style("font-size", "9px");
  svg.append("g").attr("transform", `translate(0,${{height}})`).call(d3.axisBottom(x).ticks(5));
  svg.append("text").attr("x", width/2).attr("y", height + 35)
    .attr("text-anchor", "middle").attr("font-size", "9px").text("% druggable genes");
}})();
</script>
</body>
</html>"""

    out_path = OUT_DIR / "Fig7_immune_signaling.html"
    with open(out_path, "w") as f:
        f.write(html)
    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()

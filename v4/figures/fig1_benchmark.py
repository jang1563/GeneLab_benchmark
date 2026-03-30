#!/usr/bin/env python3
"""
Fig 1 — Multi-Method Benchmark Overview
GeneLabBench v4 Publication Figure

Panels:
  A: 8-tissue × 8-method AUROC heatmap (gene features)
  B: Forest plot (RE AUROC ± 95% CI) for PCA-LR gene, 6 LOMO tissues
  C: Gene vs best-pathway AUROC comparison scatter (PCA-LR, 8 tissues)
  D: Nemenyi critical difference diagram (8 methods, all-8 Friedman)

Output: v4/figures/html/Fig1_benchmark.html (self-contained D3.js v7)
"""

import json
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EVAL_DIR = os.path.join(BASE, "evaluation")
OUT_DIR = os.path.join(BASE, "figures", "html")
os.makedirs(OUT_DIR, exist_ok=True)

# --- Load data ---
with open(os.path.join(EVAL_DIR, "M1_summary.json")) as f:
    m1 = json.load(f)
with open(os.path.join(EVAL_DIR, "META_forest_plots.json")) as f:
    meta_forest = json.load(f)
with open(os.path.join(EVAL_DIR, "META_method_ranking.json")) as f:
    meta_rank = json.load(f)

# --- Constants ---
METHODS_ORDER = ["pca_lr", "elasticnet_lr", "xgb", "rf", "knn", "mlp", "svm_rbf", "tabnet"]
METHOD_LABELS = {
    "pca_lr": "PCA-LR", "elasticnet_lr": "ElasticNet-LR", "rf": "Random Forest",
    "xgb": "XGBoost", "svm_rbf": "SVM-RBF", "knn": "kNN", "mlp": "MLP", "tabnet": "TabNet"
}
LOMO_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin"]
ALL_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin", "lung", "colon"]
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastrocnemius", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin", "lung": "Lung†", "colon": "Colon†"
}

# Sort tissues by PCA-LR gene AUROC descending for heatmap
tissue_order = sorted(ALL_TISSUES, key=lambda t: m1[t]["gene"]["pca_lr"]["auroc"], reverse=True)

# --- Panel A: Heatmap data (8 tissues × 8 methods, gene features) ---
heatmap_data = []
for tissue in tissue_order:
    for method in METHODS_ORDER:
        entry = m1[tissue]["gene"][method]
        heatmap_data.append({
            "tissue": TISSUE_LABELS[tissue],
            "method": METHOD_LABELS[method],
            "auroc": round(entry["auroc"], 3),
            "perm_p": entry["perm_p"],
            "cv": entry["cv"]
        })

# --- Panel B: Forest plot data (PCA-LR gene, 6 LOMO tissues) ---
forest_data = []
for fp in meta_forest["forest_plots"]:
    if fp["method"] == "pca_lr" and fp["feature_type"] == "gene" and fp["tissue"] in LOMO_TISSUES:
        forest_data.append({
            "tissue": TISSUE_LABELS[fp["tissue"]],
            "tissue_key": fp["tissue"],
            "pooled_re": round(fp["pooled_re"]["auroc"], 3),
            "ci_lower_re": round(fp["pooled_re"]["ci_lower"], 3),
            "ci_upper_re": round(fp["pooled_re"]["ci_upper"], 3),
            "I2": round(fp["heterogeneity"]["I2"], 1),
            "Q_pvalue": round(fp["heterogeneity"]["Q_pvalue"], 3),
            "studies": [{
                "fold": s["fold"],
                "auroc": round(s["auroc"], 3),
                "ci_lower": round(s["ci_lower"], 3),
                "ci_upper": round(s["ci_upper"], 3),
                "weight_re": round(s["weight_re"], 3)
            } for s in fp["studies"]]
        })

# Sort forest data by pooled RE AUROC
forest_data.sort(key=lambda x: x["pooled_re"], reverse=True)

# --- Panel C: Gene vs Pathway scatter (PCA-LR, all 8 tissues) ---
scatter_data = []
for tissue in ALL_TISSUES:
    gene_auroc = m1[tissue]["gene"]["pca_lr"]["auroc"]
    # Best pathway = max of hallmark and kegg
    hallmark_auroc = m1[tissue]["pathway_hallmark"]["pca_lr"]["auroc"]
    kegg_auroc = m1[tissue]["pathway_kegg"]["pca_lr"]["auroc"]
    best_pw = max(hallmark_auroc, kegg_auroc)
    best_pw_type = "KEGG" if kegg_auroc > hallmark_auroc else "Hallmark"
    scatter_data.append({
        "tissue": TISSUE_LABELS[tissue],
        "gene": round(gene_auroc, 3),
        "pathway": round(best_pw, 3),
        "pw_type": best_pw_type,
        "cv": m1[tissue]["gene"]["pca_lr"]["cv"]
    })

# --- Panel D: Nemenyi CD diagram data (use all-8 Friedman) ---
friedman = meta_rank["friedman_all8"]
cd_data = {
    "chi2": round(friedman["chi2"], 2),
    "p_value": friedman["p_value"],
    "cd": round(friedman["critical_difference"], 3),
    "k": friedman["k_methods"],
    "n": friedman["n_tissues"],
    "ranks": {METHOD_LABELS[m]: round(r, 3) for m, r in friedman["avg_ranks"].items()},
    "significant_pairs": [[METHOD_LABELS[a], METHOD_LABELS[b]] for a, b in friedman["significant_pairs"]],
    "groups": friedman["nemenyi_groups"]
}

# Sort methods by average rank for CD diagram
cd_methods_sorted = sorted(cd_data["ranks"].items(), key=lambda x: x[1])

# --- Build HTML ---
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig 1 — Multi-Method Benchmark for Spaceflight Transcriptomics Classification</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{
    font-family: Arial, Helvetica, sans-serif;
    margin: 0; padding: 20px; background: #fff;
  }}
  .figure-container {{ width: 1200px; margin: 0 auto; }}
  .figure-title {{ font-size: 14px; font-weight: bold; text-align: center; margin-bottom: 4px; }}
  .figure-subtitle {{ font-size: 10px; color: #666; text-align: center; margin-bottom: 16px; }}
  .panels {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
  .panel {{ border: 1px solid #e0e0e0; border-radius: 4px; padding: 12px; background: #fafafa; }}
  .panel-label {{ font-size: 12px; font-weight: bold; margin-bottom: 2px; }}
  .panel-desc {{ font-size: 9px; color: #666; margin-bottom: 8px; }}
  .axis text {{ font-size: 9px; }}
  .axis-label {{ font-size: 10px; }}
  .stat-text {{ font-size: 8px; fill: #444; }}
  .annotation {{ font-size: 8px; fill: #888; }}
  .tooltip {{
    position: absolute; padding: 6px 10px; background: rgba(0,0,0,0.85);
    color: #fff; border-radius: 4px; font-size: 10px; pointer-events: none;
    z-index: 100; line-height: 1.4;
  }}
  .interpretation {{
    margin-top: 16px; padding: 10px; background: #f5f5f5;
    border-left: 3px solid #999; font-size: 9px; color: #444; line-height: 1.5;
  }}
  .interpretation strong {{ color: #333; }}
  .dagger-note {{ font-size: 8px; color: #888; margin-top: 6px; }}
</style>
</head>
<body>
<div class="figure-container">
  <div class="figure-title">Fig 1 — Multi-Method Benchmark for Spaceflight Transcriptomics Classification</div>
  <div class="figure-subtitle">8 tissues × 8 methods × 4 feature types = 256 evaluations | Leave-One-Mission-Out CV (6 tissues) + 5-fold CV (2 tissues)</div>
  <div class="panels">
    <div class="panel" id="panelA">
      <div class="panel-label">A</div>
      <div class="panel-desc">AUROC Heatmap — Gene Features (8 tissues × 8 methods)</div>
      <div id="heatmap"></div>
      <div class="dagger-note">† Lung, Colon: 5-fold stratified CV (single mission); not directly comparable to LOMO-CV</div>
    </div>
    <div class="panel" id="panelB">
      <div class="panel-label">B</div>
      <div class="panel-desc">Random-Effects Forest Plot — PCA-LR, Gene Features (6 LOMO tissues)</div>
      <div id="forest"></div>
    </div>
    <div class="panel" id="panelC">
      <div class="panel-label">C</div>
      <div class="panel-desc">Gene vs Best-Pathway AUROC — PCA-LR (8 tissues)</div>
      <div id="scatter"></div>
    </div>
    <div class="panel" id="panelD">
      <div class="panel-label">D</div>
      <div class="panel-desc">Nemenyi Critical Difference Diagram — Average Ranks (8 methods, 8 tissues)</div>
      <div id="cd_diagram"></div>
    </div>
  </div>
  <div class="interpretation">
    <strong>Key findings:</strong>
    PCA-LR achieves the best overall performance (avg rank 2.1, gene AUROC 0.776 across 8 tissues),
    followed by ElasticNet-LR (rank 2.9, 0.762). SVM-RBF (rank 6.5) and TabNet (rank 7.3) perform worst
    — deep learning does not help with small-n spaceflight datasets.
    Friedman test confirms significant method differences (χ²={cd_data["chi2"]}, p={cd_data["p_value"]:.4f}).
    Pathway features (Hallmark/KEGG) improve classification in kidney (+0.245), eye (+0.126), and thymus (+0.040),
    while gene features are superior for skin and lung.
    Substantial between-mission heterogeneity exists (liver I²=51.9%, thymus I²=82.5%).
  </div>
</div>

<div class="tooltip" id="tooltip" style="display:none;"></div>

<script>
// === DATA ===
const heatmapData = {json.dumps(heatmap_data)};
const forestData = {json.dumps(forest_data)};
const scatterData = {json.dumps(scatter_data)};
const cdData = {json.dumps(cd_data)};
const cdMethodsSorted = {json.dumps(cd_methods_sorted)};

// === Okabe-Ito palette ===
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

// =============================================
// PANEL A: AUROC Heatmap
// =============================================
(function() {{
  const margin = {{top: 5, right: 80, bottom: 60, left: 100}};
  const cellW = 58, cellH = 32;
  const tissues = [...new Set(heatmapData.map(d => d.tissue))];
  const methods = [...new Set(heatmapData.map(d => d.method))];
  const width = methods.length * cellW;
  const height = tissues.length * cellH;

  const svg = d3.select("#heatmap").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(methods).range([0, width]).padding(0.04);
  const y = d3.scaleBand().domain(tissues).range([0, height]).padding(0.04);

  // Color scale: 0.3=light, 1.0=dark blue
  const color = d3.scaleSequential()
    .domain([0.30, 1.0])
    .interpolator(d3.interpolateBlues);

  // Cells
  svg.selectAll("rect.cell")
    .data(heatmapData).enter().append("rect")
    .attr("class", "cell")
    .attr("x", d => x(d.method))
    .attr("y", d => y(d.tissue))
    .attr("width", x.bandwidth())
    .attr("height", y.bandwidth())
    .attr("fill", d => color(Math.max(0.30, Math.min(1.0, d.auroc))))
    .attr("stroke", "#fff").attr("stroke-width", 1)
    .attr("rx", 2)
    .on("mouseover", function(event, d) {{
      tooltip.style("display", "block")
        .html(`<b>${{d.tissue}}</b> × <b>${{d.method}}</b><br>AUROC: ${{d.auroc.toFixed(3)}}<br>p_perm: ${{d.perm_p.toFixed(3)}}<br>CV: ${{d.cv}}`);
    }})
    .on("mousemove", function(event) {{
      tooltip.style("left", (event.pageX + 12) + "px").style("top", (event.pageY - 10) + "px");
    }})
    .on("mouseout", function() {{ tooltip.style("display", "none"); }});

  // AUROC text in cells
  svg.selectAll("text.cell-val")
    .data(heatmapData).enter().append("text")
    .attr("class", "cell-val")
    .attr("x", d => x(d.method) + x.bandwidth() / 2)
    .attr("y", d => y(d.tissue) + y.bandwidth() / 2)
    .attr("text-anchor", "middle").attr("dominant-baseline", "central")
    .attr("font-size", "8px").attr("font-weight", "bold")
    .attr("fill", d => d.auroc > 0.7 ? "#fff" : "#333")
    .text(d => {{
      let t = d.auroc.toFixed(2);
      if (d.perm_p < 0.01) t += "**";
      else if (d.perm_p < 0.05) t += "*";
      return t;
    }});

  // Axes
  svg.append("g").attr("class", "axis")
    .attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).tickSize(0))
    .selectAll("text")
    .attr("transform", "rotate(-40)").style("text-anchor", "end")
    .attr("dx", "-0.5em").attr("dy", "0.3em").style("font-size", "9px");

  svg.append("g").attr("class", "axis")
    .call(d3.axisLeft(y).tickSize(0))
    .selectAll("text").style("font-size", "9px");

  // Color legend
  const legendW = 12, legendH = height;
  const legendX = width + 15;
  const defs = svg.append("defs");
  const gradient = defs.append("linearGradient")
    .attr("id", "heatGrad").attr("x1", "0%").attr("y1", "100%").attr("x2", "0%").attr("y2", "0%");
  gradient.append("stop").attr("offset", "0%").attr("stop-color", color(0.30));
  gradient.append("stop").attr("offset", "100%").attr("stop-color", color(1.0));

  svg.append("rect").attr("x", legendX).attr("y", 0)
    .attr("width", legendW).attr("height", legendH)
    .style("fill", "url(#heatGrad)").attr("stroke", "#ccc");

  const legendScale = d3.scaleLinear().domain([0.30, 1.0]).range([legendH, 0]);
  svg.append("g").attr("transform", `translate(${{legendX + legendW}},0)`)
    .call(d3.axisRight(legendScale).ticks(5).tickFormat(d3.format(".1f")))
    .selectAll("text").style("font-size", "8px");

  svg.append("text").attr("x", legendX + legendW / 2).attr("y", -4)
    .attr("text-anchor", "middle").attr("font-size", "8px").text("AUROC");

  // Significance note
  svg.append("text").attr("x", 0).attr("y", height + 52)
    .attr("font-size", "8px").attr("fill", "#666")
    .text("*p < 0.05  **p < 0.01 (permutation test, n=1000)");
}})();

// =============================================
// PANEL B: Forest Plot
// =============================================
(function() {{
  const margin = {{top: 10, right: 90, bottom: 35, left: 100}};
  const width = 380, rowH = 22;

  // Flatten: one row per study + one row per pooled
  let rows = [];
  forestData.forEach(td => {{
    td.studies.forEach(s => {{
      rows.push({{ tissue: td.tissue, fold: s.fold, auroc: s.auroc,
                   ci_lo: s.ci_lower, ci_hi: s.ci_upper, weight: s.weight_re, type: "study" }});
    }});
    rows.push({{ tissue: td.tissue, fold: "Pooled (RE)", auroc: td.pooled_re,
                 ci_lo: td.ci_lower_re, ci_hi: td.ci_upper_re, I2: td.I2, Q_p: td.Q_pvalue, type: "pooled" }});
    rows.push({{ tissue: "", fold: "", auroc: null, type: "spacer" }});
  }});

  const height = rows.length * rowH;
  const svg = d3.select("#forest").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleLinear().domain([0, 1.1]).range([0, width]);

  // Reference line at 0.5
  svg.append("line").attr("x1", x(0.5)).attr("x2", x(0.5))
    .attr("y1", 0).attr("y2", height)
    .attr("stroke", "#ccc").attr("stroke-dasharray", "3,3");

  rows.forEach((row, i) => {{
    const cy = i * rowH + rowH / 2;
    if (row.type === "spacer") return;

    // Label
    if (row.type === "pooled") {{
      svg.append("text").attr("x", -5).attr("y", cy + 4)
        .attr("text-anchor", "end").attr("font-size", "9px").attr("font-weight", "bold")
        .text(row.fold);
    }} else {{
      svg.append("text").attr("x", -5).attr("y", cy + 4)
        .attr("text-anchor", "end").attr("font-size", "8px").attr("fill", "#555")
        .text(row.fold);
    }}

    // Tissue label (only on first study row)
    if (row.type === "study") {{
      const tissueRows = rows.filter(r => r.tissue === row.tissue && r.type !== "spacer");
      if (tissueRows[0] === row) {{
        svg.append("text").attr("x", -75).attr("y", cy + 4)
          .attr("text-anchor", "start").attr("font-size", "9px").attr("font-weight", "bold")
          .attr("fill", tissueColors[row.tissue] || "#333")
          .text(row.tissue);
      }}
    }}

    if (row.auroc === null) return;

    // CI line
    const ciLo = Math.max(0, row.ci_lo);
    const ciHi = Math.min(1.1, row.ci_hi);
    svg.append("line")
      .attr("x1", x(ciLo)).attr("x2", x(ciHi))
      .attr("y1", cy).attr("y2", cy)
      .attr("stroke", row.type === "pooled" ? "#333" : "#888")
      .attr("stroke-width", row.type === "pooled" ? 2 : 1);

    if (row.type === "pooled") {{
      // Diamond for pooled
      const dw = 6, dh = 8;
      const cx = x(row.auroc);
      svg.append("polygon")
        .attr("points", `${{cx}},${{cy-dh}} ${{cx+dw}},${{cy}} ${{cx}},${{cy+dh}} ${{cx-dw}},${{cy}}`)
        .attr("fill", "#333");
      // I² annotation
      svg.append("text").attr("x", width + 5).attr("y", cy + 4)
        .attr("font-size", "8px").attr("fill", "#444")
        .text(`I²=${{row.I2}}%`);
    }} else {{
      // Square for study, size ~ weight
      const size = Math.max(3, Math.sqrt(row.weight) * 20);
      svg.append("rect")
        .attr("x", x(row.auroc) - size/2).attr("y", cy - size/2)
        .attr("width", size).attr("height", size)
        .attr("fill", "#555");
    }}

    // AUROC value on right
    svg.append("text").attr("x", width + 50).attr("y", cy + 4)
      .attr("text-anchor", "end").attr("font-size", "8px")
      .attr("fill", row.type === "pooled" ? "#000" : "#666")
      .attr("font-weight", row.type === "pooled" ? "bold" : "normal")
      .text(`${{row.auroc.toFixed(2)}} [${{row.ci_lo.toFixed(2)}}, ${{row.ci_hi.toFixed(2)}}]`);
  }});

  // X axis
  svg.append("g").attr("class", "axis").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).ticks(6).tickFormat(d3.format(".1f")));
  svg.append("text").attr("x", width/2).attr("y", height + 30)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("AUROC");
}})();

// =============================================
// PANEL C: Gene vs Pathway Scatter
// =============================================
(function() {{
  const margin = {{top: 15, right: 20, bottom: 45, left: 50}};
  const width = 350, height = 300;

  const svg = d3.select("#scatter").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const xScale = d3.scaleLinear().domain([0.4, 1.0]).range([0, width]);
  const yScale = d3.scaleLinear().domain([0.4, 1.0]).range([height, 0]);

  // Diagonal (equivalence line)
  svg.append("line")
    .attr("x1", xScale(0.4)).attr("y1", yScale(0.4))
    .attr("x2", xScale(1.0)).attr("y2", yScale(1.0))
    .attr("stroke", "#ccc").attr("stroke-dasharray", "5,5").attr("stroke-width", 1);

  // Label regions
  svg.append("text").attr("x", xScale(0.88)).attr("y", yScale(0.76))
    .attr("font-size", "8px").attr("fill", "#aaa").attr("transform", "rotate(-45," + xScale(0.88) + "," + yScale(0.76) + ")")
    .text("Gene better →");
  svg.append("text").attr("x", xScale(0.55)).attr("y", yScale(0.68))
    .attr("font-size", "8px").attr("fill", "#aaa").attr("transform", "rotate(-45," + xScale(0.55) + "," + yScale(0.68) + ")")
    .text("← Pathway better");

  // Points
  svg.selectAll("circle.point")
    .data(scatterData).enter().append("circle")
    .attr("cx", d => xScale(d.gene))
    .attr("cy", d => yScale(d.pathway))
    .attr("r", 7)
    .attr("fill", d => tissueColors[d.tissue] || "#333")
    .attr("stroke", d => d.cv === "5fold_stratified" ? "#333" : "none")
    .attr("stroke-width", d => d.cv === "5fold_stratified" ? 1.5 : 0)
    .attr("stroke-dasharray", d => d.cv === "5fold_stratified" ? "2,2" : "none")
    .attr("opacity", 0.85)
    .on("mouseover", function(event, d) {{
      tooltip.style("display", "block")
        .html(`<b>${{d.tissue}}</b><br>Gene: ${{d.gene.toFixed(3)}}<br>Best Pathway (${{d.pw_type}}): ${{d.pathway.toFixed(3)}}<br>Δ: ${{(d.pathway - d.gene).toFixed(3)}}<br>CV: ${{d.cv}}`);
    }})
    .on("mousemove", function(event) {{
      tooltip.style("left", (event.pageX + 12) + "px").style("top", (event.pageY - 10) + "px");
    }})
    .on("mouseout", function() {{ tooltip.style("display", "none"); }});

  // Labels
  scatterData.forEach(d => {{
    const labelOffsets = {{
      "Liver": [8, -8], "Gastrocnemius": [8, 12], "Kidney": [8, -8],
      "Thymus": [-10, -10], "Eye": [8, 10], "Skin": [8, -8],
      "Lung†": [8, -8], "Colon†": [8, 10]
    }};
    const off = labelOffsets[d.tissue] || [8, -8];
    svg.append("text")
      .attr("x", xScale(d.gene) + off[0])
      .attr("y", yScale(d.pathway) + off[1])
      .attr("font-size", "8px").attr("fill", tissueColors[d.tissue] || "#333")
      .text(d.tissue);
  }});

  // Axes
  svg.append("g").attr("class", "axis").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(xScale).ticks(6).tickFormat(d3.format(".2f")));
  svg.append("text").attr("x", width/2).attr("y", height + 35)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("Gene AUROC (PCA-LR)");

  svg.append("g").attr("class", "axis")
    .call(d3.axisLeft(yScale).ticks(6).tickFormat(d3.format(".2f")));
  svg.append("text").attr("x", -height/2).attr("y", -38)
    .attr("text-anchor", "middle").attr("font-size", "10px")
    .attr("transform", "rotate(-90)").text("Best Pathway AUROC (PCA-LR)");

  // Legend for dashed border = 5-fold
  svg.append("circle").attr("cx", width - 100).attr("cy", 10).attr("r", 5)
    .attr("fill", "#ccc").attr("stroke", "#333").attr("stroke-width", 1.5).attr("stroke-dasharray", "2,2");
  svg.append("text").attr("x", width - 90).attr("y", 14)
    .attr("font-size", "8px").attr("fill", "#555").text("5-fold CV");
  svg.append("circle").attr("cx", width - 100).attr("cy", 26).attr("r", 5)
    .attr("fill", "#ccc");
  svg.append("text").attr("x", width - 90).attr("y", 30)
    .attr("font-size", "8px").attr("fill", "#555").text("LOMO-CV");
}})();

// =============================================
// PANEL D: Nemenyi Critical Difference Diagram
// =============================================
(function() {{
  const margin = {{top: 30, right: 30, bottom: 80, left: 30}};
  const width = 500, height = 140;

  const svg = d3.select("#cd_diagram").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const minRank = 1, maxRank = 8;
  const x = d3.scaleLinear().domain([minRank, maxRank]).range([0, width]);

  // Axis at top
  svg.append("g").attr("class", "axis")
    .call(d3.axisTop(x).ticks(8).tickFormat(d => d));
  svg.append("text").attr("x", width/2).attr("y", -22)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("Average Rank");

  // CD bar
  const cdWidth = x(minRank + cdData.cd) - x(minRank);
  svg.append("rect").attr("x", x(minRank) - 2).attr("y", -16)
    .attr("width", cdWidth).attr("height", 6)
    .attr("fill", OI.vermilion).attr("opacity", 0.7);
  svg.append("text").attr("x", x(minRank) + cdWidth / 2).attr("y", -18)
    .attr("text-anchor", "middle").attr("font-size", "8px").attr("fill", OI.vermilion)
    .text(`CD = ${{cdData.cd}}`);

  // Methods: place on left and right sides
  const methods = cdMethodsSorted;
  const leftMethods = methods.slice(0, 4);   // best 4
  const rightMethods = methods.slice(4, 8);  // worst 4

  const leftY0 = 30, rightY0 = 30, rowGap = 22;

  // Draw left methods (best ranks, labels on left)
  leftMethods.forEach((m, i) => {{
    const rank = m[1];
    const yPos = leftY0 + i * rowGap;
    // Horizontal connector from rank position to left
    svg.append("line").attr("x1", 0).attr("x2", x(rank))
      .attr("y1", yPos).attr("y2", yPos)
      .attr("stroke", "#333").attr("stroke-width", 1);
    // Tick down from axis
    svg.append("line").attr("x1", x(rank)).attr("x2", x(rank))
      .attr("y1", 0).attr("y2", yPos)
      .attr("stroke", "#333").attr("stroke-width", 1);
    // Label
    svg.append("text").attr("x", -5).attr("y", yPos + 4)
      .attr("text-anchor", "end").attr("font-size", "9px")
      .text(`${{m[0]}} (${{rank.toFixed(1)}})`);
  }});

  // Draw right methods (worst ranks, labels on right)
  rightMethods.forEach((m, i) => {{
    const rank = m[1];
    const yPos = rightY0 + i * rowGap;
    // Horizontal connector from rank position to right
    svg.append("line").attr("x1", x(rank)).attr("x2", width)
      .attr("y1", yPos).attr("y2", yPos)
      .attr("stroke", "#333").attr("stroke-width", 1);
    // Tick down from axis
    svg.append("line").attr("x1", x(rank)).attr("x2", x(rank))
      .attr("y1", 0).attr("y2", yPos)
      .attr("stroke", "#333").attr("stroke-width", 1);
    // Label
    svg.append("text").attr("x", width + 5).attr("y", yPos + 4)
      .attr("text-anchor", "start").attr("font-size", "9px")
      .text(`${{m[0]}} (${{rank.toFixed(1)}})`);
  }});

  // Draw Nemenyi groups as thick bars connecting non-significant methods
  const groups = cdData.groups;
  const groupColors = [OI.blue, OI.green, OI.orange];
  groups.forEach((g, gi) => {{
    const methodNames = g.methods.map(m => {{
      // Map back from key to label
      const mapping = {json.dumps({v: k for k, v in METHOD_LABELS.items()})};
      return cdData.ranks[Object.keys(cdData.ranks).find(k => {{
        const key = mapping[k];
        return key !== undefined;
      }})];
    }});
    // Get rank range for this group
    const ranks = g.avg_ranks;
    const minR = Math.min(...ranks);
    const maxR = Math.max(...ranks);
    const barY = height - 30 + gi * 14;

    svg.append("line")
      .attr("x1", x(minR)).attr("x2", x(maxR))
      .attr("y1", barY).attr("y2", barY)
      .attr("stroke", groupColors[gi % groupColors.length])
      .attr("stroke-width", 3)
      .attr("opacity", 0.7);
    svg.append("text").attr("x", x((minR + maxR) / 2)).attr("y", barY - 4)
      .attr("text-anchor", "middle").attr("font-size", "8px")
      .attr("fill", groupColors[gi % groupColors.length])
      .text(`Group ${{g.group}}`);
  }});

  // Stat annotation
  svg.append("text").attr("x", width / 2).attr("y", height + 30)
    .attr("text-anchor", "middle").attr("font-size", "9px").attr("fill", "#444")
    .text(`Friedman χ² = ${{cdData.chi2}}, p = ${{cdData.p_value.toFixed(4)}} | ${{cdData.k}} methods, ${{cdData.n}} tissues | CD = ${{cdData.cd}} (α = 0.05)`);
}})();
</script>
</body>
</html>"""

# --- Write output ---
out_path = os.path.join(OUT_DIR, "Fig1_benchmark.html")
with open(out_path, "w") as f:
    f.write(html)
print(f"Written: {out_path}")
print(f"  Panel A: {len(heatmap_data)} cells (8×8)")
print(f"  Panel B: {len(forest_data)} tissue forest plots")
print(f"  Panel C: {len(scatter_data)} tissue scatter points")
print(f"  Panel D: CD={cd_data['cd']}, Friedman p={cd_data['p_value']:.4f}")

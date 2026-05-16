#!/usr/bin/env python3
"""
Fig S1 — Per-Tissue Detailed Results
8 panels (one per tissue), all methods bar chart with CI + significance

Output: v4/figures/html/FigS1_pertissue.html
"""

import json
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EVAL_DIR = os.path.join(BASE, "evaluation")
OUT_DIR = os.path.join(BASE, "figures", "html")
os.makedirs(OUT_DIR, exist_ok=True)

with open(os.path.join(EVAL_DIR, "M1_summary.json")) as f:
    m1 = json.load(f)

ALL_TISSUES = ["thymus", "colon", "lung", "skin", "kidney", "liver", "gastrocnemius", "eye"]
METHODS = ["pca_lr", "elasticnet_lr", "xgb", "rf", "knn", "mlp", "svm_rbf", "tabnet"]
METHOD_LABELS = {
    "pca_lr": "PCA-LR", "elasticnet_lr": "EN-LR", "rf": "RF",
    "xgb": "XGB", "svm_rbf": "SVM", "knn": "kNN", "mlp": "MLP", "tabnet": "TabNet"
}
FEATURES = ["gene", "pathway_hallmark", "pathway_kegg", "combined"]
FEAT_LABELS = {"gene": "Gene", "pathway_hallmark": "Hallmark", "pathway_kegg": "KEGG", "combined": "Combined"}
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastrocnemius", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin", "lung": "Lung", "colon": "Colon"
}

# Build per-tissue data
tissue_data = {}
for tissue in ALL_TISSUES:
    tissue_data[tissue] = []
    for feat in FEATURES:
        for method in METHODS:
            entry = m1[tissue][feat][method]
            tissue_data[tissue].append({
                "method": METHOD_LABELS[method],
                "feature": FEAT_LABELS[feat],
                "auroc": round(entry["auroc"], 3),
                "std": round(entry["std"], 3),
                "ci_lower": round(entry["ci_lower"], 3),
                "perm_p": round(entry["perm_p"], 4),
                "cv": entry["cv"]
            })

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig S1 — Per-Tissue Detailed Results</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ font-family: Arial, Helvetica, sans-serif; margin: 0; padding: 20px; background: #fff; }}
  .figure-container {{ width: 1200px; margin: 0 auto; }}
  .figure-title {{ font-size: 14px; font-weight: bold; text-align: center; margin-bottom: 4px; }}
  .figure-subtitle {{ font-size: 10px; color: #666; text-align: center; margin-bottom: 16px; }}
  .panels {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
  .panel {{ border: 1px solid #e0e0e0; border-radius: 4px; padding: 10px; background: #fafafa; }}
  .panel-label {{ font-size: 11px; font-weight: bold; margin-bottom: 4px; }}
  .axis text {{ font-size: 8px; }}
  .tooltip {{
    position: absolute; padding: 6px 10px; background: rgba(0,0,0,0.85);
    color: #fff; border-radius: 4px; font-size: 10px; pointer-events: none; z-index: 100; line-height: 1.4;
  }}
</style>
</head>
<body>
<div class="figure-container">
  <div class="figure-title">Fig S1 — Per-Tissue Classification Performance (All Methods × All Features)</div>
  <div class="figure-subtitle">8 methods × 4 feature types per tissue | * p < 0.05, ** p < 0.01 (permutation test)</div>
  <div class="panels" id="panels"></div>
</div>
<div class="tooltip" id="tooltip" style="display:none;"></div>

<script>
const OI = {{
  black: "#000000", orange: "#E69F00", skyblue: "#56B4E9", green: "#009E73",
  yellow: "#F0E442", blue: "#0072B2", vermilion: "#D55E00", pink: "#CC79A7"
}};
const featColors = {{"Gene": OI.blue, "Hallmark": OI.orange, "KEGG": OI.green, "Combined": OI.pink}};
const tooltip = d3.select("#tooltip");
const tissueData = {json.dumps(tissue_data)};
const tissueLabels = {json.dumps(TISSUE_LABELS)};
const tissueOrder = {json.dumps(ALL_TISSUES)};

tissueOrder.forEach((tissue, idx) => {{
  const data = tissueData[tissue];
  const label = tissueLabels[tissue];
  const cv = data[0].cv;

  const panel = d3.select("#panels").append("div").attr("class", "panel");
  panel.append("div").attr("class", "panel-label").text(`${{String.fromCharCode(65 + idx)}}. ${{label}} (${{cv}})`);

  const margin = {{top: 10, right: 10, bottom: 45, left: 40}};
  const methods = [...new Set(data.map(d => d.method))];
  const features = [...new Set(data.map(d => d.feature))];
  const barW = 10;
  const groupW = methods.length * (features.length * barW + 4);
  const width = methods.length * (features.length * barW + 12) + 20;
  const height = 180;

  const svg = panel.append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x0 = d3.scaleBand().domain(methods).range([0, width]).padding(0.15);
  const x1 = d3.scaleBand().domain(features).range([0, x0.bandwidth()]).padding(0.05);
  const y = d3.scaleLinear().domain([0, 1.05]).range([height, 0]);

  // 0.5 line
  svg.append("line").attr("x1", 0).attr("x2", width).attr("y1", y(0.5)).attr("y2", y(0.5))
    .attr("stroke", "#eee").attr("stroke-dasharray", "3,3");

  data.forEach(d => {{
    const bx = x0(d.method) + x1(d.feature);
    svg.append("rect")
      .attr("x", bx).attr("y", y(d.auroc))
      .attr("width", x1.bandwidth())
      .attr("height", Math.max(0, height - y(d.auroc)))
      .attr("fill", featColors[d.feature]).attr("opacity", 0.8)
      .on("mouseover", function(event) {{
        tooltip.style("display", "block")
          .html(`<b>${{d.method}}</b> (${{d.feature}})<br>AUROC: ${{d.auroc}} ± ${{d.std}}<br>CI lower: ${{d.ci_lower}}<br>p: ${{d.perm_p}}`);
      }})
      .on("mousemove", function(event) {{
        tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
      }})
      .on("mouseout", () => tooltip.style("display", "none"));

    // Significance
    if (d.perm_p < 0.05) {{
      svg.append("text")
        .attr("x", bx + x1.bandwidth()/2).attr("y", y(d.auroc) - 2)
        .attr("text-anchor", "middle").attr("font-size", "7px").attr("fill", OI.vermilion)
        .text(d.perm_p < 0.01 ? "**" : "*");
    }}
  }});

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x0).tickSize(0))
    .selectAll("text").style("font-size", "7px").attr("transform", "rotate(-30)").style("text-anchor", "end");
  svg.append("g").call(d3.axisLeft(y).ticks(5).tickFormat(d3.format(".1f")));

  // Feature legend (only first panel)
  if (idx === 0) {{
    let lx = width - 80, ly = 5;
    Object.entries(featColors).forEach(([feat, col]) => {{
      svg.append("rect").attr("x", lx).attr("y", ly).attr("width", 8).attr("height", 8).attr("fill", col).attr("opacity", 0.8);
      svg.append("text").attr("x", lx + 12).attr("y", ly + 8).attr("font-size", "7px").text(feat);
      ly += 12;
    }});
  }}
}});
</script>
</body>
</html>"""

out_path = os.path.join(OUT_DIR, "FigS1_pertissue.html")
with open(out_path, "w") as f:
    f.write(html)
print(f"Written: {out_path}")
print(f"  {len(ALL_TISSUES)} tissue panels × {len(METHODS)} methods × {len(FEATURES)} features")

#!/usr/bin/env python3
"""
Fig S5 — Negative Controls & Batch Detection
A: NC1 shuffled labels — AUROC violin/box per method
B: NC2 random features — AUROC violin/box per method
C: D3 mission ID prediction — macro-F1 heatmap (gene vs pathway)
D: D5 hardware type prediction — macro-F1 bar chart

Requires: NC1_shuffled_*.json, NC2_random_*.json, D3_v4_*.json, D5_v4_*.json
         in v4/evaluation/

Output: v4/figures/html/FigS5_controls.html
"""

import json
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EVAL_DIR = os.path.join(BASE, "evaluation")
OUT_DIR = os.path.join(BASE, "figures", "html")
os.makedirs(OUT_DIR, exist_ok=True)

METHODS = ["pca_lr", "elasticnet_lr", "xgb", "rf", "knn", "mlp", "svm_rbf", "tabnet"]
METHOD_LABELS = {
    "pca_lr": "PCA-LR", "elasticnet_lr": "EN-LR", "rf": "RF",
    "xgb": "XGB", "svm_rbf": "SVM", "knn": "kNN", "mlp": "MLP", "tabnet": "TabNet"
}
ALL_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin", "lung", "colon"]
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastro", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin", "lung": "Lung", "colon": "Colon"
}

# --- Load NC1 data ---
nc1_data = {}
for tissue in ALL_TISSUES:
    for method in METHODS:
        fname = f"NC1_shuffled_{tissue}_{method}.json"
        fpath = os.path.join(EVAL_DIR, fname)
        if os.path.exists(fpath):
            with open(fpath) as f:
                d = json.load(f)
            key = METHOD_LABELS.get(method, method)
            if key not in nc1_data:
                nc1_data[key] = []
            nc1_data[key].append({
                "tissue": TISSUE_LABELS.get(tissue, tissue),
                "auroc": round(d.get("mean_auroc", 0.5), 4),
                "folds": [round(f["auroc"], 4) for f in d.get("folds", [])]
            })

# --- Load NC2 data ---
nc2_data = {}
for tissue in ALL_TISSUES:
    for method in METHODS:
        fname = f"NC2_random_{tissue}_{method}.json"
        fpath = os.path.join(EVAL_DIR, fname)
        if os.path.exists(fpath):
            with open(fpath) as f:
                d = json.load(f)
            key = METHOD_LABELS.get(method, method)
            if key not in nc2_data:
                nc2_data[key] = []
            nc2_data[key].append({
                "tissue": TISSUE_LABELS.get(tissue, tissue),
                "auroc": round(d.get("mean_auroc", 0.5), 4),
                "folds": [round(f["auroc"], 4) for f in d.get("folds", [])]
            })

# --- Load D3 data ---
d3_data = {}
for tissue in ALL_TISSUES:
    for method in METHODS:
        for feat in ["gene", "pathway_hallmark"]:
            fname = f"D3_v4_{tissue}_{method}_{feat}.json"
            fpath = os.path.join(EVAL_DIR, fname)
            if os.path.exists(fpath):
                with open(fpath) as f:
                    d = json.load(f)
                tl = TISSUE_LABELS.get(tissue, tissue)
                if tl not in d3_data:
                    d3_data[tl] = {}
                ml = METHOD_LABELS.get(method, method)
                if ml not in d3_data[tl]:
                    d3_data[tl][ml] = {}
                d3_data[tl][ml][feat] = round(d.get("macro_f1", 0), 4)

# --- Load D5 data ---
d5_data = {}
for tissue in ALL_TISSUES:
    for method in METHODS:
        for feat in ["gene", "pathway_hallmark"]:
            fname = f"D5_v4_{tissue}_{method}_{feat}.json"
            fpath = os.path.join(EVAL_DIR, fname)
            if os.path.exists(fpath):
                with open(fpath) as f:
                    d = json.load(f)
                tl = TISSUE_LABELS.get(tissue, tissue)
                if tl not in d5_data:
                    d5_data[tl] = {}
                ml = METHOD_LABELS.get(method, method)
                if ml not in d5_data[tl]:
                    d5_data[tl][ml] = {}
                d5_data[tl][ml][feat] = round(d.get("macro_f1", 0), 4)

# Check data availability
n_nc1 = sum(len(v) for v in nc1_data.values())
n_nc2 = sum(len(v) for v in nc2_data.values())
n_d3 = sum(len(m) for m in d3_data.values())
n_d5 = sum(len(m) for m in d5_data.values())

if n_nc1 + n_nc2 + n_d3 + n_d5 == 0:
    print("ERROR: No NC1/NC2/D3/D5 JSON files found in v4/evaluation/")
    print("       Run negative control scripts on HPC first, then sync results.")
    sys.exit(1)

method_labels = [METHOD_LABELS[m] for m in METHODS]

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig S5 — Negative Controls & Batch Detection</title>
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
  .note {{ font-size: 8px; color: #888; margin-top: 8px; }}
</style>
</head>
<body>
<div class="figure-container">
  <div class="figure-title">Fig S5 — Negative Controls & Batch Detection</div>
  <div class="figure-subtitle">NC1: shuffled labels, NC2: random features, D3: mission ID, D5: hardware type</div>
  <div class="panels" id="panels"></div>
  <div class="note">
    A-B: All AUROC values should cluster around 0.5 (chance). Any method with AUROC &gt; 0.65 indicates data leakage.
    C-D: High macro-F1 for gene features confirms batch effects; low F1 for pathway features confirms batch-resistance.
  </div>
</div>
<div class="tooltip" id="tooltip" style="display:none;"></div>

<script>
const OI = {{
  black: "#000000", orange: "#E69F00", skyblue: "#56B4E9", green: "#009E73",
  yellow: "#F0E442", blue: "#0072B2", vermilion: "#D55E00", pink: "#CC79A7"
}};
const tooltip = d3.select("#tooltip");
const nc1Data = {json.dumps(nc1_data)};
const nc2Data = {json.dumps(nc2_data)};
const d3Data = {json.dumps(d3_data)};
const d5Data = {json.dumps(d5_data)};
const methodLabels = {json.dumps(method_labels)};

// Helper: draw box/strip plot for NC data
function drawNCPanel(parentId, panelLabel, data, color) {{
  const panel = d3.select("#panels").append("div").attr("class", "panel");
  panel.append("div").attr("class", "panel-label").text(panelLabel);

  const methods = Object.keys(data);
  if (methods.length === 0) {{
    panel.append("div").style("font-size", "9px").style("color", "#999")
      .text("No data available. Run NC scripts on HPC first.");
    return;
  }}

  const margin = {{top: 10, right: 10, bottom: 40, left: 40}};
  const width = 450, height = 200;

  const svg = panel.append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(methods).range([0, width]).padding(0.2);
  const y = d3.scaleLinear().domain([0, 1]).range([height, 0]);

  // 0.5 reference line
  svg.append("line").attr("x1", 0).attr("x2", width).attr("y1", y(0.5)).attr("y2", y(0.5))
    .attr("stroke", "#ccc").attr("stroke-dasharray", "4,3");

  methods.forEach(method => {{
    const vals = data[method].map(d => d.auroc);
    const cx = x(method) + x.bandwidth() / 2;

    // Jittered points
    vals.forEach((v, i) => {{
      const jitter = (Math.random() - 0.5) * x.bandwidth() * 0.5;
      svg.append("circle")
        .attr("cx", cx + jitter).attr("cy", y(v)).attr("r", 3)
        .attr("fill", color).attr("opacity", 0.6)
        .on("mouseover", function(event) {{
          const entry = data[method][i];
          tooltip.style("display", "block")
            .html(`<b>${{method}}</b> (${{entry.tissue}})<br>AUROC: ${{entry.auroc}}`);
        }})
        .on("mousemove", function(event) {{
          tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
        }})
        .on("mouseout", () => tooltip.style("display", "none"));
    }});

    // Mean line
    const mean = d3.mean(vals);
    svg.append("line")
      .attr("x1", x(method) + 2).attr("x2", x(method) + x.bandwidth() - 2)
      .attr("y1", y(mean)).attr("y2", y(mean))
      .attr("stroke", "#333").attr("stroke-width", 2);
  }});

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).tickSize(0))
    .selectAll("text").attr("transform", "rotate(-30)").style("text-anchor", "end").style("font-size", "8px");
  svg.append("g").call(d3.axisLeft(y).ticks(5).tickFormat(d3.format(".1f")));
  svg.append("text").attr("transform", "rotate(-90)").attr("x", -height/2).attr("y", -30)
    .attr("text-anchor", "middle").attr("font-size", "9px").text("AUROC");
}}

// Panel A: NC1
drawNCPanel("panels", "A. NC1 — Shuffled Labels (expected AUROC \\u2248 0.5)", nc1Data, OI.blue);

// Panel B: NC2
drawNCPanel("panels", "B. NC2 — Random Features (expected AUROC \\u2248 0.5)", nc2Data, OI.orange);

// Panel C: D3 mission ID
(() => {{
  const panel = d3.select("#panels").append("div").attr("class", "panel");
  panel.append("div").attr("class", "panel-label").text("C. D3 — Mission ID Prediction (macro-F1)");

  const tissues = Object.keys(d3Data);
  if (tissues.length === 0) {{
    panel.append("div").style("font-size", "9px").style("color", "#999")
      .text("No data available. Run D3 script on HPC first.");
    return;
  }}

  const features = ["gene", "pathway_hallmark"];
  const featLabels = {{"gene": "Gene", "pathway_hallmark": "Pathway"}};
  const margin = {{top: 10, right: 10, bottom: 45, left: 60}};
  const cellW = 30, cellH = 20;
  const methods = methodLabels.filter(m => tissues.some(t => d3Data[t] && d3Data[t][m]));

  features.forEach((feat, fidx) => {{
    panel.append("div").style("font-size", "9px").style("font-weight", "bold")
      .style("margin-top", fidx > 0 ? "8px" : "0").text(featLabels[feat]);

    const width = methods.length * cellW;
    const height = tissues.length * cellH;

    const svg = panel.append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

    const x = d3.scaleBand().domain(methods).range([0, width]).padding(0.05);
    const yScale = d3.scaleBand().domain(tissues).range([0, height]).padding(0.05);
    const color = d3.scaleSequential(d3.interpolateYlOrRd).domain([0, 1]);

    tissues.forEach(tissue => {{
      methods.forEach(method => {{
        const val = (d3Data[tissue] && d3Data[tissue][method]) ? d3Data[tissue][method][feat] : null;
        if (val === null || val === undefined) return;

        svg.append("rect")
          .attr("x", x(method)).attr("y", yScale(tissue))
          .attr("width", x.bandwidth()).attr("height", yScale.bandwidth())
          .attr("fill", color(val)).attr("rx", 2)
          .on("mouseover", function(event) {{
            tooltip.style("display", "block")
              .html(`<b>${{tissue}} / ${{method}}</b><br>F1=${{val}}`);
          }})
          .on("mousemove", function(event) {{
            tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
          }})
          .on("mouseout", () => tooltip.style("display", "none"));

        svg.append("text")
          .attr("x", x(method) + x.bandwidth()/2).attr("y", yScale(tissue) + yScale.bandwidth()/2 + 3)
          .attr("text-anchor", "middle").attr("font-size", "6px")
          .attr("fill", val > 0.5 ? "#fff" : "#333")
          .text(val.toFixed(2));
      }});
    }});

    svg.append("g").attr("transform", `translate(0,${{height}})`)
      .call(d3.axisBottom(x).tickSize(0))
      .selectAll("text").attr("transform", "rotate(-30)").style("text-anchor", "end").style("font-size", "7px");
    svg.append("g").call(d3.axisLeft(yScale).tickSize(0))
      .selectAll("text").style("font-size", "7px");
  }});
}})();

// Panel D: D5 hardware
(() => {{
  const panel = d3.select("#panels").append("div").attr("class", "panel");
  panel.append("div").attr("class", "panel-label").text("D. D5 — Hardware Type Prediction (RR vs MHU)");

  const tissues = Object.keys(d5Data);
  if (tissues.length === 0) {{
    panel.append("div").style("font-size", "9px").style("color", "#999")
      .text("No data available. Run D5 script on HPC first.");
    return;
  }}

  const margin = {{top: 10, right: 10, bottom: 40, left: 60}};
  const methods = methodLabels.filter(m => tissues.some(t => d5Data[t] && d5Data[t][m]));
  const barW = 12;
  const width = methods.length * (barW * 2 + 8) * tissues.length + 40;
  const height = 180;

  const svg = panel.append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const xGroup = d3.scaleBand().domain(methods).range([0, width]).padding(0.2);
  const xFeat = d3.scaleBand().domain(["gene", "pathway_hallmark"]).range([0, xGroup.bandwidth()]).padding(0.1);
  const y = d3.scaleLinear().domain([0, 1.05]).range([height, 0]);
  const featColors = {{"gene": OI.blue, "pathway_hallmark": OI.orange}};

  tissues.forEach((tissue, tidx) => {{
    methods.forEach(method => {{
      ["gene", "pathway_hallmark"].forEach(feat => {{
        const val = (d5Data[tissue] && d5Data[tissue][method]) ? d5Data[tissue][method][feat] : null;
        if (val === null || val === undefined) return;

        const bx = xGroup(method) + xFeat(feat);
        svg.append("rect")
          .attr("x", bx).attr("y", y(val))
          .attr("width", xFeat.bandwidth()).attr("height", Math.max(0, height - y(val)))
          .attr("fill", featColors[feat]).attr("opacity", 0.8).attr("rx", 1)
          .on("mouseover", function(event) {{
            tooltip.style("display", "block")
              .html(`<b>${{tissue}} / ${{method}}</b><br>${{feat}}<br>F1=${{val}}`);
          }})
          .on("mousemove", function(event) {{
            tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
          }})
          .on("mouseout", () => tooltip.style("display", "none"));
      }});
    }});
  }});

  // 0.5 line
  svg.append("line").attr("x1", 0).attr("x2", width).attr("y1", y(0.5)).attr("y2", y(0.5))
    .attr("stroke", "#ccc").attr("stroke-dasharray", "3,3");

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(xGroup).tickSize(0))
    .selectAll("text").attr("transform", "rotate(-30)").style("text-anchor", "end").style("font-size", "7px");
  svg.append("g").call(d3.axisLeft(y).ticks(5).tickFormat(d3.format(".1f")));

  // Legend
  let lx = width - 80, ly = 5;
  Object.entries(featColors).forEach(([feat, col]) => {{
    svg.append("rect").attr("x", lx).attr("y", ly).attr("width", 8).attr("height", 8).attr("fill", col).attr("opacity", 0.8);
    svg.append("text").attr("x", lx + 12).attr("y", ly + 8).attr("font-size", "7px")
      .text(feat === "gene" ? "Gene" : "Pathway");
    ly += 12;
  }});
}})();
</script>
</body>
</html>"""

out_path = os.path.join(OUT_DIR, "FigS5_controls.html")
with open(out_path, "w") as f:
    f.write(html)
print(f"Written: {out_path}")
print(f"  NC1: {n_nc1} entries, NC2: {n_nc2} entries, D3: {n_d3} entries, D5: {n_d5} entries")

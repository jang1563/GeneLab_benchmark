#!/usr/bin/env python3
"""
Fig S2 — DeLong Pairwise AUROC Comparison Heatmaps
8 tissues × C(8,2)=28 method pairs, BH-FDR corrected p-values

Output: v4/figures/html/FigS2_delong.html
"""

import json
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EVAL_DIR = os.path.join(BASE, "evaluation")
OUT_DIR = os.path.join(BASE, "figures", "html")
os.makedirs(OUT_DIR, exist_ok=True)

with open(os.path.join(EVAL_DIR, "META_method_ranking.json")) as f:
    meta = json.load(f)

METHODS = ["pca_lr", "elasticnet_lr", "xgb", "rf", "knn", "mlp", "svm_rbf", "tabnet"]
METHOD_LABELS = {
    "pca_lr": "PCA-LR", "elasticnet_lr": "EN-LR", "rf": "RF",
    "xgb": "XGB", "svm_rbf": "SVM", "knn": "kNN", "mlp": "MLP", "tabnet": "TabNet"
}
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastro", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin", "lung": "Lung", "colon": "Colon"
}

# Build DeLong data
delong_data = {}
delong_raw = meta.get("delong_pairwise", {})
for tissue, pairs in delong_raw.items():
    tissue_label = TISSUE_LABELS.get(tissue, tissue)
    delong_data[tissue_label] = []
    for pair_key, val in pairs.items():
        methods = pair_key.split("_vs_")
        if len(methods) == 2:
            ma, mb = methods
        else:
            # Handle complex method names (e.g., elasticnet_lr_vs_svm_rbf)
            parts = pair_key.split("_vs_")
            ma = parts[0]
            mb = parts[1] if len(parts) > 1 else ""
        delong_data[tissue_label].append({
            "method_a": METHOD_LABELS.get(ma, ma),
            "method_b": METHOD_LABELS.get(mb, mb),
            "z": round(val.get("z", 0), 3),
            "p_raw": round(val.get("p_raw", 1), 4),
            "p_fdr": round(val.get("p_fdr", 1), 4),
            "auc_diff": round(val.get("auc_diff", 0), 4)
        })

# Also get Wilcoxon data
wilcoxon_data = {}
wilcoxon_raw = meta.get("wilcoxon_pairwise", {})
for tissue, pairs in wilcoxon_raw.items():
    tissue_label = TISSUE_LABELS.get(tissue, tissue)
    wilcoxon_data[tissue_label] = []
    for pair_key, val in pairs.items():
        methods = pair_key.split("_vs_")
        if len(methods) >= 2:
            ma = methods[0]
            mb = methods[1]
        else:
            continue
        wilcoxon_data[tissue_label].append({
            "method_a": METHOD_LABELS.get(ma, ma),
            "method_b": METHOD_LABELS.get(mb, mb),
            "stat": val.get("stat"),
            "p_raw": round(val.get("p_raw", 1), 4),
            "p_fdr": round(val.get("p_fdr", 1), 4)
        })

method_labels = [METHOD_LABELS[m] for m in METHODS]

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig S2 — DeLong Pairwise AUROC Comparison</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ font-family: Arial, Helvetica, sans-serif; margin: 0; padding: 20px; background: #fff; }}
  .figure-container {{ width: 1200px; margin: 0 auto; }}
  .figure-title {{ font-size: 14px; font-weight: bold; text-align: center; margin-bottom: 4px; }}
  .figure-subtitle {{ font-size: 10px; color: #666; text-align: center; margin-bottom: 16px; }}
  .panels {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
  .panel {{ border: 1px solid #e0e0e0; border-radius: 4px; padding: 10px; background: #fafafa; }}
  .panel-label {{ font-size: 10px; font-weight: bold; margin-bottom: 4px; }}
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
  <div class="figure-title">Fig S2 — DeLong Pairwise AUROC Comparison (Gene Features)</div>
  <div class="figure-subtitle">8 tissues × C(8,2)=28 method pairs | BH-FDR corrected p-values</div>
  <div class="panels" id="panels"></div>
  <div class="note">
    Bold border = p_FDR &lt; 0.05. Color = AUROC difference (row − column): blue = row better, red = column better.
    DeLong test on concatenated LOMO/5-fold predictions.
  </div>
</div>
<div class="tooltip" id="tooltip" style="display:none;"></div>

<script>
const OI = {{orange: "#E69F00", blue: "#0072B2", vermilion: "#D55E00"}};
const tooltip = d3.select("#tooltip");
const delongData = {json.dumps(delong_data)};
const methodLabels = {json.dumps(method_labels)};
const tissues = Object.keys(delongData);

tissues.forEach((tissue, idx) => {{
  const data = delongData[tissue];
  const panel = d3.select("#panels").append("div").attr("class", "panel");
  panel.append("div").attr("class", "panel-label").text(`${{tissue}}`);

  const margin = {{top: 5, right: 5, bottom: 50, left: 55}};
  const cellS = 38;
  const width = methodLabels.length * cellS;
  const height = methodLabels.length * cellS;

  const svg = panel.append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(methodLabels).range([0, width]).padding(0.05);
  const y = d3.scaleBand().domain(methodLabels).range([0, height]).padding(0.05);

  const color = d3.scaleDiverging()
    .domain([-0.3, 0, 0.3])
    .interpolator(d3.interpolateRdBu);

  // Build lookup
  const lookup = {{}};
  data.forEach(d => {{
    lookup[d.method_a + "_" + d.method_b] = d;
    // Reverse with negated diff
    lookup[d.method_b + "_" + d.method_a] = {{
      ...d,
      method_a: d.method_b, method_b: d.method_a,
      auc_diff: -d.auc_diff, z: -d.z
    }};
  }});

  methodLabels.forEach(ma => {{
    methodLabels.forEach(mb => {{
      if (ma === mb) {{
        svg.append("rect")
          .attr("x", x(mb)).attr("y", y(ma))
          .attr("width", x.bandwidth()).attr("height", y.bandwidth())
          .attr("fill", "#f0f0f0").attr("rx", 2);
        return;
      }}
      const entry = lookup[ma + "_" + mb];
      if (!entry) return;

      svg.append("rect")
        .attr("x", x(mb)).attr("y", y(ma))
        .attr("width", x.bandwidth()).attr("height", y.bandwidth())
        .attr("fill", color(entry.auc_diff))
        .attr("stroke", entry.p_fdr < 0.05 ? "#333" : "#eee")
        .attr("stroke-width", entry.p_fdr < 0.05 ? 2 : 0.5)
        .attr("rx", 2)
        .on("mouseover", function(event) {{
          tooltip.style("display", "block")
            .html(`<b>${{entry.method_a}} vs ${{entry.method_b}}</b><br>ΔAUROC: ${{entry.auc_diff > 0 ? "+" : ""}}${{entry.auc_diff}}<br>z: ${{entry.z}}<br>p_raw: ${{entry.p_raw}}<br>p_FDR: ${{entry.p_fdr}}`);
        }})
        .on("mousemove", function(event) {{
          tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
        }})
        .on("mouseout", () => tooltip.style("display", "none"));

      // p-value text
      if (entry.p_fdr < 0.05) {{
        svg.append("text")
          .attr("x", x(mb) + x.bandwidth()/2)
          .attr("y", y(ma) + y.bandwidth()/2 + 4)
          .attr("text-anchor", "middle").attr("font-size", "7px")
          .attr("font-weight", "bold")
          .attr("fill", Math.abs(entry.auc_diff) > 0.15 ? "#fff" : "#333")
          .text(entry.p_fdr < 0.001 ? "***" : entry.p_fdr < 0.01 ? "**" : "*");
      }}
    }});
  }});

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).tickSize(0))
    .selectAll("text").attr("transform", "rotate(-40)")
    .style("text-anchor", "end").style("font-size", "7px");
  svg.append("g").call(d3.axisLeft(y).tickSize(0))
    .selectAll("text").style("font-size", "7px");
}});
</script>
</body>
</html>"""

out_path = os.path.join(OUT_DIR, "FigS2_delong.html")
with open(out_path, "w") as f:
    f.write(html)
print(f"Written: {out_path}")
print(f"  {len(delong_data)} tissue panels with pairwise DeLong comparisons")

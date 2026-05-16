#!/usr/bin/env python3
"""
Fig 2 — Ablation Studies and Performance Sensitivity
GeneLabBench v4 Publication Figure

Panels:
  A: Feature count ablation (top-K genes) — line chart per tissue, 2 methods
  B: PCA component ablation — line chart per tissue
  C: Sample size learning curves — 4 largest tissues, 2 methods
  D: Bootstrap CI convergence — CI width vs n_bootstrap

Output: v4/figures/html/Fig2_ablation.html (self-contained D3.js v7)
"""

import json
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EVAL_DIR = os.path.join(BASE, "evaluation")
OUT_DIR = os.path.join(BASE, "figures", "html")
os.makedirs(OUT_DIR, exist_ok=True)

with open(os.path.join(EVAL_DIR, "ABL_summary.json")) as f:
    abl = json.load(f)

# --- Constants ---
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastrocnemius", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin", "lung": "Lung", "colon": "Colon"
}

# --- Panel A: Feature count ---
feat_data = {}
s1 = abl["S1_feature_count"]
K_VALUES = ["100", "500", "1000", "2000", "5000", "10000", "all"]
for tissue in s1:
    feat_data[tissue] = {}
    for method in s1[tissue]:
        feat_data[tissue][method] = []
        for k in K_VALUES:
            if k in s1[tissue][method]:
                entry = s1[tissue][method][k]
                feat_data[tissue][method].append({
                    "k": entry["top_k_actual"] if k != "all" else entry["top_k_actual"],
                    "k_label": k,
                    "auroc": round(entry["mean_auroc"], 4),
                    "std": round(entry["std_auroc"], 4)
                })

# --- Panel B: PCA components ---
pca_data = {}
s2 = abl["S2_pca_components"]
NC_VALUES = ["5", "10", "20", "50", "100", "200"]
for tissue in s2:
    pca_data[tissue] = []
    for nc in NC_VALUES:
        if nc in s2[tissue]:
            entry = s2[tissue][nc]
            pca_data[tissue].append({
                "nc": int(nc),
                "auroc": round(entry["mean_auroc"], 4),
                "std": round(entry["std_auroc"], 4),
                "var_explained": round(entry.get("mean_explained_var", 0), 4)
            })

# --- Panel C: Sample size ---
sample_data = {}
s3 = abl["S3_sample_size"]
FRACS = ["0.2", "0.4", "0.6", "0.8", "1.0"]
for tissue in s3:
    sample_data[tissue] = {}
    for method in s3[tissue]:
        sample_data[tissue][method] = []
        for frac in FRACS:
            if frac in s3[tissue][method]:
                entry = s3[tissue][method][frac]
                sample_data[tissue][method].append({
                    "frac": float(frac),
                    "auroc": round(entry["mean_auroc"], 4),
                    "std": round(entry["std_auroc"], 4)
                })

# --- Panel D: Bootstrap stability ---
boot_data = {}
s4 = abl["S4_bootstrap_stability"]
BOOT_VALUES = ["100", "500", "1000", "2000", "5000"]
for tissue in s4:
    boot_data[tissue] = []
    for nb in BOOT_VALUES:
        if nb in s4[tissue]:
            entry = s4[tissue][nb]
            boot_data[tissue].append({
                "n_boot": int(nb),
                "auroc": round(entry["mean_auroc"], 4),
                "ci_width": round(entry["mean_ci_width"], 4)
            })

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig 2 — Ablation Studies and Performance Sensitivity</title>
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
  .axis text {{ font-size: 9px; }}
  .tooltip {{
    position: absolute; padding: 6px 10px; background: rgba(0,0,0,0.85);
    color: #fff; border-radius: 4px; font-size: 10px; pointer-events: none; z-index: 100;
  }}
  .interpretation {{
    margin-top: 16px; padding: 10px; background: #f5f5f5;
    border-left: 3px solid #999; font-size: 9px; color: #444; line-height: 1.5;
  }}
  .interpretation strong {{ color: #333; }}
  .legend {{ font-size: 9px; }}
</style>
</head>
<body>
<div class="figure-container">
  <div class="figure-title">Fig 2 — Ablation Studies and Performance Sensitivity</div>
  <div class="figure-subtitle">Feature count, PCA dimensionality, sample size, and bootstrap convergence analysis</div>
  <div class="panels">
    <div class="panel"><div class="panel-label">A</div>
      <div class="panel-desc">Feature Count Ablation — Top-K Genes by Variance (PCA-LR & ElasticNet-LR)</div>
      <div id="feat_chart"></div>
    </div>
    <div class="panel"><div class="panel-label">B</div>
      <div class="panel-desc">PCA Component Ablation — PCA-LR (8 tissues)</div>
      <div id="pca_chart"></div>
    </div>
    <div class="panel"><div class="panel-label">C</div>
      <div class="panel-desc">Sample Size Learning Curves — Subsampled Training Data (4 tissues)</div>
      <div id="sample_chart"></div>
    </div>
    <div class="panel"><div class="panel-label">D</div>
      <div class="panel-desc">Bootstrap CI Width Convergence — PCA-LR Gene (8 tissues)</div>
      <div id="boot_chart"></div>
    </div>
  </div>
  <div class="interpretation">
    <strong>Key findings:</strong>
    (A) Performance peaks at K=1000–2000 genes for most tissues, then plateaus or slightly declines with all genes (~20K).
    (B) PCA with 20–50 components is optimal; more components do not consistently improve performance.
    (C) Learning curves show liver is most robust to subsampling; smaller tissues (gastrocnemius) are more sensitive.
    (D) Bootstrap CI width stabilizes by n=1000, justifying the choice of n_bootstrap=2000 for all evaluations.
  </div>
</div>
<div class="tooltip" id="tooltip" style="display:none;"></div>

<script>
const OI = {{
  black: "#000000", orange: "#E69F00", skyblue: "#56B4E9", green: "#009E73",
  yellow: "#F0E442", blue: "#0072B2", vermilion: "#D55E00", pink: "#CC79A7"
}};
const tissueColors = {{
  "liver": OI.blue, "gastrocnemius": OI.vermilion, "kidney": OI.green,
  "thymus": OI.orange, "eye": OI.pink, "skin": OI.skyblue,
  "lung": "#888888", "colon": "#555555"
}};
const tissueLabels = {json.dumps(TISSUE_LABELS)};
const tooltip = d3.select("#tooltip");

// =============================================
// PANEL A: Feature Count Ablation
// =============================================
(function() {{
  const data = {json.dumps(feat_data)};
  const margin = {{top: 15, right: 140, bottom: 45, left: 55}};
  const width = 370, height = 280;

  const svg = d3.select("#feat_chart").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const kLabels = ["100", "500", "1K", "2K", "5K", "10K", "All"];
  const kNums = [100, 500, 1000, 2000, 5000, 10000, 25000];
  const x = d3.scaleLog().domain([80, 30000]).range([0, width]);
  const y = d3.scaleLinear().domain([0.3, 1.0]).range([height, 0]);

  // Grid
  svg.append("line").attr("x1", 0).attr("x2", width).attr("y1", y(0.5)).attr("y2", y(0.5))
    .attr("stroke", "#eee").attr("stroke-dasharray", "3,3");

  const tissues = Object.keys(data);
  const methodDash = {{"pca_lr": "", "elasticnet_lr": "5,3"}};

  tissues.forEach(tissue => {{
    Object.keys(data[tissue]).forEach(method => {{
      const pts = data[tissue][method];
      if (pts.length === 0) return;

      const line = d3.line()
        .x(d => x(d.k))
        .y(d => y(d.auroc));

      svg.append("path")
        .datum(pts)
        .attr("d", line)
        .attr("fill", "none")
        .attr("stroke", tissueColors[tissue])
        .attr("stroke-width", method === "pca_lr" ? 2 : 1.5)
        .attr("stroke-dasharray", methodDash[method] || "")
        .attr("opacity", 0.8);

      // Dots
      svg.selectAll(`.dot-${{tissue}}-${{method}}`)
        .data(pts).enter().append("circle")
        .attr("cx", d => x(d.k)).attr("cy", d => y(d.auroc))
        .attr("r", 3)
        .attr("fill", tissueColors[tissue])
        .attr("opacity", 0.8)
        .on("mouseover", function(event, d) {{
          tooltip.style("display", "block")
            .html(`<b>${{tissueLabels[tissue]}}</b> (${{method === "pca_lr" ? "PCA-LR" : "EN-LR"}})<br>K=${{d.k_label}} (${{d.k}})<br>AUROC: ${{d.auroc.toFixed(3)}} ± ${{d.std.toFixed(3)}}`);
        }})
        .on("mousemove", function(event) {{
          tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
        }})
        .on("mouseout", () => tooltip.style("display", "none"));
    }});
  }});

  // X axis
  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).tickValues(kNums).tickFormat(d => {{
      if (d >= 1000) return d >= 25000 ? "All" : (d/1000)+"K";
      return d;
    }}));
  svg.append("text").attr("x", width/2).attr("y", height+35)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("Top-K Genes by Variance");

  // Y axis
  svg.append("g").call(d3.axisLeft(y).ticks(7).tickFormat(d3.format(".1f")));
  svg.append("text").attr("x", -height/2).attr("y", -40)
    .attr("text-anchor", "middle").attr("font-size", "10px")
    .attr("transform", "rotate(-90)").text("Mean AUROC");

  // Legend
  let ly = 0;
  tissues.forEach(t => {{
    svg.append("line").attr("x1", width+10).attr("x2", width+30)
      .attr("y1", ly+4).attr("y2", ly+4)
      .attr("stroke", tissueColors[t]).attr("stroke-width", 2);
    svg.append("text").attr("x", width+33).attr("y", ly+8)
      .attr("font-size", "8px").text(tissueLabels[t]);
    ly += 14;
  }});
  ly += 8;
  svg.append("line").attr("x1", width+10).attr("x2", width+30)
    .attr("y1", ly+4).attr("y2", ly+4)
    .attr("stroke", "#333").attr("stroke-width", 2);
  svg.append("text").attr("x", width+33).attr("y", ly+8)
    .attr("font-size", "8px").text("PCA-LR (solid)");
  ly += 14;
  svg.append("line").attr("x1", width+10).attr("x2", width+30)
    .attr("y1", ly+4).attr("y2", ly+4)
    .attr("stroke", "#333").attr("stroke-width", 1.5).attr("stroke-dasharray", "5,3");
  svg.append("text").attr("x", width+33).attr("y", ly+8)
    .attr("font-size", "8px").text("EN-LR (dashed)");
}})();

// =============================================
// PANEL B: PCA Components Ablation
// =============================================
(function() {{
  const data = {json.dumps(pca_data)};
  const margin = {{top: 15, right: 140, bottom: 45, left: 55}};
  const width = 370, height = 280;

  const svg = d3.select("#pca_chart").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleLog().domain([4, 250]).range([0, width]);
  const y = d3.scaleLinear().domain([0.3, 1.0]).range([height, 0]);

  svg.append("line").attr("x1", 0).attr("x2", width).attr("y1", y(0.5)).attr("y2", y(0.5))
    .attr("stroke", "#eee").attr("stroke-dasharray", "3,3");

  // Highlight n=50 (Phase 1 default)
  svg.append("line").attr("x1", x(50)).attr("x2", x(50))
    .attr("y1", 0).attr("y2", height)
    .attr("stroke", "#ddd").attr("stroke-dasharray", "3,3");
  svg.append("text").attr("x", x(50)+3).attr("y", 10)
    .attr("font-size", "7px").attr("fill", "#999").text("default (n=50)");

  const tissues = Object.keys(data);
  tissues.forEach(tissue => {{
    const pts = data[tissue];
    if (pts.length === 0) return;

    const line = d3.line().x(d => x(d.nc)).y(d => y(d.auroc));
    svg.append("path").datum(pts).attr("d", line)
      .attr("fill", "none").attr("stroke", tissueColors[tissue])
      .attr("stroke-width", 2).attr("opacity", 0.8);

    svg.selectAll(`.pca-dot-${{tissue}}`)
      .data(pts).enter().append("circle")
      .attr("cx", d => x(d.nc)).attr("cy", d => y(d.auroc))
      .attr("r", 3).attr("fill", tissueColors[tissue]).attr("opacity", 0.8)
      .on("mouseover", function(event, d) {{
        tooltip.style("display", "block")
          .html(`<b>${{tissueLabels[tissue]}}</b><br>n_components=${{d.nc}}<br>AUROC: ${{d.auroc.toFixed(3)}} ± ${{d.std.toFixed(3)}}<br>Var explained: ${{(d.var_explained*100).toFixed(1)}}%`);
      }})
      .on("mousemove", function(event) {{
        tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
      }})
      .on("mouseout", () => tooltip.style("display", "none"));
  }});

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).tickValues([5,10,20,50,100,200]).tickFormat(d => d));
  svg.append("text").attr("x", width/2).attr("y", height+35)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("Number of PCA Components");

  svg.append("g").call(d3.axisLeft(y).ticks(7).tickFormat(d3.format(".1f")));
  svg.append("text").attr("x", -height/2).attr("y", -40)
    .attr("text-anchor", "middle").attr("font-size", "10px")
    .attr("transform", "rotate(-90)").text("Mean AUROC (PCA-LR)");

  // Legend
  let ly = 0;
  tissues.forEach(t => {{
    svg.append("line").attr("x1", width+10).attr("x2", width+30)
      .attr("y1", ly+4).attr("y2", ly+4)
      .attr("stroke", tissueColors[t]).attr("stroke-width", 2);
    svg.append("text").attr("x", width+33).attr("y", ly+8)
      .attr("font-size", "8px").text(tissueLabels[t]);
    ly += 14;
  }});
}})();

// =============================================
// PANEL C: Sample Size Learning Curves
// =============================================
(function() {{
  const data = {json.dumps(sample_data)};
  const margin = {{top: 15, right: 140, bottom: 45, left: 55}};
  const width = 370, height = 280;

  const svg = d3.select("#sample_chart").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleLinear().domain([0.15, 1.05]).range([0, width]);
  const y = d3.scaleLinear().domain([0.3, 1.0]).range([height, 0]);

  svg.append("line").attr("x1", 0).attr("x2", width).attr("y1", y(0.5)).attr("y2", y(0.5))
    .attr("stroke", "#eee").attr("stroke-dasharray", "3,3");

  const methodDash = {{"pca_lr": "", "elasticnet_lr": "5,3"}};
  const tissues = Object.keys(data);

  tissues.forEach(tissue => {{
    Object.keys(data[tissue]).forEach(method => {{
      const pts = data[tissue][method];
      if (pts.length === 0) return;

      const line = d3.line().x(d => x(d.frac)).y(d => y(d.auroc));
      svg.append("path").datum(pts).attr("d", line)
        .attr("fill", "none").attr("stroke", tissueColors[tissue])
        .attr("stroke-width", method === "pca_lr" ? 2 : 1.5)
        .attr("stroke-dasharray", methodDash[method] || "")
        .attr("opacity", 0.8);

      // Error bars + dots
      pts.forEach(d => {{
        if (d.std > 0) {{
          svg.append("line")
            .attr("x1", x(d.frac)).attr("x2", x(d.frac))
            .attr("y1", y(Math.max(0, d.auroc - d.std)))
            .attr("y2", y(Math.min(1, d.auroc + d.std)))
            .attr("stroke", tissueColors[tissue]).attr("stroke-width", 1).attr("opacity", 0.4);
        }}
        svg.append("circle")
          .attr("cx", x(d.frac)).attr("cy", y(d.auroc))
          .attr("r", 3).attr("fill", tissueColors[tissue]).attr("opacity", 0.8);
      }});
    }});
  }});

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).tickValues([0.2,0.4,0.6,0.8,1.0]).tickFormat(d => (d*100)+"%"));
  svg.append("text").attr("x", width/2).attr("y", height+35)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("Training Data Fraction");

  svg.append("g").call(d3.axisLeft(y).ticks(7).tickFormat(d3.format(".1f")));
  svg.append("text").attr("x", -height/2).attr("y", -40)
    .attr("text-anchor", "middle").attr("font-size", "10px")
    .attr("transform", "rotate(-90)").text("Mean AUROC");

  // Legend
  let ly = 0;
  tissues.forEach(t => {{
    svg.append("line").attr("x1", width+10).attr("x2", width+30)
      .attr("y1", ly+4).attr("y2", ly+4)
      .attr("stroke", tissueColors[t]).attr("stroke-width", 2);
    svg.append("text").attr("x", width+33).attr("y", ly+8)
      .attr("font-size", "8px").text(tissueLabels[t]);
    ly += 14;
  }});
  ly += 8;
  svg.append("line").attr("x1", width+10).attr("x2", width+30)
    .attr("y1", ly+4).attr("y2", ly+4)
    .attr("stroke", "#333").attr("stroke-width", 2);
  svg.append("text").attr("x", width+33).attr("y", ly+8)
    .attr("font-size", "8px").text("PCA-LR (solid)");
  ly += 14;
  svg.append("line").attr("x1", width+10).attr("x2", width+30)
    .attr("y1", ly+4).attr("y2", ly+4)
    .attr("stroke", "#333").attr("stroke-width", 1.5).attr("stroke-dasharray", "5,3");
  svg.append("text").attr("x", width+33).attr("y", ly+8)
    .attr("font-size", "8px").text("EN-LR (dashed)");
}})();

// =============================================
// PANEL D: Bootstrap CI Convergence
// =============================================
(function() {{
  const data = {json.dumps(boot_data)};
  const margin = {{top: 15, right: 140, bottom: 45, left: 55}};
  const width = 370, height = 280;

  const svg = d3.select("#boot_chart").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleLog().domain([80, 6000]).range([0, width]);
  const y = d3.scaleLinear().domain([0, 0.7]).range([height, 0]);

  // Highlight n=2000 (chosen default)
  svg.append("line").attr("x1", x(2000)).attr("x2", x(2000))
    .attr("y1", 0).attr("y2", height)
    .attr("stroke", "#ddd").attr("stroke-dasharray", "3,3");
  svg.append("text").attr("x", x(2000)+3).attr("y", 10)
    .attr("font-size", "7px").attr("fill", "#999").text("chosen (n=2000)");

  const tissues = Object.keys(data);
  tissues.forEach(tissue => {{
    const pts = data[tissue];
    if (pts.length === 0) return;

    const line = d3.line().x(d => x(d.n_boot)).y(d => y(d.ci_width));
    svg.append("path").datum(pts).attr("d", line)
      .attr("fill", "none").attr("stroke", tissueColors[tissue])
      .attr("stroke-width", 2).attr("opacity", 0.8);

    svg.selectAll(`.boot-dot-${{tissue}}`)
      .data(pts).enter().append("circle")
      .attr("cx", d => x(d.n_boot)).attr("cy", d => y(d.ci_width))
      .attr("r", 3).attr("fill", tissueColors[tissue]).attr("opacity", 0.8)
      .on("mouseover", function(event, d) {{
        tooltip.style("display", "block")
          .html(`<b>${{tissueLabels[tissue]}}</b><br>n_bootstrap=${{d.n_boot}}<br>AUROC: ${{d.auroc.toFixed(4)}}<br>CI width: ${{d.ci_width.toFixed(4)}}`);
      }})
      .on("mousemove", function(event) {{
        tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
      }})
      .on("mouseout", () => tooltip.style("display", "none"));
  }});

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).tickValues([100,500,1000,2000,5000]).tickFormat(d => d >= 1000 ? (d/1000)+"K" : d));
  svg.append("text").attr("x", width/2).attr("y", height+35)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("Number of Bootstrap Resamples");

  svg.append("g").call(d3.axisLeft(y).ticks(7).tickFormat(d3.format(".2f")));
  svg.append("text").attr("x", -height/2).attr("y", -40)
    .attr("text-anchor", "middle").attr("font-size", "10px")
    .attr("transform", "rotate(-90)").text("Mean 95% CI Width");

  // Legend
  let ly = 0;
  tissues.forEach(t => {{
    svg.append("line").attr("x1", width+10).attr("x2", width+30)
      .attr("y1", ly+4).attr("y2", ly+4)
      .attr("stroke", tissueColors[t]).attr("stroke-width", 2);
    svg.append("text").attr("x", width+33).attr("y", ly+8)
      .attr("font-size", "8px").text(tissueLabels[t]);
    ly += 14;
  }});
}})();
</script>
</body>
</html>"""

out_path = os.path.join(OUT_DIR, "Fig2_ablation.html")
with open(out_path, "w") as f:
    f.write(html)
print(f"Written: {out_path}")
print(f"  Panel A: {len(feat_data)} tissues × 2 methods × 7 K-levels")
print(f"  Panel B: {len(pca_data)} tissues × 6 PCA levels")
print(f"  Panel C: {len(sample_data)} tissues × 2 methods × 5 fractions")
print(f"  Panel D: {len(boot_data)} tissues × 5 bootstrap levels")

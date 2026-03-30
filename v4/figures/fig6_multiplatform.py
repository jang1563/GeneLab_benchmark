#!/usr/bin/env python3
"""
Fig 6 — Multi-Platform Extension (scRNA, Spatial, Power Analysis)
GeneLabBench v4 Publication Figure

Panels:
  A: scRNA LOAO AUROC (RRRM-2, top cell types per tissue)
  B: Spatial Visium brain summary (negative result)
  C: Power analysis recommendations (sample size curves)
  D: Summary schematic table (GeneLabBench overview)

Output: v4/figures/html/Fig6_multiplatform.html
"""

import json
import os
import math

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
V3_EVAL = os.path.join(os.path.dirname(BASE), "v3", "evaluation")
V4_EVAL = os.path.join(BASE, "evaluation")
OUT_DIR = os.path.join(BASE, "figures", "html")
os.makedirs(OUT_DIR, exist_ok=True)

with open(os.path.join(V3_EVAL, "F5C_rrrm2_loao.json")) as f:
    f5c = json.load(f)
with open(os.path.join(V4_EVAL, "META_power_analysis.json")) as f:
    power = json.load(f)

# --- Panel A: scRNA LOAO ---
scrna_data = []
for tissue_key, cell_types in f5c["results"].items():
    tissue_label = {"pbmc": "PBMC", "spleen": "Spleen", "femur_bm": "Bone Marrow (Femur)",
                    "humerus_bm": "Bone Marrow (Humerus)"}.get(tissue_key, tissue_key)
    for ct, metrics in cell_types.items():
        auroc = metrics.get("auroc")
        if auroc is None or (isinstance(auroc, float) and math.isnan(auroc)):
            continue
        scrna_data.append({
            "tissue": tissue_label,
            "cell_type": ct.replace("_", " ").title(),
            "auroc": round(auroc, 3),
            "ci_lo": round(metrics.get("ci_low", 0), 3),
            "ci_hi": round(metrics.get("ci_high", 1), 3),
            "p": round(metrics.get("p_perm", 1), 4),
            "n_cells": metrics.get("n_cells_total", 0)
        })

# Sort by AUROC descending, take top per tissue
scrna_data.sort(key=lambda x: -x["auroc"])
# Keep top 5 per tissue for display
tissue_counts = {}
scrna_top = []
for d in scrna_data:
    tc = tissue_counts.get(d["tissue"], 0)
    if tc < 5:
        scrna_top.append(d)
        tissue_counts[d["tissue"]] = tc + 1

# --- Panel B: Spatial Visium (hardcoded negative result) ---
spatial_data = {
    "tissue": "Brain (RR-3 Visium)",
    "section_auroc": 0.139,
    "section_ci": [0.0, 0.429],
    "section_p": 0.900,
    "animal_auroc": 0.444,
    "animal_p": 0.400,
    "bulk_auroc": 0.000,
    "bulk_note": "Perfect inversion — genuine overfitting at n=6",
    "n_deg": 1,
    "pc1_effect": "42.5% = slide effect (not condition)",
    "conclusion": "Brain is a genuine negative tissue for spaceflight classification"
}

# --- Panel C: Power analysis ---
power_curves = []
for tissue, sim in power["simulation"].items():
    for entry in sim.get("power_curve_vary_n", []):
        power_curves.append({
            "tissue": tissue.replace("gastrocnemius", "Gastro").replace("liver", "Liver")
                    .replace("kidney", "Kidney").replace("thymus", "Thymus")
                    .replace("eye", "Eye").replace("skin", "Skin"),
            "n_per_group": entry["n_per_group"],
            "n_missions": entry["n_missions"],
            "power": round(entry["power"], 3)
        })

power_recommendations = power.get("recommendations", {})
analytical = power.get("analytical", {})
analytical_data = []
for tissue, val in analytical.items():
    label = tissue.replace("gastrocnemius", "Gastro").replace("liver", "Liver") \
                .replace("kidney", "Kidney").replace("thymus", "Thymus") \
                .replace("eye", "Eye").replace("skin", "Skin")
    analytical_data.append({
        "tissue": label,
        "auroc": round(val["observed_pooled_auroc"], 3),
        "n80": val["required_n_per_group"]["power_80"],
        "n90": val["required_n_per_group"]["power_90"],
        "n95": val["required_n_per_group"]["power_95"]
    })

# --- Panel D: Overview summary ---
overview = [
    {"category": "Bulk RNA-seq", "scope": "8 tissues × 8 methods × 4 features", "n_eval": 256, "key_result": "PCA-LR best (0.776 mean)"},
    {"category": "Ablation", "scope": "Feature count, PCA, sample size, bootstrap", "n_eval": 608, "key_result": "K=2000 optimal, PCA 20-50"},
    {"category": "Meta-analysis", "scope": "Forest plots, Friedman, DeLong, power", "n_eval": "—", "key_result": "Friedman p=0.015"},
    {"category": "SHAP", "scope": "8 tissues × 2 methods = 16 analyses", "n_eval": 16, "key_result": "5 consensus genes, NPAS2/PER2"},
    {"category": "WGCNA", "scope": "6 tissues, 60 modules, 15 preservation pairs", "n_eval": "—", "key_result": "liver→kidney Z=35.65"},
    {"category": "Cross-tissue", "scope": "7×7 transfer (Method A+C)", "n_eval": 42, "key_result": "liver→kidney 0.73"},
    {"category": "Foundation models", "scope": "UCE + scFoundation, 7 tissues", "n_eval": 14, "key_result": "Both < PCA-LR"},
    {"category": "scRNA-seq", "scope": "RRRM-2: 4 tissues, ~50 cell types", "n_eval": "~50", "key_result": "PBMC NK 0.845"},
    {"category": "Spatial", "scope": "RR-3 brain Visium", "n_eval": 1, "key_result": "NEGATIVE (AUROC=0.14)"},
    {"category": "Cross-species", "scope": "Drosophila ↔ 6 mouse tissues (KEGG)", "n_eval": 6, "key_result": "Weak negative r"},
    {"category": "Radiation", "scope": "3 tissues × 3 contrasts", "n_eval": 9, "key_result": "Near chance"},
]

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig 6 — Multi-Platform Extension</title>
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
    color: #fff; border-radius: 4px; font-size: 10px; pointer-events: none; z-index: 100; line-height: 1.4;
  }}
  .interpretation {{
    margin-top: 16px; padding: 10px; background: #f5f5f5;
    border-left: 3px solid #999; font-size: 9px; color: #444; line-height: 1.5;
  }}
  .interpretation strong {{ color: #333; }}
  table.data-table {{ border-collapse: collapse; width: 100%; font-size: 8px; margin-top: 6px; }}
  table.data-table th, table.data-table td {{ border: 1px solid #ddd; padding: 2px 5px; text-align: left; }}
  table.data-table th {{ background: #f0f0f0; font-weight: bold; }}
  .sig {{ color: #D55E00; font-weight: bold; }}
  .negative {{ color: #888; font-style: italic; }}
</style>
</head>
<body>
<div class="figure-container">
  <div class="figure-title">Fig 6 — Multi-Platform Extensions and Power Analysis</div>
  <div class="figure-subtitle">scRNA-seq (RRRM-2), spatial Visium, power recommendations, and GeneLabBench overview</div>
  <div class="panels">
    <div class="panel"><div class="panel-label">A</div>
      <div class="panel-desc">scRNA-seq LOAO Classification — RRRM-2 (Top 5 Cell Types per Tissue)</div>
      <div id="scrna"></div>
    </div>
    <div class="panel"><div class="panel-label">B</div>
      <div class="panel-desc">Spatial Visium — Brain (RR-3, NEGATIVE result)</div>
      <div id="spatial"></div>
    </div>
    <div class="panel"><div class="panel-label">C</div>
      <div class="panel-desc">Power Analysis — Required Samples per Group (6 LOMO tissues)</div>
      <div id="power"></div>
    </div>
    <div class="panel"><div class="panel-label">D</div>
      <div class="panel-desc">GeneLabBench Summary — Scope and Key Results</div>
      <div id="overview"></div>
    </div>
  </div>
  <div class="interpretation">
    <strong>Key findings:</strong>
    (A) PBMC NK cells (0.845) and T cells (0.752) show strongest scRNA-seq spaceflight signal.
    Bone marrow shows no signal (all cell types near chance).
    (B) Brain Visium is a genuine negative: section-level AUROC=0.14, PC1=slide effect.
    (C) Power analysis: thymus and gastrocnemius are well-powered; kidney needs 35+ samples per group.
    (D) GeneLabBench encompasses 1000+ evaluations across bulk, single-cell, spatial, and cross-species platforms.
  </div>
</div>
<div class="tooltip" id="tooltip" style="display:none;"></div>

<script>
const OI = {{
  black: "#000000", orange: "#E69F00", skyblue: "#56B4E9", green: "#009E73",
  yellow: "#F0E442", blue: "#0072B2", vermilion: "#D55E00", pink: "#CC79A7"
}};
const tooltip = d3.select("#tooltip");
const tissueColors = {{
  "PBMC": OI.blue, "Spleen": OI.green,
  "Bone Marrow (Femur)": OI.orange, "Bone Marrow (Humerus)": OI.pink
}};

// =============================================
// PANEL A: scRNA LOAO
// =============================================
(function() {{
  const data = {json.dumps(scrna_top)};
  const margin = {{top: 10, right: 30, bottom: 40, left: 160}};
  const rowH = 16;
  const width = 320, height = data.length * rowH;

  const svg = d3.select("#scrna").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleLinear().domain([0, 1.0]).range([0, width]);

  // 0.5 reference
  svg.append("line").attr("x1", x(0.5)).attr("x2", x(0.5)).attr("y1", 0).attr("y2", height)
    .attr("stroke", "#ccc").attr("stroke-dasharray", "3,3");

  data.forEach((d, i) => {{
    const cy = i * rowH + rowH / 2;
    const color = tissueColors[d.tissue] || "#888";

    // CI line
    svg.append("line")
      .attr("x1", x(Math.max(0, d.ci_lo))).attr("x2", x(Math.min(1, d.ci_hi)))
      .attr("y1", cy).attr("y2", cy)
      .attr("stroke", color).attr("stroke-width", 1).attr("opacity", 0.5);
    // Point
    svg.append("circle")
      .attr("cx", x(d.auroc)).attr("cy", cy).attr("r", 3)
      .attr("fill", d.p < 0.05 ? color : "#ddd")
      .attr("stroke", color).attr("stroke-width", 0.5);
    // Label
    svg.append("text").attr("x", -5).attr("y", cy + 4)
      .attr("text-anchor", "end").attr("font-size", "7px")
      .attr("fill", color)
      .text(`${{d.tissue}}: ${{d.cell_type}}`);
    // Value
    svg.append("text").attr("x", x(Math.min(1, d.ci_hi)) + 3).attr("y", cy + 3)
      .attr("font-size", "7px").attr("fill", d.p < 0.05 ? "#333" : "#aaa")
      .text(`${{d.auroc.toFixed(2)}}${{d.p < 0.01 ? "**" : d.p < 0.05 ? "*" : ""}}`);
  }});

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).ticks(5).tickFormat(d3.format(".1f")));
  svg.append("text").attr("x", width/2).attr("y", height + 32)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("AUROC (LOAO)");
}})();

// =============================================
// PANEL B: Spatial Visium (Negative)
// =============================================
(function() {{
  const data = {json.dumps(spatial_data)};
  const container = d3.select("#spatial");

  const margin = {{top: 10, right: 20, bottom: 10, left: 20}};
  const svg = container.append("svg")
    .attr("width", 460).attr("height", 180)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  // Negative result box
  svg.append("rect").attr("x", 0).attr("y", 0).attr("width", 420).attr("height", 150)
    .attr("fill", "#fff5f5").attr("stroke", "#e0c0c0").attr("rx", 6);

  svg.append("text").attr("x", 210).attr("y", 25)
    .attr("text-anchor", "middle").attr("font-size", "12px").attr("font-weight", "bold")
    .attr("fill", OI.vermilion).text("NEGATIVE RESULT");

  const lines = [
    `Section-level AUROC = ${{data.section_auroc}} [95% CI: ${{data.section_ci[0]}}, ${{data.section_ci[1]}}], p = ${{data.section_p}}`,
    `Animal-level AUROC = ${{data.animal_auroc}}, p = ${{data.animal_p}}`,
    `Companion bulk AUROC = ${{data.bulk_auroc}} (${{data.bulk_note}})`,
    `DEGs at padj<0.05: ${{data.n_deg}}`,
    `PC1 (${{data.pc1_effect}})`,
    `→ ${{data.conclusion}}`
  ];
  lines.forEach((line, i) => {{
    const isBold = i === lines.length - 1;
    svg.append("text").attr("x", 20).attr("y", 48 + i * 17)
      .attr("font-size", "9px")
      .attr("fill", isBold ? "#333" : "#555")
      .attr("font-weight", isBold ? "bold" : "normal")
      .text(line);
  }});
}})();

// =============================================
// PANEL C: Power Analysis
// =============================================
(function() {{
  const curves = {json.dumps(power_curves)};
  const analytical = {json.dumps(analytical_data)};
  const margin = {{top: 15, right: 130, bottom: 45, left: 55}};
  const width = 320, height = 220;

  const svg = d3.select("#power").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom + 80)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleLinear().domain([0, 55]).range([0, width]);
  const y = d3.scaleLinear().domain([0, 1.0]).range([height, 0]);

  // 0.8 power reference
  svg.append("line").attr("x1", 0).attr("x2", width).attr("y1", y(0.8)).attr("y2", y(0.8))
    .attr("stroke", "#ddd").attr("stroke-dasharray", "3,3");
  svg.append("text").attr("x", width + 3).attr("y", y(0.8) + 4)
    .attr("font-size", "7px").attr("fill", "#aaa").text("80% power");

  const tissueColorMap = {{
    "Liver": OI.blue, "Gastro": OI.vermilion, "Kidney": OI.green,
    "Thymus": OI.orange, "Eye": OI.pink, "Skin": OI.skyblue
  }};

  const tissues = [...new Set(curves.map(d => d.tissue))];
  tissues.forEach(tissue => {{
    const pts = curves.filter(d => d.tissue === tissue);
    if (pts.length === 0) return;

    const line = d3.line().x(d => x(d.n_per_group)).y(d => y(d.power));
    svg.append("path").datum(pts).attr("d", line)
      .attr("fill", "none").attr("stroke", tissueColorMap[tissue] || "#333")
      .attr("stroke-width", 2).attr("opacity", 0.8);
  }});

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).ticks(6));
  svg.append("text").attr("x", width/2).attr("y", height + 35)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("Samples per Group");
  svg.append("g").call(d3.axisLeft(y).ticks(5).tickFormat(d3.format(".0%")));
  svg.append("text").attr("x", -height/2).attr("y", -40)
    .attr("text-anchor", "middle").attr("font-size", "10px")
    .attr("transform", "rotate(-90)").text("Statistical Power");

  // Legend
  let ly = 0;
  tissues.forEach(t => {{
    svg.append("line").attr("x1", width+10).attr("x2", width+30)
      .attr("y1", ly+4).attr("y2", ly+4)
      .attr("stroke", tissueColorMap[t] || "#333").attr("stroke-width", 2);
    svg.append("text").attr("x", width+33).attr("y", ly+8)
      .attr("font-size", "8px").text(t);
    ly += 14;
  }});

  // Analytical table below
  const tableY = height + 50;
  svg.append("text").attr("x", 0).attr("y", tableY)
    .attr("font-size", "9px").attr("font-weight", "bold").text("Required n per group:");
  analytical.forEach((d, i) => {{
    svg.append("text").attr("x", 0).attr("y", tableY + 14 + i * 12)
      .attr("font-size", "8px").attr("fill", tissueColorMap[d.tissue] || "#333")
      .text(`${{d.tissue}} (AUROC=${{d.auroc}}): n=${{d.n80}} (80%), n=${{d.n90}} (90%)`);
  }});
}})();

// =============================================
// PANEL D: GeneLabBench Overview Table
// =============================================
(function() {{
  const data = {json.dumps(overview)};
  const container = d3.select("#overview");

  let html = '<table class="data-table"><tr><th>Category</th><th>Scope</th><th>N Eval</th><th>Key Result</th></tr>';
  data.forEach(d => {{
    html += `<tr><td><b>${{d.category}}</b></td><td>${{d.scope}}</td><td>${{d.n_eval}}</td><td>${{d.key_result}}</td></tr>`;
  }});
  html += '</table>';
  html += '<div style="font-size:8px; color:#666; margin-top:8px;">Total: 1000+ evaluations across 8 tissues, 8 methods, 4 features, single-cell, spatial, and cross-species platforms.</div>';

  container.append("div").html(html);
}})();
</script>
</body>
</html>"""

out_path = os.path.join(OUT_DIR, "Fig6_multiplatform.html")
with open(out_path, "w") as f:
    f.write(html)
print(f"Written: {out_path}")
print(f"  Panel A: {len(scrna_top)} scRNA cell type entries")
print(f"  Panel C: {len(power_curves)} power curve points")
print(f"  Panel D: {len(overview)} overview rows")

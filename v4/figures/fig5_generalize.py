#!/usr/bin/env python3
"""
Fig 5 — Generalizability (Cross-tissue, FM comparison, Cross-species, Radiation)
GeneLabBench v4 Publication Figure

Panels:
  A: 7×7 cross-tissue transfer heatmap (Method A gene)
  B: FM vs PCA-LR bar chart (UCE + scFoundation, 7 tissues)
  C: Cross-species correlation (Drosophila vs mouse tissues)
  D: Radiation analog classification (brain, spleen, skin)

Output: v4/figures/html/Fig5_generalize.html
"""

import json
import os
import math

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
V3_EVAL = os.path.join(os.path.dirname(BASE), "v3", "evaluation")
V4_EVAL = os.path.join(BASE, "evaluation")
OUT_DIR = os.path.join(BASE, "figures", "html")
os.makedirs(OUT_DIR, exist_ok=True)

with open(os.path.join(V3_EVAL, "B_ext_transfer_matrix.json")) as f:
    transfer = json.load(f)
with open(os.path.join(V3_EVAL, "FM_scfoundation.json")) as f:
    fm_scf = json.load(f)
with open(os.path.join(V3_EVAL, "FM_uce.json")) as f:
    fm_uce = json.load(f)
with open(os.path.join(V3_EVAL, "E4_multispecies_nes.json")) as f:
    e4 = json.load(f)
with open(os.path.join(V3_EVAL, "R1_radiation_classification.json")) as f:
    r1 = json.load(f)
with open(os.path.join(V4_EVAL, "M1_summary.json")) as f:
    m1 = json.load(f)

TISSUE_LABELS = {
    "liver": "Liver", "kidney": "Kidney", "thymus": "Thymus",
    "gastrocnemius": "Gastro", "eye": "Eye", "lung": "Lung", "colon": "Colon", "skin": "Skin"
}

# --- Panel A: Cross-tissue transfer ---
transfer_tissues = transfer["tissues"]
transfer_labels = [TISSUE_LABELS.get(t, t) for t in transfer_tissues]
matrix_a = transfer["method_a_matrix"]
# Flatten to array of objects
transfer_data = []
for i, t_train in enumerate(transfer_tissues):
    for j, t_test in enumerate(transfer_tissues):
        val = matrix_a[i][j]
        transfer_data.append({
            "train": TISSUE_LABELS.get(t_train, t_train),
            "test": TISSUE_LABELS.get(t_test, t_test),
            "auroc": round(val, 3) if val is not None else None
        })

# --- Panel B: FM comparison ---
fm_data = []
for tissue in ["liver", "gastrocnemius", "kidney", "thymus", "eye", "lung", "colon"]:
    # PCA-LR from v4 M1
    pca_auroc = m1[tissue]["gene"]["pca_lr"]["auroc"]
    pca_p = m1[tissue]["gene"]["pca_lr"]["perm_p"]

    # scFoundation
    scf_entry = next((r for r in fm_scf["results"] if r["tissue"] == tissue), None)
    scf_auroc = scf_entry["lomo_auroc"] if scf_entry else None
    scf_p = scf_entry.get("perm_p") if scf_entry else None

    # UCE
    uce_entry = next((r for r in fm_uce["results"] if r["tissue"] == tissue), None)
    uce_auroc = uce_entry["lomo_auroc"] if uce_entry else None
    uce_p = uce_entry.get("perm_p") if uce_entry else None

    fm_data.append({
        "tissue": TISSUE_LABELS.get(tissue, tissue),
        "pca_lr": round(pca_auroc, 3),
        "pca_p": round(pca_p, 4) if pca_p else None,
        "scf": round(scf_auroc, 3) if scf_auroc else None,
        "scf_p": round(scf_p, 4) if scf_p else None,
        "uce": round(uce_auroc, 3) if uce_auroc else None,
        "uce_p": round(uce_p, 4) if uce_p else None
    })

# --- Panel C: Cross-species (Drosophila vs mouse tissues) ---
cross_species = []
for pc in e4["pairwise_concordance"]:
    a, b = pc["species_a"], pc["species_b"]
    # Drosophila vs mouse pairs
    if "drosophila" in a.lower() or "drosophila" in b.lower():
        mouse_tissue = b if "drosophila" in a.lower() else a
        mouse_label = mouse_tissue.replace("mouse_", "")
        cross_species.append({
            "mouse_tissue": TISSUE_LABELS.get(mouse_label, mouse_label),
            "r": round(pc["spearman_r"], 3),
            "p": round(pc["spearman_p"], 4),
            "n_common": pc["n_common_pathways"],
            "ci": [round(pc["bootstrap_ci"][0], 3), round(pc["bootstrap_ci"][1], 3)]
        })

# Mouse intra-species pairs
mouse_pairs = []
for pc in e4["pairwise_concordance"]:
    a, b = pc["species_a"], pc["species_b"]
    if "mouse" in a and "mouse" in b:
        ta = a.replace("mouse_", "")
        tb = b.replace("mouse_", "")
        mouse_pairs.append({
            "tissue_a": TISSUE_LABELS.get(ta, ta),
            "tissue_b": TISSUE_LABELS.get(tb, tb),
            "r": round(pc["spearman_r"], 3),
            "p": round(pc["spearman_p"], 4)
        })

# --- Panel D: Radiation analog ---
rad_data = []
for tissue in r1:
    for entry in r1[tissue]:
        rad_data.append({
            "tissue": tissue.replace("_rad", ""),
            "contrast": entry["contrast"],
            "auroc": round(entry["auroc"], 3),
            "ci_lo": round(entry.get("ci_lower", 0), 3),
            "ci_hi": round(entry.get("ci_upper", 1), 3),
            "p": round(entry.get("perm_p", 1), 3),
            "n": entry.get("n", 0)
        })

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig 5 — Generalizability</title>
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
</style>
</head>
<body>
<div class="figure-container">
  <div class="figure-title">Fig 5 — Generalizability of Spaceflight Transcriptomic Signatures</div>
  <div class="figure-subtitle">Cross-tissue transfer, foundation model comparison, cross-species concordance, and radiation analogs</div>
  <div class="panels">
    <div class="panel"><div class="panel-label">A</div>
      <div class="panel-desc">Cross-Tissue Transfer — Method A (Gene Features), 7×7 Tissues</div>
      <div id="transfer"></div>
    </div>
    <div class="panel"><div class="panel-label">B</div>
      <div class="panel-desc">Foundation Model vs PCA-LR — AUROC Comparison (7 tissues)</div>
      <div id="fm_compare"></div>
    </div>
    <div class="panel"><div class="panel-label">C</div>
      <div class="panel-desc">Cross-Species KEGG Concordance — Drosophila vs Mouse Tissues</div>
      <div id="cross_species"></div>
    </div>
    <div class="panel"><div class="panel-label">D</div>
      <div class="panel-desc">Radiation Analog Classification (Brain, Spleen, Skin)</div>
      <div id="radiation"></div>
    </div>
  </div>
  <div class="interpretation">
    <strong>Key findings:</strong>
    (A) Cross-tissue transfer is heterogeneous: liver→kidney best (0.73), while most pairs are near chance.
    (B) Both foundation models (scFoundation, UCE) underperform PCA-LR on all tissues — pre-trained cell
    atlas knowledge does not capture spaceflight-specific perturbations.
    (C) Drosophila-mouse KEGG concordance is weak and often negative (r = −0.19 to −0.59), suggesting
    spaceflight responses are species-specific at the pathway level.
    (D) Radiation analog classification is near chance for all tissues, indicating radiation alone does not
    recapitulate the full spaceflight transcriptomic signature.
  </div>
</div>
<div class="tooltip" id="tooltip" style="display:none;"></div>

<script>
const OI = {{
  black: "#000000", orange: "#E69F00", skyblue: "#56B4E9", green: "#009E73",
  yellow: "#F0E442", blue: "#0072B2", vermilion: "#D55E00", pink: "#CC79A7"
}};
const tooltip = d3.select("#tooltip");

// =============================================
// PANEL A: Cross-Tissue Transfer Heatmap
// =============================================
(function() {{
  const data = {json.dumps(transfer_data)};
  const labels = {json.dumps(transfer_labels)};
  const margin = {{top: 5, right: 60, bottom: 60, left: 70}};
  const cellS = 52;
  const width = labels.length * cellS, height = labels.length * cellS;

  const svg = d3.select("#transfer").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(labels).range([0, width]).padding(0.05);
  const y = d3.scaleBand().domain(labels).range([0, height]).padding(0.05);
  const color = d3.scaleSequential().domain([0.3, 0.85]).interpolator(d3.interpolateBlues);

  svg.selectAll("rect.cell")
    .data(data.filter(d => d.auroc !== null)).enter().append("rect")
    .attr("x", d => x(d.test)).attr("y", d => y(d.train))
    .attr("width", x.bandwidth()).attr("height", y.bandwidth())
    .attr("fill", d => d.train === d.test ? "#f0f0f0" : color(Math.max(0.3, d.auroc)))
    .attr("stroke", "#fff").attr("stroke-width", 1).attr("rx", 2)
    .on("mouseover", function(event, d) {{
      tooltip.style("display", "block")
        .html(`Train: <b>${{d.train}}</b> → Test: <b>${{d.test}}</b><br>AUROC: ${{d.auroc}}`);
    }})
    .on("mousemove", function(event) {{
      tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
    }})
    .on("mouseout", () => tooltip.style("display", "none"));

  // AUROC text
  svg.selectAll("text.val")
    .data(data.filter(d => d.auroc !== null && d.train !== d.test)).enter().append("text")
    .attr("x", d => x(d.test) + x.bandwidth()/2)
    .attr("y", d => y(d.train) + y.bandwidth()/2 + 4)
    .attr("text-anchor", "middle").attr("font-size", "7px")
    .attr("fill", d => d.auroc > 0.65 ? "#fff" : "#333")
    .text(d => d.auroc.toFixed(2));

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).tickSize(0))
    .selectAll("text").attr("transform", "rotate(-40)")
    .style("text-anchor", "end").style("font-size", "8px");
  svg.append("g").call(d3.axisLeft(y).tickSize(0))
    .selectAll("text").style("font-size", "8px");

  svg.append("text").attr("x", width/2).attr("y", height + 52)
    .attr("text-anchor", "middle").attr("font-size", "9px").text("Test Tissue →");
  svg.append("text").attr("x", -height/2).attr("y", -55)
    .attr("text-anchor", "middle").attr("font-size", "9px")
    .attr("transform", "rotate(-90)").text("← Train Tissue");
}})();

// =============================================
// PANEL B: FM vs PCA-LR
// =============================================
(function() {{
  const data = {json.dumps(fm_data)};
  const margin = {{top: 15, right: 20, bottom: 50, left: 50}};
  const width = 420, height = 280;
  const barW = 14;

  const svg = d3.select("#fm_compare").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(data.map(d => d.tissue)).range([0, width]).padding(0.25);
  const y = d3.scaleLinear().domain([0, 1.0]).range([height, 0]);

  // Reference line at 0.5
  svg.append("line").attr("x1", 0).attr("x2", width).attr("y1", y(0.5)).attr("y2", y(0.5))
    .attr("stroke", "#ccc").attr("stroke-dasharray", "3,3");

  data.forEach(d => {{
    const cx = x(d.tissue) + x.bandwidth() / 2;
    // PCA-LR bar
    svg.append("rect")
      .attr("x", cx - barW * 1.5 - 1).attr("y", y(d.pca_lr))
      .attr("width", barW).attr("height", height - y(d.pca_lr))
      .attr("fill", OI.blue).attr("opacity", 0.8);
    // scFoundation bar
    if (d.scf !== null) {{
      svg.append("rect")
        .attr("x", cx - barW * 0.5).attr("y", y(d.scf))
        .attr("width", barW).attr("height", height - y(d.scf))
        .attr("fill", OI.vermilion).attr("opacity", 0.8);
    }}
    // UCE bar
    if (d.uce !== null) {{
      svg.append("rect")
        .attr("x", cx + barW * 0.5 + 1).attr("y", y(d.uce))
        .attr("width", barW).attr("height", height - y(d.uce))
        .attr("fill", OI.green).attr("opacity", 0.8);
    }}
    // Significance stars
    if (d.scf_p && d.scf_p < 0.05) {{
      svg.append("text").attr("x", cx).attr("y", y(d.scf) - 3)
        .attr("text-anchor", "middle").attr("font-size", "9px").attr("fill", OI.vermilion).text("*");
    }}
  }});

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).tickSize(0))
    .selectAll("text").attr("transform", "rotate(-30)")
    .style("text-anchor", "end").style("font-size", "9px");
  svg.append("g").call(d3.axisLeft(y).ticks(5).tickFormat(d3.format(".1f")));
  svg.append("text").attr("x", -height/2).attr("y", -38)
    .attr("text-anchor", "middle").attr("font-size", "10px")
    .attr("transform", "rotate(-90)").text("AUROC");

  // Legend
  const lx = width - 130, ly = 10;
  [["PCA-LR", OI.blue], ["scFoundation", OI.vermilion], ["UCE", OI.green]].forEach((item, i) => {{
    svg.append("rect").attr("x", lx).attr("y", ly + i*16).attr("width", 12).attr("height", 10).attr("fill", item[1]).attr("opacity", 0.8);
    svg.append("text").attr("x", lx + 16).attr("y", ly + i*16 + 9).attr("font-size", "9px").text(item[0]);
  }});
}})();

// =============================================
// PANEL C: Cross-species
// =============================================
(function() {{
  const drosData = {json.dumps(cross_species)};
  const margin = {{top: 15, right: 20, bottom: 45, left: 100}};
  const width = 350, height = 220;

  const svg = d3.select("#cross_species").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const y = d3.scaleBand().domain(drosData.map(d => d.mouse_tissue)).range([0, height]).padding(0.3);
  const x = d3.scaleLinear().domain([-0.7, 0.5]).range([0, width]);

  // Zero line
  svg.append("line").attr("x1", x(0)).attr("x2", x(0)).attr("y1", 0).attr("y2", height)
    .attr("stroke", "#ccc").attr("stroke-dasharray", "3,3");

  drosData.forEach(d => {{
    const cy = y(d.mouse_tissue) + y.bandwidth() / 2;
    // CI line
    svg.append("line")
      .attr("x1", x(d.ci[0])).attr("x2", x(d.ci[1]))
      .attr("y1", cy).attr("y2", cy)
      .attr("stroke", "#888").attr("stroke-width", 1);
    // Point
    svg.append("circle")
      .attr("cx", x(d.r)).attr("cy", cy)
      .attr("r", 5)
      .attr("fill", d.p < 0.05 ? OI.vermilion : "#aaa")
      .attr("stroke", "#333").attr("stroke-width", 0.5);
    // Value
    svg.append("text").attr("x", x(d.ci[1]) + 5).attr("y", cy + 4)
      .attr("font-size", "8px").attr("fill", d.p < 0.05 ? OI.vermilion : "#666")
      .text(`r=${{d.r}} ${{d.p < 0.01 ? "**" : d.p < 0.05 ? "*" : ""}}`);
  }});

  svg.append("g").call(d3.axisLeft(y).tickSize(0)).selectAll("text").style("font-size", "9px");
  svg.append("g").attr("transform", `translate(0,${{height}})`).call(d3.axisBottom(x).ticks(6));
  svg.append("text").attr("x", width/2).attr("y", height + 35)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("Spearman ρ (KEGG NES)");
  svg.append("text").attr("x", 0).attr("y", -5)
    .attr("font-size", "9px").attr("fill", OI.vermilion).text("Drosophila ↔ Mouse");
}})();

// =============================================
// PANEL D: Radiation Analog
// =============================================
(function() {{
  const data = {json.dumps(rad_data)};
  const margin = {{top: 15, right: 20, bottom: 45, left: 130}};
  const width = 300, height = 260;

  const svg = d3.select("#radiation").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const y = d3.scaleBand().domain(data.map(d => d.tissue + ": " + d.contrast.split(" ")[0]))
    .range([0, height]).padding(0.2);
  const x = d3.scaleLinear().domain([0, 1.0]).range([0, width]);

  // 0.5 reference
  svg.append("line").attr("x1", x(0.5)).attr("x2", x(0.5)).attr("y1", 0).attr("y2", height)
    .attr("stroke", "#ccc").attr("stroke-dasharray", "3,3");

  data.forEach(d => {{
    const label = d.tissue + ": " + d.contrast.split(" ")[0];
    const cy = y(label) + y.bandwidth() / 2;

    // CI
    svg.append("line")
      .attr("x1", x(d.ci_lo)).attr("x2", x(d.ci_hi))
      .attr("y1", cy).attr("y2", cy)
      .attr("stroke", "#888").attr("stroke-width", 1);
    // Point
    svg.append("circle")
      .attr("cx", x(d.auroc)).attr("cy", cy)
      .attr("r", 5)
      .attr("fill", d.p < 0.05 ? OI.vermilion : "#aaa");
    // Value
    svg.append("text").attr("x", x(d.ci_hi) + 5).attr("y", cy + 4)
      .attr("font-size", "8px").text(d.auroc.toFixed(2));
  }});

  svg.append("g").call(d3.axisLeft(y).tickSize(0)).selectAll("text").style("font-size", "7px");
  svg.append("g").attr("transform", `translate(0,${{height}})`).call(d3.axisBottom(x).ticks(5));
  svg.append("text").attr("x", width/2).attr("y", height + 35)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("AUROC");
}})();
</script>
</body>
</html>"""

out_path = os.path.join(OUT_DIR, "Fig5_generalize.html")
with open(out_path, "w") as f:
    f.write(html)
print(f"Written: {out_path}")
print(f"  Panel A: {len(transfer_data)} transfer pairs")
print(f"  Panel B: {len(fm_data)} tissues × 3 methods")
print(f"  Panel C: {len(cross_species)} Drosophila-mouse pairs")
print(f"  Panel D: {len(rad_data)} radiation entries")

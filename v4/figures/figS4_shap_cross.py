#!/usr/bin/env python3
"""
Fig S4 — SHAP Cross-Method Agreement & Interaction Analysis
A: Method agreement (Jaccard + Spearman) per tissue
B: SHAP top-100 rank comparison heatmap
C: SHAP direction bias (up vs down in FLT)
D: Gene-gene interaction top pairs (tree methods)

Output: v4/figures/html/FigS4_shap_cross.html
"""

import json
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EVAL_DIR = os.path.join(BASE, "evaluation")
OUT_DIR = os.path.join(BASE, "figures", "html")
os.makedirs(OUT_DIR, exist_ok=True)

# Load consensus
with open(os.path.join(EVAL_DIR, "SHAP_consensus.json")) as f:
    consensus = json.load(f)

# Load individual SHAP files to build cross-method heatmap
ALL_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin", "lung", "colon"]
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastro", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin", "lung": "Lung", "colon": "Colon"
}

# Discover available SHAP files
shap_files = {}
for fname in os.listdir(EVAL_DIR):
    if fname.startswith("SHAP_") and not fname.startswith("SHAP_interactions") \
       and not fname.startswith("SHAP_consensus") and not fname.startswith("SHAP_WGCNA") \
       and fname.endswith(".json"):
        with open(os.path.join(EVAL_DIR, fname)) as f:
            sd = json.load(f)
        tissue = sd["tissue"]
        method = sd["method"]
        if tissue not in shap_files:
            shap_files[tissue] = {}
        shap_files[tissue][method] = {
            "top_100": sd["top_100_genes"],
            "method_label": sd.get("method_label", method)
        }

# Build cross-method Jaccard matrix per tissue
cross_method_data = {}
for tissue in ALL_TISSUES:
    if tissue not in shap_files or len(shap_files[tissue]) < 2:
        continue
    methods = sorted(shap_files[tissue].keys())
    method_labels = [shap_files[tissue][m]["method_label"] for m in methods]
    jaccard_matrix = []
    for ma in methods:
        row = []
        for mb in methods:
            if ma == mb:
                row.append(1.0)
            else:
                sa = set(shap_files[tissue][ma]["top_100"])
                sb = set(shap_files[tissue][mb]["top_100"])
                inter = len(sa & sb)
                union = len(sa | sb)
                row.append(round(inter / union, 3) if union > 0 else 0)
        jaccard_matrix.append(row)
    cross_method_data[tissue] = {
        "methods": method_labels,
        "jaccard": jaccard_matrix
    }

# Method agreement from consensus
agreement_data = consensus.get("method_agreement", {})

# SHAP directions from consensus
direction_data = consensus.get("shap_directions", {})

# Load interaction files
interaction_files = []
for fname in sorted(os.listdir(EVAL_DIR)):
    if fname.startswith("SHAP_interactions_") and fname.endswith(".json"):
        with open(os.path.join(EVAL_DIR, fname)) as f:
            idata = json.load(f)
        interaction_files.append({
            "tissue": TISSUE_LABELS.get(idata["tissue"], idata["tissue"]),
            "method": idata["method"],
            "top_20_interactions": idata.get("top_50_interactions", [])[:20],
            "top_20_main": idata.get("top_20_main_effects", {})
        })

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig S4 — SHAP Cross-Method Agreement</title>
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
  <div class="figure-title">Fig S4 — SHAP Cross-Method Agreement & Gene-Gene Interactions</div>
  <div class="figure-subtitle">Top-100 gene overlap, direction consistency, and interaction effects</div>
  <div class="panels" id="panels"></div>
  <div class="note">
    A: Jaccard index and Spearman &rho; for top-100 SHAP genes between method pairs.
    B: Pairwise Jaccard heatmaps per tissue (2 methods per tissue, top-2 by AUROC).
    C: SHAP direction bias — fraction of top-100 genes upregulated (SHAP &gt; 0) vs downregulated in flight.
    D: Top gene-gene interactions from TreeSHAP (XGBoost/RF only).
  </div>
</div>
<div class="tooltip" id="tooltip" style="display:none;"></div>

<script>
const OI = {{
  black: "#000000", orange: "#E69F00", skyblue: "#56B4E9", green: "#009E73",
  yellow: "#F0E442", blue: "#0072B2", vermilion: "#D55E00", pink: "#CC79A7"
}};
const tooltip = d3.select("#tooltip");
const agreementData = {json.dumps(agreement_data)};
const directionData = {json.dumps(direction_data)};
const crossMethodData = {json.dumps(cross_method_data)};
const interactionFiles = {json.dumps(interaction_files)};
const tissueLabels = {json.dumps(TISSUE_LABELS)};

// --- Panel A: Method Agreement Summary ---
(() => {{
  const panel = d3.select("#panels").append("div").attr("class", "panel");
  panel.append("div").attr("class", "panel-label").text("A. Cross-Method Agreement (Jaccard & Spearman)");

  const tissues = Object.keys(agreementData);
  const margin = {{top: 10, right: 20, bottom: 45, left: 80}};
  const barH = 22;
  const height = tissues.length * barH * 2 + 10;
  const width = 400;

  const svg = panel.append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const xJac = d3.scaleLinear().domain([0, 0.5]).range([0, width/2 - 10]);
  const xRho = d3.scaleLinear().domain([0, 0.5]).range([0, width/2 - 10]);

  tissues.forEach((tissue, i) => {{
    const ag = agreementData[tissue];
    const tLabel = tissueLabels[tissue] || tissue;
    const y = i * barH * 2;

    // Tissue label
    svg.append("text").attr("x", -5).attr("y", y + barH).attr("text-anchor", "end")
      .attr("font-size", "8px").attr("dominant-baseline", "middle")
      .text(`${{tLabel}} (${{ag.method_a}} vs ${{ag.method_b}})`);

    // Jaccard bar
    svg.append("rect").attr("x", 0).attr("y", y).attr("width", xJac(ag.jaccard_top100))
      .attr("height", barH - 2).attr("fill", OI.blue).attr("opacity", 0.7).attr("rx", 2);
    svg.append("text").attr("x", xJac(ag.jaccard_top100) + 3).attr("y", y + barH/2)
      .attr("font-size", "7px").attr("dominant-baseline", "middle")
      .text(`J=${{ag.jaccard_top100.toFixed(2)}} (${{ag.n_overlap_top100}} genes)`);

    // Spearman bar
    const rhoVal = Math.abs(ag.spearman_rho);
    svg.append("rect").attr("x", width/2 + 5).attr("y", y).attr("width", xRho(rhoVal))
      .attr("height", barH - 2).attr("fill", OI.orange).attr("opacity", 0.7).attr("rx", 2);
    svg.append("text").attr("x", width/2 + 5 + xRho(rhoVal) + 3).attr("y", y + barH/2)
      .attr("font-size", "7px").attr("dominant-baseline", "middle")
      .text(`\\u03C1=${{ag.spearman_rho.toFixed(3)}}`);
  }});

  // Column headers
  svg.append("text").attr("x", width/4).attr("y", -2).attr("text-anchor", "middle")
    .attr("font-size", "8px").attr("font-weight", "bold").text("Jaccard (top-100)");
  svg.append("text").attr("x", width*3/4).attr("y", -2).attr("text-anchor", "middle")
    .attr("font-size", "8px").attr("font-weight", "bold").text("Spearman \\u03C1 (all genes)");
}})();

// --- Panel B: Jaccard heatmaps per tissue (only tissues with 2 methods) ---
(() => {{
  const panel = d3.select("#panels").append("div").attr("class", "panel");
  panel.append("div").attr("class", "panel-label").text("B. Method Pair Jaccard Heatmaps");

  const tissues = Object.keys(crossMethodData);
  const cellS = 35;
  const margin = {{top: 5, right: 5, bottom: 40, left: 55}};

  tissues.forEach((tissue, tidx) => {{
    const td = crossMethodData[tissue];
    const methods = td.methods;
    const tLabel = tissueLabels[tissue] || tissue;

    panel.append("div").style("font-size", "9px").style("font-weight", "bold")
      .style("margin-top", tidx > 0 ? "8px" : "0").text(tLabel);

    const width = methods.length * cellS;
    const height = methods.length * cellS;

    const svg = panel.append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

    const x = d3.scaleBand().domain(methods).range([0, width]).padding(0.05);
    const y = d3.scaleBand().domain(methods).range([0, height]).padding(0.05);
    const color = d3.scaleSequential(d3.interpolateBlues).domain([0, 1]);

    methods.forEach((ma, i) => {{
      methods.forEach((mb, j) => {{
        const val = td.jaccard[i][j];
        svg.append("rect")
          .attr("x", x(mb)).attr("y", y(ma))
          .attr("width", x.bandwidth()).attr("height", y.bandwidth())
          .attr("fill", i === j ? "#f0f0f0" : color(val))
          .attr("stroke", "#eee").attr("rx", 2)
          .on("mouseover", function(event) {{
            tooltip.style("display", "block")
              .html(`<b>${{ma}} vs ${{mb}}</b><br>Jaccard=${{val}}`);
          }})
          .on("mousemove", function(event) {{
            tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
          }})
          .on("mouseout", () => tooltip.style("display", "none"));

        if (i !== j) {{
          svg.append("text")
            .attr("x", x(mb) + x.bandwidth()/2).attr("y", y(ma) + y.bandwidth()/2 + 3)
            .attr("text-anchor", "middle").attr("font-size", "7px")
            .attr("fill", val > 0.5 ? "#fff" : "#333")
            .text(val.toFixed(2));
        }}
      }});
    }});

    svg.append("g").attr("transform", `translate(0,${{height}})`)
      .call(d3.axisBottom(x).tickSize(0))
      .selectAll("text").attr("transform", "rotate(-30)").style("text-anchor", "end").style("font-size", "7px");
    svg.append("g").call(d3.axisLeft(y).tickSize(0))
      .selectAll("text").style("font-size", "7px");
  }});
}})();

// --- Panel C: SHAP direction bias ---
(() => {{
  const panel = d3.select("#panels").append("div").attr("class", "panel");
  panel.append("div").attr("class", "panel-label").text("C. SHAP Direction Bias (Top-100 Genes)");

  const entries = Object.entries(directionData).map(([key, val]) => ({{
    label: `${{(tissueLabels[val.tissue] || val.tissue)}} (${{val.method}})`,
    tissue: val.tissue,
    method: val.method,
    n_up: val.n_positive_in_top100,
    n_down: val.n_negative_in_top100,
    frac_up: val.n_positive_in_top100 / 100
  }}));
  entries.sort((a, b) => b.frac_up - a.frac_up);

  const margin = {{top: 10, right: 10, bottom: 30, left: 120}};
  const barH = 14;
  const width = 350;
  const height = entries.length * barH;

  const svg = panel.append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const xScale = d3.scaleLinear().domain([0, 100]).range([0, width]);
  const yBand = d3.scaleBand().domain(entries.map(e => e.label)).range([0, height]).padding(0.1);

  entries.forEach(e => {{
    // Up bar (positive SHAP = upregulated in FLT)
    svg.append("rect")
      .attr("x", 0).attr("y", yBand(e.label))
      .attr("width", xScale(e.n_up)).attr("height", yBand.bandwidth())
      .attr("fill", OI.vermilion).attr("opacity", 0.7);

    // Down bar
    svg.append("rect")
      .attr("x", xScale(e.n_up)).attr("y", yBand(e.label))
      .attr("width", xScale(e.n_down)).attr("height", yBand.bandwidth())
      .attr("fill", OI.blue).attr("opacity", 0.7);

    // Count labels
    svg.append("text")
      .attr("x", xScale(e.n_up / 2)).attr("y", yBand(e.label) + yBand.bandwidth()/2 + 3)
      .attr("text-anchor", "middle").attr("font-size", "6px").attr("fill", "#fff")
      .text(e.n_up > 15 ? e.n_up : "");
    svg.append("text")
      .attr("x", xScale(e.n_up + e.n_down / 2)).attr("y", yBand(e.label) + yBand.bandwidth()/2 + 3)
      .attr("text-anchor", "middle").attr("font-size", "6px").attr("fill", "#fff")
      .text(e.n_down > 15 ? e.n_down : "");
  }});

  svg.append("g").call(d3.axisLeft(yBand).tickSize(0))
    .selectAll("text").style("font-size", "7px");
  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(xScale).ticks(5))
    .selectAll("text").style("font-size", "7px");

  // Legend
  svg.append("rect").attr("x", width - 100).attr("y", -8).attr("width", 8).attr("height", 8).attr("fill", OI.vermilion).attr("opacity", 0.7);
  svg.append("text").attr("x", width - 88).attr("y", -1).attr("font-size", "7px").text("Up in FLT");
  svg.append("rect").attr("x", width - 45).attr("y", -8).attr("width", 8).attr("height", 8).attr("fill", OI.blue).attr("opacity", 0.7);
  svg.append("text").attr("x", width - 33).attr("y", -1).attr("font-size", "7px").text("Down");
}})();

// --- Panel D: Gene-gene interactions (TreeSHAP) ---
(() => {{
  const panel = d3.select("#panels").append("div").attr("class", "panel");
  panel.append("div").attr("class", "panel-label").text("D. Top Gene-Gene Interactions (TreeSHAP)");

  if (interactionFiles.length === 0) {{
    panel.append("div").style("font-size", "9px").style("color", "#999")
      .text("No interaction data available.");
    return;
  }}

  interactionFiles.forEach((ifile, fidx) => {{
    panel.append("div").style("font-size", "9px").style("font-weight", "bold")
      .style("margin-top", fidx > 0 ? "10px" : "0")
      .text(`${{ifile.tissue}} — ${{ifile.method.toUpperCase()}}`);

    const top10 = ifile.top_20_interactions.slice(0, 10);
    const margin = {{top: 5, right: 10, bottom: 25, left: 120}};
    const barH = 14;
    const width = 300;
    const height = top10.length * barH;

    const svg = panel.append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
      .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

    const xMax = d3.max(top10, d => d.interaction_strength) * 1.1;
    const xScale = d3.scaleLinear().domain([0, xMax]).range([0, width]);
    const yBand = d3.scaleBand()
      .domain(top10.map(d => `${{d.gene_a}} \\u00D7 ${{d.gene_b}}`))
      .range([0, height]).padding(0.1);

    top10.forEach(d => {{
      const label = `${{d.gene_a}} \\u00D7 ${{d.gene_b}}`;
      svg.append("rect")
        .attr("x", 0).attr("y", yBand(label))
        .attr("width", xScale(d.interaction_strength)).attr("height", yBand.bandwidth())
        .attr("fill", OI.green).attr("opacity", 0.7).attr("rx", 1)
        .on("mouseover", function(event) {{
          tooltip.style("display", "block")
            .html(`<b>${{d.gene_a}} \\u00D7 ${{d.gene_b}}</b><br>Interaction: ${{d.interaction_strength.toFixed(4)}}`);
        }})
        .on("mousemove", function(event) {{
          tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
        }})
        .on("mouseout", () => tooltip.style("display", "none"));

      svg.append("text")
        .attr("x", xScale(d.interaction_strength) + 3).attr("y", yBand(label) + yBand.bandwidth()/2 + 3)
        .attr("font-size", "6px").text(d.interaction_strength.toFixed(4));
    }});

    svg.append("g").call(d3.axisLeft(yBand).tickSize(0))
      .selectAll("text").style("font-size", "6px");
    svg.append("g").attr("transform", `translate(0,${{height}})`)
      .call(d3.axisBottom(xScale).ticks(4).tickFormat(d3.format(".3f")))
      .selectAll("text").style("font-size", "7px");
  }});
}})();
</script>
</body>
</html>"""

out_path = os.path.join(OUT_DIR, "FigS4_shap_cross.html")
with open(out_path, "w") as f:
    f.write(html)
print(f"Written: {out_path}")
print(f"  {len(agreement_data)} tissue agreements, {len(direction_data)} direction entries, {len(interaction_files)} interaction files")

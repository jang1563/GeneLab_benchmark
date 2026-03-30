#!/usr/bin/env python3
"""
Fig 3 — Consensus Spaceflight Gene Signature
GeneLabBench v4 Publication Figure

Panels:
  A: UpSet-style intersection plot (SHAP top-100 overlap across tissues)
  B: Consensus + tissue-specific genes × tissue heatmap (SHAP importance)
  C: Method agreement (Jaccard & Spearman per tissue)
  D: Cell 2020 validation + SHAP direction summary

Output: v4/figures/html/Fig3_consensus.html (self-contained D3.js v7)
"""

import json
import os
import glob

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EVAL_DIR = os.path.join(BASE, "evaluation")
OUT_DIR = os.path.join(BASE, "figures", "html")
os.makedirs(OUT_DIR, exist_ok=True)

with open(os.path.join(EVAL_DIR, "SHAP_consensus.json")) as f:
    consensus = json.load(f)

# Load per-tissue SHAP files for the gene × tissue heatmap
shap_files = glob.glob(os.path.join(EVAL_DIR, "SHAP_*_*.json"))
shap_files = [f for f in shap_files if "consensus" not in f and "WGCNA" not in f and "interactions" not in f]

per_tissue_shap = {}  # tissue -> {gene: importance}
per_tissue_top100 = {}  # tissue -> set of top-100 genes
for sf in shap_files:
    with open(sf) as f:
        d = json.load(f)
    tissue = d["tissue"]
    method = d["method"]
    # Use first method per tissue (or accumulate)
    if tissue not in per_tissue_shap:
        per_tissue_shap[tissue] = {}
        per_tissue_top100[tissue] = set()
    # Merge: take max importance per gene across methods
    for gene, imp in d["gene_importance"].items():
        gene_upper = gene.upper()
        if gene_upper not in per_tissue_shap[tissue] or imp > per_tissue_shap[tissue][gene_upper]:
            per_tissue_shap[tissue][gene_upper] = imp
    per_tissue_top100[tissue].update([g.upper() for g in d["top_100_genes"]])

TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin", "lung", "colon"]
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastrocnemius", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin", "lung": "Lung", "colon": "Colon"
}

# --- Panel A: UpSet-style data ---
# Count how many tissues each gene appears in (top-100 union)
gene_tissue_membership = {}
for tissue in TISSUES:
    if tissue in per_tissue_top100:
        for gene in per_tissue_top100[tissue]:
            if gene not in gene_tissue_membership:
                gene_tissue_membership[gene] = set()
            gene_tissue_membership[gene].add(tissue)

# Intersection size distribution
from collections import Counter
intersection_sizes = Counter()
for gene, tissues_set in gene_tissue_membership.items():
    k = len(tissues_set)
    intersection_sizes[k] += 1

upset_bars = [{"n_tissues": k, "n_genes": intersection_sizes[k]} for k in sorted(intersection_sizes.keys())]

# Top intersection sets (genes in 2+ tissues)
multi_tissue_genes = []
for gene, tissues_set in gene_tissue_membership.items():
    if len(tissues_set) >= 2:
        multi_tissue_genes.append({
            "gene": gene,
            "n_tissues": len(tissues_set),
            "tissues": sorted(tissues_set)
        })
multi_tissue_genes.sort(key=lambda x: (-x["n_tissues"], x["gene"]))

# --- Panel B: Gene × tissue heatmap ---
# Select genes: 5 consensus + top tissue-specific (up to 5 per tissue)
heatmap_genes = []
for g in consensus["consensus_genes"]:
    heatmap_genes.append({"gene": g["gene"], "type": "consensus", "n_tissues": g["n_tissues"]})

for tissue in TISSUES:
    ts_genes = consensus["tissue_specific_genes"].get(tissue, [])
    for g in ts_genes[:3]:  # Top 3 per tissue
        if g["gene"] not in [hg["gene"] for hg in heatmap_genes]:
            heatmap_genes.append({"gene": g["gene"], "type": f"specific_{tissue}", "n_tissues": 1})

# Build heatmap matrix
heatmap_data = []
for ginfo in heatmap_genes:
    gene_upper = ginfo["gene"].upper()
    for tissue in TISSUES:
        imp = per_tissue_shap.get(tissue, {}).get(gene_upper, 0)
        in_top100 = gene_upper in per_tissue_top100.get(tissue, set())
        heatmap_data.append({
            "gene": ginfo["gene"],
            "tissue": TISSUE_LABELS[tissue],
            "importance": round(imp, 6),
            "in_top100": in_top100,
            "gene_type": ginfo["type"]
        })

# --- Panel C: Method agreement ---
agreement_data = []
for tissue in TISSUES:
    if tissue in consensus["method_agreement"]:
        ma = consensus["method_agreement"][tissue]
        agreement_data.append({
            "tissue": TISSUE_LABELS[tissue],
            "method_a": ma["method_a"],
            "method_b": ma["method_b"],
            "jaccard": round(ma["jaccard_top100"], 3),
            "overlap": ma["n_overlap_top100"],
            "spearman": round(ma["spearman_rho"], 3)
        })

# --- Panel D: Cell 2020 validation ---
validation_data = []
for key, val in consensus["validation"].items():
    if "cell2020" in key:
        parts = key.split("_")
        # tissue_method_cell2020
        tissue = parts[0]
        method = "_".join(parts[1:-1])
        validation_data.append({
            "tissue": TISSUE_LABELS.get(tissue, tissue),
            "method": method,
            "n_overlap": val.get("n_overlap", 0),
            "n_reference": val.get("n_reference", 0),
            "p_value": round(val.get("p_value", 1.0), 4)
        })

# SHAP direction summary (up/down counts)
direction_data = []
for key, sd in consensus["shap_directions"].items():
    direction_data.append({
        "key": key,
        "tissue": sd["tissue"],
        "method": sd.get("method", key.split("_", 1)[-1] if "_" in key else ""),
        "n_up": sd["n_positive_in_top100"],
        "n_down": sd["n_negative_in_top100"],
        "top_up": [g[0] if isinstance(g, list) else g for g in sd.get("top10_upregulated", [])[:5]],
        "top_down": [g[0] if isinstance(g, list) else g for g in sd.get("top10_downregulated", [])[:5]]
    })

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig 3 — Consensus Spaceflight Gene Signature</title>
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
  table.data-table {{ border-collapse: collapse; width: 100%; font-size: 9px; margin-top: 8px; }}
  table.data-table th, table.data-table td {{ border: 1px solid #ddd; padding: 3px 6px; text-align: left; }}
  table.data-table th {{ background: #f0f0f0; font-weight: bold; }}
  .sig {{ color: #D55E00; font-weight: bold; }}
</style>
</head>
<body>
<div class="figure-container">
  <div class="figure-title">Fig 3 — Consensus Spaceflight Gene Signature from SHAP Interpretability</div>
  <div class="figure-subtitle">SHAP top-100 genes across 8 tissues × 2 methods per tissue (16 SHAP analyses)</div>
  <div class="panels">
    <div class="panel"><div class="panel-label">A</div>
      <div class="panel-desc">SHAP Top-100 Gene Overlap — Intersection Size Distribution</div>
      <div id="upset"></div>
    </div>
    <div class="panel"><div class="panel-label">B</div>
      <div class="panel-desc">Gene × Tissue SHAP Importance Heatmap (Consensus + Tissue-Specific)</div>
      <div id="gene_heatmap"></div>
    </div>
    <div class="panel"><div class="panel-label">C</div>
      <div class="panel-desc">Cross-Method Agreement — Jaccard & Spearman per Tissue</div>
      <div id="agreement"></div>
    </div>
    <div class="panel"><div class="panel-label">D</div>
      <div class="panel-desc">SHAP Direction Summary & Cell 2020 Validation</div>
      <div id="validation"></div>
    </div>
  </div>
  <div class="interpretation">
    <strong>Key findings:</strong>
    5 consensus genes appear in ≥3 tissues and ≥2 methods: <b>MUP22</b> (4 tissues), <b>NPAS2</b> (3 tissues, circadian),
    <b>PER2</b> (3 tissues, circadian), <b>GM24497</b> (3 tissues), <b>THRSP</b> (3 tissues, lipid metabolism).
    Cross-method agreement varies: skin/lung show highest Jaccard overlap (0.25), while thymus is lowest (0.005).
    Cell 2020 validation shows limited overlap (1-2/66 genes), reflecting that SHAP importance captures predictive
    features (not necessarily fold-change magnitude), and spaceflight signatures are tissue-specific.
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
  "Lung": "#888888", "Colon": "#555555"
}};
const tooltip = d3.select("#tooltip");

// =============================================
// PANEL A: UpSet-style intersection bar chart
// =============================================
(function() {{
  const upsetBars = {json.dumps(upset_bars)};
  const multiGenes = {json.dumps(multi_tissue_genes[:20])};
  const margin = {{top: 15, right: 20, bottom: 45, left: 55}};
  const width = 400, height = 180;

  const svg = d3.select("#upset").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom + 120)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(upsetBars.map(d => d.n_tissues)).range([0, width]).padding(0.3);
  const y = d3.scaleLinear().domain([0, d3.max(upsetBars, d => d.n_genes) * 1.1]).range([height, 0]);

  svg.selectAll("rect.bar")
    .data(upsetBars).enter().append("rect")
    .attr("x", d => x(d.n_tissues))
    .attr("y", d => y(d.n_genes))
    .attr("width", x.bandwidth())
    .attr("height", d => height - y(d.n_genes))
    .attr("fill", (d, i) => d.n_tissues >= 3 ? OI.vermilion : d.n_tissues === 2 ? OI.orange : OI.skyblue)
    .attr("opacity", 0.8);

  svg.selectAll("text.bar-label")
    .data(upsetBars).enter().append("text")
    .attr("x", d => x(d.n_tissues) + x.bandwidth()/2)
    .attr("y", d => y(d.n_genes) - 4)
    .attr("text-anchor", "middle").attr("font-size", "9px").attr("font-weight", "bold")
    .text(d => d.n_genes);

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x));
  svg.append("text").attr("x", width/2).attr("y", height+35)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("Number of Tissues");

  svg.append("g").call(d3.axisLeft(y).ticks(5));
  svg.append("text").attr("x", -height/2).attr("y", -40)
    .attr("text-anchor", "middle").attr("font-size", "10px")
    .attr("transform", "rotate(-90)").text("Number of Genes");

  // Multi-tissue gene table below
  const tableY = height + 55;
  svg.append("text").attr("x", 0).attr("y", tableY)
    .attr("font-size", "9px").attr("font-weight", "bold").text("Multi-tissue genes (≥2 tissues):");

  multiGenes.slice(0, 12).forEach((g, i) => {{
    svg.append("text").attr("x", 0).attr("y", tableY + 14 + i * 12)
      .attr("font-size", "8px")
      .attr("fill", g.n_tissues >= 3 ? OI.vermilion : "#444")
      .text(`${{g.gene}} (${{g.n_tissues}}): ${{g.tissues.join(", ")}}`);
  }});
}})();

// =============================================
// PANEL B: Gene × Tissue SHAP Heatmap
// =============================================
(function() {{
  const data = {json.dumps(heatmap_data)};
  const margin = {{top: 5, right: 15, bottom: 60, left: 120}};
  const genes = [...new Set(data.map(d => d.gene))];
  const tissues = [...new Set(data.map(d => d.tissue))];
  const cellW = 46, cellH = 16;
  const width = tissues.length * cellW;
  const height = genes.length * cellH;

  const svg = d3.select("#gene_heatmap").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(tissues).range([0, width]).padding(0.05);
  const y = d3.scaleBand().domain(genes).range([0, height]).padding(0.05);

  // Max importance for color scale
  const maxImp = d3.max(data.filter(d => d.importance > 0), d => d.importance);
  const color = d3.scaleSequential().domain([0, maxImp * 0.5]).interpolator(d3.interpolateYlOrRd);

  svg.selectAll("rect.cell")
    .data(data).enter().append("rect")
    .attr("x", d => x(d.tissue))
    .attr("y", d => y(d.gene))
    .attr("width", x.bandwidth())
    .attr("height", y.bandwidth())
    .attr("fill", d => d.importance > 0 ? color(d.importance) : "#f5f5f5")
    .attr("stroke", d => d.in_top100 ? "#333" : "#eee")
    .attr("stroke-width", d => d.in_top100 ? 1.5 : 0.5)
    .attr("rx", 1)
    .on("mouseover", function(event, d) {{
      tooltip.style("display", "block")
        .html(`<b>${{d.gene}}</b> × ${{d.tissue}}<br>Importance: ${{d.importance.toExponential(2)}}<br>In top-100: ${{d.in_top100 ? "Yes" : "No"}}<br>Type: ${{d.gene_type}}`);
    }})
    .on("mousemove", function(event) {{
      tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
    }})
    .on("mouseout", () => tooltip.style("display", "none"));

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).tickSize(0))
    .selectAll("text").attr("transform", "rotate(-40)")
    .style("text-anchor", "end").attr("dx", "-0.5em").style("font-size", "8px");

  svg.append("g").call(d3.axisLeft(y).tickSize(0))
    .selectAll("text").style("font-size", "7px")
    .attr("fill", function(d) {{
      const ginfo = data.find(dd => dd.gene === d);
      return ginfo && ginfo.gene_type === "consensus" ? OI.vermilion : "#333";
    }});

  // Note
  svg.append("text").attr("x", 0).attr("y", height + 50)
    .attr("font-size", "7px").attr("fill", "#888")
    .text("Bold border = SHAP top-100. Red labels = consensus genes.");
}})();

// =============================================
// PANEL C: Method Agreement
// =============================================
(function() {{
  const data = {json.dumps(agreement_data)};
  const margin = {{top: 15, right: 30, bottom: 45, left: 100}};
  const width = 380, height = 250;

  const svg = d3.select("#agreement").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const y = d3.scaleBand().domain(data.map(d => d.tissue)).range([0, height]).padding(0.25);
  const x = d3.scaleLinear().domain([0, 0.45]).range([0, width]);

  // Jaccard bars
  svg.selectAll("rect.jaccard")
    .data(data).enter().append("rect")
    .attr("x", 0).attr("y", d => y(d.tissue))
    .attr("width", d => x(d.jaccard))
    .attr("height", y.bandwidth() / 2)
    .attr("fill", OI.blue).attr("opacity", 0.8);

  // Spearman bars
  svg.selectAll("rect.spearman")
    .data(data).enter().append("rect")
    .attr("x", 0).attr("y", d => y(d.tissue) + y.bandwidth() / 2)
    .attr("width", d => x(Math.max(0, d.spearman)))
    .attr("height", y.bandwidth() / 2)
    .attr("fill", OI.orange).attr("opacity", 0.8);

  // Value labels
  data.forEach(d => {{
    svg.append("text")
      .attr("x", x(d.jaccard) + 3)
      .attr("y", y(d.tissue) + y.bandwidth() / 4 + 3)
      .attr("font-size", "8px").attr("fill", OI.blue)
      .text(`J=${{d.jaccard}} (n=${{d.overlap}})`);
    svg.append("text")
      .attr("x", x(Math.max(0, d.spearman)) + 3)
      .attr("y", y(d.tissue) + y.bandwidth() * 3/4 + 3)
      .attr("font-size", "8px").attr("fill", OI.orange)
      .text(`ρ=${{d.spearman}}`);
  }});

  svg.append("g").call(d3.axisLeft(y).tickSize(0))
    .selectAll("text").style("font-size", "9px");
  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).ticks(5));
  svg.append("text").attr("x", width/2).attr("y", height+35)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("Agreement Score");

  // Legend
  svg.append("rect").attr("x", width-120).attr("y", 5).attr("width", 12).attr("height", 8).attr("fill", OI.blue).attr("opacity", 0.8);
  svg.append("text").attr("x", width-105).attr("y", 13).attr("font-size", "8px").text("Jaccard (top-100)");
  svg.append("rect").attr("x", width-120).attr("y", 18).attr("width", 12).attr("height", 8).attr("fill", OI.orange).attr("opacity", 0.8);
  svg.append("text").attr("x", width-105).attr("y", 26).attr("font-size", "8px").text("Spearman ρ (rank)");
}})();

// =============================================
// PANEL D: Validation & Directions
// =============================================
(function() {{
  const validData = {json.dumps(validation_data)};
  const dirData = {json.dumps(direction_data)};

  const container = d3.select("#validation");

  // Direction summary table
  container.append("div").attr("style", "font-size:9px; font-weight:bold; margin-bottom:4px;")
    .text("SHAP Direction (FLT ↑ vs GC ↑):");

  let html = '<table class="data-table"><tr><th>Tissue</th><th>Method</th><th>↑FLT</th><th>↓FLT</th><th>Top Up</th></tr>';
  dirData.forEach(d => {{
    html += `<tr><td>${{d.tissue}}</td><td>${{d.method}}</td><td>${{d.n_up}}</td><td>${{d.n_down}}</td><td style="font-size:8px">${{d.top_up.slice(0,3).join(", ")}}</td></tr>`;
  }});
  html += '</table>';
  container.append("div").html(html);

  // Cell 2020 validation table
  container.append("div").attr("style", "font-size:9px; font-weight:bold; margin-top:12px; margin-bottom:4px;")
    .text("da Silveira Cell 2020 Validation (67 genes):");

  // Filter to show only entries with overlap > 0
  const sigValid = validData.filter(d => d.n_overlap > 0);
  let vhtml = '<table class="data-table"><tr><th>Tissue</th><th>Method</th><th>Overlap</th><th>p (hypergeom)</th></tr>';
  if (sigValid.length > 0) {{
    sigValid.forEach(d => {{
      const cls = d.p_value < 0.05 ? ' class="sig"' : '';
      vhtml += `<tr><td>${{d.tissue}}</td><td>${{d.method}}</td><td>${{d.n_overlap}}/${{d.n_reference}}</td><td${{cls}}>${{d.p_value.toFixed(4)}}</td></tr>`;
    }});
  }} else {{
    vhtml += '<tr><td colspan="4">No overlaps found across all tissue×method combinations</td></tr>';
  }}
  vhtml += '</table>';
  container.append("div").html(vhtml);

  container.append("div").attr("style", "font-size:8px; color:#888; margin-top:8px;")
    .text("Note: Limited Cell 2020 overlap expected — da Silveira genes are mostly from liver/muscle DEGs, while SHAP captures predictive importance (not fold-change magnitude).");
}})();
</script>
</body>
</html>"""

out_path = os.path.join(OUT_DIR, "Fig3_consensus.html")
with open(out_path, "w") as f:
    f.write(html)
print(f"Written: {out_path}")
print(f"  Panel A: {len(upset_bars)} bars, {len(multi_tissue_genes)} multi-tissue genes")
print(f"  Panel B: {len(heatmap_genes)} genes × 8 tissues")
print(f"  Panel C: {len(agreement_data)} tissues")
print(f"  Panel D: {len(direction_data)} SHAP directions, {len(validation_data)} validations")

#!/usr/bin/env python3
"""
Fig 4 — Network Biology (WGCNA + SHAP Integration + PPI)
GeneLabBench v4 Publication Figure

Panels:
  A: Module-trait correlation heatmap (6 tissues × significant modules)
  B: Preservation Zsummary lower-triangle heatmap (6×6 tissues)
  C: SHAP→WGCNA enrichment (dot plot) + eye crystallin module highlight
  D: PPI enrichment summary (top terms per tissue)

Output: v4/figures/html/Fig4_network.html (self-contained D3.js v7)
"""

import json
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EVAL_DIR = os.path.join(BASE, "evaluation")
OUT_DIR = os.path.join(BASE, "figures", "html")
os.makedirs(OUT_DIR, exist_ok=True)

WGCNA_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin"]
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastro", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin"
}

# --- Load WGCNA module data ---
wgcna_modules = {}
module_trait_data = []
for tissue in WGCNA_TISSUES:
    fpath = os.path.join(EVAL_DIR, f"WGCNA_{tissue}_modules.json")
    with open(fpath) as f:
        d = json.load(f)
    wgcna_modules[tissue] = d
    # condition_modules has modules with significant trait correlations
    for cm in d.get("condition_modules", []):
        module_trait_data.append({
            "tissue": TISSUE_LABELS[tissue],
            "module": cm["module"],
            "r": round(cm["r"], 3),
            "p": cm["p"],
            "direction": cm["direction"],
            "n_genes": d["modules"].get(cm["module"], {}).get("n_genes", 0)
        })
    # Also add non-significant modules for context (all modules with their trait corr)
    for mod_color, mod_info in d["modules"].items():
        # Check if already added
        existing = [mt for mt in module_trait_data
                    if mt["tissue"] == TISSUE_LABELS[tissue] and mt["module"] == mod_color]
        if not existing:
            # Get trait correlation from significant_trait_associations
            sta = mod_info.get("significant_trait_associations", [])
            r_val = sta[0]["r"] if sta else 0
            p_val = sta[0]["p"] if sta else 1.0
            module_trait_data.append({
                "tissue": TISSUE_LABELS[tissue],
                "module": mod_color,
                "r": round(r_val, 3),
                "p": round(p_val, 4),
                "direction": "FLT_higher" if r_val > 0 else "GC_higher" if r_val < 0 else "none",
                "n_genes": mod_info.get("n_genes", 0)
            })

# --- Load preservation data ---
with open(os.path.join(EVAL_DIR, "WGCNA_preservation.json")) as f:
    pres = json.load(f)

pres_matrix = []
for ps in pres["pair_summary"]:
    pair = ps["pair"]
    t1, t2 = pair.split("_vs_")
    z = ps["mean_Zsummary"]
    pres_matrix.append({
        "ref": TISSUE_LABELS.get(t1, t1),
        "test": TISSUE_LABELS.get(t2, t2),
        "Zsummary": round(z, 2) if z is not None else 0,
        "n_high": ps["n_high_pres"],
        "n_moderate": ps["n_moderate_pres"],
        "n_not": ps["n_not_pres"]
    })

# --- Load SHAP-WGCNA integration ---
with open(os.path.join(EVAL_DIR, "SHAP_WGCNA_integration.json")) as f:
    swi = json.load(f)

shap_wgcna_data = []
for tissue in swi.get("tissues", {}):
    for em in swi["tissues"][tissue].get("enriched_modules", []):
        shap_wgcna_data.append({
            "tissue": TISSUE_LABELS.get(tissue, tissue),
            "module": em["module"],
            "n_shap": em["n_shap_in_module"],
            "n_module": em["n_module_genes"],
            "p_value": em["p_value"],
            "ratio": round(em["enrichment_ratio"], 2),
            "genes": em["overlap_genes"][:8]
        })

# --- Load PPI enrichment ---
with open(os.path.join(EVAL_DIR, "PPI_enrichment.json")) as f:
    ppi = json.load(f)

ppi_data = []
for tissue_name in ppi.get("tissues", {}):
    t = ppi["tissues"][tissue_name]
    density = t.get("network", {}).get("density", 0) if isinstance(t.get("network"), dict) else 0
    terms = t.get("enrichment", {}).get("terms", []) if isinstance(t.get("enrichment"), dict) else []
    if not terms:
        terms = t.get("enriched_terms", [])  # fallback field name
    n_sig = t.get("n_significant_terms", len(terms))
    top_terms = []
    for term in terms[:5]:
        top_terms.append({
            "desc": term.get("description", ""),
            "p_fdr": term.get("p_fdr", 1.0),
            "genes": term.get("genes_in_term", [])[:5]
        })
    ppi_data.append({
        "tissue": TISSUE_LABELS.get(tissue_name, tissue_name),
        "density": round(density, 3) if density else 0,
        "n_terms": n_sig,
        "top_terms": top_terms
    })

# Also add tissues with separate data (lung, colon, skin may be in separate structure)
ALL_PPI_TISSUES = [t["tissue"] for t in ppi_data]
for tissue in ["Skin", "Lung", "Colon"]:
    t_key = [k for k in ppi.get("tissues", {}) if TISSUE_LABELS.get(k) == tissue or k == tissue.lower()]
    if t_key and tissue not in ALL_PPI_TISSUES:
        # Already handled above
        pass

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig 4 — Network Biology</title>
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
  table.data-table {{ border-collapse: collapse; width: 100%; font-size: 9px; margin-top: 4px; }}
  table.data-table th, table.data-table td {{ border: 1px solid #ddd; padding: 3px 5px; text-align: left; }}
  table.data-table th {{ background: #f0f0f0; font-weight: bold; }}
  .sig {{ color: #D55E00; font-weight: bold; }}
</style>
</head>
<body>
<div class="figure-container">
  <div class="figure-title">Fig 4 — Network Biology: WGCNA Co-Expression and PPI Enrichment</div>
  <div class="figure-subtitle">6 tissues × WGCNA modules + SHAP integration + STRING PPI (v12)</div>
  <div class="panels">
    <div class="panel"><div class="panel-label">A</div>
      <div class="panel-desc">Module-Trait Correlation Heatmap (6 tissues, condition: FLT vs GC)</div>
      <div id="trait_heatmap"></div>
    </div>
    <div class="panel"><div class="panel-label">B</div>
      <div class="panel-desc">Module Preservation — Mean Zsummary (15 tissue pairs)</div>
      <div id="preservation"></div>
    </div>
    <div class="panel"><div class="panel-label">C</div>
      <div class="panel-desc">SHAP → WGCNA Module Enrichment</div>
      <div id="shap_wgcna"></div>
    </div>
    <div class="panel"><div class="panel-label">D</div>
      <div class="panel-desc">STRING PPI Enrichment Summary</div>
      <div id="ppi_summary"></div>
    </div>
  </div>
  <div class="interpretation">
    <strong>Key findings:</strong>
    Gastrocnemius saddlebrown module (r=+0.496, FLT↑) and thymus black module (r=−0.479, GC↑) show
    the strongest condition associations. Liver→kidney preservation is strongest (Z=35.65), suggesting
    shared core transcriptomic programs. SHAP-enriched modules include eye green (crystallins: CRYAA,
    CRYBA1, CRYBB2, CRYGB, CRYGC, CRYGS) and gastrocnemius saddlebrown (7.5× enrichment).
    STRING PPI: skin has highest network density (0.53) with 39 enriched terms.
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
// PANEL A: Module-Trait Correlation Heatmap
// =============================================
(function() {{
  const data = {json.dumps(module_trait_data)};
  const margin = {{top: 5, right: 60, bottom: 60, left: 80}};

  const tissues = [...new Set(data.map(d => d.tissue))];
  const modules = [...new Set(data.map(d => d.module))];
  const cellW = 26, cellH = 26;
  const width = modules.length * cellW;
  const height = tissues.length * cellH;

  const svg = d3.select("#trait_heatmap").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(modules).range([0, width]).padding(0.05);
  const y = d3.scaleBand().domain(tissues).range([0, height]).padding(0.05);

  const color = d3.scaleDiverging()
    .domain([-0.5, 0, 0.5])
    .interpolator(d3.interpolateRdBu);

  svg.selectAll("rect.cell")
    .data(data).enter().append("rect")
    .attr("x", d => x(d.module))
    .attr("y", d => y(d.tissue))
    .attr("width", x.bandwidth())
    .attr("height", y.bandwidth())
    .attr("fill", d => d.r !== 0 ? color(-d.r) : "#f5f5f5")
    .attr("stroke", d => d.p < 0.05 ? "#333" : "#eee")
    .attr("stroke-width", d => d.p < 0.05 ? 2 : 0.5)
    .attr("rx", 2)
    .on("mouseover", function(event, d) {{
      tooltip.style("display", "block")
        .html(`<b>${{d.tissue}}</b> — ${{d.module}}<br>r = ${{d.r}}<br>p = ${{d.p.toExponential(2)}}<br>n_genes = ${{d.n_genes}}<br>Direction: ${{d.direction}}`);
    }})
    .on("mousemove", function(event) {{
      tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
    }})
    .on("mouseout", () => tooltip.style("display", "none"));

  // Significance stars
  svg.selectAll("text.sig")
    .data(data.filter(d => d.p < 0.05)).enter().append("text")
    .attr("x", d => x(d.module) + x.bandwidth()/2)
    .attr("y", d => y(d.tissue) + y.bandwidth()/2 + 4)
    .attr("text-anchor", "middle").attr("font-size", "10px").attr("font-weight", "bold")
    .attr("fill", d => Math.abs(d.r) > 0.3 ? "#fff" : "#333")
    .text(d => d.p < 0.01 ? "**" : "*");

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).tickSize(0))
    .selectAll("text").attr("transform", "rotate(-55)")
    .style("text-anchor", "end").attr("dx", "-0.3em").style("font-size", "7px");

  svg.append("g").call(d3.axisLeft(y).tickSize(0))
    .selectAll("text").style("font-size", "9px");

  // Legend
  const lx = width + 10;
  [[-0.5, "GC ↑ (r < 0)"], [0, "Neutral"], [0.5, "FLT ↑ (r > 0)"]].forEach((item, i) => {{
    svg.append("rect").attr("x", lx).attr("y", i * 18)
      .attr("width", 12).attr("height", 12)
      .attr("fill", color(-item[0]));
    svg.append("text").attr("x", lx + 16).attr("y", i * 18 + 10)
      .attr("font-size", "7px").text(item[1]);
  }});
  svg.append("text").attr("x", lx).attr("y", 68)
    .attr("font-size", "7px").attr("fill", "#666").text("Bold border = p<0.05");
}})();

// =============================================
// PANEL B: Preservation Zsummary Heatmap
// =============================================
(function() {{
  const data = {json.dumps(pres_matrix)};
  const tissues = ["Liver", "Gastro", "Kidney", "Thymus", "Eye", "Skin"];
  const margin = {{top: 5, right: 60, bottom: 60, left: 70}};
  const cellS = 55;
  const width = tissues.length * cellS;
  const height = tissues.length * cellS;

  const svg = d3.select("#preservation").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(tissues).range([0, width]).padding(0.05);
  const y = d3.scaleBand().domain(tissues).range([0, height]).padding(0.05);

  const color = d3.scaleThreshold()
    .domain([2, 10])
    .range(["#e0e0e0", "#FFD700", "#2ca02c"]);

  // Build lookup
  const lookup = {{}};
  data.forEach(d => {{
    lookup[d.ref + "_" + d.test] = d;
  }});

  // Draw lower triangle
  for (let i = 0; i < tissues.length; i++) {{
    for (let j = 0; j < tissues.length; j++) {{
      if (i === j) continue;
      const key = tissues[i] + "_" + tissues[j];
      const entry = lookup[key];
      if (!entry) continue;

      svg.append("rect")
        .attr("x", x(tissues[j]))
        .attr("y", y(tissues[i]))
        .attr("width", x.bandwidth())
        .attr("height", y.bandwidth())
        .attr("fill", color(entry.Zsummary))
        .attr("stroke", "#fff").attr("stroke-width", 1)
        .attr("rx", 3)
        .on("mouseover", function(event) {{
          tooltip.style("display", "block")
            .html(`<b>${{entry.ref}} → ${{entry.test}}</b><br>Mean Z = ${{entry.Zsummary}}<br>High: ${{entry.n_high}}, Mod: ${{entry.n_moderate}}, Low: ${{entry.n_not}}`);
        }})
        .on("mousemove", function(event) {{
          tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
        }})
        .on("mouseout", () => tooltip.style("display", "none"));

      // Z value text
      svg.append("text")
        .attr("x", x(tissues[j]) + x.bandwidth()/2)
        .attr("y", y(tissues[i]) + y.bandwidth()/2 + 4)
        .attr("text-anchor", "middle").attr("font-size", "8px")
        .attr("font-weight", entry.Zsummary > 10 ? "bold" : "normal")
        .attr("fill", entry.Zsummary > 10 ? "#fff" : "#333")
        .text(entry.Zsummary.toFixed(1));
    }}
  }}

  // Diagonal
  tissues.forEach(t => {{
    svg.append("rect")
      .attr("x", x(t)).attr("y", y(t))
      .attr("width", x.bandwidth()).attr("height", y.bandwidth())
      .attr("fill", "#f0f0f0").attr("stroke", "#ddd").attr("rx", 3);
    svg.append("text")
      .attr("x", x(t) + x.bandwidth()/2)
      .attr("y", y(t) + y.bandwidth()/2 + 4)
      .attr("text-anchor", "middle").attr("font-size", "8px").attr("fill", "#999")
      .text("—");
  }});

  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).tickSize(0))
    .selectAll("text").style("font-size", "9px");
  svg.append("g").call(d3.axisLeft(y).tickSize(0))
    .selectAll("text").style("font-size", "9px");

  // Legend
  const lx = width + 10;
  [["#e0e0e0", "Z < 2 (not preserved)"], ["#FFD700", "2 ≤ Z < 10 (moderate)"], ["#2ca02c", "Z ≥ 10 (high)"]].forEach((item, i) => {{
    svg.append("rect").attr("x", lx).attr("y", i * 16).attr("width", 12).attr("height", 12).attr("fill", item[0]);
    svg.append("text").attr("x", lx + 16).attr("y", i * 16 + 10).attr("font-size", "7px").text(item[1]);
  }});
}})();

// =============================================
// PANEL C: SHAP → WGCNA Enrichment
// =============================================
(function() {{
  const data = {json.dumps(shap_wgcna_data)};
  const container = d3.select("#shap_wgcna");

  if (data.length === 0) {{
    container.append("div").attr("style", "font-size:10px; color:#888; padding:20px;")
      .text("No significant SHAP→WGCNA module enrichments found.");
    return;
  }}

  const margin = {{top: 15, right: 20, bottom: 40, left: 120}};
  const barH = 28;
  const width = 380;
  const height = data.length * barH;

  const svg = container.append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom + 100)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleLinear()
    .domain([0, d3.max(data, d => d.ratio) * 1.2])
    .range([0, width]);
  const y = d3.scaleBand()
    .domain(data.map(d => d.tissue + " " + d.module))
    .range([0, height]).padding(0.2);

  svg.selectAll("rect.bar")
    .data(data).enter().append("rect")
    .attr("x", 0)
    .attr("y", d => y(d.tissue + " " + d.module))
    .attr("width", d => x(d.ratio))
    .attr("height", y.bandwidth())
    .attr("fill", d => d.p_value < 0.001 ? OI.vermilion : OI.orange)
    .attr("opacity", 0.8)
    .on("mouseover", function(event, d) {{
      tooltip.style("display", "block")
        .html(`<b>${{d.tissue}} — ${{d.module}}</b><br>Enrichment: ${{d.ratio}}×<br>SHAP genes in module: ${{d.n_shap}}/${{d.n_module}}<br>p = ${{d.p_value.toExponential(2)}}<br>Genes: ${{d.genes.join(", ")}}`);
    }})
    .on("mousemove", function(event) {{
      tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
    }})
    .on("mouseout", () => tooltip.style("display", "none"));

  // Labels
  data.forEach(d => {{
    const label = d.tissue + " " + d.module;
    svg.append("text")
      .attr("x", x(d.ratio) + 4)
      .attr("y", y(label) + y.bandwidth()/2 + 4)
      .attr("font-size", "8px").attr("fill", "#333")
      .text(`${{d.ratio}}× (${{d.n_shap}}/${{d.n_module}})`);
  }});

  svg.append("g").call(d3.axisLeft(y).tickSize(0))
    .selectAll("text").style("font-size", "8px");
  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).ticks(5));
  svg.append("text").attr("x", width/2).attr("y", height + 32)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("Enrichment Ratio");

  // Eye crystallin highlight
  const eyeEntry = data.find(d => d.tissue === "Eye" && d.module === "green");
  if (eyeEntry) {{
    const cy = height + 55;
    svg.append("text").attr("x", 0).attr("y", cy)
      .attr("font-size", "9px").attr("font-weight", "bold").attr("fill", OI.pink)
      .text("Eye green module — Crystallin genes:");
    svg.append("text").attr("x", 0).attr("y", cy + 14)
      .attr("font-size", "8px").attr("fill", "#444")
      .text(eyeEntry.genes.join(", "));
  }}
}})();

// =============================================
// PANEL D: PPI Enrichment Summary
// =============================================
(function() {{
  const data = {json.dumps(ppi_data)};
  const container = d3.select("#ppi_summary");

  // Bar chart: n_terms per tissue
  const margin = {{top: 15, right: 20, bottom: 45, left: 70}};
  const width = 300, height = 200;

  const svg = container.append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const sortedData = data.filter(d => d.n_terms > 0).sort((a, b) => b.n_terms - a.n_terms);
  const y = d3.scaleBand().domain(sortedData.map(d => d.tissue)).range([0, height]).padding(0.3);
  const x = d3.scaleLinear().domain([0, d3.max(sortedData, d => d.n_terms) * 1.1]).range([0, width]);

  svg.selectAll("rect.bar")
    .data(sortedData).enter().append("rect")
    .attr("x", 0).attr("y", d => y(d.tissue))
    .attr("width", d => x(d.n_terms))
    .attr("height", y.bandwidth())
    .attr("fill", OI.blue).attr("opacity", 0.8);

  sortedData.forEach(d => {{
    svg.append("text")
      .attr("x", x(d.n_terms) + 4)
      .attr("y", y(d.tissue) + y.bandwidth()/2 + 4)
      .attr("font-size", "8px").attr("fill", "#333")
      .text(`${{d.n_terms}} terms (d=${{d.density}})`);
  }});

  svg.append("g").call(d3.axisLeft(y).tickSize(0))
    .selectAll("text").style("font-size", "9px");
  svg.append("g").attr("transform", `translate(0,${{height}})`)
    .call(d3.axisBottom(x).ticks(5));
  svg.append("text").attr("x", width/2).attr("y", height + 35)
    .attr("text-anchor", "middle").attr("font-size", "10px").text("Enriched Terms (FDR < 0.05)");

  // Top terms table
  let html = '<div style="margin-top:12px; font-size:9px; font-weight:bold;">Top Enriched GO/KEGG Terms:</div>';
  html += '<table class="data-table"><tr><th>Tissue</th><th>Term</th><th>p_FDR</th><th>Genes</th></tr>';
  sortedData.forEach(d => {{
    d.top_terms.slice(0, 2).forEach(t => {{
      html += `<tr><td>${{d.tissue}}</td><td>${{t.desc}}</td><td>${{t.p_fdr < 0.001 ? t.p_fdr.toExponential(1) : t.p_fdr.toFixed(3)}}</td><td style="font-size:8px">${{t.genes.slice(0,4).join(", ")}}</td></tr>`;
    }});
  }});
  html += '</table>';
  container.append("div").html(html);
}})();
</script>
</body>
</html>"""

out_path = os.path.join(OUT_DIR, "Fig4_network.html")
with open(out_path, "w") as f:
    f.write(html)
print(f"Written: {out_path}")
print(f"  Panel A: {len(module_trait_data)} module-trait entries across 6 tissues")
print(f"  Panel B: {len(pres_matrix)} preservation pairs")
print(f"  Panel C: {len(shap_wgcna_data)} enriched modules")
print(f"  Panel D: {len(ppi_data)} tissues with PPI data")

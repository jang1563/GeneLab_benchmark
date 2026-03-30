#!/usr/bin/env python3
"""
Fig S3 — WGCNA Soft Threshold Selection & Module Overview
6 tissues: soft threshold R² vs power curves + module size bar charts

Output: v4/figures/html/FigS3_wgcna_overview.html
"""

import csv
import json
import os

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WGCNA_DIR = os.path.join(BASE, "wgcna_outputs")
EVAL_DIR = os.path.join(BASE, "evaluation")
OUT_DIR = os.path.join(BASE, "figures", "html")
os.makedirs(OUT_DIR, exist_ok=True)

TISSUES = ["liver", "kidney", "thymus", "skin", "eye", "gastrocnemius"]
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastrocnemius", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin"
}

# Load data per tissue
tissue_data = {}
for tissue in TISSUES:
    tdir = os.path.join(WGCNA_DIR, tissue)

    # Summary
    with open(os.path.join(tdir, "summary.json")) as f:
        summary = json.load(f)

    # Soft threshold
    sft_rows = []
    with open(os.path.join(tdir, "soft_threshold.csv")) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sft_rows.append({
                "power": int(row["Power"]),
                "r2": round(float(row["SFT.R.sq"]), 4),
                "slope": round(float(row["slope"]), 4),
                "mean_k": round(float(row["mean.k"]), 2)
            })

    # Module sizes
    mod_rows = []
    with open(os.path.join(tdir, "module_sizes.csv")) as f:
        reader = csv.DictReader(f)
        for row in reader:
            mod_rows.append({
                "color": row["module_color"],
                "n_genes": int(row["n_genes"])
            })

    # Trait correlations from WGCNA_*_modules.json
    modfile = os.path.join(EVAL_DIR, f"WGCNA_{tissue}_modules.json")
    cond_mods = []
    if os.path.exists(modfile):
        with open(modfile) as f:
            mdata = json.load(f)
        cond_mods = mdata.get("condition_modules", []) or []

    tissue_data[tissue] = {
        "label": TISSUE_LABELS[tissue],
        "beta": summary["soft_threshold_beta"],
        "r2": round(summary["r2_achieved"], 4),
        "n_samples": summary["n_samples"],
        "n_genes": summary["n_genes_input"],
        "n_modules": summary["n_modules_final"],
        "n_grey": summary["n_grey_genes"],
        "sft": sft_rows,
        "modules": mod_rows,
        "condition_modules": cond_mods
    }

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig S3 — WGCNA Soft Threshold & Module Overview</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ font-family: Arial, Helvetica, sans-serif; margin: 0; padding: 20px; background: #fff; }}
  .figure-container {{ width: 1200px; margin: 0 auto; }}
  .figure-title {{ font-size: 14px; font-weight: bold; text-align: center; margin-bottom: 4px; }}
  .figure-subtitle {{ font-size: 10px; color: #666; text-align: center; margin-bottom: 16px; }}
  .tissue-row {{ display: flex; gap: 12px; margin-bottom: 16px; border: 1px solid #e0e0e0; border-radius: 4px; padding: 10px; background: #fafafa; }}
  .panel {{ flex: 1; }}
  .panel-label {{ font-size: 10px; font-weight: bold; margin-bottom: 4px; }}
  .tissue-info {{ font-size: 8px; color: #666; margin-bottom: 4px; }}
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
  <div class="figure-title">Fig S3 — WGCNA Network Construction: Soft Threshold Selection & Module Summary</div>
  <div class="figure-subtitle">6 tissues | Signed hybrid network, R&sup2; &ge; 0.80 threshold, minModuleSize=30</div>
  <div id="panels"></div>
  <div class="note">
    Left: Scale-free topology fit (R&sup2;) vs soft threshold power. Dashed line = R&sup2;=0.80 cutoff. Red dot = selected &beta;.
    Right: Module gene counts. Color = WGCNA module color. Border = significant condition association (p&lt;0.05).
  </div>
</div>
<div class="tooltip" id="tooltip" style="display:none;"></div>

<script>
const tooltip = d3.select("#tooltip");
const data = {json.dumps(tissue_data)};
const tissueOrder = {json.dumps(TISSUES)};

// Module color map for common WGCNA colors
const wgcnaColorMap = {{
  "turquoise": "#40E0D0", "blue": "#4169E1", "brown": "#8B4513",
  "yellow": "#FFD700", "green": "#228B22", "red": "#DC143C",
  "black": "#333333", "pink": "#FF69B4", "magenta": "#FF00FF",
  "purple": "#800080", "greenyellow": "#ADFF2F", "tan": "#D2B48C",
  "salmon": "#FA8072", "cyan": "#00CED1", "midnightblue": "#191970",
  "lightcyan": "#E0FFFF", "grey": "#C0C0C0", "lightgreen": "#90EE90",
  "lightyellow": "#FFFFE0", "royalblue": "#4169E1", "darkred": "#8B0000",
  "darkgreen": "#006400", "darkturquoise": "#00CED1", "darkgrey": "#A9A9A9",
  "orange": "#FF8C00", "darkorange": "#FF8C00", "white": "#F0F0F0",
  "saddlebrown": "#8B4513", "steelblue": "#4682B4", "paleturquoise": "#AFEEEE",
  "violet": "#EE82EE", "darkolivegreen": "#556B2F", "orangered": "#FF4500",
  "darkmagenta": "#8B008B", "sienna": "#A0522D", "skyblue": "#87CEEB",
  "mediumpurple": "#9370DB", "plum": "#DDA0DD", "yellowgreen": "#9ACD32",
  "lightsteelblue": "#B0C4DE", "floralwhite": "#FFFAF0", "honeydew": "#F0FFF0",
  "ivory": "#FFFFF0", "navajowhite": "#FFDEAD"
}};

tissueOrder.forEach((tissue, idx) => {{
  const td = data[tissue];
  const row = d3.select("#panels").append("div").attr("class", "tissue-row");

  // --- LEFT: Soft threshold R² curve ---
  const leftPanel = row.append("div").attr("class", "panel");
  leftPanel.append("div").attr("class", "panel-label")
    .text(`${{String.fromCharCode(65 + idx)}}. ${{td.label}} (n=${{td.n_samples}}, ${{td.n_genes}} genes)`);
  leftPanel.append("div").attr("class", "tissue-info")
    .text(`Selected \\u03B2=${{td.beta}}, R\\u00B2=${{td.r2}}, ${{td.n_modules}} modules`);

  const margin = {{top: 5, right: 30, bottom: 30, left: 35}};
  const w = 260, h = 120;

  const svg1 = leftPanel.append("svg")
    .attr("width", w + margin.left + margin.right)
    .attr("height", h + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const sft = td.sft;
  const xScale = d3.scaleLinear().domain([0, 30]).range([0, w]);
  const yScale = d3.scaleLinear().domain([0, 1]).range([h, 0]);

  // R²=0.80 threshold line
  svg1.append("line")
    .attr("x1", 0).attr("x2", w)
    .attr("y1", yScale(0.80)).attr("y2", yScale(0.80))
    .attr("stroke", "#E69F00").attr("stroke-dasharray", "4,3").attr("stroke-width", 1);

  // R² curve
  const line = d3.line().x(d => xScale(d.power)).y(d => yScale(d.r2));
  svg1.append("path").datum(sft)
    .attr("d", line).attr("fill", "none").attr("stroke", "#0072B2").attr("stroke-width", 1.5);

  // Points
  sft.forEach(d => {{
    svg1.append("circle")
      .attr("cx", xScale(d.power)).attr("cy", yScale(d.r2))
      .attr("r", d.power === td.beta ? 5 : 2.5)
      .attr("fill", d.power === td.beta ? "#D55E00" : "#0072B2")
      .attr("stroke", d.power === td.beta ? "#D55E00" : "none")
      .attr("stroke-width", d.power === td.beta ? 2 : 0)
      .on("mouseover", function(event) {{
        tooltip.style("display", "block")
          .html(`<b>\\u03B2=${{d.power}}</b><br>R\\u00B2=${{d.r2}}<br>slope=${{d.slope}}<br>mean k=${{d.mean_k}}`);
      }})
      .on("mousemove", function(event) {{
        tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
      }})
      .on("mouseout", () => tooltip.style("display", "none"));
  }});

  // Selected beta label
  const selPt = sft.find(d => d.power === td.beta);
  if (selPt) {{
    svg1.append("text")
      .attr("x", xScale(selPt.power) + 6).attr("y", yScale(selPt.r2) - 6)
      .attr("font-size", "7px").attr("fill", "#D55E00").attr("font-weight", "bold")
      .text(`\\u03B2=${{td.beta}}`);
  }}

  svg1.append("g").attr("transform", `translate(0,${{h}})`)
    .call(d3.axisBottom(xScale).ticks(6).tickSize(3))
    .selectAll("text").style("font-size", "7px");
  svg1.append("g").call(d3.axisLeft(yScale).ticks(5).tickFormat(d3.format(".1f")).tickSize(3))
    .selectAll("text").style("font-size", "7px");
  svg1.append("text").attr("x", w/2).attr("y", h + 22).attr("text-anchor", "middle")
    .attr("font-size", "8px").text("Soft threshold power (\\u03B2)");
  svg1.append("text").attr("transform", "rotate(-90)").attr("x", -h/2).attr("y", -25)
    .attr("text-anchor", "middle").attr("font-size", "8px").text("Scale-free R\\u00B2");

  // --- RIGHT: Module sizes bar chart ---
  const rightPanel = row.append("div").attr("class", "panel");
  rightPanel.append("div").attr("class", "panel-label").text("Module gene counts");

  const mods = td.modules.filter(m => m.color !== "grey");
  const modW = Math.min(400, Math.max(200, mods.length * 25));
  const modH = 120;

  const svg2 = rightPanel.append("svg")
    .attr("width", modW + margin.left + margin.right + 20)
    .attr("height", modH + margin.top + margin.bottom + 10)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const xMod = d3.scaleBand().domain(mods.map(m => m.color)).range([0, modW]).padding(0.15);
  const yMod = d3.scaleLinear().domain([0, d3.max(mods, m => m.n_genes) * 1.1]).range([modH, 0]);

  // Condition-associated module lookup
  const sigMods = new Set((td.condition_modules || []).filter(c => c.p < 0.05).map(c => c.module));

  mods.forEach(m => {{
    svg2.append("rect")
      .attr("x", xMod(m.color)).attr("y", yMod(m.n_genes))
      .attr("width", xMod.bandwidth())
      .attr("height", Math.max(0, modH - yMod(m.n_genes)))
      .attr("fill", wgcnaColorMap[m.color] || "#999")
      .attr("stroke", sigMods.has(m.color) ? "#333" : "#ccc")
      .attr("stroke-width", sigMods.has(m.color) ? 2.5 : 0.5)
      .attr("rx", 1)
      .on("mouseover", function(event) {{
        const cm = (td.condition_modules || []).find(c => c.module === m.color);
        let extra = "";
        if (cm) extra = `<br>r=${{cm.r.toFixed(3)}}, p=${{cm.p.toFixed(4)}}`;
        tooltip.style("display", "block")
          .html(`<b>${{m.color}}</b><br>${{m.n_genes}} genes${{extra}}`);
      }})
      .on("mousemove", function(event) {{
        tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-10)+"px");
      }})
      .on("mouseout", () => tooltip.style("display", "none"));

    // Gene count label
    if (xMod.bandwidth() > 12) {{
      svg2.append("text")
        .attr("x", xMod(m.color) + xMod.bandwidth()/2)
        .attr("y", yMod(m.n_genes) - 2)
        .attr("text-anchor", "middle").attr("font-size", "6px")
        .text(m.n_genes);
    }}
  }});

  // Grey gene count annotation
  const greyMod = td.modules.find(m => m.color === "grey");
  if (greyMod) {{
    svg2.append("text")
      .attr("x", modW).attr("y", 8)
      .attr("text-anchor", "end").attr("font-size", "7px").attr("fill", "#999")
      .text(`grey: ${{greyMod.n_genes}}`);
  }}

  svg2.append("g").attr("transform", `translate(0,${{modH}})`)
    .call(d3.axisBottom(xMod).tickSize(0))
    .selectAll("text").attr("transform", "rotate(-45)").style("text-anchor", "end").style("font-size", "6px");
  svg2.append("g").call(d3.axisLeft(yMod).ticks(4).tickSize(3))
    .selectAll("text").style("font-size", "7px");
}});
</script>
</body>
</html>"""

out_path = os.path.join(OUT_DIR, "FigS3_wgcna_overview.html")
with open(out_path, "w") as f:
    f.write(html)
print(f"Written: {out_path}")
print(f"  {len(TISSUES)} tissue panels with soft threshold + module overview")

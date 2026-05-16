#!/usr/bin/env python3
"""
rrrm1_fig4_integrated.py — Fig4: RRRM-1 scRNA-seq Summary (4 panels)

Panel A: Cell-type composition FLT vs GC (4 tissues, grouped bars)
Panel B: Cell-type pathway NES heatmap (blood, top pathways by variance)
Panel C: LOAO AUROC per cell type vs v1.0 bulk baseline
Panel D: F2-D cross-species concordance bar chart

Usage:
  python3 v2/scripts/rrrm1_fig4_integrated.py

Output:
  v2/figures/Fig4_scrna_summary.html
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent

F2A_JSON = REPO_ROOT / "v2/evaluation/F2A_composition.json"
F2B_JSON = REPO_ROOT / "v2/evaluation/F2B_pseudobulk_fgsea.json"
F2C_JSON = REPO_ROOT / "v2/evaluation/F2C_loao_classifier.json"
F2D_JSON = REPO_ROOT / "v2/evaluation/F2D_crossspecies.json"
FIG_DIR = REPO_ROOT / "v2/figures"

# Okabe-Ito palette
PALETTE = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999",
    "#000000", "#88CCEE",
]


def load_json(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def build_panel_a_data(f2a: dict) -> str:
    """Panel A: mean composition per condition per tissue."""
    results = f2a["results"]
    data = []
    for tissue in ["blood", "eye", "muscle", "skin"]:
        if tissue not in results:
            continue
        for ct, v in results[tissue].items():
            data.append({
                "tissue": tissue,
                "cell_type": ct.replace("_", " "),
                "flt_mean": round(v["flt_mean"] * 100, 1),
                "gc_mean": round(v["gc_mean"] * 100, 1),
                "delta": round(v["delta"] * 100, 1),
                "p_raw": v["p_raw"],
                "padj": v["padj"],
            })
    return json.dumps(data)


def build_panel_b_data(f2b: dict) -> str:
    """Panel B: NES heatmap for blood cell types, top 15 pathways by variance."""
    blood = f2b["results"].get("blood", {})
    muscle = f2b["results"].get("muscle", {})

    # Use blood — collect NES per cell type per pathway
    # We need the full NES CSVs which are stored locally now
    blood_dir = REPO_ROOT / "v2/processed/F2B_blood"
    nes_dict = {}
    for csv_path in sorted(blood_dir.glob("*_fgsea_hallmark.csv")):
        ct = csv_path.stem.replace("_fgsea_hallmark", "").replace("_", " ")
        df = pd.read_csv(csv_path)
        nes_dict[ct] = df.set_index("pathway")["NES"].to_dict()

    if not nes_dict:
        return json.dumps({"cell_types": [], "pathways": [], "nes": []})

    # Build matrix
    cts = sorted(nes_dict.keys())
    all_pws = set()
    for v in nes_dict.values():
        all_pws.update(v.keys())
    all_pws = sorted(all_pws)

    # Select top 15 by variance across cell types
    variance = {}
    for pw in all_pws:
        vals = [nes_dict[ct].get(pw, 0) for ct in cts]
        variance[pw] = np.var(vals)
    top_pws = sorted(variance, key=lambda x: -variance[x])[:15]

    # Build NES matrix
    nes_matrix = []
    for pw in top_pws:
        row = [round(nes_dict[ct].get(pw, 0), 3) for ct in cts]
        nes_matrix.append(row)

    return json.dumps({
        "cell_types": cts,
        "pathways": [pw.replace("HALLMARK_", "").replace("_", " ").title() for pw in top_pws],
        "pathways_raw": top_pws,
        "nes": nes_matrix,
    })


def build_panel_c_data(f2c: dict) -> str:
    """Panel C: AUROC per cell type, all tissues."""
    data = []
    bulk = f2c.get("v1_bulk_auroc", {})
    for tissue in ["blood", "eye", "muscle", "skin"]:
        if tissue not in f2c["results"]:
            continue
        for ct, v in f2c["results"][tissue].items():
            sig = v["p_perm"] < 0.05
            data.append({
                "tissue": tissue,
                "cell_type": ct.replace("_", " "),
                "auroc": round(v["auroc"], 3),
                "ci_low": round(v["ci_low"], 3),
                "ci_high": round(v["ci_high"], 3),
                "p_perm": v["p_perm"],
                "sig": sig,
                "n_cells": v["n_cells_total"],
            })
        # Add bulk baseline
        if tissue in bulk:
            data.append({
                "tissue": tissue,
                "cell_type": "v1.0 bulk",
                "auroc": bulk[tissue],
                "ci_low": None,
                "ci_high": None,
                "p_perm": None,
                "sig": True,
                "n_cells": None,
                "is_bulk": True,
            })
    return json.dumps(data)


def build_panel_d_data(f2d: dict) -> str:
    """Panel D: F2-D concordance bar chart."""
    data = []
    e1_r = f2d["e1_baseline"]["spearman_r"]

    for m in f2d["matched_pairs"]:
        if m.get("skipped"):
            continue
        data.append({
            "label": m["matched_lineage"],
            "r": round(m["spearman_r"], 3),
            "ci_low": round(m["ci_low"], 3),
            "ci_high": round(m["ci_high"], 3),
            "n": m["n_pathways"],
            "type": "matched",
            "exceeds_e1": m.get("exceeds_e1", False),
        })

    for s in f2d["sensitivity"]:
        if s.get("skipped"):
            continue
        data.append({
            "label": s["matched_lineage"],
            "r": round(s["spearman_r"], 3),
            "ci_low": round(s["ci_low"], 3),
            "ci_high": round(s["ci_high"], 3),
            "n": s["n_pathways"],
            "type": "sensitivity",
        })

    return json.dumps({"bars": data, "e1_r": round(e1_r, 3)})


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Fig4: RRRM-1 scRNA-seq Summary ===\n")

    f2a = load_json(F2A_JSON)
    f2b = load_json(F2B_JSON)
    f2c = load_json(F2C_JSON)
    f2d = load_json(F2D_JSON)

    panel_a = build_panel_a_data(f2a)
    panel_b = build_panel_b_data(f2b)
    panel_c = build_panel_c_data(f2c)
    panel_d = build_panel_d_data(f2d)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig4: RRRM-1 scRNA-seq Benchmark Summary</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; font-size: 12px; background: #fff; margin: 16px; }}
  .row {{ display: flex; flex-wrap: wrap; gap: 16px; margin-bottom: 16px; }}
  .panel {{ background: #fff; border: 1px solid #e0e0e0; padding: 10px; }}
  .panel-label {{ font-size: 14px; font-weight: bold; margin-bottom: 4px; }}
  .panel-title {{ font-size: 11px; font-weight: bold; margin-bottom: 2px; }}
  .panel-sub {{ font-size: 8px; color: #666; margin-bottom: 8px; }}
  .axis text {{ font-size: 8px; }}
  .axis line, .axis path {{ stroke: #333; }}
  .tooltip {{
    position: absolute; background: rgba(255,255,255,0.95);
    border: 1px solid #ccc; padding: 5px; border-radius: 3px;
    font-size: 9px; pointer-events: none; display: none;
  }}
  .ref-line {{ stroke: #D55E00; stroke-width: 1.5; stroke-dasharray: 5,3; }}
  .hm-cell {{ stroke: #fff; stroke-width: 0.5; }}
  #download-btn {{
    margin-top: 10px; padding: 5px 12px; cursor: pointer;
    font-size: 11px; border: 1px solid #555; border-radius: 3px; background: #f5f5f5;
  }}
</style>
</head>
<body>
<h2 style="font-size:14px;margin-bottom:4px">Figure 4: RRRM-1 Single-Cell RNA-seq Benchmark</h2>
<p style="font-size:9px;color:#555;margin-top:0">
  (A) Cell-type composition FLT vs GC. (B) Blood cell-type Hallmark NES heatmap.
  (C) LOAO classifier AUROC by cell type. (D) Cross-species cell-type pathway concordance.
</p>

<div class="tooltip" id="tooltip"></div>

<div class="row">
  <div class="panel" id="panelA"></div>
  <div class="panel" id="panelB"></div>
</div>
<div class="row">
  <div class="panel" id="panelC"></div>
  <div class="panel" id="panelD"></div>
</div>

<button id="download-btn">Download all panels (SVG)</button>

<script>
const palette = {json.dumps(PALETTE)};
const tooltip = d3.select("#tooltip");

// ═══ Panel A: Composition (grouped bar per tissue) ═══════════════════
(function() {{
  const data = {panel_a};
  const margin = {{top: 24, right: 10, bottom: 70, left: 50}};
  const W = 560, H = 340;
  const w = W - margin.left - margin.right;
  const h = H - margin.top - margin.bottom;

  const div = d3.select("#panelA");
  div.append("div").attr("class","panel-label").text("A");
  div.append("div").attr("class","panel-title").text("Cell-type composition (FLT vs GC)");
  div.append("div").attr("class","panel-sub").text("Mean proportion (%). Wilcoxon rank-sum, BH FDR. n=4 FLT vs n=4 GC per tissue.");

  const svg = div.append("svg").attr("width",W).attr("height",H)
    .append("g").attr("transform",`translate(${{margin.left}},${{margin.top}})`);

  // Group by tissue
  const tissues = [...new Set(data.map(d=>d.tissue))];
  const cellTypes = [...new Set(data.map(d=>d.cell_type))];

  // Stacked bar per tissue: show top 5 cell types per tissue by mean proportion
  tissues.forEach((tissue, ti) => {{
    const tData = data.filter(d => d.tissue === tissue)
      .sort((a,b) => (b.flt_mean + b.gc_mean)/2 - (a.flt_mean + a.gc_mean)/2)
      .slice(0, 6);

    const panelW = w / tissues.length - 10;
    const offsetX = ti * (panelW + 10);

    const x = d3.scaleBand()
      .domain(tData.map(d=>d.cell_type))
      .range([0, panelW]).padding(0.15);
    const y = d3.scaleLinear().domain([0, 100]).range([h, 0]);

    const g = svg.append("g").attr("transform", `translate(${{offsetX}},0)`);

    // FLT bars (left half)
    g.selectAll(".flt").data(tData).enter().append("rect")
      .attr("x", d => x(d.cell_type))
      .attr("y", d => y(d.flt_mean))
      .attr("width", x.bandwidth()/2)
      .attr("height", d => h - y(d.flt_mean))
      .attr("fill", "#D55E00").attr("opacity", 0.8);

    // GC bars (right half)
    g.selectAll(".gc").data(tData).enter().append("rect")
      .attr("x", d => x(d.cell_type) + x.bandwidth()/2)
      .attr("y", d => y(d.gc_mean))
      .attr("width", x.bandwidth()/2)
      .attr("height", d => h - y(d.gc_mean))
      .attr("fill", "#0072B2").attr("opacity", 0.8);

    // X axis
    g.append("g").attr("transform",`translate(0,${{h}})`).call(d3.axisBottom(x).tickSize(2))
      .selectAll("text").attr("transform","rotate(-40)").style("text-anchor","end").style("font-size","7px");

    // Tissue label
    g.append("text").attr("x", panelW/2).attr("y", -6)
      .attr("text-anchor","middle").attr("font-size","9px").attr("font-weight","bold")
      .text(tissue.charAt(0).toUpperCase() + tissue.slice(1));
  }});

  // Y axis
  svg.append("g").attr("class","axis").call(d3.axisLeft(d3.scaleLinear().domain([0,100]).range([h,0])).ticks(5).tickFormat(d=>d+"%").tickSize(3));
  svg.append("text").attr("transform","rotate(-90)").attr("x",-h/2).attr("y",-38)
    .attr("text-anchor","middle").attr("font-size","9px").text("Mean proportion (%)");

  // Legend
  const lg = svg.append("g").attr("transform",`translate(${{w-90}},${{-10}})`);
  lg.append("rect").attr("width",10).attr("height",10).attr("fill","#D55E00").attr("opacity",0.8);
  lg.append("text").attr("x",14).attr("y",9).attr("font-size","8px").text("FLT");
  lg.append("rect").attr("x",50).attr("width",10).attr("height",10).attr("fill","#0072B2").attr("opacity",0.8);
  lg.append("text").attr("x",64).attr("y",9).attr("font-size","8px").text("GC");
}})();

// ═══ Panel B: NES Heatmap (blood cell types) ═════════════════════════
(function() {{
  const data = {panel_b};
  const margin = {{top: 24, right: 30, bottom: 10, left: 180}};
  const cellW = 55, cellH = 18;
  const nC = data.cell_types.length, nR = data.pathways.length;
  const W = margin.left + cellW * nC + margin.right;
  const H = margin.top + cellH * nR + margin.bottom + 60;

  const div = d3.select("#panelB");
  div.append("div").attr("class","panel-label").text("B");
  div.append("div").attr("class","panel-title").text("Blood cell-type Hallmark NES");
  div.append("div").attr("class","panel-sub").text("Top 15 pathways by variance across cell types. Pseudo-bulk fGSEA (FLT vs GC).");

  if (!data.pathways.length) {{ div.append("p").text("No data"); return; }}

  const svg = div.append("svg").attr("width", W).attr("height", H)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top + 50}})`);

  const color = d3.scaleDiverging(d3.interpolateRdBu).domain([2.5, 0, -2.5]);

  data.pathways.forEach((pw, i) => {{
    data.cell_types.forEach((ct, j) => {{
      const nes = data.nes[i][j];
      svg.append("rect")
        .attr("x", j*cellW).attr("y", i*cellH)
        .attr("width", cellW).attr("height", cellH)
        .attr("class","hm-cell")
        .attr("fill", color(nes));
      svg.append("text")
        .attr("x", j*cellW + cellW/2).attr("y", i*cellH + cellH/2 + 3)
        .attr("text-anchor","middle").attr("font-size","7px")
        .attr("fill", Math.abs(nes) > 1.5 ? "#fff" : "#333")
        .text(nes.toFixed(1));
    }});
  }});

  // Row labels
  data.pathways.forEach((pw, i) => {{
    svg.append("text").attr("x",-6).attr("y", i*cellH + cellH/2 + 3)
      .attr("text-anchor","end").attr("font-size","8px").text(pw);
  }});

  // Col labels
  data.cell_types.forEach((ct, j) => {{
    svg.append("text")
      .attr("x", j*cellW + cellW/2).attr("y",-6)
      .attr("text-anchor","end").attr("font-size","8px")
      .attr("transform", `rotate(-40, ${{j*cellW+cellW/2}}, -6)`)
      .text(ct);
  }});
}})();

// ═══ Panel C: AUROC bar chart ════════════════════════════════════════
(function() {{
  const data = {panel_c};
  const margin = {{top: 24, right: 20, bottom: 80, left: 50}};
  const W = 560, H = 340;
  const w = W - margin.left - margin.right;
  const h = H - margin.top - margin.bottom;

  const div = d3.select("#panelC");
  div.append("div").attr("class","panel-label").text("C");
  div.append("div").attr("class","panel-title").text("LOAO classifier AUROC per cell type");
  div.append("div").attr("class","panel-sub").text("PCA(50)-LR, leave-one-animal-out. Dashed orange = v1.0 bulk baseline. Error bars = 95% bootstrap CI.");

  const svg = div.append("svg").attr("width",W).attr("height",H)
    .append("g").attr("transform",`translate(${{margin.left}},${{margin.top}})`);

  // Focus on key tissues: blood + muscle (skip eye/skin for clarity)
  const tissues = ["blood", "muscle"];
  const tissuePalette = {{"blood": "#56B4E9", "muscle": "#009E73"}};
  const bulkRef = {{}};

  tissues.forEach((tissue, ti) => {{
    const tData = data.filter(d => d.tissue === tissue && !d.is_bulk)
      .sort((a,b) => b.auroc - a.auroc);
    const bulkRow = data.find(d => d.tissue === tissue && d.is_bulk);
    if (bulkRow) bulkRef[tissue] = bulkRow.auroc;

    const panelW = w / tissues.length - 10;
    const offsetX = ti * (panelW + 10);

    const x = d3.scaleBand().domain(tData.map(d=>d.cell_type)).range([0, panelW]).padding(0.2);
    const y = d3.scaleLinear().domain([0, 1]).range([h, 0]);

    const g = svg.append("g").attr("transform", `translate(${{offsetX}},0)`);

    // Bars
    g.selectAll(".bar").data(tData).enter().append("rect")
      .attr("x", d => x(d.cell_type))
      .attr("y", d => y(d.auroc))
      .attr("width", x.bandwidth())
      .attr("height", d => h - y(d.auroc))
      .attr("fill", d => d.sig ? tissuePalette[tissue] : "#ccc")
      .attr("opacity", 0.8);

    // CI error bars
    tData.forEach(d => {{
      if (d.ci_low && d.ci_high) {{
        const cx = x(d.cell_type) + x.bandwidth()/2;
        g.append("line").attr("x1",cx).attr("x2",cx)
          .attr("y1",y(d.ci_low)).attr("y2",y(d.ci_high))
          .attr("stroke","#333").attr("stroke-width",0.8);
      }}
    }});

    // Bulk baseline
    if (bulkRef[tissue]) {{
      g.append("line").attr("class","ref-line")
        .attr("x1",0).attr("x2",panelW)
        .attr("y1",y(bulkRef[tissue])).attr("y2",y(bulkRef[tissue]));
      g.append("text").attr("x",panelW-2).attr("y",y(bulkRef[tissue])-4)
        .attr("text-anchor","end").attr("font-size","7px").attr("fill","#D55E00")
        .text(`bulk ${{bulkRef[tissue].toFixed(3)}}`);
    }}

    // Chance line
    g.append("line").attr("x1",0).attr("x2",panelW)
      .attr("y1",y(0.5)).attr("y2",y(0.5))
      .attr("stroke","#999").attr("stroke-dasharray","2,2").attr("stroke-width",0.5);

    // X axis
    g.append("g").attr("transform",`translate(0,${{h}})`).call(d3.axisBottom(x).tickSize(2))
      .selectAll("text").attr("transform","rotate(-40)").style("text-anchor","end").style("font-size","7px");

    // Tissue label
    g.append("text").attr("x", panelW/2).attr("y", -6)
      .attr("text-anchor","middle").attr("font-size","9px").attr("font-weight","bold")
      .text(tissue.charAt(0).toUpperCase() + tissue.slice(1));
  }});

  // Y axis
  svg.append("g").attr("class","axis").call(d3.axisLeft(d3.scaleLinear().domain([0,1]).range([h,0])).ticks(5).tickSize(3));
  svg.append("text").attr("transform","rotate(-90)").attr("x",-h/2).attr("y",-38)
    .attr("text-anchor","middle").attr("font-size","9px").text("AUROC");
}})();

// ═══ Panel D: Cross-species concordance bar chart ════════════════════
(function() {{
  const data = {panel_d};
  const margin = {{top: 24, right: 20, bottom: 80, left: 50}};
  const W = 350, H = 340;
  const w = W - margin.left - margin.right;
  const h = H - margin.top - margin.bottom;

  const div = d3.select("#panelD");
  div.append("div").attr("class","panel-label").text("D");
  div.append("div").attr("class","panel-title").text("Cross-species cell-type NES concordance");
  div.append("div").attr("class","panel-sub").text("Spearman r: RRRM-1 mouse blood vs I4 human PBMC. Dashed = E1 bulk baseline.");

  const svg = div.append("svg").attr("width",W).attr("height",H)
    .append("g").attr("transform",`translate(${{margin.left}},${{margin.top}})`);

  const bars = data.bars;
  const e1_r = data.e1_r;

  const x = d3.scaleBand().domain(bars.map(d=>d.label)).range([0,w]).padding(0.25);
  const y = d3.scaleLinear().domain([-0.5, 1]).range([h, 0]);

  // Zero line
  svg.append("line").attr("x1",0).attr("x2",w)
    .attr("y1",y(0)).attr("y2",y(0)).attr("stroke","#999").attr("stroke-width",0.5);

  // E1 baseline
  svg.append("line").attr("class","ref-line")
    .attr("x1",0).attr("x2",w).attr("y1",y(e1_r)).attr("y2",y(e1_r));
  svg.append("text").attr("x",w-2).attr("y",y(e1_r)-4)
    .attr("text-anchor","end").attr("font-size","7px").attr("fill","#D55E00")
    .text(`E1 r = ${{e1_r}}`);

  // Bars
  svg.selectAll(".bar").data(bars).enter().append("rect")
    .attr("x", d => x(d.label))
    .attr("y", d => d.r >= 0 ? y(d.r) : y(0))
    .attr("width", x.bandwidth())
    .attr("height", d => Math.abs(y(0) - y(d.r)))
    .attr("fill", d => d.type === "matched" ? "#0072B2" : "#88CCEE")
    .attr("opacity", 0.8);

  // CI error bars
  bars.forEach(d => {{
    const cx = x(d.label) + x.bandwidth()/2;
    if (d.ci_low != null && d.ci_high != null) {{
      svg.append("line").attr("x1",cx).attr("x2",cx)
        .attr("y1",y(d.ci_low)).attr("y2",y(d.ci_high))
        .attr("stroke","#333").attr("stroke-width",0.8);
      svg.append("line").attr("x1",cx-3).attr("x2",cx+3)
        .attr("y1",y(d.ci_low)).attr("y2",y(d.ci_low)).attr("stroke","#333");
      svg.append("line").attr("x1",cx-3).attr("x2",cx+3)
        .attr("y1",y(d.ci_high)).attr("y2",y(d.ci_high)).attr("stroke","#333");
    }}
    // n label
    const yPos = d.r >= 0 ? y(Math.max(d.ci_high||d.r, d.r)) - 10 : y(Math.min(d.ci_low||d.r, d.r)) + 14;
    svg.append("text").attr("x",cx).attr("y",yPos)
      .attr("text-anchor","middle").attr("font-size","7px").attr("fill","#555")
      .text(`n=${{d.n}}`);
  }});

  // Axes
  svg.append("g").attr("class","axis").attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(x).tickSize(2))
    .selectAll("text").attr("transform","rotate(-30)").style("text-anchor","end").style("font-size","7px");
  svg.append("g").attr("class","axis").call(d3.axisLeft(y).ticks(6).tickSize(3));
  svg.append("text").attr("transform","rotate(-90)").attr("x",-h/2).attr("y",-38)
    .attr("text-anchor","middle").attr("font-size","9px").text("Spearman r");

  // Legend
  const lg = svg.append("g").attr("transform",`translate(10,10)`);
  lg.append("rect").attr("width",10).attr("height",10).attr("fill","#0072B2").attr("opacity",0.8);
  lg.append("text").attr("x",14).attr("y",9).attr("font-size","7px").text("Matched pair");
  lg.append("rect").attr("y",14).attr("width",10).attr("height",10).attr("fill","#88CCEE").attr("opacity",0.8);
  lg.append("text").attr("x",14).attr("y",23).attr("font-size","7px").text("Sensitivity");
}})();

// ─── Download ────────────────────────────────────────────────────────
document.getElementById("download-btn").addEventListener("click", function() {{
  const svgs = document.querySelectorAll("svg");
  const serializer = new XMLSerializer();
  let combined = '<?xml version="1.0"?>\\n';
  svgs.forEach((s, i) => {{
    combined += `<!-- Panel ${{["A","B","C","D"][i]}} -->\\n` + serializer.serializeToString(s) + "\\n";
  }});
  const link = document.createElement("a");
  link.download = "Fig4_scrna_summary.svg";
  link.href = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(combined);
  link.click();
}});
</script>
</body>
</html>"""

    out_path = FIG_DIR / "Fig4_scrna_summary.html"
    out_path.write_text(html)
    print(f"  Fig4 saved: {out_path}")
    print("\n=== Fig4 Complete ===")


if __name__ == "__main__":
    main()

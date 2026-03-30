#!/usr/bin/env python
"""Fig S7: Metabolic Flux Detail.

Extended heatmap of subsystem-level flux changes (FLT − GC) across 6 LOMO tissues.
Shows top-30 metabolic subsystems ranked by absolute mean flux difference.
"""
import json
import os
from pathlib import Path

import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent.parent
V5_EVAL_DIR = BASE_DIR / "v5" / "evaluation"
OUT_DIR = BASE_DIR / "v5" / "figures" / "html"

FLUX_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin"]
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastrocnemius", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin",
}


def load_flux_data():
    """Load metabolic flux results and build subsystem heatmap."""
    tissue_data = {}
    all_subsystems = {}  # subsystem → max abs_mean_diff across tissues

    for tissue in FLUX_TISSUES:
        path = V5_EVAL_DIR / f"metabolic_flux_{tissue}.json"
        if not path.exists():
            continue
        with open(path) as f:
            d = json.load(f)
        if d.get("status") != "success":
            continue

        subs = d.get("all_subsystems", d.get("top_subsystems", {}))
        tissue_data[tissue] = {
            "subsystems": subs,
            "coverage": d.get("gene_coverage_pct", 0),
            "obj_diff": d.get("fba_objective_diff", 0),
            "n_flt": d.get("n_samples_flt", 0),
            "n_gc": d.get("n_samples_gc", 0),
        }

        for ss, info in subs.items():
            adiff = abs(info.get("mean_flux_diff", 0))
            if ss not in all_subsystems or adiff > all_subsystems[ss]:
                all_subsystems[ss] = adiff

    # Top-40 subsystems by max absolute diff across any tissue
    top_ss = sorted(all_subsystems.items(), key=lambda x: -x[1])[:40]
    top_ss_names = [s[0] for s in top_ss]

    return tissue_data, top_ss_names


def build_heatmap_data(tissue_data, subsystems):
    """Build heatmap cells and summary stats."""
    cells = []
    summary_data = []

    for tissue in FLUX_TISSUES:
        td = tissue_data.get(tissue)
        if not td:
            continue
        subs = td["subsystems"]
        label = TISSUE_LABELS[tissue]

        summary_data.append({
            "tissue": label,
            "coverage": td["coverage"],
            "obj_diff": round(td["obj_diff"], 4),
            "n_flt": td["n_flt"],
            "n_gc": td["n_gc"],
        })

        for ss in subsystems:
            info = subs.get(ss, {})
            mfd = info.get("mean_flux_diff", 0)
            cells.append({
                "tissue": label,
                "subsystem": ss,
                "mean_flux_diff": round(float(mfd), 4),
                "abs_mean_diff": round(float(info.get("abs_mean_diff", 0)), 4),
                "n_reactions": info.get("n_reactions", 0),
                "n_increased": info.get("n_increased", 0),
                "n_decreased": info.get("n_decreased", 0),
            })

    return cells, summary_data


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    tissue_data, subsystems = load_flux_data()
    cells, summary_data = build_heatmap_data(tissue_data, subsystems)

    print(f"Fig S7: {len(cells)} cells ({len(subsystems)} subsystems × {len(tissue_data)} tissues)")
    print(f"Summary: {len(summary_data)} tissues")

    tissues_list = [TISSUE_LABELS[t] for t in FLUX_TISSUES if t in tissue_data]

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig S7: Metabolic Flux Subsystem Detail</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
body {{ font-family: Arial, Helvetica, sans-serif; margin: 0; padding: 20px; background: #fff; }}
.figure-container {{ width: 1400px; margin: 0 auto; }}
.figure-title {{ font-size: 14px; font-weight: bold; text-align: center; margin-bottom: 4px; }}
.figure-subtitle {{ font-size: 10px; color: #666; text-align: center; margin-bottom: 16px; }}
.panels {{ display: grid; grid-template-columns: 1fr; gap: 16px; }}
.panel {{ border: 1px solid #e0e0e0; border-radius: 4px; padding: 12px; background: #fafafa; }}
.panel-label {{ font-size: 12px; font-weight: bold; margin-bottom: 2px; }}
.panel-desc {{ font-size: 9px; color: #666; margin-bottom: 8px; }}
.tooltip {{ position: absolute; background: rgba(0,0,0,0.85); color: #fff; padding: 6px 10px;
            border-radius: 3px; font-size: 9px; pointer-events: none; z-index: 100; max-width: 300px; }}
.interpretation {{ margin-top: 16px; padding: 12px; background: #f0f7ff; border-radius: 4px;
                   font-size: 10px; line-height: 1.5; }}
.summary-table {{ font-size: 9px; border-collapse: collapse; margin: 8px auto; }}
.summary-table th, .summary-table td {{ border: 1px solid #ddd; padding: 4px 8px; text-align: center; }}
.summary-table th {{ background: #f0f0f0; font-weight: bold; }}
</style>
</head>
<body>
<div class="figure-container">
  <div class="figure-title">Fig S7: E-Flux Metabolic Subsystem Flux Changes During Spaceflight</div>
  <div class="figure-subtitle">iMM1865 mouse GEM constrained by tissue-specific gene expression (FLT vs GC). Top-40 subsystems by max |flux difference|.</div>

  <div class="panels">
    <div class="panel" id="panelA">
      <div class="panel-label">A</div>
      <div class="panel-desc">Subsystem flux difference heatmap (mean flux FLT − GC, clipped to [−500, +500])</div>
      <div id="flux-heatmap"></div>
    </div>
    <div class="panel" id="panelB">
      <div class="panel-label">B</div>
      <div class="panel-desc">Gene-to-model coverage and FBA objective summary</div>
      <div id="summary-table-container"></div>
    </div>
  </div>

  <div class="interpretation">
    <strong>Key findings:</strong>
    E-Flux metabolic modeling reveals tissue-specific flux reprogramming during spaceflight.
    Vitamin B12, glutamate, and citric acid cycle subsystems show the largest flux changes across tissues.
    Most tissues converge to identical FBA biomass objectives (798.81), indicating that spaceflight-induced
    expression changes redistribute metabolic flux rather than reducing total metabolic capacity.
    Skin is the exception (FBA objective FLT=782 vs GC=787), suggesting genuine metabolic stress.
    Gene coverage ranges from 84.6% (gastrocnemius) to 92.9% (thymus), indicating robust model-data alignment.
  </div>
</div>
<div class="tooltip" id="tooltip" style="display:none;"></div>

<script>
const OI = {{
  black: "#000000", orange: "#E69F00", skyblue: "#56B4E9", green: "#009E73",
  yellow: "#F0E442", blue: "#0072B2", vermilion: "#D55E00", pink: "#CC79A7"
}};
const tooltip = d3.select("#tooltip");

// Panel A: Subsystem flux heatmap
(function() {{
  const data = {json.dumps(cells)};
  const tissues = {json.dumps(tissues_list)};
  const subsystems = {json.dumps(subsystems)};

  const margin = {{top: 10, right: 80, bottom: 60, left: 280}};
  const cellW = 80, cellH = 18;
  const width = tissues.length * cellW, height = subsystems.length * cellH;

  const svg = d3.select("#flux-heatmap").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(tissues).range([0, width]).padding(0.05);
  const y = d3.scaleBand().domain(subsystems).range([0, height]).padding(0.03);

  // Clip extreme values for color scale
  const clipVal = 500;
  const color = d3.scaleDiverging(d3.interpolateRdBu)
    .domain([clipVal, 0, -clipVal]);

  svg.selectAll("rect").data(data).join("rect")
    .attr("x", d => x(d.tissue)).attr("y", d => y(d.subsystem))
    .attr("width", x.bandwidth()).attr("height", y.bandwidth())
    .attr("fill", d => {{
      const v = Math.max(-clipVal, Math.min(clipVal, d.mean_flux_diff));
      return Math.abs(d.mean_flux_diff) < 0.001 ? "#f5f5f5" : color(v);
    }})
    .attr("stroke", "#eee").attr("stroke-width", 0.3)
    .on("mouseover", (e, d) => {{
      tooltip.style("display", "block")
        .html(`<b>${{d.subsystem}}</b> in ${{d.tissue}}<br>
               Mean flux diff: ${{d.mean_flux_diff}}<br>
               |Mean diff|: ${{d.abs_mean_diff}}<br>
               Reactions: ${{d.n_reactions}} (${{d.n_increased}}↑ ${{d.n_decreased}}↓)`)
        .style("left", (e.pageX+10)+"px").style("top", (e.pageY-10)+"px");
    }}).on("mouseout", () => tooltip.style("display", "none"));

  svg.append("g").attr("transform", `translate(0,${{height}})`).call(d3.axisBottom(x))
    .selectAll("text").attr("transform", "rotate(-30)").style("text-anchor", "end").style("font-size", "9px");
  svg.append("g").call(d3.axisLeft(y)).selectAll("text").style("font-size", "7px");

  // Color legend
  const legendW = 15, legendH = height;
  const legendG = svg.append("g").attr("transform", `translate(${{width + 20}}, 0)`);
  const legendScale = d3.scaleLinear().domain([-clipVal, clipVal]).range([legendH, 0]);
  const legendAxis = d3.axisRight(legendScale).ticks(5);

  const gradient = svg.append("defs").append("linearGradient")
    .attr("id", "flux-gradient").attr("x1", "0%").attr("y1", "100%").attr("x2", "0%").attr("y2", "0%");
  gradient.append("stop").attr("offset", "0%").attr("stop-color", d3.interpolateRdBu(0));
  gradient.append("stop").attr("offset", "50%").attr("stop-color", d3.interpolateRdBu(0.5));
  gradient.append("stop").attr("offset", "100%").attr("stop-color", d3.interpolateRdBu(1));

  legendG.append("rect").attr("width", legendW).attr("height", legendH)
    .style("fill", "url(#flux-gradient)");
  legendG.append("g").attr("transform", `translate(${{legendW}},0)`).call(legendAxis)
    .selectAll("text").style("font-size", "7px");
  legendG.append("text").attr("transform", `rotate(-90) translate(${{-legendH/2}}, -5)`)
    .attr("text-anchor", "middle").attr("font-size", "8px").text("FLT − GC flux");
}})();

// Panel B: Summary table
(function() {{
  const summary = {json.dumps(summary_data)};
  let html = '<table class="summary-table"><tr><th>Tissue</th><th>n FLT</th><th>n GC</th><th>Gene Coverage (%)</th><th>FBA Obj Diff (FLT−GC)</th></tr>';
  summary.forEach(s => {{
    html += `<tr><td>${{s.tissue}}</td><td>${{s.n_flt}}</td><td>${{s.n_gc}}</td><td>${{s.coverage}}%</td><td>${{s.obj_diff}}</td></tr>`;
  }});
  html += '</table>';
  document.getElementById("summary-table-container").innerHTML = html;
}})();
</script>
</body>
</html>"""

    out_path = OUT_DIR / "FigS7_metabolic_detail.html"
    with open(out_path, "w") as f:
        f.write(html)
    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()

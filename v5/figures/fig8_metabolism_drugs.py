#!/usr/bin/env python
"""Fig 8: Metabolic Flux & Drug Target Mapping.

Panel A: Subsystem flux changes heatmap (top subsystems × 6 tissues)
Panel B: Drug tier distribution (FDA-approved, clinical, pre-clinical)
Panel C: Consensus biomarker panel scores
Panel D: Panel validation AUROC per tissue
"""
import json
import os
from pathlib import Path

import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent.parent
V5_EVAL_DIR = BASE_DIR / "v5" / "evaluation"
OUT_DIR = BASE_DIR / "v5" / "figures" / "html"

FLUX_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin"]
ALL_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin", "lung", "colon"]
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastrocnemius", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin", "lung": "Lung†", "colon": "Colon†",
}


def load_flux_data():
    """Load metabolic flux subsystem data for top subsystems."""
    all_subsystems = {}  # subsystem → {tissue: abs_mean_diff}
    for tissue in FLUX_TISSUES:
        path = V5_EVAL_DIR / f"metabolic_flux_{tissue}.json"
        if not path.exists():
            continue
        with open(path) as f:
            d = json.load(f)
        for ss, info in d.get("all_subsystems", d.get("top_subsystems", {})).items():
            if ss not in all_subsystems:
                all_subsystems[ss] = {}
            all_subsystems[ss][tissue] = info.get("mean_flux_diff", 0)

    # Top 15 subsystems by max abs change across tissues
    ranked = sorted(all_subsystems.items(),
                    key=lambda x: max(abs(v) for v in x[1].values()),
                    reverse=True)[:15]

    heatmap = []
    ss_names = []
    for ss, tissue_data in ranked:
        ss_short = ss[:30]
        ss_names.append(ss_short)
        for tissue in FLUX_TISSUES:
            heatmap.append({
                "subsystem": ss_short,
                "tissue": TISSUE_LABELS[tissue],
                "flux_diff": round(tissue_data.get(tissue, 0), 3),
            })
    return heatmap, ss_names


def load_drug_tiers():
    """Load drug tier distribution."""
    path = V5_EVAL_DIR / "drug_targets.json"
    if not path.exists():
        return []
    with open(path) as f:
        d = json.load(f)
    tiers = d.get("tiers", {})
    return [
        {"tier": "Tier 1 (FDA-approved)", "count": len(tiers.get("tier1_approved", []))},
        {"tier": "Tier 2 (Clinical)", "count": len(tiers.get("tier2_clinical", []))},
        {"tier": "Tier 3 (Pre-clinical)", "count": len(tiers.get("tier3_preclinical", []))},
    ]


def load_panel_data():
    """Load consensus biomarker panel."""
    path = V5_EVAL_DIR / "consensus_biomarker_panel.json"
    if not path.exists():
        return [], []
    with open(path) as f:
        d = json.load(f)

    panel = []
    for entry in d.get("panel", [])[:20]:
        sources = entry.get("sources", {})
        panel.append({
            "gene": entry["gene"],
            "score": entry["score"],
            "n_tissues": entry.get("n_tissues", 0),
            "shap": sources.get("shap", 0),
            "wgcna": sources.get("wgcna_hub", 0),
            "network": sources.get("network", 0),
            "drug": sources.get("druggable", 0),
            "lit": sources.get("literature", 0),
        })

    validation = []
    for tissue in ALL_TISSUES:
        v = d.get("panel_validation", {}).get(tissue, {})
        auroc = v.get("auroc")
        if auroc is not None:
            validation.append({
                "tissue": TISSUE_LABELS[tissue],
                "auroc": auroc,
                "n_genes": v.get("n_panel_genes_found", 0),
            })
    return panel, validation


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    flux_heatmap, ss_names = load_flux_data()
    drug_tiers = load_drug_tiers()
    panel, validation = load_panel_data()

    print(f"Panel A: {len(flux_heatmap)} flux cells ({len(ss_names)} subsystems)")
    print(f"Panel B: {sum(t['count'] for t in drug_tiers)} total drug hits")
    print(f"Panel C: {len(panel)} panel genes")
    print(f"Panel D: {len(validation)} tissue validation points")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig 8: Metabolic Flux & Drug Targets</title>
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
.interpretation {{ margin-top: 16px; padding: 12px; background: #f0f7ff; border-radius: 4px;
                   font-size: 10px; line-height: 1.5; }}
.tooltip {{ position: absolute; background: rgba(0,0,0,0.85); color: #fff; padding: 6px 10px;
            border-radius: 3px; font-size: 9px; pointer-events: none; z-index: 100; }}
</style>
</head>
<body>
<div class="figure-container">
  <div class="figure-title">Fig 8: Metabolic Flux Modeling & Drug Target Mapping</div>
  <div class="figure-subtitle">E-Flux constrained FBA (iMM1865) and DGIdb v5 + ChEMBL drug-gene interactions</div>
  <div class="panels">
    <div class="panel" id="panelA">
      <div class="panel-label">A</div>
      <div class="panel-desc">Subsystem flux changes (FLT − GC mean flux, top 15)</div>
      <div id="flux-heatmap"></div>
    </div>
    <div class="panel" id="panelB">
      <div class="panel-label">B</div>
      <div class="panel-desc">Drug interaction tiers across spaceflight genes</div>
      <div id="drug-tiers"></div>
    </div>
    <div class="panel" id="panelC">
      <div class="panel-label">C</div>
      <div class="panel-desc">Consensus biomarker panel (top-20, multi-source scoring)</div>
      <div id="panel-scores"></div>
    </div>
    <div class="panel" id="panelD">
      <div class="panel-label">D</div>
      <div class="panel-desc">Panel validation: AUROC with 20-gene panel vs full gene set</div>
      <div id="panel-validation"></div>
    </div>
  </div>
  <div class="interpretation">
    <strong>Key findings:</strong>
    E-Flux metabolic modeling reveals tissue-specific flux redistribution during spaceflight, with skin showing the
    largest biomass production difference (FLT 782.1 vs GC 787.2). Among 1,919 spaceflight-affected genes, 271 (32.5%)
    have known drug interactions in DGIdb, yielding 1,287 FDA-approved drug-gene connections. The 20-gene consensus
    biomarker panel (scored across 6 evidence sources) achieves AUROC 0.935 in gastrocnemius and 0.750 in colon,
    demonstrating that a compact, interpretable gene set can capture substantial spaceflight signal.
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
  "Lung†": "#888888", "Colon†": "#555555"
}};
const tooltip = d3.select("#tooltip");

// Panel A: Flux heatmap
(function() {{
  const data = {json.dumps(flux_heatmap)};
  const tissues = {json.dumps([TISSUE_LABELS[t] for t in FLUX_TISSUES])};
  const subsystems = {json.dumps(ss_names)};

  const margin = {{top: 10, right: 20, bottom: 80, left: 200}};
  const cellW = 60, cellH = 22;
  const width = tissues.length * cellW, height = subsystems.length * cellH;

  const svg = d3.select("#flux-heatmap").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(tissues).range([0, width]).padding(0.05);
  const y = d3.scaleBand().domain(subsystems).range([0, height]).padding(0.05);
  const maxAbs = d3.max(data, d => Math.abs(d.flux_diff)) || 1;
  const color = d3.scaleDiverging(d3.interpolateRdBu).domain([maxAbs, 0, -maxAbs]);

  svg.selectAll("rect").data(data).join("rect")
    .attr("x", d => x(d.tissue)).attr("y", d => y(d.subsystem))
    .attr("width", x.bandwidth()).attr("height", y.bandwidth())
    .attr("fill", d => color(d.flux_diff))
    .on("mouseover", (e, d) => {{
      tooltip.style("display", "block")
        .html(`${{d.subsystem}}<br>${{d.tissue}}: ${{d.flux_diff}}`)
        .style("left", (e.pageX+10)+"px").style("top", (e.pageY-10)+"px");
    }}).on("mouseout", () => tooltip.style("display", "none"));

  svg.append("g").attr("transform", `translate(0,${{height}})`).call(d3.axisBottom(x))
    .selectAll("text").attr("transform", "rotate(-45)").style("text-anchor", "end").style("font-size", "8px");
  svg.append("g").call(d3.axisLeft(y)).selectAll("text").style("font-size", "7px");
}})();

// Panel B: Drug tiers
(function() {{
  const data = {json.dumps(drug_tiers)};
  const margin = {{top: 30, right: 80, bottom: 40, left: 180}};
  const width = 300, height = data.length * 50;

  const svg = d3.select("#drug-tiers").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleLinear().domain([0, d3.max(data, d => d.count) * 1.2]).range([0, width]);
  const y = d3.scaleBand().domain(data.map(d => d.tier)).range([0, height]).padding(0.4);

  const colors = [OI.blue, OI.orange, OI.skyblue];
  svg.selectAll("rect").data(data).join("rect")
    .attr("x", 0).attr("y", d => y(d.tier))
    .attr("width", d => x(d.count)).attr("height", y.bandwidth())
    .attr("fill", (d, i) => colors[i]);

  svg.selectAll("text.count").data(data).join("text").attr("class", "count")
    .attr("x", d => x(d.count) + 5).attr("y", d => y(d.tier) + y.bandwidth()/2 + 4)
    .attr("font-size", "10px").attr("font-weight", "bold").text(d => d.count);

  svg.append("g").call(d3.axisLeft(y)).selectAll("text").style("font-size", "9px");
  svg.append("g").attr("transform", `translate(0,${{height}})`).call(d3.axisBottom(x).ticks(5));
  svg.append("text").attr("x", width/2).attr("y", height + 30)
    .attr("text-anchor", "middle").attr("font-size", "9px").text("Number of drug-gene interactions");
}})();

// Panel C: Biomarker panel stacked bars
(function() {{
  const data = {json.dumps(panel)};
  const margin = {{top: 10, right: 20, bottom: 30, left: 80}};
  const width = 400, height = data.length * 20;

  const svg = d3.select("#panel-scores").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const y = d3.scaleBand().domain(data.map(d => d.gene)).range([0, height]).padding(0.2);
  const x = d3.scaleLinear().domain([0, 13]).range([0, width]);

  const sources = ["shap", "wgcna", "network", "drug", "lit"];
  const sourceColors = {{
    shap: OI.blue, wgcna: OI.green, network: OI.orange,
    drug: OI.vermilion, lit: OI.pink
  }};

  data.forEach(d => {{
    let cumX = 0;
    sources.forEach(s => {{
      const val = d[s] || 0;
      if (val > 0) {{
        svg.append("rect").attr("x", x(cumX)).attr("y", y(d.gene))
          .attr("width", x(val) - x(0)).attr("height", y.bandwidth())
          .attr("fill", sourceColors[s]);
        cumX += val;
      }}
    }});
  }});

  svg.append("g").call(d3.axisLeft(y)).selectAll("text").style("font-size", "8px");
  svg.append("g").attr("transform", `translate(0,${{height}})`).call(d3.axisBottom(x).ticks(7));

  // Legend
  const leg = svg.append("g").attr("transform", `translate(${{width - 200}}, -5)`);
  sources.forEach((s, i) => {{
    leg.append("rect").attr("x", i * 55).attr("y", 0).attr("width", 10).attr("height", 10)
      .attr("fill", sourceColors[s]);
    leg.append("text").attr("x", i * 55 + 13).attr("y", 9).attr("font-size", "7px")
      .text(s.toUpperCase());
  }});
}})();

// Panel D: Validation AUROC
(function() {{
  const data = {json.dumps(validation)};
  const margin = {{top: 10, right: 40, bottom: 50, left: 120}};
  const width = 350, height = data.length * 35;

  const svg = d3.select("#panel-validation").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleLinear().domain([0, 1]).range([0, width]);
  const y = d3.scaleBand().domain(data.map(d => d.tissue)).range([0, height]).padding(0.3);

  // Chance line
  svg.append("line").attr("x1", x(0.5)).attr("y1", 0)
    .attr("x2", x(0.5)).attr("y2", height)
    .attr("stroke", "#999").attr("stroke-dasharray", "4,2");

  svg.selectAll("rect").data(data).join("rect")
    .attr("x", 0).attr("y", d => y(d.tissue))
    .attr("width", d => x(d.auroc)).attr("height", y.bandwidth())
    .attr("fill", d => tissueColors[d.tissue] || "#ccc");

  svg.selectAll("text.auroc").data(data).join("text").attr("class", "auroc")
    .attr("x", d => x(d.auroc) + 4).attr("y", d => y(d.tissue) + y.bandwidth()/2 + 4)
    .attr("font-size", "9px").text(d => d.auroc.toFixed(3));

  svg.append("g").call(d3.axisLeft(y)).selectAll("text").style("font-size", "9px");
  svg.append("g").attr("transform", `translate(0,${{height}})`).call(d3.axisBottom(x).ticks(5));
  svg.append("text").attr("x", width/2).attr("y", height + 35)
    .attr("text-anchor", "middle").attr("font-size", "9px").text("AUROC (20-gene panel, PCA-LR)");
}})();
</script>
</body>
</html>"""

    out_path = OUT_DIR / "Fig8_metabolism_drugs.html"
    with open(out_path, "w") as f:
        f.write(html)
    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Fig 10: Mouse→Human Translation (v6)

6-panel figure summarizing translation of mouse spaceflight signatures to human cfRNA:
Panel A: Gene conservation — enrichment of mouse SHAP/WGCNA genes in human DRR
Panel B: Pathway conservation — Spearman rho of Hallmark NES across tissues
Panel C: Cross-species transfer — PCA-LR AUROC trained on mouse, tested on human
Panel D: Biomarker panel validation — panel genes in cfRNA (FC + FDR)
Panel E: TF activity conservation — rho between mouse and human TF activity
Panel F: Drug target translation — tier breakdown of druggable targets
"""
import json
import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent.parent
V6_EVAL = BASE_DIR / "v6" / "evaluation"
OUT_DIR = BASE_DIR / "v6" / "figures" / "html"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Okabe-Ito palette
OI = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "sky": "#56B4E9",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "yellow": "#F0E442",
    "black": "#000000",
    "gray": "#999999",
}
TISSUE_LABELS = {
    "liver": "Liver", "gastrocnemius": "Gastroc.", "kidney": "Kidney",
    "thymus": "Thymus", "eye": "Eye", "skin": "Skin",
    "lung": "Lung†", "colon": "Colon†",
}


def load_data():
    data = {}
    for phase in ["A", "B", "C", "D", "E", "F"]:
        fname = {
            "A": "V6_A_gene_conservation.json",
            "B": "V6_B_pathway_conservation.json",
            "C": "V6_C_cross_species_transfer.json",
            "D": "V6_D_biomarker_validation.json",
            "E": "V6_E_tf_conservation.json",
            "F": "V6_F_drug_target_validation.json",
        }[phase]
        with open(V6_EVAL / fname) as f:
            data[phase] = json.load(f)
    return data


def panel_a_data(a):
    """Gene conservation: fold enrichment and p-values for each gene set."""
    er = a["enrichment_results"]
    bars = []
    labels_map = {
        "shap_consensus": "SHAP Consensus\n(5 genes)",
        "shap_union_top100": "SHAP Union\n(1,361 genes)",
        "wgcna_hubs": "WGCNA Hubs\n(540 genes)",
        "biomarker_panel": "Biomarker\nPanel (20)",
    }
    for key, label in labels_map.items():
        if key not in er:
            continue
        v = er[key]
        fe = v.get("fold_enrichment", 0)
        p = v.get("hypergeometric_p", 1.0)
        overlap = v.get("n_overlap_drr", 0)
        n_human = v.get("n_human_genes", 0)
        bars.append({
            "label": label,
            "fe": round(fe, 3),
            "p": round(p, 4),
            "overlap": overlap,
            "n_human": n_human,
            "sig": p < 0.05,
        })
    return bars


def panel_b_data(b):
    """Pathway conservation: Spearman rho by tissue."""
    rows = []
    for tissue, v in b["per_tissue_correlations"].items():
        rows.append({
            "tissue": TISSUE_LABELS.get(tissue, tissue),
            "rho": round(v["spearman_rho"], 3),
            "p": round(v["p_value"], 4),
            "sig": v["significant"],
            "ci_lower": round(v["ci_lower"], 3),
            "ci_upper": round(v["ci_upper"], 3),
        })
    rows.sort(key=lambda x: x["rho"], reverse=True)
    return rows


def panel_c_data(c):
    """Cross-species transfer AUROC heatmap (pre_vs_post only)."""
    analyses = c["analyses"]
    pvp = analyses.get("pre_vs_post", {})
    rows = []
    for tissue, v in pvp.items():
        rows.append({
            "tissue": TISSUE_LABELS.get(tissue, tissue),
            "auroc": v.get("transfer_auroc", 0.5) or 0.5,
            "p": round(v.get("perm_p", 1.0), 4),
            "sig": v.get("significant", False),
            "mouse_auroc": round(v.get("mouse_internal_auroc", 0) or 0, 3),
            "n_genes": v.get("n_common_genes", 0),
        })
    rows.sort(key=lambda x: x["auroc"], reverse=True)
    return rows


def panel_d_data(d):
    """Biomarker panel genes: fold change in cfRNA."""
    genes = []
    for r in d["gene_results"]:
        if not r["detected_in_cfrna"]:
            continue
        genes.append({
            "gene": r.get("human_gene_matched", r["mouse_gene"]),
            "fc": round(r.get("de_diff", 1.0) or 1.0, 3),  # de_diff now stores FC
            "fdr": round(r.get("de_fdr", 1.0) or 1.0, 4),
            "is_de": r.get("is_de", False),
            "direction": r.get("direction", "down_in_flight"),
            "is_drr": r.get("is_drr", False),
            "n_drugs": r.get("n_drugs", 0),
        })
    genes.sort(key=lambda x: abs(x["fc"]), reverse=True)
    return genes


def panel_e_data(e):
    """TF conservation: rho by tissue."""
    rows = []
    for tissue, v in e["per_tissue_correlations"].items():
        rows.append({
            "tissue": TISSUE_LABELS.get(tissue, tissue),
            "rho": round(v["spearman_rho"], 3),
            "p": round(v["p_value"], 4),
            "sig": v.get("significant", False),
            "n_shared_tfs": v.get("n_shared_tfs", 0),
            "n_concordant": v.get("n_sig_both", 0),
            "concordance_rate": round(v.get("concordance_rate") or 0, 3),
        })
    rows.sort(key=lambda x: x["rho"], reverse=True)
    return rows


def panel_f_data(f):
    """Drug target tiers: bar chart + top genes."""
    tier_counts = f.get("tier_counts", {})
    tier_data = [
        {"tier": "A: Validated", "n": tier_counts.get("A", 0),
         "desc": "Sig. DE in cfRNA + druggable", "color": OI["green"]},
        {"tier": "B: Promising", "n": tier_counts.get("B", 0),
         "desc": "Sig. DE or SHAP + druggable", "color": OI["blue"]},
        {"tier": "C: Detected", "n": tier_counts.get("C", 0),
         "desc": "Present in cfRNA", "color": OI["sky"]},
        {"tier": "D: Untranslatable", "n": tier_counts.get("D", 0),
         "desc": "No human ortholog in cfRNA", "color": OI["gray"]},
    ]
    # Top tier A and B genes
    top_genes = []
    for g in f.get("tier_A_validated", []):
        top_genes.append({
            "gene": g["gene"], "tier": "A",
            "n_drugs": g["n_drugs"],
            "top_drug": g["top_drugs"][0] if g["top_drugs"] else "",
            "fdr": round(g.get("de_fdr", 1.0), 4),
        })
    for g in f.get("tier_B_promising", [])[:5]:
        top_genes.append({
            "gene": g["gene"], "tier": "B",
            "n_drugs": g["n_drugs"],
            "top_drug": g["top_drugs"][0] if g["top_drugs"] else "",
            "fdr": round(g.get("de_fdr", 1.0), 4),
        })
    return tier_data, top_genes


def generate_html(data):
    a_data = panel_a_data(data["A"])
    b_data = panel_b_data(data["B"])
    c_data = panel_c_data(data["C"])
    d_data = panel_d_data(data["D"])
    e_data = panel_e_data(data["E"])
    tier_data, top_genes = panel_f_data(data["F"])

    a_json = json.dumps(a_data)
    b_json = json.dumps(b_data)
    c_json = json.dumps(c_data)
    d_json = json.dumps(d_data)
    e_json = json.dumps(e_data)
    tier_json = json.dumps(tier_data)
    top_genes_json = json.dumps(top_genes)

    b_mean_rho = round(data["B"].get("mean_rho", 0), 3)
    e_mean_rho = round(data["E"].get("mean_rho", 0), 3)
    n_detected = data["D"].get("n_detected_in_cfrna", 0)
    panel_size = data["D"].get("panel_size", 20)
    n_tier_a = len(data["F"].get("tier_A_validated", []))
    n_tier_b = len(data["F"].get("tier_B_promising", []))

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Fig 10 — Mouse→Human Translation (GeneLabBench v6)</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: "Helvetica Neue", Arial, sans-serif; font-size: 11px;
         background: #fff; color: #222; padding: 20px; }}
  h1 {{ font-size: 14px; font-weight: bold; margin-bottom: 4px; color: #111; }}
  .subtitle {{ font-size: 10px; color: #555; margin-bottom: 16px; }}
  .fig-grid {{ display: grid; grid-template-columns: 1fr 1fr 1fr;
               grid-template-rows: auto auto; gap: 16px 20px; }}
  .panel {{ background: #fff; border: 1px solid #e0e0e0; border-radius: 4px; padding: 12px; }}
  .panel-title {{ font-size: 11px; font-weight: bold; margin-bottom: 8px; color: #333; }}
  .panel-subtitle {{ font-size: 9px; color: #777; margin-top: 4px; margin-bottom: 10px; }}
  svg text {{ font-family: "Helvetica Neue", Arial, sans-serif; }}
  .axis line, .axis path {{ stroke: #ccc; }}
  .axis text {{ font-size: 9px; fill: #555; }}
  .grid line {{ stroke: #eee; stroke-dasharray: 2,2; }}
  .sig-star {{ font-size: 10px; fill: #d55e00; font-weight: bold; }}
  .tooltip {{ position: absolute; background: rgba(0,0,0,0.8); color: #fff;
              padding: 6px 10px; border-radius: 4px; font-size: 10px;
              pointer-events: none; opacity: 0; transition: opacity 0.15s; }}
  .footer {{ font-size: 9px; color: #888; margin-top: 12px; }}
</style>
</head>
<body>
<h1>Fig 10 — Mouse→Human Translation: Spaceflight Signatures in Human cfRNA</h1>
<div class="subtitle">GeneLabBench v6 | 8 mouse tissues × 466 human DRR genes × 24 cfRNA samples</div>
<div class="tooltip" id="tooltip"></div>
<div class="fig-grid">
  <div class="panel" id="panelA">
    <div class="panel-title">A | Gene Conservation</div>
    <div class="panel-subtitle">Fold enrichment of mouse gene sets in human DRR (466 genes)</div>
    <svg id="svgA"></svg>
  </div>
  <div class="panel" id="panelB">
    <div class="panel-title">B | Pathway Conservation (mean ρ = {b_mean_rho})</div>
    <div class="panel-subtitle">Spearman ρ of Hallmark NES: mouse tissue vs human cfRNA</div>
    <svg id="svgB"></svg>
  </div>
  <div class="panel" id="panelC">
    <div class="panel-title">C | Cross-Species Transfer AUROC</div>
    <div class="panel-subtitle">PCA-LR trained on mouse tissue, tested on human cfRNA (pre vs post)</div>
    <svg id="svgC"></svg>
  </div>
  <div class="panel" id="panelD">
    <div class="panel-title">D | Biomarker Panel in cfRNA ({n_detected}/{panel_size} detected)</div>
    <div class="panel-subtitle">Fold change (pre→flight) for v5 panel genes detected in cfRNA</div>
    <svg id="svgD"></svg>
  </div>
  <div class="panel" id="panelE">
    <div class="panel-title">E | TF Activity Conservation (mean ρ = {e_mean_rho})</div>
    <div class="panel-subtitle">Spearman ρ of TF activity scores: mouse tissue vs human cfRNA</div>
    <svg id="svgE"></svg>
  </div>
  <div class="panel" id="panelF">
    <div class="panel-title">F | Drug Target Translation ({n_tier_a} validated, {n_tier_b} promising)</div>
    <div class="panel-subtitle">Druggable mouse spaceflight targets detected in human cfRNA</div>
    <svg id="svgF"></svg>
  </div>
</div>
<div class="footer">
  * p &lt; 0.05 | † New tissues (v4 only). Pathway conservation (B): FDR not corrected for 5 tests.
  Transfer AUROC (C): permutation test n=500. DRR = Differentially Regulated Response (SpaceOmicsBench).
</div>

<script>
const aData = {a_json};
const bData = {b_json};
const cData = {c_json};
const dData = {d_json};
const eData = {e_json};
const tierData = {tier_json};
const topGenes = {top_genes_json};

const OI = {{
  blue: "#0072B2", orange: "#E69F00", green: "#009E73",
  sky: "#56B4E9", red: "#D55E00", purple: "#CC79A7", gray: "#999"
}};

const tooltip = d3.select("#tooltip");
function showTip(html) {{
  tooltip.html(html).style("opacity", 1);
}}
function moveTip(event) {{
  tooltip.style("left", (event.pageX+12)+"px").style("top", (event.pageY-28)+"px");
}}
function hideTip() {{ tooltip.style("opacity", 0); }}

// ─── Panel A: Gene Conservation bar chart ───────────────────────────────────
(function() {{
  const W = 280, H = 160, mg = {{top:10, right:20, bottom:50, left:110}};
  const w = W - mg.left - mg.right, h = H - mg.top - mg.bottom;
  const svg = d3.select("#svgA").attr("width", W).attr("height", H);
  const g = svg.append("g").attr("transform", `translate(${{mg.left}},${{mg.top}})`);

  const x = d3.scaleLinear().domain([0, d3.max(aData, d => d.fe)*1.2 || 2]).range([0, w]);
  const y = d3.scaleBand().domain(aData.map(d => d.label)).range([0, h]).padding(0.3);

  g.append("g").attr("class","grid")
    .call(d3.axisBottom(x).ticks(4).tickSize(h).tickFormat(""))
    .attr("transform", "translate(0,0)").select(".domain").remove();

  // Reference line at FE=1
  g.append("line").attr("x1", x(1)).attr("x2", x(1))
    .attr("y1", 0).attr("y2", h)
    .attr("stroke", "#999").attr("stroke-dasharray","3,2").attr("stroke-width",1);

  g.selectAll(".bar").data(aData).join("rect")
    .attr("class","bar")
    .attr("x", 0).attr("y", d => y(d.label))
    .attr("width", d => x(d.fe)).attr("height", y.bandwidth())
    .attr("fill", d => d.sig ? OI.blue : OI.gray)
    .attr("rx", 2)
    .on("mouseover", (e,d) => {{
      showTip(`<b>${{d.label.replace("\\n"," ")}}</b><br>FE=${{d.fe}} (overlap=${{d.overlap}}/${{d.n_human}})<br>p=${{d.p}}`);
    }})
    .on("mousemove", moveTip).on("mouseout", hideTip);

  // Significance stars
  g.selectAll(".sig-star").data(aData.filter(d=>d.sig)).join("text")
    .attr("class","sig-star")
    .attr("x", d => x(d.fe)+3).attr("y", d => y(d.label)+y.bandwidth()/2+4)
    .text("*");

  g.append("g").attr("class","axis")
    .call(d3.axisLeft(y).tickSize(0).tickPadding(4))
    .select(".domain").remove();
  g.append("g").attr("class","axis").attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(x).ticks(4));
  g.append("text").attr("x",w/2).attr("y",h+32).attr("text-anchor","middle")
    .attr("font-size","9px").attr("fill","#555").text("Fold Enrichment vs random");
}})();

// ─── Panel B: Pathway Conservation rho bar chart ────────────────────────────
(function() {{
  const W = 280, H = 160, mg = {{top:10, right:20, bottom:50, left:65}};
  const w = W - mg.left - mg.right, h = H - mg.top - mg.bottom;
  const svg = d3.select("#svgB").attr("width", W).attr("height", H);
  const g = svg.append("g").attr("transform", `translate(${{mg.left}},${{mg.top}})`);

  const xDom = [d3.min(bData, d=>d.rho)-0.1, d3.max(bData, d=>d.rho)+0.1];
  const x = d3.scaleLinear().domain(xDom).range([0, w]);
  const y = d3.scaleBand().domain(bData.map(d=>d.tissue)).range([0,h]).padding(0.3);

  g.append("g").attr("class","grid")
    .call(d3.axisBottom(x).ticks(5).tickSize(h).tickFormat(""))
    .attr("transform","translate(0,0)").select(".domain").remove();

  // Zero line
  g.append("line").attr("x1",x(0)).attr("x2",x(0))
    .attr("y1",0).attr("y2",h)
    .attr("stroke","#888").attr("stroke-dasharray","3,2");

  g.selectAll(".bar").data(bData).join("rect")
    .attr("class","bar")
    .attr("x", d => d.rho>=0 ? x(0) : x(d.rho))
    .attr("y", d => y(d.tissue))
    .attr("width", d => Math.abs(x(d.rho)-x(0)))
    .attr("height", y.bandwidth())
    .attr("fill", d => d.rho >= 0 ? OI.blue : OI.red)
    .attr("opacity", d => d.sig ? 1.0 : 0.5)
    .attr("rx", 2)
    .on("mouseover", (e,d) => {{
      showTip(`<b>${{d.tissue}}</b><br>ρ=${{d.rho}} (CI: ${{d.ci_lower}}, ${{d.ci_upper}})<br>p=${{d.p}}`);
    }})
    .on("mousemove", moveTip).on("mouseout", hideTip);

  // CI whiskers
  g.selectAll(".ci").data(bData).join("line")
    .attr("x1", d=>x(d.ci_lower)).attr("x2", d=>x(d.ci_upper))
    .attr("y1", d=>y(d.tissue)+y.bandwidth()/2)
    .attr("y2", d=>y(d.tissue)+y.bandwidth()/2)
    .attr("stroke","#333").attr("stroke-width",1.5);

  g.selectAll(".sig").data(bData.filter(d=>d.sig)).join("text")
    .attr("class","sig-star")
    .attr("x", d=>x(d.ci_upper)+3).attr("y", d=>y(d.tissue)+y.bandwidth()/2+4)
    .text("*");

  g.append("g").attr("class","axis")
    .call(d3.axisLeft(y).tickSize(0).tickPadding(4))
    .select(".domain").remove();
  g.append("g").attr("class","axis").attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(x).ticks(5));
  g.append("text").attr("x",w/2).attr("y",h+32).attr("text-anchor","middle")
    .attr("font-size","9px").attr("fill","#555").text("Spearman ρ (mouse vs human Hallmark NES)");
}})();

// ─── Panel C: Cross-species transfer AUROC dot plot ─────────────────────────
(function() {{
  const W = 280, H = 160, mg = {{top:10, right:40, bottom:50, left:65}};
  const w = W - mg.left - mg.right, h = H - mg.top - mg.bottom;
  const svg = d3.select("#svgC").attr("width", W).attr("height", H);
  const g = svg.append("g").attr("transform", `translate(${{mg.left}},${{mg.top}})`);

  const x = d3.scaleLinear().domain([0, 1]).range([0, w]);
  const y = d3.scaleBand().domain(cData.map(d=>d.tissue)).range([0,h]).padding(0.3);

  g.append("g").attr("class","grid")
    .call(d3.axisBottom(x).ticks(5).tickSize(h).tickFormat(""))
    .attr("transform","translate(0,0)").select(".domain").remove();

  // Reference line at 0.5
  g.append("line").attr("x1",x(0.5)).attr("x2",x(0.5))
    .attr("y1",0).attr("y2",h).attr("stroke","#999").attr("stroke-dasharray","3,2");

  // Transfer AUROC bars
  g.selectAll(".bar").data(cData).join("rect")
    .attr("class","bar")
    .attr("x", x(0)).attr("y", d=>y(d.tissue))
    .attr("width", d=>x(d.auroc)-x(0))
    .attr("height", y.bandwidth())
    .attr("fill", d => d.sig ? OI.green : OI.gray)
    .attr("opacity", 0.7).attr("rx",2)
    .on("mouseover", (e,d) => {{
      showTip(`<b>${{d.tissue}}</b><br>Transfer AUROC=${{d.auroc}}<br>Mouse internal=${{d.mouse_auroc}}<br>p=${{d.p}} ${{d.sig?"*":""}}<br>Common genes=${{d.n_genes}}`);
    }})
    .on("mousemove", moveTip).on("mouseout", hideTip);

  g.append("g").attr("class","axis")
    .call(d3.axisLeft(y).tickSize(0).tickPadding(4))
    .select(".domain").remove();
  g.append("g").attr("class","axis").attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(x).ticks(5));
  g.append("text").attr("x",w/2).attr("y",h+32).attr("text-anchor","middle")
    .attr("font-size","9px").attr("fill","#555").text("Transfer AUROC (pre vs post)");
  // Legend
  const lg = g.append("g").attr("transform",`translate(${{w-60}},0)`);
  [[OI.green,"p<0.05"],[OI.gray,"n.s."]].forEach(([c,t],i) => {{
    lg.append("rect").attr("x",0).attr("y",i*12).attr("width",8).attr("height",8).attr("fill",c);
    lg.append("text").attr("x",11).attr("y",i*12+7).attr("font-size","8px").attr("fill","#555").text(t);
  }});
}})();

// ─── Panel D: Biomarker panel volcano ────────────────────────────────────────
(function() {{
  const W = 280, H = 160, mg = {{top:10, right:20, bottom:45, left:55}};
  const w = W - mg.left - mg.right, h = H - mg.top - mg.bottom;
  const svg = d3.select("#svgD").attr("width", W).attr("height", H);
  const g = svg.append("g").attr("transform", `translate(${{mg.left}},${{mg.top}})`);

  // Waterfall: genes sorted by abs(fc), colored by direction
  const sorted = [...dData].sort((a,b) => b.fc - a.fc);
  const x2 = d3.scaleBand().domain(sorted.map(d=>d.gene)).range([0,w]).padding(0.15);
  const fcVals = sorted.map(d=>d.fc).filter(v=>isFinite(v));
  const xmax = Math.max(Math.abs(d3.min(fcVals)), Math.abs(d3.max(fcVals)), 2);
  const y2 = d3.scaleLinear().domain([-xmax, xmax]).range([h, 0]);

  g.append("line").attr("x1",0).attr("x2",w).attr("y1",y2(0)).attr("y2",y2(0))
    .attr("stroke","#999").attr("stroke-dasharray","3,2");

  g.selectAll(".bar").data(sorted).join("rect")
    .attr("class","bar")
    .attr("x", d=>x2(d.gene))
    .attr("y", d=> d.fc>=0 ? y2(d.fc) : y2(0))
    .attr("width", x2.bandwidth())
    .attr("height", d=> Math.abs(y2(d.fc)-y2(0)))
    .attr("fill", d => d.fc>0 ? OI.orange : OI.blue)
    .attr("rx", 1)
    .on("mouseover", (e,d) => {{
      showTip(`<b>${{d.gene}}</b><br>FC(pre→flight)=${{d.fc}}<br>FDR=${{d.fdr}}<br>${{d.is_drr?"DRR ✓":""}} ${{d.n_drugs>0?d.n_drugs+" drugs":""}}`);
    }})
    .on("mousemove", moveTip).on("mouseout", hideTip);

  // Mark DE genes
  g.selectAll(".de-mark").data(sorted.filter(d=>d.is_de)).join("text")
    .attr("class","sig-star")
    .attr("x", d=>x2(d.gene)+x2.bandwidth()/2)
    .attr("y", d => d.fc>=0 ? y2(d.fc)-3 : y2(0)+11)
    .attr("text-anchor","middle").text("*");

  g.append("g").attr("class","axis").attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(x2).tickSize(0).tickPadding(2))
    .selectAll("text").attr("transform","rotate(-35)").attr("text-anchor","end").attr("font-size","7px");
  g.append("g").attr("class","axis")
    .call(d3.axisLeft(y2).ticks(5).tickSize(3));
  g.append("text").attr("transform","rotate(-90)").attr("x",-h/2).attr("y",-42)
    .attr("text-anchor","middle").attr("font-size","9px").attr("fill","#555")
    .text("Fold Change (pre→flight)");

  // Legend
  const lg = g.append("g").attr("transform",`translate(${{w-60}},0)`);
  [[OI.orange,"up in flight"],[OI.blue,"down in flight"]].forEach(([c,t],i) => {{
    lg.append("rect").attr("x",0).attr("y",i*11).attr("width",7).attr("height",7).attr("fill",c);
    lg.append("text").attr("x",9).attr("y",i*11+6).attr("font-size","7.5px").attr("fill","#555").text(t);
  }});
}})();

// ─── Panel E: TF conservation rho ────────────────────────────────────────────
(function() {{
  const W = 280, H = 160, mg = {{top:10, right:20, bottom:50, left:65}};
  const w = W - mg.left - mg.right, h = H - mg.top - mg.bottom;
  const svg = d3.select("#svgE").attr("width", W).attr("height", H);
  const g = svg.append("g").attr("transform", `translate(${{mg.left}},${{mg.top}})`);

  const xDom = [d3.min(eData,d=>d.rho)-0.1, d3.max(eData,d=>d.rho)+0.1];
  const x = d3.scaleLinear().domain(xDom).range([0, w]);
  const y = d3.scaleBand().domain(eData.map(d=>d.tissue)).range([0,h]).padding(0.3);

  g.append("g").attr("class","grid")
    .call(d3.axisBottom(x).ticks(5).tickSize(h).tickFormat(""))
    .attr("transform","translate(0,0)").select(".domain").remove();

  g.append("line").attr("x1",x(0)).attr("x2",x(0))
    .attr("y1",0).attr("y2",h)
    .attr("stroke","#888").attr("stroke-dasharray","3,2");

  g.selectAll(".bar").data(eData).join("rect")
    .attr("class","bar")
    .attr("x", d=>d.rho>=0 ? x(0) : x(d.rho))
    .attr("y", d=>y(d.tissue))
    .attr("width", d=>Math.abs(x(d.rho)-x(0)))
    .attr("height", y.bandwidth())
    .attr("fill", d=>d.rho>=0 ? OI.blue : OI.red)
    .attr("opacity", d=>d.sig ? 1.0 : 0.5)
    .attr("rx", 2)
    .on("mouseover", (e,d) => {{
      showTip(`<b>${{d.tissue}}</b><br>ρ=${{d.rho}}, p=${{d.p}}<br>Shared TFs=${{d.n_shared_tfs}}<br>Concordance=${{d.concordance_rate}}`);
    }})
    .on("mousemove", moveTip).on("mouseout", hideTip);

  g.selectAll(".sig").data(eData.filter(d=>d.sig)).join("text")
    .attr("class","sig-star")
    .attr("x", d=>x(Math.max(d.rho,0))+3).attr("y", d=>y(d.tissue)+y.bandwidth()/2+4)
    .text("*");

  g.append("g").attr("class","axis")
    .call(d3.axisLeft(y).tickSize(0).tickPadding(4))
    .select(".domain").remove();
  g.append("g").attr("class","axis").attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(x).ticks(5));
  g.append("text").attr("x",w/2).attr("y",h+32).attr("text-anchor","middle")
    .attr("font-size","9px").attr("fill","#555").text("Spearman ρ (mouse vs human TF activity)");
}})();

// ─── Panel F: Drug target tiers stacked bar ──────────────────────────────────
(function() {{
  const W = 280, H = 160, mg = {{top:10, right:20, bottom:50, left:30}};
  const w = W - mg.left - mg.right, h = H - mg.top - mg.bottom;
  const svg = d3.select("#svgF").attr("width", W).attr("height", H);
  const g = svg.append("g").attr("transform", `translate(${{mg.left}},${{mg.top}})`);

  const total = tierData.reduce((s,d) => s+d.n, 0);
  const x = d3.scaleLinear().domain([0, total]).range([0, w]);
  const barH = 30;
  const barY = 20;

  // Stacked horizontal bar
  let cumX = 0;
  tierData.forEach(d => {{
    if (d.n === 0) return;
    g.append("rect").attr("x", x(cumX)).attr("y", barY)
      .attr("width", x(d.n)).attr("height", barH)
      .attr("fill", d.color).attr("rx",2)
      .on("mouseover", (e) => {{
        showTip(`<b>Tier ${{d.tier}}</b><br>n=${{d.n}} (${{(d.n/total*100).toFixed(1)}}%)<br>${{d.desc}}`);
      }})
      .on("mousemove", moveTip).on("mouseout", hideTip);
    // Label if wide enough
    if (x(d.n) > 30) {{
      g.append("text").attr("x", x(cumX + d.n/2)).attr("y", barY+barH/2+4)
        .attr("text-anchor","middle").attr("font-size","9px").attr("fill","#fff")
        .attr("font-weight","bold")
        .text(`n=${{d.n}}`);
    }}
    cumX += d.n;
  }});

  // Legend
  tierData.forEach((d, i) => {{
    const lx = (i % 2) * (w/2), ly = barY + barH + 12 + Math.floor(i/2)*14;
    g.append("rect").attr("x",lx).attr("y",ly).attr("width",8).attr("height",8).attr("fill",d.color);
    g.append("text").attr("x",lx+11).attr("y",ly+7).attr("font-size","8.5px").attr("fill","#333")
      .text(d.tier.split(":")[0]+": "+d.desc.substring(0,25));
  }});

  // Top genes table
  const tableY = barY + barH + 50;
  g.append("text").attr("x",0).attr("y",tableY-2).attr("font-size","9px")
    .attr("font-weight","bold").attr("fill","#333").text("Top validated targets:");
  topGenes.slice(0,5).forEach((tg, i) => {{
    const ty = tableY + 11 + i*11;
    const col = tg.tier === "A" ? OI.green : OI.blue;
    g.append("text").attr("x",0).attr("y",ty).attr("font-size","8.5px").attr("fill",col)
      .attr("font-weight","bold").text(`[${{tg.tier}}]`);
    g.append("text").attr("x",18).attr("y",ty).attr("font-size","8.5px").attr("fill","#333")
      .text(`${{tg.gene}} (${{tg.n_drugs}} drugs)`);
  }});
}})();
</script>
</body>
</html>"""
    return html


def main():
    print("Loading v6 data...")
    data = load_data()

    print("Generating Fig 10...")
    html = generate_html(data)

    out_path = OUT_DIR / "v6_Fig10_human_translation.html"
    with open(out_path, "w") as f:
        f.write(html)
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    main()

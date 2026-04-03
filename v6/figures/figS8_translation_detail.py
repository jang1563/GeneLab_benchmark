#!/usr/bin/env python3
"""
Fig S8: Mouse→Human Translation — Supplementary Details (v6)

Panel A: Per-tissue gene set overlap heatmap (tissues × gene sets)
Panel B: Top concordant/discordant Hallmark pathways (cross-tissue)
Panel C: Transfer AUROC: pre_vs_post vs pre_vs_post_recovery
Panel D: Biomarker panel gene-level detail (FC + DRR + drugs)
Panel E: TF concordance rates by tissue
Panel F: Full drug target list (Tier A + B)
"""
import json
import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent.parent
V6_EVAL = BASE_DIR / "v6" / "evaluation"
OUT_DIR = BASE_DIR / "v6" / "figures" / "html"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OI = {
    "blue": "#0072B2", "orange": "#E69F00", "green": "#009E73",
    "sky": "#56B4E9", "red": "#D55E00", "purple": "#CC79A7",
    "yellow": "#F0E442", "gray": "#999999",
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
    """Per-tissue enrichment table for heatmap."""
    tissues = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin", "lung", "colon"]
    per_tissue = a.get("per_tissue_enrichment", {})
    cells = []
    for tissue in tissues:
        if tissue not in per_tissue:
            continue
        v = per_tissue[tissue]
        cells.append({
            "tissue": TISSUE_LABELS.get(tissue, tissue),
            "n_mouse": v.get("n_mouse_genes", 0),
            "n_human": v.get("n_human_genes", 0),
            "n_overlap": v.get("n_overlap_drr", 0),
            "fe": round(v.get("fold_enrichment", 0), 3),
            "p": round(v.get("hypergeometric_p", 1.0), 4),
            "sig": v.get("significant", False),
            "n_unmapped": v.get("n_unmapped", 0),
        })
    cells.sort(key=lambda x: x["fe"], reverse=True)
    return cells


def panel_b_data(b):
    """Top concordant and discordant Hallmark pathways across tissues."""
    concordant_counts = {}
    discordant_counts = {}

    for tissue, v in b["per_tissue_correlations"].items():
        for pw_info in v.get("top_concordant", [])[:5]:
            pw = pw_info.get("pathway", "")
            concordant_counts[pw] = concordant_counts.get(pw, 0) + 1
        for pw_info in v.get("top_discordant", [])[:5]:
            pw = pw_info.get("pathway", "")
            discordant_counts[pw] = discordant_counts.get(pw, 0) + 1

    # Get all pathways with their per-tissue NES
    all_pathways = set(list(concordant_counts.keys())[:5] + list(discordant_counts.keys())[:5])

    # Per-tissue NES for top pathways
    pathway_tissue_nes = {}
    for tissue, v in b["per_tissue_correlations"].items():
        t_label = TISSUE_LABELS.get(tissue, tissue)
        for pw_info in (v.get("top_concordant", []) + v.get("top_discordant", [])):
            pw = pw_info.get("pathway", "")
            if pw in all_pathways:
                if pw not in pathway_tissue_nes:
                    pathway_tissue_nes[pw] = {}
                m_nes = pw_info.get("mouse_nes", 0)
                h_nes = pw_info.get("human_nes", 0)
                pathway_tissue_nes[pw][t_label] = {
                    "mouse_nes": round(m_nes, 3),
                    "human_nes": round(h_nes, 3),
                }

    top_concordant = sorted(concordant_counts.items(), key=lambda x: -x[1])[:6]
    top_discordant = sorted(discordant_counts.items(), key=lambda x: -x[1])[:4]

    return top_concordant, top_discordant, pathway_tissue_nes


def panel_c_data(c):
    """Transfer AUROC comparison: pre_vs_post vs pre_vs_post_recovery."""
    pvp = c["analyses"].get("pre_vs_post", {})
    pvpr = c["analyses"].get("pre_vs_post_recovery", {})

    rows = []
    all_tissues = set(list(pvp.keys()) + list(pvpr.keys()))
    for tissue in all_tissues:
        a1 = pvp.get(tissue, {})
        a2 = pvpr.get(tissue, {})
        rows.append({
            "tissue": TISSUE_LABELS.get(tissue, tissue),
            "pvp": a1.get("transfer_auroc", 0.5) or 0.5,
            "pvpr": a2.get("transfer_auroc", 0.5) or 0.5,
            "pvp_p": round(a1.get("perm_p", 1.0), 4),
            "pvpr_p": round(a2.get("perm_p", 1.0), 4),
            "n_genes": a1.get("n_common_genes", 0),
            "mouse_auroc": round(a1.get("mouse_internal_auroc", 0) or 0, 3),
        })
    rows.sort(key=lambda x: x["pvp"], reverse=True)
    return rows


def panel_d_data(d):
    """Full biomarker panel gene detail for supplementary table."""
    rows = []
    for r in d["gene_results"]:
        rows.append({
            "rank": r.get("rank", 0),
            "mouse_gene": r.get("mouse_gene", ""),
            "human_gene": r.get("human_gene_matched") or "-",
            "score": r.get("panel_score", 0),
            "detected": r.get("detected_in_cfrna", False),
            "fc": round(r.get("de_diff", 1.0) or 1.0, 3),
            "fdr": round(r.get("de_fdr", 1.0) or 1.0, 4),
            "is_de": r.get("is_de", False),
            "is_drr": r.get("is_drr", False),
            "direction": r.get("direction", "-"),
            "n_drugs": r.get("n_drugs", 0),
            "n_tissues": len(r.get("tissues", [])),
        })
    return rows


def panel_e_data(e):
    """TF conservation detail: n_sig_both, concordance_rate by tissue."""
    rows = []
    for tissue, v in e["per_tissue_correlations"].items():
        rows.append({
            "tissue": TISSUE_LABELS.get(tissue, tissue),
            "rho": round(v["spearman_rho"], 3),
            "p": round(v["p_value"], 4),
            "n_shared": v.get("n_shared_tfs", 0),
            "n_sig_both": v.get("n_sig_both", 0),
            "concordance": round(v.get("concordance_rate") or 0, 3),
            "sig": v.get("significant", False),
        })
    rows.sort(key=lambda x: x["rho"], reverse=True)
    return rows


def panel_f_data(f):
    """Full tier A + B drug target list."""
    genes = []
    for tier, key in [("A", "tier_A_validated"), ("B", "tier_B_promising")]:
        for g in f.get(key, []):
            genes.append({
                "tier": tier,
                "gene": g["gene"],
                "n_drugs": g["n_drugs"],
                "top_drugs": ", ".join(g.get("top_drugs", [])[:3]),
                "fda": g.get("has_fda_approved", False),
                "fdr": round(g.get("de_fdr", 1.0), 4),
                "direction": g.get("direction", "-"),
            })
    return genes


def generate_html(data):
    a_cells = panel_a_data(data["A"])
    top_conc, top_disc, pw_nes = panel_b_data(data["B"])
    c_rows = panel_c_data(data["C"])
    d_rows = panel_d_data(data["D"])
    e_rows = panel_e_data(data["E"])
    f_genes = panel_f_data(data["F"])

    # For panel B: pathway correlation heatmap data
    tissue_list = list(data["B"]["per_tissue_correlations"].keys())
    # Create NES comparison table for top pathways
    top_pws = [p for p, _ in top_conc[:4]] + [p for p, _ in top_disc[:3]]
    pw_cells = []
    for pw in top_pws:
        for tissue in tissue_list:
            t_label = TISSUE_LABELS.get(tissue, tissue)
            nes_info = pw_nes.get(pw, {}).get(t_label, {})
            if nes_info:
                pw_cells.append({
                    "pathway": pw.replace("_", " ").title()[:28],
                    "tissue": t_label,
                    "mouse_nes": nes_info.get("mouse_nes", 0),
                    "human_nes": nes_info.get("human_nes", 0),
                    "concordant": (nes_info.get("mouse_nes", 0) > 0) == (nes_info.get("human_nes", 0) > 0),
                })

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Fig S8 — Mouse→Human Translation: Supplementary Details (GeneLabBench v6)</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: "Helvetica Neue", Arial, sans-serif; font-size: 11px;
         background: #fff; color: #222; padding: 20px; }}
  h1 {{ font-size: 13px; font-weight: bold; margin-bottom: 4px; }}
  .subtitle {{ font-size: 10px; color: #555; margin-bottom: 14px; }}
  .fig-grid {{ display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 14px 18px; }}
  .panel {{ background: #fff; border: 1px solid #e0e0e0; border-radius: 4px; padding: 10px; }}
  .panel-title {{ font-size: 11px; font-weight: bold; margin-bottom: 6px; color: #333; }}
  .panel-subtitle {{ font-size: 9px; color: #777; margin-top: 2px; margin-bottom: 8px; }}
  svg text {{ font-family: "Helvetica Neue", Arial, sans-serif; }}
  .axis line, .axis path {{ stroke: #ccc; }}
  .axis text {{ font-size: 8.5px; fill: #555; }}
  .grid line {{ stroke: #eee; stroke-dasharray: 2,2; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 9px; }}
  th {{ background: #f5f5f5; font-weight: bold; padding: 3px 5px; border: 1px solid #ddd;
        text-align: left; position: sticky; top: 0; }}
  td {{ padding: 2px 5px; border: 1px solid #eee; }}
  tr:hover {{ background: #f9f9f9; }}
  .tier-a {{ color: #009E73; font-weight: bold; }}
  .tier-b {{ color: #0072B2; font-weight: bold; }}
  .detected {{ color: #009E73; }}
  .not-detected {{ color: #999; }}
  .sig {{ color: #D55E00; font-weight: bold; }}
  .tooltip {{ position: absolute; background: rgba(0,0,0,0.8); color: #fff;
              padding: 5px 8px; border-radius: 3px; font-size: 9px;
              pointer-events: none; opacity: 0; transition: opacity 0.15s; }}
  .footer {{ font-size: 9px; color: #888; margin-top: 10px; }}
  .scroll-box {{ max-height: 200px; overflow-y: auto; }}
</style>
</head>
<body>
<h1>Fig S8 — Mouse→Human Translation: Supplementary Details</h1>
<div class="subtitle">GeneLabBench v6 | Extended analysis of spaceflight signature conservation</div>
<div class="tooltip" id="tooltip"></div>
<div class="fig-grid">
  <div class="panel" id="panelA">
    <div class="panel-title">S8A | Per-Tissue Gene Overlap with Human DRR</div>
    <div class="panel-subtitle">SHAP top-100 genes mapped to human orthologs, overlap with 466 DRR genes</div>
    <svg id="svgA"></svg>
  </div>
  <div class="panel" id="panelB">
    <div class="panel-title">S8B | Top Concordant Hallmark Pathways</div>
    <div class="panel-subtitle">Pathways most consistently co-dysregulated in mouse + human</div>
    <svg id="svgB"></svg>
  </div>
  <div class="panel" id="panelC">
    <div class="panel-title">S8C | Transfer AUROC: Two Comparisons</div>
    <div class="panel-subtitle">Pre vs Post (black) and Pre vs Post+Recovery (gray)</div>
    <svg id="svgC"></svg>
  </div>
  <div class="panel" id="panelD">
    <div class="panel-title">S8D | Biomarker Panel in Human cfRNA</div>
    <div class="panel-subtitle">All 20 panel genes with FC (pre→flight) and detection status</div>
    <div class="scroll-box">
    <table id="tableD">
      <thead><tr>
        <th>Gene</th><th>Human</th><th>Score</th><th>Tissues</th><th>FC</th><th>FDR</th><th>Drugs</th>
      </tr></thead>
      <tbody id="tbodyD"></tbody>
    </table>
    </div>
  </div>
  <div class="panel" id="panelE">
    <div class="panel-title">S8E | TF Activity Concordance by Tissue</div>
    <div class="panel-subtitle">Number of TFs significant in both mouse and human</div>
    <svg id="svgE"></svg>
  </div>
  <div class="panel" id="panelF">
    <div class="panel-title">S8F | Drug Target Summary (Tier A + B)</div>
    <div class="panel-subtitle">Mouse spaceflight drug targets translated to human cfRNA</div>
    <div class="scroll-box">
    <table id="tableF">
      <thead><tr>
        <th>Tier</th><th>Gene</th><th>Drugs</th><th>Top Drug</th><th>FDR</th><th>Dir.</th>
      </tr></thead>
      <tbody id="tbodyF"></tbody>
    </table>
    </div>
  </div>
</div>
<div class="footer">
  DRR = Differentially Regulated Response (SpaceOmicsBench cfRNA). TF activity from decoupler-py + CollecTRI.
  * p &lt; 0.05. FC = fold change pre→flight (EdgeR). † New tissues added in v4.
</div>

<script>
const aCells = {json.dumps(a_cells)};
const bConcordant = {json.dumps(top_conc)};
const bDiscordant = {json.dumps(top_disc)};
const pwCells = {json.dumps(pw_cells)};
const cRows = {json.dumps(c_rows)};
const dRows = {json.dumps(d_rows)};
const eRows = {json.dumps(e_rows)};
const fGenes = {json.dumps(f_genes)};

const OI = {{
  blue: "#0072B2", orange: "#E69F00", green: "#009E73",
  sky: "#56B4E9", red: "#D55E00", gray: "#999"
}};
const tooltip = d3.select("#tooltip");
function showTip(h) {{ tooltip.html(h).style("opacity",1); }}
function moveTip(e) {{ tooltip.style("left",(e.pageX+12)+"px").style("top",(e.pageY-28)+"px"); }}
function hideTip() {{ tooltip.style("opacity",0); }}

// ─── Panel A: Per-tissue overlap bar ────────────────────────────────────────
(function() {{
  const W=280, H=200, mg={{top:10,right:50,bottom:20,left:65}};
  const w=W-mg.left-mg.right, h=H-mg.top-mg.bottom;
  const svg=d3.select("#svgA").attr("width",W).attr("height",H);
  const g=svg.append("g").attr("transform",`translate(${{mg.left}},${{mg.top}})`);

  const x=d3.scaleLinear().domain([0,d3.max(aCells,d=>d.fe)*1.2||3]).range([0,w]);
  const y=d3.scaleBand().domain(aCells.map(d=>d.tissue)).range([0,h]).padding(0.25);

  g.append("line").attr("x1",x(1)).attr("x2",x(1)).attr("y1",0).attr("y2",h)
    .attr("stroke","#999").attr("stroke-dasharray","3,2");

  g.selectAll(".bar").data(aCells).join("rect")
    .attr("class","bar").attr("x",0).attr("y",d=>y(d.tissue))
    .attr("width",d=>x(d.fe)).attr("height",y.bandwidth())
    .attr("fill",d=>d.sig?OI.blue:OI.gray).attr("rx",2)
    .on("mouseover",(e,d)=>{{showTip(`<b>${{d.tissue}}</b><br>FE=${{d.fe}} (overlap=${{d.n_overlap}}/${{d.n_human}})<br>p=${{d.p}}<br>Unmapped=${{d.n_unmapped}}`); }})
    .on("mousemove",moveTip).on("mouseout",hideTip);

  // Overlap labels
  g.selectAll(".olabel").data(aCells).join("text")
    .attr("x",d=>x(d.fe)+4).attr("y",d=>y(d.tissue)+y.bandwidth()/2+3)
    .attr("font-size","8px").attr("fill","#555")
    .text(d=>d.n_overlap>0?`n=${{d.n_overlap}}`:"");

  g.append("g").attr("class","axis")
    .call(d3.axisLeft(y).tickSize(0).tickPadding(4)).select(".domain").remove();
  g.append("g").attr("class","axis").attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(x).ticks(4));
  g.append("text").attr("x",w/2).attr("y",h+16).attr("text-anchor","middle")
    .attr("font-size","8.5px").attr("fill","#555").text("Fold Enrichment in DRR");
}})();

// ─── Panel B: Top concordant pathways NES scatter ────────────────────────────
(function() {{
  const W=280, H=200, mg={{top:10,right:20,bottom:30,left:100}};
  const w=W-mg.left-mg.right, h=H-mg.top-mg.bottom;
  const svg=d3.select("#svgB").attr("width",W).attr("height",H);
  const g=svg.append("g").attr("transform",`translate(${{mg.left}},${{mg.top}})`);

  // Bar chart: count of tissues where pathway is concordant
  const allPws = [...bConcordant, ...bDiscordant];
  if (allPws.length === 0) {{
    g.append("text").attr("x",w/2).attr("y",h/2).attr("text-anchor","middle")
      .attr("font-size","10px").attr("fill","#999").text("No pathway data");
    return;
  }}

  const maxCount = d3.max(allPws, d=>d[1]);
  const x=d3.scaleLinear().domain([0,maxCount]).range([0,w]);
  const y=d3.scaleBand().domain(allPws.map(d=>d[0].replace(/_/g," ").substring(0,22)))
    .range([0,h]).padding(0.25);

  g.selectAll(".bar").data(allPws).join("rect")
    .attr("class","bar").attr("x",0)
    .attr("y",d=>y(d[0].replace(/_/g," ").substring(0,22)))
    .attr("width",d=>x(d[1])).attr("height",y.bandwidth())
    .attr("fill",(d,i)=>i<bConcordant.length?OI.blue:OI.orange)
    .attr("rx",2)
    .on("mouseover",(e,d)=>{{showTip(`<b>${{d[0]}}</b><br>${{d[1]}} tissues concordant`); }})
    .on("mousemove",moveTip).on("mouseout",hideTip);

  g.selectAll(".cnt").data(allPws).join("text")
    .attr("x",d=>x(d[1])+3).attr("y",d=>y(d[0].replace(/_/g," ").substring(0,22))+y.bandwidth()/2+3)
    .attr("font-size","8px").attr("fill","#555").text(d=>d[1]);

  g.append("g").attr("class","axis")
    .call(d3.axisLeft(y).tickSize(0).tickPadding(3)).select(".domain").remove();
  g.append("g").attr("class","axis").attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(x).ticks(maxCount).tickFormat(d3.format("d")));
  g.append("text").attr("x",w/2).attr("y",h+22).attr("text-anchor","middle")
    .attr("font-size","8.5px").attr("fill","#555").text("No. tissues (concordant=blue, discordant=orange)");
}})();

// ─── Panel C: Transfer AUROC grouped bar ────────────────────────────────────
(function() {{
  const W=280, H=200, mg={{top:10,right:20,bottom:30,left:65}};
  const w=W-mg.left-mg.right, h=H-mg.top-mg.bottom;
  const svg=d3.select("#svgC").attr("width",W).attr("height",H);
  const g=svg.append("g").attr("transform",`translate(${{mg.left}},${{mg.top}})`);

  const x0=d3.scaleBand().domain(cRows.map(d=>d.tissue)).range([0,w]).padding(0.2);
  const x1=d3.scaleBand().domain(["pvp","pvpr"]).range([0,x0.bandwidth()]).padding(0.05);
  const y=d3.scaleLinear().domain([0,1]).range([h,0]);

  // Reference 0.5 line
  g.append("line").attr("x1",0).attr("x2",w).attr("y1",y(0.5)).attr("y2",y(0.5))
    .attr("stroke","#999").attr("stroke-dasharray","3,2");

  cRows.forEach(d => {{
    const xBase = x0(d.tissue);
    [{{"key":"pvp","val":d.pvp,"p":d.pvp_p,"col":OI.blue}},
     {{"key":"pvpr","val":d.pvpr,"p":d.pvpr_p,"col":OI.sky}}].forEach(dd => {{
      g.append("rect")
        .attr("x",xBase+x1(dd.key)).attr("y",y(dd.val))
        .attr("width",x1.bandwidth()).attr("height",Math.max(0,y(0)-y(dd.val)))
        .attr("fill",dd.col).attr("rx",1)
        .on("mouseover",(e)=>{{showTip(`<b>${{d.tissue}}</b><br>AUROC=${{dd.val.toFixed(3)}}<br>p=${{dd.p}}<br>Mouse=${{d.mouse_auroc}}`); }})
        .on("mousemove",moveTip).on("mouseout",hideTip);
    }});
  }});

  g.append("g").attr("class","axis").attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(x0).tickSize(0).tickPadding(3))
    .selectAll("text").attr("transform","rotate(-30)").attr("text-anchor","end").attr("font-size","8px");
  g.append("g").attr("class","axis")
    .call(d3.axisLeft(y).ticks(5));

  const lg=g.append("g").attr("transform",`translate(${{w-70}},2)`);
  [[OI.blue,"Pre vs Post"],[OI.sky,"Pre vs Post+Rec"]].forEach(([c,t],i)=>{{
    lg.append("rect").attr("x",0).attr("y",i*12).attr("width",8).attr("height",8).attr("fill",c);
    lg.append("text").attr("x",10).attr("y",i*12+7).attr("font-size","8px").attr("fill","#555").text(t);
  }});
}})();

// ─── Panel D: Biomarker table ────────────────────────────────────────────────
(function() {{
  const tbody = d3.select("#tbodyD");
  dRows.forEach(r => {{
    const tr = tbody.append("tr");
    tr.append("td").style("font-weight","bold")
      .style("color", r.detected ? "#333" : "#bbb")
      .text(r.mouse_gene);
    tr.append("td").attr("class", r.detected ? "detected" : "not-detected")
      .text(r.detected ? r.human_gene : "—");
    tr.append("td").text(r.score);
    tr.append("td").text(r.n_tissues);
    tr.append("td").style("color", r.fc > 1.5 ? OI.orange : r.fc < 0.67 ? OI.blue : "#555")
      .text(r.detected ? r.fc.toFixed(2) : "—");
    tr.append("td").attr("class", r.is_de ? "sig" : "")
      .text(r.detected ? r.fdr.toFixed(4) : "—");
    tr.append("td").text(r.n_drugs || "—");
  }});
}})();

// ─── Panel E: TF concordance grouped bar ────────────────────────────────────
(function() {{
  const W=280, H=200, mg={{top:10,right:20,bottom:20,left:65}};
  const w=W-mg.left-mg.right, h=H-mg.top-mg.bottom;
  const svg=d3.select("#svgE").attr("width",W).attr("height",H);
  const g=svg.append("g").attr("transform",`translate(${{mg.left}},${{mg.top}})`);

  const x=d3.scaleBand().domain(eRows.map(d=>d.tissue)).range([0,w]).padding(0.25);
  const x1=d3.scaleBand().domain(["rho","conc"]).range([0,x.bandwidth()]).padding(0.1);
  const maxRho = Math.max(Math.abs(d3.min(eRows,d=>d.rho)), d3.max(eRows,d=>d.rho));
  const y=d3.scaleLinear().domain([-maxRho*1.2, maxRho*1.2]).range([h,0]);
  const yRight=d3.scaleLinear().domain([0,1]).range([h,0]);

  g.append("line").attr("x1",0).attr("x2",w).attr("y1",y(0)).attr("y2",y(0))
    .attr("stroke","#999").attr("stroke-dasharray","3,2");

  eRows.forEach(d => {{
    const xBase=x(d.tissue);
    // Rho bar
    g.append("rect")
      .attr("x",xBase+x1("rho")).attr("y",d.rho>=0?y(d.rho):y(0))
      .attr("width",x1.bandwidth()).attr("height",Math.abs(y(d.rho)-y(0)))
      .attr("fill",d.rho>=0?OI.blue:OI.red)
      .attr("opacity",d.sig?1:0.5).attr("rx",1)
      .on("mouseover",(e)=>{{showTip(`<b>${{d.tissue}}</b><br>ρ=${{d.rho}}, p=${{d.p}}<br>Shared TFs=${{d.n_shared}}<br>Sig in both=${{d.n_sig_both}}<br>Concordance=${{d.concordance}}`); }})
      .on("mousemove",moveTip).on("mouseout",hideTip);
    // Concordance dot (right axis scale)
    const cy=yRight(d.concordance);
    g.append("circle")
      .attr("cx",xBase+x.bandwidth()/2).attr("cy",cy).attr("r",3)
      .attr("fill",OI.orange).attr("stroke","#fff").attr("stroke-width",0.5)
      .on("mouseover",(e)=>{{showTip(`<b>${{d.tissue}}</b><br>Concordance rate=${{d.concordance}}`); }})
      .on("mousemove",moveTip).on("mouseout",hideTip);
  }});

  g.append("g").attr("class","axis").attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(x).tickSize(0).tickPadding(3))
    .selectAll("text").attr("transform","rotate(-25)").attr("text-anchor","end").attr("font-size","8px");
  g.append("g").attr("class","axis").call(d3.axisLeft(y).ticks(5));

  // Legend
  const lg=g.append("g").attr("transform",`translate(${{w-70}},2)`);
  [[OI.blue,"Spearman ρ"],[OI.orange,"Concordance"]].forEach(([c,t],i)=>{{
    lg.append("rect").attr("x",0).attr("y",i*12).attr("width",8).attr("height",8).attr("fill",c);
    lg.append("text").attr("x",10).attr("y",i*12+7).attr("font-size","8px").attr("fill","#555").text(t);
  }});
}})();

// ─── Panel F: Drug target table ──────────────────────────────────────────────
(function() {{
  const tbody = d3.select("#tbodyF");
  fGenes.forEach(g => {{
    const tr = tbody.append("tr");
    tr.append("td").attr("class",`tier-${{g.tier.toLowerCase()}}`).text(g.tier);
    tr.append("td").style("font-weight","bold").text(g.gene);
    tr.append("td").text(g.n_drugs);
    tr.append("td").style("font-size","8px").style("color","#555").text(g.top_drugs || "—");
    tr.append("td").attr("class",g.fdr<0.05?"sig":"").text(g.fdr);
    tr.append("td").style("color",g.direction==="up_in_flight"?OI.orange:OI.blue)
      .text(g.direction==="up_in_flight"?"↑":"↓");
  }});
}})();
</script>
</body>
</html>"""
    return html


def main():
    print("Loading v6 data...")
    data = load_data()

    print("Generating Fig S8...")
    html = generate_html(data)

    out_path = OUT_DIR / "v6_FigS8_translation_detail.html"
    with open(out_path, "w") as f:
        f.write(html)
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    main()

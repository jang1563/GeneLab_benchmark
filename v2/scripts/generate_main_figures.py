#!/usr/bin/env python3
"""
generate_main_figures.py
Generates 3 publication-ready integrated main figures for GeneLab v2 paper.

Output:
  v2/figures/Fig1_temporal.html       — Mouse temporal dynamics (T1/T2/T3)
  v2/figures/Fig2_crossspecies.html   — Cross-species conservation (E1/E2/E3)
  v2/figures/Fig3_pbmc_celltype.html  — PBMC cell-type resolution (F1)

Style: D3.js v7, Okabe-Ito palette, Nature Methods typography
"""

import json, math
import numpy as np
import pandas as pd
from pathlib import Path

BASE = Path("v2")
EVAL = BASE / "evaluation"
PROC = BASE / "processed"
FIGS = BASE / "figures"
FIGS.mkdir(exist_ok=True)

D3_CDN = "https://d3js.org/d3.v7.min.js"

# ── Okabe-Ito palette ─────────────────────────────────────────────────────────
OI = {
    "blue":       "#0072B2",
    "orange":     "#E69F00",
    "green":      "#009E73",
    "sky":        "#56B4E9",
    "vermillion": "#D55E00",
    "pink":       "#CC79A7",
    "yellow":     "#F0E442",
    "black":      "#000000",
    "gray":       "#999999",
    "lgray":      "#cccccc",
}

HTML_HEADER = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{title}</title>
<script src="{d3}"></script>
<style>
  body {{ font-family: Arial, sans-serif; font-size: 12px; margin: 20px; background: #fff; }}
  .fig-title {{ font-size: 13px; font-weight: bold; fill: #111; }}
  .panel-label {{ font-size: 13px; font-weight: bold; fill: #111; }}
  .axis-label {{ font-size: 10px; fill: #333; }}
  .tick text {{ font-size: 8.5px; }}
  .ref-line {{ stroke: #999; stroke-dasharray: 4,3; stroke-width: 1px; }}
  .anno-text {{ font-size: 8.5px; fill: #555; }}
  #tooltip {{ position: absolute; background: rgba(0,0,0,0.8); color: #fff;
    padding: 5px 9px; border-radius: 4px; font-size: 9.5px; pointer-events: none; display: none; }}
</style>
</head>
<body>
<div id="tooltip"></div>
"""

HTML_FOOTER = """</body>\n</html>\n"""


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1: Mouse Temporal Dynamics (T1 / T2 / T3)
# ══════════════════════════════════════════════════════════════════════════════

def make_fig1():
    with open(EVAL / "T_temporal_summary.json") as f:
        T = json.load(f)

    # ── Panel A: T1 – ISS-T vs LAR preservation artifact (RR-8 liver) ────────
    t1b = T["T1"]["T1b_RR-8_liver"]
    t1_bars = [
        {"cond": "FLT", "gene": t1b["conditions"]["FLT_gene"]["auroc"],
         "gene_lo": t1b["conditions"]["FLT_gene"]["ci_low"],
         "gene_hi": t1b["conditions"]["FLT_gene"]["ci_high"],
         "path": t1b["conditions"]["FLT_pathway"]["auroc"],
         "path_lo": t1b["conditions"]["FLT_pathway"]["ci_low"],
         "path_hi": t1b["conditions"]["FLT_pathway"]["ci_high"]},
        {"cond": "GC", "gene": t1b["conditions"]["GC_gene"]["auroc"],
         "gene_lo": t1b["conditions"]["GC_gene"]["ci_low"],
         "gene_hi": t1b["conditions"]["GC_gene"]["ci_high"],
         "path": t1b["conditions"]["GC_pathway"]["auroc"],
         "path_lo": t1b["conditions"]["GC_pathway"]["ci_low"],
         "path_hi": t1b["conditions"]["GC_pathway"]["ci_high"]},
        {"cond": "BSL", "gene": t1b["conditions"]["BSL_gene"]["auroc"],
         "gene_lo": t1b["conditions"]["BSL_gene"]["ci_low"],
         "gene_hi": t1b["conditions"]["BSL_gene"]["ci_high"],
         "path": t1b["conditions"]["BSL_pathway"]["auroc"],
         "path_lo": t1b["conditions"]["BSL_pathway"]["ci_low"],
         "path_hi": t1b["conditions"]["BSL_pathway"]["ci_high"]},
    ]

    # ── Panel B: T2 – post-flight recovery (RR-6 pathway scatter) ────────────
    pr6 = T["T2"]["T2_RR-6"]["pathway_recovery"]
    rr6_ratio = round(T["T2"]["T2_RR-6"]["pca_recovery_ratio"], 3)
    rr8_ratio = round(T["T2"]["T2_RR-8"]["pca_recovery_ratio"], 3)

    t2_pts = []
    for pw, v in pr6.items():
        df = v.get("delta_flight", 0)
        dr = v.get("delta_return", 0)
        rf = v.get("recovery_fraction")
        rev = v.get("direction_reversed", False)
        if abs(df) >= 0.10:
            if rev:
                cat = "overshoot"
            elif rf is not None and rf > 0.1:
                cat = "recover"
            else:
                cat = "persist"
            t2_pts.append({
                "pw": pw.replace("HALLMARK_", "").replace("_", " ").title()[:22],
                "df": round(df, 3), "dr": round(dr, 3),
                "rf": round(rf, 3) if rf is not None else None,
                "cat": cat
            })

    t2_data = {
        "pts": t2_pts,
        "rr6_ratio": rr6_ratio,
        "rr8_ratio": rr8_ratio,
    }

    # ── Panel C: T3 – age amplification ──────────────────────────────────────
    sub = T["T3"]["subtasks"]
    t3_bars = [
        {"age": "OLD", "feat": "gene",
         "auroc": round(sub["T3d_OLD_gene"]["auroc"], 3),
         "ci_lo": round(sub["T3d_OLD_gene"]["ci_low"], 3),
         "ci_hi": round(sub["T3d_OLD_gene"]["ci_high"], 3)},
        {"age": "YNG", "feat": "gene",
         "auroc": round(sub["T3d_YNG_gene"]["auroc"], 3),
         "ci_lo": round(sub["T3d_YNG_gene"]["ci_low"], 3),
         "ci_hi": round(sub["T3d_YNG_gene"]["ci_high"], 3)},
        {"age": "OLD", "feat": "pathway",
         "auroc": round(sub["T3d_OLD_pathway"]["auroc"], 3),
         "ci_lo": round(sub["T3d_OLD_pathway"].get("ci_low", sub["T3d_OLD_pathway"]["auroc"] - 0.07), 3),
         "ci_hi": round(sub["T3d_OLD_pathway"].get("ci_high", min(1, sub["T3d_OLD_pathway"]["auroc"] + 0.07)), 3)},
        {"age": "YNG", "feat": "pathway",
         "auroc": round(sub["T3d_YNG_pathway"]["auroc"], 3),
         "ci_lo": round(sub["T3d_YNG_pathway"].get("ci_low", sub["T3d_YNG_pathway"]["auroc"] - 0.10), 3),
         "ci_hi": round(sub["T3d_YNG_pathway"].get("ci_high", min(1, sub["T3d_YNG_pathway"]["auroc"] + 0.10)), 3)},
    ]

    fig_data = {
        "t1": t1_bars,
        "t2": t2_data,
        "t3": t3_bars,
        "colors": OI,
    }

    html = HTML_HEADER.format(title="Fig1: Mouse Temporal Dynamics", d3=D3_CDN)
    html += f"""
<svg id="fig1"></svg>
<script>
const D = {json.dumps(fig_data)};
const C = D.colors;

// ── Layout ────────────────────────────────────────────────────────────────────
const totalW = 920, totalH = 480;
const mT = 50, mB = 40;
const svg = d3.select("#fig1").attr("width", totalW).attr("height", totalH);
const tip = d3.select("#tooltip");

// Panel offsets: A(0-280), B(300-610), C(630-880)
const pA = {{x: 55, y: mT, w: 230, h: 360}};
const pB = {{x: 320, y: mT, w: 280, h: 360}};
const pC = {{x: 650, y: mT, w: 240, h: 360}};

// Panel labels
[["A", pA.x-40, mT-22], ["B", pB.x-20, mT-22], ["C", pC.x-20, mT-22]].forEach(([l,x,y]) => {{
  svg.append("text").attr("class","panel-label").attr("x",x).attr("y",y).text(l);
}});

// Figure title
svg.append("text").attr("class","fig-title").attr("x",totalW/2).attr("y",18)
  .attr("text-anchor","middle").text("Fig. 1 · Mouse Spaceflight Transcriptomics: Temporal Dynamics");

// ── Panel A: T1 — ISS-T vs LAR AUROC (RR-8 liver) ────────────────────────────
(function() {{
  const gA = svg.append("g").attr("transform", `translate(${{pA.x}},${{pA.y}})`);

  gA.append("text").attr("x", pA.w/2).attr("y", -8)
    .attr("text-anchor","middle").attr("font-size","10px").attr("font-weight","bold")
    .text("T1: Preservation artifact");
  gA.append("text").attr("x", pA.w/2).attr("y", 6)
    .attr("text-anchor","middle").attr("font-size","8.5px").attr("fill","#555")
    .text("ISS-T vs LAR AUROC (RR-8 liver)");

  const conds = D.t1.map(d=>d.cond);
  const yBand = d3.scaleBand().domain(conds).range([0, pA.h]).padding(0.3);
  const xA = d3.scaleLinear().domain([0.4, 1.05]).range([0, pA.w]);
  const bw = yBand.bandwidth() * 0.42;

  // Gridlines
  [0.5, 0.6, 0.7, 0.8, 0.9, 1.0].forEach(v => {{
    gA.append("line").attr("class","ref-line")
      .attr("x1",xA(v)).attr("x2",xA(v)).attr("y1",0).attr("y2",pA.h)
      .attr("stroke-dasharray", v===0.5?"4,3":"2,2").attr("opacity",0.5);
  }});

  // Bars for gene (top) and pathway (bottom) within each condition
  D.t1.forEach(d => {{
    const y0 = yBand(d.cond);
    // Gene bar
    gA.append("rect").attr("x",xA(0.5)).attr("y",y0)
      .attr("width", xA(d.gene)-xA(0.5)).attr("height",bw)
      .attr("fill", C.blue).attr("opacity",0.9)
      .on("mouseover", ev => tip.style("display","block").html(
        `<b>${{d.cond}} Gene</b><br>AUROC = ${{d.gene.toFixed(3)}}<br>95% CI: [${{d.gene_lo.toFixed(3)}}, ${{d.gene_hi.toFixed(3)}}]`))
      .on("mousemove", ev => tip.style("left",(ev.pageX+12)+"px").style("top",(ev.pageY-20)+"px"))
      .on("mouseout", () => tip.style("display","none"));
    // Gene CI
    gA.append("line").attr("x1",xA(d.gene_lo)).attr("x2",xA(d.gene_hi))
      .attr("y1",y0+bw/2).attr("y2",y0+bw/2)
      .attr("stroke","#333").attr("stroke-width",1.5);
    // Pathway bar
    gA.append("rect").attr("x",xA(0.5)).attr("y",y0+bw+2)
      .attr("width", xA(d.path)-xA(0.5)).attr("height",bw)
      .attr("fill", C.orange).attr("opacity",0.9)
      .on("mouseover", ev => tip.style("display","block").html(
        `<b>${{d.cond}} Pathway</b><br>AUROC = ${{d.path.toFixed(3)}}<br>95% CI: [${{d.path_lo.toFixed(3)}}, ${{d.path_hi.toFixed(3)}}]`))
      .on("mousemove", ev => tip.style("left",(ev.pageX+12)+"px").style("top",(ev.pageY-20)+"px"))
      .on("mouseout", () => tip.style("display","none"));
    gA.append("line").attr("x1",xA(d.path_lo)).attr("x2",xA(d.path_hi))
      .attr("y1",y0+bw+2+bw/2).attr("y2",y0+bw+2+bw/2)
      .attr("stroke","#333").attr("stroke-width",1.5);
    // Cond label
    gA.append("text").attr("x",-5).attr("y",y0+bw+1)
      .attr("text-anchor","end").attr("font-size","9px").attr("dominant-baseline","middle")
      .text(d.cond);
  }});

  // AUROC ref line at 0.5
  gA.append("line").attr("x1",xA(0.5)).attr("x2",xA(0.5))
    .attr("y1",0).attr("y2",pA.h).attr("stroke","#333").attr("stroke-width",1.5);

  // Axis
  gA.append("g").attr("transform",`translate(0,${{pA.h}})`).call(d3.axisBottom(xA).tickValues([0.5,0.6,0.7,0.8,0.9,1.0]).tickFormat(d3.format(".1f"))).selectAll("text").attr("font-size","8px");
  gA.append("text").attr("class","axis-label").attr("x",pA.w/2).attr("y",pA.h+28)
    .attr("text-anchor","middle").text("AUROC (ISS-T vs LAR classifier)");

  // Key annotation: FLT ≈ GC
  gA.append("text").attr("x",pA.w+2).attr("y",yBand("FLT")+bw+1)
    .attr("font-size","8px").attr("fill","#555")
    .text(`${{D.t1[0].gene.toFixed(2)}}`);
  gA.append("text").attr("x",pA.w+2).attr("y",yBand("GC")+bw+1)
    .attr("font-size","8px").attr("fill","#555")
    .text(`${{D.t1[1].gene.toFixed(2)}}`);

  // FLT≈GC brace annotation
  const fltY = yBand("FLT") + bw;
  const gcY  = yBand("GC")  + bw;
  gA.append("text").attr("x",pA.w-4).attr("y",(fltY+gcY)/2-2)
    .attr("text-anchor","end").attr("font-size","8px").attr("fill",C.vermillion)
    .text("Δ≈0 → artifact");

  // Legend
  const legY = pA.h + 42;
  [[C.blue,"Gene"],[C.orange,"Pathway"]].forEach(([col,lbl], i) => {{
    gA.append("rect").attr("x",i*65).attr("y",legY).attr("width",10).attr("height",10).attr("fill",col);
    gA.append("text").attr("x",i*65+13).attr("y",legY+8).attr("font-size","8.5px").text(lbl);
  }});
}})();

// ── Panel B: T2 — Recovery scatter (RR-6 delta_flight vs delta_return) ────────
(function() {{
  const gB = svg.append("g").attr("transform", `translate(${{pB.x}},${{pB.y}})`);

  gB.append("text").attr("x", pB.w/2).attr("y", -8)
    .attr("text-anchor","middle").attr("font-size","10px").attr("font-weight","bold")
    .text("T2: Post-flight recovery");
  gB.append("text").attr("x", pB.w/2).attr("y", 6)
    .attr("text-anchor","middle").attr("font-size","8.5px").attr("fill","#555")
    .text("RR-6 liver: Δ_flight vs Δ_return (NES)");

  const pts = D.t2.pts;
  const allD = pts.flatMap(p => [p.df, p.dr]);
  const ext = 0.12;
  const dMin = d3.min(allD) - ext, dMax = d3.max(allD) + ext;
  const scX = d3.scaleLinear().domain([dMin, dMax]).range([0, pB.w]);
  const scY = d3.scaleLinear().domain([dMin, dMax]).range([pB.h, 0]);

  // Zones
  // y=0 line (perfect recovery)
  gB.append("line").attr("x1",0).attr("x2",pB.w)
    .attr("y1",scY(0)).attr("y2",scY(0)).attr("class","ref-line");
  // y=x line (no recovery)
  gB.append("line").attr("x1",scX(dMin)).attr("x2",scX(dMax))
    .attr("y1",scY(dMin)).attr("y2",scY(dMax))
    .attr("stroke","#bbb").attr("stroke-dasharray","3,3").attr("stroke-width",1);
  // x=0 line
  gB.append("line").attr("x1",scX(0)).attr("x2",scX(0))
    .attr("y1",0).attr("y2",pB.h).attr("class","ref-line");

  // Zone labels
  gB.append("text").attr("x",scX(0.07)).attr("y",scY(0.15))
    .attr("font-size","8px").attr("fill",C.orange).text("← overshoot");
  gB.append("text").attr("x",scX(-0.38)).attr("y",scY(-0.18))
    .attr("font-size","8px").attr("fill",C.vermillion).text("persist ↓");
  gB.append("text").attr("x",scX(-0.3)).attr("y",scY(-0.07))
    .attr("font-size","8px").attr("fill",C.green).text("recover →");

  // Color by category
  const catColor = {{
    "overshoot": C.orange,
    "recover":   C.green,
    "persist":   C.vermillion,
  }};

  // Points
  gB.selectAll(".pt").data(pts).join("circle")
    .attr("cx", d => scX(d.df)).attr("cy", d => scY(d.dr))
    .attr("r", 5).attr("fill", d => catColor[d.cat]).attr("opacity",0.85)
    .attr("stroke","#fff").attr("stroke-width",0.5)
    .on("mouseover", (ev,d) => tip.style("display","block").html(
      `<b>${{d.pw}}</b><br>Δ_flight: ${{d.df}}<br>Δ_return: ${{d.dr}}<br>Recovery: ${{d.rf !== null ? (d.rf*100).toFixed(0)+"%" : "—"}}<br>Category: ${{d.cat}}`))
    .on("mousemove", ev => tip.style("left",(ev.pageX+12)+"px").style("top",(ev.pageY-20)+"px"))
    .on("mouseout", () => tip.style("display","none"));

  // Axes
  gB.append("g").attr("transform",`translate(0,${{pB.h}})`).call(d3.axisBottom(scX).ticks(5)).selectAll("text").attr("font-size","8px");
  gB.append("g").call(d3.axisLeft(scY).ticks(5)).selectAll("text").attr("font-size","8px");
  gB.append("text").attr("class","axis-label").attr("x",pB.w/2).attr("y",pB.h+28)
    .attr("text-anchor","middle").text("Δ_flight (flight-induced NES change)");
  gB.append("text").attr("class","axis-label").attr("transform","rotate(-90)")
    .attr("x",-pB.h/2).attr("y",-36).attr("text-anchor","middle")
    .text("Δ_return (post-return NES change)");

  // Recovery ratios annotation
  const ratioY = pB.h + 42;
  gB.append("text").attr("x",0).attr("y",ratioY+8).attr("font-size","8.5px")
    .text(`PCA recovery: RR-6 ${{(D.t2.rr6_ratio*100).toFixed(0)}}%  |  RR-8 ${{(D.t2.rr8_ratio*100).toFixed(0)}}%`);

  // Legend
  const legY = pB.h + 52;
  Object.entries(catColor).forEach(([cat, col], i) => {{
    gB.append("circle").attr("cx",i*80+6).attr("cy",legY+5).attr("r",5).attr("fill",col);
    gB.append("text").attr("x",i*80+13).attr("y",legY+9).attr("font-size","8.5px")
      .text(cat.charAt(0).toUpperCase()+cat.slice(1));
  }});
}})();

// ── Panel C: T3 — Age amplification ──────────────────────────────────────────
(function() {{
  const gC = svg.append("g").attr("transform", `translate(${{pC.x}},${{pC.y}})`);

  gC.append("text").attr("x", pC.w/2).attr("y", -8)
    .attr("text-anchor","middle").attr("font-size","10px").attr("font-weight","bold")
    .text("T3: Age amplification");
  gC.append("text").attr("x", pC.w/2).attr("y", 6)
    .attr("text-anchor","middle").attr("font-size","8.5px").attr("fill","#555")
    .text("Spaceflight AUROC: OLD vs YNG (RR-8)");

  const groups = ["gene", "pathway"];
  const ages   = ["OLD", "YNG"];
  const yBand  = d3.scaleBand().domain(groups).range([0, pC.h]).padding(0.28);
  const xC     = d3.scaleLinear().domain([0.4, 1.05]).range([0, pC.w]);
  const bw     = yBand.bandwidth() * 0.42;
  const ageColor = {{OLD: C.blue, YNG: C.sky}};

  // Gridlines
  [0.5, 0.6, 0.7, 0.8, 0.9, 1.0].forEach(v => {{
    gC.append("line").attr("class","ref-line").attr("opacity",0.5)
      .attr("x1",xC(v)).attr("x2",xC(v)).attr("y1",0).attr("y2",pC.h);
  }});

  D.t3.forEach(d => {{
    const feat = d.feat;
    const y0 = yBand(feat) + (d.age === "OLD" ? 0 : bw + 2);
    gC.append("rect").attr("x",xC(0.5)).attr("y",y0)
      .attr("width", xC(d.auroc)-xC(0.5)).attr("height", bw)
      .attr("fill", ageColor[d.age]).attr("opacity",0.9)
      .on("mouseover", ev => tip.style("display","block").html(
        `<b>${{d.age}} ${{feat}}</b><br>AUROC = ${{d.auroc.toFixed(3)}}<br>95% CI: [${{d.ci_lo.toFixed(3)}}, ${{d.ci_hi.toFixed(3)}}]`))
      .on("mousemove", ev => tip.style("left",(ev.pageX+12)+"px").style("top",(ev.pageY-20)+"px"))
      .on("mouseout", () => tip.style("display","none"));
    // CI whisker
    gC.append("line").attr("x1",xC(d.ci_lo)).attr("x2",xC(d.ci_hi))
      .attr("y1",y0+bw/2).attr("y2",y0+bw/2).attr("stroke","#333").attr("stroke-width",1.5);
  }});

  // Feature labels
  groups.forEach(g => {{
    gC.append("text").attr("x",-5).attr("y",yBand(g)+bw+1)
      .attr("text-anchor","end").attr("dominant-baseline","middle")
      .attr("font-size","9px").text(g.charAt(0).toUpperCase()+g.slice(1));
  }});

  gC.append("line").attr("x1",xC(0.5)).attr("x2",xC(0.5))
    .attr("y1",0).attr("y2",pC.h).attr("stroke","#333").attr("stroke-width",1.5);
  gC.append("g").attr("transform",`translate(0,${{pC.h}})`).call(d3.axisBottom(xC).tickValues([0.5,0.6,0.7,0.8,0.9,1.0]).tickFormat(d3.format(".1f"))).selectAll("text").attr("font-size","8px");
  gC.append("text").attr("class","axis-label").attr("x",pC.w/2).attr("y",pC.h+28)
    .attr("text-anchor","middle").text("Spaceflight AUROC (FLT vs GC)");

  // Delta annotation
  const oldGene = D.t3.find(d=>d.age==="OLD"&&d.feat==="gene").auroc;
  const yngGene = D.t3.find(d=>d.age==="YNG"&&d.feat==="gene").auroc;
  gC.append("text").attr("x",pC.w-4).attr("y",yBand("gene")+bw+1)
    .attr("text-anchor","end").attr("font-size","8px").attr("fill",C.vermillion)
    .text(`Δ = +${{(oldGene-yngGene).toFixed(3)}}`);

  // Legend
  const legY = pC.h + 42;
  [["OLD", C.blue], ["YNG", C.sky]].forEach(([lbl,col], i) => {{
    gC.append("rect").attr("x",i*60).attr("y",legY).attr("width",10).attr("height",10).attr("fill",col);
    gC.append("text").attr("x",i*60+13).attr("y",legY+8).attr("font-size","8.5px").text(lbl);
  }});
}})();
</script>
"""
    html += HTML_FOOTER

    out = FIGS / "Fig1_temporal.html"
    out.write_text(html, encoding="utf-8")
    print(f"Saved: {out}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2: Cross-Species Conservation (E1 / E2 / E3)
# ══════════════════════════════════════════════════════════════════════════════

def make_fig2():
    with open(EVAL / "E1_crossspecies_nes.json") as f:
        E1 = json.load(f)
    with open(EVAL / "E2_mission_conservation.json") as f:
        E2 = json.load(f)
    with open(EVAL / "E3_cfrna_origin.json") as f:
        E3 = json.load(f)

    def short(p):
        return p.replace("HALLMARK_", "").replace("_", " ").title()

    # E1 scatter points
    e1_pts = [{"pw": short(d["pathway"]),
               "x": round(d["mouse_nes_mean"], 3),
               "y": round(d["human_nes"], 3),
               "conc": d["concordant"]}
              for d in E1["pathways"]]

    # E1 r value
    e1_res = E1["results"]["mission_averaged"]

    # E2 bar data
    e2_bars = [
        {"label": "I4 (3d)",    "dur": 3,   "r": round(E2["results"]["I4_vs_mouse"]["spearman_r"],3),
         "ci_lo": round(E2["results"]["I4_vs_mouse"]["ci_low"],3),
         "ci_hi": round(E2["results"]["I4_vs_mouse"]["ci_high"],3),
         "sig": E2["results"]["I4_vs_mouse"]["p_permutation"] < 0.05},
        {"label": "JAXA (120d)","dur": 120,  "r": round(E2["results"]["JAXA_vs_mouse"]["spearman_r"],3),
         "ci_lo": round(E2["results"]["JAXA_vs_mouse"]["ci_low"],3),
         "ci_hi": round(E2["results"]["JAXA_vs_mouse"]["ci_high"],3),
         "sig": E2["results"]["JAXA_vs_mouse"]["p_permutation"] < 0.05},
    ]

    # E3 bar data (cell types)
    CT_ORDER = ["CD14+ Monocyte","CD16+ Monocyte","Dendritic Cell","Natural Killer Cell",
                "B Cell","CD4+ T Cell","CD8+ T Cell","Other T Cell","Other","PBMC Pseudobulk"]
    e3_bars = []
    for ct in CT_ORDER:
        if ct in E3["results"]:
            r = E3["results"][ct]
            e3_bars.append({"ct": ct, "r": r["r"], "ci_lo": r["ci_lo"],
                             "ci_hi": r["ci_hi"], "p": r["p"], "n": r["n"],
                             "sig": r["p"] < 0.05})

    fig_data = {
        "e1_pts": e1_pts,
        "e1_r": round(e1_res["spearman_r"], 3),
        "e1_p": e1_res["p_permutation"],
        "e1_n": e1_res["n_pathways"],
        "e2": e2_bars,
        "e3": e3_bars,
        "colors": OI,
    }

    html = HTML_HEADER.format(title="Fig2: Cross-Species Conservation", d3=D3_CDN)
    html += f"""
<svg id="fig2"></svg>
<script>
const D = {json.dumps(fig_data)};
const C = D.colors;

const totalW = 940, totalH = 440;
const mT = 50, mB = 40;
const svg = d3.select("#fig2").attr("width", totalW).attr("height", totalH);
const tip = d3.select("#tooltip");

const pA = {{x: 55, y: mT, w: 270, h: 330}};
const pB = {{x: 380, y: mT, w: 165, h: 330}};
const pC = {{x: 590, y: mT, w: 320, h: 330}};

// Panel labels
[["A", pA.x-40, mT-22], ["B", pB.x-22, mT-22], ["C", pC.x-22, mT-22]].forEach(([l,x,y]) => {{
  svg.append("text").attr("class","panel-label").attr("x",x).attr("y",y).text(l);
}});
svg.append("text").attr("class","fig-title").attr("x",totalW/2).attr("y",18)
  .attr("text-anchor","middle").text("Fig. 2 · Cross-Species Pathway Conservation in Spaceflight");

// ── Panel A: E1 – JAXA cfRNA vs. mouse liver NES scatter ─────────────────────
(function() {{
  const gA = svg.append("g").attr("transform",`translate(${{pA.x}},${{pA.y}})`);
  gA.append("text").attr("x",pA.w/2).attr("y",-8)
    .attr("text-anchor","middle").attr("font-size","10px").attr("font-weight","bold")
    .text("E1: Human cfRNA vs. mouse liver NES");
  gA.append("text").attr("x",pA.w/2).attr("y",6)
    .attr("text-anchor","middle").attr("font-size","8.5px").attr("fill","#555")
    .text(`JAXA CFE (120d) · ${{D.e1_n}} Hallmark pathways · r = ${{D.e1_r.toFixed(3)}}*`);

  const allX = D.e1_pts.map(d=>d.x), allY = D.e1_pts.map(d=>d.y);
  const xPad = 0.2, yPad = 0.2;
  const scX = d3.scaleLinear().domain([d3.min(allX)-xPad, d3.max(allX)+xPad]).range([0, pA.w]);
  const scY = d3.scaleLinear().domain([d3.min(allY)-yPad, d3.max(allY)+yPad]).range([pA.h, 0]);

  // Ref lines
  [[scX(0), scX(0), 0, pA.h], [0, pA.w, scY(0), scY(0)]].forEach(([x1,x2,y1,y2]) => {{
    gA.append("line").attr("class","ref-line").attr("x1",x1).attr("x2",x2).attr("y1",y1).attr("y2",y2);
  }});

  // OLS fit line
  const meanX = d3.mean(allX), meanY = d3.mean(allY);
  const slope = d3.sum(allX.map((x,i)=>(x-meanX)*(allY[i]-meanY))) /
                d3.sum(allX.map(x=>(x-meanX)**2));
  const b = meanY - slope*meanX;
  const x0 = d3.min(allX)-xPad, x1 = d3.max(allX)+xPad;
  gA.append("line").attr("x1",scX(x0)).attr("x2",scX(x1))
    .attr("y1",scY(slope*x0+b)).attr("y2",scY(slope*x1+b))
    .attr("stroke",C.blue).attr("stroke-width",1.5).attr("opacity",0.6);

  // Dots
  gA.selectAll(".dot").data(D.e1_pts).join("circle")
    .attr("cx", d=>scX(d.x)).attr("cy", d=>scY(d.y)).attr("r",5)
    .attr("fill", d => d.conc ? C.blue : C.vermillion)
    .attr("opacity",0.8).attr("stroke","#fff").attr("stroke-width",0.5)
    .on("mouseover", (ev,d) => tip.style("display","block").html(
      `<b>${{d.pw}}</b><br>Mouse NES: ${{d.x.toFixed(3)}}<br>Human NES: ${{d.y.toFixed(3)}}<br>${{d.conc ? "Concordant" : "Discordant"}}`))
    .on("mousemove", ev => tip.style("left",(ev.pageX+12)+"px").style("top",(ev.pageY-20)+"px"))
    .on("mouseout", () => tip.style("display","none"));

  // Axes
  gA.append("g").attr("transform",`translate(0,${{pA.h}})`).call(d3.axisBottom(scX).ticks(5)).selectAll("text").attr("font-size","8px");
  gA.append("g").call(d3.axisLeft(scY).ticks(5)).selectAll("text").attr("font-size","8px");
  gA.append("text").attr("class","axis-label").attr("x",pA.w/2).attr("y",pA.h+28)
    .attr("text-anchor","middle").text("Mouse liver NES (mission-averaged)");
  gA.append("text").attr("class","axis-label").attr("transform","rotate(-90)")
    .attr("x",-pA.h/2).attr("y",-38).attr("text-anchor","middle").text("Human cfRNA NES (JAXA CFE)");

  // r annotation
  gA.append("text").attr("x",4).attr("y",12).attr("font-size","9px").attr("fill",C.blue)
    .text(`r = ${{D.e1_r.toFixed(3)}}* (p = ${{D.e1_p.toFixed(3)}})`);

  // Legend
  [[C.blue,"Concordant"],[C.vermillion,"Discordant"]].forEach(([col,lbl],i) => {{
    gA.append("circle").attr("cx",i*90+6).attr("cy",pA.h+46).attr("r",5).attr("fill",col);
    gA.append("text").attr("x",i*90+13).attr("y",pA.h+50).attr("font-size","8.5px").text(lbl);
  }});
}})();

// ── Panel B: E2 – Duration effect ────────────────────────────────────────────
(function() {{
  const gB = svg.append("g").attr("transform",`translate(${{pB.x}},${{pB.y}})`);
  gB.append("text").attr("x",pB.w/2).attr("y",-8)
    .attr("text-anchor","middle").attr("font-size","10px").attr("font-weight","bold")
    .text("E2: Duration effect");
  gB.append("text").attr("x",pB.w/2).attr("y",6)
    .attr("text-anchor","middle").attr("font-size","8.5px").attr("fill","#555")
    .text("cfRNA vs. mouse Spearman r");

  const missions = D.e2.map(d=>d.label);
  const yB = d3.scaleBand().domain(missions).range([0, pB.h*0.5]).padding(0.3);
  const xB = d3.scaleLinear().domain([-0.5, 0.7]).range([0, pB.w]);
  const xRef = xB(0);

  // Zero line
  gB.append("line").attr("x1",xRef).attr("x2",xRef).attr("y1",0).attr("y2",pB.h*0.5)
    .attr("stroke","#333").attr("stroke-width",1.2);

  D.e2.forEach(d => {{
    const y0 = yB(d.label);
    const bw = yB.bandwidth();
    const col = d.sig ? C.blue : C.lgray;
    gB.append("rect")
      .attr("x", d.r >= 0 ? xRef : xB(d.r)).attr("y",y0)
      .attr("width", Math.abs(xB(d.r)-xRef)).attr("height",bw)
      .attr("fill",col).attr("opacity",0.9)
      .on("mouseover", ev => tip.style("display","block").html(
        `<b>${{d.label}}</b><br>r = ${{d.r.toFixed(3)}}<br>95% CI: [${{d.ci_lo.toFixed(3)}}, ${{d.ci_hi.toFixed(3)}}]<br>${{d.sig ? "p < 0.05 *" : "NS"}}`))
      .on("mousemove", ev => tip.style("left",(ev.pageX+12)+"px").style("top",(ev.pageY-20)+"px"))
      .on("mouseout", () => tip.style("display","none"));
    // CI
    gB.append("line").attr("x1",xB(d.ci_lo)).attr("x2",xB(d.ci_hi))
      .attr("y1",y0+bw/2).attr("y2",y0+bw/2).attr("stroke","#333").attr("stroke-width",1.5);
    // Label
    gB.append("text").attr("x",-3).attr("y",y0+bw/2+3)
      .attr("text-anchor","end").attr("font-size","8.5px").text(d.label);
    // r value
    const labelX = d.r >= 0 ? xB(d.ci_hi)+3 : xB(d.ci_lo)-3;
    gB.append("text").attr("x",labelX).attr("y",y0+bw/2+3)
      .attr("text-anchor", d.r >= 0 ? "start" : "end").attr("font-size","8.5px")
      .attr("fill", d.sig ? "#000" : "#888")
      .text(`${{d.r >= 0 ? "+" : ""}}${{d.r.toFixed(3)}}${{d.sig ? "*" : ""}}`);
  }});

  // Duration arrow annotation
  const arrowY = pB.h * 0.5 + 20;
  gB.append("text").attr("x",pB.w/2).attr("y",arrowY)
    .attr("text-anchor","middle").attr("font-size","8.5px").attr("fill",C.blue)
    .text("Longer mission → stronger conservation");

  gB.append("g").attr("transform",`translate(0,${{pB.h*0.5}})`).call(d3.axisBottom(xB).tickValues([-0.4,-0.2,0,0.2,0.4,0.6])).selectAll("text").attr("font-size","8px");
  gB.append("text").attr("class","axis-label").attr("x",pB.w/2).attr("y",pB.h*0.5+30)
    .attr("text-anchor","middle").text("Spearman r (vs. mouse NES)");

  // NS legend
  gB.append("rect").attr("x",0).attr("y",pB.h*0.5+46).attr("width",10).attr("height",10).attr("fill",C.lgray);
  gB.append("text").attr("x",13).attr("y",pB.h*0.5+54).attr("font-size","8.5px").text("NS");
  gB.append("rect").attr("x",38).attr("y",pB.h*0.5+46).attr("width",10).attr("height",10).attr("fill",C.blue);
  gB.append("text").attr("x",51).attr("y",pB.h*0.5+54).attr("font-size","8.5px").text("p<0.05");
}})();

// ── Panel C: E3 – cfRNA cellular origin ──────────────────────────────────────
(function() {{
  const gC = svg.append("g").attr("transform",`translate(${{pC.x}},${{pC.y}})`);
  gC.append("text").attr("x",pC.w/2).attr("y",-8)
    .attr("text-anchor","middle").attr("font-size","10px").attr("font-weight","bold")
    .text("E3: cfRNA cellular origin");
  gC.append("text").attr("x",pC.w/2).attr("y",6)
    .attr("text-anchor","middle").attr("font-size","8.5px").attr("fill","#555")
    .text("I4 cfRNA vs. PBMC cell-type NES (R+1 vs pre)");

  const cts = D.e3.map(d=>d.ct);
  const yC = d3.scaleBand().domain(cts).range([0, pC.h]).padding(0.25);
  const xC = d3.scaleLinear().domain([-0.75, 0.2]).range([0, pC.w]);
  const xRef = xC(0);

  // Zero ref
  gC.append("line").attr("x1",xRef).attr("x2",xRef).attr("y1",0).attr("y2",pC.h)
    .attr("stroke","#333").attr("stroke-width",1.2);

  D.e3.forEach(d => {{
    const y0 = yC(d.ct);
    const bw = yC.bandwidth();
    const col = d.sig ? C.pink : "#BBBBDD";
    gC.append("rect")
      .attr("x", d.r >= 0 ? xRef : xC(d.r)).attr("y",y0)
      .attr("width", Math.abs(xC(d.r)-xRef)).attr("height",bw)
      .attr("fill",col).attr("opacity",0.9)
      .on("mouseover", ev => tip.style("display","block").html(
        `<b>${{d.ct}}</b><br>r = ${{d.r.toFixed(3)}}<br>95% CI: [${{d.ci_lo.toFixed(3)}}, ${{d.ci_hi.toFixed(3)}}]<br>p = ${{d.p.toFixed(4)}}${{d.sig ? " *" : ""}}<br>n = ${{d.n}} pathways`))
      .on("mousemove", ev => tip.style("left",(ev.pageX+12)+"px").style("top",(ev.pageY-20)+"px"))
      .on("mouseout", () => tip.style("display","none"));
    // CI
    gC.append("line").attr("x1",xC(d.ci_lo)).attr("x2",xC(d.ci_hi))
      .attr("y1",y0+bw/2).attr("y2",y0+bw/2).attr("stroke","#333").attr("stroke-width",1.5);
    // Label
    gC.append("text").attr("x",-3).attr("y",y0+bw/2+3)
      .attr("text-anchor","end").attr("font-size","8.5px").text(d.ct);
    // r value
    gC.append("text").attr("x",xC(d.ci_lo)-3).attr("y",y0+bw/2+3)
      .attr("text-anchor","end").attr("font-size","8px")
      .attr("fill", d.sig ? "#000" : "#888")
      .text(`${{d.r.toFixed(3)}}${{d.sig ? "*" : ""}}`);
  }});

  gC.append("g").attr("transform",`translate(0,${{pC.h}})`).call(d3.axisBottom(xC).ticks(6)).selectAll("text").attr("font-size","8px");
  gC.append("text").attr("class","axis-label").attr("x",pC.w/2).attr("y",pC.h+28)
    .attr("text-anchor","middle").text("Spearman r (cell-type NES vs. cfRNA NES)");

  // Anti-corr annotation
  gC.append("text").attr("x",xRef-4).attr("y",12)
    .attr("text-anchor","end").attr("font-size","8px").attr("fill",C.vermillion)
    .text("← All negative (anti-correlated)");

  // Legend
  [[C.pink,"p<0.05"],["#BBBBDD","NS"]].forEach(([col,lbl],i) => {{
    gC.append("rect").attr("x",i*65).attr("y",pC.h+44).attr("width",10).attr("height",10).attr("fill",col);
    gC.append("text").attr("x",i*65+13).attr("y",pC.h+52).attr("font-size","8.5px").text(lbl);
  }});
}})();
</script>
"""
    html += HTML_FOOTER

    out = FIGS / "Fig2_crossspecies.html"
    out.write_text(html, encoding="utf-8")
    print(f"Saved: {out}")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3: PBMC Cell-Type Resolution (F1)
# ══════════════════════════════════════════════════════════════════════════════

def make_fig3():
    ct_df = pd.read_csv(PROC / "F1_scrna/i4_snrnaseq_celltype_fgsea.csv")
    tp_df = pd.read_csv(PROC / "F1_scrna/i4_snrnaseq_temporal_fgsea.csv")

    CT_ORDER = ["CD14+ Monocyte","CD16+ Monocyte","Dendritic Cell","Natural Killer Cell",
                "B Cell","CD4+ T Cell","CD8+ T Cell","Other T Cell","Other","PBMC Pseudobulk"]
    SHEET_ORDER = ["I4-FP1","I4-LP1","I4-LP3","I4-RP1","I4-RP3"]
    SHEET_LABEL = {"I4-FP1":"FP1 (R+1 vs pre)","I4-LP1":"LP1 (R+45 vs pre)",
                   "I4-LP3":"LP3 (R+45/82 vs pre)","I4-RP1":"RP1 (R+45 vs R+1)",
                   "I4-RP3":"RP3 (R+45/82 vs R+1)"}

    # NES heatmap: all sig pathways + top by mean |NES| to fill up to 14 rows
    sig_counts = (ct_df[ct_df["padj"] < 0.05]
                  .groupby("pathway")["Cell_Type"].count()
                  .sort_values(ascending=False))
    sig_pws = list(sig_counts.index)  # all pathways with ≥1 sig cell type

    # Add top pathways by mean |NES| to reach 14 rows
    nes_all = ct_df.pivot_table(index="pathway", columns="Cell_Type", values="NES").fillna(0)
    mean_abs_nes = nes_all.abs().mean(axis=1).sort_values(ascending=False)
    extra_pws = [p for p in mean_abs_nes.index if p not in sig_pws]
    top_pws = list(sig_pws) + extra_pws[: max(0, 14 - len(sig_pws))]
    # Sort: sig first (by sig_count), then by mean |NES|
    top_pws = top_pws[:14]

    nes_piv = (ct_df[ct_df["pathway"].isin(top_pws)]
               .pivot_table(index="pathway", columns="Cell_Type", values="NES")
               .reindex(columns=CT_ORDER, fill_value=0))
    nes_piv = nes_piv.reindex(top_pws)

    def short(p):
        return p.replace("HALLMARK_","").replace("_"," ").title()

    cells_nes = []
    for ri, pw in enumerate(top_pws):
        for ci, ct in enumerate(CT_ORDER):
            if ct in nes_piv.columns:
                v = nes_piv.loc[pw, ct]
                # sig info
                sub = ct_df[(ct_df["pathway"]==pw) & (ct_df["Cell_Type"]==ct)]
                is_sig = bool(len(sub) > 0 and sub["padj"].values[0] < 0.05)
                cells_nes.append({"row": ri, "col": ci,
                                  "nes": round(float(v), 3) if not np.isnan(v) else 0,
                                  "sig": is_sig})
            else:
                cells_nes.append({"row": ri, "col": ci, "nes": 0, "sig": False})

    # Temporal sig count heatmap
    sig_pivot = (tp_df[tp_df["padj"] < 0.05]
                 .groupby(["sheet","Cell_Type"])["pathway"]
                 .count()
                 .unstack(fill_value=0)
                 .reindex(index=SHEET_ORDER, columns=CT_ORDER, fill_value=0))

    cells_sig = []
    for ri, sheet in enumerate(SHEET_ORDER):
        for ci, ct in enumerate(CT_ORDER):
            cells_sig.append({"row": ri, "col": ci,
                               "n": int(sig_pivot.loc[sheet, ct]) if sheet in sig_pivot.index else 0})

    nes_vals = [c["nes"] for c in cells_nes if c["nes"] != 0]
    abs_max = max(abs(min(nes_vals)), abs(max(nes_vals))) if nes_vals else 3.0

    fig_data = {
        "pw_labels": [short(p) for p in top_pws],
        "ct_labels": CT_ORDER,
        "cells_nes": cells_nes,
        "abs_max_nes": round(abs_max, 2),
        "sheet_labels": [SHEET_LABEL[s] for s in SHEET_ORDER],
        "cells_sig": cells_sig,
        "max_sig": int(sig_pivot.values.max()),
        "colors": OI,
    }

    html = HTML_HEADER.format(title="Fig3: PBMC Cell-Type Resolution", d3=D3_CDN)
    html += f"""
<svg id="fig3"></svg>
<script>
const D = {json.dumps(fig_data)};
const C = D.colors;

// Use 40px cell width to fit both panels side by side
const cWA = 40, cWB = 40;
const nCtA = D.ct_labels.length, nPw = D.pw_labels.length;
const nCtB = D.ct_labels.length, nSh = D.sheet_labels.length;
const mT = 55, mLa = 200, mLb = 22, gap = 70, mR = 20;
const pAh = nPw * 28, pBh = nSh * 46;
const pAw = nCtA * cWA, pBw = nCtB * cWB;
const totalW = mLa + pAw + gap + mLb + pBw + mR;
const totalH = mT + Math.max(pAh, pBh) + 90;

const svg = d3.select("#fig3").attr("width", totalW).attr("height", totalH);
const tip = d3.select("#tooltip");

const pA = {{x: mLa, y: mT, w: pAw, h: pAh}};
const pB_x = mLa + pAw + gap + mLb;
const pB = {{x: pB_x, y: mT, w: pBw, h: pBh}};

[["A", pA.x-160, mT-22], ["B", pB.x-mLb, mT-22]].forEach(([l,x,y]) => {{
  svg.append("text").attr("class","panel-label").attr("x",x).attr("y",y).text(l);
}});
svg.attr("width", totalW).attr("height", totalH);
svg.append("text").attr("class","fig-title").attr("x",totalW/2).attr("y",18)
  .attr("text-anchor","middle").text("Fig. 3 · I4 PBMC Cell-Type-Resolved Pathway Response");

// ── Panel A: NES heatmap (FP1: R+1 vs pre-flight) ─────────────────────────────
(function() {{
  const gA = svg.append("g").attr("transform",`translate(${{pA.x}},${{pA.y}})`);
  const cW = cWA, cH = 28;
  const nR = D.pw_labels.length, nC = D.ct_labels.length;

  gA.append("text").attr("x",nC*cW/2).attr("y",-22)
    .attr("text-anchor","middle").attr("font-size","10px").attr("font-weight","bold")
    .text("A: Pathway NES at R+1 (FP1, R+1 vs. pre-flight)");
  gA.append("text").attr("x",nC*cW/2).attr("y",-8)
    .attr("text-anchor","middle").attr("font-size","8.5px").attr("fill","#555")
    .text("Hallmark fGSEA NES · asterisk = padj < 0.05");

  // Diverging color scale: blue (negative) → white → red (positive)
  const colScale = d3.scaleDiverging(d3.interpolateRdBu)
    .domain([D.abs_max_nes, 0, -D.abs_max_nes]);

  // Column headers (rotated)
  gA.selectAll(".colH").data(D.ct_labels).join("text")
    .attr("x",(d,i)=>i*cW+cW/2).attr("y",-4)
    .attr("text-anchor","end").attr("font-size","8px")
    .attr("transform",(d,i)=>`rotate(-38,${{i*cW+cW/2}},-4)`)
    .text(d=>d);

  // Row labels
  gA.selectAll(".rowH").data(D.pw_labels).join("text")
    .attr("x",-6).attr("y",(d,i)=>i*cH+cH/2+3.5)
    .attr("text-anchor","end").attr("font-size","8.5px").text(d=>d);

  // Cells
  const cellG = gA.selectAll(".cG").data(D.cells_nes).join("g")
    .attr("transform",d=>`translate(${{d.col*cW}},${{d.row*cH}})`);
  cellG.append("rect").attr("width",cW).attr("height",cH)
    .attr("fill",d=>d.nes===0?"#f5f5f5":colScale(d.nes))
    .attr("stroke","#fff").attr("stroke-width",0.8)
    .on("mouseover",(ev,d)=>tip.style("display","block").html(
      `<b>${{D.pw_labels[d.row]}}</b><br><b>${{D.ct_labels[d.col]}}</b><br>NES = ${{d.nes.toFixed(3)}}${{d.sig?"  **":""}}`))
    .on("mousemove",ev=>tip.style("left",(ev.pageX+12)+"px").style("top",(ev.pageY-20)+"px"))
    .on("mouseout",()=>tip.style("display","none"));
  // NES text
  cellG.filter(d=>d.nes!==0).append("text").attr("x",cW/2).attr("y",cH/2+3.5)
    .attr("text-anchor","middle").attr("font-size","7.5px")
    .attr("fill",d=>Math.abs(d.nes)>2.2?"#fff":"#222")
    .text(d=>d.nes.toFixed(1));
  // Sig asterisk
  cellG.filter(d=>d.sig).append("text").attr("x",cW-3).attr("y",9)
    .attr("text-anchor","end").attr("font-size","9px").attr("fill","#fff")
    .attr("font-weight","bold").text("*");

  // Legend
  const legG = gA.append("g").attr("transform",`translate(0,${{nR*cH+18}})`);
  const legW = 110, legH = 10;
  const defs = svg.append("defs");
  const grad = defs.append("linearGradient").attr("id","nesGrad").attr("x1","0%").attr("x2","100%");
  d3.range(0,1.01,0.1).forEach(t=>{{
    grad.append("stop").attr("offset",`${{t*100}}%`)
      .attr("stop-color", colScale(D.abs_max_nes - t*2*D.abs_max_nes));
  }});
  legG.append("rect").attr("width",legW).attr("height",legH).style("fill","url(#nesGrad)");
  legG.append("text").attr("y",-3).attr("font-size","8.5px").text("NES");
  [[-D.abs_max_nes, 0], [0, legW/2], [D.abs_max_nes, legW]].forEach(([v,x])=>{{
    legG.append("line").attr("x1",x).attr("x2",x).attr("y1",legH).attr("y2",legH+4).attr("stroke","#666");
    legG.append("text").attr("x",x).attr("y",legH+12).attr("text-anchor","middle").attr("font-size","7.5px").text(v.toFixed(1));
  }});
}})();

// ── Panel B: Temporal sig count heatmap ───────────────────────────────────────
(function() {{
  const gB = svg.append("g").attr("transform",`translate(${{pB.x}},${{pB.y}})`);
  const cW = cWB, cH = 46;
  const nR = D.sheet_labels.length, nC = D.ct_labels.length;
  const maxN = D.max_sig;

  gB.append("text").attr("x",nC*cW/2).attr("y",-22)
    .attr("text-anchor","middle").attr("font-size","10px").attr("font-weight","bold")
    .text("B: Temporal pathway response (sig. pathways per cell type)");
  gB.append("text").attr("x",nC*cW/2).attr("y",-8)
    .attr("text-anchor","middle").attr("font-size","8.5px").attr("fill","#555")
    .text("Hallmark fGSEA padj < 0.05 · 5 timepoint comparisons × 10 cell types");

  // Color scale: 0=white, max=deep blue
  const colorN = d3.scaleSequential().domain([0,maxN]).interpolator(d3.interpolateBlues);

  // Phase separator y position (after row 2 = LP3/FP3)
  const sepY = 3*cH;

  // Column headers
  gB.selectAll(".colH").data(D.ct_labels).join("text")
    .attr("x",(d,i)=>i*cW+cW/2).attr("y",-4)
    .attr("text-anchor","end").attr("font-size","8px")
    .attr("transform",(d,i)=>`rotate(-38,${{i*cW+cW/2}},-4)`)
    .text(d=>d);

  // Phase labels (left)
  gB.append("text").attr("x",-6).attr("y",1.5*cH).attr("text-anchor","end")
    .attr("font-size","8.5px").attr("fill",C.blue).attr("dominant-baseline","middle")
    .text("vs. pre-flight");
  gB.append("text").attr("x",-6).attr("y",sepY+cH).attr("text-anchor","end")
    .attr("font-size","8.5px").attr("fill",C.orange).attr("dominant-baseline","middle")
    .text("recovery");

  // Divider line
  gB.append("line").attr("x1",-4).attr("x2",nC*cW+4)
    .attr("y1",sepY-2).attr("y2",sepY-2).attr("stroke","#aaa").attr("stroke-dasharray","4,3");

  // Row labels
  gB.selectAll(".rowH").data(D.sheet_labels).join("text")
    .attr("x",-5).attr("y",(d,i)=>i*cH+cH/2+3.5)
    .attr("text-anchor","end").attr("font-size","8.5px").text(d=>d);

  // Cells
  const cellG = gB.selectAll(".cG").data(D.cells_sig).join("g")
    .attr("transform",d=>`translate(${{d.col*cW}},${{d.row*cH}})`);
  cellG.append("rect").attr("width",cW).attr("height",cH)
    .attr("fill",d=>d.n===0?"#f5f5f5":colorN(d.n))
    .attr("stroke","#fff").attr("stroke-width",0.8)
    .on("mouseover",(ev,d)=>tip.style("display","block").html(
      `<b>${{D.ct_labels[d.col]}}</b><br><b>${{D.sheet_labels[d.row]}}</b><br>Sig pathways: ${{d.n}}`))
    .on("mousemove",ev=>tip.style("left",(ev.pageX+12)+"px").style("top",(ev.pageY-20)+"px"))
    .on("mouseout",()=>tip.style("display","none"));
  cellG.filter(d=>d.n>0).append("text").attr("x",cW/2).attr("y",cH/2+4)
    .attr("text-anchor","middle").attr("font-size","9px")
    .attr("fill",d=>d.n>=maxN*0.55?"#fff":"#222").text(d=>d.n);

  // Legend
  const legG = gB.append("g").attr("transform",`translate(0,${{nR*cH+18}})`);
  const legW = 110, legH = 10;
  const grad2 = svg.select("defs").append("linearGradient").attr("id","sigGrad").attr("x1","0%").attr("x2","100%");
  d3.range(0,1.01,0.1).forEach(t=>{{
    grad2.append("stop").attr("offset",`${{t*100}}%`).attr("stop-color",colorN(t*maxN));
  }});
  legG.append("rect").attr("width",legW).attr("height",legH).style("fill","url(#sigGrad)");
  legG.append("text").attr("y",-3).attr("font-size","8.5px").text("Sig pathways (padj<0.05)");
  [0, Math.round(maxN/2), maxN].forEach(v=>{{
    const x = v/maxN*legW;
    legG.append("line").attr("x1",x).attr("x2",x).attr("y1",legH).attr("y2",legH+4).attr("stroke","#666");
    legG.append("text").attr("x",x).attr("y",legH+12).attr("text-anchor","middle").attr("font-size","7.5px").text(v);
  }});
}})();
</script>
"""
    html += HTML_FOOTER

    out = FIGS / "Fig3_pbmc_celltype.html"
    out.write_text(html, encoding="utf-8")
    print(f"Saved: {out}")


# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import os
    os.chdir(Path(__file__).parent.parent.parent)  # → GeneLab_benchmark/

    print("Generating main figures...")
    make_fig1()
    make_fig2()
    make_fig3()
    print("\nDone. All 3 figures saved to v2/figures/")

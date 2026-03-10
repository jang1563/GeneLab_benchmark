#!/usr/bin/env python3
"""
E3: cfRNA cellular origin analysis
Q: Which I4 PBMC cell type's pathway response best predicts plasma cfRNA NES?
Data: E2 cfRNA fGSEA (48 pathways) + F1 snRNA-seq cell-type fGSEA (10 cell types × 34 pathways)
Output: v2/evaluation/E3_cfrna_origin.json, v2/figures/E3_cfrna_origin.html
"""

import json
import math
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

np.random.seed(42)

# ── Load data ──────────────────────────────────────────────────────────────
cfrna_df = pd.read_csv("v2/processed/E2_human/i4_cfrna_hallmark_fgsea.csv")
ct_df    = pd.read_csv("v2/processed/F1_scrna/i4_snrnaseq_celltype_fgsea.csv")

cfrna_nes = cfrna_df.set_index("pathway")["NES"]
ct_pivot  = ct_df.pivot_table(index="pathway", columns="Cell_Type", values="NES")

common_pw = sorted(set(cfrna_nes.index) & set(ct_pivot.index))
print(f"Common pathways: {len(common_pw)}")

# ── Spearman correlation with bootstrap CI ─────────────────────────────────
def spearman_with_ci(x, y, n_boot=1000, n_perm=10000, seed=42):
    rng = np.random.default_rng(seed)
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    n = len(x)
    if n < 8:
        return None
    r, _ = stats.spearmanr(x, y)
    boots = [stats.spearmanr(
        rng.choice(x, n, replace=True),
        rng.choice(y, n, replace=True)
    ).statistic for _ in range(n_boot)]
    ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])
    perms = [stats.spearmanr(rng.permutation(x), y).statistic for _ in range(n_perm)]
    p_perm = (np.sum(np.abs(perms) >= abs(r)) + 1) / (n_perm + 1)
    return {"r": round(float(r), 4), "ci_lo": round(float(ci_lo), 4),
            "ci_hi": round(float(ci_hi), 4), "p": round(float(p_perm), 4), "n": n}

# ── Run per cell type ──────────────────────────────────────────────────────
results = {}
cfrna_vals = cfrna_nes.reindex(common_pw).values.astype(float)

CT_ORDER = ["CD14+ Monocyte", "CD16+ Monocyte", "Dendritic Cell",
            "Natural Killer Cell", "B Cell", "CD4+ T Cell",
            "CD8+ T Cell", "Other T Cell", "Other", "PBMC Pseudobulk"]

print("\nE3: cfRNA vs. cell-type NES correlation (I4 FP1, R+1 vs. pre-flight)")
print(f"{'Cell Type':<25} {'r':>6}  {'95% CI':>16}  {'p':>8}  {'n':>3}")
print("-" * 70)

for ct in CT_ORDER:
    if ct not in ct_pivot.columns:
        continue
    ct_vals = ct_pivot[ct].reindex(common_pw).values.astype(float)
    res = spearman_with_ci(ct_vals, cfrna_vals)
    if res:
        results[ct] = res
        sig = "*" if res["p"] < 0.05 else " "
        print(f"  {ct:<23} {res['r']:>+.3f}  [{res['ci_lo']:+.3f}, {res['ci_hi']:+.3f}]  "
              f"p={res['p']:.4f}{sig}  n={res['n']}")

# Save JSON
eval_path = Path("v2/evaluation/E3_cfrna_origin.json")
eval_path.parent.mkdir(exist_ok=True)
eval_path.write_text(json.dumps({"results": results,
                                  "n_common_pathways": len(common_pw),
                                  "comparison": "I4_FP1_R1_vs_preflight"}, indent=2))
print(f"\nSaved: {eval_path}")

# ── HTML Figure ────────────────────────────────────────────────────────────
# All r are negative → use highest |r| (strongest anti-correlation)
best_ct = max(results, key=lambda c: abs(results[c]["r"]))
best_r  = results[best_ct]["r"]
print(f"\nStrongest anti-correlation: {best_ct} (r = {best_r:+.3f})")

# Scatter data for best cell type
ct_best_vals = ct_pivot[best_ct].reindex(common_pw)
cfrna_plot   = cfrna_nes.reindex(common_pw)

def short(p):
    return p.replace("HALLMARK_", "").replace("_", " ").title()

scatter_pts = []
for pw in common_pw:
    x = ct_best_vals[pw]
    y = cfrna_plot[pw]
    if not (math.isnan(float(x)) or math.isnan(float(y))):
        scatter_pts.append({"pw": short(pw), "x": round(float(x),3), "y": round(float(y),3)})

bar_data = [{"ct": ct, **results[ct]} for ct in CT_ORDER if ct in results]

fig_data = {
    "bar":      bar_data,
    "scatter":  scatter_pts,
    "best_ct":  best_ct,
    "best_r":   best_r,
    "n_common": len(common_pw),
}

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>E3: I4 cfRNA Cellular Origin</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; font-size: 12px; margin: 20px; background: #fff; }}
  h2   {{ font-size: 14px; font-weight: bold; margin-bottom: 4px; }}
  .subtitle {{ font-size: 10px; color: #555; margin-bottom: 16px; }}
  .bar-best {{ fill: #0072B2; }}
  .bar-other {{ fill: #88BBDD; }}
  .bar-neg  {{ fill: #CC79A7; }}
  .ci-line  {{ stroke: #333; stroke-width: 1.5px; }}
  .dot      {{ fill: #0072B2; opacity: 0.8; stroke: #fff; stroke-width: 0.5px; }}
  .dot:hover {{ opacity: 1; stroke: #333; stroke-width: 1px; cursor: pointer; }}
  .ref-line {{ stroke: #999; stroke-dasharray: 4,3; stroke-width: 1px; }}
  .fit-line {{ stroke: #0072B2; stroke-width: 1.5px; }}
  .panel-title {{ font-size: 11px; font-weight: bold; fill: #222; }}
  #tooltip {{
    position: absolute; background: rgba(0,0,0,0.78); color: #fff;
    padding: 5px 9px; border-radius: 4px; font-size: 10px;
    pointer-events: none; display: none;
  }}
</style>
</head>
<body>
<h2>E3 · I4 PBMC snRNA-seq vs. cfRNA: Pathway-Level Correlation by Cell Type</h2>
<div class="subtitle">
  Q: Which PBMC cell type's transcriptome best predicts plasma cfRNA pathway signals?<br>
  I4-FP1 (R+1 vs. pre-flight) · {len(common_pw)} common Hallmark pathways · Spearman r ± 95% CI
</div>
<div id="tooltip"></div>
<svg id="main"></svg>

<script>
const D = {json.dumps(fig_data)};

const barW = 420, barH = 220;
const scW  = 300, scH  = 260;
const mL = 160, mT = 40, mB = 50, mR = 20;
const sep = 50;
const totalW = mL + barW + sep + scW + mR + 30;
const totalH = mT + Math.max(barH, scH) + mB + 20;

const svg = d3.select("#main").attr("width", totalW).attr("height", totalH);
const tip = d3.select("#tooltip");

// ── Panel A: horizontal bar chart ────────────────────────────────────────
const bgA = svg.append("g").attr("transform", `translate(${{mL}},${{mT}})`);
svg.append("text").attr("class","panel-title").attr("x", 8).attr("y", mT-8).text("A");

const cts  = D.bar.map(d => d.ct);
const yA   = d3.scaleBand().domain(cts).range([0, barH]).padding(0.3);
const xA   = d3.scaleLinear().domain([-0.6, 0.8]).range([0, barW]);
const xRef = xA(0);

// Zero line
bgA.append("line").attr("class","ref-line")
   .attr("x1", xRef).attr("x2", xRef).attr("y1", 0).attr("y2", barH);

// Bars
bgA.selectAll(".bar")
  .data(D.bar).join("rect")
    .attr("class", d => d.r >= 0 ? (d.ct === D.best_ct ? "bar-best" : "bar-other") : "bar-neg")
    .attr("x", d => d.r >= 0 ? xRef : xA(d.r))
    .attr("y", d => yA(d.ct))
    .attr("width", d => Math.abs(xA(d.r) - xRef))
    .attr("height", yA.bandwidth())
    .on("mouseover", (ev, d) => {{
      tip.style("display","block")
         .html(`<b>${{d.ct}}</b><br>r = ${{d.r.toFixed(3)}}<br>95% CI: [${{d.ci_lo.toFixed(3)}}, ${{d.ci_hi.toFixed(3)}}]<br>p = ${{d.p.toFixed(4)}}<br>n = ${{d.n}} pathways`);
    }})
    .on("mousemove", ev => tip.style("left",(ev.pageX+12)+"px").style("top",(ev.pageY-20)+"px"))
    .on("mouseout",  ()  => tip.style("display","none"));

// CI lines
bgA.selectAll(".ci").data(D.bar).join("line")
    .attr("class","ci-line")
    .attr("x1", d => xA(d.ci_lo)).attr("x2", d => xA(d.ci_hi))
    .attr("y1", d => yA(d.ct) + yA.bandwidth()/2)
    .attr("y2", d => yA(d.ct) + yA.bandwidth()/2);

// r value labels
bgA.selectAll(".rlab").data(D.bar).join("text")
    .attr("x", d => d.r >= 0 ? xA(d.ci_hi) + 4 : xA(d.ci_lo) - 4)
    .attr("y", d => yA(d.ct) + yA.bandwidth()/2 + 3.5)
    .attr("text-anchor", d => d.r >= 0 ? "start" : "end")
    .attr("font-size","9px")
    .attr("fill", d => d.p < 0.05 ? "#000" : "#888")
    .text(d => `${{d.r >= 0 ? "+" : ""}}${{d.r.toFixed(3)}}${{d.p < 0.05 ? "*" : ""}}`);

// Axes
bgA.append("g").call(d3.axisLeft(yA).tickSize(3)).selectAll("text").attr("font-size","9px");
bgA.append("g").attr("transform",`translate(0,${{barH}})`).call(d3.axisBottom(xA).ticks(6))
   .selectAll("text").attr("font-size","8px");
bgA.append("text").attr("x", barW/2).attr("y", barH+38).attr("text-anchor","middle")
   .attr("font-size","10px").text("Spearman r (cell-type NES vs. cfRNA NES)");

// ── Panel B: scatter for best cell type ──────────────────────────────────
const offX = mL + barW + sep;
const scG  = svg.append("g").attr("transform", `translate(${{offX}}, ${{mT}})`);
svg.append("text").attr("class","panel-title").attr("x", offX - 12).attr("y", mT-8).text("B");

const allX = D.scatter.map(d=>d.x), allY = D.scatter.map(d=>d.y);
const xPad = 0.3, yPad = 0.3;
const xS = d3.scaleLinear().domain([d3.min(allX)-xPad, d3.max(allX)+xPad]).range([0, scW]);
const yS = d3.scaleLinear().domain([d3.min(allY)-yPad, d3.max(allY)+yPad]).range([scH, 0]);

// Ref lines at 0
scG.append("line").attr("class","ref-line")
   .attr("x1",xS(0)).attr("x2",xS(0)).attr("y1",0).attr("y2",scH);
scG.append("line").attr("class","ref-line")
   .attr("x1",0).attr("x2",scW).attr("y1",yS(0)).attr("y2",yS(0));

// Fit line (OLS for visual)
const n = D.scatter.length;
const meanX = d3.mean(allX), meanY = d3.mean(allY);
const slope = d3.sum(allX.map((x,i)=>(x-meanX)*(allY[i]-meanY))) /
              d3.sum(allX.map(x=>(x-meanX)**2));
const intercept = meanY - slope*meanX;
const x0 = d3.min(allX)-xPad, x1 = d3.max(allX)+xPad;
scG.append("line").attr("class","fit-line")
   .attr("x1",xS(x0)).attr("x2",xS(x1))
   .attr("y1",yS(slope*x0+intercept)).attr("y2",yS(slope*x1+intercept));

// Dots
scG.selectAll(".dot").data(D.scatter).join("circle")
    .attr("class","dot").attr("r", 5)
    .attr("cx", d => xS(d.x)).attr("cy", d => yS(d.y))
    .on("mouseover", (ev, d) => {{
      tip.style("display","block")
         .html(`<b>${{d.pw}}</b><br>Cell-type NES: ${{d.x.toFixed(3)}}<br>cfRNA NES: ${{d.y.toFixed(3)}}`);
    }})
    .on("mousemove", ev => tip.style("left",(ev.pageX+12)+"px").style("top",(ev.pageY-20)+"px"))
    .on("mouseout",  ()  => tip.style("display","none"));

// Axes
scG.append("g").attr("transform",`translate(0,${{scH}})`).call(d3.axisBottom(xS).ticks(5))
   .selectAll("text").attr("font-size","8px");
scG.append("g").call(d3.axisLeft(yS).ticks(5)).selectAll("text").attr("font-size","8px");

scG.append("text").attr("x", scW/2).attr("y", scH+38).attr("text-anchor","middle")
   .attr("font-size","9px").text(`${{D.best_ct}} NES`);
scG.append("text").attr("transform","rotate(-90)").attr("x",-scH/2).attr("y",-38)
   .attr("text-anchor","middle").attr("font-size","9px").text("I4 cfRNA NES");

// r annotation
const bestRes = D.bar.find(d=>d.ct===D.best_ct);
scG.append("text").attr("x",4).attr("y",12).attr("font-size","9px").attr("fill","#0072B2")
   .text(`r = ${{D.best_r.toFixed(3)}}${{bestRes.p<0.05?"*":""}}  (n=${{bestRes.n}})`);
</script>
</body>
</html>
"""

out_html = Path("v2/figures/E3_cfrna_origin.html")
out_html.write_text(html, encoding="utf-8")
print(f"Saved: {out_html}")

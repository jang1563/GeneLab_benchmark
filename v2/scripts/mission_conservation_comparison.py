#!/usr/bin/env python3
"""
mission_conservation_comparison.py — GeneLab_benchmark v2, Category E2

Cross-mission duration effect on NES conservation:
  Compares how well short-duration (I4, 3d) vs. long-duration (JAXA, 120d)
  human spaceflight cfRNA NES correlates with mouse liver NES.

Inputs:
  I4 NES:   v2/processed/E2_human/i4_cfrna_hallmark_fgsea.csv
  JAXA NES: v2/processed/E1_human/human_cfrna_hallmark_fgsea.csv
  Mouse NES: v2/evaluation/E1_crossspecies_nes.json (mission-averaged)

Outputs:
  v2/evaluation/E2_mission_conservation.json
  v2/figures/E2_duration_conservation.html  (2-panel: scatter + bar chart)

Usage:
  python3 v2/scripts/mission_conservation_comparison.py
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

# ── Paths ─────────────────────────────────────────────────────────────────────

REPO_ROOT       = Path(__file__).resolve().parent.parent.parent
I4_FGSEA_CSV    = REPO_ROOT / "v2/processed/E2_human/i4_cfrna_hallmark_fgsea.csv"
JAXA_FGSEA_CSV  = REPO_ROOT / "v2/processed/E1_human/human_cfrna_hallmark_fgsea.csv"
E1_JSON         = REPO_ROOT / "v2/evaluation/E1_crossspecies_nes.json"
OUT_DIR         = REPO_ROOT / "v2/evaluation"
FIG_DIR         = REPO_ROOT / "v2/figures"

OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)


# ── Utilities ─────────────────────────────────────────────────────────────────

def normalize_pathway_name(name: str) -> str:
    """Normalize Hallmark pathway name to HALLMARK_XXX format."""
    name = str(name).strip().upper().replace(" ", "_").replace("-", "_")
    if not name.startswith("HALLMARK_"):
        name = "HALLMARK_" + name
    return name


def load_nes_from_fgsea_csv(csv_path: Path, mission_label: str) -> dict:
    """Load NES values from fGSEA CSV output."""
    df = pd.read_csv(csv_path)
    nes_map = {}
    for _, row in df.iterrows():
        pname = normalize_pathway_name(str(row["pathway"]))
        nes_val = row.get("NES", np.nan)
        if not (isinstance(nes_val, float) and np.isnan(nes_val)):
            nes_map[pname] = float(nes_val)
    print(f"  {mission_label}: {len(nes_map)} pathways loaded")
    return nes_map


def spearman_with_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 1000,
                     n_perm: int = 10000, seed: int = 42) -> dict:
    """Spearman r with bootstrap 95% CI and permutation p-value."""
    r, p_scipy = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    n = len(x)

    # Bootstrap CI
    boot_r = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            r_b, _ = stats.spearmanr(x[idx], y[idx])
            boot_r.append(r_b)
        except Exception:
            pass
    boot_r = np.array(boot_r)

    # Permutation p
    perm_r = []
    for _ in range(n_perm):
        y_perm = rng.permutation(y)
        r_p, _ = stats.spearmanr(x, y_perm)
        perm_r.append(r_p)
    perm_r = np.array(perm_r)
    p_perm = float(np.mean(np.abs(perm_r) >= np.abs(r)))
    p_perm = max(p_perm, 1.0 / n_perm)

    sign_agree = float(np.mean(np.sign(x) == np.sign(y)))

    return {
        "spearman_r":             float(r),
        "ci_low":                 float(np.percentile(boot_r, 2.5)),
        "ci_high":                float(np.percentile(boot_r, 97.5)),
        "p_scipy":                float(p_scipy),
        "p_permutation":          p_perm,
        "sign_agreement_fraction": sign_agree,
        "n_pathways":             int(len(x)),
    }


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=== E2: Cross-Mission Duration Conservation ===")

    # ── 1. Load NES values ─────────────────────────────────────────────────────
    print("\n[1] Loading NES values...")

    i4_nes   = load_nes_from_fgsea_csv(I4_FGSEA_CSV,   "I4 (3d)")
    jaxa_nes = load_nes_from_fgsea_csv(JAXA_FGSEA_CSV, "JAXA (120d)")

    # Mouse NES (mission-averaged) from E1 JSON
    with open(E1_JSON) as f:
        e1 = json.load(f)
    mouse_nes = {
        normalize_pathway_name(p["pathway"]): p["mouse_nes_mean"]
        for p in e1["pathways"]
    }
    print(f"  Mouse (mission-averaged): {len(mouse_nes)} pathways")

    # JAXA results already in E1 JSON
    jaxa_r_e1   = e1["results"]["mission_averaged"]["spearman_r"]
    jaxa_ci_low = e1["results"]["mission_averaged"]["ci_low"]
    jaxa_ci_hi  = e1["results"]["mission_averaged"]["ci_high"]
    jaxa_p_perm = e1["results"]["mission_averaged"]["p_permutation"]
    jaxa_n      = e1["results"]["mission_averaged"]["n_pathways"]

    # ── 2. Align I4 vs mouse ───────────────────────────────────────────────────
    print("\n[2] Aligning I4 vs mouse NES...")
    common_i4 = sorted(set(i4_nes.keys()) & set(mouse_nes.keys()))
    print(f"  Common pathways (I4 ∩ mouse): {len(common_i4)}")

    i4_vec    = np.array([i4_nes[p]   for p in common_i4])
    mouse_i4  = np.array([mouse_nes[p] for p in common_i4])

    # ── 3. Align JAXA vs mouse (for scatter, may differ from E1 due to I4 pathway subset) ──
    common_jaxa = sorted(set(jaxa_nes.keys()) & set(mouse_nes.keys()))
    jaxa_vec    = np.array([jaxa_nes[p]  for p in common_jaxa])
    mouse_jaxa  = np.array([mouse_nes[p] for p in common_jaxa])

    # ── 4. Compute I4 Spearman r ───────────────────────────────────────────────
    print("\n[3] Computing I4 vs mouse Spearman r...")
    i4_stats = spearman_with_ci(i4_vec, mouse_i4)
    print(f"  I4:   r = {i4_stats['spearman_r']:.3f} "
          f"(95% CI: {i4_stats['ci_low']:.3f}–{i4_stats['ci_high']:.3f}), "
          f"p = {i4_stats['p_permutation']:.3f}")
    print(f"  JAXA: r = {jaxa_r_e1:.3f} "
          f"(95% CI: {jaxa_ci_low:.3f}–{jaxa_ci_hi:.3f}), "
          f"p = {jaxa_p_perm:.3f}  [from E1]")

    # ── 5. Pathway-level details ───────────────────────────────────────────────
    pathway_details_i4 = [
        {
            "pathway":       p,
            "i4_nes":        float(i4_nes[p]),
            "mouse_nes":     float(mouse_nes[p]),
            "concordant":    bool(np.sign(i4_nes[p]) == np.sign(mouse_nes[p])),
        }
        for p in common_i4
    ]

    # ── 6. Save JSON results ───────────────────────────────────────────────────
    print("\n[4] Saving JSON results...")
    out_json = OUT_DIR / "E2_mission_conservation.json"
    result = {
        "task": "E2",
        "description": "Cross-mission duration effect on NES conservation (I4 3d vs. JAXA 120d)",
        "missions": {
            "I4": {
                "label": "Inspiration4",
                "duration_days": 3,
                "crew_size": 4,
                "data_source": "cfrna_crossmission_r1.csv (log2FoldChange_I4)",
                "n_genes_ranked": len(i4_nes),
            },
            "JAXA": {
                "label": "JAXA CFE (OSD-530)",
                "duration_days": 120,
                "crew_size": 6,
                "data_source": "human_cfrna_hallmark_fgsea.csv",
                "n_genes_ranked": 26845,
            },
        },
        "mouse_reference": {
            "tissue": "liver",
            "missions": ["MHU-2", "RR-1", "RR-3", "RR-6", "RR-8", "RR-9"],
            "averaging": "arithmetic mean across missions",
        },
        "results": {
            "I4_vs_mouse": {**i4_stats, "n_pathways": len(common_i4)},
            "JAXA_vs_mouse": {
                "spearman_r": jaxa_r_e1,
                "ci_low":     jaxa_ci_low,
                "ci_high":    jaxa_ci_hi,
                "p_permutation": jaxa_p_perm,
                "n_pathways": jaxa_n,
                "source": "E1_crossspecies_nes.json",
            },
        },
        "pathways_i4": pathway_details_i4,
    }

    with open(out_json, "w") as f:
        json.dump(result, f, indent=2)
    print(f"  Saved: {out_json}")

    # ── 7. Generate 2-panel HTML figure ───────────────────────────────────────
    print("\n[5] Generating 2-panel HTML figure...")
    fig_path = FIG_DIR / "E2_duration_conservation.html"
    _make_figure(
        i4_vec=i4_vec,
        mouse_i4=mouse_i4,
        pathways_i4=common_i4,
        i4_stats=i4_stats,
        jaxa_r=jaxa_r_e1,
        jaxa_ci_low=jaxa_ci_low,
        jaxa_ci_hi=jaxa_ci_hi,
        jaxa_p=jaxa_p_perm,
        out_path=fig_path,
    )
    print(f"  Saved: {fig_path}")

    # ── Summary ────────────────────────────────────────────────────────────────
    print("\n=== E2 Summary ===")
    print(f"  I4  (3d):   r = {i4_stats['spearman_r']:.3f}  "
          f"(CI: {i4_stats['ci_low']:.3f}–{i4_stats['ci_high']:.3f}), "
          f"p = {i4_stats['p_permutation']:.3f}")
    print(f"  JAXA(120d): r = {jaxa_r_e1:.3f}  "
          f"(CI: {jaxa_ci_low:.3f}–{jaxa_ci_hi:.3f}), "
          f"p = {jaxa_p_perm:.3f}")
    delta = jaxa_r_e1 - i4_stats["spearman_r"]
    print(f"  Δr (JAXA − I4) = {delta:.3f}")


# ── Figure ────────────────────────────────────────────────────────────────────

def _make_figure(i4_vec, mouse_i4, pathways_i4, i4_stats,
                 jaxa_r, jaxa_ci_low, jaxa_ci_hi, jaxa_p, out_path):
    """Generate self-contained 2-panel HTML (D3.js v7, Okabe-Ito)."""

    # Panel 1 data: I4 scatter
    scatter_data = [
        {
            "pathway": p,
            "i4_nes": float(h),
            "mouse_nes": float(m),
            "concordant": bool(np.sign(h) == np.sign(m)),
        }
        for p, h, m in zip(pathways_i4, i4_vec, mouse_i4)
    ]

    # Panel 2 data: duration bar chart
    bar_data = [
        {
            "label": "I4 (3d)",
            "r": i4_stats["spearman_r"],
            "ci_low": i4_stats["ci_low"],
            "ci_high": i4_stats["ci_high"],
            "p": i4_stats["p_permutation"],
            "duration": 3,
        },
        {
            "label": "JAXA (120d)",
            "r": jaxa_r,
            "ci_low": jaxa_ci_low,
            "ci_high": jaxa_ci_hi,
            "p": jaxa_p,
            "duration": 120,
        },
    ]

    i4_r = i4_stats["spearman_r"]
    i4_ci_low = i4_stats["ci_low"]
    i4_ci_hi  = i4_stats["ci_high"]
    i4_p = i4_stats["p_permutation"]
    i4_n  = i4_stats["n_pathways"]
    jaxa_n = 50  # from E1

    i4_p_str   = f"p = {i4_p:.3f}"   if i4_p   >= 0.001 else "p < 0.001"
    jaxa_p_str = f"p = {jaxa_p:.3f}" if jaxa_p >= 0.001 else "p < 0.001"

    scatter_json = json.dumps(scatter_data)
    bar_json     = json.dumps(bar_data)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>E2 Cross-Mission Duration Conservation</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; font-size: 12px; background: #fff; margin: 20px; }}
  .axis text {{ font-size: 9px; }}
  .axis line, .axis path {{ stroke: #333; }}
  .gridline {{ stroke: #e0e0e0; stroke-width: 0.5; }}
  .dot {{ opacity: 0.8; cursor: pointer; }}
  .dot:hover {{ opacity: 1; stroke: #333; stroke-width: 1; }}
  .tooltip {{
    position: absolute; background: rgba(255,255,255,0.95);
    border: 1px solid #ccc; padding: 6px; border-radius: 3px;
    font-size: 10px; pointer-events: none; display: none;
  }}
  .panel-title {{ font-size: 12px; font-weight: bold; }}
  .panel-sub {{ font-size: 9px; fill: #555; }}
  svg {{ display: inline-block; vertical-align: top; }}
  .bar {{ opacity: 0.85; }}
  .bar:hover {{ opacity: 1; }}
  .errbar {{ stroke: #333; stroke-width: 1.5; }}
  .zeroline {{ stroke: #999; stroke-width: 0.8; stroke-dasharray: 4,3; }}
  .sig-star {{ font-size: 11px; fill: #333; }}
</style>
</head>
<body>
<h2 style="font-size:14px;margin-bottom:4px">E2: Cross-Mission Duration Effect on cfRNA–Mouse NES Conservation</h2>
<div id="tooltip" class="tooltip"></div>

<script>
const scatterData = {scatter_json};
const barData     = {bar_json};

// ─── Okabe-Ito palette ───────────────────────────────────────────────────────
const OKI = {{
  concordant:   "#0072B2",  // blue
  discordant:   "#E69F00",  // orange
  i4_bar:       "#56B4E9",  // sky blue
  jaxa_bar:     "#0072B2",  // deep blue
  zero:         "#999999",
}};

const tooltip = d3.select("#tooltip");

// ─── Panel 1: I4 scatter ─────────────────────────────────────────────────────
(function() {{
  const W = 380, H = 360;
  const margin = {{top:50, right:20, bottom:55, left:60}};
  const w = W - margin.left - margin.right;
  const h = H - margin.top  - margin.bottom;

  const svg = d3.select("body")
    .append("svg").attr("width", W).attr("height", H);
  const g = svg.append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const allX = scatterData.map(d => d.i4_nes);
  const allY = scatterData.map(d => d.mouse_nes);
  const xExt = d3.extent(allX), yExt = d3.extent(allY);
  const pad = 0.3;

  const xScale = d3.scaleLinear().domain([xExt[0]-pad, xExt[1]+pad]).range([0, w]).nice();
  const yScale = d3.scaleLinear().domain([yExt[0]-pad, yExt[1]+pad]).range([h, 0]).nice();

  // Gridlines
  g.append("g").attr("class","gridline")
    .call(d3.axisLeft(yScale).ticks(6).tickSize(-w).tickFormat(""))
    .selectAll("line").style("stroke","#e8e8e8");
  g.append("g").attr("class","gridline")
    .attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(xScale).ticks(6).tickSize(-h).tickFormat(""))
    .selectAll("line").style("stroke","#e8e8e8");

  // Zero lines
  g.append("line").attr("class","zeroline")
    .attr("x1",xScale(0)).attr("x2",xScale(0))
    .attr("y1",0).attr("y2",h);
  g.append("line").attr("class","zeroline")
    .attr("x1",0).attr("x2",w)
    .attr("y1",yScale(0)).attr("y2",yScale(0));

  // Axes
  g.append("g").attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(xScale).ticks(6));
  g.append("g").call(d3.axisLeft(yScale).ticks(6));

  // Axis labels
  g.append("text").attr("text-anchor","middle")
    .attr("x", w/2).attr("y", h+40)
    .style("font-size","10px")
    .text("I4 cfRNA NES (R+1 vs. pre-flight)");
  g.append("text").attr("text-anchor","middle")
    .attr("transform","rotate(-90)")
    .attr("x", -h/2).attr("y", -45)
    .style("font-size","10px")
    .text("Mouse Liver NES (mission-averaged)");

  // Dots
  g.selectAll(".dot")
    .data(scatterData)
    .enter().append("circle")
      .attr("class","dot")
      .attr("cx", d => xScale(d.i4_nes))
      .attr("cy", d => yScale(d.mouse_nes))
      .attr("r", 5)
      .attr("fill", d => d.concordant ? OKI.concordant : OKI.discordant)
      .on("mouseover", function(event, d) {{
        tooltip.style("display","block")
          .html(`<b>${{d.pathway.replace("HALLMARK_","").replace(/_/g," ")}}</b><br>
                 I4 NES: ${{d.i4_nes.toFixed(3)}}<br>
                 Mouse NES: ${{d.mouse_nes.toFixed(3)}}<br>
                 Sign: ${{d.concordant ? "Concordant" : "Discordant"}}`);
      }})
      .on("mousemove", function(event) {{
        tooltip.style("left",(event.pageX+10)+"px").style("top",(event.pageY-20)+"px");
      }})
      .on("mouseout", () => tooltip.style("display","none"));

  // Legend
  const leg = g.append("g").attr("transform",`translate(${{w-130}},5)`);
  leg.append("circle").attr("r",5).attr("fill",OKI.concordant).attr("cy",0);
  leg.append("text").attr("x",9).attr("y",4).style("font-size","8px").text("Concordant sign");
  leg.append("circle").attr("r",5).attr("fill",OKI.discordant).attr("cy",14);
  leg.append("text").attr("x",9).attr("y",18).style("font-size","8px").text("Discordant sign");

  // Title
  svg.append("text").attr("class","panel-title")
    .attr("x", W/2).attr("y", 18).attr("text-anchor","middle")
    .text("I4 (3-day) vs. Mouse Liver NES");
  svg.append("text").attr("class","panel-sub")
    .attr("x", W/2).attr("y", 32).attr("text-anchor","middle")
    .text(`r = {i4_r:.3f} (CI: {i4_ci_low:.3f}–{i4_ci_hi:.3f}), {i4_p_str}, n = {i4_n}`);
}})();

// ─── Panel 2: Duration bar chart ─────────────────────────────────────────────
(function() {{
  const W = 300, H = 360;
  const margin = {{top:50, right:30, bottom:80, left:70}};
  const w = W - margin.left - margin.right;
  const h = H - margin.top  - margin.bottom;

  const svg = d3.select("body")
    .append("svg").attr("width", W).attr("height", H);
  const g = svg.append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const labels = barData.map(d => d.label);
  const xScale = d3.scaleBand().domain(labels).range([0, w]).padding(0.45);

  const allR = barData.flatMap(d => [d.ci_low, d.ci_high]);
  const yMin = Math.min(0, d3.min(allR)) - 0.05;
  const yMax = Math.max(0, d3.max(allR)) + 0.05;
  const yScale = d3.scaleLinear().domain([yMin, yMax]).range([h, 0]).nice();

  // Gridlines
  g.append("g").attr("class","gridline")
    .call(d3.axisLeft(yScale).ticks(6).tickSize(-w).tickFormat(""))
    .selectAll("line").style("stroke","#e8e8e8");

  // Zero line
  g.append("line").attr("class","zeroline")
    .attr("x1",0).attr("x2",w)
    .attr("y1",yScale(0)).attr("y2",yScale(0));

  // Axes
  g.append("g").attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(xScale))
    .selectAll("text").style("font-size","10px");
  g.append("g").call(d3.axisLeft(yScale).ticks(6));

  // Y-axis label
  g.append("text").attr("text-anchor","middle")
    .attr("transform","rotate(-90)")
    .attr("x", -h/2).attr("y", -55)
    .style("font-size","10px")
    .text("Spearman r (vs. Mouse Liver NES)");

  // Bars
  const colors = [OKI.i4_bar, OKI.jaxa_bar];
  g.selectAll(".bar")
    .data(barData)
    .enter().append("rect")
      .attr("class","bar")
      .attr("x", d => xScale(d.label))
      .attr("y", d => d.r >= 0 ? yScale(d.r) : yScale(0))
      .attr("width", xScale.bandwidth())
      .attr("height", d => Math.abs(yScale(d.r) - yScale(0)))
      .attr("fill", (d,i) => colors[i]);

  // Error bars (95% CI)
  const errW = 6;
  barData.forEach(d => {{
    const cx = xScale(d.label) + xScale.bandwidth()/2;
    g.append("line").attr("class","errbar")
      .attr("x1",cx).attr("x2",cx)
      .attr("y1",yScale(d.ci_low)).attr("y2",yScale(d.ci_high));
    g.append("line").attr("class","errbar")
      .attr("x1",cx-errW).attr("x2",cx+errW)
      .attr("y1",yScale(d.ci_low)).attr("y2",yScale(d.ci_low));
    g.append("line").attr("class","errbar")
      .attr("x1",cx-errW).attr("x2",cx+errW)
      .attr("y1",yScale(d.ci_high)).attr("y2",yScale(d.ci_high));
  }});

  // p-value annotations
  barData.forEach(d => {{
    const cx = xScale(d.label) + xScale.bandwidth()/2;
    const yTop = yScale(Math.max(d.r, d.ci_high)) - 8;
    const pStr = d.p < 0.05 ? (d.p < 0.01 ? "**" : "*") : "ns";
    g.append("text").attr("class","sig-star")
      .attr("x", cx).attr("y", yTop)
      .attr("text-anchor","middle")
      .text(pStr);
  }});

  // Duration labels below bars
  g.append("g").attr("transform",`translate(0,${{h+2}})`)
    .selectAll("text")
    .data(barData)
    .enter().append("text")
      .attr("x", d => xScale(d.label) + xScale.bandwidth()/2)
      .attr("y", 28)
      .attr("text-anchor","middle")
      .style("font-size","9px")
      .style("fill","#555")
      .text(d => `${{d.duration}}d`);

  // Title
  svg.append("text").attr("class","panel-title")
    .attr("x", W/2).attr("y", 18).attr("text-anchor","middle")
    .text("Duration vs. Mouse NES Conservation");
  svg.append("text").attr("class","panel-sub")
    .attr("x", W/2).attr("y", 33).attr("text-anchor","middle")
    .text("* p<0.05, ** p<0.01, ns = not significant");

  // Download button
  svg.append("foreignObject").attr("x", W-60).attr("y", H-22).attr("width",55).attr("height",20)
    .append("xhtml:button")
    .style("font-size","8px").style("cursor","pointer")
    .text("Save SVG")
    .on("click", function() {{
      const svgEl = svg.node();
      const serializer = new XMLSerializer();
      const blob = new Blob([serializer.serializeToString(svgEl)], {{type:"image/svg+xml"}});
      const a = document.createElement("a");
      a.href = URL.createObjectURL(blob);
      a.download = "E2_duration_bar.svg";
      a.click();
    }});
}})();
</script>
</body>
</html>
"""
    with open(out_path, "w") as f:
        f.write(html)


if __name__ == "__main__":
    main()

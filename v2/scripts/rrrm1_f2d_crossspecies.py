#!/usr/bin/env python3
"""
rrrm1_f2d_crossspecies.py — F2-D: Cross-species Pathway Concordance

RRRM-1 mouse blood cell-type fGSEA NES vs I4 human PBMC cell-type fGSEA NES.

Question: Does cell-type matching improve cross-species pathway concordance
          beyond E1 bulk baseline (r=0.352)?

Approach:
  1. Load RRRM-1 blood per-cell-type NES (from F2-B output)
  2. Load I4 PBMC per-cell-type NES (from F1 output)
  3. Compute Spearman r for matched cell-type pairs
  4. Compute all-pairs matrix (5 RRRM-1 × 10 I4)
  5. Bootstrap 95% CI + permutation p-value (same as E1)
  6. Compare monocyte pair r vs E1 bulk baseline r=0.352

Sign convention:
  RRRM-1: FLT vs GC → positive NES = spaceflight enriched
  I4:     FP1_R1 vs preflight → positive NES = post-flight enriched
  Both represent spaceflight effect → directly comparable.

Usage:
  python3 v2/scripts/rrrm1_f2d_crossspecies.py

Output:
  v2/evaluation/F2D_crossspecies.json
  v2/figures/F2D_crossspecies_scatter.html
"""

import json
import sys
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────

REPO_ROOT = Path(__file__).resolve().parent.parent.parent  # GeneLab_benchmark/

RRRM1_BLOOD_DIR = REPO_ROOT / "v2/processed/F2B_blood"
I4_FGSEA_CSV = REPO_ROOT / "v2/processed/F1_scrna/i4_snrnaseq_celltype_fgsea.csv"
E1_JSON = REPO_ROOT / "v2/evaluation/E1_crossspecies_nes.json"

OUT_DIR = REPO_ROOT / "v2/evaluation"
FIG_DIR = REPO_ROOT / "v2/figures"

# Cell type mapping: (I4 human, RRRM-1 mouse, lineage label)
MATCHED_PAIRS = [
    ("CD14+ Monocyte", "monocyte_macrophage", "monocyte"),
    ("B Cell", "b_cell", "B cell"),
    ("CD4+ T Cell", "t_cell", "T cell"),
]

SENSITIVITY_PAIRS = [
    ("CD16+ Monocyte", "monocyte_macrophage", "monocyte (CD16+)"),
    ("CD8+ T Cell", "t_cell", "T cell (CD8+)"),
    ("Other T Cell", "t_cell", "T cell (Other)"),
]

MIN_PATHWAYS = 8  # Minimum pathway overlap for meaningful Spearman

# Okabe-Ito palette
CONCORDANT_COLOR = "#0072B2"
DISCORDANT_COLOR = "#E69F00"


# ── Data Loading ──────────────────────────────────────────────────────────────

def load_rrrm1_blood_nes(data_dir: Path) -> dict:
    """Load RRRM-1 blood per-cell-type NES from F2-B CSVs."""
    nes_dict = {}
    for csv_path in sorted(data_dir.glob("*_fgsea_hallmark.csv")):
        ct = csv_path.stem.replace("_fgsea_hallmark", "")
        df = pd.read_csv(csv_path)
        nes_dict[ct] = df.set_index("pathway")["NES"]
        print(f"  RRRM-1 {ct}: {len(nes_dict[ct])} pathways")
    return nes_dict


def load_i4_pbmc_nes(csv_path: Path) -> dict:
    """Load I4 PBMC per-cell-type NES from F1 fGSEA CSV."""
    df = pd.read_csv(csv_path)
    df = df[df["comparison"] == "FP1_R1_vs_preflight"]

    nes_dict = {}
    for ct in sorted(df["Cell_Type"].unique()):
        sub = df[df["Cell_Type"] == ct]
        nes = sub.set_index("pathway")["NES"]
        nes_dict[ct] = nes
        print(f"  I4 {ct}: {len(nes)} pathways")
    return nes_dict


def load_e1_baseline(json_path: Path) -> dict:
    """Load E1 baseline results."""
    with open(json_path) as f:
        e1 = json.load(f)
    r = e1["results"]["mission_averaged"]["spearman_r"]
    n = e1["results"]["mission_averaged"]["n_pathways"]
    print(f"  E1 baseline: r={r:.3f}, n={n} pathways")
    return {"spearman_r": r, "n_pathways": n}


# ── Statistics (same method as E1: cross_species_nes_comparison.py) ───────────

def compute_concordance(nes_x: pd.Series, nes_y: pd.Series,
                        name_x: str, name_y: str) -> dict:
    """Compute Spearman r with bootstrap CI and permutation p between two NES vectors."""
    common = sorted(set(nes_x.index) & set(nes_y.index))
    n = len(common)

    if n < MIN_PATHWAYS:
        return {
            "spearman_r": np.nan, "ci_low": np.nan, "ci_high": np.nan,
            "p_scipy": np.nan, "p_permutation": np.nan,
            "sign_agreement": np.nan, "n_pathways": n,
            "skipped": True, "reason": f"n_pathways={n} < {MIN_PATHWAYS}",
        }

    x = nes_x.reindex(common).values
    y = nes_y.reindex(common).values

    r, p_scipy = stats.spearmanr(x, y)

    # Bootstrap 95% CI
    rng = np.random.default_rng(42)
    boot_r = []
    for _ in range(1000):
        idx = rng.integers(0, n, size=n)
        try:
            r_b, _ = stats.spearmanr(x[idx], y[idx])
            if not np.isnan(r_b):
                boot_r.append(r_b)
        except Exception:
            pass
    if boot_r:
        ci_low = float(np.percentile(boot_r, 2.5))
        ci_high = float(np.percentile(boot_r, 97.5))
    else:
        ci_low, ci_high = np.nan, np.nan

    # Permutation p-value (two-tailed, 10K)
    perm_r = []
    for _ in range(10000):
        r_p, _ = stats.spearmanr(x, rng.permutation(y))
        perm_r.append(r_p)
    perm_r = np.array(perm_r)
    p_perm = float(max(np.mean(np.abs(perm_r) >= np.abs(r)), 1 / 10000))

    # Sign agreement
    sign_agree = float(np.mean(np.sign(x) == np.sign(y)))

    # Per-pathway details
    pathways_detail = [
        {
            "pathway": p,
            "rrrm1_nes": float(nes_x[p]),
            "i4_nes": float(nes_y[p]),
            "concordant": bool(np.sign(nes_x[p]) == np.sign(nes_y[p])),
        }
        for p in common
    ]

    return {
        "spearman_r": float(r),
        "ci_low": ci_low,
        "ci_high": ci_high,
        "p_scipy": float(p_scipy),
        "p_permutation": p_perm,
        "sign_agreement": sign_agree,
        "n_pathways": n,
        "pathways": pathways_detail,
    }


def compute_all_pairs_matrix(rrrm1_nes: dict, i4_nes: dict) -> dict:
    """Compute Spearman r for all (RRRM-1 × I4) cell-type pairs."""
    rrrm1_cts = sorted(rrrm1_nes.keys())
    i4_cts = sorted(i4_nes.keys())

    r_matrix = []
    n_matrix = []
    for i4_ct in i4_cts:
        r_row = []
        n_row = []
        for rrrm1_ct in rrrm1_cts:
            common = sorted(set(rrrm1_nes[rrrm1_ct].index) & set(i4_nes[i4_ct].index))
            n = len(common)
            if n >= MIN_PATHWAYS:
                x = rrrm1_nes[rrrm1_ct].reindex(common).values
                y = i4_nes[i4_ct].reindex(common).values
                r, _ = stats.spearmanr(x, y)
                r_row.append(round(float(r), 4))
            else:
                r_row.append(None)
            n_row.append(n)
        r_matrix.append(r_row)
        n_matrix.append(n_row)

    return {
        "rrrm1_celltypes": rrrm1_cts,
        "i4_celltypes": i4_cts,
        "spearman_r": r_matrix,
        "n_pathways": n_matrix,
    }


# ── HTML Figure ───────────────────────────────────────────────────────────────

def make_f2d_html(matched_results: list, sensitivity_results: list,
                  matrix: dict, e1_baseline: dict, out_path: Path):
    """Generate F2-D HTML with scatter + heatmap + bar chart (D3.js v7)."""

    # Primary monocyte scatter data
    mono = matched_results[0]  # CD14+ Monocyte pair
    scatter_data = json.dumps(mono.get("pathways", []))
    mono_r = mono["spearman_r"]
    mono_ci = f"{mono['ci_low']:.3f}–{mono['ci_high']:.3f}"
    mono_p = mono["p_permutation"]
    mono_n = mono["n_pathways"]
    mono_sign = mono["sign_agreement"]
    p_str = f"p = {mono_p:.3f}" if mono_p >= 0.001 else "p < 0.001"

    # Bar chart data: matched + sensitivity + E1 baseline
    bar_data = []
    for m in matched_results:
        if not m.get("skipped"):
            bar_data.append({
                "label": m["matched_lineage"],
                "r": m["spearman_r"],
                "ci_low": m["ci_low"],
                "ci_high": m["ci_high"],
                "n": m["n_pathways"],
                "type": "matched",
            })
    for s in sensitivity_results:
        if not s.get("skipped"):
            bar_data.append({
                "label": s["matched_lineage"],
                "r": s["spearman_r"],
                "ci_low": s["ci_low"],
                "ci_high": s["ci_high"],
                "n": s["n_pathways"],
                "type": "sensitivity",
            })
    bar_json = json.dumps(bar_data)
    matrix_json = json.dumps(matrix)
    e1_r = e1_baseline["spearman_r"]

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>F2-D Cross-Species Cell-Type NES Concordance</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; font-size: 12px; background: #fff; margin: 20px; }}
  h2 {{ font-size: 14px; font-weight: bold; margin-bottom: 4px; }}
  .subtitle {{ font-size: 9px; color: #555; margin-bottom: 16px; }}
  .panel {{ display: inline-block; vertical-align: top; margin: 10px; }}
  .axis text {{ font-size: 9px; }}
  .axis line, .axis path {{ stroke: #333; }}
  .gridline {{ stroke: #e0e0e0; stroke-width: 0.5; }}
  .dot {{ opacity: 0.8; cursor: pointer; }}
  .dot:hover {{ opacity: 1; stroke: #333; stroke-width: 1.5; }}
  .tooltip {{
    position: absolute; background: rgba(255,255,255,0.95);
    border: 1px solid #ccc; padding: 6px; border-radius: 3px;
    font-size: 10px; pointer-events: none; display: none;
  }}
  .heatmap-cell {{ stroke: #fff; stroke-width: 1; }}
  .heatmap-cell.matched {{ stroke: #000; stroke-width: 2; }}
  .ref-line {{ stroke: #D55E00; stroke-width: 1.5; stroke-dasharray: 5,3; }}
  #download-btn {{
    margin-top: 10px; padding: 5px 12px; cursor: pointer;
    font-size: 11px; border: 1px solid #555; border-radius: 3px; background: #f5f5f5;
  }}
</style>
</head>
<body>
<h2>F2-D: Cross-Species Cell-Type Pathway Concordance (RRRM-1 Mouse Blood vs I4 Human PBMC)</h2>
<p class="subtitle">Spearman r on matched Hallmark NES vectors. E1 bulk baseline: r = {e1_r:.3f} (mouse liver vs. JAXA cfRNA, 50 pathways).</p>

<div class="tooltip" id="tooltip"></div>

<div>
  <div class="panel" id="scatter"></div>
  <div class="panel" id="barchart"></div>
</div>
<div>
  <div class="panel" id="heatmap"></div>
</div>

<button id="download-btn">Download SVG</button>

<script>
const scatterData = {scatter_data};
const barData = {bar_json};
const matrix = {matrix_json};
const e1_r = {e1_r};
const monoR = {mono_r if not np.isnan(mono_r) else 'null'};

// ─── Panel 1: Monocyte Scatter ───────────────────────────────────────
(function() {{
  const margin = {{top: 50, right: 30, bottom: 60, left: 70}};
  const W = 420, H = 400;
  const w = W - margin.left - margin.right;
  const h = H - margin.top - margin.bottom;

  const svg = d3.select("#scatter").append("svg")
    .attr("width", W).attr("height", H)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  if (!scatterData.length) {{
    svg.append("text").attr("x", w/2).attr("y", h/2)
      .attr("text-anchor","middle").attr("font-size","11px")
      .text("No overlapping pathways for monocyte pair");
    return;
  }}

  const xExt = d3.extent(scatterData, d => d.rrrm1_nes);
  const yExt = d3.extent(scatterData, d => d.i4_nes);
  const pad = 0.3;
  const xS = d3.scaleLinear().domain([xExt[0]-pad, xExt[1]+pad]).range([0, w]);
  const yS = d3.scaleLinear().domain([yExt[0]-pad, yExt[1]+pad]).range([h, 0]);

  // Gridlines
  svg.selectAll(".gx").data(xS.ticks(6)).enter().append("line")
    .attr("class","gridline").attr("x1",d=>xS(d)).attr("x2",d=>xS(d)).attr("y1",0).attr("y2",h);
  svg.selectAll(".gy").data(yS.ticks(6)).enter().append("line")
    .attr("class","gridline").attr("x1",0).attr("x2",w).attr("y1",d=>yS(d)).attr("y2",d=>yS(d));

  // Zero lines
  svg.append("line").attr("x1",xS(0)).attr("x2",xS(0)).attr("y1",0).attr("y2",h)
    .attr("stroke","#999").attr("stroke-dasharray","3,2").attr("stroke-width",0.8);
  svg.append("line").attr("x1",0).attr("x2",w).attr("y1",yS(0)).attr("y2",yS(0))
    .attr("stroke","#999").attr("stroke-dasharray","3,2").attr("stroke-width",0.8);

  // Axes
  svg.append("g").attr("class","axis").attr("transform",`translate(0,${{h}})`).call(d3.axisBottom(xS).ticks(6).tickSize(3));
  svg.append("g").attr("class","axis").call(d3.axisLeft(yS).ticks(6).tickSize(3));

  // Axis labels
  svg.append("text").attr("text-anchor","middle").attr("x",w/2).attr("y",h+45)
    .attr("font-size","10px").text("RRRM-1 Mouse monocyte/macrophage NES (FLT vs GC)");
  svg.append("text").attr("text-anchor","middle").attr("transform","rotate(-90)")
    .attr("x",-h/2).attr("y",-55).attr("font-size","10px")
    .text("I4 Human CD14+ Monocyte NES (flight vs pre)");

  // Dots
  const tooltip = d3.select("#tooltip");
  svg.selectAll(".dot").data(scatterData).enter().append("circle")
    .attr("class","dot").attr("cx",d=>xS(d.rrrm1_nes)).attr("cy",d=>yS(d.i4_nes)).attr("r",5)
    .attr("fill",d => d.concordant ? "{CONCORDANT_COLOR}" : "{DISCORDANT_COLOR}")
    .on("mouseover", function(event, d) {{
      tooltip.style("display","block")
        .html(`<b>${{d.pathway.replace("HALLMARK_","").replace(/_/g," ")}}</b><br>
               RRRM-1: ${{d.rrrm1_nes.toFixed(3)}}<br>I4: ${{d.i4_nes.toFixed(3)}}`);
    }})
    .on("mousemove", function(event) {{
      tooltip.style("left",(event.pageX+10)+"px").style("top",(event.pageY-20)+"px");
    }})
    .on("mouseout", () => tooltip.style("display","none"));

  // Title + stats
  svg.append("text").attr("text-anchor","middle").attr("x",w/2).attr("y",-30)
    .attr("font-size","11px").attr("font-weight","bold")
    .text("CD14+ Monocyte ↔ monocyte/macrophage");

  const rStr = monoR !== null ? monoR.toFixed(3) : "N/A";
  svg.append("text").attr("text-anchor","middle").attr("x",w/2).attr("y",-14)
    .attr("font-size","9px").attr("fill","#555")
    .text(`r = ${{rStr}} (95% CI: {mono_ci}), {p_str}, n = {mono_n} pathways`);

  // Legend
  const lg = svg.append("g").attr("transform",`translate(${{w-130}},${{h-50}})`);
  lg.append("circle").attr("cx",6).attr("cy",0).attr("r",5).attr("fill","{CONCORDANT_COLOR}");
  lg.append("text").attr("x",14).attr("y",4).attr("font-size","9px").text("Concordant");
  lg.append("circle").attr("cx",6).attr("cy",16).attr("r",5).attr("fill","{DISCORDANT_COLOR}");
  lg.append("text").attr("x",14).attr("y",20).attr("font-size","9px").text("Discordant");
}})();

// ─── Panel 2: Matched Pairs Bar Chart ────────────────────────────────
(function() {{
  const margin = {{top: 50, right: 20, bottom: 80, left: 55}};
  const W = 350, H = 400;
  const w = W - margin.left - margin.right;
  const h = H - margin.top - margin.bottom;

  const svg = d3.select("#barchart").append("svg")
    .attr("width", W).attr("height", H)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const x = d3.scaleBand().domain(barData.map(d=>d.label)).range([0,w]).padding(0.3);
  const y = d3.scaleLinear().domain([-1, 1]).range([h, 0]);

  // Axes
  svg.append("g").attr("class","axis").attr("transform",`translate(0,${{h}})`)
    .call(d3.axisBottom(x).tickSize(2))
    .selectAll("text").attr("transform","rotate(-30)").style("text-anchor","end").style("font-size","8px");
  svg.append("g").attr("class","axis").call(d3.axisLeft(y).ticks(8).tickSize(3));

  // Zero line
  svg.append("line").attr("x1",0).attr("x2",w)
    .attr("y1",y(0)).attr("y2",y(0)).attr("stroke","#999").attr("stroke-width",0.5);

  // E1 baseline reference
  svg.append("line").attr("class","ref-line")
    .attr("x1",0).attr("x2",w).attr("y1",y(e1_r)).attr("y2",y(e1_r));
  svg.append("text").attr("x",w-2).attr("y",y(e1_r)-4)
    .attr("text-anchor","end").attr("font-size","8px").attr("fill","#D55E00")
    .text(`E1 bulk r = ${{e1_r.toFixed(3)}}`);

  // Bars
  svg.selectAll(".bar").data(barData).enter().append("rect")
    .attr("x", d => x(d.label))
    .attr("y", d => d.r >= 0 ? y(d.r) : y(0))
    .attr("width", x.bandwidth())
    .attr("height", d => Math.abs(y(0) - y(d.r)))
    .attr("fill", d => d.type === "matched" ? "{CONCORDANT_COLOR}" : "#88CCEE")
    .attr("opacity", 0.8);

  // CI error bars
  barData.forEach(d => {{
    if (d.ci_low != null && d.ci_high != null) {{
      const cx = x(d.label) + x.bandwidth() / 2;
      svg.append("line").attr("x1",cx).attr("x2",cx)
        .attr("y1",y(d.ci_low)).attr("y2",y(d.ci_high))
        .attr("stroke","#333").attr("stroke-width",1);
      svg.append("line").attr("x1",cx-3).attr("x2",cx+3)
        .attr("y1",y(d.ci_low)).attr("y2",y(d.ci_low)).attr("stroke","#333");
      svg.append("line").attr("x1",cx-3).attr("x2",cx+3)
        .attr("y1",y(d.ci_high)).attr("y2",y(d.ci_high)).attr("stroke","#333");
    }}
  }});

  // n labels above bars
  barData.forEach(d => {{
    const cx = x(d.label) + x.bandwidth() / 2;
    const yPos = d.r >= 0 ? y(d.ci_high || d.r) - 10 : y(d.ci_low || d.r) + 14;
    svg.append("text").attr("x",cx).attr("y",yPos)
      .attr("text-anchor","middle").attr("font-size","7.5px").attr("fill","#555")
      .text(`n=${{d.n}}`);
  }});

  svg.append("text").attr("text-anchor","middle").attr("x",w/2).attr("y",-30)
    .attr("font-size","11px").attr("font-weight","bold")
    .text("Cell-Type Matched Concordance");
  svg.append("text").attr("text-anchor","middle").attr("transform","rotate(-90)")
    .attr("x",-h/2).attr("y",-40).attr("font-size","10px").text("Spearman r");

  // Legend
  const lg = svg.append("g").attr("transform",`translate(10,${{h-50}})`);
  lg.append("rect").attr("width",10).attr("height",10).attr("fill","{CONCORDANT_COLOR}").attr("opacity",0.8);
  lg.append("text").attr("x",14).attr("y",9).attr("font-size","8px").text("Matched");
  lg.append("rect").attr("y",14).attr("width",10).attr("height",10).attr("fill","#88CCEE").attr("opacity",0.8);
  lg.append("text").attr("x",14).attr("y",23).attr("font-size","8px").text("Sensitivity");
}})();

// ─── Panel 3: All-Pairs Heatmap ──────────────────────────────────────
(function() {{
  const margin = {{top: 80, right: 30, bottom: 20, left: 180}};
  const cellW = 60, cellH = 22;
  const nCols = matrix.rrrm1_celltypes.length;
  const nRows = matrix.i4_celltypes.length;
  const W = margin.left + cellW * nCols + margin.right;
  const H = margin.top + cellH * nRows + margin.bottom;

  const svg = d3.select("#heatmap").append("svg")
    .attr("width", W).attr("height", H)
    .append("g").attr("transform", `translate(${{margin.left}},${{margin.top}})`);

  const colorScale = d3.scaleDiverging(d3.interpolateRdBu).domain([1, 0, -1]);

  // Matched pair set for highlighting
  const matchedSet = new Set([
    "CD14+ Monocyte|monocyte_macrophage",
    "B Cell|b_cell",
    "CD4+ T Cell|t_cell",
  ]);

  // Cells
  matrix.i4_celltypes.forEach((i4ct, i) => {{
    matrix.rrrm1_celltypes.forEach((rct, j) => {{
      const r = matrix.spearman_r[i][j];
      const n = matrix.n_pathways[i][j];
      const isMatched = matchedSet.has(`${{i4ct}}|${{rct}}`);

      const cell = svg.append("rect")
        .attr("x", j * cellW).attr("y", i * cellH)
        .attr("width", cellW).attr("height", cellH)
        .attr("class", "heatmap-cell" + (isMatched ? " matched" : ""));

      if (r !== null && n >= {MIN_PATHWAYS}) {{
        cell.attr("fill", colorScale(r));
        svg.append("text")
          .attr("x", j * cellW + cellW/2).attr("y", i * cellH + cellH/2 + 3)
          .attr("text-anchor","middle").attr("font-size","8px")
          .attr("fill", Math.abs(r) > 0.5 ? "#fff" : "#333")
          .text(r.toFixed(2));
      }} else {{
        cell.attr("fill", "#f0f0f0");
        svg.append("text")
          .attr("x", j * cellW + cellW/2).attr("y", i * cellH + cellH/2 + 3)
          .attr("text-anchor","middle").attr("font-size","8px").attr("fill","#999")
          .text("×");
      }}
    }});
  }});

  // Row labels (I4 cell types)
  matrix.i4_celltypes.forEach((ct, i) => {{
    svg.append("text").attr("x",-6).attr("y", i * cellH + cellH/2 + 3)
      .attr("text-anchor","end").attr("font-size","9px").text(ct);
  }});

  // Column labels (RRRM-1 cell types)
  matrix.rrrm1_celltypes.forEach((ct, j) => {{
    svg.append("text")
      .attr("x", j * cellW + cellW/2).attr("y",-6)
      .attr("text-anchor","end").attr("font-size","9px")
      .attr("transform", `rotate(-40, ${{j * cellW + cellW/2}}, -6)`)
      .text(ct.replace(/_/g, " "));
  }});

  // Title
  svg.append("text").attr("x", (nCols * cellW) / 2).attr("y", -55)
    .attr("text-anchor","middle").attr("font-size","11px").attr("font-weight","bold")
    .text("All-Pairs Spearman r (I4 Human × RRRM-1 Mouse Blood)");

  // Color legend
  const lgX = nCols * cellW + 10;
  const lgScale = d3.scaleLinear().domain([-1, 1]).range([80, 0]);
  const lgAxis = d3.axisRight(lgScale).ticks(5).tickSize(3);
  const lgG = svg.append("g").attr("transform", `translate(${{lgX}}, 20)`);

  const gradient = svg.append("defs").append("linearGradient")
    .attr("id","hm-grad").attr("x1","0%").attr("y1","100%").attr("x2","0%").attr("y2","0%");
  gradient.append("stop").attr("offset","0%").attr("stop-color", colorScale(-1));
  gradient.append("stop").attr("offset","50%").attr("stop-color", colorScale(0));
  gradient.append("stop").attr("offset","100%").attr("stop-color", colorScale(1));

  lgG.append("rect").attr("width",12).attr("height",80).attr("fill","url(#hm-grad)");
  lgG.append("g").attr("transform","translate(12,0)").call(lgAxis).selectAll("text").attr("font-size","7px");
}})();

// ─── Download ────────────────────────────────────────────────────────
document.getElementById("download-btn").addEventListener("click", function() {{
  const svgs = document.querySelectorAll("svg");
  const serializer = new XMLSerializer();
  let combined = "";
  svgs.forEach((s, i) => {{
    combined += `<!-- Panel ${{i+1}} -->\\n` + serializer.serializeToString(s) + "\\n";
  }});
  const link = document.createElement("a");
  link.download = "F2D_crossspecies.svg";
  link.href = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(combined);
  link.click();
}});
</script>
</body>
</html>"""

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html)
    print(f"  HTML saved: {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("=== F2-D: Cross-Species Cell-Type Pathway Concordance ===\n")

    # 1. Load data
    print("[1] Loading RRRM-1 blood NES...")
    rrrm1_nes = load_rrrm1_blood_nes(RRRM1_BLOOD_DIR)

    print("\n[2] Loading I4 PBMC NES...")
    i4_nes = load_i4_pbmc_nes(I4_FGSEA_CSV)

    print("\n[3] Loading E1 baseline...")
    e1_baseline = load_e1_baseline(E1_JSON)

    # 2. Matched pairs
    print("\n[4] Computing matched pair concordance...")
    matched_results = []
    for i4_ct, rrrm1_ct, lineage in MATCHED_PAIRS:
        if i4_ct not in i4_nes:
            print(f"  WARNING: {i4_ct} not in I4 data, skipping")
            continue
        if rrrm1_ct not in rrrm1_nes:
            print(f"  WARNING: {rrrm1_ct} not in RRRM-1 data, skipping")
            continue

        print(f"\n  {i4_ct} ↔ {rrrm1_ct} ({lineage}):")
        result = compute_concordance(rrrm1_nes[rrrm1_ct], i4_nes[i4_ct],
                                     rrrm1_ct, i4_ct)
        result["i4_celltype"] = i4_ct
        result["rrrm1_celltype"] = rrrm1_ct
        result["matched_lineage"] = lineage

        if result.get("skipped"):
            print(f"    SKIPPED: {result['reason']}")
        else:
            r = result["spearman_r"]
            ci = f"[{result['ci_low']:.3f}, {result['ci_high']:.3f}]"
            exceeds = r > e1_baseline["spearman_r"]
            result["exceeds_e1"] = exceeds
            tag = "EXCEEDS E1" if exceeds else "below E1"
            print(f"    r = {r:.3f} {ci}, p_perm = {result['p_permutation']:.4f}, "
                  f"n = {result['n_pathways']}, sign_agree = {result['sign_agreement']:.0%} "
                  f"→ {tag}")

        matched_results.append(result)

    # 3. Sensitivity pairs
    print("\n[5] Computing sensitivity pairs...")
    sensitivity_results = []
    for i4_ct, rrrm1_ct, lineage in SENSITIVITY_PAIRS:
        if i4_ct not in i4_nes or rrrm1_ct not in rrrm1_nes:
            continue
        print(f"\n  {i4_ct} ↔ {rrrm1_ct} ({lineage}):")
        result = compute_concordance(rrrm1_nes[rrrm1_ct], i4_nes[i4_ct],
                                     rrrm1_ct, i4_ct)
        result["i4_celltype"] = i4_ct
        result["rrrm1_celltype"] = rrrm1_ct
        result["matched_lineage"] = lineage

        if result.get("skipped"):
            print(f"    SKIPPED: {result['reason']}")
        else:
            r = result["spearman_r"]
            print(f"    r = {r:.3f}, p_perm = {result['p_permutation']:.4f}, n = {result['n_pathways']}")

        sensitivity_results.append(result)

    # 4. All-pairs matrix
    print("\n[6] Computing all-pairs matrix...")
    matrix = compute_all_pairs_matrix(rrrm1_nes, i4_nes)
    print(f"  {len(matrix['i4_celltypes'])} × {len(matrix['rrrm1_celltypes'])} pairs")

    # 5. Summary
    print("\n=== Summary ===")
    print(f"E1 bulk baseline: r = {e1_baseline['spearman_r']:.3f} (n = {e1_baseline['n_pathways']})")
    for m in matched_results:
        if not m.get("skipped"):
            tag = ">" if m.get("exceeds_e1") else "<"
            print(f"  {m['matched_lineage']}: r = {m['spearman_r']:.3f} "
                  f"(n = {m['n_pathways']}) {tag} E1")

    # 6. Save JSON
    output = {
        "task": "F2-D",
        "description": "RRRM-1 mouse blood cell-type NES vs I4 human PBMC cell-type NES",
        "e1_baseline": {
            "spearman_r": e1_baseline["spearman_r"],
            "n_pathways": e1_baseline["n_pathways"],
            "description": "E1: mouse liver bulk vs JAXA cfRNA (different tissues)",
        },
        "matched_pairs": matched_results,
        "sensitivity": sensitivity_results,
        "all_pairs_matrix": matrix,
    }

    out_json = OUT_DIR / "F2D_crossspecies.json"
    with open(out_json, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\n  JSON saved: {out_json}")

    # 7. HTML figure
    print("\n[7] Generating HTML figure...")
    make_f2d_html(matched_results, sensitivity_results, matrix, e1_baseline,
                  FIG_DIR / "F2D_crossspecies_scatter.html")

    print("\n=== F2-D Complete ===")


if __name__ == "__main__":
    main()

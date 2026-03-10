#!/usr/bin/env python3
"""
cross_species_nes_comparison.py — GeneLab_benchmark v2, Category E1

Cross-species NES conservation analysis:
  Mouse liver fGSEA Hallmark NES (6 missions)
  vs.
  Human cfRNA fGSEA Hallmark NES (JAXA CFE OSD-530, 6 astronauts)

Approach:
  1. Load pre-computed human cfRNA fGSEA NES from R script output
     (run_human_cfrna_fgsea.R must have been run first to generate
      v2/processed/E1_human/human_cfrna_hallmark_fgsea.csv)
  2. Load mouse liver fGSEA Hallmark NES per mission, compute mission average
  3. Spearman r between human NES and mouse NES vectors (50 pathways)
  4. Bootstrap 95% CI, permutation p-value
  5. Save JSON results + HTML scatter plot

No ortholog mapping required (Hallmark gene set names are species-independent).

Sign convention:
  Both human and mouse: positive NES = pathway enriched in upregulated genes
  (flight vs. ground in mouse; in-flight vs. pre-flight in human cfRNA)

Usage:
  python3 v2/scripts/cross_species_nes_comparison.py
  python3 v2/scripts/cross_species_nes_comparison.py --out-dir v2/evaluation

Output:
  v2/evaluation/E1_crossspecies_nes.json
  v2/figures/E1_crossspecies_scatter.html
"""

import argparse
import json
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────

REPO_ROOT = Path(__file__).resolve().parent.parent.parent  # GeneLab_benchmark/

HUMAN_FGSEA_CSV = REPO_ROOT / "v2/processed/E1_human/human_cfrna_hallmark_fgsea.csv"
MOUSE_FGSEA_DIR = REPO_ROOT / "processed/fgsea/liver"
OUT_DIR = REPO_ROOT / "v2/evaluation"
FIG_DIR = REPO_ROOT / "v2/figures"

MISSIONS = ["MHU-2", "RR-1", "RR-3", "RR-6", "RR-8", "RR-9"]


# ── Step 1: Load human cfRNA ranking ──────────────────────────────────────────

def load_human_cfrna_ranking(path: Path) -> pd.Series:
    """Load human cfRNA gene ranking for fGSEA (in-flight vs. pre-flight)."""
    df = pd.read_csv(path)
    # Use edge_pre_vs_flight_diff: mean_flight_normalized - mean_pre_normalized
    # Positive = up in flight, negative = down in flight
    rank_col = "edge_pre_vs_flight_diff"
    if rank_col not in df.columns:
        raise ValueError(f"Column '{rank_col}' not found. Available: {list(df.columns)}")

    rnk = df.set_index("gene")[rank_col].dropna()
    # Remove genes with 0 variance (no ranking signal)
    rnk = rnk[rnk != 0]
    # Sort descending for prerank
    rnk = rnk.sort_values(ascending=False)
    print(f"  Human cfRNA ranking: {len(rnk):,} genes loaded")
    return rnk


# ── Step 2: Run fGSEA on human cfRNA ──────────────────────────────────────────

def run_human_fgsea(rnk: pd.Series, out_dir: Path) -> pd.DataFrame:
    """Run gseapy.prerank against MSigDB Hallmark for Homo sapiens."""
    try:
        import gseapy
    except ImportError:
        print("ERROR: gseapy not installed. Run: pip install gseapy")
        sys.exit(1)

    out_dir.mkdir(parents=True, exist_ok=True)
    prerank_out = out_dir / "human_cfrna_hallmark_prerank"

    print("  Running gseapy.prerank (Hallmark, hsapiens)...")
    print("  (Fetching MSigDB Hallmark gene sets from internet if not cached...)")

    res = gseapy.prerank(
        rnk=rnk,
        gene_sets="MSigDB_Hallmark_2020",  # fetched from gseapy's built-in MSigDB
        outdir=str(prerank_out),
        processes=4,
        permutation_num=1000,
        ascending=False,
        seed=42,
        verbose=False,
        min_size=15,
        max_size=500,
    )
    df = res.res2d
    print(f"  Human fGSEA complete: {len(df)} pathways")
    return df


def run_human_fgsea_alternative(rnk: pd.Series) -> pd.DataFrame:
    """
    Alternative: use gseapy with H: Hallmark MSigDB gene sets fetched as dict.
    Falls back to this if 'MSigDB_Hallmark_2020' key fails.
    """
    import gseapy

    print("  Fetching Hallmark gene sets via gseapy.get_library (alternative)...")
    hallmark = gseapy.get_library("MSigDB_Hallmark_2020", organism="Human")

    res = gseapy.prerank(
        rnk=rnk,
        gene_sets=hallmark,
        processes=4,
        permutation_num=1000,
        ascending=False,
        seed=42,
        verbose=False,
        min_size=15,
        max_size=500,
    )
    return res.res2d


# ── Step 2b: Load pre-computed human fGSEA (from R script) ────────────────────

def load_human_fgsea_from_r(csv_path: Path) -> pd.DataFrame:
    """Load pre-computed human cfRNA fGSEA results from run_human_cfrna_fgsea.R."""
    if not csv_path.exists():
        raise FileNotFoundError(
            f"Human fGSEA CSV not found: {csv_path}\n"
            "Run: Rscript v2/scripts/run_human_cfrna_fgsea.R"
        )
    df = pd.read_csv(csv_path)
    print(f"  Human fGSEA loaded: {len(df)} pathways from {csv_path.name}")
    return df


# ── Step 3: Load mouse liver NES ──────────────────────────────────────────────

def load_mouse_liver_nes(fgsea_dir: Path, missions: list) -> pd.DataFrame:
    """Load per-mission Hallmark NES values for mouse liver, return DataFrame."""
    records = []
    for mission in missions:
        fpath = fgsea_dir / f"{mission}_fgsea_hallmark.csv"
        if not fpath.exists():
            print(f"  WARNING: {fpath.name} not found, skipping")
            continue
        df = pd.read_csv(fpath)
        for _, row in df.iterrows():
            records.append({
                "pathway": row["pathway"],
                "NES": row["NES"],
                "padj": row["padj"],
                "mission": mission,
            })

    nes_df = pd.DataFrame(records)
    print(f"  Mouse liver NES: {len(missions)} missions, {nes_df['pathway'].nunique()} pathways")
    return nes_df


def compute_mouse_mission_average(nes_df: pd.DataFrame) -> pd.Series:
    """Compute mission-averaged NES per pathway."""
    avg = nes_df.groupby("pathway")["NES"].mean()
    return avg


# ── Step 4: Normalize pathway names ───────────────────────────────────────────

def normalize_pathway_name(name: str) -> str:
    """Normalize Hallmark pathway name for matching between mouse/human.

    Handles multiple naming conventions:
    - Mouse R fgsea: "HALLMARK_ADIPOGENESIS"
    - gseapy MSigDB_Hallmark_2020: "HALLMARK ADIPOGENESIS" or "ADIPOGENESIS"
    - gseapy enrichr format: "Adipogenesis" (title case with spaces)
    """
    name = name.upper().strip()
    # Replace spaces with underscores
    name = name.replace(" ", "_")
    # Remove any leading/trailing underscores
    name = name.strip("_")
    # Add HALLMARK_ prefix if not present
    if not name.startswith("HALLMARK_"):
        name = "HALLMARK_" + name
    return name


def extract_human_nes(human_df: pd.DataFrame) -> dict:
    """Extract NES values from gseapy result DataFrame, handling column name variations."""
    nes_map = {}

    # Detect column names (gseapy versions differ)
    term_col = None
    nes_col = None
    for c in human_df.columns:
        c_lower = c.lower()
        if c_lower in ("term", "name", "pathway"):
            term_col = c
        elif c_lower == "nes":
            nes_col = c

    if term_col is None or nes_col is None:
        print(f"  WARNING: Could not find Term/NES columns. Available: {list(human_df.columns)}")
        # Try positional approach for common gseapy output format
        # Standard columns: Term, ES, NES, NOM p-val, FDR q-val, ...
        if len(human_df.columns) >= 3:
            term_col = human_df.columns[0]
            # Find NES by position (usually column 2 or 3)
            for i, c in enumerate(human_df.columns):
                if "nes" in c.lower():
                    nes_col = c
                    break

    if term_col is None or nes_col is None:
        raise ValueError(f"Cannot find Term/NES columns in: {list(human_df.columns)}")

    for _, row in human_df.iterrows():
        pname = normalize_pathway_name(str(row[term_col]))
        try:
            nes_val = float(row[nes_col])
            if not np.isnan(nes_val):
                nes_map[pname] = nes_val
        except (ValueError, TypeError):
            pass

    return nes_map


def align_nes_vectors(human_df: pd.DataFrame, mouse_avg: pd.Series) -> tuple:
    """Align human and mouse NES vectors by pathway name."""
    human_nes = extract_human_nes(human_df)
    mouse_nes = {normalize_pathway_name(k): v for k, v in mouse_avg.items()}

    # Intersect
    common = sorted(set(human_nes.keys()) & set(mouse_nes.keys()))
    print(f"  Common pathways: {len(common)}")

    if len(common) == 0:
        print("  ERROR: No common pathways found!")
        print("  Human pathway names (sample):", list(human_nes.keys())[:5])
        print("  Mouse pathway names (sample):", list(mouse_nes.keys())[:5])
        raise ValueError("No common pathways between human and mouse NES")

    h_vec = np.array([human_nes[p] for p in common])
    m_vec = np.array([mouse_nes[p] for p in common])

    return h_vec, m_vec, common


# ── Step 5: Statistics ─────────────────────────────────────────────────────────

def spearman_with_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 1000, seed: int = 42):
    """Spearman r with bootstrap 95% CI and permutation p-value."""
    r, p_scipy = stats.spearmanr(x, y)

    rng = np.random.default_rng(seed)

    # Bootstrap CI
    boot_r = []
    n = len(x)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        try:
            r_b, _ = stats.spearmanr(x[idx], y[idx])
            boot_r.append(r_b)
        except Exception:
            pass
    boot_r = np.array(boot_r)
    ci_low = float(np.percentile(boot_r, 2.5))
    ci_high = float(np.percentile(boot_r, 97.5))

    # Permutation p (two-tailed)
    n_perm = 10000
    perm_r = []
    for _ in range(n_perm):
        y_perm = rng.permutation(y)
        r_p, _ = stats.spearmanr(x, y_perm)
        perm_r.append(r_p)
    perm_r = np.array(perm_r)
    p_perm = float(np.mean(np.abs(perm_r) >= np.abs(r)))

    # Sign agreement
    sign_agree = float(np.mean(np.sign(x) == np.sign(y)))

    return {
        "spearman_r": float(r),
        "ci_low": ci_low,
        "ci_high": ci_high,
        "p_scipy": float(p_scipy),
        "p_permutation": max(p_perm, 1 / n_perm),  # avoid p=0
        "sign_agreement_fraction": sign_agree,
        "n_pathways": len(x),
    }


def per_mission_correlations(human_df: pd.DataFrame, nes_df: pd.DataFrame) -> dict:
    """Compute Spearman r for each mouse mission vs. human NES."""
    results = {}
    missions = nes_df["mission"].unique()

    human_nes_map = extract_human_nes(human_df)

    for mission in sorted(missions):
        m_df = nes_df[nes_df["mission"] == mission]
        m_map = {normalize_pathway_name(r["pathway"]): r["NES"] for _, r in m_df.iterrows()}

        common = sorted(set(human_nes_map.keys()) & set(m_map.keys()))
        if len(common) < 10:
            continue

        h_vec = np.array([human_nes_map[p] for p in common])
        m_vec = np.array([m_map[p] for p in common])
        r, p = stats.spearmanr(h_vec, m_vec)
        results[mission] = {
            "spearman_r": float(r),
            "p_scipy": float(p),
            "n_pathways": len(common),
        }
    return results


# ── Step 6: HTML scatter plot ──────────────────────────────────────────────────

def make_scatter_html(
    h_vec: np.ndarray,
    m_vec: np.ndarray,
    pathways: list,
    stats_dict: dict,
    out_path: Path,
):
    """Generate self-contained HTML scatter plot (D3.js v7, Okabe-Ito)."""
    points = [
        {"pathway": p, "human_nes": float(h), "mouse_nes": float(m)}
        for p, h, m in zip(pathways, h_vec, m_vec)
    ]
    points_json = json.dumps(points)

    r = stats_dict["spearman_r"]
    ci_low = stats_dict["ci_low"]
    ci_high = stats_dict["ci_high"]
    p_perm = stats_dict["p_permutation"]
    sign_agree = stats_dict["sign_agreement_fraction"]
    n = stats_dict["n_pathways"]

    p_str = f"p = {p_perm:.3f}" if p_perm >= 0.001 else f"p < 0.001"
    subtitle = (
        f"Spearman r = {r:.3f} (95% CI: {ci_low:.3f}–{ci_high:.3f}), "
        f"{p_str}, sign agreement = {sign_agree:.0%}, n = {n} pathways"
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>E1 Cross-Species NES Conservation</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; font-size: 12px; background: #fff; margin: 20px; }}
  .axis text {{ font-size: 9px; }}
  .axis line, .axis path {{ stroke: #333; }}
  .gridline {{ stroke: #e0e0e0; stroke-width: 0.5; }}
  .dot {{ opacity: 0.8; }}
  .dot:hover {{ opacity: 1; stroke: #333; stroke-width: 1; }}
  .label-text {{ font-size: 7.5px; fill: #333; }}
  .tooltip {{
    position: absolute; background: rgba(255,255,255,0.95);
    border: 1px solid #ccc; padding: 6px; border-radius: 3px;
    font-size: 10px; pointer-events: none; display: none;
  }}
  .title {{ font-size: 12px; font-weight: bold; }}
  .subtitle {{ font-size: 9px; fill: #555; }}
  #download-btn {{
    margin-top: 10px; padding: 5px 12px; cursor: pointer;
    font-size: 11px; border: 1px solid #555; border-radius: 3px;
    background: #f5f5f5;
  }}
</style>
</head>
<body>
<div class="tooltip" id="tooltip"></div>
<div id="chart"></div>
<button id="download-btn">Download SVG</button>
<script>
const data = {points_json};

const margin = {{top: 50, right: 30, bottom: 60, left: 70}};
const width = 500 - margin.left - margin.right;
const height = 480 - margin.top - margin.bottom;

const svg = d3.select("#chart").append("svg")
  .attr("id", "main-svg")
  .attr("width", width + margin.left + margin.right)
  .attr("height", height + margin.top + margin.bottom)
  .append("g")
  .attr("transform", `translate(${{margin.left}},${{margin.top}})`);

// Axes
const xExtent = d3.extent(data, d => d.human_nes);
const yExtent = d3.extent(data, d => d.mouse_nes);
const pad = 0.3;
const xScale = d3.scaleLinear()
  .domain([xExtent[0] - pad, xExtent[1] + pad]).range([0, width]);
const yScale = d3.scaleLinear()
  .domain([yExtent[0] - pad, yExtent[1] + pad]).range([height, 0]);

// Gridlines
svg.append("g").attr("class", "gridlines")
  .selectAll("line.gridline-x")
  .data(xScale.ticks(6)).enter().append("line")
  .attr("class", "gridline")
  .attr("x1", d => xScale(d)).attr("x2", d => xScale(d))
  .attr("y1", 0).attr("y2", height);

svg.append("g").attr("class", "gridlines")
  .selectAll("line.gridline-y")
  .data(yScale.ticks(6)).enter().append("line")
  .attr("class", "gridline")
  .attr("x1", 0).attr("x2", width)
  .attr("y1", d => yScale(d)).attr("y2", d => yScale(d));

// Zero lines
svg.append("line")
  .attr("x1", xScale(0)).attr("x2", xScale(0))
  .attr("y1", 0).attr("y2", height)
  .attr("stroke", "#999").attr("stroke-width", 0.8).attr("stroke-dasharray", "3,2");
svg.append("line")
  .attr("x1", 0).attr("x2", width)
  .attr("y1", yScale(0)).attr("y2", yScale(0))
  .attr("stroke", "#999").attr("stroke-width", 0.8).attr("stroke-dasharray", "3,2");

// Axes
const xAxis = d3.axisBottom(xScale).ticks(6).tickSize(3);
const yAxis = d3.axisLeft(yScale).ticks(6).tickSize(3);

svg.append("g").attr("class", "axis").attr("transform", `translate(0,${{height}})`).call(xAxis);
svg.append("g").attr("class", "axis").call(yAxis);

// Axis labels
svg.append("text")
  .attr("text-anchor", "middle")
  .attr("x", width / 2).attr("y", height + 45)
  .attr("font-size", "10px")
  .text("Human cfRNA NES (JAXA CFE, in-flight vs. pre-flight)");

svg.append("text")
  .attr("text-anchor", "middle")
  .attr("transform", "rotate(-90)")
  .attr("x", -height / 2).attr("y", -55)
  .attr("font-size", "10px")
  .text("Mouse Liver NES (mission average, 6 missions)");

// Color: Okabe-Ito blue for concordant, orange for discordant
const concordColor = "#0072B2";
const discordColor = "#E69F00";

// Dots
const tooltip = d3.select("#tooltip");
svg.selectAll(".dot")
  .data(data).enter().append("circle")
  .attr("class", "dot")
  .attr("cx", d => xScale(d.human_nes))
  .attr("cy", d => yScale(d.mouse_nes))
  .attr("r", 5)
  .attr("fill", d => Math.sign(d.human_nes) === Math.sign(d.mouse_nes) ? concordColor : discordColor)
  .on("mouseover", function(event, d) {{
    tooltip.style("display", "block")
      .html(`<b>${{d.pathway.replace("HALLMARK_", "").replace(/_/g, " ")}}</b><br>
             Human NES: ${{d.human_nes.toFixed(3)}}<br>
             Mouse NES: ${{d.mouse_nes.toFixed(3)}}`);
  }})
  .on("mousemove", function(event) {{
    tooltip.style("left", (event.pageX + 10) + "px").style("top", (event.pageY - 20) + "px");
  }})
  .on("mouseout", () => tooltip.style("display", "none"));

// Title
svg.append("text").attr("class", "title")
  .attr("text-anchor", "middle").attr("x", width / 2).attr("y", -30)
  .text("E1: Cross-Species Hallmark NES Conservation");

svg.append("text").attr("class", "subtitle")
  .attr("text-anchor", "middle").attr("x", width / 2).attr("y", -14)
  .text("{subtitle}");

// Legend
const legend = svg.append("g").attr("transform", `translate(${{width - 130}}, ${{height - 60}})`);
legend.append("circle").attr("cx", 6).attr("cy", 0).attr("r", 5).attr("fill", concordColor);
legend.append("text").attr("x", 14).attr("y", 4).attr("font-size", "9px").text("Concordant sign");
legend.append("circle").attr("cx", 6).attr("cy", 18).attr("r", 5).attr("fill", discordColor);
legend.append("text").attr("x", 14).attr("y", 22).attr("font-size", "9px").text("Discordant sign");

// Download SVG
document.getElementById("download-btn").addEventListener("click", function() {{
  const svgEl = document.getElementById("main-svg");
  const serializer = new XMLSerializer();
  const source = '<?xml version="1.0"?>\\n' + serializer.serializeToString(svgEl);
  const link = document.createElement("a");
  link.download = "E1_crossspecies_scatter.svg";
  link.href = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(source);
  link.click();
}});
</script>
</body>
</html>"""

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html)
    print(f"  Scatter plot saved: {out_path}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="E1 Cross-Species NES Comparison")
    parser.add_argument("--out-dir", default=str(OUT_DIR), help="Output directory")
    parser.add_argument("--fig-dir", default=str(FIG_DIR), help="Figure directory")
    parser.add_argument("--n-boot", type=int, default=1000, help="Bootstrap iterations")
    parser.add_argument("--n-perm", type=int, default=10000, help="Permutation iterations")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    fig_dir = Path(args.fig_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    print("=== E1: Cross-Species NES Comparison ===")
    print(f"Human fGSEA CSV: {HUMAN_FGSEA_CSV}")
    print(f"Mouse fGSEA: {MOUSE_FGSEA_DIR}")

    # 1. Load pre-computed human cfRNA fGSEA results (from run_human_cfrna_fgsea.R)
    print("\n[1] Loading human cfRNA fGSEA NES (pre-computed by R script)...")
    human_df = load_human_fgsea_from_r(HUMAN_FGSEA_CSV)
    print(f"  Columns: {list(human_df.columns)}")
    print(human_df[["pathway", "NES", "padj"]].head(5).to_string())

    # 3. Load mouse liver NES
    print("\n[3] Loading mouse liver fGSEA NES...")
    nes_df = load_mouse_liver_nes(MOUSE_FGSEA_DIR, MISSIONS)
    mouse_avg = compute_mouse_mission_average(nes_df)

    # 4. Align and compute Spearman r
    print("\n[4] Aligning NES vectors...")
    h_vec, m_vec, pathways = align_nes_vectors(human_df, mouse_avg)

    print("\n[5] Computing Spearman correlation (mission average)...")
    stat_result = spearman_with_ci(h_vec, m_vec, n_boot=args.n_boot)
    r = stat_result["spearman_r"]
    ci_low = stat_result["ci_low"]
    ci_high = stat_result["ci_high"]
    p_perm = stat_result["p_permutation"]
    sign_agree = stat_result["sign_agreement_fraction"]

    print(f"  Spearman r = {r:.3f} ({ci_low:.3f}–{ci_high:.3f})")
    print(f"  Permutation p = {p_perm:.4f}")
    print(f"  Sign agreement = {sign_agree:.1%}")

    # 5. Per-mission correlation
    print("\n[6] Per-mission correlations...")
    per_mission = per_mission_correlations(human_df, nes_df)
    for m, v in sorted(per_mission.items()):
        print(f"  {m}: r = {v['spearman_r']:.3f} (p = {v['p_scipy']:.3f})")

    # 6. Top concordant and discordant pathways
    sign_match = [
        p for p, h, m in zip(pathways, h_vec, m_vec)
        if np.sign(h) == np.sign(m) and abs(h) > 0.5 and abs(m) > 0.5
    ]
    sign_discord = [
        p for p, h, m in zip(pathways, h_vec, m_vec)
        if np.sign(h) != np.sign(m) and abs(h) > 0.5 and abs(m) > 0.5
    ]
    print(f"\n  Concordant pathways (|NES|>0.5 both): {len(sign_match)}")
    for p in sign_match[:10]:
        idx = pathways.index(p)
        print(f"    {p}: human={h_vec[idx]:.2f}, mouse={m_vec[idx]:.2f}")
    print(f"\n  Discordant pathways (|NES|>0.5 both): {len(sign_discord)}")
    for p in sign_discord[:5]:
        idx = pathways.index(p)
        print(f"    {p}: human={h_vec[idx]:.2f}, mouse={m_vec[idx]:.2f}")

    # 7. Save results JSON
    pathway_data = [
        {
            "pathway": p,
            "human_nes": float(h),
            "mouse_nes_mean": float(m),
            "concordant": bool(np.sign(h) == np.sign(m)),
        }
        for p, h, m in zip(pathways, h_vec, m_vec)
    ]

    output = {
        "task": "E1",
        "description": "Mouse liver vs. human cfRNA Hallmark NES Spearman correlation",
        "human_data": {
            "source": "JAXA CFE cfRNA (OSD-530)",
            "ranking": "edge_pre_vs_flight_diff (mean_flight - mean_pre, normalized)",
            "n_pathways": len(human_df),
            "file": HUMAN_FGSEA_CSV.name,
        },
        "mouse_data": {
            "tissue": "liver",
            "missions": MISSIONS,
            "db": "hallmark",
            "averaging": "arithmetic mean across missions",
        },
        "results": {
            "mission_averaged": stat_result,
            "per_mission": per_mission,
        },
        "pathways": pathway_data,
        "concordant_pathways_high": sign_match,
        "discordant_pathways_high": sign_discord,
    }

    out_json = out_dir / "E1_crossspecies_nes.json"
    out_json.write_text(json.dumps(output, indent=2))
    print(f"\n  Results saved: {out_json}")

    # 8. HTML scatter plot
    print("\n[7] Generating scatter plot...")
    make_scatter_html(h_vec, m_vec, pathways, stat_result, fig_dir / "E1_crossspecies_scatter.html")

    print("\n=== E1 Complete ===")
    print(f"  Spearman r = {r:.3f} (95% CI: {ci_low:.3f}–{ci_high:.3f}), p = {p_perm:.4f}")
    print(f"  Sign agreement = {sign_agree:.1%}")
    print(f"  Results: {out_json}")
    print(f"  Figure: {fig_dir / 'E1_crossspecies_scatter.html'}")


if __name__ == "__main__":
    main()

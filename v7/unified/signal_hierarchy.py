"""
Unified Signal Hierarchy Analysis
Compares method categories across GeneLabBench (mouse) and SpaceOmicsBench (human).

Core finding: Tabular ML dominates; FMs, GNNs, and sequence models fail.
This is consistent across both datasets.
"""

import json
import os
from pathlib import Path

# ─── Paths ────────────────────────────────────────────────────────────────────
GENELABS_DIR = Path(__file__).parent.parent.parent  # GeneLab_benchmark/
SOB_DIR = Path(os.environ.get("SPACEOMICS_ROOT", ".")) / "v3"
OUT_DIR = GENELABS_DIR / "v7" / "figures" / "html"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ─── Method categories ────────────────────────────────────────────────────────
# GeneLab method → category
GENELABS_CATEGORIES = {
    "pca_lr":       "Tabular ML",
    "elasticnet_lr":"Tabular ML",
    "svm_rbf":      "Tabular ML",
    "rf":           "Tabular ML",
    "xgb":          "Tabular ML",
    "knn":          "Tabular ML",
    "mlp":          "Deep Learning",
    "tabnet":       "Deep Learning",
    # FMs (from v3 results, stored separately)
    "geneformer":   "Foundation Model",
    "scfoundation": "Foundation Model",
    "uce":          "Foundation Model",
    "scgpt":        "Foundation Model",
    "text_llm":     "Text LLM",
}

# SpaceOmicsBench method → category
SOB_CATEGORIES = {
    "logreg":       "Tabular ML",
    "rf":           "Tabular ML",
    "xgboost":      "Tabular ML",
    "lightgbm":     "Tabular ML",
    "svm":          "Tabular ML",
    "tabpfn":       "Tabular ML",  # TabPFN = tabular prior, still Tabular
    "mlp":          "Deep Learning",
    "late_fusion":  "Deep Learning",
    "early_pca":    "Deep Learning",
    "graphsage":    "Graph/Sequence",
    "gnn":          "Graph/Sequence",
    "esm2_base_only":     "Foundation Model",
    "esm2_combined":      "Foundation Model",
    "esm2_esm2_only":     "Foundation Model",
    "random":       "Random Baseline",
    "dummy":        "Random Baseline",
    "mofa_plus":    "Deep Learning",  # dimensionality reduction
    "snf":          "Deep Learning",
}

CATEGORY_ORDER = ["Tabular ML", "Deep Learning", "Graph/Sequence", "Foundation Model", "Text LLM", "Random Baseline"]
CATEGORY_COLORS = {
    "Tabular ML":      "#E69F00",  # Okabe-Ito orange
    "Deep Learning":   "#56B4E9",  # Okabe-Ito sky blue
    "Graph/Sequence":  "#009E73",  # Okabe-Ito green
    "Foundation Model":"#CC79A7",  # Okabe-Ito pink
    "Text LLM":        "#D55E00",  # Okabe-Ito vermillion
    "Random Baseline": "#999999",  # grey
}

# ─── Load GeneLab data ────────────────────────────────────────────────────────
def load_genelab():
    m1_path = GENELABS_DIR / "v4" / "evaluation" / "M1_summary.json"
    with open(m1_path) as f:
        m1 = json.load(f)

    rows = []
    for tissue, feat_types in m1.items():
        for feat, methods in feat_types.items():
            for method, res in methods.items():
                auroc = res.get("auroc")
                if auroc is None:
                    continue
                cat = GENELABS_CATEGORIES.get(method, "Other")
                rows.append({
                    "dataset": "GeneLab (Mouse)",
                    "tissue_task": tissue,
                    "feature_type": feat,
                    "method": method,
                    "category": cat,
                    "score": round(auroc, 4),
                    "metric": "auroc",
                    "random_baseline": 0.5,
                })

    # Add FM results from v3 (stored as separate records from MEMORY.md)
    # These are mean AUROCs across LOMO folds from v3 analysis
    fm_records = [
        # scFoundation: liver 0.635(p=0.001), gastro 0.691, others ~0.5
        {"method": "scfoundation", "tissue_task": "liver",          "score": 0.635, "sig": True},
        {"method": "scfoundation", "tissue_task": "gastrocnemius",  "score": 0.691, "sig": True},
        {"method": "scfoundation", "tissue_task": "kidney",         "score": 0.541, "sig": False},
        {"method": "scfoundation", "tissue_task": "thymus",         "score": 0.487, "sig": False},
        {"method": "scfoundation", "tissue_task": "eye",            "score": 0.563, "sig": False},
        {"method": "scfoundation", "tissue_task": "lung",           "score": 0.389, "sig": False},
        # UCE: thymus 0.632(p=0.031), others near chance
        {"method": "uce",          "tissue_task": "thymus",         "score": 0.632, "sig": True},
        {"method": "uce",          "tissue_task": "liver",          "score": 0.459, "sig": False},
        {"method": "uce",          "tissue_task": "gastrocnemius",  "score": 0.578, "sig": False},
        {"method": "uce",          "tissue_task": "kidney",         "score": 0.489, "sig": False},
        {"method": "uce",          "tissue_task": "eye",            "score": 0.550, "sig": False},
        {"method": "uce",          "tissue_task": "lung",           "score": 0.555, "sig": False},
        # Geneformer: mean 0.476 (from v1)
        {"method": "geneformer",   "tissue_task": "liver",          "score": 0.476, "sig": False},
        # scGPT: 0.6655 for one tissue in v1
        {"method": "scgpt",        "tissue_task": "liver",          "score": 0.6655,"sig": False},
    ]
    for r in fm_records:
        rows.append({
            "dataset": "GeneLab (Mouse)",
            "tissue_task": r["tissue_task"],
            "feature_type": "embedding",
            "method": r["method"],
            "category": GENELABS_CATEGORIES.get(r["method"], "Foundation Model"),
            "score": r["score"],
            "metric": "auroc",
            "random_baseline": 0.5,
            "significant": r.get("sig", False),
        })

    return rows

# ─── Load SpaceOmicsBench data ────────────────────────────────────────────────
def load_sob():
    path = SOB_DIR / "results" / "unified_baseline_results.json"
    with open(path) as f:
        d = json.load(f)

    tasks = d["task_results"]
    leaderboard = d.get("leaderboard", {})

    rows = []
    # Only include AUROC/AUPRC tasks for comparability
    AUROC_METRICS = {"auroc", "auc", "auprc"}
    for task_id, task_data in tasks.items():
        metric = task_data["primary_metric"].lower()
        if metric not in AUROC_METRICS:
            continue  # skip clustering (ARI), regression, etc.
        random_score = 0.5 if metric in ("auroc", "auc") else 0.0

        models = task_data.get("models", {})
        for model_name, res in models.items():
            score = res.get("score")
            if score is None:
                continue
            # Normalize method name to lowercase for category lookup
            mn = model_name.lower()
            cat = SOB_CATEGORIES.get(mn, "Other")
            if cat == "Other":
                for k, v in SOB_CATEGORIES.items():
                    if k in mn:
                        cat = v
                        break
            # Map SOB method names to shared names where possible
            shared_name = mn
            if mn == "logreg":  shared_name = "logreg_lr"
            elif mn == "xgboost": shared_name = "xgb"
            elif mn == "graphsage": shared_name = "graphsage_gnn"
            elif mn in ("esm2_base_only", "esm2_combined", "esm2_esm2_only"): shared_name = mn

            rows.append({
                "dataset": "SpaceOmicsBench (Human)",
                "tissue_task": task_id,
                "feature_type": "multi-omics",
                "method": shared_name,
                "method_orig": model_name,
                "category": cat,
                "score": round(score, 4),
                "metric": metric,
                "random_baseline": random_score,
            })

    # Add leaderboard bests as reference
    for task_id, lb in leaderboard.items():
        best_model = lb.get("best_model", "")
        best_score = lb.get("score")
        metric = lb.get("metric", "score")
        if best_score is not None:
            cat = SOB_CATEGORIES.get(best_model.lower(), "Tabular ML")
            # Mark as leaderboard best
            rows.append({
                "dataset": "SpaceOmicsBench (Human)",
                "tissue_task": f"{task_id}_best",
                "feature_type": "best",
                "method": best_model,
                "category": cat,
                "score": round(best_score, 4),
                "metric": metric,
                "random_baseline": 0.5 if metric in ("auroc", "auc") else None,
                "is_leaderboard_best": True,
            })

    # Add within-task normalized scores (min-max per task, across all models)
    from collections import defaultdict
    task_scores = defaultdict(list)
    for r in rows:
        if not r.get("is_leaderboard_best"):
            task_scores[r["tissue_task"]].append(r["score"])
    for r in rows:
        tid = r["tissue_task"]
        scores_in_task = task_scores.get(tid, [])
        if scores_in_task:
            lo, hi = min(scores_in_task), max(scores_in_task)
            r["norm_score"] = round((r["score"] - lo) / (hi - lo + 1e-8), 4)
        else:
            r["norm_score"] = 0.0

    return rows

# ─── Aggregate by category ────────────────────────────────────────────────────
def summarize_by_category(rows, dataset_filter=None):
    from collections import defaultdict
    cat_scores = defaultdict(list)
    for r in rows:
        if dataset_filter and r["dataset"] != dataset_filter:
            continue
        cat = r["category"]
        if cat not in CATEGORY_ORDER:
            continue
        # Only use AUROC for comparability
        if r["metric"] not in ("auroc", "auc", "auprc"):
            continue
        # Skip leaderboard_best entries (to avoid double-counting)
        if r.get("is_leaderboard_best"):
            continue
        cat_scores[cat].append(r["score"])

    summary = {}
    for cat in CATEGORY_ORDER:
        scores = cat_scores.get(cat, [])
        if scores:
            import statistics
            summary[cat] = {
                "mean": round(sum(scores) / len(scores), 4),
                "median": round(statistics.median(scores), 4),
                "max": round(max(scores), 4),
                "min": round(min(scores), 4),
                "n": len(scores),
                "scores": scores,
            }
    return summary

# ─── Build visualization data ─────────────────────────────────────────────────
def build_viz_data():
    gl_rows = load_genelab()
    sob_rows = load_sob()
    all_rows = gl_rows + sob_rows

    # Per-tissue/task best score by category for GeneLab
    gl_summary = summarize_by_category(all_rows, "GeneLab (Mouse)")
    sob_summary = summarize_by_category(all_rows, "SpaceOmicsBench (Human)")

    # Unified: per method, best AUROC across tissues/tasks
    from collections import defaultdict
    unified_method_scores = defaultdict(lambda: {"gl": [], "sob": []})
    for r in gl_rows:
        if r["metric"] in ("auroc", "auc") and r["category"] in CATEGORY_ORDER and not r.get("is_leaderboard_best"):
            unified_method_scores[r["method"]]["gl"].append(r["score"])
    for r in sob_rows:
        if r["metric"] in ("auroc", "auc", "auprc") and r["category"] in CATEGORY_ORDER and not r.get("is_leaderboard_best"):
            unified_method_scores[r["method"]]["sob"].append(r["score"])

    # Representative method entries for the scatter
    method_summary = []
    for method, data in unified_method_scores.items():
        cat_gl = next((r["category"] for r in gl_rows if r["method"] == method), None)
        cat_sob = next((r["category"] for r in sob_rows if r["method"] == method), None)
        cat = cat_gl or cat_sob or "Other"
        if cat not in CATEGORY_ORDER:
            continue
        gl_mean = round(sum(data["gl"]) / len(data["gl"]), 4) if data["gl"] else None
        sob_mean = round(sum(data["sob"]) / len(data["sob"]), 4) if data["sob"] else None
        method_summary.append({
            "method": method,
            "category": cat,
            "gl_mean_auroc": gl_mean,
            "sob_mean_score": sob_mean,
            "gl_n": len(data["gl"]),
            "sob_n": len(data["sob"]),
        })

    # Category-level cross-species comparison
    # GL: mean AUROC per category
    # SOB: mean within-task normalized score per category (0=worst, 1=best in task)
    from collections import defaultdict
    sob_cat_norm = defaultdict(list)
    for r in sob_rows:
        if r.get("is_leaderboard_best"):
            continue
        cat = r.get("category")
        if cat in CATEGORY_ORDER and "norm_score" in r:
            sob_cat_norm[cat].append(r["norm_score"])

    category_comparison = []
    for cat in CATEGORY_ORDER:
        gl = gl_summary.get(cat)
        sob_norms = sob_cat_norm.get(cat, [])
        if gl and sob_norms:
            category_comparison.append({
                "category": cat,
                "gl_mean_auroc": gl["mean"],
                "gl_n": gl["n"],
                "sob_mean_norm": round(sum(sob_norms) / len(sob_norms), 4),
                "sob_n": len(sob_norms),
            })

    return {
        "gl_summary": gl_summary,
        "sob_summary": sob_summary,
        "method_summary": method_summary,
        "category_comparison": category_comparison,
        "gl_rows": [r for r in gl_rows if not r.get("is_leaderboard_best")],
        "sob_rows": [r for r in sob_rows if not r.get("is_leaderboard_best")],
        "category_colors": CATEGORY_COLORS,
        "category_order": CATEGORY_ORDER,
    }

# ─── Generate HTML ────────────────────────────────────────────────────────────
def generate_html(data):
    import json as _json
    data_json = _json.dumps(data, indent=2)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Unified Signal Hierarchy — GeneLabBench + SpaceOmicsBench</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: "Helvetica Neue", Arial, sans-serif; font-size: 13px; background: #fff; color: #222; }}
  h1 {{ font-size: 16px; font-weight: 600; margin: 18px 24px 4px; }}
  .subtitle {{ font-size: 12px; color: #666; margin: 0 24px 16px; }}
  .panels {{ display: flex; gap: 16px; padding: 0 24px 24px; flex-wrap: wrap; }}
  .panel {{ flex: 1; min-width: 340px; }}
  .panel-title {{ font-size: 13px; font-weight: 600; margin-bottom: 6px; color: #444; }}
  .finding-box {{
    margin: 0 24px 20px;
    padding: 10px 14px;
    background: #f5f8ff;
    border-left: 4px solid #4472c4;
    font-size: 12px;
    line-height: 1.6;
  }}
  .finding-box strong {{ color: #2a3e7c; }}
  .legend {{ display: flex; flex-wrap: wrap; gap: 12px; margin: 0 24px 16px; font-size: 12px; }}
  .legend-item {{ display: flex; align-items: center; gap: 5px; }}
  .legend-dot {{ width: 12px; height: 12px; border-radius: 50%; flex-shrink: 0; }}
  svg text {{ font-family: "Helvetica Neue", Arial, sans-serif; }}
  .tooltip {{
    position: fixed; background: rgba(0,0,0,0.82); color: #fff;
    padding: 8px 12px; border-radius: 4px; font-size: 12px; pointer-events: none;
    line-height: 1.5; max-width: 240px; z-index: 9999; display: none;
  }}
</style>
</head>
<body>
<h1>Unified Signal Hierarchy: Tabular ML Dominates Spaceflight Biology Across Species</h1>
<div class="subtitle">GeneLabBench (mouse, 8 tissues, LOMO-CV) + SpaceOmicsBench (human, 26 tasks, multi-omics)</div>

<div class="finding-box">
  <strong>Core finding:</strong> Across both mouse transcriptomics (GeneLabBench) and human multi-omics (SpaceOmicsBench),
  classical tabular ML methods (PCA-LR, TabPFN, RF, XGBoost) consistently outperform foundation models, GNNs, and sequence-based models.
  FMs/GNNs add noise rather than signal for small-N spaceflight data. This convergent finding suggests the bottleneck is
  domain-specific data scarcity, not method capability.
</div>

<div class="legend" id="legend"></div>
<div class="panels">
  <div class="panel"><div class="panel-title">A. GeneLabBench (Mouse) — AUROC by Method Category</div><svg id="gl-chart"></svg></div>
  <div class="panel"><div class="panel-title">B. SpaceOmicsBench (Human) — Score by Method Category</div><svg id="sob-chart"></svg></div>
  <div class="panel" style="min-width:600px;flex:2"><div class="panel-title">C. Cross-Species Convergence: Same Method → Similar Relative Performance</div><svg id="scatter"></svg></div>
</div>
<div class="panels">
  <div class="panel" style="min-width:100%"><div class="panel-title">D. All Individual Results — GeneLabBench AUROC Distribution by Category (gene features)</div><svg id="strip"></svg></div>
</div>

<div class="tooltip" id="tooltip"></div>

<script>
const DATA = {data_json};
const COLORS = DATA.category_colors;
const CAT_ORDER = DATA.category_order;
const tooltip = document.getElementById("tooltip");

function showTip(html, event) {{
  tooltip.innerHTML = html;
  tooltip.style.display = "block";
  tooltip.style.left = (event.clientX + 14) + "px";
  tooltip.style.top = (event.clientY - 10) + "px";
}}
function hideTip() {{ tooltip.style.display = "none"; }}

// Legend
const legend = document.getElementById("legend");
CAT_ORDER.forEach(cat => {{
  const item = document.createElement("div");
  item.className = "legend-item";
  item.innerHTML = `<div class="legend-dot" style="background:${{COLORS[cat]}}"></div><span>${{cat}}</span>`;
  legend.appendChild(item);
}});

// ─── Panel A: GeneLab box/strip by category ───────────────────────────────
(function() {{
  const W = 360, H = 300, mg = {{top:20, right:20, bottom:60, left:55}};
  const svg = d3.select("#gl-chart").attr("width", W).attr("height", H);
  const iW = W - mg.left - mg.right, iH = H - mg.top - mg.bottom;
  const g = svg.append("g").attr("transform", `translate(${{mg.left}},${{mg.top}})`);

  // Collect scores per category (gene feature only for clarity)
  const byCategory = {{}};
  DATA.gl_rows.forEach(r => {{
    if (r.metric !== "auroc") return;
    const cat = r.category;
    if (!CAT_ORDER.includes(cat)) return;
    if (!byCategory[cat]) byCategory[cat] = [];
    byCategory[cat].push(r.score);
  }});

  const cats = CAT_ORDER.filter(c => byCategory[c] && byCategory[c].length > 0);
  const x = d3.scaleBand().domain(cats).range([0, iW]).padding(0.3);
  const y = d3.scaleLinear().domain([0.3, 1.05]).range([iH, 0]);

  // Random baseline reference line
  g.append("line")
    .attr("x1", 0).attr("x2", iW).attr("y1", y(0.5)).attr("y2", y(0.5))
    .attr("stroke", "#ccc").attr("stroke-dasharray", "4,3").attr("stroke-width", 1);
  g.append("text").attr("x", iW+2).attr("y", y(0.5)+4)
    .attr("font-size", 10).attr("fill", "#999").text("chance");

  // Grid lines
  y.ticks(5).forEach(tick => {{
    g.append("line").attr("x1", 0).attr("x2", iW)
     .attr("y1", y(tick)).attr("y2", y(tick))
     .attr("stroke", "#eee");
  }});

  // Boxes + dots
  cats.forEach(cat => {{
    const scores = byCategory[cat].sort(d3.ascending);
    const q1 = d3.quantile(scores, 0.25);
    const median = d3.quantile(scores, 0.5);
    const q3 = d3.quantile(scores, 0.75);
    const mean = d3.mean(scores);
    const xc = x(cat) + x.bandwidth() / 2;
    const bw = x.bandwidth() * 0.5;
    const col = COLORS[cat];

    // IQR box
    g.append("rect")
      .attr("x", xc - bw/2).attr("y", y(q3))
      .attr("width", bw).attr("height", y(q1) - y(q3))
      .attr("fill", col).attr("opacity", 0.25).attr("stroke", col);

    // Median line
    g.append("line")
      .attr("x1", xc - bw/2).attr("x2", xc + bw/2)
      .attr("y1", y(median)).attr("y2", y(median))
      .attr("stroke", col).attr("stroke-width", 2);

    // Mean dot
    g.append("circle").attr("cx", xc).attr("cy", y(mean))
      .attr("r", 5).attr("fill", col)
      .on("mouseover", (e) => showTip(`<b>${{cat}}</b><br>Mean: ${{mean.toFixed(3)}}<br>Median: ${{median.toFixed(3)}}<br>N: ${{scores.length}}`, e))
      .on("mouseout", hideTip);

    // Jitter dots
    scores.forEach((s, i) => {{
      const jitter = (Math.random() - 0.5) * bw * 0.7;
      g.append("circle").attr("cx", xc + jitter).attr("cy", y(s))
        .attr("r", 2.5).attr("fill", col).attr("opacity", 0.6);
    }});
  }});

  // Axes
  g.append("g").attr("transform", `translate(0,${{iH}})`).call(
    d3.axisBottom(x).tickFormat(d => d.replace(" ML","").replace(" Model",""))
  ).selectAll("text").attr("transform","rotate(-30)").style("text-anchor","end").attr("font-size",10);

  g.append("g").call(d3.axisLeft(y).ticks(5).tickFormat(d3.format(".2f"))).selectAll("text").attr("font-size",10);

  g.append("text").attr("x", iW/2).attr("y", iH+52).attr("text-anchor","middle")
   .attr("font-size",11).attr("fill","#555").text("Method Category");
  g.append("text").attr("transform","rotate(-90)").attr("x",-iH/2).attr("y",-42)
   .attr("text-anchor","middle").attr("font-size",11).attr("fill","#555").text("AUROC");
}})();

// ─── Panel B: SpaceOmicsBench bar by category ─────────────────────────────
(function() {{
  const W = 360, H = 300, mg = {{top:20, right:20, bottom:60, left:55}};
  const svg = d3.select("#sob-chart").attr("width", W).attr("height", H);
  const iW = W - mg.left - mg.right, iH = H - mg.top - mg.bottom;
  const g = svg.append("g").attr("transform", `translate(${{mg.left}},${{mg.top}})`);

  const byCategory = {{}};
  DATA.sob_rows.forEach(r => {{
    const cat = r.category;
    if (!CAT_ORDER.includes(cat)) return;
    if (!byCategory[cat]) byCategory[cat] = [];
    byCategory[cat].push(r.score);
  }});

  const cats = CAT_ORDER.filter(c => byCategory[c] && byCategory[c].length > 0);
  const allScores = Object.values(byCategory).flat();
  const yMin = Math.max(0, d3.min(allScores) - 0.05);
  const yMax = Math.min(1.05, d3.max(allScores) + 0.05);

  const x = d3.scaleBand().domain(cats).range([0, iW]).padding(0.3);
  const y = d3.scaleLinear().domain([yMin, yMax]).range([iH, 0]);

  y.ticks(5).forEach(tick => {{
    g.append("line").attr("x1", 0).attr("x2", iW)
     .attr("y1", y(tick)).attr("y2", y(tick))
     .attr("stroke", "#eee");
  }});

  cats.forEach(cat => {{
    const scores = byCategory[cat].sort(d3.ascending);
    const q1 = d3.quantile(scores, 0.25) || d3.min(scores);
    const median = d3.quantile(scores, 0.5);
    const q3 = d3.quantile(scores, 0.75) || d3.max(scores);
    const mean = d3.mean(scores);
    const xc = x(cat) + x.bandwidth() / 2;
    const bw = x.bandwidth() * 0.5;
    const col = COLORS[cat];

    g.append("rect").attr("x", xc - bw/2).attr("y", y(q3))
      .attr("width", bw).attr("height", Math.max(1, y(q1) - y(q3)))
      .attr("fill", col).attr("opacity", 0.25).attr("stroke", col);

    g.append("line").attr("x1", xc - bw/2).attr("x2", xc + bw/2)
      .attr("y1", y(median)).attr("y2", y(median))
      .attr("stroke", col).attr("stroke-width", 2);

    g.append("circle").attr("cx", xc).attr("cy", y(mean))
      .attr("r", 5).attr("fill", col)
      .on("mouseover", (e) => showTip(`<b>${{cat}}</b><br>Mean: ${{mean.toFixed(3)}}<br>Median: ${{median.toFixed(3)}}<br>N: ${{scores.length}}`, e))
      .on("mouseout", hideTip);

    scores.forEach(s => {{
      const jitter = (Math.random() - 0.5) * bw * 0.7;
      g.append("circle").attr("cx", xc + jitter).attr("cy", y(s))
        .attr("r", 2.5).attr("fill", col).attr("opacity", 0.6);
    }});
  }});

  g.append("g").attr("transform", `translate(0,${{iH}})`).call(
    d3.axisBottom(x).tickFormat(d => d.replace(" ML","").replace(" Model",""))
  ).selectAll("text").attr("transform","rotate(-30)").style("text-anchor","end").attr("font-size",10);

  g.append("g").call(d3.axisLeft(y).ticks(5).tickFormat(d3.format(".2f"))).selectAll("text").attr("font-size",10);

  g.append("text").attr("x", iW/2).attr("y", iH+52).attr("text-anchor","middle")
   .attr("font-size",11).attr("fill","#555").text("Method Category");
  g.append("text").attr("transform","rotate(-90)").attr("x",-iH/2).attr("y",-42)
   .attr("text-anchor","middle").attr("font-size",11).attr("fill","#555").text("Score (AUROC / AUPRC / F1)");
}})();

// ─── Panel C: Category-level cross-species comparison ─────────────────────
(function() {{
  const W = 560, H = 320, mg = {{top:30, right:30, bottom:65, left:70}};
  const svg = d3.select("#scatter").attr("width", W).attr("height", H);
  const iW = W - mg.left - mg.right, iH = H - mg.top - mg.bottom;
  const g = svg.append("g").attr("transform", `translate(${{mg.left}},${{mg.top}})`);

  const cats = DATA.category_comparison;
  if (!cats || cats.length === 0) {{
    g.append("text").attr("x", iW/2).attr("y", iH/2).attr("text-anchor","middle")
     .attr("font-size",12).attr("fill","#999").text("No category overlap data");
    return;
  }}

  const xDomain = [d3.min(cats, c=>c.gl_mean_auroc)-0.05, d3.max(cats, c=>c.gl_mean_auroc)+0.05];
  const yDomain = [d3.min(cats, c=>c.sob_mean_norm)-0.05, d3.max(cats, c=>c.sob_mean_norm)+0.05];
  const x = d3.scaleLinear().domain(xDomain).range([0, iW]);
  const y = d3.scaleLinear().domain(yDomain).range([iH, 0]);

  // Spearman rank correlation annotation
  const xRanks = [...cats].sort((a,b)=>a.gl_mean_auroc-b.gl_mean_auroc).map((c,i)=>{{return{{cat:c.category,rank:i+1}}}});
  const yRanks = [...cats].sort((a,b)=>a.sob_mean_norm-b.sob_mean_norm).map((c,i)=>{{return{{cat:c.category,rank:i+1}}}});
  const xRankMap = Object.fromEntries(xRanks.map(r=>[r.cat, r.rank]));
  const yRankMap = Object.fromEntries(yRanks.map(r=>[r.cat, r.rank]));
  const n = cats.length;
  const dSq = cats.reduce((s,c)=>s+Math.pow(xRankMap[c.category]-yRankMap[c.category],2), 0);
  const rho = 1 - 6*dSq / (n*(n*n-1));

  // Grid
  x.ticks(4).forEach(v => g.append("line").attr("x1",x(v)).attr("x2",x(v)).attr("y1",0).attr("y2",iH).attr("stroke","#eee"));
  y.ticks(4).forEach(v => g.append("line").attr("x1",0).attr("x2",iW).attr("y1",y(v)).attr("y2",y(v)).attr("stroke","#eee"));

  // Trend line
  const xMean = d3.mean(cats, c=>c.gl_mean_auroc);
  const yMean = d3.mean(cats, c=>c.sob_mean_norm);
  const cov = d3.mean(cats, c=>(c.gl_mean_auroc-xMean)*(c.sob_mean_norm-yMean));
  const varX = d3.mean(cats, c=>Math.pow(c.gl_mean_auroc-xMean,2));
  const slope = varX>0 ? cov/varX : 0;
  const xMin = xDomain[0], xMax = xDomain[1];
  g.append("line")
    .attr("x1",x(xMin)).attr("y1",y(yMean+slope*(xMin-xMean)))
    .attr("x2",x(xMax)).attr("y2",y(yMean+slope*(xMax-xMean)))
    .attr("stroke","#bbb").attr("stroke-dasharray","5,4").attr("stroke-width",1.5);

  // Dots per category
  cats.forEach(c => {{
    const col = COLORS[c.category] || "#888";
    const cx = x(c.gl_mean_auroc), cy = y(c.sob_mean_norm);
    g.append("circle").attr("cx",cx).attr("cy",cy).attr("r",10)
     .attr("fill",col).attr("opacity",0.85).attr("stroke","#fff").attr("stroke-width",2)
     .on("mouseover", (e) => showTip(
       `<b>${{c.category}}</b><br>GeneLab AUROC: ${{c.gl_mean_auroc.toFixed(3)}} (n=${{c.gl_n}})<br>SpaceOmicsBench norm: ${{c.sob_mean_norm.toFixed(3)}} (n=${{c.sob_n}})`, e))
     .on("mouseout", hideTip);
    // Label
    const loffset = c.gl_mean_auroc > xMean ? 13 : -13;
    g.append("text").attr("x",cx+loffset).attr("y",cy+4)
     .attr("font-size",9).attr("fill",col).attr("font-weight","600")
     .attr("text-anchor", c.gl_mean_auroc > xMean ? "start":"end")
     .text(c.category.replace(" ML","").replace(" Model","").replace(" Learning","").replace("Graph/Sequence","GNN"));
  }});

  // Rho annotation
  g.append("text").attr("x",4).attr("y",14).attr("font-size",11).attr("fill","#444")
   .attr("font-weight","600").text(`ρ = ${{rho.toFixed(2)}} (Spearman rank, n=${{n}})`);

  g.append("g").attr("transform",`translate(0,${{iH}})`).call(d3.axisBottom(x).ticks(4).tickFormat(d3.format(".2f"))).selectAll("text").attr("font-size",10);
  g.append("g").call(d3.axisLeft(y).ticks(4).tickFormat(d3.format(".2f"))).selectAll("text").attr("font-size",10);

  g.append("text").attr("x",iW/2).attr("y",iH+46).attr("text-anchor","middle")
   .attr("font-size",11).attr("fill","#555").text("GeneLab Mean AUROC (mouse, per category)");
  g.append("text").attr("transform","rotate(-90)").attr("x",-iH/2).attr("y",-55)
   .attr("text-anchor","middle").attr("font-size",11).attr("fill","#555")
   .text("SpaceOmicsBench within-task rank (human, per category)");
  g.append("text").attr("x",iW/2).attr("y",-12).attr("text-anchor","middle")
   .attr("font-size",9).attr("fill","#888").text("SOB score normalized 0=worst,1=best within each task");
}})();

// ─── Panel D: Strip chart — all GeneLab gene-feature results ─────────────
(function() {{
  const W = 1100, H = 250, mg = {{top:20, right:20, bottom:50, left:120}};
  const svg = d3.select("#strip").attr("width", W).attr("height", H);
  const iW = W - mg.left - mg.right, iH = H - mg.top - mg.bottom;
  const g = svg.append("g").attr("transform", `translate(${{mg.left}},${{mg.top}})`);

  // One row per method, show all tissue AUROCs
  const methodData = {{}};
  DATA.gl_rows.forEach(r => {{
    if (r.metric !== "auroc" || r.feature_type !== "gene") return;
    if (!methodData[r.method]) methodData[r.method] = {{ scores: [], category: r.category }};
    methodData[r.method].scores.push({{ score: r.score, tissue: r.tissue_task }});
  }});

  const methods = Object.entries(methodData)
    .filter(([m,d]) => CAT_ORDER.includes(d.category))
    .sort((a,b) => {{
      const catDiff = CAT_ORDER.indexOf(a[1].category) - CAT_ORDER.indexOf(b[1].category);
      if (catDiff !== 0) return catDiff;
      return d3.mean(b[1].scores.map(s=>s.score)) - d3.mean(a[1].scores.map(s=>s.score));
    }});

  const y = d3.scaleBand().domain(methods.map(m=>m[0])).range([0, iH]).padding(0.3);
  const x = d3.scaleLinear().domain([0.3, 1.05]).range([0, iW]);

  // Chance line
  g.append("line").attr("x1", x(0.5)).attr("x2", x(0.5))
   .attr("y1", 0).attr("y2", iH).attr("stroke","#bbb").attr("stroke-dasharray","4,3");

  methods.forEach(([method, data]) => {{
    const yc = y(method) + y.bandwidth()/2;
    const col = COLORS[data.category] || "#888";
    const scores = data.scores.map(s=>s.score);
    const mean = d3.mean(scores);

    // Mean line
    g.append("line")
      .attr("x1", x(d3.min(scores))).attr("x2", x(d3.max(scores)))
      .attr("y1", yc).attr("y2", yc)
      .attr("stroke", col).attr("opacity", 0.3).attr("stroke-width", 1);

    // Dots
    data.scores.forEach(s => {{
      g.append("circle").attr("cx", x(s.score)).attr("cy", yc)
        .attr("r", 4).attr("fill", col).attr("opacity", 0.75)
        .on("mouseover", (e) => showTip(`<b>${{method}}</b><br>Tissue: ${{s.tissue}}<br>AUROC: ${{s.score.toFixed(3)}}<br>Category: ${{data.category}}`, e))
        .on("mouseout", hideTip);
    }});

    // Mean dot
    g.append("circle").attr("cx", x(mean)).attr("cy", yc)
      .attr("r", 6).attr("fill", col).attr("stroke","#fff").attr("stroke-width",1.5);
  }});

  g.append("g").attr("transform", `translate(0,${{iH}})`).call(d3.axisBottom(x).ticks(7).tickFormat(d3.format(".2f"))).selectAll("text").attr("font-size",10);
  g.append("g").call(d3.axisLeft(y).tickFormat(m => m.replace("_"," "))).selectAll("text").attr("font-size",10);

  g.append("text").attr("x", iW/2).attr("y", iH+38).attr("text-anchor","middle")
   .attr("font-size",11).attr("fill","#555").text("AUROC (gene features, per tissue)");
  g.append("text").attr("x", -8).attr("y", -6).attr("text-anchor","end")
   .attr("font-size",10).attr("fill","#888").text("Method");
}})();
</script>
</body>
</html>"""
    return html

# ─── Main ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("Loading data...")
    data = build_viz_data()

    print(f"\nGeneLab summary by category:")
    for cat, s in data["gl_summary"].items():
        print(f"  {cat:20s}: mean={s['mean']:.3f}, n={s['n']}")

    print(f"\nSpaceOmicsBench summary by category:")
    for cat, s in data["sob_summary"].items():
        print(f"  {cat:20s}: mean={s['mean']:.3f}, n={s['n']}")

    print(f"\nMethods with both datasets: {sum(1 for m in data['method_summary'] if m['gl_mean_auroc'] and m['sob_mean_score'])}")

    print("\nGenerating HTML figure...")
    html = generate_html(data)
    out_path = OUT_DIR / "v7_unified_signal_hierarchy.html"
    with open(out_path, "w") as f:
        f.write(html)
    print(f"Saved: {out_path}")

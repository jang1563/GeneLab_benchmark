#!/usr/bin/env python3
"""
rrrm1_f2c_loao_classifier.py  —  F2-C: Per-cell-type LOAO spaceflight classifier

Leave-One-Animal-Out (LOAO) cross-validation of PCA-LR binary classifier
(FLT vs GC) for each tissue × cell_type combination.

LOAO structure (8 animals total: 4 FLT + 4 GC):
  - 8 folds, each fold leaves out 1 animal
  - Train on 7 animals, test on 1 → aggregate predictions → single AUROC

Usage:
    python3 rrrm1_f2c_loao_classifier.py --all
    python3 rrrm1_f2c_loao_classifier.py --tissue muscle --verbose

Inputs:
    *_hardened.h5ad  (obs: condition, animal_id, broad_celltype; log1p normalized HVGs)

Outputs:
    v2/evaluation/F2C_loao_classifier.json
    v2/figures/F2C_loao_auroc.html
"""

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────
SCRATCH = Path("/athena/masonlab/scratch/users/jak4013/rrrm1_scrna")
HARDENED_DIR = SCRATCH / "downstream_initial" / "hardening" / "objects"
BASE_DIR = Path("/home/fs01/jak4013/rrrm1_scrna")
EVAL_DIR = BASE_DIR / "evaluation"
FIG_DIR = BASE_DIR / "figures"

# v1.0 bulk AUROC reference (tissue → mean LOMO AUROC from A-series)
V1_BULK_AUROC = {
    "eye":    0.700,   # A6 eye (approx from v1 results)
    "muscle": 0.758,   # gastrocnemius baseline
    "skin":   0.821,   # A5 skin LOMO mean
    "blood":  None,    # no v1 bulk blood (PBMC only)
}

TISSUE_OSD = {
    "blood":  "OSD-918",
    "eye":    "OSD-920",
    "muscle": "OSD-924",
    "skin":   "OSD-934",
}

# Minimum cells per condition (total across all animals) to include a cell type
MIN_CELLS_PER_CONDITION = 20
# Minimum PCA components
MIN_PCA = 10
MAX_PCA = 50
# Variance percentile to filter (keep top X% variable genes; 0.25 = top 75%)
VARIANCE_FILTER_PERCENTILE = 0.25


# ── Statistical utilities ──────────────────────────────────────────────────

def bootstrap_auroc_ci(y_true: np.ndarray, y_score: np.ndarray,
                       n_boot: int = 2000, seed: int = 42) -> tuple:
    """Percentile-bootstrap 95% CI for AUROC."""
    rng = np.random.default_rng(seed)
    n = len(y_true)
    boot_aurocs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yt = y_true[idx]
        ys = y_score[idx]
        if len(np.unique(yt)) < 2:
            continue
        try:
            boot_aurocs.append(roc_auc_score(yt, ys))
        except Exception:
            pass
    if len(boot_aurocs) < 50:
        return (float("nan"), float("nan"))
    return (float(np.percentile(boot_aurocs, 2.5)),
            float(np.percentile(boot_aurocs, 97.5)))


def permutation_test_auroc(y_true: np.ndarray, y_score: np.ndarray,
                           n_perm: int = 5000, seed: int = 42) -> float:
    """One-sided permutation p-value: P(AUROC_perm >= AUROC_obs)."""
    rng = np.random.default_rng(seed)
    try:
        observed = roc_auc_score(y_true, y_score)
    except Exception:
        return float("nan")
    count_ge = 0
    for _ in range(n_perm):
        y_perm = rng.permutation(y_true)
        try:
            perm_auroc = roc_auc_score(y_perm, y_score)
            if perm_auroc >= observed:
                count_ge += 1
        except Exception:
            pass
    return float((count_ge + 1) / (n_perm + 1))


# ── PCA-LR classifier ─────────────────────────────────────────────────────

def pca_lr_predict(X_train: np.ndarray, y_train: np.ndarray,
                   X_test: np.ndarray) -> np.ndarray:
    """
    StandardScaler → variance filter → PCA(50) → LogisticRegression(balanced).
    Returns predicted probabilities for the positive class (FLT=1).
    """
    # Variance filter on train set only
    var = np.var(X_train, axis=0)
    threshold = np.percentile(var, VARIANCE_FILTER_PERCENTILE * 100)
    keep = var > threshold
    if keep.sum() < 50:
        keep = np.ones(X_train.shape[1], dtype=bool)

    X_tr = X_train[:, keep]
    X_te = X_test[:, keep]

    # Standardize
    scaler = StandardScaler()
    X_tr = scaler.fit_transform(X_tr)
    X_te = scaler.transform(X_te)

    # PCA
    n_comp = min(MAX_PCA, X_tr.shape[0] - 1, X_tr.shape[1])
    n_comp = max(n_comp, MIN_PCA)
    if n_comp < 2:
        return np.full(len(X_te), 0.5)

    pca = PCA(n_components=n_comp, random_state=42)
    X_tr = pca.fit_transform(X_tr)
    X_te = pca.transform(X_te)

    # Logistic Regression
    clf = LogisticRegression(
        solver="lbfgs",
        class_weight="balanced",
        max_iter=2000,
        C=1.0,
        random_state=42,
    )
    try:
        clf.fit(X_tr, y_train)
        y_score = clf.predict_proba(X_te)[:, 1]
    except Exception as e:
        y_score = np.full(len(X_te), 0.5)
    return y_score


# ── LOAO cross-validation ─────────────────────────────────────────────────

def loao_cv(adata_sub: sc.AnnData, verbose: bool = False) -> dict:
    """
    Leave-One-Animal-Out CV on a cell-type subset.
    Returns dict with AUROC, CI, p-value, and per-fold details.
    """
    obs = adata_sub.obs
    animals = obs["animal_id"].unique()

    # Dense feature matrix
    X = adata_sub.X
    if sp.issparse(X):
        X = X.toarray()
    X = X.astype(np.float32)

    y = (obs["condition"].values == "FLT").astype(int)

    # Check class balance
    n_flt = y.sum()
    n_gc = len(y) - n_flt
    if n_flt < MIN_CELLS_PER_CONDITION or n_gc < MIN_CELLS_PER_CONDITION:
        return None

    all_y_true = []
    all_y_score = []
    fold_details = []

    for animal_out in animals:
        train_mask = (obs["animal_id"] != animal_out).values
        test_mask = ~train_mask

        X_train, X_test = X[train_mask], X[test_mask]
        y_train, y_test = y[train_mask], y[test_mask]

        if len(X_test) == 0 or len(np.unique(y_train)) < 2:
            if verbose:
                print(f"      Fold {animal_out}: skipped (no test cells or single class in train)")
            continue

        y_score = pca_lr_predict(X_train, y_train, X_test)

        all_y_true.extend(y_test.tolist())
        all_y_score.extend(y_score.tolist())

        # fold_auroc is NaN when test set has only one class (expected in LOAO:
        # each fold holds out 1 animal = 1 condition). Overall AUROC pools all
        # holdout predictions, restoring both classes.
        fold_auroc = float("nan")
        if len(np.unique(y_test)) == 2:
            try:
                fold_auroc = roc_auc_score(y_test, y_score)
            except Exception:
                pass

        fold_details.append({
            "animal_out": animal_out,
            "condition_out": obs[obs["animal_id"] == animal_out]["condition"].iloc[0],
            "n_test": int(len(y_test)),
            "n_train": int(len(y_train)),
            "fold_auroc": fold_auroc,
        })

        if verbose:
            print(f"      Fold {animal_out}: n_test={len(y_test)} AUROC={fold_auroc:.3f}")

    if len(all_y_true) < 10 or len(np.unique(all_y_true)) < 2:
        return None

    all_y_true = np.array(all_y_true)
    all_y_score = np.array(all_y_score)

    try:
        auroc = float(roc_auc_score(all_y_true, all_y_score))
    except Exception:
        return None

    ci_low, ci_high = bootstrap_auroc_ci(all_y_true, all_y_score)
    p_perm = permutation_test_auroc(all_y_true, all_y_score)

    return {
        "auroc": auroc,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "p_perm": p_perm,
        "n_animals": int(len(animals)),
        "n_cells_flt": int(n_flt),
        "n_cells_gc": int(n_gc),
        "n_cells_total": int(len(y)),
        "fold_details": fold_details,
    }


# ── Per-tissue analysis ────────────────────────────────────────────────────

def run_tissue(tissue: str, verbose: bool = False) -> dict:
    path = HARDENED_DIR / f"RRRM1_{tissue}_hardened.h5ad"
    if not path.exists():
        raise FileNotFoundError(f"Not found: {path}")

    print(f"\n[F2-C] {tissue}")
    adata = sc.read_h5ad(path)

    required = ["condition", "animal_id", "broad_celltype"]
    missing = [c for c in required if c not in adata.obs.columns]
    if missing:
        raise ValueError(f"[{tissue}] Missing obs columns: {missing}")

    print(f"  {adata.n_obs} cells, {adata.n_vars} features")
    print(f"  Conditions: {dict(adata.obs['condition'].value_counts())}")
    print(f"  Animals: {adata.obs['animal_id'].nunique()}")
    print(f"  Cell types: {sorted(adata.obs['broad_celltype'].unique())}")

    tissue_results = {}
    cell_types = [ct for ct in adata.obs["broad_celltype"].unique() if ct != "unknown"]

    for ct in sorted(cell_types):
        sub = adata[adata.obs["broad_celltype"] == ct]
        n_flt = (sub.obs["condition"] == "FLT").sum()
        n_gc = (sub.obs["condition"] == "GC").sum()

        if n_flt < MIN_CELLS_PER_CONDITION or n_gc < MIN_CELLS_PER_CONDITION:
            print(f"  SKIP {ct}: n_FLT={n_flt} n_GC={n_gc}")
            continue

        print(f"  {ct}: {len(sub)} cells (FLT={n_flt}, GC={n_gc})", end=" ... ", flush=True)
        result = loao_cv(sub, verbose=verbose)

        if result is None:
            print("insufficient data")
            continue

        print(f"AUROC={result['auroc']:.3f} [{result['ci_low']:.3f}, {result['ci_high']:.3f}] "
              f"p={result['p_perm']:.3f}")
        tissue_results[ct] = result

    return tissue_results


# ── HTML figure ────────────────────────────────────────────────────────────

def make_auroc_html(all_results: dict) -> str:
    """Grouped bar chart per tissue, cell type × AUROC with CI error bars."""
    # Flatten for D3
    plot_data = []
    for tissue, ct_dict in all_results.items():
        v1_auroc = V1_BULK_AUROC.get(tissue)
        for ct, res in ct_dict.items():
            plot_data.append({
                "tissue": tissue,
                "cell_type": ct,
                "auroc": res["auroc"],
                "ci_low": res["ci_low"],
                "ci_high": res["ci_high"],
                "p_perm": res["p_perm"],
                "n_cells": res["n_cells_total"],
                "significant": res["p_perm"] < 0.05,
            })
        # Add v1.0 bulk as reference line
        if v1_auroc is not None:
            plot_data.append({
                "tissue": tissue,
                "cell_type": "_bulk_v1",
                "auroc": v1_auroc,
                "ci_low": v1_auroc,
                "ci_high": v1_auroc,
                "p_perm": None,
                "n_cells": 0,
                "significant": False,
                "is_bulk": True,
            })

    data_json = json.dumps(plot_data)
    tissues = list(all_results.keys())

    return f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>F2-C RRRM-1 LOAO Classifier AUROC per Cell Type</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; font-size: 11px; margin: 20px; background: #fff; }}
  h2 {{ font-size: 13px; font-weight: bold; }}
  .panel {{ display: inline-block; vertical-align: top; margin: 15px; }}
  .panel-title {{ font-weight: bold; font-size: 12px; margin-bottom: 4px; }}
  .axis text {{ font-size: 9px; }}
  .bar-sig {{ fill: #0072B2; }}
  .bar-ns {{ fill: #999999; }}
  .ref-line {{ stroke: #D55E00; stroke-width: 1.5; stroke-dasharray: 4,3; }}
</style>
</head>
<body>
<h2>F2-C: RRRM-1 scRNA-seq LOAO Classifier AUROC (FLT vs GC, per cell type)</h2>
<p style="font-size:10px;color:#666">
  Blue = p &lt; 0.05 (permutation). Orange dashed = v1.0 bulk AUROC reference.
  Error bars = 95% bootstrap CI.
</p>
<div id="charts"></div>
<script>
const allData = {data_json};
const tissues = {json.dumps(tissues)};

const W = 320, H = 200;
const margin = {{top: 20, right: 20, bottom: 90, left: 55}};
const innerW = W - margin.left - margin.right;
const innerH = H - margin.top - margin.bottom;

const container = d3.select("#charts");

tissues.forEach(tissue => {{
  const tData = allData.filter(d => d.tissue === tissue && !d.is_bulk);
  const bulkRef = allData.find(d => d.tissue === tissue && d.is_bulk);
  if (tData.length === 0) return;

  const cellTypes = tData.map(d => d.cell_type);
  const x = d3.scaleBand().domain(cellTypes).range([0, innerW]).padding(0.3);
  const y = d3.scaleLinear().domain([0, 1]).range([innerH, 0]);

  const div = container.append("div").attr("class","panel");
  div.append("div").attr("class","panel-title")
    .text(tissue.charAt(0).toUpperCase() + tissue.slice(1));
  const svg = div.append("svg").attr("width",W).attr("height",H);
  const g = svg.append("g").attr("transform",`translate(${{margin.left}},${{margin.top}})`);

  // Chance line at 0.5
  g.append("line")
    .attr("x1",0).attr("x2",innerW).attr("y1",y(0.5)).attr("y2",y(0.5))
    .attr("stroke","#ccc").attr("stroke-dasharray","3,2");

  // v1.0 bulk reference line
  if (bulkRef) {{
    g.append("line")
      .attr("class","ref-line")
      .attr("x1",0).attr("x2",innerW)
      .attr("y1",y(bulkRef.auroc)).attr("y2",y(bulkRef.auroc));
    g.append("text")
      .attr("x",innerW-2).attr("y",y(bulkRef.auroc)-3)
      .attr("text-anchor","end").style("font-size","8px").style("fill","#D55E00")
      .text("bulk v1");
  }}

  // Bars
  g.selectAll(".bar")
    .data(tData)
    .join("rect")
    .attr("class", d => d.significant ? "bar-sig" : "bar-ns")
    .attr("x", d => x(d.cell_type))
    .attr("y", d => y(Math.max(d.auroc, 0.5)))
    .attr("width", x.bandwidth())
    .attr("height", d => Math.abs(y(d.auroc) - y(0.5)))
    .attr("fill", d => d.significant ? "#0072B2" : "#999");

  // Error bars (CI)
  g.selectAll(".ci")
    .data(tData)
    .join("line")
    .attr("x1", d => x(d.cell_type) + x.bandwidth()/2)
    .attr("x2", d => x(d.cell_type) + x.bandwidth()/2)
    .attr("y1", d => y(isNaN(d.ci_high) ? d.auroc : d.ci_high))
    .attr("y2", d => y(isNaN(d.ci_low) ? d.auroc : d.ci_low))
    .attr("stroke", "#333").attr("stroke-width", 1);

  // AUROC labels
  g.selectAll(".lbl")
    .data(tData)
    .join("text")
    .attr("x", d => x(d.cell_type) + x.bandwidth()/2)
    .attr("y", d => y(d.auroc) - 3)
    .attr("text-anchor","middle").style("font-size","8px")
    .text(d => d.auroc.toFixed(2));

  // Axes
  g.append("g").attr("transform",`translate(0,${{innerH}})`).call(
    d3.axisBottom(x).tickSize(2)
  ).selectAll("text")
    .attr("transform","rotate(-40)").style("text-anchor","end").style("font-size","8px");

  g.append("g").call(
    d3.axisLeft(y).ticks(5).tickFormat(d3.format(".2f")).tickSize(2)
  );

  // Y label
  svg.append("text").attr("transform","rotate(-90)")
    .attr("x",-(H/2)).attr("y",14).attr("text-anchor","middle")
    .style("font-size","9px").text("AUROC");
}});
</script>
</body>
</html>
"""


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="F2-C LOAO classifier per cell type")
    parser.add_argument("--tissue", choices=list(TISSUE_OSD.keys()))
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--verbose", action="store_true", help="Print per-fold AUROC")
    args = parser.parse_args()

    if not args.tissue and not args.all:
        parser.error("Specify --tissue or --all")

    EVAL_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    tissues = list(TISSUE_OSD.keys()) if args.all else [args.tissue]

    all_results = {}
    for tissue in tissues:
        try:
            res = run_tissue(tissue, verbose=args.verbose)
            all_results[tissue] = res
        except FileNotFoundError as e:
            print(f"  SKIP {tissue}: {e}")
        except ValueError as e:
            print(f"  ERROR {tissue}: {e}")

    if not all_results:
        print("No tissues processed.")
        sys.exit(1)

    # Summary
    print("\n=== F2-C Summary ===")
    print(f"{'Tissue':<10} {'Cell type':<35} {'AUROC':>6} {'CI':>16} {'p':>8}")
    print("-" * 75)
    for tissue, ct_dict in sorted(all_results.items()):
        v1 = V1_BULK_AUROC.get(tissue, "N/A")
        for ct, res in sorted(ct_dict.items(), key=lambda x: -x[1]["auroc"]):
            sig = "*" if res["p_perm"] < 0.05 else " "
            ci_str = f"[{res['ci_low']:.2f}, {res['ci_high']:.2f}]"
            print(f"{tissue:<10} {ct:<35} {res['auroc']:.3f}{sig} {ci_str:>16} {res['p_perm']:.3f}")
        print(f"{'':>10} {'v1.0 bulk baseline':>35} {str(v1):>6}")
        print()

    # Save JSON
    out_json = EVAL_DIR / "F2C_loao_classifier.json"
    output = {
        "task": "F2-C",
        "method": "LOAO PCA(50)-LR(balanced), bootstrap CI n=2000, permutation p n=5000",
        "results": all_results,
        "v1_bulk_auroc": {k: v for k, v in V1_BULK_AUROC.items() if v is not None},
    }
    with open(out_json, "w") as f:
        json.dump(output, f, indent=2)
    print(f"Saved: {out_json}")

    # Save HTML
    out_html = FIG_DIR / "F2C_loao_auroc.html"
    out_html.write_text(make_auroc_html(all_results))
    print(f"Saved: {out_html}")


if __name__ == "__main__":
    main()

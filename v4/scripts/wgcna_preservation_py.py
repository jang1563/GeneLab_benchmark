#!/usr/bin/env python3
"""
wgcna_preservation_py.py — Module preservation analysis across 15 tissue pairs.

Python implementation of WGCNA modulePreservation() core statistics.
For each (reference, test) tissue pair, computes:
  - Z_density: Z-score of intra-module co-expression density in test tissue
  - Z_connectivity: Z-score of intra-module mean adjacency in test tissue
  - Zsummary: composite preservation score (avg of Z_density, Z_connectivity)
  - medianRank: rank of module preservation vs random same-size gene sets

Interpretation (same as R WGCNA):
  Zsummary > 10  → strongly preserved
  2 < Zsummary < 10 → moderately preserved
  Zsummary < 2  → not preserved (or underpowered)

Input:  v4/wgcna_inputs/{tissue}_expr.csv    (genes × samples)
        v4/wgcna_outputs/{tissue}/module_assignments.csv

Output: v4/evaluation/WGCNA_preservation.json

Usage:
  python wgcna_preservation_py.py
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from itertools import combinations

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR    = Path(__file__).resolve().parent.parent.parent
INPUTS_DIR  = BASE_DIR / "v4" / "wgcna_inputs"
OUTPUTS_DIR = BASE_DIR / "v4" / "wgcna_outputs"
EVAL_DIR    = BASE_DIR / "v4" / "evaluation"

LOMO_TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin"]
N_PERM       = 100    # permutations per module for Z-score
RANDOM_SEED  = 42
MIN_MODULE_GENES = 10  # minimum genes in module for preservation test


def load_expr(tissue):
    """Load expression matrix: returns DataFrame (genes × samples)."""
    expr_path = INPUTS_DIR / f"{tissue}_expr.csv"
    expr = pd.read_csv(expr_path, index_col=0)  # genes × samples
    return expr


def load_modules(tissue):
    """Load module assignments. Returns dict: gene → module_color."""
    mod_path = OUTPUTS_DIR / tissue / "module_assignments.csv"
    df = pd.read_csv(mod_path)
    return dict(zip(df["gene"], df["module_color"]))


def compute_mean_cor(gene_list, expr_genes_x_samples):
    """Mean pairwise Pearson correlation of genes in expr (genes × samples).

    Returns the upper-triangle mean (excluding diagonal).
    """
    valid = [g for g in gene_list if g in expr_genes_x_samples.index]
    if len(valid) < 2:
        return np.nan

    sub = expr_genes_x_samples.loc[valid].values.astype(np.float64)
    # Compute correlation matrix
    C = np.corrcoef(sub)  # n_valid × n_valid
    n = len(valid)
    if n < 2:
        return np.nan
    # Mean of upper triangle (excluding diagonal)
    triu_vals = C[np.triu_indices(n, k=1)]
    return float(np.nanmean(triu_vals))


def compute_mean_adj(gene_list, expr_genes_x_samples, beta=8):
    """Mean pairwise adjacency (soft-thresholded) for gene_list in expr."""
    valid = [g for g in gene_list if g in expr_genes_x_samples.index]
    if len(valid) < 2:
        return np.nan

    sub = expr_genes_x_samples.loc[valid].values.astype(np.float64)
    C = np.corrcoef(sub)
    adj = np.power(np.clip(0.5 + 0.5 * C, 0.0, 1.0), beta)
    n = len(valid)
    np.fill_diagonal(adj, 0.0)
    triu_vals = adj[np.triu_indices(n, k=1)]
    return float(np.nanmean(triu_vals))


def preservation_stats(module_genes, test_expr, all_test_genes, beta=8, n_perm=100, seed=42):
    """Compute preservation Z-scores for one module in test tissue.

    module_genes: list of gene names in the module
    test_expr: genes × samples DataFrame of test tissue
    all_test_genes: list of all genes available in test tissue (for random sampling)
    beta: soft threshold for adjacency

    Returns: {Z_density, Z_connectivity, Zsummary, medianRank_density, n_module_genes_available}
    """
    rng = np.random.default_rng(seed)

    # Genes available in test tissue
    avail = [g for g in module_genes if g in test_expr.index]
    k = len(avail)

    if k < 2:
        return {
            "Z_density": np.nan, "Z_connectivity": np.nan, "Zsummary": np.nan,
            "medianRank_density": np.nan, "n_genes_in_ref": len(module_genes),
            "n_genes_in_test": k, "obs_density": np.nan, "obs_connectivity": np.nan
        }

    # Observed density and connectivity
    obs_density = compute_mean_cor(avail, test_expr)
    obs_connectivity = compute_mean_adj(avail, test_expr, beta=beta)

    # Null distribution: n_perm random gene sets of same size k
    background_genes = [g for g in all_test_genes if g not in set(avail)]
    if len(background_genes) < k:
        # Not enough background — use all available genes
        background_genes = list(test_expr.index)

    null_density = []
    null_connectivity = []
    for _ in range(n_perm):
        rand_genes = rng.choice(background_genes, size=k, replace=False).tolist()
        nd = compute_mean_cor(rand_genes, test_expr)
        na = compute_mean_adj(rand_genes, test_expr, beta=beta)
        null_density.append(nd)
        null_connectivity.append(na)

    null_density = np.array([v for v in null_density if not np.isnan(v)])
    null_connectivity = np.array([v for v in null_connectivity if not np.isnan(v)])

    def z_score(obs, null_arr):
        if len(null_arr) < 2:
            return np.nan
        mu, sd = null_arr.mean(), null_arr.std()
        if sd < 1e-10:
            return np.nan
        return float((obs - mu) / sd)

    Z_d = z_score(obs_density, null_density)
    Z_c = z_score(obs_connectivity, null_connectivity)

    # Zsummary: average of available Z scores
    z_vals = [z for z in [Z_d, Z_c] if not np.isnan(z)]
    Zsummary = float(np.mean(z_vals)) if z_vals else np.nan

    # medianRank for density: rank of obs in null_density (lower = better preserved)
    if len(null_density) > 0 and not np.isnan(obs_density):
        rank = float(np.sum(null_density >= obs_density))  # how many nulls ≥ obs
    else:
        rank = np.nan

    return {
        "Z_density": round(Z_d, 3) if not np.isnan(Z_d) else None,
        "Z_connectivity": round(Z_c, 3) if not np.isnan(Z_c) else None,
        "Zsummary": round(Zsummary, 3) if not np.isnan(Zsummary) else None,
        "medianRank_density": round(rank, 1) if not np.isnan(rank) else None,
        "n_genes_in_ref": len(module_genes),
        "n_genes_in_test": k,
        "obs_density": round(obs_density, 4) if not np.isnan(obs_density) else None,
        "obs_connectivity": round(obs_connectivity, 4) if not np.isnan(obs_connectivity) else None,
    }


def run_preservation_pair(ref_tissue, test_tissue,
                           ref_modules, test_expr,
                           all_test_genes, beta=8):
    """Compute preservation for all modules in ref_tissue in test_tissue."""
    # Group ref genes by module
    modules = {}
    for gene, mod in ref_modules.items():
        if mod == "grey":
            continue
        modules.setdefault(mod, []).append(gene)

    pair_results = {}
    for mod_color, genes in sorted(modules.items()):
        if len(genes) < MIN_MODULE_GENES:
            continue

        seed = abs(hash(f"{ref_tissue}_{test_tissue}_{mod_color}")) % (2**31)
        stats = preservation_stats(genes, test_expr, all_test_genes,
                                    beta=beta, n_perm=N_PERM, seed=seed)
        pair_results[mod_color] = stats

        zsumm = stats["Zsummary"]
        flag = "**" if (zsumm is not None and zsumm > 10) else \
               "*"  if (zsumm is not None and zsumm > 2)  else ""
        print(f"    {mod_color:20s} k={stats['n_genes_in_ref']:4d} "
              f"avail={stats['n_genes_in_test']:4d} "
              f"Z={zsumm:.2f}{flag}" if zsumm is not None
              else f"    {mod_color:20s} k={stats['n_genes_in_ref']:4d} avail={stats['n_genes_in_test']:4d} Z=NA")

    return pair_results


def main():
    EVAL_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading expression matrices for 6 tissues...")
    expr_data   = {}
    module_data = {}
    for tissue in LOMO_TISSUES:
        expr_path = INPUTS_DIR / f"{tissue}_expr.csv"
        mod_path  = OUTPUTS_DIR / tissue / "module_assignments.csv"
        if not expr_path.exists():
            print(f"  WARNING: {tissue} expression not found, skipping")
            continue
        if not mod_path.exists():
            print(f"  WARNING: {tissue} module assignments not found, skipping")
            continue
        expr_data[tissue]   = pd.read_csv(expr_path, index_col=0)  # genes × samples
        module_data[tissue] = load_modules(tissue)
        n_genes = len(expr_data[tissue])
        n_mods  = len(set(v for v in module_data[tissue].values() if v != "grey"))
        print(f"  {tissue}: {n_genes} genes, {n_mods} modules")

    available_tissues = [t for t in LOMO_TISSUES if t in expr_data]
    pairs = list(combinations(available_tissues, 2))
    print(f"\nComputing preservation for C({len(available_tissues)},2)={len(pairs)} tissue pairs "
          f"(n_perm={N_PERM})...")

    all_results = {}
    pair_summary = []

    for ref_tissue, test_tissue in pairs:
        print(f"\n{'─'*60}")
        print(f"Reference: {ref_tissue}  →  Test: {test_tissue}")

        ref_modules  = module_data[ref_tissue]
        test_expr    = expr_data[test_tissue]
        all_test_genes = list(test_expr.index)

        # Use test tissue's dominant beta (proxy: median common beta)
        beta = 9  # default; could read from summary.json

        result = run_preservation_pair(ref_tissue, test_tissue,
                                        ref_modules, test_expr,
                                        all_test_genes, beta=beta)

        key = f"{ref_tissue}_vs_{test_tissue}"
        all_results[key] = result

        # Summary stats
        zs = [v["Zsummary"] for v in result.values() if v["Zsummary"] is not None]
        n_high   = sum(1 for z in zs if z > 10)
        n_mod    = sum(1 for z in zs if 2 < z <= 10)
        n_notpres = sum(1 for z in zs if z <= 2)
        mean_z   = float(np.mean(zs)) if zs else np.nan
        pair_summary.append({
            "pair": key, "n_modules": len(result),
            "n_high_pres": n_high, "n_moderate_pres": n_mod, "n_not_pres": n_notpres,
            "mean_Zsummary": round(mean_z, 2) if not np.isnan(mean_z) else None
        })
        print(f"  Summary: {len(result)} modules tested | "
              f"high(Z>10)={n_high}, moderate(2<Z≤10)={n_mod}, not(Z≤2)={n_notpres}, "
              f"mean_Z={mean_z:.2f}" if not np.isnan(mean_z) else "  Summary: no valid Z scores")

    # Save
    output = {
        "metadata": {
            "n_perm": N_PERM,
            "min_module_size": MIN_MODULE_GENES,
            "tissues": available_tissues,
            "n_pairs": len(pairs),
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "note": "Python implementation: Z_density + Z_connectivity → Zsummary"
        },
        "pair_summary": pair_summary,
        "preservation": all_results
    }

    out_path = EVAL_DIR / "WGCNA_preservation.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=lambda x: None if (isinstance(x, float) and np.isnan(x)) else x)
    print(f"\nSaved: {out_path}")

    # Print summary table
    print(f"\n{'═'*70}")
    print(f"{'Pair':<35} {'Modules':>8} {'High':>6} {'Mod':>6} {'Mean Z':>8}")
    print(f"{'─'*70}")
    for s in pair_summary:
        print(f"{s['pair']:<35} {s['n_modules']:>8} {s['n_high_pres']:>6} "
              f"{s['n_moderate_pres']:>6} {str(s['mean_Zsummary']):>8}")


if __name__ == "__main__":
    main()

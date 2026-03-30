#!/usr/bin/env python3
"""
wgcna_analysis_py.py — Per-tissue WGCNA co-expression network analysis.

Pure Python/numpy/scipy implementation (no R dependency).
Algorithm matches WGCNA R package: signed hybrid adjacency → TOM →
hierarchical clustering → dynamic tree cut → module merging → eigengenes → kME.

Key optimizations vs naïve implementation:
  - TOM via matrix multiply (BLAS-accelerated, ~2-3s for n=5000)
  - kME via vectorized Pearson (matrix ops, fast)
  - soft threshold from correlation matrix computed once

Input:  v4/wgcna_inputs/{tissue}_expr.csv   (genes × samples)
        v4/wgcna_inputs/{tissue}_traits.csv  (samples × traits)
Output: v4/wgcna_outputs/{tissue}/
          module_assignments.csv
          eigengenes.csv
          trait_correlations_r.csv / trait_correlations_p.csv
          soft_threshold.csv
          module_sizes.csv
          summary.json

Usage:
  python wgcna_analysis_py.py --tissue liver
  python wgcna_analysis_py.py --tissue gastrocnemius
"""

import argparse
import json
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime
from scipy.stats import pearsonr, linregress
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR    = Path(__file__).resolve().parent.parent.parent
INPUTS_DIR  = BASE_DIR / "v4" / "wgcna_inputs"
OUTPUTS_DIR = BASE_DIR / "v4" / "wgcna_outputs"

# ── Parameters (DD-24) ────────────────────────────────────────────────────────
R2_THRESHOLD     = 0.80
MERGE_CUT_HEIGHT = 0.25
DEEP_SPLIT       = 2
SMALL_N_TISSUES  = {"gastrocnemius", "eye"}

# Standard WGCNA color palette
WGCNA_COLORS = [
    "turquoise", "blue", "brown", "yellow", "green", "red", "black",
    "pink", "magenta", "purple", "greenyellow", "tan", "salmon",
    "cyan", "midnightblue", "lightcyan", "grey60", "lightgreen",
    "lightyellow", "royalblue", "darkred", "darkgreen", "darkturquoise",
    "darkgrey", "orange", "darkorange", "white", "skyblue", "saddlebrown",
    "steelblue", "paleturquoise",
]


# ── Data loading ──────────────────────────────────────────────────────────────

def load_data(tissue):
    """Load expression (genes×samples) and traits (samples×traits)."""
    expr_path   = INPUTS_DIR / f"{tissue}_expr.csv"
    traits_path = INPUTS_DIR / f"{tissue}_traits.csv"

    if not expr_path.exists():
        raise FileNotFoundError(f"Expression file not found: {expr_path}")
    if not traits_path.exists():
        raise FileNotFoundError(f"Traits file not found: {traits_path}")

    # Expression: genes × samples → transpose to samples × genes for WGCNA
    expr = pd.read_csv(expr_path, index_col=0)
    print(f"  Loaded expression: {expr.shape[0]} genes × {expr.shape[1]} samples")
    datExpr = expr.T.astype(float)

    # Traits: samples × traits
    traits = pd.read_csv(traits_path, index_col=0)
    traits = traits.apply(pd.to_numeric, errors="coerce")

    # Align samples
    common = sorted(set(datExpr.index) & set(traits.index))
    if len(common) < 5:
        raise ValueError(f"Too few common samples: {len(common)}")
    datExpr = datExpr.loc[common]
    traits  = traits.loc[common]
    print(f"  Aligned: {len(common)} samples")
    return datExpr, traits


# ── Soft threshold selection ──────────────────────────────────────────────────

def pick_soft_threshold(X, powers=None, is_small_n=False):
    """Compute scale-free topology R² for a range of soft threshold powers.

    X: samples × genes (numpy array, standardized)
    Returns: (chosen_beta, r2_achieved, r2_threshold_met, sft_df)
    """
    if powers is None:
        powers = list(range(1, 21)) + list(range(22, 31, 2))

    print(f"  Computing gene-gene correlation matrix ({X.shape[1]} genes)...")
    cor = np.corrcoef(X.T)  # genes × genes
    np.fill_diagonal(cor, 0.0)

    sft_results = []
    for beta in powers:
        # Signed hybrid adjacency: ((1 + cor) / 2) ^ beta
        adj = np.power(np.clip(0.5 + 0.5 * cor, 0.0, 1.0), beta)
        np.fill_diagonal(adj, 0.0)
        k = adj.sum(axis=1)  # node strength (degree for soft networks)

        # Scale-free topology fit: log(p(k)) ~ log(k)
        k_pos = k[k > 1e-10]
        if len(k_pos) < 10:
            sft_results.append({"Power": beta, "SFT.R.sq": 0.0, "slope": 0.0, "mean.k": float(k.mean())})
            continue

        n_bins = min(100, max(10, len(k_pos) // 5))
        hist, edges = np.histogram(k_pos, bins=n_bins)
        centers = (edges[:-1] + edges[1:]) / 2.0
        valid = hist > 0

        if valid.sum() < 4:
            sft_results.append({"Power": beta, "SFT.R.sq": 0.0, "slope": 0.0, "mean.k": float(k.mean())})
            continue

        log_c = np.log10(centers[valid])
        log_h = np.log10(hist[valid].astype(float))
        sl, ic, r, p, se = linregress(log_c, log_h)
        r2 = r**2 if sl < 0 else 0.0
        sft_results.append({"Power": beta, "SFT.R.sq": round(r2, 4), "slope": round(sl, 4),
                             "mean.k": round(float(k.mean()), 4)})

    sft_df = pd.DataFrame(sft_results)

    # Select first power where R² ≥ threshold AND slope < 0
    good = sft_df[(sft_df["SFT.R.sq"] >= R2_THRESHOLD) & (sft_df["slope"] < 0)]
    if len(good) > 0:
        chosen_beta = int(good.iloc[0]["Power"])
        r2_achieved = float(good.iloc[0]["SFT.R.sq"])
        r2_met = True
        print(f"  Chosen beta={chosen_beta} (R²={r2_achieved:.3f} >= {R2_THRESHOLD})")
    else:
        # Fallback: beta=6 for small-n (n<50), else 8
        chosen_beta = 6 if is_small_n else 8
        r2_achieved = float(sft_df["SFT.R.sq"].max())
        r2_met = False
        print(f"  WARNING: R² never reached {R2_THRESHOLD} (max={r2_achieved:.3f}). "
              f"Using fallback beta={chosen_beta}")

    return chosen_beta, r2_achieved, r2_met, sft_df, cor


# ── TOM computation (vectorized) ──────────────────────────────────────────────

def compute_adjacency(cor, beta):
    """Signed hybrid adjacency matrix."""
    adj = np.power(np.clip(0.5 + 0.5 * cor, 0.0, 1.0), beta, dtype=np.float32)
    np.fill_diagonal(adj, 0.0)
    return adj


def compute_tom(adj):
    """Vectorized Topological Overlap Matrix.

    Uses matrix multiply approximation for numerator:
      TOM(i,j) = (sum_u a[i,u]*a[j,u] + a[i,j]) / (min(k_i,k_j) + 1 - a[i,j])
    where the exact formula uses sum_u min(a[i,u],a[j,u]) but the matrix
    multiply approximation (inner product) gives equivalent clustering results
    for soft-thresholded signed networks and is BLAS-accelerated (~2s for n=5000).

    adj: float32 numpy array (n_genes × n_genes)
    Returns: float32 TOM matrix
    """
    print(f"  TOM computation: {adj.shape[0]}×{adj.shape[0]} (matrix multiply, float32)...")
    # BLAS matrix multiply: L[i,j] = sum_u adj[i,u]*adj[j,u]
    L = adj @ adj  # shape: (n, n), float32, BLAS-accelerated

    k = adj.sum(axis=1)  # node strengths, shape (n,)

    numerator = L + adj  # (n, n)
    denominator = np.minimum(k[:, None], k[None, :]) + 1.0 - adj  # (n, n)

    # Avoid division by zero
    np.clip(denominator, 1e-10, None, out=denominator)

    TOM = numerator / denominator
    np.fill_diagonal(TOM, 0.0)
    np.clip(TOM, 0.0, 1.0, out=TOM)
    return TOM


# ── Module detection ──────────────────────────────────────────────────────────

def detect_modules(TOM, min_module_size, deep_split=2):
    """Hierarchical clustering + dynamic tree cut.

    TOM: float32 (n×n) TOM matrix
    Returns: labels array (int, 0=unassigned)
    """
    dissTOM = 1.0 - TOM
    np.fill_diagonal(dissTOM, 0.0)
    np.clip(dissTOM, 0.0, 1.0, out=dissTOM)

    print(f"  Hierarchical clustering (average linkage)...")
    # squareform converts n×n symmetric matrix to condensed form
    Z = linkage(squareform(dissTOM.astype(np.float64)), method="average")

    print(f"  Dynamic tree cut (deepSplit={deep_split}, minModuleSize={min_module_size})...")
    try:
        from dynamicTreeCut import cutreeHybrid
        cut = cutreeHybrid(Z, dissTOM.astype(np.float64),
                           minClusterSize=min_module_size,
                           deepSplit=deep_split,
                           verbose=0)
        labels = cut["labels"]
        print(f"  cutreeHybrid: {len(set(labels)) - (1 if 0 in labels else 0)} modules "
              f"({sum(labels == 0)} grey)")
    except ImportError:
        # Fallback: scipy fcluster
        from scipy.cluster.hierarchy import fcluster
        labels = fcluster(Z, t=0.75, criterion="distance")
        print(f"  fcluster fallback: {len(set(labels))} modules")

    return labels


def labels_to_colors(labels):
    """Assign WGCNA color names to integer cluster labels (0 → grey)."""
    unique = sorted(set(labels))
    color_map = {0: "grey"}
    color_idx = 0
    for lbl in unique:
        if lbl == 0:
            continue
        color_map[lbl] = WGCNA_COLORS[color_idx % len(WGCNA_COLORS)]
        color_idx += 1
    return np.array([color_map.get(lbl, "grey") for lbl in labels])


# ── Module eigengenes ─────────────────────────────────────────────────────────

def compute_eigengenes(datExpr, module_colors):
    """Compute PC1 of each module's gene expression (module eigengene).

    datExpr: samples × genes DataFrame
    module_colors: array of color strings, one per gene
    Returns: DataFrame (samples × modules), columns named "ME{color}"
    """
    modules = sorted(set(module_colors))
    MEs = pd.DataFrame(index=datExpr.index)

    for mod in modules:
        mask = np.array(module_colors) == mod
        mod_expr = datExpr.values[:, mask]  # samples × module_genes

        if mod_expr.shape[1] < 2:
            MEs[f"ME{mod}"] = mod_expr[:, 0] if mod_expr.shape[1] == 1 else 0.0
            continue

        try:
            pca = PCA(n_components=1, random_state=42)
            pc1 = pca.fit_transform(mod_expr)[:, 0]
            # Flip sign so PC1 correlates positively with average module expression
            if np.corrcoef(pc1, mod_expr.mean(axis=1))[0, 1] < 0:
                pc1 = -pc1
            MEs[f"ME{mod}"] = pc1
        except Exception:
            MEs[f"ME{mod}"] = mod_expr.mean(axis=1)

    return MEs


def merge_close_modules(datExpr, module_colors, cut_height=0.25, max_iter=20):
    """Merge modules whose eigengene Pearson correlation > 1 - cut_height.

    Matches R WGCNA mergeCloseModules() behavior.
    """
    module_colors = list(module_colors)
    modules = [m for m in sorted(set(module_colors)) if m != "grey"]

    if len(modules) <= 1:
        return np.array(module_colors)

    n_before = len(modules)
    for iteration in range(max_iter):
        MEs = compute_eigengenes(datExpr, module_colors)
        me_cols = [f"ME{m}" for m in modules if f"ME{m}" in MEs.columns]

        if len(me_cols) <= 1:
            break

        me_matrix = MEs[me_cols].values.T  # modules × samples
        me_cor = np.corrcoef(me_matrix)
        np.fill_diagonal(me_cor, -np.inf)

        max_cor_val = float(np.max(me_cor))
        if max_cor_val < (1.0 - cut_height):
            break

        i, j = np.unravel_index(np.argmax(me_cor), me_cor.shape)
        col_i = me_cols[i].replace("ME", "")
        col_j = me_cols[j].replace("ME", "")

        # Merge smaller → larger module (keep the color with more genes)
        n_i = sum(1 for c in module_colors if c == col_i)
        n_j = sum(1 for c in module_colors if c == col_j)
        keep = col_i if n_i >= n_j else col_j
        drop = col_j if n_i >= n_j else col_i

        module_colors = [keep if c == drop else c for c in module_colors]
        modules = [m for m in sorted(set(module_colors)) if m != "grey"]
        print(f"    Iter {iteration+1}: merged {drop} → {keep} "
              f"(r={max_cor_val:.3f}, remaining={len(modules)} modules)")

    n_after = len(modules)
    print(f"  Merging: {n_before} → {n_after} modules (cut_height={cut_height})")
    return np.array(module_colors)


# ── Module membership (kME) ───────────────────────────────────────────────────

def compute_kme(datExpr, MEs):
    """Vectorized kME: Pearson correlation of each gene with each ME.

    datExpr: samples × genes
    MEs: samples × modules
    Returns: DataFrame (genes × modules), columns named "kME{color}"
    """
    n_genes   = datExpr.shape[1]
    gene_names = datExpr.columns.tolist()
    me_cols   = MEs.columns.tolist()

    # Stack genes and MEs, compute pairwise Pearson in one matrix op
    # X: genes × samples, Y: modules × samples
    X = datExpr.values.T.astype(np.float64)  # n_genes × n_samples
    Y = MEs.values.T.astype(np.float64)      # n_modules × n_samples

    # Center
    X = X - X.mean(axis=1, keepdims=True)
    Y = Y - Y.mean(axis=1, keepdims=True)

    # Normalize
    X_norm = np.linalg.norm(X, axis=1, keepdims=True)
    Y_norm = np.linalg.norm(Y, axis=1, keepdims=True)
    X_norm = np.where(X_norm > 1e-10, X_norm, 1.0)
    Y_norm = np.where(Y_norm > 1e-10, Y_norm, 1.0)
    X = X / X_norm
    Y = Y / Y_norm

    kme_values = X @ Y.T  # n_genes × n_modules (vectorized Pearson)
    np.clip(kme_values, -1.0, 1.0, out=kme_values)

    # Rename columns from ME{color} to kME{color}
    kme_colnames = [c.replace("ME", "kME") for c in me_cols]
    return pd.DataFrame(kme_values, index=gene_names, columns=kme_colnames)


# ── Module-trait correlations ─────────────────────────────────────────────────

def compute_trait_cors(MEs, traits):
    """Pearson correlation between MEs and traits (pairwise complete obs)."""
    # Filter traits: numeric, >0 variance, ≥5 non-NA
    valid_traits = []
    for col in traits.columns:
        vals = pd.to_numeric(traits[col], errors="coerce")
        if vals.notna().sum() >= 5 and vals.std() > 0:
            valid_traits.append(col)

    if not valid_traits:
        return pd.DataFrame(), pd.DataFrame()

    me_cols = MEs.columns.tolist()
    r_df = pd.DataFrame(index=me_cols, columns=valid_traits, dtype=float)
    p_df = pd.DataFrame(index=me_cols, columns=valid_traits, dtype=float)

    for me_col in me_cols:
        me_v = pd.to_numeric(MEs[me_col], errors="coerce")
        for trait_col in valid_traits:
            tr_v = pd.to_numeric(traits[trait_col], errors="coerce")
            valid = me_v.notna() & tr_v.notna()
            if valid.sum() < 5:
                r_df.loc[me_col, trait_col] = np.nan
                p_df.loc[me_col, trait_col] = np.nan
                continue
            r, p = pearsonr(me_v[valid].values, tr_v[valid].values)
            r_df.loc[me_col, trait_col] = round(r, 4)
            p_df.loc[me_col, trait_col] = round(p, 6)

    return r_df, p_df


# ── Save results ──────────────────────────────────────────────────────────────

def save_results(tissue, datExpr, module_colors, MEs, trait_r, trait_p,
                 kme_df, beta, r2_achieved, r2_met, min_module_size,
                 sft_df, n_raw_modules):
    output_path = OUTPUTS_DIR / tissue
    output_path.mkdir(parents=True, exist_ok=True)

    gene_names = datExpr.columns.tolist()

    # 1. Module assignments + kME
    mod_df = pd.DataFrame({"gene": gene_names, "module_color": module_colors})

    # kME for assigned module
    kme_assigned = []
    for gene, mod in zip(gene_names, module_colors):
        kme_col = f"kME{mod}"
        val = float(kme_df.loc[gene, kme_col]) if (kme_df is not None and
              kme_col in kme_df.columns and gene in kme_df.index) else np.nan
        kme_assigned.append(val)
    mod_df["kME_assigned"] = kme_assigned

    # Full kME matrix
    if kme_df is not None:
        kme_reset = kme_df.reset_index().rename(columns={"index": "gene"})
        mod_df = mod_df.merge(kme_reset, on="gene", how="left")

    mod_df.to_csv(output_path / "module_assignments.csv", index=False)
    print(f"  Wrote: module_assignments.csv ({len(mod_df)} genes)")

    # 2. Eigengenes
    MEs.to_csv(output_path / "eigengenes.csv")

    # 3. Trait correlations
    n_sig = 0
    if not trait_r.empty:
        trait_r.reset_index().rename(columns={"index": "module"}).to_csv(
            output_path / "trait_correlations_r.csv", index=False)
        trait_p.reset_index().rename(columns={"index": "module"}).to_csv(
            output_path / "trait_correlations_p.csv", index=False)
        n_sig = int((trait_p < 0.05).sum().sum())
        print(f"  Wrote: trait_correlations (r + p). Significant (p<0.05): {n_sig}")
    else:
        print(f"  No valid traits for correlation")

    # 4. Module sizes
    sizes = pd.Series(module_colors).value_counts().reset_index()
    sizes.columns = ["module_color", "n_genes"]
    sizes.to_csv(output_path / "module_sizes.csv", index=False)

    # 5. Soft threshold
    if sft_df is not None:
        sft_df.to_csv(output_path / "soft_threshold.csv", index=False)

    # 6. Summary JSON
    non_grey = sorted(set(m for m in module_colors if m != "grey"))
    n_grey   = int(sum(1 for c in module_colors if c == "grey"))

    summary = {
        "tissue": tissue,
        "n_samples": int(datExpr.shape[0]),
        "n_genes_input": int(datExpr.shape[1]),
        "soft_threshold_beta": int(beta),
        "r2_achieved": round(float(r2_achieved), 4) if r2_achieved is not None else None,
        "r2_threshold_met": bool(r2_met),
        "min_module_size": int(min_module_size),
        "n_modules_raw": int(n_raw_modules),
        "n_modules_final": len(non_grey),
        "n_grey_genes": n_grey,
        "modules": non_grey,
        "n_genes_per_module": {
            mod: int(sum(1 for c in module_colors if c == mod))
            for mod in non_grey
        },
        "n_traits": int(len(trait_r.columns)) if not trait_r.empty else 0,
        "n_significant_module_trait": int(n_sig),
        "network_type": "signed hybrid",
        "cor_function": "pearson",
        "merge_cut_height": MERGE_CUT_HEIGHT,
        "deep_split": DEEP_SPLIT,
        "implementation": "manual_numpy_scipy",
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }

    with open(output_path / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"  Wrote: summary.json")
    print(f"\n  Modules: {len(non_grey)} (+ {n_grey} grey genes)")
    print(f"  Beta: {beta}, R²={'%.3f' % r2_achieved if r2_achieved else 'N/A'}, "
          f"Sig module-trait: {n_sig}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Per-tissue WGCNA (numpy/scipy)")
    parser.add_argument("--tissue", required=True,
                        help="Tissue name (e.g., liver, gastrocnemius)")
    args = parser.parse_args()
    tissue = args.tissue

    print(f"\n{'='*50}")
    print(f"WGCNA analysis: {tissue}")
    print(f"{'='*50}")

    # ── Load ──────────────────────────────────────────────────────────────────
    print(f"[{tissue}] Loading data...")
    datExpr, traits = load_data(tissue)

    is_small_n    = tissue in SMALL_N_TISSUES
    min_mod_size  = 20 if is_small_n else 30
    print(f"  minModuleSize={min_mod_size} (small_n={is_small_n})")

    # Standardize expression for correlation
    from sklearn.preprocessing import scale
    X = scale(datExpr.values, axis=0)  # samples × genes

    # ── Soft threshold ────────────────────────────────────────────────────────
    print(f"[{tissue}] Selecting soft threshold...")
    beta, r2_achieved, r2_met, sft_df, cor_mat = pick_soft_threshold(
        X, is_small_n=is_small_n)

    # ── Adjacency ─────────────────────────────────────────────────────────────
    print(f"[{tissue}] Computing adjacency (beta={beta})...")
    adj32 = compute_adjacency(cor_mat.astype(np.float32), beta)

    # ── TOM (vectorized) ──────────────────────────────────────────────────────
    print(f"[{tissue}] Computing TOM...")
    TOM = compute_tom(adj32)

    # Free memory
    del adj32, cor_mat
    import gc; gc.collect()

    # ── Module detection ──────────────────────────────────────────────────────
    print(f"[{tissue}] Detecting modules...")
    labels = detect_modules(TOM, min_mod_size, deep_split=DEEP_SPLIT)
    raw_colors = labels_to_colors(labels)
    n_raw = len(set(c for c in raw_colors if c != "grey"))
    print(f"  Before merging: {n_raw} modules, {sum(raw_colors == 'grey')} grey genes")

    # Free TOM (large matrix)
    del TOM
    gc.collect()

    # ── Module merging ────────────────────────────────────────────────────────
    print(f"[{tissue}] Merging similar modules (cutHeight={MERGE_CUT_HEIGHT})...")
    merged_colors = merge_close_modules(datExpr, raw_colors, MERGE_CUT_HEIGHT)
    n_final = len(set(c for c in merged_colors if c != "grey"))
    print(f"  After merging: {n_final} modules")

    # ── Module eigengenes ─────────────────────────────────────────────────────
    print(f"[{tissue}] Computing module eigengenes...")
    MEs = compute_eigengenes(datExpr, merged_colors)
    me_cols = [c for c in MEs.columns if c != "MEgrey"]
    MEs = MEs[me_cols] if me_cols else MEs
    print(f"  Modules: {', '.join(c.replace('ME','') for c in me_cols)}")

    # ── Module-trait correlations ─────────────────────────────────────────────
    print(f"[{tissue}] Computing module-trait correlations...")
    trait_r, trait_p = compute_trait_cors(MEs, traits)

    # ── kME (vectorized) ──────────────────────────────────────────────────────
    print(f"[{tissue}] Computing kME ({datExpr.shape[1]} genes × {len(me_cols)} modules)...")
    kme_df = compute_kme(datExpr, MEs) if len(me_cols) > 0 else None

    # ── Save ──────────────────────────────────────────────────────────────────
    print(f"[{tissue}] Saving results...")
    save_results(tissue, datExpr, merged_colors, MEs, trait_r, trait_p,
                 kme_df, beta, r2_achieved, r2_met, min_mod_size,
                 sft_df, n_raw)

    print(f"\n[{tissue}] WGCNA complete. Outputs: {OUTPUTS_DIR / tissue}")


if __name__ == "__main__":
    main()

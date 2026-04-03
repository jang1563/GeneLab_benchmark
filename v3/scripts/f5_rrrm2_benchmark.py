#!/usr/bin/env python3
"""F5 RRRM-2 scRNA-seq benchmark: composition, fGSEA, LOAO classifier.

Tasks:
  F5a — Cell-type composition shifts (FLT vs GC, Wilcoxon)
  F5b — Pseudo-bulk fGSEA per cell type (Hallmark 50)
  F5c — Leave-One-Animal-Out classifier per cell type (PCA-LR)
  F5d — Cross-mission PBMC concordance (RRRM-1 vs RRRM-2 NES)
  F5e — Bone marrow cell-type benchmark (novel tissue analysis)

Usage:
  python f5_rrrm2_benchmark.py [--tasks F5a,F5b,F5c] [--data-dir DIR]
"""
import argparse
import gc
import json
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Paths ──────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent  # v3/
EVAL_DIR = BASE_DIR / "evaluation"
FIG_DIR = BASE_DIR / "figures"
EVAL_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Default HPC data directory
DEFAULT_DATA_DIR = Path(
    os.environ.get("SCRATCH_DIR", ".")) / "huggingface" / "benchmark" / (
    "GeneLab_benchmark/v3/data/rrrm2"
)

# RRRM-2 tissue → GLDS mapping
TISSUE_GLDS = {
    "femur_bm": "402",
    "humerus_bm": "403",
    "pbmc": "404",
    "spleen": "405",
}

# RRRM-1 F2-B NES results for cross-mission concordance (F5d)
RRRM1_EVAL_DIR = Path(__file__).resolve().parent.parent.parent / "v2" / "evaluation"

# Thresholds
MIN_CELLS_PER_CONDITION = 20
MIN_ANIMALS_PER_GROUP = 2
MAX_CELLS_PER_CELLTYPE = 50000  # Subsample to prevent OOM on dense conversion
N_BOOTSTRAP = 2000
N_PERMUTATION = 5000
VARIANCE_FILTER_PERCENTILE = 0.25  # keep top 75%
MIN_PCA = 10
MAX_PCA = 50
RANDOM_SEED = 42


class NumpyEncoder(json.JSONEncoder):
    """JSON encoder that handles numpy types."""
    def default(self, obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def save_json(data: dict, path: Path):
    """Save dict to JSON with numpy-safe encoding."""
    with open(path, "w") as f:
        json.dump(data, f, indent=2, cls=NumpyEncoder)
    print(f"Saved: {path}")


# ── Data Loading ───────────────────────────────────────────────────────
def load_h5ad(data_dir: Path, glds: str):
    """Load h5ad for a GLDS dataset."""
    import anndata as ad

    path = data_dir / f"GLDS_{glds}.h5ad"
    if not path.exists():
        raise FileNotFoundError(f"h5ad not found: {path}")
    print(f"  Loading {path.name}...")
    adata = ad.read_h5ad(path)
    print(f"  Shape: {adata.shape}")
    return adata


def harmonize_metadata(adata, tissue: str):
    """Map RRRM-2 metadata columns to standard names.

    RRRM-2 uses:
      - exp: Flight / Ground  →  condition: FLT / GC
      - predicted.id  →  broad_celltype
      - orig.ident  →  animal_id  (e.g. FO1, GY4)
      - age: Old / Young  (kept as-is for F5e analysis)
    """
    obs = adata.obs.copy()

    # Condition mapping (strip whitespace, handle Categorical)
    if "exp" in obs.columns:
        obs["exp"] = obs["exp"].astype(str).str.strip()
        cond_map = {"Flight": "FLT", "Ground": "GC"}
        obs["condition"] = obs["exp"].map(cond_map).astype(str)
        unmapped = obs["condition"] == "nan"
        if unmapped.sum() > 0:
            bad_vals = obs.loc[unmapped, "exp"].unique()
            print(f"  WARNING: {unmapped.sum()} cells with unmapped exp values: {bad_vals}")
            obs.loc[unmapped, "condition"] = "UNKNOWN"
    else:
        raise KeyError(f"Missing 'exp' column in {tissue}")

    # Cell type
    ct_col_found = False
    for ct_col in ["predicted.id", "celltype", "cell_type", "broad_celltype"]:
        if ct_col in obs.columns:
            obs["broad_celltype"] = obs[ct_col].astype(str)
            ct_col_found = True
            break
    if not ct_col_found:
        print(f"  WARNING: No cell type column found for {tissue}")
        obs["broad_celltype"] = "unknown"
        obs["_no_celltype"] = True  # Flag for downstream tasks

    # Animal ID
    if "orig.ident" in obs.columns:
        obs["animal_id"] = obs["orig.ident"].astype(str)
    else:
        raise KeyError(f"Missing 'orig.ident' column in {tissue}")

    # Age (keep original)
    if "age" in obs.columns:
        obs["age_group"] = obs["age"].astype(str)

    adata.obs = obs
    return adata


# ── F5a: Composition ──────────────────────────────────────────────────
def compute_proportions(adata) -> pd.DataFrame:
    """Compute per-animal cell-type proportions."""
    obs = adata.obs
    animals = obs["animal_id"].unique()
    cell_types = sorted(obs["broad_celltype"].unique())

    rows = []
    for animal in animals:
        mask = obs["animal_id"] == animal
        sub = obs[mask]
        total = len(sub)
        cond = sub["condition"].iloc[0]
        row = {"animal_id": animal, "condition": cond, "n_cells": total}
        for ct in cell_types:
            row[ct] = (sub["broad_celltype"] == ct).sum() / total
        rows.append(row)

    return pd.DataFrame(rows).set_index("animal_id")


def run_f5a(data_dir: Path) -> dict:
    """F5a: Cell-type composition shifts (FLT vs GC)."""
    print("\n=== F5a: Composition Analysis ===")
    all_results = {}
    all_pvals = []  # For BH correction across all tissue × cell_type pairs

    for tissue, glds in TISSUE_GLDS.items():
        print(f"\n--- {tissue} (GLDS-{glds}) ---")
        try:
            adata = load_h5ad(data_dir, glds)
            adata = harmonize_metadata(adata, tissue)
        except (FileNotFoundError, KeyError) as e:
            print(f"  SKIP: {e}")
            continue

        prop_df = compute_proportions(adata)
        cell_types = [c for c in prop_df.columns
                      if c not in ("condition", "n_cells")]

        flt_animals = prop_df[prop_df["condition"] == "FLT"]
        gc_animals = prop_df[prop_df["condition"] == "GC"]
        print(f"  Animals: {len(flt_animals)} FLT, {len(gc_animals)} GC")

        tissue_results = {}
        for ct in cell_types:
            flt_vals = flt_animals[ct].values
            gc_vals = gc_animals[ct].values

            if len(flt_vals) < MIN_ANIMALS_PER_GROUP or len(gc_vals) < MIN_ANIMALS_PER_GROUP:
                continue

            try:
                stat, p = stats.mannwhitneyu(flt_vals, gc_vals, alternative="two-sided")
            except ValueError:
                p = 1.0

            result = {
                "flt_mean": float(np.mean(flt_vals)),
                "flt_std": float(np.std(flt_vals, ddof=1)) if len(flt_vals) > 1 else 0.0,
                "gc_mean": float(np.mean(gc_vals)),
                "gc_std": float(np.std(gc_vals, ddof=1)) if len(gc_vals) > 1 else 0.0,
                "delta": float(np.mean(flt_vals) - np.mean(gc_vals)),
                "p_raw": float(p),
                "n_flt": int(len(flt_vals)),
                "n_gc": int(len(gc_vals)),
            }
            tissue_results[ct] = result
            all_pvals.append((tissue, ct, p))

        all_results[tissue] = tissue_results
        del adata; gc.collect()

    # BH correction across all pairs
    if all_pvals:
        tissues_ct = [(t, ct) for t, ct, _ in all_pvals]
        pvals = [p for _, _, p in all_pvals]
        _, padj_arr, _, _ = multipletests(pvals, method="fdr_bh")
        for (t, ct), padj in zip(tissues_ct, padj_arr):
            all_results[t][ct]["padj"] = float(padj)

    return {"task": "F5-A", "description": "RRRM-2 cell-type composition", "results": all_results}


# ── F5b: Pseudo-bulk fGSEA ────────────────────────────────────────────
def aggregate_pseudobulk(adata, animal_col="animal_id") -> tuple:
    """Aggregate raw counts per animal → pseudo-bulk DataFrame.

    Uses .raw.X (raw counts) if available, else .X.
    Sparse-efficient: sums per-animal without full dense conversion.
    Returns: (bulk_df [animals × genes], conditions Series)
    """
    if adata.raw is not None:
        X = adata.raw.X
        genes = adata.raw.var_names.tolist()
    else:
        X = adata.X
        genes = adata.var_names.tolist()

    is_sparse = sp.issparse(X)
    obs = adata.obs
    animals = obs[animal_col].unique()

    bulk_data = []
    conditions = {}
    for animal in animals:
        mask = (obs[animal_col] == animal).values
        indices = np.where(mask)[0]
        if is_sparse:
            # Sparse row slicing + sum → stays sparse until final 1D array
            summed = np.asarray(X[indices].sum(axis=0)).flatten().astype(np.float32)
        else:
            summed = X[indices].sum(axis=0).astype(np.float32)
            if isinstance(summed, np.matrix):
                summed = np.asarray(summed).flatten()
        bulk_data.append(summed)
        conditions[animal] = obs.loc[mask, "condition"].iloc[0]

    bulk_df = pd.DataFrame(bulk_data, index=animals, columns=genes)
    cond_series = pd.Series(conditions)
    return bulk_df, cond_series


def rank_genes_ttest(bulk_df: pd.DataFrame, cond_series: pd.Series) -> pd.Series:
    """Welch t-test ranking (FLT vs GC). Returns t-statistics per gene (vectorized)."""
    flt_idx = cond_series[cond_series == "FLT"].index
    gc_idx = cond_series[cond_series == "GC"].index

    flt_arr = bulk_df.loc[flt_idx].values.astype(float)
    gc_arr = bulk_df.loc[gc_idx].values.astype(float)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        t_stats, _ = stats.ttest_ind(flt_arr, gc_arr, equal_var=False, axis=0)

    t_stats = np.where(np.isfinite(t_stats), t_stats, 0.0)
    return pd.Series(t_stats, index=bulk_df.columns)


def convert_ensmusg_to_symbol(ranking: pd.Series, var_source) -> pd.Series:
    """Convert ENSMUSG IDs to gene symbols for fGSEA compatibility.

    var_source: object with .var and .var_names (adata or adata.raw).
    """
    if len(ranking) == 0:
        return ranking
    if not str(ranking.index[0]).startswith("ENSMUSG"):
        return ranking

    # Check if var has gene_symbols or gene_name column
    sym_col = None
    for col in ["gene_symbols", "gene_name", "gene_short_name"]:
        if col in var_source.var.columns:
            sym_col = col
            break

    if sym_col is None:
        print("  WARNING: No gene symbol column found, using raw IDs")
        return ranking

    sym_map = dict(zip(var_source.var_names.astype(str),
                       var_source.var[sym_col].astype(str)))
    ranking.index = ranking.index.map(lambda x: sym_map.get(x, x))
    ranking = ranking[~ranking.index.duplicated(keep="first")]
    # Drop unmapped (still ENSMUSG)
    ranking = ranking[~ranking.index.str.startswith("ENSMUSG")]
    return ranking


def run_preranked_fgsea(ranking: pd.Series) -> pd.DataFrame:
    """Run gseapy preranked fGSEA with MSigDB Hallmark."""
    import gseapy as gp

    # Convert to uppercase for gseapy MSigDB matching
    ranking_upper = ranking.copy()
    ranking_upper.index = ranking_upper.index.str.upper()
    ranking_upper = ranking_upper[~ranking_upper.index.duplicated(keep="first")]
    ranking_upper = ranking_upper.dropna()
    ranking_upper = ranking_upper[ranking_upper != 0]

    if len(ranking_upper) < 100:
        print(f"  WARNING: Only {len(ranking_upper)} ranked genes, skipping fGSEA")
        return pd.DataFrame()

    try:
        res = gp.prerank(
            rnk=ranking_upper,
            gene_sets="MSigDB_Hallmark_2020",
            organism="Human",
            outdir=None,
            min_size=15,
            max_size=500,
            permutation_num=1000,
            seed=RANDOM_SEED,
            verbose=False,
        )
        result_df = res.res2d
    except Exception as e:
        print(f"  fGSEA error: {e}")
        return pd.DataFrame()

    # Normalize column names
    col_map = {}
    for col in result_df.columns:
        lc = col.lower()
        if lc in ("nes",):
            col_map[col] = "NES"
        elif "pval" in lc or "nom p" in lc:
            col_map[col] = "pval"
        elif "fdr" in lc:
            col_map[col] = "fdr"
        elif lc in ("term", "name"):
            col_map[col] = "pathway"
    result_df = result_df.rename(columns=col_map)

    return result_df


def run_f5b(data_dir: Path) -> dict:
    """F5b: Pseudo-bulk fGSEA per cell type."""
    print("\n=== F5b: Pseudo-bulk fGSEA ===")
    all_results = {}

    for tissue, glds in TISSUE_GLDS.items():
        print(f"\n--- {tissue} (GLDS-{glds}) ---")
        try:
            adata = load_h5ad(data_dir, glds)
            adata = harmonize_metadata(adata, tissue)
        except (FileNotFoundError, KeyError) as e:
            print(f"  SKIP: {e}")
            continue

        # Skip tissues with no cell type annotation entirely
        has_celltype = "_no_celltype" not in adata.obs.columns
        if not has_celltype:
            print(f"  No cell type annotation — skipping (run annotation first)")
            del adata; gc.collect()
            continue

        cell_types = sorted(adata.obs["broad_celltype"].unique())
        tissue_results = {}

        for ct in cell_types:
            ct_safe = ct.replace(" ", "_").replace("/", "_")
            mask = adata.obs["broad_celltype"] == ct
            sub = adata[mask].copy()

            n_flt = (sub.obs["condition"] == "FLT").sum()
            n_gc = (sub.obs["condition"] == "GC").sum()

            if n_flt < MIN_CELLS_PER_CONDITION or n_gc < MIN_CELLS_PER_CONDITION:
                print(f"  {ct}: SKIP (FLT={n_flt}, GC={n_gc})")
                del sub; continue

            # Count animals per condition
            flt_animals = sub.obs[sub.obs["condition"] == "FLT"]["animal_id"].nunique()
            gc_animals = sub.obs[sub.obs["condition"] == "GC"]["animal_id"].nunique()

            if flt_animals < MIN_ANIMALS_PER_GROUP or gc_animals < MIN_ANIMALS_PER_GROUP:
                print(f"  {ct}: SKIP (FLT animals={flt_animals}, GC animals={gc_animals})")
                del sub; continue

            print(f"  {ct}: {n_flt} FLT, {n_gc} GC ({flt_animals}+{gc_animals} animals)")

            # Pseudo-bulk aggregation (uses .raw if available)
            bulk_df, cond_series = aggregate_pseudobulk(sub)
            del sub; gc.collect()
            ranking = rank_genes_ttest(bulk_df, cond_series)
            del bulk_df; gc.collect()
            # Symbol mapping: use .raw.var if available, else .var
            sym_source = adata.raw if adata.raw is not None else adata
            ranking = convert_ensmusg_to_symbol(ranking, sym_source)

            # fGSEA
            fgsea_df = run_preranked_fgsea(ranking)
            del ranking
            if fgsea_df.empty:
                tissue_results[ct] = {
                    "n_cells": int(n_flt + n_gc),
                    "n_animals_flt": int(flt_animals),
                    "n_animals_gc": int(gc_animals),
                    "n_pathways": 0,
                    "status": "fgsea_failed",
                }
                continue

            # Save per-cell-type CSV
            csv_dir = EVAL_DIR / "F5B" / tissue
            csv_dir.mkdir(parents=True, exist_ok=True)
            fgsea_df.to_csv(csv_dir / f"{ct_safe}_fgsea_hallmark.csv", index=False)

            # Extract NES for top pathways
            nes_values = pd.Series(dtype=float)
            if "NES" in fgsea_df.columns:
                # Find pathway name column
                pw_col = "pathway" if "pathway" in fgsea_df.columns else None
                if pw_col is None:
                    for c in fgsea_df.columns:
                        if c.lower() in ("term", "name", "gene_set"):
                            pw_col = c
                            break
                if pw_col is not None:
                    nes_values = fgsea_df.set_index(pw_col)["NES"].astype(float)

            tissue_results[ct] = {
                "n_cells": int(n_flt + n_gc),
                "n_animals_flt": int(flt_animals),
                "n_animals_gc": int(gc_animals),
                "n_pathways": len(fgsea_df),
                "top_enriched": nes_values.nlargest(3).index.tolist() if len(nes_values) > 0 else [],
                "top_depleted": nes_values.nsmallest(3).index.tolist() if len(nes_values) > 0 else [],
            }

        all_results[tissue] = tissue_results
        del adata; gc.collect()

    return {"task": "F5-B", "description": "RRRM-2 pseudo-bulk fGSEA", "results": all_results}


# ── F5c: LOAO Classifier ──────────────────────────────────────────────
def pca_lr_predict(X_train, y_train, X_test):
    """Variance filter → StandardScaler → PCA → LogisticRegression."""
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    # Variance filter on raw data (top 75%) — before scaling
    variances = np.var(X_train, axis=0)
    threshold = np.percentile(variances, VARIANCE_FILTER_PERCENTILE * 100)
    keep = variances > threshold
    if keep.sum() < MIN_PCA:
        keep = np.ones(X_train.shape[1], dtype=bool)
    X_train = X_train[:, keep]
    X_test = X_test[:, keep]

    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)

    # PCA
    n_components = min(MAX_PCA, X_train_s.shape[0] - 1, X_train_s.shape[1])
    n_components = max(n_components, MIN_PCA)
    n_components = min(n_components, X_train_s.shape[1], X_train_s.shape[0] - 1)

    if n_components < 1:
        return np.full(X_test.shape[0], 0.5)

    pca = PCA(n_components=n_components, random_state=RANDOM_SEED)
    X_train_p = pca.fit_transform(X_train_s)
    X_test_p = pca.transform(X_test_s)

    # Logistic Regression
    lr = LogisticRegression(
        class_weight="balanced",
        max_iter=2000,
        random_state=RANDOM_SEED,
    )
    lr.fit(X_train_p, y_train)
    return lr.predict_proba(X_test_p)[:, 1]


def loao_cv(adata_sub) -> dict:
    """Leave-One-Animal-Out cross-validation."""
    from sklearn.metrics import roc_auc_score

    X = adata_sub.X
    if sp.issparse(X):
        X = X.toarray()
    X = X.astype(np.float32)

    obs = adata_sub.obs
    y = (obs["condition"] == "FLT").astype(int).values
    animals = obs["animal_id"].values
    unique_animals = np.unique(animals)

    all_y_true = []
    all_y_score = []
    fold_details = []

    for animal_out in unique_animals:
        test_mask = animals == animal_out
        train_mask = ~test_mask

        X_train, y_train = X[train_mask], y[train_mask]
        X_test, y_test = X[test_mask], y[test_mask]

        if len(np.unique(y_train)) < 2:
            continue

        y_score = pca_lr_predict(X_train, y_train, X_test)
        all_y_true.extend(y_test.tolist())
        all_y_score.extend(y_score.tolist())

        cond_out = obs.loc[test_mask, "condition"].iloc[0]
        fold_auroc = float("nan")
        if len(np.unique(y_test)) == 2:
            fold_auroc = float(roc_auc_score(y_test, y_score))

        fold_details.append({
            "animal_out": str(animal_out),
            "condition_out": cond_out,
            "n_test": int(test_mask.sum()),
            "n_train": int(train_mask.sum()),
            "fold_auroc": fold_auroc,
        })

    all_y_true = np.array(all_y_true)
    all_y_score = np.array(all_y_score)

    if len(np.unique(all_y_true)) < 2 or len(all_y_true) == 0:
        return {"auroc": float("nan"), "status": "insufficient_data"}

    auroc = float(roc_auc_score(all_y_true, all_y_score))

    # Bootstrap CI
    rng = np.random.RandomState(RANDOM_SEED)
    boot_aurocs = []
    for _ in range(N_BOOTSTRAP):
        idx = rng.choice(len(all_y_true), size=len(all_y_true), replace=True)
        if len(np.unique(all_y_true[idx])) < 2:
            continue
        boot_aurocs.append(roc_auc_score(all_y_true[idx], all_y_score[idx]))
    ci_low = float(np.percentile(boot_aurocs, 2.5)) if boot_aurocs else float("nan")
    ci_high = float(np.percentile(boot_aurocs, 97.5)) if boot_aurocs else float("nan")

    # Permutation test
    perm_count = 0
    for _ in range(N_PERMUTATION):
        perm_y = rng.permutation(all_y_true)
        if len(np.unique(perm_y)) < 2:
            continue
        perm_auroc = roc_auc_score(perm_y, all_y_score)
        if perm_auroc >= auroc:
            perm_count += 1
    p_perm = (perm_count + 1) / (N_PERMUTATION + 1)

    return {
        "auroc": auroc,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "p_perm": float(p_perm),
        "n_animals": len(unique_animals),
        "n_cells_flt": int((obs["condition"] == "FLT").sum()),
        "n_cells_gc": int((obs["condition"] == "GC").sum()),
        "n_cells_total": len(obs),
        "fold_details": fold_details,
    }


def run_f5c(data_dir: Path) -> dict:
    """F5c: Leave-One-Animal-Out classifier per cell type."""
    print("\n=== F5c: LOAO Classifier ===")
    all_results = {}

    for tissue, glds in TISSUE_GLDS.items():
        print(f"\n--- {tissue} (GLDS-{glds}) ---")
        try:
            adata = load_h5ad(data_dir, glds)
            adata = harmonize_metadata(adata, tissue)
        except (FileNotFoundError, KeyError) as e:
            print(f"  SKIP: {e}")
            continue

        # Skip tissues with no cell type annotation entirely
        has_celltype = "_no_celltype" not in adata.obs.columns
        if not has_celltype:
            print(f"  No cell type annotation — skipping (run annotation first)")
            del adata; gc.collect()
            continue

        cell_types = sorted(adata.obs["broad_celltype"].unique())
        tissue_results = {}

        for ct in cell_types:
            mask = adata.obs["broad_celltype"] == ct
            sub = adata[mask].copy()

            n_flt = (sub.obs["condition"] == "FLT").sum()
            n_gc = (sub.obs["condition"] == "GC").sum()

            if n_flt < MIN_CELLS_PER_CONDITION or n_gc < MIN_CELLS_PER_CONDITION:
                print(f"  {ct}: SKIP (FLT={n_flt}, GC={n_gc})")
                del sub; continue

            # Check annotation bias (>90% cells from one condition)
            flt_frac = n_flt / (n_flt + n_gc)
            if flt_frac > 0.9 or flt_frac < 0.1:
                print(f"  {ct}: FLAG annotation bias (FLT={flt_frac:.1%})")
                tissue_results[ct] = {
                    "status": "annotation_bias",
                    "flt_fraction": float(flt_frac),
                    "n_flt": int(n_flt),
                    "n_gc": int(n_gc),
                }
                del sub; continue

            # Subsample if too many cells (prevent OOM on dense conversion)
            total_cells = n_flt + n_gc
            if total_cells > MAX_CELLS_PER_CELLTYPE:
                rng = np.random.default_rng(RANDOM_SEED)
                idx = rng.choice(total_cells, size=MAX_CELLS_PER_CELLTYPE, replace=False)
                sub = sub[idx].copy()
                print(f"  {ct}: subsampled {total_cells} → {MAX_CELLS_PER_CELLTYPE} cells")

            print(f"  {ct}: {n_flt} FLT + {n_gc} GC cells, LOAO...")
            result = loao_cv(sub)
            del sub; gc.collect()
            result["n_cells_original"] = int(total_cells)
            tissue_results[ct] = result

            if "auroc" in result and not np.isnan(result["auroc"]):
                print(f"    AUROC={result['auroc']:.3f} [{result.get('ci_low', 0):.3f}, {result.get('ci_high', 0):.3f}] p={result.get('p_perm', 1):.4f}")

        all_results[tissue] = tissue_results
        del adata; gc.collect()

    return {
        "task": "F5-C",
        "method": "LOAO PCA(50)-LR(balanced), bootstrap CI n=2000, permutation p n=5000",
        "description": "RRRM-2 cell-type LOAO classifier",
        "results": all_results,
    }


# ── F5d: Cross-Mission PBMC Concordance ───────────────────────────────
def load_rrrm1_nes(tissue: str = "blood") -> dict:
    """Load RRRM-1 (v2) F2-B NES results for a tissue."""
    nes_dir = RRRM1_EVAL_DIR / "F2B" / tissue
    if not nes_dir.exists():
        print(f"  RRRM-1 F2-B NES dir not found: {nes_dir}")
        return {}

    nes_by_ct = {}
    for csv_path in sorted(nes_dir.glob("*_fgsea_hallmark.csv")):
        ct = csv_path.stem.replace("_fgsea_hallmark", "")
        df = pd.read_csv(csv_path)

        # Find NES column
        nes_col = None
        for c in df.columns:
            if c.upper() == "NES":
                nes_col = c
                break
        if nes_col is None:
            continue

        # Find pathway column
        pw_col = df.columns[0]
        for c in df.columns:
            if "term" in c.lower() or "pathway" in c.lower() or "name" in c.lower():
                pw_col = c
                break

        nes_by_ct[ct] = dict(zip(df[pw_col], df[nes_col].astype(float)))

    return nes_by_ct


def run_f5d(data_dir: Path) -> dict:
    """F5d: Cross-mission PBMC concordance (RRRM-1 blood vs RRRM-2 PBMC NES).

    Requires F5b to have been run first (reads F5B/pbmc/*.csv).
    """
    print("\n=== F5d: Cross-Mission PBMC Concordance ===")

    # Load RRRM-1 blood NES
    rrrm1_nes = load_rrrm1_nes("blood")
    if not rrrm1_nes:
        print("  No RRRM-1 F2-B blood NES found, skipping F5d")
        return {"task": "F5-D", "status": "no_rrrm1_data"}

    # Load RRRM-2 PBMC NES from F5b CSV outputs
    rrrm2_nes_dir = EVAL_DIR / "F5B" / "pbmc"
    rrrm2_nes = {}
    if rrrm2_nes_dir.exists():
        for csv_path in sorted(rrrm2_nes_dir.glob("*_fgsea_hallmark.csv")):
            ct = csv_path.stem.replace("_fgsea_hallmark", "")
            df = pd.read_csv(csv_path)
            nes_col = None
            for c in df.columns:
                if c.upper() == "NES":
                    nes_col = c
                    break
            if nes_col is None:
                continue
            pw_col = df.columns[0]
            for c in df.columns:
                if "term" in c.lower() or "pathway" in c.lower() or "name" in c.lower():
                    pw_col = c
                    break
            rrrm2_nes[ct] = dict(zip(df[pw_col], df[nes_col].astype(float)))

    if not rrrm2_nes:
        print("  No RRRM-2 PBMC NES found. Run F5b first.")
        return {"task": "F5-D", "status": "no_rrrm2_nes"}

    # Compute pairwise Spearman correlations
    concordance = []
    for ct1 in sorted(rrrm1_nes.keys()):
        for ct2 in sorted(rrrm2_nes.keys()):
            nes1 = rrrm1_nes[ct1]
            nes2 = rrrm2_nes[ct2]
            common = sorted(set(nes1.keys()) & set(nes2.keys()))

            if len(common) < 5:
                continue

            x = np.array([nes1[pw] for pw in common])
            y = np.array([nes2[pw] for pw in common])
            r, p = stats.spearmanr(x, y)

            concordance.append({
                "rrrm1_celltype": ct1,
                "rrrm2_celltype": ct2,
                "n_common_pathways": len(common),
                "spearman_r": float(r),
                "spearman_p": float(p),
            })

    return {
        "task": "F5-D",
        "description": "RRRM-1 blood vs RRRM-2 PBMC NES concordance",
        "concordance": concordance,
    }


# ── F5e: Bone Marrow Benchmark ────────────────────────────────────────
def run_f5e(data_dir: Path) -> dict:
    """F5e: Bone marrow cell-type benchmark (femur + humerus, novel tissues)."""
    print("\n=== F5e: Bone Marrow Benchmark ===")
    bm_results = {}

    for tissue in ["femur_bm", "humerus_bm"]:
        glds = TISSUE_GLDS[tissue]
        print(f"\n--- {tissue} (GLDS-{glds}) ---")
        try:
            adata = load_h5ad(data_dir, glds)
            adata = harmonize_metadata(adata, tissue)
        except (FileNotFoundError, KeyError) as e:
            print(f"  SKIP: {e}")
            continue

        cell_types = sorted(adata.obs["broad_celltype"].unique())
        n_cells = len(adata)
        n_flt = (adata.obs["condition"] == "FLT").sum()
        n_gc = (adata.obs["condition"] == "GC").sum()

        tissue_summary = {
            "n_cells_total": int(n_cells),
            "n_flt": int(n_flt),
            "n_gc": int(n_gc),
            "n_cell_types": len(cell_types),
            "cell_types": cell_types,
        }

        # LOAO per cell type (same as F5c but focused on BM)
        ct_results = {}
        for ct in cell_types:
            mask = adata.obs["broad_celltype"] == ct
            sub = adata[mask].copy()

            ct_n_flt = (sub.obs["condition"] == "FLT").sum()
            ct_n_gc = (sub.obs["condition"] == "GC").sum()

            if ct_n_flt < MIN_CELLS_PER_CONDITION or ct_n_gc < MIN_CELLS_PER_CONDITION:
                ct_results[ct] = {"status": "insufficient_cells",
                                  "n_flt": int(ct_n_flt), "n_gc": int(ct_n_gc)}
                continue

            flt_frac = ct_n_flt / (ct_n_flt + ct_n_gc)
            if flt_frac > 0.9 or flt_frac < 0.1:
                ct_results[ct] = {"status": "annotation_bias",
                                  "flt_fraction": float(flt_frac)}
                continue

            print(f"  {ct}: {ct_n_flt} FLT + {ct_n_gc} GC, LOAO...")
            result = loao_cv(sub)
            ct_results[ct] = result

            if "auroc" in result and not np.isnan(result["auroc"]):
                print(f"    AUROC={result['auroc']:.3f}")

        tissue_summary["cell_type_results"] = ct_results
        bm_results[tissue] = tissue_summary
        del adata; gc.collect()

    # Cross-site concordance (femur vs humerus)
    if "femur_bm" in bm_results and "humerus_bm" in bm_results:
        femur_ct = set(bm_results["femur_bm"].get("cell_type_results", {}).keys())
        humerus_ct = set(bm_results["humerus_bm"].get("cell_type_results", {}).keys())
        common_ct = femur_ct & humerus_ct

        cross_site = {}
        for ct in sorted(common_ct):
            f_res = bm_results["femur_bm"]["cell_type_results"][ct]
            h_res = bm_results["humerus_bm"]["cell_type_results"][ct]
            if "auroc" in f_res and "auroc" in h_res:
                cross_site[ct] = {
                    "femur_auroc": f_res["auroc"],
                    "humerus_auroc": h_res["auroc"],
                    "delta": f_res["auroc"] - h_res["auroc"],
                }
        bm_results["cross_site_concordance"] = cross_site

    return {
        "task": "F5-E",
        "description": "RRRM-2 bone marrow cell-type benchmark",
        "results": bm_results,
    }


# ── Main ───────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="F5 RRRM-2 benchmark")
    parser.add_argument("--tasks", default="F5a,F5b,F5c,F5d,F5e",
                        help="Comma-separated task list")
    parser.add_argument("--data-dir", type=str, default=str(DEFAULT_DATA_DIR),
                        help="Directory with GLDS_*.h5ad files")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    tasks = [t.strip().upper() for t in args.tasks.split(",")]
    print(f"Data dir: {data_dir}")
    print(f"Tasks: {tasks}")

    results = {}

    if "F5A" in tasks:
        r = run_f5a(data_dir)
        results["F5a"] = r
        save_json(r, EVAL_DIR / "F5A_rrrm2_composition.json")

    if "F5B" in tasks:
        r = run_f5b(data_dir)
        results["F5b"] = r
        save_json(r, EVAL_DIR / "F5B_rrrm2_fgsea.json")

    if "F5C" in tasks:
        r = run_f5c(data_dir)
        results["F5c"] = r
        save_json(r, EVAL_DIR / "F5C_rrrm2_loao.json")

    if "F5D" in tasks:
        r = run_f5d(data_dir)
        results["F5d"] = r
        save_json(r, EVAL_DIR / "F5D_rrrm2_cross_mission.json")

    if "F5E" in tasks:
        r = run_f5e(data_dir)
        results["F5e"] = r
        save_json(r, EVAL_DIR / "F5E_rrrm2_bone_marrow.json")

    print("\n=== F5 RRRM-2 Benchmark Complete ===")


if __name__ == "__main__":
    main()

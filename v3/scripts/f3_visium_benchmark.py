#!/usr/bin/env python3
"""
f3_visium_benchmark.py — OSD-352 RR-3 Brain Visium Benchmark (F3a-d)

Tasks:
  F3a: Section-level LOAO classification (FLT vs GC)
  F3b: Spatially variable genes (Moran's I) FLT vs GC comparison
  F3c: Brain region-specific AUROC (requires Seurat cluster metadata)
  F3d: Spatial pseudo-bulk vs companion bulk RNA-seq comparison

Data: Masarapu et al. 2024 Nature Communications
  OSD-352: 3 FLT + 3 GC animals, 12 Visium sections

Statistical notes:
  - n=6 animals (3+3) → C(6,3) = 20 exact permutations
  - Minimum possible p-value = 1/20 = 0.05
  - Effect sizes (Cohen's d) reported alongside p-values
"""

import sys
import json
import argparse
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats, sparse
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score

# ── Constants ──────────────────────────────────────────────────────────────
SEED = 42
N_BOOTSTRAP = 2000
MAX_PCA = 5           # Conservative due to small n
MIN_PCA = 2
VARIANCE_FILTER_PCTL = 0.50  # Keep top 50% by variance (conservative)

BASE_DIR = Path(__file__).resolve().parent.parent  # v3/
DATA_DIR = BASE_DIR / "data" / "spatial" / "osd352"
EVAL_DIR = BASE_DIR / "evaluation"
FIG_DIR = BASE_DIR / "figures"

for d in [EVAL_DIR, FIG_DIR]:
    d.mkdir(parents=True, exist_ok=True)


# ── Section → Animal/Condition Mapping ─────────────────────────────────────
# OSD-352: 3 slides (158, 159, 304), 4 capture areas each (A1, B1, C1, D1)
# From Masarapu et al. 2024 supplementary:
#   Slide 158: sections from 2 animals (likely mixed FLT/GC on same slide)
#   Slide 159: sections from 2 animals
#   Slide 304: sections from 2 animals
#
# Mapping will be loaded from seurat_metadata.csv or manually set.
# Each capture area = 1 section from 1 animal's brain.
# With 12 sections / 6 animals = 2 sections per animal.

def load_section_condition_mapping(data_dir: Path) -> dict:
    """Load section → (animal_id, condition) mapping.

    Returns dict: section_name → {"animal": str, "condition": "FLT"|"GC"}
    """
    mapping_path = data_dir / "section_mapping.json"
    if mapping_path.exists():
        with open(mapping_path) as f:
            return json.load(f)

    seurat_meta_path = data_dir / "seurat_metadata.csv"
    if seurat_meta_path.exists():
        meta = pd.read_csv(seurat_meta_path, index_col=0)
        # Derive mapping from Seurat metadata
        mapping = derive_mapping_from_seurat(meta)
        # Cache it
        with open(mapping_path, "w") as f:
            json.dump(mapping, f, indent=2)
        return mapping

    print("ERROR: No section_mapping.json or seurat_metadata.csv found!")
    print("Please create section_mapping.json with format:")
    print('  {"Sample_158_A1": {"animal": "F1", "condition": "FLT"}, ...}')
    sys.exit(1)


def derive_mapping_from_seurat(meta: pd.DataFrame) -> dict:
    """Derive section→animal→condition from Seurat metadata.

    Seurat barcodes typically contain section identifiers.
    """
    mapping = {}

    # Try to find sample/orig.ident column
    sample_col = None
    for col in ["orig.ident", "sample", "library_id"]:
        if col in meta.columns:
            sample_col = col
            break

    if sample_col is None:
        print("  WARNING: No sample column found in Seurat metadata")
        return mapping

    # Try to find condition column
    cond_col = None
    for col in meta.columns:
        vals = set(meta[col].astype(str).str.lower().unique())
        # Look for columns with FLT/GC or Flight/Ground values
        if vals & {"flt", "gc", "flight", "ground", "space.flight",
                    "ground.control", "spaceflight", "groundcontrol"}:
            cond_col = col
            break

    if cond_col is None:
        print("  WARNING: No condition column found in Seurat metadata")
        return mapping

    # Build mapping: for each section (orig.ident value),
    # what is the condition?
    for section_id in meta[sample_col].unique():
        section_meta = meta[meta[sample_col] == section_id]
        conditions = section_meta[cond_col].unique()
        if len(conditions) == 1:
            cond_raw = str(conditions[0]).strip()
            # Normalize condition
            cond_lower = cond_raw.lower()
            if cond_lower in ("flt", "flight", "space.flight", "spaceflight"):
                condition = "FLT"
            elif cond_lower in ("gc", "ground", "ground.control",
                                "groundcontrol"):
                condition = "GC"
            else:
                condition = cond_raw

            # Try to extract animal ID
            # Pattern: section_id might be like "Sample_158_A1"
            # Animal ID might be embedded or from another column
            animal_id = section_id  # Default: use section as animal proxy

            # Check for animal column
            for acol in ["animal", "animal_id", "mouse_id", "replicate"]:
                if acol in section_meta.columns:
                    animals = section_meta[acol].unique()
                    if len(animals) == 1:
                        animal_id = str(animals[0])
                    break

            mapping[section_id] = {
                "animal": animal_id,
                "condition": condition
            }

    print(f"  Derived mapping for {len(mapping)} sections:")
    for s, info in sorted(mapping.items()):
        print(f"    {s}: animal={info['animal']}, condition={info['condition']}")
    return mapping


# ── Utility Functions ──────────────────────────────────────────────────────

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def save_json(data: dict, path: Path):
    with open(path, "w") as f:
        json.dump(data, f, indent=2, cls=NumpyEncoder)
    print(f"  Saved: {path}")


def cohens_d(group1: np.ndarray, group2: np.ndarray) -> float:
    """Compute Cohen's d effect size."""
    n1, n2 = len(group1), len(group2)
    if n1 < 2 or n2 < 2:
        return float("nan")
    var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)
    pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))
    if pooled_std == 0:
        return float("nan")
    return (np.mean(group1) - np.mean(group2)) / pooled_std


# ── PCA-LR Classification ─────────────────────────────────────────────────

def pca_lr_predict(X_train, y_train, X_test):
    """Variance filter → Scale → PCA → LR prediction.

    Adapted from f5_rrrm2_benchmark.py for small n.
    """
    n_train, n_features = X_train.shape

    # Variance filter (keep top 50%)
    variances = np.var(X_train, axis=0)
    threshold = np.percentile(variances, VARIANCE_FILTER_PCTL * 100)
    mask = variances > threshold
    if mask.sum() < MIN_PCA:
        mask = np.ones(n_features, dtype=bool)  # Keep all
    X_train = X_train[:, mask]
    X_test = X_test[:, mask]
    n_features = X_train.shape[1]

    # Scale
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    # PCA
    n_components = min(MAX_PCA, n_train - 1, n_features)
    n_components = max(n_components, MIN_PCA)
    n_components = min(n_components, n_features, n_train - 1)
    if n_components < 1:
        return np.full(X_test.shape[0], 0.5)

    pca = PCA(n_components=n_components, random_state=SEED)
    X_train_pca = pca.fit_transform(X_train)
    X_test_pca = pca.transform(X_test)

    # Logistic Regression
    lr = LogisticRegression(class_weight="balanced", max_iter=2000,
                            random_state=SEED)
    try:
        lr.fit(X_train_pca, y_train)
        return lr.predict_proba(X_test_pca)[:, 1]
    except Exception:
        return np.full(X_test.shape[0], 0.5)


def _run_loao_cv(X, y, animals):
    """Run leave-one-animal-out CV (internal helper).

    Returns:
        (all_y_true, all_y_score) as np.ndarray, or (None, None) if failed
    """
    unique_animals = np.unique(animals)
    all_y_true = []
    all_y_score = []

    for animal_out in unique_animals:
        test_mask = animals == animal_out
        train_mask = ~test_mask
        y_train, y_test = y[train_mask], y[test_mask]
        X_train, X_test = X[train_mask], X[test_mask]

        if len(np.unique(y_train)) < 2:
            continue

        y_score = pca_lr_predict(X_train, y_train, X_test)
        all_y_true.extend(y_test)
        all_y_score.extend(y_score)

    all_y_true = np.array(all_y_true)
    all_y_score = np.array(all_y_score)

    if len(all_y_true) == 0 or len(np.unique(all_y_true)) < 2:
        return None, None

    return all_y_true, all_y_score


def loao_classification(X, y, animals, section_names=None):
    """Leave-One-Animal-Out CV with exact permutation test (full re-fit).

    Args:
        X: feature matrix (n_samples × n_features)
        y: labels (0=GC, 1=FLT)
        animals: animal IDs per sample
        section_names: optional section names for fold details

    Returns:
        dict with auroc, ci, p-values, fold details
    """
    rng = np.random.default_rng(SEED)
    unique_animals = np.unique(animals)
    n_animals = len(unique_animals)

    # Run observed LOAO
    all_y_true, all_y_score = _run_loao_cv(X, y, animals)
    if all_y_true is None:
        return {"auroc": float("nan"), "status": "insufficient_classes"}

    auroc = roc_auc_score(all_y_true, all_y_score)

    # Build fold details (observed labels)
    fold_details = []
    for animal_out in unique_animals:
        test_mask = animals == animal_out
        cond_out = "FLT" if y[test_mask][0] == 1 else "GC"
        detail = {
            "animal_out": str(animal_out),
            "condition_out": cond_out,
            "n_test": int(test_mask.sum()),
            "n_train": int((~test_mask).sum()),
        }
        if section_names is not None:
            detail["sections_out"] = [section_names[i] for i in
                                       np.where(test_mask)[0]]
        fold_details.append(detail)

    # Bootstrap CI
    boot_aurocs = []
    for _ in range(N_BOOTSTRAP):
        idx = rng.choice(len(all_y_true), size=len(all_y_true), replace=True)
        if len(np.unique(all_y_true[idx])) >= 2:
            boot_aurocs.append(roc_auc_score(all_y_true[idx],
                                              all_y_score[idx]))
    ci_low = np.percentile(boot_aurocs, 2.5) if boot_aurocs else float("nan")
    ci_high = np.percentile(boot_aurocs, 97.5) if boot_aurocs else float("nan")

    # Exact permutation test (full re-fit, animal-block)
    # C(n_animals, n_flt) unique assignments — C(6,3) = 20
    n_flt = sum(1 for a in unique_animals
                if y[animals == a][0] == 1)

    perm_count = 0
    n_perms = 0
    for combo in combinations(range(n_animals), n_flt):
        perm_y = np.zeros_like(y)
        for i, a in enumerate(unique_animals):
            mask = animals == a
            if i in combo:
                perm_y[mask] = 1
            else:
                perm_y[mask] = 0

        if len(np.unique(perm_y)) < 2:
            continue
        n_perms += 1

        # Full re-fit LOAO with permuted labels
        perm_true, perm_score = _run_loao_cv(X, perm_y, animals)
        if perm_true is None:
            continue
        perm_auroc = roc_auc_score(perm_true, perm_score)
        if perm_auroc >= auroc:
            perm_count += 1

    p_exact = perm_count / n_perms if n_perms > 0 else float("nan")

    # Cohen's d on PC1 scores
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    n_comp = min(MAX_PCA, X.shape[0] - 1, X.shape[1])
    if n_comp >= 1:
        pca = PCA(n_components=n_comp, random_state=SEED)
        X_pca = pca.fit_transform(X_scaled)
        pc1_flt = X_pca[y == 1, 0]
        pc1_gc = X_pca[y == 0, 0]
        effect_size = cohens_d(pc1_flt, pc1_gc)
    else:
        effect_size = float("nan")

    return {
        "auroc": auroc,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "p_exact": p_exact,
        "n_permutations": n_perms,
        "effect_size_cohens_d": effect_size,
        "n_animals": n_animals,
        "n_flt": n_flt,
        "n_gc": n_animals - n_flt,
        "n_samples": len(y),
        "fold_details": fold_details
    }


# ── F3a: Section-Level Classification ──────────────────────────────────────

def run_f3a(data_dir: Path, mapping: dict) -> dict:
    """F3a: Section-level and animal-level LOAO classification."""
    print("\n" + "=" * 60)
    print("F3a: Section-Level LOAO Classification")
    print("=" * 60)

    processed_dir = data_dir / "processed"
    X = np.load(processed_dir / "pseudobulk_matrix.npy")
    gene_names = np.load(processed_dir / "gene_names.npy", allow_pickle=True)
    with open(processed_dir / "section_names.json") as f:
        section_names = json.load(f)

    print(f"  Loaded: {X.shape} (sections × genes)")

    # Map sections to animals and conditions
    animals = []
    conditions = []
    valid_sections = []
    valid_indices = []

    for i, sec in enumerate(section_names):
        if sec in mapping:
            animals.append(mapping[sec]["animal"])
            conditions.append(mapping[sec]["condition"])
            valid_sections.append(sec)
            valid_indices.append(i)
        else:
            print(f"  WARNING: {sec} not in mapping, skipping")

    X = X[valid_indices]
    animals = np.array(animals)
    y = np.array([1 if c == "FLT" else 0 for c in conditions])
    section_names = valid_sections

    print(f"  FLT sections: {sum(y == 1)}, GC sections: {sum(y == 0)}")
    print(f"  Animals: {np.unique(animals)}")

    # Section-level LOAO (animal-block permutation)
    print("\n  [Section-level LOAO]")
    section_result = loao_classification(X, y, animals, section_names)
    print(f"  AUROC = {section_result['auroc']:.3f} "
          f"[{section_result['ci_low']:.3f}, {section_result['ci_high']:.3f}]")
    print(f"  Exact perm p = {section_result['p_exact']:.3f} "
          f"(n_perms={section_result['n_permutations']})")
    print(f"  Cohen's d = {section_result['effect_size_cohens_d']:.3f}")

    # Animal-level LOAO (average 2 sections per animal → n=6)
    print("\n  [Animal-level LOAO]")
    unique_animals = np.unique(animals)
    n_genes = X.shape[1]
    X_animal = np.zeros((len(unique_animals), n_genes), dtype=np.float32)
    y_animal = np.zeros(len(unique_animals), dtype=int)
    animal_names = []

    for i, a in enumerate(unique_animals):
        mask = animals == a
        X_animal[i] = X[mask].mean(axis=0)
        y_animal[i] = y[mask][0]  # Same condition for same animal
        animal_names.append(a)

    animal_result = loao_classification(X_animal, y_animal,
                                        np.array(animal_names))
    print(f"  AUROC = {animal_result['auroc']:.3f} "
          f"[{animal_result['ci_low']:.3f}, {animal_result['ci_high']:.3f}]")
    print(f"  Exact perm p = {animal_result['p_exact']:.3f}")

    result = {
        "task": "F3a",
        "description": "OSD-352 RR-3 Brain Visium section-level classification",
        "dataset": "OSD-352",
        "mission": "RR-3",
        "tissue": "brain",
        "method": f"PCA(≤{MAX_PCA})-LR(balanced), LOAO, "
                  f"exact permutation (C(6,3)=20)",
        "section_level": section_result,
        "animal_level": animal_result,
        "n_genes": int(n_genes),
        "section_mapping": mapping
    }

    with open(processed_dir / "section_info.json") as f:
        spot_info = json.load(f)
    result["spot_counts"] = spot_info

    return result


# ── F3b: Spatially Variable Genes ──────────────────────────────────────────

def precompute_spatial_neighbors(coords: np.ndarray, k: int = 6) -> np.ndarray:
    """Build k-NN index once for a section's spatial coordinates.

    Args:
        coords: spatial coordinates (n_spots, 2)
        k: number of nearest neighbors

    Returns:
        neighbor_indices: (n_spots, k) indices of k nearest neighbors
    """
    from scipy.spatial import cKDTree
    tree = cKDTree(coords)
    _, indices = tree.query(coords, k=k + 1)  # +1 for self
    return indices[:, 1:]  # Remove self


def compute_morans_i(values: np.ndarray,
                     neighbor_indices: np.ndarray) -> tuple:
    """Compute Moran's I using precomputed k-NN neighbors (vectorized).

    Args:
        values: gene expression values (n_spots,)
        neighbor_indices: k-NN indices (n_spots, k) from
                          precompute_spatial_neighbors()

    Returns:
        (morans_i, z_score, p_value)
    """
    n = len(values)
    k = neighbor_indices.shape[1]
    if n < 10:
        return float("nan"), float("nan"), float("nan")

    x = values - np.mean(values)
    ss = np.sum(x ** 2)
    if ss == 0:
        return 0.0, 0.0, 1.0

    # Vectorized Moran's I
    W = n * k  # total weight (binary k-NN)
    neighbor_x = x[neighbor_indices]  # (n, k)
    numerator = np.sum(x * neighbor_x.sum(axis=1))

    I = (n / W) * (numerator / ss)

    # Expected value and variance under null (normality assumption)
    E_I = -1.0 / (n - 1)
    var_I = 1.0 / (n - 1) - E_I ** 2
    if var_I <= 0:
        return I, 0.0, 1.0

    z = (I - E_I) / np.sqrt(var_I)
    p = 2 * (1 - stats.norm.cdf(abs(z)))  # Two-tailed

    return float(I), float(z), float(p)


def run_f3b(data_dir: Path, mapping: dict) -> dict:
    """F3b: Spatially variable genes FLT vs GC comparison."""
    print("\n" + "=" * 60)
    print("F3b: Spatially Variable Genes (Moran's I)")
    print("=" * 60)

    import anndata as ad
    processed_dir = data_dir / "processed"
    merged_path = processed_dir / "osd352_visium_merged.h5ad"

    if not merged_path.exists():
        print("  ERROR: merged h5ad not found. Run process_osd352_visium.py first.")
        return {"task": "F3b", "status": "data_not_found"}

    adata = ad.read_h5ad(merged_path)
    print(f"  Loaded merged: {adata.n_obs} spots × {adata.n_vars} genes")

    # Select top variable genes for Moran's I (computing for all genes is slow)
    # Use HVGs as a proxy
    sc_available = True
    try:
        import scanpy as sc
        adata_hvg = adata.copy()
        sc.pp.highly_variable_genes(adata_hvg, n_top_genes=2000)
        hvg_genes = adata_hvg.var_names[adata_hvg.var["highly_variable"]]
        print(f"  Selected {len(hvg_genes)} HVGs for Moran's I analysis")
    except Exception as e:
        print(f"  WARNING: HVG selection failed ({e}), using top 1000 by variance")
        X_dense = adata.X.toarray() if sparse.issparse(adata.X) else adata.X
        var = np.var(X_dense, axis=0)
        top_idx = np.argsort(var)[-1000:]
        hvg_genes = adata.var_names[top_idx]
        sc_available = False

    # Compute Moran's I per section per gene
    sections = adata.obs["section"].unique()
    flt_sections = [s for s in sections if s in mapping
                    and mapping[s]["condition"] == "FLT"]
    gc_sections = [s for s in sections if s in mapping
                   and mapping[s]["condition"] == "GC"]

    print(f"  FLT sections: {len(flt_sections)}, GC sections: {len(gc_sections)}")

    # For each section, compute Moran's I for HVGs
    section_svg_counts = {}  # section → n_significant_SVGs
    gene_morans = {}  # gene → {section: morans_I}

    for sec in sorted(sections):
        if sec not in mapping:
            continue
        sec_mask = adata.obs["section"] == sec
        adata_sec = adata[sec_mask]

        if "spatial" not in adata_sec.obsm:
            print(f"  WARNING: {sec} has no spatial coordinates, skipping")
            continue

        coords = adata_sec.obsm["spatial"]
        n_spots = adata_sec.n_obs
        print(f"  {sec} ({mapping[sec]['condition']}): {n_spots} spots, "
              f"computing Moran's I for {len(hvg_genes)} genes...")

        # Precompute spatial neighbors once per section
        neighbor_indices = precompute_spatial_neighbors(coords, k=6)

        X_sec = adata_sec[:, hvg_genes].X
        if sparse.issparse(X_sec):
            X_sec = X_sec.toarray()

        n_sig = 0
        for g_idx, gene in enumerate(hvg_genes):
            values = X_sec[:, g_idx].flatten()
            if np.std(values) == 0:
                continue

            I, z, p = compute_morans_i(values, neighbor_indices)

            if gene not in gene_morans:
                gene_morans[gene] = {}
            gene_morans[gene][sec] = {"I": I, "z": z, "p": p}

            if p < 0.05:
                n_sig += 1

        section_svg_counts[sec] = n_sig
        print(f"    → {n_sig} SVGs (p < 0.05)")

    # Approach A: Compare SVG counts FLT vs GC
    flt_svg_counts = [section_svg_counts.get(s, 0) for s in flt_sections]
    gc_svg_counts = [section_svg_counts.get(s, 0) for s in gc_sections]

    if flt_svg_counts and gc_svg_counts:
        svg_stat, svg_p = stats.mannwhitneyu(flt_svg_counts, gc_svg_counts,
                                              alternative="two-sided")
        svg_effect = cohens_d(np.array(flt_svg_counts),
                              np.array(gc_svg_counts))
    else:
        svg_stat, svg_p, svg_effect = float("nan"), float("nan"), float("nan")

    # Approach B: Per-gene comparison (exploratory)
    # For each gene, compare Moran's I across FLT vs GC sections
    per_gene_results = []
    for gene, sec_data in gene_morans.items():
        flt_i = [sec_data[s]["I"] for s in flt_sections if s in sec_data]
        gc_i = [sec_data[s]["I"] for s in gc_sections if s in sec_data]

        if len(flt_i) >= 2 and len(gc_i) >= 2:
            delta_i = np.mean(flt_i) - np.mean(gc_i)
            try:
                _, p = stats.mannwhitneyu(flt_i, gc_i,
                                           alternative="two-sided")
            except ValueError:
                p = 1.0
            per_gene_results.append({
                "gene": gene,
                "flt_mean_I": np.mean(flt_i),
                "gc_mean_I": np.mean(gc_i),
                "delta_I": delta_i,
                "p": p
            })

    per_gene_results.sort(key=lambda x: abs(x["delta_I"]), reverse=True)
    n_nominal_sig = sum(1 for r in per_gene_results if r["p"] < 0.05)

    # FLT-specific SVGs: significant in most FLT but not GC sections
    flt_specific = []
    gc_specific = []
    for gene, sec_data in gene_morans.items():
        flt_sig = sum(1 for s in flt_sections
                      if s in sec_data and sec_data[s]["p"] < 0.05)
        gc_sig = sum(1 for s in gc_sections
                     if s in sec_data and sec_data[s]["p"] < 0.05)
        n_flt = sum(1 for s in flt_sections if s in sec_data)
        n_gc = sum(1 for s in gc_sections if s in sec_data)

        if n_flt > 0 and n_gc > 0:
            flt_frac = flt_sig / n_flt
            gc_frac = gc_sig / n_gc
            if flt_frac >= 0.5 and gc_frac < 0.2:
                flt_specific.append(gene)
            if gc_frac >= 0.5 and flt_frac < 0.2:
                gc_specific.append(gene)

    # Jaccard overlap
    all_flt_svg = set()
    all_gc_svg = set()
    for gene, sec_data in gene_morans.items():
        if any(sec_data.get(s, {}).get("p", 1) < 0.05 for s in flt_sections):
            all_flt_svg.add(gene)
        if any(sec_data.get(s, {}).get("p", 1) < 0.05 for s in gc_sections):
            all_gc_svg.add(gene)

    union = all_flt_svg | all_gc_svg
    intersection = all_flt_svg & all_gc_svg
    jaccard = len(intersection) / len(union) if union else 0.0

    result = {
        "task": "F3b",
        "description": "Spatially variable genes: Moran's I FLT vs GC",
        "n_genes_tested": len(hvg_genes),
        "svg_count_comparison": {
            "flt_svg_counts": flt_svg_counts,
            "gc_svg_counts": gc_svg_counts,
            "flt_mean": np.mean(flt_svg_counts) if flt_svg_counts else 0,
            "gc_mean": np.mean(gc_svg_counts) if gc_svg_counts else 0,
            "mannwhitney_p": svg_p,
            "cohens_d": svg_effect
        },
        "per_gene_exploratory": {
            "n_nominal_sig_p05": n_nominal_sig,
            "n_total_tested": len(per_gene_results),
            "top_10_by_delta_I": per_gene_results[:10]
        },
        "condition_specific_svgs": {
            "flt_specific": flt_specific[:20],
            "gc_specific": gc_specific[:20],
            "n_flt_specific": len(flt_specific),
            "n_gc_specific": len(gc_specific)
        },
        "svg_overlap": {
            "n_flt_union": len(all_flt_svg),
            "n_gc_union": len(all_gc_svg),
            "n_intersection": len(intersection),
            "jaccard": jaccard
        },
        "section_svg_counts": section_svg_counts
    }

    print(f"\n  SVG counts: FLT mean={np.mean(flt_svg_counts):.0f}, "
          f"GC mean={np.mean(gc_svg_counts):.0f}, "
          f"p={svg_p:.3f}, d={svg_effect:.2f}")
    print(f"  Per-gene: {n_nominal_sig}/{len(per_gene_results)} nominal p<0.05")
    print(f"  FLT-specific SVGs: {len(flt_specific)}, "
          f"GC-specific: {len(gc_specific)}")
    print(f"  Jaccard overlap: {jaccard:.3f}")

    return result


# ── F3c: Brain Region Analysis ─────────────────────────────────────────────

def run_f3c(data_dir: Path, mapping: dict) -> dict:
    """F3c: Brain region-specific AUROC.

    Requires Seurat cluster metadata for region annotation.
    """
    print("\n" + "=" * 60)
    print("F3c: Brain Region-Specific AUROC")
    print("=" * 60)

    seurat_meta_path = data_dir / "seurat_metadata.csv"
    if not seurat_meta_path.exists():
        print("  SKIPPED: seurat_metadata.csv not found.")
        print("  F3c requires Seurat cluster annotations for brain region mapping.")
        return {
            "task": "F3c",
            "status": "skipped_no_seurat_metadata",
            "description": "Requires Seurat RDS → metadata extraction first"
        }

    import anndata as ad
    processed_dir = data_dir / "processed"
    merged_path = processed_dir / "osd352_visium_merged.h5ad"

    if not merged_path.exists():
        return {"task": "F3c", "status": "data_not_found"}

    adata = ad.read_h5ad(merged_path)
    meta = pd.read_csv(seurat_meta_path, index_col=0)

    # Find cluster column in Seurat metadata
    cluster_col = None
    for col in meta.columns:
        if "cluster" in col.lower() or col == "seurat_clusters":
            cluster_col = col
            break

    if cluster_col is None:
        return {"task": "F3c", "status": "no_cluster_column"}

    print(f"  Cluster column: {cluster_col}")
    print(f"  Unique clusters: {sorted(meta[cluster_col].unique())}")

    # Merge cluster labels into adata
    # Seurat barcode format may differ from SpaceRanger
    # Try direct merge by index
    common_barcodes = adata.obs.index.intersection(meta.index)
    if len(common_barcodes) < adata.n_obs * 0.5:
        # Try stripping suffixes
        meta_stripped = meta.copy()
        meta_stripped.index = meta_stripped.index.str.replace(
            r"_\d+$", "", regex=True)
        common_barcodes = adata.obs.index.intersection(meta_stripped.index)

    if len(common_barcodes) < 100:
        print(f"  WARNING: Only {len(common_barcodes)} matching barcodes. "
              f"Barcode format mismatch?")
        return {"task": "F3c", "status": "barcode_mismatch",
                "n_matched": len(common_barcodes)}

    # Assign clusters
    adata = adata[common_barcodes].copy()
    adata.obs["cluster"] = meta.loc[common_barcodes, cluster_col].values

    # Group clusters into brain regions (to be refined with paper)
    # For now, use raw clusters
    clusters = adata.obs["cluster"].unique()
    print(f"  {len(clusters)} clusters, {adata.n_obs} spots with annotations")

    # Per-cluster pseudo-bulk → LOAO
    region_results = {}
    for cluster in sorted(clusters):
        cluster_mask = adata.obs["cluster"] == cluster
        adata_cluster = adata[cluster_mask]

        # Build per-section pseudo-bulk for this cluster
        sections_with_cluster = adata_cluster.obs["section"].unique()
        if len(sections_with_cluster) < 4:
            region_results[str(cluster)] = {
                "status": "too_few_sections",
                "n_sections": len(sections_with_cluster)
            }
            continue

        # Build feature matrix
        section_names_c = []
        animals_c = []
        y_c = []
        X_list = []

        for sec in sorted(sections_with_cluster):
            if sec not in mapping:
                continue
            sec_mask = (adata_cluster.obs["section"] == sec)
            adata_sec = adata_cluster[sec_mask]

            n_spots = adata_sec.n_obs
            if n_spots < 5:
                continue

            X_sec = adata_sec.X
            if sparse.issparse(X_sec):
                X_sec = X_sec.toarray()
            pb = np.mean(X_sec, axis=0).astype(np.float32)

            X_list.append(pb)
            section_names_c.append(sec)
            animals_c.append(mapping[sec]["animal"])
            y_c.append(1 if mapping[sec]["condition"] == "FLT" else 0)

        if len(X_list) < 4 or len(set(y_c)) < 2:
            region_results[str(cluster)] = {
                "status": "insufficient_data",
                "n_sections": len(X_list)
            }
            continue

        X_c = np.array(X_list)
        y_c = np.array(y_c)
        animals_c = np.array(animals_c)

        result = loao_classification(X_c, y_c, animals_c)
        result["n_spots_total"] = int(cluster_mask.sum())
        result["n_sections"] = len(X_list)
        region_results[str(cluster)] = result

        auroc = result.get("auroc", float("nan"))
        print(f"  Cluster {cluster}: AUROC={auroc:.3f}, "
              f"n_spots={result['n_spots_total']}, "
              f"n_sections={result['n_sections']}")

    return {
        "task": "F3c",
        "description": "Brain region/cluster-specific AUROC",
        "n_clusters": len(clusters),
        "method": f"Per-cluster pseudo-bulk → PCA(≤{MAX_PCA})-LR LOAO, "
                  f"exact permutation",
        "cluster_results": region_results
    }


# ── F3d: Spatial vs Bulk Comparison ────────────────────────────────────────

def run_f3d(data_dir: Path, mapping: dict) -> dict:
    """F3d: Spatial pseudo-bulk vs companion bulk RNA-seq comparison."""
    print("\n" + "=" * 60)
    print("F3d: Spatial vs Companion Bulk Comparison")
    print("=" * 60)

    processed_dir = data_dir / "processed"

    # Load spatial pseudo-bulk (section-level)
    X_spatial = np.load(processed_dir / "pseudobulk_matrix.npy")
    gene_names_spatial = np.load(processed_dir / "gene_names.npy",
                                  allow_pickle=True)
    with open(processed_dir / "section_names.json") as f:
        section_names = json.load(f)

    # Load companion bulk RNA-seq
    bulk_counts_path = data_dir / "GLDS-352_rna_seq_Normalized_Counts.csv"
    bulk_sample_path = data_dir / "GLDS-352_rna_seq_SampleTable.csv"

    if not bulk_counts_path.exists():
        return {"task": "F3d", "status": "bulk_data_not_found"}

    # Bulk counts: genes (ENSMUSG) × samples
    bulk_df = pd.read_csv(bulk_counts_path, index_col=0)
    bulk_samples = pd.read_csv(bulk_sample_path, index_col=0)

    print(f"  Bulk: {bulk_df.shape} (genes × samples)")
    print(f"  Bulk samples: {list(bulk_df.columns)}")
    # Print condition info (column name may vary across OSDR datasets)
    cond_col = None
    for col in bulk_samples.columns:
        if col.lower() in ("condition", "group", "factor_value",
                           "factor.value", "sample_type"):
            cond_col = col
            break
    if cond_col:
        print(f"  Bulk conditions ({cond_col}): "
              f"{dict(bulk_samples[cond_col])}")
    else:
        print(f"  Bulk sample columns: {list(bulk_samples.columns)}")

    # Map bulk samples to condition and animal
    # Sample names: RR3_BRN_{condition}_{animal}
    bulk_y = []
    bulk_animals = []
    for sample in bulk_df.columns:
        parts = sample.split("_")
        if "FLT" in parts:
            bulk_y.append(1)
        elif "GC" in parts:
            bulk_y.append(0)
        else:
            bulk_y.append(-1)
        bulk_animals.append(parts[-1])  # F1, F2, F7, G7, G8, G9

    bulk_y = np.array(bulk_y)
    bulk_animals = np.array(bulk_animals)

    # Remove samples with unknown condition
    valid = bulk_y >= 0
    bulk_df = bulk_df.iloc[:, valid]
    bulk_y = bulk_y[valid]
    bulk_animals = bulk_animals[valid]

    # Transpose: samples × genes
    X_bulk = bulk_df.values.T.astype(np.float32)
    gene_names_bulk = bulk_df.index.values

    # log1p transform if not already (check if values suggest raw counts)
    if X_bulk.max() > 100:
        print("  Bulk values > 100, applying log1p transform")
        X_bulk = np.log1p(X_bulk)

    print(f"  Bulk matrix: {X_bulk.shape} (samples × genes)")

    # Bulk LOAO
    print("\n  [Bulk LOAO]")
    bulk_result = loao_classification(X_bulk, bulk_y, bulk_animals)
    print(f"  Bulk AUROC = {bulk_result['auroc']:.3f} "
          f"[{bulk_result['ci_low']:.3f}, {bulk_result['ci_high']:.3f}]")
    print(f"  Exact perm p = {bulk_result['p_exact']:.3f}")

    # Spatial animal-level LOAO (for fair comparison with bulk n=6)
    # Average sections per animal
    unique_animals = np.unique(
        [mapping[s]["animal"] for s in section_names if s in mapping])
    X_spatial_animal = []
    y_spatial_animal = []
    animals_spatial = []

    for a in unique_animals:
        sec_indices = [i for i, s in enumerate(section_names)
                       if s in mapping and mapping[s]["animal"] == a]
        if sec_indices:
            X_spatial_animal.append(X_spatial[sec_indices].mean(axis=0))
            y_spatial_animal.append(
                1 if mapping[section_names[sec_indices[0]]]["condition"] == "FLT"
                else 0)
            animals_spatial.append(a)

    X_spatial_animal = np.array(X_spatial_animal)
    y_spatial_animal = np.array(y_spatial_animal)
    animals_spatial = np.array(animals_spatial)

    print("\n  [Spatial animal-level LOAO]")
    spatial_result = loao_classification(X_spatial_animal, y_spatial_animal,
                                          animals_spatial)
    print(f"  Spatial AUROC = {spatial_result['auroc']:.3f} "
          f"[{spatial_result['ci_low']:.3f}, {spatial_result['ci_high']:.3f}]")
    print(f"  Exact perm p = {spatial_result['p_exact']:.3f}")

    # Delta AUROC
    delta_auroc = spatial_result["auroc"] - bulk_result["auroc"]
    print(f"\n  Δ AUROC (spatial - bulk) = {delta_auroc:+.3f}")

    # DEG overlap
    bulk_deg_path = data_dir / "GLDS-352_rna_seq_differential_expression.csv"
    deg_overlap = None
    if bulk_deg_path.exists():
        bulk_deg = pd.read_csv(bulk_deg_path, index_col=0)
        # Find padj column
        padj_col = None
        for col in bulk_deg.columns:
            if "adj" in col.lower() and "p" in col.lower():
                padj_col = col
                break
        if padj_col is None:
            for col in bulk_deg.columns:
                if "fdr" in col.lower():
                    padj_col = col
                    break

        if padj_col is not None:
            bulk_deg_genes = set(bulk_deg.index[bulk_deg[padj_col] < 0.05])
            print(f"  Bulk DEGs (padj<0.05): {len(bulk_deg_genes)}")
            deg_overlap = {
                "bulk_deg_count": len(bulk_deg_genes),
                "padj_column_used": padj_col
            }
        else:
            print(f"  WARNING: No padj column found. Columns: {list(bulk_deg.columns)}")

    result = {
        "task": "F3d",
        "description": "Spatial pseudo-bulk vs companion bulk RNA-seq comparison",
        "bulk_auroc": bulk_result,
        "spatial_auroc": spatial_result,
        "delta_auroc": delta_auroc,
        "comparison": {
            "spatial_n_animals": len(unique_animals),
            "bulk_n_samples": len(bulk_y),
            "same_mission": True,
            "same_cv_design": "LOAO, animal-level",
            "apple_to_apple": True
        }
    }
    if deg_overlap:
        result["deg_overlap"] = deg_overlap

    return result


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="OSD-352 RR-3 Brain Visium Benchmark (F3a-d)")
    parser.add_argument("--data-dir", type=str, default=str(DATA_DIR))
    parser.add_argument("--tasks", type=str, default="F3a,F3b,F3c,F3d",
                        help="Comma-separated tasks to run")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    tasks = [t.strip().upper() for t in args.tasks.split(",")]

    print("=" * 60)
    print("OSD-352 RR-3 Brain Visium Benchmark")
    print("=" * 60)

    # Load section mapping
    print("\nLoading section → animal/condition mapping...")
    mapping = load_section_condition_mapping(data_dir)

    results = {}

    if "F3A" in tasks:
        results["F3a"] = run_f3a(data_dir, mapping)
        save_json(results["F3a"], EVAL_DIR / "F3a_visium_classification.json")

    if "F3B" in tasks:
        results["F3b"] = run_f3b(data_dir, mapping)
        save_json(results["F3b"], EVAL_DIR / "F3b_visium_svg.json")

    if "F3C" in tasks:
        results["F3c"] = run_f3c(data_dir, mapping)
        save_json(results["F3c"], EVAL_DIR / "F3c_visium_regions.json")

    if "F3D" in tasks:
        results["F3d"] = run_f3d(data_dir, mapping)
        save_json(results["F3d"],
                  EVAL_DIR / "F3d_visium_cross_resolution.json")

    print("\n" + "=" * 60)
    print("DONE — All tasks complete")
    print("=" * 60)


if __name__ == "__main__":
    main()

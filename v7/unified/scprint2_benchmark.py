"""
scPRINT-2 Benchmark — v7 GeneLabBench
=======================================
Evaluate scPRINT-2 (Gonzalez-Blas et al., 2025) zero-shot embeddings
on mouse spaceflight bulk RNA-seq classification.

scPRINT-2: trained on 350M cells, 16 organisms (including mouse),
single-cell foundation model. We treat each bulk RNA-seq sample as a
"pseudo-cell" and extract zero-shot embeddings → LR classifier.

Scientific question: Does multi-organism (350M cell) pretraining
finally beat PCA-LR on spaceflight bulk RNA-seq?

Protocol:
  1. Load bulk RNA-seq (log2-normalized, ~17-28K genes per tissue)
  2. Convert to AnnData (pseudo-cells): each sample = 1 cell
  3. Preprocess via scdataloader.Preprocessor (maps to scPRINT gene vocab)
  4. Extract zero-shot embeddings via Python API (Embedder class)
  5. LOMO-CV: LR on embeddings (512-dim), compare to PCA-LR baseline

Usage:
  python scprint2_benchmark.py --tissue liver
  python scprint2_benchmark.py --tissue all

Requires:
  pip install scprint2
  Default checkpoint: medium-v1.5.ckpt
  Legacy local filename medium-v1.5_fixed.ckpt is also accepted
  LaminDB+bionty instance initialized with mouse/human organisms

Output:
  v7/evaluation/SCPRINT2_{tissue}.json
"""

import argparse
import json
import re
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR   = Path(__file__).resolve().parent.parent.parent
PROC_DIR   = BASE_DIR / "processed" / "A_detection"
M1_SUMMARY = BASE_DIR / "v4" / "evaluation" / "M1_summary.json"
OUT_DIR    = BASE_DIR / "v7" / "evaluation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Reuse from gnn_wgcna.py
TISSUE_MISSIONS = {
    "liver":         ["RR-1", "RR-3", "RR-6", "RR-8", "RR-9", "MHU-2"],
    "gastrocnemius": ["RR-1", "RR-5", "RR-9"],
    "kidney":        ["RR-1", "RR-3", "RR-7"],
    "thymus":        ["MHU-1", "MHU-2", "RR-6", "RR-9"],
    "eye":           ["RR-1", "RR-3", "OSD-397"],
    "skin":          ["MHU-2", "RR-6", "RR-7"],
}

MISSION_ALIASES = {"TBD": "OSD-397"}

EXPR_FILE = {
    "liver":         "liver_all_missions_log2_norm_limma_rbe.csv",
    "gastrocnemius": "gastrocnemius_all_missions_log2_norm_limma_rbe.csv",
    "kidney":        "kidney_all_missions_log2_norm_limma_rbe.csv",
    "thymus":        "thymus_all_missions_log2_norm_limma_rbe.csv",
    "eye":           "eye_all_missions_log2_norm_limma_rbe.csv",
    "skin":          "skin_all_missions_log2_norm.csv",
}
META_FILE = {t: f"{t}_all_missions_metadata.csv" for t in TISSUE_MISSIONS}

LABEL_MAP = {
    "FLT": 1, "Flight": 1,
    "GC": 0, "Ground Control": 0, "Ground": 0,
    "Basal": 0, "BC": 0, "VC": 0, "Vivarium": 0,
}
EXCLUDE_LABELS = {"AG", "Artificial Gravity"}


# ── Data Loading ──────────────────────────────────────────────────────────────
def load_tissue_data(tissue):
    """Load expression + metadata + labels."""
    tdir = PROC_DIR / tissue
    expr = pd.read_csv(tdir / EXPR_FILE[tissue], index_col=0)   # samples × genes
    meta = pd.read_csv(tdir / META_FILE[tissue])

    # Align metadata index with expression index
    if "sample_name" in meta.columns:
        meta = meta.set_index("sample_name")
    elif "Unnamed: 0" in meta.columns:
        meta = meta.set_index("Unnamed: 0")
    common = expr.index.intersection(meta.index)
    expr = expr.loc[common]
    meta = meta.loc[common]

    # Find label column
    for candidate in ["label", "condition", "Label", "Condition"]:
        if candidate in meta.columns:
            label_col = candidate
            break
    else:
        label_col = meta.columns[0]
    y_raw = meta[label_col].astype(str).str.strip()
    valid = y_raw.apply(lambda v: v not in EXCLUDE_LABELS and v in LABEL_MAP)
    expr = expr.loc[valid]
    meta = meta.loc[valid]
    y = y_raw[valid].map(LABEL_MAP).values.astype(np.int32)

    if "mission" not in meta.columns:
        meta = meta.copy()
        meta["mission"] = "unknown"
    else:
        meta = meta.copy()
        meta["mission"] = meta["mission"].replace(MISSION_ALIASES)

    # Drop any non-numeric columns from expression (e.g. skin has mission/osd_id)
    numeric_cols = expr.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) < len(expr.columns):
        dropped = set(expr.columns) - set(numeric_cols)
        print(f"  Dropping non-numeric columns from expression: {dropped}")
        expr = expr[numeric_cols]

    return expr, y, meta


# ── scPRINT AnnData Preparation ───────────────────────────────────────────────
def expr_to_anndata(expr_df, y, meta_df, tissue):
    """
    Convert bulk RNA-seq DataFrame to AnnData for scPRINT.
    expr_df: samples × genes (log2-normalized)
    Returns: anndata.AnnData with pseudo-count X, organism annotation.
    """
    import anndata as ad

    # scPRINT expects count-like data. Convert log2-norm to pseudo-counts.
    X = np.round(np.power(2.0, expr_df.values.astype(np.float32)) - 1.0)
    X = np.clip(X, 0, None)  # ensure non-negative

    obs = meta_df.copy()
    obs["condition_label"] = y
    obs["organism_ontology_term_id"] = "NCBITaxon:10090"
    obs["is_primary_data"] = True  # avoid "too many secondary" filter

    var = pd.DataFrame(index=expr_df.columns)
    var.index.name = "gene_id"

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names = expr_df.columns.tolist()
    print(f"  AnnData: {adata.n_obs} samples × {adata.n_vars} genes")
    return adata


def resolve_scprint2_checkpoint(ckpt_arg):
    """Resolve a scPRINT-2 checkpoint path."""
    requested = Path(ckpt_arg).expanduser()
    model_dir = BASE_DIR / "v7" / "models" / "scprint2"
    candidates = []

    def add_candidate(path):
        if path not in candidates:
            candidates.append(path)

    add_candidate(requested)
    if not requested.is_absolute():
        add_candidate(BASE_DIR / requested)
        add_candidate(model_dir / requested.name)

    if requested.name == "medium-v1.5.ckpt":
        add_candidate(model_dir / "medium-v1.5.ckpt")
        add_candidate(model_dir / "medium-v1.5_fixed.ckpt")
    elif requested.name == "medium-v1.5_fixed.ckpt":
        add_candidate(model_dir / "medium-v1.5_fixed.ckpt")
        add_candidate(model_dir / "medium-v1.5.ckpt")

    for ckpt_path in candidates:
        if ckpt_path.exists():
            return ckpt_path.resolve()

    raise FileNotFoundError(
        f"Checkpoint not found: {ckpt_arg}. "
        "Run v7/scripts/hpc_setup_v7.sh or place a checkpoint under "
        "v7/models/scprint2/."
    )


def load_scprint_model(ckpt_path, device="cuda"):
    """Load scPRINT-2 model from a fixed checkpoint using Python API."""
    import torch
    from scprint2 import scPRINT2

    print(f"  Loading scPRINT-2 model from {Path(ckpt_path).name}...")
    t0 = time.time()

    model = scPRINT2.load_from_checkpoint(
        str(ckpt_path),
        precpt_gene_emb=None,
        gene_pos_file=None,
        map_location=device if device != "cuda" else "cpu",
        strict=False,
    )
    model.eval()
    if device == "cuda":
        import torch
        if torch.cuda.is_available():
            model = model.cuda()
    print(f"  Model loaded ({time.time()-t0:.1f}s)")
    return model


def extract_scprint_embeddings(adata, model):
    """Extract embeddings using the Python API (Embedder + Preprocessor)."""
    from scdataloader import Preprocessor
    from scprint2.tasks.cell_emb import Embedder

    print(f"  Extracting scPRINT-2 embeddings (n={adata.n_obs})...")
    t0 = time.time()

    # Preprocess: maps genes to scPRINT vocabulary, normalizes
    preprocessor = Preprocessor(
        do_postp=False,
        force_preprocess=True,
        skip_validate=True,
        drop_non_primary=False,
        min_dataset_size=1,
        min_valid_genes_id=0,
        min_nnz_genes=0,
        maxdropamount=100,
    )
    adata = preprocessor(adata)
    print(f"  After preprocessing: {adata.n_obs} × {adata.n_vars} genes")

    # Extract embeddings
    embedder = Embedder(
        batch_size=64,
        num_workers=0,
        max_len=2000,
        doclass=False,
        doplot=False,
        how="random expr",
        pred_embedding=["cell_type_ontology_term_id"],
    )
    adata, metrics = embedder(model, adata)

    # Find embedding key
    emb_key = None
    for key in ["scprint_emb", "X_scprint", "X_scprint2"]:
        if key in adata.obsm:
            emb_key = key
            break
    if emb_key is None:
        for key in adata.obsm.keys():
            arr = np.asarray(adata.obsm[key])
            if arr.ndim == 2 and arr.shape[0] == adata.n_obs and arr.shape[1] > 1:
                emb_key = key
                break
    if emb_key is None:
        raise KeyError(f"No embedding found in adata.obsm. Keys: {list(adata.obsm.keys())}")

    emb = np.asarray(adata.obsm[emb_key], dtype=np.float32)
    print(f"  Embeddings: {emb.shape} (key={emb_key}, {time.time()-t0:.1f}s)")
    return emb


def cache_slug_for_checkpoint(ckpt_path):
    """Generate a stable cache key for a checkpoint path."""
    stem = Path(checkpoint_record_value(ckpt_path)).stem
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("_")
    return slug or "default"


def embedding_cache_file(cache_dir, tissue, ckpt_path):
    """Return the checkpoint-specific embedding cache path."""
    cache_dir = Path(cache_dir)
    return cache_dir / f"SCPRINT2_emb_{tissue}__{cache_slug_for_checkpoint(ckpt_path)}.npy"


def checkpoint_record_value(ckpt_path):
    """Serialize checkpoint paths portably for results and cache matching."""
    if not ckpt_path:
        return ""
    path = Path(ckpt_path).expanduser()
    canonical_name = "medium-v1.5.ckpt" if path.name == "medium-v1.5_fixed.ckpt" else path.name
    parts = path.parts
    if "v7" in parts:
        v7_idx = parts.index("v7")
        return str(Path(*parts[v7_idx:]).with_name(canonical_name))
    return str(path.with_name(canonical_name))


def result_matches_request(result, *, tissue, ckpt_path, n_bootstrap, n_perm, seed):
    """Check whether an existing result matches the current invocation."""
    expected = {
        "method": "scPRINT-2",
        "tissue": tissue,
        "ckpt_path": checkpoint_record_value(ckpt_path),
        "n_bootstrap": n_bootstrap,
        "n_perm": n_perm,
        "seed": seed,
    }
    for key, value in expected.items():
        if key == "ckpt_path":
            if checkpoint_record_value(result.get(key, "")) != value:
                return False
            continue
        if result.get(key) != value:
            return False
    return True


# ── LOMO-CV ────────────────────────────────────────────────────────────────────
def lomo_cv_on_embeddings(embeddings, y, meta, missions, n_bootstrap, n_perm, seed):
    """Run LOMO-CV with LR on embeddings. Returns fold results."""
    rng = np.random.default_rng(seed)
    fold_results = []

    for test_mission in missions:
        test_mask  = (meta["mission"].values == test_mission)
        train_mask = ~test_mask

        if test_mask.sum() == 0:
            continue
        if len(np.unique(y[test_mask])) < 2:
            print(f"    [SKIP] {test_mission}: single class in test")
            continue
        if train_mask.sum() < 4:
            print(f"    [SKIP] {test_mission}: too few train samples")
            continue

        emb_train = embeddings[train_mask]
        emb_test  = embeddings[test_mask]
        y_train   = y[train_mask]
        y_test    = y[test_mask]

        # Standardize embeddings (fit on train only)
        scaler = StandardScaler()
        emb_train_s = scaler.fit_transform(emb_train)
        emb_test_s  = scaler.transform(emb_test)

        # LR with L2 regularization (matches PCA-LR pipeline)
        lr = LogisticRegression(
            C=1.0, max_iter=1000, solver="lbfgs",
            class_weight="balanced", random_state=seed
        )
        lr.fit(emb_train_s, y_train)
        y_score = lr.predict_proba(emb_test_s)[:, 1]
        auroc   = roc_auc_score(y_test, y_score)

        # Bootstrap CI
        boot_aurocs = []
        for _ in range(n_bootstrap):
            idx = rng.integers(0, len(y_test), size=len(y_test))
            if len(np.unique(y_test[idx])) < 2:
                continue
            try:
                boot_aurocs.append(roc_auc_score(y_test[idx], y_score[idx]))
            except Exception:
                continue
        ci_lo = np.percentile(boot_aurocs, 2.5)  if boot_aurocs else np.nan
        ci_hi = np.percentile(boot_aurocs, 97.5) if boot_aurocs else np.nan

        print(f"    {test_mission}: AUROC={auroc:.3f} [{ci_lo:.3f}, {ci_hi:.3f}]")
        fold_results.append({
            "test_mission": test_mission,
            "auroc": round(float(auroc), 4),
            "ci_lower": round(float(ci_lo), 4),
            "ci_upper": round(float(ci_hi), 4),
            "n_train": int(train_mask.sum()),
            "n_test":  int(test_mask.sum()),
        })

    # Mean + permutation test
    if not fold_results:
        return fold_results, None, None

    mean_auroc = float(np.mean([f["auroc"] for f in fold_results]))

    # Permutation test: shuffle y, retrain LR on embeddings
    # Fast since LR is cheap (embeddings are fixed)
    perm_rng   = np.random.default_rng(seed + 7777)
    perm_means = []
    for perm_i in range(n_perm):
        y_perm = perm_rng.permutation(y)
        fold_aurocs_p = []

        for test_mission in missions:
            test_mask  = (meta["mission"].values == test_mission)
            train_mask = ~test_mask
            if test_mask.sum() == 0 or len(np.unique(y_perm[test_mask])) < 2:
                continue

            emb_train_p = scaler.fit_transform(embeddings[train_mask])
            emb_test_p  = scaler.transform(embeddings[test_mask])
            y_train_p   = y_perm[train_mask]
            y_test_p    = y_perm[test_mask]

            lr_p = LogisticRegression(
                C=1.0, max_iter=500, solver="liblinear",
                class_weight="balanced", random_state=seed + perm_i
            )
            try:
                lr_p.fit(emb_train_p, y_train_p)
                y_score_p = lr_p.predict_proba(emb_test_p)[:, 1]
                fold_aurocs_p.append(roc_auc_score(y_test_p, y_score_p))
            except Exception:
                pass

        if fold_aurocs_p:
            perm_means.append(np.mean(fold_aurocs_p))

    perm_p = (np.sum(np.array(perm_means) >= mean_auroc) + 1) / (len(perm_means) + 1)
    return fold_results, mean_auroc, perm_p


def bootstrap_mean_ci_from_folds(fold_results, seed):
    """Estimate a CI on the mean AUROC by resampling folds."""
    if not fold_results:
        return None, None

    values = np.array([f["auroc"] for f in fold_results], dtype=np.float32)
    if len(values) == 1:
        return float(values[0]), float(values[0])

    rng = np.random.default_rng(seed)
    boot_means = []
    for _ in range(5000):
        idx = rng.integers(0, len(values), size=len(values))
        boot_means.append(float(np.mean(values[idx])))

    return float(np.percentile(boot_means, 2.5)), float(np.percentile(boot_means, 97.5))


def load_pcalr_baseline(tissue):
    try:
        with open(M1_SUMMARY) as f:
            d = json.load(f)
        return d[tissue]["gene"]["pca_lr"]["auroc"]
    except Exception:
        return None


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tissue", default="liver",
                        help="Tissue name or 'all'")
    parser.add_argument("--ckpt-path", default="medium-v1.5.ckpt",
                        help="Local checkpoint path; legacy _fixed filenames are also accepted")
    parser.add_argument("--model-id", default=None,
                        help="Deprecated alias for --ckpt-path")
    parser.add_argument("--n-bootstrap", type=int, default=1000)
    parser.add_argument("--n-perm",      type=int, default=100)
    parser.add_argument("--seed",        type=int, default=42)
    parser.add_argument("--embedding-cache-dir", default=None,
                        help="Directory to cache per-tissue embeddings (optional)")
    args = parser.parse_args()

    if args.model_id and args.ckpt_path == "medium-v1.5.ckpt":
        args.ckpt_path = args.model_id

    print("Resolving scPRINT-2 checkpoint...")
    ckpt_path = resolve_scprint2_checkpoint(args.ckpt_path)
    print(f"  Checkpoint: {ckpt_path}")

    # Determine device
    import torch
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"  Device: {device}")

    # Load model once (shared across tissues)
    model = load_scprint_model(ckpt_path, device=device)

    # ── Process tissues ───────────────────────────────────────────────────────
    tissues = list(TISSUE_MISSIONS.keys()) if args.tissue == "all" else [args.tissue]

    for tissue in tissues:
        out_file = OUT_DIR / f"SCPRINT2_{tissue}.json"
        if out_file.exists():
            with open(out_file) as f:
                existing = json.load(f)
            if result_matches_request(
                existing,
                tissue=tissue,
                ckpt_path=ckpt_path,
                n_bootstrap=args.n_bootstrap,
                n_perm=args.n_perm,
                seed=args.seed,
            ):
                print(f"[SKIP] {out_file.name} already exists for this parameter set")
                continue
            print(f"[RECOMPUTE] {out_file.name} exists but parameters changed")

        print(f"\n{'='*60}")
        print(f"Tissue: {tissue}")
        print(f"{'='*60}")

        # Load expression data
        expr, y, meta = load_tissue_data(tissue)
        missions = TISSUE_MISSIONS[tissue]

        # Check embedding cache
        cache_dir = Path(args.embedding_cache_dir) if args.embedding_cache_dir else OUT_DIR / "embeddings_cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_file = embedding_cache_file(cache_dir, tissue, ckpt_path)

        if cache_file.exists():
            print(f"  Loading cached embeddings: {cache_file.name}")
            embeddings = np.load(cache_file)
            print(f"  Embeddings: {embeddings.shape}")
        else:
            # Convert to AnnData
            adata = expr_to_anndata(expr, y, meta, tissue)

            # Extract scPRINT embeddings via Python API
            embeddings = extract_scprint_embeddings(adata, model)

            np.save(cache_file, embeddings)
            print(f"  Cached embeddings to: {cache_file.name}")

        emb_dim = embeddings.shape[1]
        print(f"  Embedding dim: {emb_dim}")

        # LOMO-CV on embeddings
        print(f"  LOMO-CV ({len(missions)} missions)...")
        fold_results, mean_auroc, perm_p = lomo_cv_on_embeddings(
            embeddings, y, meta, missions,
            args.n_bootstrap, args.n_perm, args.seed
        )

        if mean_auroc is not None:
            print(f"  Mean AUROC: {mean_auroc:.4f}, perm_p: {perm_p:.3f}")

        mean_ci_lo, mean_ci_hi = bootstrap_mean_ci_from_folds(
            fold_results, args.seed + 2026
        )
        pca_lr_auroc = load_pcalr_baseline(tissue)

        result = {
            "method":       "scPRINT-2",
            "ckpt_path":    checkpoint_record_value(ckpt_path),
            "tissue":       tissue,
            "embedding_dim": emb_dim,
            "mean_auroc":   round(float(mean_auroc), 4) if mean_auroc is not None else None,
            "mean_ci_lower": round(float(mean_ci_lo), 4) if mean_ci_lo is not None else None,
            "mean_ci_upper": round(float(mean_ci_hi), 4) if mean_ci_hi is not None else None,
            "perm_pvalue":  round(float(perm_p), 4)    if perm_p is not None else None,
            "pca_lr_auroc": pca_lr_auroc,
            "delta_vs_pca_lr": round(float(mean_auroc - pca_lr_auroc), 4)
                               if (mean_auroc is not None and pca_lr_auroc is not None) else None,
            "n_folds":      len(fold_results),
            "folds":        fold_results,
            "seed":         args.seed,
            "n_bootstrap":  args.n_bootstrap,
            "n_perm":       args.n_perm,
        }

        with open(out_file, "w") as f:
            json.dump(result, f, indent=2)
        print(f"  Saved: {out_file.name}")

    # Print summary table
    print("\n" + "="*70)
    print(f"{'Tissue':<14} {'scPRINT-2':<10} {'perm_p':<10} {'PCA-LR':<10} {'Delta'}")
    print("-"*70)
    for tissue in tissues:
        out_file = OUT_DIR / f"SCPRINT2_{tissue}.json"
        if out_file.exists():
            with open(out_file) as f:
                r = json.load(f)
            auroc  = r.get("mean_auroc")
            perm_p = r.get("perm_pvalue")
            pcalr  = r.get("pca_lr_auroc")
            sig = "*" if perm_p and perm_p < 0.05 else ""
            delta = f"{auroc - pcalr:+.3f}" if (auroc is not None and pcalr is not None) else "N/A"
            print(f"{tissue:<14} {str(auroc):<10} {(str(perm_p)+sig):<10} {str(pcalr):<10} {delta}")

    if len(tissues) > 1:
        summary = {
            "method": "scPRINT-2",
            "ckpt_path": checkpoint_record_value(ckpt_path),
            "tissues": tissues,
            "results": [],
        }
        for tissue in tissues:
            out_file = OUT_DIR / f"SCPRINT2_{tissue}.json"
            if out_file.exists():
                with open(out_file) as f:
                    summary["results"].append(json.load(f))
        summary_path = OUT_DIR / "SCPRINT2_summary.json"
        with open(summary_path, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"\nSummary saved: {summary_path}")


if __name__ == "__main__":
    main()

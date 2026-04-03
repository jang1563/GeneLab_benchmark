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
  3. Map ENSMUSG → scPRINT gene vocabulary (ENSMUSG supported for mouse)
  4. Extract zero-shot embeddings via the official `scprint2 embed` workflow
  5. LOMO-CV: LR on embeddings (512-dim), compare to PCA-LR baseline
  6. Also test: fine-tuned LR on all-but-test-mission embeddings

Usage:
  python scprint2_benchmark.py --tissue liver
  python scprint2_benchmark.py --tissue all

Requires:
  pip install scprint2
  Default checkpoint: medium-v1.5.ckpt (from jkobject/scPRINT)

Output:
  v7/evaluation/SCPRINT2_{tissue}.json
"""

import argparse
import json
import subprocess
import tempfile
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
    "eye":           ["RR-1", "RR-3", "TBD"],
    "skin":          ["MHU-2", "RR-6", "RR-7"],
}

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

    return expr, y, meta


# ── scPRINT AnnData Preparation ───────────────────────────────────────────────
def expr_to_anndata(expr_df, y, meta_df, tissue):
    """
    Convert bulk RNA-seq DataFrame to AnnData for scPRINT.
    expr_df: samples × genes (log2-normalized)
    Returns: anndata.AnnData with:
      - X = raw-like expression (we use log2-norm directly since scPRINT can handle)
      - obs: metadata including condition, mission
      - var_names: ENSMUSG IDs (scPRINT mouse vocab)
    """
    import anndata as ad

    # scPRINT expects count-like data in X, but we have log2-normalized.
    # Use 2^(expr) - 1 to approximate counts (pseudo-counts), then round.
    # scPRINT's encoder handles log normalization internally.
    # Alternative: provide log2-norm directly and specify log_norm=True.
    # We'll use the log2-normalized directly with normalization flag.
    X = expr_df.values.astype(np.float32)  # (n_samples, n_genes)

    obs = meta_df.copy()
    obs["condition_label"] = y
    obs["n_counts"] = X.sum(axis=1)
    obs["organism_ontology_term_id"] = "NCBITaxon:10090"
    obs["tissue"] = tissue

    var = pd.DataFrame(index=expr_df.columns)
    var.index.name = "gene_id"
    # scPRINT uses ENSEMBL IDs — our gene IDs may be ENSMUSG (correct for mouse)
    # or symbol. Keep as-is for mouse organism.
    var["gene_name"] = var.index.tolist()

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names = expr_df.columns.tolist()
    print(f"  AnnData: {adata.n_obs} samples × {adata.n_vars} genes")
    return adata


def resolve_scprint2_checkpoint(ckpt_arg):
    """Resolve a scPRINT-2 checkpoint path, downloading from HF if needed."""
    ckpt_path = Path(ckpt_arg).expanduser()
    if ckpt_path.exists():
        return ckpt_path.resolve()

    filename = ckpt_path.name if ckpt_path.suffix == ".ckpt" else "medium-v1.5.ckpt"
    repo_candidates = ["jkobject/scPRINT"]

    # Allow repo-qualified inputs like jkobject/scPRINT/v2-medium.ckpt.
    arg_parts = ckpt_arg.replace("\\", "/").split("/")
    if len(arg_parts) >= 3 and arg_parts[-1].endswith(".ckpt"):
        repo_candidates = ["/".join(arg_parts[:-1])] + repo_candidates
        filename = arg_parts[-1]

    try:
        from huggingface_hub import hf_hub_download
    except ImportError as exc:
        raise FileNotFoundError(
            f"Checkpoint not found locally: {ckpt_arg}. "
            "Install huggingface_hub or pre-download the checkpoint."
        ) from exc

    local_dir = OUT_DIR / "models" / "scprint2"
    local_dir.mkdir(parents=True, exist_ok=True)
    last_error = None
    for repo_id in repo_candidates:
        try:
            downloaded = hf_hub_download(repo_id=repo_id, filename=filename, local_dir=local_dir)
            return Path(downloaded).resolve()
        except Exception as exc:  # pragma: no cover - depends on external network/auth
            last_error = exc

    raise FileNotFoundError(
        f"Could not resolve scPRINT-2 checkpoint from {ckpt_arg}. "
        f"Last error: {last_error}"
    )


def pick_embedding_key(adata_emb):
    """Best-effort lookup for the embedding matrix written by scprint2."""
    preferred = ["X_scprint", "X_scprint2", "X_scPRINT", "X_scPRINT2"]
    for key in preferred:
        if key in adata_emb.obsm:
            return key

    for key in adata_emb.obsm.keys():
        if "scprint" in key.lower():
            return key

    for key, value in adata_emb.obsm.items():
        arr = np.asarray(value)
        if arr.ndim == 2 and arr.shape[0] == adata_emb.n_obs and arr.shape[1] > 1:
            return key

    raise KeyError("Could not find a 2D embedding matrix in AnnData.obsm")


def extract_scprint_embeddings(adata, ckpt_path):
    """Run the official scprint2 CLI and read back the embedding matrix."""
    import anndata as ad

    print(f"  Extracting scPRINT-2 embeddings (n={adata.n_obs})...")
    t0 = time.time()

    with tempfile.TemporaryDirectory(prefix="scprint2_embed_", dir=str(OUT_DIR)) as tmpdir:
        tmpdir = Path(tmpdir)
        in_path = tmpdir / "input.h5ad"
        out_name = "embedded.h5ad"
        out_path = tmpdir / out_name
        adata.write_h5ad(in_path)

        cmd = [
            "scprint2", "embed",
            "--adata", str(in_path),
            "--ckpt_path", str(ckpt_path),
            "--species", "NCBITaxon:10090",
            "--output_filename", out_name,
        ]
        proc = subprocess.run(
            cmd,
            cwd=tmpdir,
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                "scprint2 embed failed.\n"
                f"Command: {' '.join(cmd)}\n"
                f"stdout:\n{proc.stdout}\n"
                f"stderr:\n{proc.stderr}"
            )
        if not out_path.exists():
            raise FileNotFoundError(
                f"scprint2 embed completed but did not create {out_path}"
            )

        adata_emb = ad.read_h5ad(out_path)
        key = pick_embedding_key(adata_emb)
        emb = np.asarray(adata_emb.obsm[key], dtype=np.float32)

    print(f"  Embeddings shape: {emb.shape} ({time.time()-t0:.1f}s)")
    return emb


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
                        help="Local checkpoint path or HF-style repo/file reference")
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

    try:
        subprocess.run(
            ["scprint2", "--help"],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            "The scprint2 CLI is not available in this environment. "
            "Run `pip install scprint2` in the active conda env first."
        ) from exc

    # ── Process tissues ───────────────────────────────────────────────────────
    tissues = list(TISSUE_MISSIONS.keys()) if args.tissue == "all" else [args.tissue]

    for tissue in tissues:
        out_file = OUT_DIR / f"SCPRINT2_{tissue}.json"
        if out_file.exists():
            print(f"[SKIP] {out_file.name} already exists")
            continue

        print(f"\n{'='*60}")
        print(f"Tissue: {tissue}")
        print(f"{'='*60}")

        # Load expression data
        expr, y, meta = load_tissue_data(tissue)
        missions = TISSUE_MISSIONS[tissue]

        # Check embedding cache
        cache_dir = Path(args.embedding_cache_dir) if args.embedding_cache_dir else OUT_DIR / "embeddings_cache"
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_file = cache_dir / f"SCPRINT2_emb_{tissue}.npy"

        if cache_file.exists():
            print(f"  Loading cached embeddings: {cache_file.name}")
            embeddings = np.load(cache_file)
            print(f"  Embeddings: {embeddings.shape}")
        else:
            # Convert to AnnData
            adata = expr_to_anndata(expr, y, meta, tissue)

            # Extract scPRINT embeddings
            embeddings = extract_scprint_embeddings(adata, ckpt_path)

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
            "ckpt_path":    str(ckpt_path),
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
            "ckpt_path": str(ckpt_path),
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

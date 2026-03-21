#!/usr/bin/env python3
"""
uce_eval.py — Evaluate UCE (Universal Cell Embeddings) on benchmark tissues.

UCE (Rosen et al., 2024) provides species-agnostic cell embeddings.
We treat each bulk RNA-seq sample as a "pseudo-cell" and extract embeddings,
then use PCA-LR LOMO cross-validation for spaceflight classification.

Requires:
    pip install uce-model anndata scanpy
    Pre-trained model: auto-downloads (4-layer) or manual (33-layer from Figshare)

Input: log2_norm CSVs → reverse to approximate counts → AnnData → UCE → embeddings → PCA-LR

Usage:
    python v3/scripts/uce_eval.py --tissue liver           # Single tissue
    python v3/scripts/uce_eval.py --tissue liver --fast     # Skip bootstrap/perm
    python v3/scripts/uce_eval.py --all                     # All tissues

Output:
    v3/evaluation/FM_uce.json
"""

import argparse
import gc
import json
import os
import shlex
import shutil
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent.parent  # GeneLab_benchmark/
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
DATA_DIR = BASE_DIR / "data" / "mouse"
OUTPUT_DIR = BASE_DIR / "v3" / "evaluation"

# ── Config ───────────────────────────────────────────────────────────────────
FLIGHT_LABELS = {"Flight"}
GROUND_LABELS = {"GC", "VC", "Ground"}
VARIANCE_PERCENTILE = 0.25
N_BOOTSTRAP = 2000
N_PERMUTATIONS = 10000

# Tissues with missions
TISSUE_MISSIONS = {
    "liver": ["RR-1", "RR-3", "RR-6", "RR-8", "RR-9", "MHU-2"],
    "kidney": ["RR-1", "RR-3", "RR-7"],
    "thymus": ["RR-6", "MHU-1", "MHU-2", "RR-9"],
    "gastrocnemius": ["RR-1", "RR-5", "RR-9"],
    "eye": ["RR-1", "RR-3", "TBD"],
    "lung": ["RR-6"],
    "colon": ["RR-6"],
}


# ── Runtime helpers ──────────────────────────────────────────────────────────

def resolve_uce_runner(override=None):
    """Resolve the installed UCE inference entrypoint.

    The Cayuga env exposes `uce-eval-single-anndata` as a console script.
    Older notes referenced `eval_single_anndata`, and module-style invocation
    may exist in some installs. We support all three to make the workflow
    portable across env rebuilds.
    """
    override = override or os.environ.get("UCE_EVAL_CMD")
    checked = []

    if override:
        runner = shlex.split(override)
        exe = runner[0]
        checked.append(override)
        if shutil.which(exe) or Path(exe).exists():
            return runner
        raise FileNotFoundError(
            f"Configured UCE entrypoint not found: {override}"
        )

    for name in ("uce-eval-single-anndata", "eval_single_anndata",
                 "uce_eval_single_anndata"):
        checked.append(name)
        exe = shutil.which(name)
        if exe:
            return [exe]

    import importlib.util

    for module_name in ("eval_single_anndata", "uce_eval_single_anndata"):
        checked.append(f"{sys.executable} -m {module_name}")
        if importlib.util.find_spec(module_name):
            return [sys.executable, "-m", module_name]

    raise RuntimeError(
        "Could not locate a UCE inference entrypoint. "
        f"Checked: {', '.join(checked)}"
    )


# ── Gene symbol mapping ─────────────────────────────────────────────────────

def build_ensmusg_to_symbol(tissue="liver"):
    """Build ENSMUSG → gene symbol mapping.

    Priority:
      1. Pre-cached TSV file (v3/data/ensmusg_to_symbol.tsv)
      2. DGE file from data/mouse/
    """
    # Try pre-cached file first
    cache_path = BASE_DIR / "v3" / "data" / "ensmusg_to_symbol.tsv"
    if cache_path.exists():
        df = pd.read_csv(cache_path, sep="\t")
        mapping = dict(zip(df["ENSEMBL"], df["SYMBOL"]))
        print(f"  Gene mapping: {len(mapping)} genes from {cache_path.name}")
        return mapping

    # Try DGE files
    import glob
    for pattern in [
        str(DATA_DIR / f"{tissue}/*/GLDS-*_rna_seq_differential_expression*.csv"),
        str(DATA_DIR / "*/*/GLDS-*_rna_seq_differential_expression*.csv"),
    ]:
        for f in glob.iglob(pattern):
            try:
                dge = pd.read_csv(f, usecols=["ENSEMBL", "SYMBOL"], nrows=60000)
                dge = dge.dropna(subset=["ENSEMBL", "SYMBOL"])
                mapping = dict(zip(dge["ENSEMBL"], dge["SYMBOL"]))
                if len(mapping) > 1000:
                    print(f"  Gene mapping: {len(mapping)} genes from {Path(f).name}")
                    return mapping
            except Exception:
                continue

    raise RuntimeError("No gene mapping found. Create v3/data/ensmusg_to_symbol.tsv")


# ── Data loading ─────────────────────────────────────────────────────────────

def load_tissue_for_uce(tissue, ensmusg_to_symbol):
    """Load tissue data and prepare for UCE embedding.

    Returns (adata, labels, missions) where:
      - adata: AnnData with approximate raw counts + gene symbols
      - labels: binary series (1=Flight, 0=Ground)
      - missions: series of mission assignments
    """
    import anndata as ad

    tissue_dir = PROCESSED_DIR / tissue

    # Load all_missions files (v1 tissues)
    counts_path = tissue_dir / f"{tissue}_all_missions_log2_norm.csv"
    meta_path = tissue_dir / f"{tissue}_all_missions_metadata.csv"

    if counts_path.exists():
        counts = pd.read_csv(counts_path, index_col=0)
        meta = pd.read_csv(meta_path, index_col=0)
    else:
        # Fallback: load per-mission files and concatenate (v3 tissues)
        import glob
        norm_files = sorted(glob.glob(str(tissue_dir / f"{tissue}_*_log2_norm.csv")))
        meta_files = sorted(glob.glob(str(tissue_dir / f"{tissue}_*_metadata.csv")))
        if not norm_files:
            print(f"  [ERROR] No data files found for {tissue}")
            return None, None, None
        counts_parts, meta_parts = [], []
        for nf, mf in zip(norm_files, meta_files):
            counts_parts.append(pd.read_csv(nf, index_col=0))
            meta_parts.append(pd.read_csv(mf, index_col=0))
        counts = pd.concat(counts_parts, axis=0)
        meta = pd.concat(meta_parts, axis=0)
        # Intersect columns (gene sets may differ across missions)
        common_cols = counts_parts[0].columns
        for part in counts_parts[1:]:
            common_cols = common_cols.intersection(part.columns)
        counts = counts[common_cols]
        print(f"  Loaded {len(norm_files)} per-mission files for {tissue}")

    # Drop non-gene columns
    gene_cols = [c for c in counts.columns if str(c).startswith("ENSMUSG")]
    counts = counts[gene_cols]
    counts = counts.apply(pd.to_numeric, errors="coerce").fillna(0)

    # Filter REMOVE
    if "REMOVE" in meta.columns:
        keep = meta["REMOVE"] != True
        meta = meta[keep]
        counts = counts.loc[counts.index.intersection(meta.index)]

    # Align
    common = counts.index.intersection(meta.index)
    counts = counts.loc[common]
    meta = meta.loc[common]

    # Binary labels
    label_col = "label"
    if label_col not in meta.columns:
        for c in ["label_raw", "condition"]:
            if c in meta.columns:
                label_col = c
                break

    labels = pd.Series(np.nan, index=meta.index)
    labels[meta[label_col].isin(FLIGHT_LABELS)] = 1
    labels[meta[label_col].isin(GROUND_LABELS)] = 0
    valid = ~labels.isna()

    counts = counts[valid]
    meta = meta[valid]
    labels = labels[valid]

    n_flt = int((labels == 1).sum())
    n_gnd = int((labels == 0).sum())
    print(f"  {tissue}: {len(labels)} samples ({n_flt} Flight, {n_gnd} Ground)")

    # Map ENSMUSG → gene symbols
    symbol_cols = []
    ensmusg_cols = []
    for col in counts.columns:
        sym = ensmusg_to_symbol.get(col)
        if sym:
            symbol_cols.append(sym)
            ensmusg_cols.append(col)

    counts_sym = counts[ensmusg_cols].copy()
    counts_sym.columns = symbol_cols

    # Remove duplicate symbols (keep first)
    counts_sym = counts_sym.loc[:, ~counts_sym.columns.duplicated()]
    print(f"  Mapped {len(counts_sym.columns)} genes to symbols")

    # Reverse log2(x+1) → approximate raw counts
    # log2_norm = log2(normalized + 1) → normalized = 2^val - 1
    # UCE expects raw counts, but normalized is close enough for bulk pseudo-cells
    approx_counts = np.power(2, counts_sym.values) - 1
    approx_counts = np.round(approx_counts).astype(np.float32)
    approx_counts = np.clip(approx_counts, 0, None)

    # Create AnnData
    adata = ad.AnnData(
        X=approx_counts,
        obs=pd.DataFrame({"label": labels.values,
                          "mission": meta["mission"].values},
                         index=counts_sym.index),
        var=pd.DataFrame(index=counts_sym.columns)
    )
    adata.var_names_make_unique()

    missions = meta["mission"]
    return adata, labels, missions


# ── UCE embedding extraction ────────────────────────────────────────────────

def extract_uce_embeddings(adata, model_dir=None, batch_size=25, nlayers=4,
                           uce_cli=None):
    """Run UCE to extract cell embeddings.

    Args:
        adata: AnnData with raw counts and gene symbols in var_names
        model_dir: Path to save/load model (None = auto)
        batch_size: Batch size for inference
        nlayers: 4 (default, lighter) or 33 (paper model)

    Returns:
        embeddings: np.ndarray (n_samples, embedding_dim)
    """
    import subprocess
    import tempfile

    tmp_dir = Path(tempfile.mkdtemp(prefix="uce_eval_"))
    runner = resolve_uce_runner(uce_cli)

    try:
        # Save AnnData to temp file
        adata_path = tmp_dir / "input.h5ad"
        adata.write_h5ad(adata_path)

        output_dir = tmp_dir / "output"
        output_dir.mkdir()
        output_dir_arg = str(output_dir)
        if not output_dir_arg.endswith(os.sep):
            output_dir_arg += os.sep

        cmd = runner + [
        "--adata_path", str(adata_path),
        "--dir", output_dir_arg,
        "--species", "mouse",
        "--batch_size", str(batch_size),
        "--nlayers", str(nlayers),
        ]

        if model_dir:
            cmd.extend(["--model_loc", str(model_dir)])

        print(f"  Running UCE via: {' '.join(cmd[:3])} ...")
        print(f"  UCE config: nlayers={nlayers}, batch={batch_size}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)

        if result.returncode != 0:
            stderr_tail = (result.stderr or "").strip()[-1000:]
            stdout_tail = (result.stdout or "").strip()[-1000:]
            raise RuntimeError(
                "UCE failed. "
                f"stderr tail: {stderr_tail or '[empty]'} "
                f"stdout tail: {stdout_tail or '[empty]'}"
            )

        # Load embeddings from output.
        #
        # Older UCE installs concatenate `args.dir + filename` instead of using
        # a path join, so without a trailing slash the file may be written next
        # to the directory as e.g. `outputinput_uce_adata.h5ad`. We pass the
        # corrected dir string above and still search the full temp workspace to
        # stay compatible with both behaviors.
        import anndata as ad
        expected_name = f"{adata_path.stem}_uce_adata.h5ad"
        candidate_paths = []
        for candidate in (
            output_dir / expected_name,
            tmp_dir / f"output{expected_name}",
        ):
            if candidate.exists():
                candidate_paths.append(candidate)

        for pattern in ("*_uce_adata.h5ad", "*.h5ad"):
            for path in sorted(tmp_dir.rglob(pattern)):
                if path == adata_path or path in candidate_paths:
                    continue
                candidate_paths.append(path)

        if not candidate_paths:
            raise RuntimeError(f"No output h5ad found under {tmp_dir}")

        output_path = candidate_paths[0]
        print(f"  UCE output file: {output_path}")
        adata_out = ad.read_h5ad(output_path)
        if "X_uce" in adata_out.obsm:
            embeddings = adata_out.obsm["X_uce"]
        else:
            for key in adata_out.obsm:
                if "uce" in key.lower() or "embed" in key.lower():
                    embeddings = adata_out.obsm[key]
                    break
            else:
                raise RuntimeError(
                    f"No UCE embeddings found. obsm keys: {list(adata_out.obsm)}"
                )

        print(f"  UCE embeddings: {embeddings.shape}")
        return embeddings
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# ── LOMO Classification ──────────────────────────────────────────────────────

def lomo_cv(embeddings, labels, missions, fast=False):
    """Leave-One-Mission-Out cross-validation with PCA-LR.

    Returns dict with overall AUROC + per-mission breakdown.
    """
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler

    labels_arr = labels.values.astype(int)
    missions_arr = missions.values
    unique_missions = sorted(set(missions_arr))

    if len(unique_missions) < 2:
        print("  [WARN] Only 1 mission — using 5-fold CV instead of LOMO")
        return _kfold_cv(embeddings, labels_arr, fast=fast)

    all_y_true = []
    all_y_score = []
    per_mission = {}

    for held_out in unique_missions:
        test_mask = missions_arr == held_out
        train_mask = ~test_mask

        X_train = embeddings[train_mask]
        X_test = embeddings[test_mask]
        y_train = labels_arr[train_mask]
        y_test = labels_arr[test_mask]

        if len(np.unique(y_train)) < 2 or len(np.unique(y_test)) < 2:
            continue

        # StandardScaler + PCA + LR
        scaler = StandardScaler()
        X_tr = scaler.fit_transform(X_train)
        X_te = scaler.transform(X_test)

        n_comps = min(20, X_tr.shape[0] - 1, X_tr.shape[1])
        if n_comps < 2:
            continue

        pca = PCA(n_components=n_comps, random_state=42)
        X_tr = pca.fit_transform(X_tr)
        X_te = pca.transform(X_te)

        clf = LogisticRegression(C=1.0, class_weight="balanced",
                                 max_iter=1000, random_state=42)
        clf.fit(X_tr, y_train)
        y_score = clf.predict_proba(X_te)[:, 1]

        all_y_true.extend(y_test)
        all_y_score.extend(y_score)

        try:
            fold_auroc = float(roc_auc_score(y_test, y_score))
        except Exception:
            fold_auroc = None

        per_mission[held_out] = {
            "n_samples": int(test_mask.sum()),
            "auroc": fold_auroc,
        }

    if not all_y_true:
        return {"lomo_auroc": None, "per_mission": per_mission}

    y_true = np.array(all_y_true)
    y_score = np.array(all_y_score)
    overall_auroc = float(roc_auc_score(y_true, y_score))

    # Bootstrap CI
    ci_low, ci_high = np.nan, np.nan
    if not fast:
        rng = np.random.RandomState(42)
        boots = []
        for _ in range(N_BOOTSTRAP):
            idx = rng.randint(0, len(y_true), size=len(y_true))
            yt, ys = y_true[idx], y_score[idx]
            if len(np.unique(yt)) < 2:
                continue
            try:
                boots.append(roc_auc_score(yt, ys))
            except Exception:
                continue
        if boots:
            ci_low = float(np.percentile(boots, 2.5))
            ci_high = float(np.percentile(boots, 97.5))

    # Permutation test
    perm_p = np.nan
    if not fast:
        rng = np.random.RandomState(42)
        count = 0
        for _ in range(N_PERMUTATIONS):
            perm_y = rng.permutation(y_true)
            if len(np.unique(perm_y)) < 2:
                continue
            try:
                if roc_auc_score(perm_y, y_score) >= overall_auroc:
                    count += 1
            except Exception:
                continue
        perm_p = (count + 1) / (N_PERMUTATIONS + 1)

    return {
        "lomo_auroc": overall_auroc,
        "ci_lower": None if np.isnan(ci_low) else ci_low,
        "ci_upper": None if np.isnan(ci_high) else ci_high,
        "perm_p": None if np.isnan(perm_p) else perm_p,
        "n_folds": len(unique_missions),
        "per_mission": per_mission,
    }


def _kfold_cv(embeddings, labels_arr, k=5, fast=False):
    """Fallback: stratified k-fold CV for single-mission tissues."""
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler

    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=42)
    all_y_true, all_y_score = [], []

    for train_idx, test_idx in skf.split(embeddings, labels_arr):
        scaler = StandardScaler()
        X_tr = scaler.fit_transform(embeddings[train_idx])
        X_te = scaler.transform(embeddings[test_idx])

        n_comps = min(20, X_tr.shape[0] - 1, X_tr.shape[1])
        if n_comps < 2:
            continue

        pca = PCA(n_components=n_comps, random_state=42)
        X_tr = pca.fit_transform(X_tr)
        X_te = pca.transform(X_te)

        clf = LogisticRegression(C=1.0, class_weight="balanced",
                                 max_iter=1000, random_state=42)
        clf.fit(X_tr, labels_arr[train_idx])
        y_score = clf.predict_proba(X_te)[:, 1]
        all_y_true.extend(labels_arr[test_idx])
        all_y_score.extend(y_score)

    if not all_y_true:
        return {"lomo_auroc": None, "n_folds": k, "cv_type": "stratified_kfold",
                "per_mission": {}}

    y_true = np.array(all_y_true)
    y_score = np.array(all_y_score)
    auroc = float(roc_auc_score(y_true, y_score))

    return {"lomo_auroc": auroc, "n_folds": k, "cv_type": "stratified_kfold",
            "per_mission": {}}


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    global N_BOOTSTRAP, N_PERMUTATIONS

    parser = argparse.ArgumentParser(description="UCE FM evaluation")
    parser.add_argument("--tissue", help="Evaluate single tissue")
    parser.add_argument("--all", action="store_true", help="All tissues")
    parser.add_argument("--fast", action="store_true", help="Skip bootstrap/perm")
    parser.add_argument("--batch-size", type=int, default=25)
    parser.add_argument("--nlayers", type=int, default=4, choices=[4, 33])
    parser.add_argument("--model-dir", help="Path to UCE model weights")
    parser.add_argument("--uce-cli",
                        help="Path to the installed UCE CLI entrypoint")
    parser.add_argument("--embeddings-only", action="store_true",
                        help="Only extract embeddings, skip classification")
    args = parser.parse_args()

    if args.fast:
        N_BOOTSTRAP = 0
        N_PERMUTATIONS = 0

    print("=" * 70)
    print("GeneLabBench v3 — UCE Foundation Model Evaluation")
    print(f"Timestamp: {datetime.now().isoformat()}")
    print(f"Model: {args.nlayers}-layer UCE")
    print("=" * 70)

    # Build gene mapping
    print("\nBuilding ENSMUSG → symbol mapping...")
    ensmusg_to_symbol = build_ensmusg_to_symbol()

    # Determine tissues
    if args.tissue:
        tissues = [args.tissue]
    elif args.all:
        tissues = list(TISSUE_MISSIONS.keys())
    else:
        tissues = ["liver"]  # Default: primary benchmark tissue

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    all_results = []

    for tissue in tissues:
        print(f"\n{'='*60}")
        print(f"  Tissue: {tissue}")
        print(f"{'='*60}")

        # Load data
        adata, labels, missions = load_tissue_for_uce(tissue, ensmusg_to_symbol)
        if adata is None:
            continue

        # Extract UCE embeddings
        try:
            embeddings = extract_uce_embeddings(
                adata, model_dir=args.model_dir,
                batch_size=args.batch_size, nlayers=args.nlayers,
                uce_cli=args.uce_cli)
        except Exception as e:
            print(f"  [ERROR] UCE embedding extraction failed: {e}")
            # Save partial result
            all_results.append({
                "model": f"UCE-{args.nlayers}L",
                "tissue": tissue,
                "error": str(e),
            })
            continue

        if args.embeddings_only:
            # Save embeddings
            emb_path = OUTPUT_DIR / f"UCE_{tissue}_embeddings.npy"
            np.save(emb_path, embeddings)
            print(f"  Saved embeddings: {emb_path}")
            continue

        # LOMO classification
        print(f"\n  LOMO Cross-Validation:")
        result = lomo_cv(embeddings, labels, missions, fast=args.fast)

        tissue_result = {
            "model": f"UCE-{args.nlayers}L",
            "tissue": tissue,
            "n_samples": len(labels),
            "n_missions": len(set(missions)),
            "embedding_dim": embeddings.shape[1],
            "n_genes_used": adata.shape[1],
            **result,
        }
        all_results.append(tissue_result)

        auroc = result.get("lomo_auroc")
        auroc_str = f"{auroc:.3f}" if auroc is not None else "N/A"
        print(f"  AUROC = {auroc_str}")
        if result.get("ci_lower") is not None:
            print(f"  95% CI = [{result['ci_lower']:.3f}, {result['ci_upper']:.3f}]")
        if result.get("perm_p") is not None:
            print(f"  Perm p = {result['perm_p']:.4f}")

        gc.collect()

    # Save results
    output = {
        "model_family": "UCE",
        "nlayers": args.nlayers,
        "results": all_results,
        "config": {
            "variance_percentile": VARIANCE_PERCENTILE,
            "n_bootstrap": N_BOOTSTRAP,
            "n_permutations": N_PERMUTATIONS,
            "batch_size": args.batch_size,
        },
        "timestamp": datetime.now().isoformat(),
    }

    out_path = OUTPUT_DIR / "FM_uce.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else str(o))
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()

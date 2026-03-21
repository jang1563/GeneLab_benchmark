#!/usr/bin/env python3
"""
scfoundation_eval.py — Evaluate scFoundation on benchmark tissues.

scFoundation (Hao et al., 2024 Nature Methods) is a 100M-param gene FM
trained on 50M+ human scRNA-seq. Supports bulk RNA-seq via --input_type bulk.

Requires human gene symbols (mouse → human ortholog mapping needed).

Pipeline:
  1. Load log2_norm data → map ENSMUSG to human gene symbols
  2. Align to scFoundation's 19,264-gene vocabulary
  3. Run scFoundation get_embedding.py → cell embeddings
  4. PCA-LR LOMO cross-validation → AUROC

Requires:
    git clone https://github.com/biomap-research/scFoundation
    Download model weights to scFoundation/models/
    Mouse→human ortholog mapping file

Usage:
    python v3/scripts/scfoundation_eval.py --tissue liver --scf-dir /path/to/scFoundation
    python v3/scripts/scfoundation_eval.py --tissue liver --fast

Output:
    v3/evaluation/FM_scfoundation.json
"""

import argparse
import gc
import json
import shutil
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

# ── Paths ────────────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).resolve().parent.parent.parent
PROCESSED_DIR = BASE_DIR / "processed" / "A_detection"
DATA_DIR = BASE_DIR / "data" / "mouse"
OUTPUT_DIR = BASE_DIR / "v3" / "evaluation"

# ── Config ───────────────────────────────────────────────────────────────────
FLIGHT_LABELS = {"Flight"}
GROUND_LABELS = {"GC", "VC", "Ground"}
N_BOOTSTRAP = 2000
N_PERMUTATIONS = 10000

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

def resolve_scfoundation_layout(scf_dir):
    """Resolve the runnable scFoundation layout.

    The upstream repo may be passed either as the repo root or as the nested
    `model/` directory. Cayuga currently has `get_embedding.py` and the model
    weight folder under `scFoundation/model/`.
    """
    scf_dir = Path(scf_dir).expanduser().resolve()

    candidate_runtime_dirs = []
    if scf_dir.name == "model":
        candidate_runtime_dirs.append(scf_dir)
        candidate_runtime_dirs.append(scf_dir.parent)
    else:
        candidate_runtime_dirs.append(scf_dir / "model")
        candidate_runtime_dirs.append(scf_dir)

    seen = set()
    runtime_dirs = []
    for candidate in candidate_runtime_dirs:
        if candidate in seen:
            continue
        seen.add(candidate)
        runtime_dirs.append(candidate)

    for runtime_dir in runtime_dirs:
        get_embedding = runtime_dir / "get_embedding.py"
        if not get_embedding.exists():
            continue

        repo_root = runtime_dir.parent if runtime_dir.name == "model" else runtime_dir

        gene_index_candidates = [
            runtime_dir / "OS_scRNA_gene_index.19264.tsv",
            repo_root / "OS_scRNA_gene_index.19264.tsv",
        ]
        model_ckpt_candidates = [
            runtime_dir / "models" / "models.ckpt",
            repo_root / "models" / "models.ckpt",
        ]
        download_hint_candidates = [
            runtime_dir / "models" / "download.txt",
            repo_root / "models" / "download.txt",
        ]

        gene_index = next((p for p in gene_index_candidates if p.exists()),
                          gene_index_candidates[0])
        model_ckpt = next((p for p in model_ckpt_candidates if p.exists()),
                          model_ckpt_candidates[0])
        download_hint = next((p for p in download_hint_candidates if p.exists()),
                             download_hint_candidates[0])

        return {
            "repo_root": repo_root,
            "runtime_dir": runtime_dir,
            "get_embedding": get_embedding,
            "gene_index": gene_index,
            "model_ckpt": model_ckpt,
            "download_hint": download_hint,
        }

    raise FileNotFoundError(
        f"Could not find get_embedding.py under {scf_dir} or {scf_dir / 'model'}"
    )


# ── Ortholog mapping ─────────────────────────────────────────────────────────

def load_mouse_to_human_orthologs(ortholog_file=None):
    """Load mouse→human ortholog mapping.

    Tries (in order):
      1. Provided ortholog_file path
      2. Pre-cached file in v3/data/
      3. BioMart query (requires internet)

    Returns dict: {mouse_symbol: human_symbol}
    """
    # Try pre-cached file
    if ortholog_file is None:
        candidates = [
            BASE_DIR / "v3" / "data" / "mouse_human_orthologs.tsv",
            BASE_DIR / "processed" / "gene_sets" / "mouse_human_orthologs.tsv",
        ]
        for f in candidates:
            if f.exists():
                ortholog_file = f
                break

    if ortholog_file and Path(ortholog_file).exists():
        df = pd.read_csv(ortholog_file, sep="\t")
        # Expected columns: mouse_symbol, human_symbol (or similar)
        mouse_col = [c for c in df.columns if "mouse" in c.lower() or "mmus" in c.lower()]
        human_col = [c for c in df.columns if "human" in c.lower() or "hsap" in c.lower()]
        if mouse_col and human_col:
            mapping = dict(zip(df[mouse_col[0]], df[human_col[0]]))
            print(f"  Ortholog mapping: {len(mapping)} pairs from {Path(ortholog_file).name}")
            return mapping

    # Fallback: build from DGE files (SYMBOL is already mouse MGI)
    # Need to create a mouse→human mapping
    print("  [WARN] No ortholog file found. Attempting BioMart query...")
    try:
        return _query_biomart_orthologs()
    except Exception as e:
        print(f"  [ERROR] BioMart query failed: {e}")
        print("  Create ortholog file: v3/data/mouse_human_orthologs.tsv")
        print("  Columns: mouse_symbol\\thuman_symbol")
        raise


def _query_biomart_orthologs():
    """Query BioMart for mouse→human 1:1 orthologs."""
    import urllib.request

    url = ("http://www.ensembl.org/biomart/martservice?"
           "query=<?xml version='1.0' encoding='UTF-8'?>"
           "<!DOCTYPE Query>"
           "<Query virtualSchemaName='default' formatter='TSV' header='1'>"
           "<Dataset name='mmusculus_gene_ensembl' interface='default'>"
           "<Attribute name='external_gene_name'/>"
           "<Attribute name='hsapiens_homolog_associated_gene_name'/>"
           "<Attribute name='hsapiens_homolog_orthology_type'/>"
           "</Dataset></Query>")

    print("  Querying Ensembl BioMart for orthologs...")
    response = urllib.request.urlopen(url, timeout=120)
    data = response.read().decode("utf-8")

    lines = data.strip().split("\n")
    mapping = {}
    for line in lines[1:]:  # skip header
        parts = line.split("\t")
        if len(parts) >= 3 and parts[2] == "ortholog_one2one":
            mouse_sym, human_sym = parts[0].strip(), parts[1].strip()
            if mouse_sym and human_sym:
                mapping[mouse_sym] = human_sym

    if len(mapping) < 1000:
        raise RuntimeError(f"Too few orthologs: {len(mapping)}")

    # Cache
    cache_dir = BASE_DIR / "v3" / "data"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / "mouse_human_orthologs.tsv"
    with open(cache_path, "w") as f:
        f.write("mouse_symbol\thuman_symbol\n")
        for m, h in sorted(mapping.items()):
            f.write(f"{m}\t{h}\n")
    print(f"  Cached {len(mapping)} orthologs to {cache_path}")

    return mapping


# ── Data loading ─────────────────────────────────────────────────────────────

def build_ensmusg_to_symbol(tissue="liver"):
    """Build ENSMUSG → mouse gene symbol mapping."""
    cache_path = BASE_DIR / "v3" / "data" / "ensmusg_to_symbol.tsv"
    if cache_path.exists():
        df = pd.read_csv(cache_path, sep="\t")
        if {"ENSEMBL", "SYMBOL"}.issubset(df.columns):
            return dict(zip(df["ENSEMBL"], df["SYMBOL"]))

    import glob
    for f in glob.iglob(str(DATA_DIR / f"{tissue}/*/GLDS-*_rna_seq_differential_expression*.csv")):
        try:
            dge = pd.read_csv(f, usecols=["ENSEMBL", "SYMBOL"], nrows=60000)
            dge = dge.dropna(subset=["ENSEMBL", "SYMBOL"])
            mapping = dict(zip(dge["ENSEMBL"], dge["SYMBOL"]))
            if len(mapping) > 1000:
                return mapping
        except Exception:
            continue

    # Fallback: any tissue
    for f in glob.iglob(str(DATA_DIR / "*/*/GLDS-*_rna_seq_differential_expression*.csv")):
        try:
            dge = pd.read_csv(f, usecols=["ENSEMBL", "SYMBOL"], nrows=60000)
            dge = dge.dropna(subset=["ENSEMBL", "SYMBOL"])
            mapping = dict(zip(dge["ENSEMBL"], dge["SYMBOL"]))
            if len(mapping) > 1000:
                return mapping
        except Exception:
            continue
    raise RuntimeError("No DGE file found for symbol mapping")


def load_tissue_for_scfoundation(tissue, ensmusg_to_symbol, mouse_to_human,
                                  scf_gene_list):
    """Load tissue data and prepare for scFoundation.

    Returns (expression_df, labels, missions) where:
      - expression_df: DataFrame aligned to scFoundation's 19264-gene vocab
      - labels: binary series
      - missions: series
    """
    tissue_dir = PROCESSED_DIR / tissue
    counts_path = tissue_dir / f"{tissue}_all_missions_log2_norm.csv"
    meta_path = tissue_dir / f"{tissue}_all_missions_metadata.csv"

    if counts_path.exists():
        counts = pd.read_csv(counts_path, index_col=0)
        meta = pd.read_csv(meta_path, index_col=0)
    else:
        # Fallback: load per-mission files and concatenate (v3 tissues)
        import glob as _glob
        norm_files = sorted(_glob.glob(str(tissue_dir / f"{tissue}_*_log2_norm.csv")))
        meta_files = sorted(_glob.glob(str(tissue_dir / f"{tissue}_*_metadata.csv")))
        if not norm_files:
            return None, None, None
        counts_parts, meta_parts = [], []
        for nf, mf in zip(norm_files, meta_files):
            counts_parts.append(pd.read_csv(nf, index_col=0))
            meta_parts.append(pd.read_csv(mf, index_col=0))
        counts = pd.concat(counts_parts, axis=0)
        meta = pd.concat(meta_parts, axis=0)
        common_cols = counts_parts[0].columns
        for part in counts_parts[1:]:
            common_cols = common_cols.intersection(part.columns)
        counts = counts[common_cols]
        print(f"  Loaded {len(norm_files)} per-mission files for {tissue}")

    gene_cols = [c for c in counts.columns if str(c).startswith("ENSMUSG")]
    counts = counts[gene_cols].apply(pd.to_numeric, errors="coerce").fillna(0)

    if "REMOVE" in meta.columns:
        keep = meta["REMOVE"] != True
        meta = meta[keep]
        counts = counts.loc[counts.index.intersection(meta.index)]

    common = counts.index.intersection(meta.index)
    counts = counts.loc[common]
    meta = meta.loc[common]

    label_col = "label"
    if label_col not in meta.columns:
        for candidate in ("label_raw", "condition", "group"):
            if candidate in meta.columns:
                label_col = candidate
                break

    labels = pd.Series(np.nan, index=meta.index)
    labels[meta[label_col].isin(FLIGHT_LABELS)] = 1
    labels[meta[label_col].isin(GROUND_LABELS)] = 0
    valid = ~labels.isna()
    counts = counts[valid]
    meta = meta[valid]
    labels = labels[valid]

    print(f"  {tissue}: {len(labels)} samples "
          f"({int((labels==1).sum())} Flt, {int((labels==0).sum())} Gnd)")

    # Map: ENSMUSG → mouse symbol → human symbol
    human_cols = {}
    for ens_col in counts.columns:
        mouse_sym = ensmusg_to_symbol.get(ens_col)
        if mouse_sym:
            human_sym = mouse_to_human.get(mouse_sym)
            if human_sym and human_sym in scf_gene_list:
                human_cols[ens_col] = human_sym

    print(f"  Mapped {len(human_cols)}/{len(counts.columns)} genes to human orthologs "
          f"in scFoundation vocabulary")

    # Build expression matrix aligned to scFoundation vocabulary
    # Data is log2(x+1). scFoundation with --pre_normalized T expects log1p(natural log).
    # Convert: ln(x+1) = log2(x+1) * ln(2)
    expr_human = pd.DataFrame(0.0, index=counts.index, columns=scf_gene_list)
    for ens_col, human_sym in human_cols.items():
        # Convert log2(x+1) to ln(x+1) = log2(x+1) * ln(2)
        expr_human[human_sym] = counts[ens_col].values * np.log(2)

    missions = meta["mission"]
    return expr_human, labels, missions


# ── scFoundation embedding extraction ────────────────────────────────────────

def extract_scfoundation_embeddings(expr_df, scf_dir, pre_normalized="T"):
    """Run scFoundation get_embedding.py to extract cell embeddings.

    Args:
        expr_df: DataFrame (samples × 19264 genes), ln(x+1) normalized
        scf_dir: Path to scFoundation repo
        pre_normalized: "T" (already log1p), "F" (raw counts)

    Returns:
        embeddings: np.ndarray (n_samples, embedding_dim)
    """
    import subprocess
    import tempfile

    layout = resolve_scfoundation_layout(scf_dir)
    get_embedding = layout["get_embedding"]
    runtime_dir = layout["runtime_dir"]
    model_ckpt = layout["model_ckpt"]
    download_hint = layout["download_hint"]

    if not model_ckpt.exists():
        hint = ""
        if download_hint.exists():
            hint = download_hint.read_text().strip()
        raise FileNotFoundError(
            f"scFoundation weights not found: {model_ckpt}. "
            f"{hint}".strip()
        )

    tmp_dir = Path(tempfile.mkdtemp(prefix="scfoundation_eval_"))

    try:
        # Save expression to CSV
        data_path = tmp_dir / "bulk_expression.csv"
        expr_df.to_csv(data_path)

        save_path = tmp_dir / "output"
        save_path.mkdir()

        cmd = [
        sys.executable, str(get_embedding),
        "--task_name", "benchmark",
        "--input_type", "bulk",
        "--output_type", "cell",
        "--pool_type", "all",
        "--tgthighres", "f1",
        "--data_path", str(data_path),
        "--save_path", str(save_path),
        "--pre_normalized", pre_normalized,
        "--version", "rde",
        ]

        print(f"  Running scFoundation from: {runtime_dir}")
        result = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=3600, cwd=str(runtime_dir))

        if result.returncode != 0:
            stderr_tail = (result.stderr or "").strip()[-1000:]
            stdout_tail = (result.stdout or "").strip()[-1000:]
            raise RuntimeError(
                "scFoundation failed. "
                f"stderr tail: {stderr_tail or '[empty]'} "
                f"stdout tail: {stdout_tail or '[empty]'}"
            )

        emb_files = list(save_path.glob("*.npy"))
        if not emb_files:
            raise RuntimeError(f"No embedding files in {save_path}")

        embeddings = np.load(emb_files[0])
        print(f"  scFoundation embeddings: {embeddings.shape}")
        return embeddings
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# ── LOMO Classification (same as UCE) ────────────────────────────────────────

def lomo_cv(embeddings, labels, missions, fast=False):
    """Leave-One-Mission-Out CV with PCA-LR."""
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

    all_y_true, all_y_score = [], []
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
        per_mission[held_out] = {"n_samples": int(test_mask.sum()),
                                 "auroc": fold_auroc}

    if not all_y_true:
        return {"lomo_auroc": None, "per_mission": per_mission}

    y_true, y_score = np.array(all_y_true), np.array(all_y_score)
    overall = float(roc_auc_score(y_true, y_score))

    ci_low = ci_high = perm_p = None
    if not fast:
        rng = np.random.RandomState(42)
        boots = []
        for _ in range(N_BOOTSTRAP):
            idx = rng.randint(0, len(y_true), size=len(y_true))
            yt, ys = y_true[idx], y_score[idx]
            if len(np.unique(yt)) >= 2:
                try:
                    boots.append(roc_auc_score(yt, ys))
                except Exception:
                    pass
        if boots:
            ci_low = float(np.percentile(boots, 2.5))
            ci_high = float(np.percentile(boots, 97.5))

        count = 0
        for _ in range(N_PERMUTATIONS):
            perm_y = rng.permutation(y_true)
            if len(np.unique(perm_y)) >= 2:
                try:
                    if roc_auc_score(perm_y, y_score) >= overall:
                        count += 1
                except Exception:
                    pass
        perm_p = (count + 1) / (N_PERMUTATIONS + 1)

    return {
        "lomo_auroc": overall, "ci_lower": ci_low, "ci_upper": ci_high,
        "perm_p": perm_p, "n_folds": len(unique_missions),
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

    parser = argparse.ArgumentParser(description="scFoundation FM evaluation")
    parser.add_argument("--tissue", default="liver")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--fast", action="store_true")
    parser.add_argument("--scf-dir", required=True,
                        help="Path to scFoundation repo")
    parser.add_argument("--ortholog-file",
                        help="Mouse→human ortholog TSV")
    args = parser.parse_args()

    if args.fast:
        N_BOOTSTRAP = 0
        N_PERMUTATIONS = 0

    scf_dir = Path(args.scf_dir)
    layout = resolve_scfoundation_layout(scf_dir)

    print("=" * 70)
    print("GeneLabBench v3 — scFoundation FM Evaluation")
    print(f"Timestamp: {datetime.now().isoformat()}")
    print(f"scFoundation dir: {layout['repo_root']}")
    print(f"scFoundation runtime dir: {layout['runtime_dir']}")
    print("=" * 70)

    # Load scFoundation gene vocabulary
    gene_index_path = layout["gene_index"]
    if not gene_index_path.exists():
        raise FileNotFoundError(f"Gene index not found: {gene_index_path}")
    gene_df = pd.read_csv(gene_index_path, sep="\t")
    gene_col = "gene_name" if "gene_name" in gene_df.columns else gene_df.columns[0]
    scf_gene_list = list(gene_df[gene_col])
    print(f"scFoundation vocabulary: {len(scf_gene_list)} genes")

    # Build mappings
    print("\nBuilding gene mappings...")
    ensmusg_to_symbol = build_ensmusg_to_symbol()
    mouse_to_human = load_mouse_to_human_orthologs(args.ortholog_file)
    print(f"  ENSMUSG→symbol: {len(ensmusg_to_symbol)}, mouse→human: {len(mouse_to_human)}")

    tissues = list(TISSUE_MISSIONS.keys()) if args.all else [args.tissue]
    all_results = []

    for tissue in tissues:
        print(f"\n{'='*60}")
        print(f"  Tissue: {tissue}")
        print(f"{'='*60}")

        expr_df, labels, missions = load_tissue_for_scfoundation(
            tissue, ensmusg_to_symbol, mouse_to_human, scf_gene_list)
        if expr_df is None:
            continue

        try:
            embeddings = extract_scfoundation_embeddings(expr_df, scf_dir)
        except Exception as e:
            print(f"  [ERROR] Embedding extraction failed: {e}")
            all_results.append({"model": "scFoundation", "tissue": tissue,
                               "error": str(e)})
            continue

        result = lomo_cv(embeddings, labels, missions, fast=args.fast)
        tissue_result = {
            "model": "scFoundation",
            "tissue": tissue,
            "n_samples": len(labels),
            "n_missions": len(set(missions)),
            "embedding_dim": embeddings.shape[1],
            "n_genes_mapped": int((expr_df != 0).any(axis=0).sum()),
            **result,
        }
        all_results.append(tissue_result)

        auroc = result.get("lomo_auroc")
        auroc_str = f"{auroc:.3f}" if auroc is not None else "N/A"
        print(f"  AUROC = {auroc_str}")
        gc.collect()

    # Save
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output = {
        "model_family": "scFoundation",
        "results": all_results,
        "timestamp": datetime.now().isoformat(),
    }
    out_path = OUTPUT_DIR / "FM_scfoundation.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2,
                  default=lambda o: float(o) if isinstance(o, np.floating) else str(o))
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()

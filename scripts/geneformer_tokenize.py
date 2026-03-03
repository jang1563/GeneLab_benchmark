"""
Geneformer Tokenization Pipeline for Mouse Bulk RNA-seq (NASA GeneLab A4 Thymus)

This script implements rank-value encoding tokenization for Geneformer,
adapted for mouse bulk RNA-seq data (ENSMUSG IDs → ENSG via ortholog mapping).

Pipeline:
  1. Download Ensembl mouse→human ortholog table (one-to-one, cached)
  2. Load Geneformer V1 token dictionary + gene median expression
  3. For each sample: expression → ortholog map → median scale → rank → token sequence
  4. Save as HuggingFace Arrow Dataset for downstream fine-tuning

Usage:
  # Tokenize all A4 LOMO folds (including held-out)
  python scripts/geneformer_tokenize.py --task A4 --model-version v1

  # Tokenize single fold (dry-run/test)
  python scripts/geneformer_tokenize.py --task A4 --fold RR-9 --dry-run

Output:
  tasks/A4_thymus_lomo/fold_{fold}/geneformer_tokens/
    train_dataset/    # HuggingFace Arrow dataset
    test_dataset/     # HuggingFace Arrow dataset
    tokenize_info.json

References:
  Theodoris et al., Nature 2023: https://doi.org/10.1038/s41586-023-06139-9
  Geneformer repo: https://huggingface.co/ctheodoris/Geneformer
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import pickle
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ─── Paths ─────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).parent.parent
DATA_DIR = ROOT / "data" / "mouse"
ORTHOLOG_CACHE = DATA_DIR / "ensembl_mouse_human_orthologs.tsv"
TASKS_DIR = ROOT / "tasks"

# Geneformer HuggingFace model ID
HF_REPO = "ctheodoris/Geneformer"

# Mouse-Geneformer local paths (Cayuga)
MOUSE_GF_BASE = Path(
    os.environ.get(
        "MOUSE_GF_BASE",
        "/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/models",
    )
)

# Model version configs
MODEL_CONFIGS = {
    "v1": {
        "model_path": "Geneformer-V1-10M",
        "token_dict_file": "geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl",
        "gene_median_file": "geneformer/gene_dictionaries_30m/gene_median_dictionary_gc30M.pkl",
        "max_length": 2048,
        "special_token": False,  # V1: no <cls>/<eos>
    },
    "v2": {
        "model_path": "Geneformer-V2-104M",
        "token_dict_file": "geneformer/token_dictionary_gc104M.pkl",
        "gene_median_file": "geneformer/gene_median_dictionary_gc104M.pkl",
        "max_length": 4096,
        "special_token": True,  # V2: add <cls>/<eos>
    },
    "mouse_gf": {
        "model_path": str(MOUSE_GF_BASE / "mouse_gf_base"),
        "token_dict_file": str(MOUSE_GF_BASE / "MLM-re_token_dictionary_v1.pkl"),
        "gene_median_file": str(MOUSE_GF_BASE / "mouse_gene_median_dictionary.pkl"),
        "max_length": 2048,
        "special_token": False,  # No special tokens (same as V1)
    },
}

def resolve_task_dir(task: str, task_dir_name: Optional[str] = None) -> Path:
    """Resolve one task directory deterministically."""
    if task_dir_name:
        task_dir = TASKS_DIR / task_dir_name
        if not task_dir.exists() or not task_dir.is_dir():
            raise FileNotFoundError(f"Task directory not found: {task_dir}")
        if not task_dir.name.startswith(f"{task}_"):
            raise ValueError(
                f"--task-dir '{task_dir.name}' does not match task '{task}'"
            )
        return task_dir

    matches = sorted(TASKS_DIR.glob(f"{task}_*"))
    if not matches:
        raise FileNotFoundError(f"No task directory found for '{task}' in {TASKS_DIR}")
    if len(matches) > 1:
        names = ", ".join(d.name for d in matches)
        raise ValueError(
            f"Ambiguous task '{task}': {names}. Use --task-dir to select one."
        )
    return matches[0]


# ─── Step 1: Ortholog Mapping ───────────────────────────────────────────────────

def download_ortholog_table(cache_path: Path, timeout: int = 300) -> pd.DataFrame:
    """
    Download mouse→human one-to-one orthologs from Ensembl BioMart.
    Returns DataFrame with columns: ensembl_gene_id, hsapiens_ensembl, pct_id
    """
    if cache_path.exists():
        log.info(f"Loading cached ortholog table from {cache_path}")
        df = pd.read_csv(cache_path, sep="\t")
        log.info(f"  Cached: {len(df):,} mouse genes, {df['hsapiens_ensembl'].notna().sum():,} with human orthologs")
        return df

    log.info("Downloading mouse→human orthologs from Ensembl BioMart...")
    query_xml = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
    <Dataset name="mmusculus_gene_ensembl" interface="default">
        <Attribute name="ensembl_gene_id"/>
        <Attribute name="hsapiens_homolog_ensembl_gene"/>
        <Attribute name="hsapiens_homolog_perc_id"/>
    </Dataset>
</Query>"""

    url = "https://www.ensembl.org/biomart/martservice?query=" + urllib.parse.quote(query_xml)
    req = urllib.request.Request(url, headers={"User-Agent": "Python/GeneLab-Benchmark"})

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.time()

    with urllib.request.urlopen(req, timeout=timeout) as resp:
        content = resp.read().decode("utf-8")

    elapsed = time.time() - start
    log.info(f"  Downloaded {len(content):,} bytes in {elapsed:.1f}s")

    # Parse the TSV
    lines = content.strip().split("\n")
    rows = [line.split("\t") for line in lines[1:] if line]  # skip header

    df = pd.DataFrame(rows, columns=["ensembl_gene_id", "hsapiens_ensembl", "pct_id"])
    df["hsapiens_ensembl"] = df["hsapiens_ensembl"].replace("", None)
    df["pct_id"] = pd.to_numeric(df["pct_id"], errors="coerce")

    # Save cache
    df.to_csv(cache_path, sep="\t", index=False)
    log.info(f"  Saved to {cache_path}")
    log.info(f"  Total: {len(df):,} mouse genes, {df['hsapiens_ensembl'].notna().sum():,} with human orthologs")

    return df


def build_ortholog_map(ortholog_df: pd.DataFrame) -> dict[str, str]:
    """
    Build ENSMUSG → ENSG mapping using highest-confidence orthologs.
    For many-to-one: keep the mouse gene with highest pct_id.
    Returns: {ENSMUSG_id: ENSG_id}
    """
    # Keep only rows with human ortholog
    df = ortholog_df.dropna(subset=["hsapiens_ensembl"]).copy()
    df = df[df["hsapiens_ensembl"].str.startswith("ENSG", na=False)]

    # For cases where multiple mouse genes map to same human gene,
    # keep the one with highest percent identity
    df = df.sort_values("pct_id", ascending=False, na_position="last")
    df = df.drop_duplicates(subset="hsapiens_ensembl", keep="first")

    # For cases where one mouse gene maps to multiple human genes,
    # keep the one with highest percent identity
    df = df.drop_duplicates(subset="ensembl_gene_id", keep="first")

    mapping = dict(zip(df["ensembl_gene_id"], df["hsapiens_ensembl"]))
    log.info(f"Ortholog map: {len(mapping):,} mouse→human gene pairs (unique 1:1)")
    return mapping


# ─── Step 2: Geneformer Vocabulary ─────────────────────────────────────────────

def load_geneformer_vocab(model_version: str = "v1",
                          mouse_gf_base: Optional[Path] = None) -> tuple[dict, dict]:
    """
    Load token dictionary and gene median expression values.
    Returns: (token_dict, gene_median_dict)
      Human Geneformer (v1/v2): {ENSG_id: token_id / median}
      Mouse-Geneformer (mouse_gf): {ENSMUSG_id: token_id / median}
    """
    config = MODEL_CONFIGS[model_version]
    log.info(f"Loading {model_version} vocabulary...")

    if model_version == "mouse_gf":
        base_dir = mouse_gf_base or MOUSE_GF_BASE
        token_file = base_dir / "MLM-re_token_dictionary_v1.pkl"
        median_file = base_dir / "mouse_gene_median_dictionary.pkl"
        missing = [str(p) for p in (token_file, median_file) if not p.exists()]
        if missing:
            raise FileNotFoundError(
                "Mouse-Geneformer dictionary files not found:\n"
                + "\n".join(f"  - {m}" for m in missing)
                + "\nSet --mouse-gf-base or MOUSE_GF_BASE to the directory containing these files."
            )
        # Load local pkl files directly (no HuggingFace download needed)
        with open(token_file, "rb") as f:
            token_dict = pickle.load(f)
        with open(median_file, "rb") as f:
            gene_median_dict = pickle.load(f)
        log.info(f"  Mouse-Geneformer base: {base_dir}")
        log.info(f"  Token vocab: {len(token_dict):,} tokens (ENSMUSG + special)")
        log.info(f"  Gene medians: {len(gene_median_dict):,} genes")
        return token_dict, gene_median_dict

    from huggingface_hub import hf_hub_download

    token_path = hf_hub_download(repo_id=HF_REPO, filename=config["token_dict_file"])
    with open(token_path, "rb") as f:
        token_dict = pickle.load(f)

    median_path = hf_hub_download(repo_id=HF_REPO, filename=config["gene_median_file"])
    with open(median_path, "rb") as f:
        gene_median_dict = pickle.load(f)

    log.info(f"  Token vocab: {len(token_dict):,} tokens (ENSG + special)")
    log.info(f"  Gene medians: {len(gene_median_dict):,} genes")

    return token_dict, gene_median_dict


# ─── Step 3: Tokenization ──────────────────────────────────────────────────────

def tokenize_sample(
    expression: dict[str, float],           # {ENSMUSG_id: expression_value}
    ortholog_map: Optional[dict[str, str]], # {ENSMUSG_id: ENSG_id}; None for mouse_gf
    token_dict: dict[str, int],             # {gene_id: token_id}
    gene_median_dict: dict[str, float],     # {gene_id: median_expression}
    max_length: int = 2048,
    special_token: bool = False,
) -> list[int]:
    """
    Convert a bulk RNA-seq expression profile to a Geneformer token sequence.

    Steps:
      1. Map genes to vocab keys (ENSMUSG→ENSG for human Geneformer; identity for mouse_gf)
      2. Keep only genes with positive expression
      3. Scale by gene median expression (makes expression comparable across genes)
      4. Rank genes by scaled expression (descending) → rank value encoding
      5. Map to token IDs
      6. Truncate to max_length

    Note: expression values should be log2(count+1) normalized (as from OSDR DESeq2)
    """
    # Step 1: Map to vocab keys
    if ortholog_map is None:
        # Mouse-Geneformer: ENSMUSG directly in token_dict (no ortholog mapping)
        ensg_expr = {g: float(v) for g, v in expression.items() if g in token_dict}
    else:
        # Human Geneformer: ENSMUSG → ENSG ortholog mapping
        ensg_expr = {}
        for ensmusg, val in expression.items():
            ensg = ortholog_map.get(ensmusg)
            if ensg and ensg in token_dict:
                ensg_expr[ensg] = float(val)

    # Step 2: Keep positive expression only
    ensg_expr = {g: v for g, v in ensg_expr.items() if v > 0}

    if not ensg_expr:
        log.warning("No expressed genes after ortholog mapping!")
        return [token_dict.get("<pad>", 0)]

    # Step 3: Scale by gene median expression
    # Geneformer uses: normalized_value = expression / median
    # This makes expression comparable across genes with different baseline levels
    genes = list(ensg_expr.keys())
    raw_expr = np.array([ensg_expr[g] for g in genes], dtype=np.float32)

    medians = np.array([
        gene_median_dict.get(g, 1.0) for g in genes
    ], dtype=np.float32)
    medians = np.where(medians == 0, 1.0, medians)  # avoid division by zero

    scaled_expr = raw_expr / medians

    # Step 4: Rank by scaled expression (descending)
    sorted_idx = np.argsort(-scaled_expr)
    ranked_genes = [genes[i] for i in sorted_idx]

    # Step 5: Map to token IDs
    token_ids = [token_dict[g] for g in ranked_genes if g in token_dict]

    # Step 6: Add special tokens (V2 only) and truncate
    if special_token:
        cls_id = token_dict.get("<cls>", None)
        eos_id = token_dict.get("<eos>", None)
        if cls_id is not None:
            token_ids = [cls_id] + token_ids
        # Reserve 1 slot for EOS
        token_ids = token_ids[: max_length - (1 if eos_id is not None else 0)]
        if eos_id is not None:
            token_ids = token_ids + [eos_id]
    else:
        token_ids = token_ids[:max_length]

    # Ensure Python int (not np.int16) for HuggingFace Dataset compatibility
    return [int(t) for t in token_ids]


def tokenize_dataframe(
    X: pd.DataFrame,  # samples × genes (ENSMUSG columns)
    y: pd.Series,     # labels (0.0/1.0)
    meta: pd.DataFrame,
    ortholog_map: dict[str, str],
    token_dict: dict[str, int],
    gene_median_dict: dict[str, float],
    model_version: str = "v1",
    label_col: str = "label",
) -> list[dict]:
    """
    Tokenize a DataFrame of bulk RNA-seq samples.
    Returns list of dicts suitable for HuggingFace Dataset.
    """
    config = MODEL_CONFIGS[model_version]
    max_length = config["max_length"]
    special_token = config["special_token"]

    gene_cols = [c for c in X.columns if str(c).startswith("ENSMUSG")]
    log.info(f"  Tokenizing {len(X)} samples with {len(gene_cols)} genes...")

    records = []
    for i, (sample_id, row) in enumerate(X.iterrows()):
        expression = {gene: row[gene] for gene in gene_cols}

        token_ids = tokenize_sample(
            expression=expression,
            ortholog_map=ortholog_map,
            token_dict=token_dict,
            gene_median_dict=gene_median_dict,
            max_length=max_length,
            special_token=special_token,
        )

        label = int(y.loc[sample_id]) if sample_id in y.index else -1

        record = {
            "input_ids": token_ids,
            "length": len(token_ids),
            "label": label,
            "sample_id": str(sample_id),
        }

        # Add metadata fields if available
        if sample_id in meta.index:
            for col in ["mission", "condition", "osd_id"]:
                if col in meta.columns:
                    record[col] = str(meta.loc[sample_id, col])

        records.append(record)

        if (i + 1) % 10 == 0:
            log.info(f"    {i + 1}/{len(X)} samples tokenized (last seq_len={len(token_ids)})")

    log.info(f"  Done. Avg sequence length: {np.mean([r['length'] for r in records]):.0f}")
    return records


# ─── Step 4: Process Folds ─────────────────────────────────────────────────────

def get_fold_dirs(task: str, fold_name: Optional[str] = None,
                  task_dir_name: Optional[str] = None) -> tuple[Path, list[Path]]:
    """Get (task_dir, fold_dirs) for a task."""
    task_dir = resolve_task_dir(task, task_dir_name=task_dir_name)
    log.info(f"Task directory: {task_dir.name}")

    if fold_name:
        if fold_name.startswith("fold_"):
            candidates = [task_dir / fold_name]
        else:
            candidates = [
                task_dir / f"fold_{fold_name}_test",
                task_dir / f"fold_{fold_name}_holdout",
                task_dir / f"fold_{fold_name}",
            ]
        fold_dirs = [d for d in candidates if d.exists() and d.is_dir()]
    else:
        fold_dirs = [d for d in task_dir.iterdir()
                     if d.is_dir() and d.name.startswith("fold_")]

    return task_dir, sorted(fold_dirs)


def tokenize_fold(
    fold_dir: Path,
    ortholog_map: dict[str, str],
    token_dict: dict[str, int],
    gene_median_dict: dict[str, float],
    model_version: str = "v1",
    overwrite: bool = False,
    dry_run: bool = False,
) -> dict:
    """
    Tokenize train + test data for a single LOMO fold.
    Saves output to fold_dir/geneformer_tokens/{model_version}/
    """
    out_dir = fold_dir / "geneformer_tokens" / model_version
    if out_dir.exists() and not overwrite and not dry_run:
        log.info(f"  Skipping {fold_dir.name} (already tokenized; use --overwrite)")
        return {"status": "skipped", "fold": fold_dir.name}

    log.info(f"\nTokenizing fold: {fold_dir.name}")

    # Load data
    train_X = pd.read_csv(fold_dir / "train_X.csv", index_col=0)
    train_y = pd.read_csv(fold_dir / "train_y.csv", index_col=0).squeeze()
    test_X = pd.read_csv(fold_dir / "test_X.csv", index_col=0)
    test_y = pd.read_csv(fold_dir / "test_y.csv", index_col=0).squeeze()

    train_meta_path = fold_dir / "train_meta.csv"
    test_meta_path = fold_dir / "test_meta.csv"
    train_meta = pd.read_csv(train_meta_path, index_col=0) if train_meta_path.exists() else pd.DataFrame()
    test_meta = pd.read_csv(test_meta_path, index_col=0) if test_meta_path.exists() else pd.DataFrame()

    log.info(f"  Train: {train_X.shape[0]} samples × {train_X.shape[1]} genes")
    log.info(f"  Test:  {test_X.shape[0]} samples × {test_X.shape[1]} genes")

    # Overlap analysis
    gene_cols = [c for c in train_X.columns if str(c).startswith("ENSMUSG")]
    if ortholog_map is None:
        # Mouse-Geneformer: ENSMUSG directly in token_dict
        n_with_ortholog = len(gene_cols)
        n_in_vocab = sum(1 for g in gene_cols if g in token_dict)
        log.info(f"  Genes in Mouse-Geneformer vocab: {n_in_vocab:,}/{len(gene_cols):,} "
                 f"({100*n_in_vocab/len(gene_cols):.1f}%)")
    else:
        n_with_ortholog = sum(1 for g in gene_cols if g in ortholog_map)
        n_in_vocab = sum(
            1 for g in gene_cols
            if g in ortholog_map and ortholog_map[g] in token_dict
        )
        log.info(f"  Genes with human ortholog: {n_with_ortholog:,}/{len(gene_cols):,}")
        log.info(f"  Genes in Geneformer vocab:  {n_in_vocab:,}/{len(gene_cols):,}")

    if dry_run:
        log.info("  [DRY RUN] Stopping here (--dry-run)")
        # Tokenize one sample to check
        sample_id = train_X.index[0]
        expression = {g: train_X.loc[sample_id, g] for g in gene_cols}
        tokens = tokenize_sample(
            expression, ortholog_map, token_dict, gene_median_dict,
            max_length=MODEL_CONFIGS[model_version]["max_length"],
            special_token=MODEL_CONFIGS[model_version]["special_token"],
        )
        log.info(f"  Sample token sequence length: {len(tokens)}")
        log.info(f"  First 10 token IDs: {tokens[:10]}")
        return {
            "status": "dry_run",
            "fold": fold_dir.name,
            "n_genes": len(gene_cols),
            "n_with_ortholog": n_with_ortholog,
            "n_in_vocab": n_in_vocab,
            "sample_seq_len": len(tokens),
        }

    # Tokenize
    train_records = tokenize_dataframe(
        train_X, train_y, train_meta,
        ortholog_map, token_dict, gene_median_dict, model_version
    )
    test_records = tokenize_dataframe(
        test_X, test_y, test_meta,
        ortholog_map, token_dict, gene_median_dict, model_version
    )

    # Save as HuggingFace Dataset (Arrow format)
    from datasets import Dataset
    out_dir.mkdir(parents=True, exist_ok=True)

    train_ds = Dataset.from_list(train_records)
    test_ds = Dataset.from_list(test_records)

    train_ds.save_to_disk(str(out_dir / "train_dataset"))
    test_ds.save_to_disk(str(out_dir / "test_dataset"))

    # Save tokenization metadata
    info = {
        "fold": fold_dir.name,
        "model_version": model_version,
        "max_length": MODEL_CONFIGS[model_version]["max_length"],
        "n_train": len(train_records),
        "n_test": len(test_records),
        "avg_train_seq_len": float(np.mean([r["length"] for r in train_records])),
        "avg_test_seq_len": float(np.mean([r["length"] for r in test_records])),
        "n_genes_input": len(gene_cols),
        "n_with_ortholog": n_with_ortholog,
        "n_in_vocab": n_in_vocab,
        "ortholog_coverage_pct": round(100 * n_in_vocab / len(gene_cols), 1) if gene_cols else 0.0,
    }

    with open(out_dir / "tokenize_info.json", "w") as f:
        json.dump(info, f, indent=2)

    log.info(f"  Saved to {out_dir}")
    log.info(f"  Train avg seq len: {info['avg_train_seq_len']:.0f}")
    log.info(f"  Test avg seq len:  {info['avg_test_seq_len']:.0f}")

    return {"status": "ok", **info}


# ─── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--task", default="A4", help="Task to tokenize (e.g. A1, A4, A6)")
    parser.add_argument("--fold", default=None, help="Specific fold name (e.g. RR-9); default=all folds")
    parser.add_argument("--model-version", default="mouse_gf", choices=["v1", "v2", "mouse_gf"],
                        help="Model version: mouse_gf (default), v1, v2")
    parser.add_argument("--overwrite", action="store_true", help="Re-tokenize even if output exists")
    parser.add_argument("--dry-run", action="store_true", help="Run overlap analysis + 1 sample, no output saved")
    parser.add_argument("--ortholog-cache", default=str(ORTHOLOG_CACHE), help="Path to cached ortholog TSV")
    parser.add_argument("--task-dir", default=None,
                        help="Explicit task directory name under tasks/ (recommended for A1 variants)")
    parser.add_argument("--mouse-gf-base", default=str(MOUSE_GF_BASE),
                        help="Base directory containing Mouse-Geneformer dictionaries")
    args = parser.parse_args()

    # Step 1: Ortholog mapping (skipped for mouse_gf)
    if args.model_version == "mouse_gf":
        ortholog_map = None
        log.info("Mouse-Geneformer: skipping ortholog mapping (ENSMUSG used directly)")
    else:
        cache_path = Path(args.ortholog_cache)
        ortholog_df = download_ortholog_table(cache_path)
        ortholog_map = build_ortholog_map(ortholog_df)

    # Step 2: Vocabulary
    token_dict, gene_median_dict = load_geneformer_vocab(
        args.model_version,
        mouse_gf_base=Path(args.mouse_gf_base),
    )

    # Step 3: Find folds
    task_dir, fold_dirs = get_fold_dirs(args.task, args.fold, task_dir_name=args.task_dir)
    if not fold_dirs:
        available = sorted(
            d.name for d in task_dir.iterdir()
            if d.is_dir() and d.name.startswith("fold_")
        )
        raise FileNotFoundError(
            f"No fold found matching '{args.fold}' in {task_dir}. "
            f"Available folds: {available}"
        )
    log.info(f"\nFounds {len(fold_dirs)} fold(s) for task {args.task}:")
    for d in fold_dirs:
        log.info(f"  {d.name}")

    # Step 4: Tokenize each fold
    results = []
    for fold_dir in fold_dirs:
        result = tokenize_fold(
            fold_dir=fold_dir,
            ortholog_map=ortholog_map,
            token_dict=token_dict,
            gene_median_dict=gene_median_dict,
            model_version=args.model_version,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
        )
        results.append(result)

    # Summary
    print("\n=== Tokenization Summary ===")
    for r in results:
        status = r.get("status", "?")
        fold = r.get("fold", "?")
        if status == "ok":
            print(f"  {fold}: OK | train_seq={r.get('avg_train_seq_len', 0):.0f} | "
                  f"test_seq={r.get('avg_test_seq_len', 0):.0f} | "
                  f"vocab_coverage={r.get('ortholog_coverage_pct', 0):.1f}%")
        elif status == "dry_run":
            print(f"  {fold}: DRY_RUN | n_genes={r.get('n_genes', 0):,} | "
                  f"in_vocab={r.get('n_in_vocab', 0):,} | "
                  f"sample_seq_len={r.get('sample_seq_len', 0)}")
        else:
            print(f"  {fold}: {status}")

    # Save overall summary
    summary_path = task_dir / f"geneformer_{args.model_version}_tokenize_summary.json"
    if not args.dry_run:
        with open(summary_path, "w") as f:
            json.dump(results, f, indent=2)
        log.info(f"\nSummary saved to {summary_path}")


if __name__ == "__main__":
    main()

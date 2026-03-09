"""
scGPT Tokenization Pipeline for Mouse Bulk RNA-seq (NASA GeneLab)

Converts mouse bulk RNA-seq expression profiles to scGPT token format using:
  1. ENSMUSG → human gene symbol via extended BioMart ortholog table
  2. scGPT GeneVocab for symbol → token ID mapping
  3. Expression binning (51 bins, matching whole_human pretrain config)
  4. Saves as PyTorch .pt files for downstream fine-tuning

IMPORTANT DIFFERENCES vs Geneformer tokenization:
  - Geneformer: rank-value encoding (gene rank by expression)
  - scGPT: gene ID tokens + binned expression values (continuous value encoding)
  - scGPT vocab: human gene symbols (not ENSMUSG, not ENSG)
  - Coverage: ~50-67% of mouse genes have human orthologs in scGPT vocab

Usage:
  python scripts/scgpt_tokenize.py --task A4 --fold RR-23 --model-version whole_human
  python scripts/scgpt_tokenize.py --task A4 --fold lomo   # all LOMO folds
  python scripts/scgpt_tokenize.py --all                    # all 6 tissues
  python scripts/scgpt_tokenize.py --task A4 --dry-run      # shape check only

Output:
  tasks/{task_dir}/fold_{fold}/scgpt_tokens/
    train_data.pt    # {'gene_ids': Tensor, 'values': Tensor, 'labels': Tensor}
    test_data.pt
    tokenize_info.json

Architecture:
  scGPT whole_human: 12L, 512d, 8h, n_bins=51, max_seq_len=1200
  Pretrained on 33M human cells (CellXGene census, June 2023)
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import torch

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ─── Paths ─────────────────────────────────────────────────────────────────────

ROOT = Path(__file__).parent.parent
TASKS_DIR = ROOT / "tasks"

# HPC paths (overridable via env vars)
MODEL_BASE = Path(
    os.environ.get(
        "SCGPT_MODEL_BASE",
        "/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/models",
    )
)
DATA_BASE = Path(
    os.environ.get(
        "SCGPT_DATA_BASE",
        "/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/data",
    )
)
HPC_TASKS_BASE = Path(
    os.environ.get(
        "SCGPT_TASKS_BASE",
        "/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/tasks",
    )
)

# Model configs
MODEL_CONFIGS = {
    "whole_human": {
        "model_dir": str(MODEL_BASE / "scgpt_whole_human"),
        "vocab_file": "vocab.json",
        "model_file": "best_model.pt",
        "args_file": "args.json",
        "n_bins": 51,
        "max_seq_len": 1200,
        "pad_token": "<pad>",
        "pad_value": -2,
        "cls_token": "<cls>",
    }
}

# Tissue → task dir name → LOMO folds (test missions)
TISSUE_TASKS = {
    "A1": ("A1_liver_lomo", ["RR-1", "RR-3", "RR-6", "RR-8", "RR-9", "MHU-2"]),
    "A2": ("A2_gastrocnemius_lomo", ["RR-1", "RR-3", "RR-5"]),
    "A3": ("A3_kidney_lomo", ["RR-1", "RR-3", "RR-5"]),
    "A4": ("A4_thymus_lomo", ["MHU-1", "MHU-2", "RR-6", "RR-9"]),
    "A5": ("A5_skin_lomo", ["RR-6", "RR-7", "MHU-2"]),
    "A6": ("A6_eye_lomo", ["RR-1", "RR-3", "RR-5"]),
}


# ─── Ortholog Mapping ──────────────────────────────────────────────────────────

def load_ortholog_map(data_base: Path) -> dict[str, str]:
    """
    Load ENSMUSG → human gene symbol mapping from extended BioMart ortholog table.
    File: data/mouse/ensembl_mouse_human_orthologs_with_symbols.tsv
    Columns: Gene stable ID | Human gene stable ID | Human gene name | %id

    Returns: {ENSMUSG: gene_symbol}
    """
    ortho_file = data_base / "mouse" / "ensembl_mouse_human_orthologs_with_symbols.tsv"
    if not ortho_file.exists():
        raise FileNotFoundError(
            f"Ortholog file not found: {ortho_file}\n"
            "Expected: data/mouse/ensembl_mouse_human_orthologs_with_symbols.tsv\n"
            "Columns: ensembl_gene_id, hsapiens_ensembl, hsapiens_gene_name, pct_id"
        )

    mapping = {}
    with open(ortho_file) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[0] and parts[2]:
                ensmusg = parts[0]
                symbol = parts[2]
                if ensmusg not in mapping:
                    mapping[ensmusg] = symbol

    log.info(f"Ortholog map: {len(mapping):,} ENSMUSG → gene_symbol pairs")
    return mapping


# ─── scGPT Vocabulary ──────────────────────────────────────────────────────────

def load_vocab(model_dir: Path, vocab_file: str = "vocab.json") -> dict:
    """Load scGPT vocab: {gene_symbol: token_id}."""
    vocab_path = model_dir / vocab_file
    if not vocab_path.exists():
        raise FileNotFoundError(f"Vocab not found: {vocab_path}")
    with open(vocab_path) as f:
        vocab = json.load(f)
    n_genes = sum(1 for k in vocab if not k.startswith("<"))
    log.info(f"scGPT vocab: {len(vocab):,} total ({n_genes:,} genes + special tokens)")
    return vocab


# ─── Tokenization ──────────────────────────────────────────────────────────────

def bin_expression(values: np.ndarray, n_bins: int = 51) -> np.ndarray:
    """
    Bin expression values into n_bins discrete bins (matching scGPT pretraining).
    scGPT uses: normalize_total(1e4) → log1p → binning
    Input: 1D array of log2(count+1) values (already log-normalized from OSDR)
    Output: integer bin indices [0, n_bins-1]
    """
    if len(values) == 0:
        return np.array([], dtype=np.int64)

    # Remove zeros (scGPT excludes zero-expression genes)
    nonzero_mask = values > 0
    if nonzero_mask.sum() == 0:
        return np.zeros(len(values), dtype=np.int64)

    # Compute bin edges from non-zero values
    nonzero_vals = values[nonzero_mask]
    vmin, vmax = nonzero_vals.min(), nonzero_vals.max()

    if vmax == vmin:
        binned = np.zeros(len(values), dtype=np.int64)
        binned[nonzero_mask] = n_bins // 2
        return binned

    # Linear binning matching scGPT's approach
    bin_edges = np.linspace(vmin, vmax, n_bins - 1)
    binned = np.zeros(len(values), dtype=np.int64)
    binned[nonzero_mask] = np.digitize(nonzero_vals, bin_edges)
    return binned


def tokenize_sample(
    expression: dict[str, float],   # {ENSMUSG: log2_expr}
    ortholog_map: dict[str, str],   # {ENSMUSG: gene_symbol}
    vocab: dict[str, int],          # {gene_symbol: token_id}
    n_bins: int = 51,
    max_seq_len: int = 1200,
    include_cls: bool = True,
    cls_token: str = "<cls>",
    pad_token: str = "<pad>",
    pad_value: int = -2,
) -> tuple[list[int], list[int]]:
    """
    Convert bulk RNA-seq sample to scGPT token + value sequences.

    Steps:
      1. ENSMUSG → gene symbol via ortholog map
      2. Filter: keep only genes in scGPT vocab with nonzero expression
      3. Sort by expression (descending) — most expressed first
      4. Bin expression values (51 bins)
      5. Truncate to max_seq_len - 1 (reserve 1 for CLS)
      6. Prepend CLS token if include_cls=True

    Returns:
      gene_ids: list of token IDs
      binned_values: list of binned expression values
    """
    # Step 1-2: Map and filter
    symbol_expr = {}
    for ensmusg, expr in expression.items():
        if ensmusg in ortholog_map and expr > 0:
            symbol = ortholog_map[ensmusg]
            if symbol in vocab:
                # Keep highest expression if multiple ENSMUSG → same symbol
                if symbol not in symbol_expr or expr > symbol_expr[symbol]:
                    symbol_expr[symbol] = expr

    if not symbol_expr:
        return [], []

    # Step 3: Sort by expression descending
    sorted_items = sorted(symbol_expr.items(), key=lambda x: x[1], reverse=True)

    # Truncate to fit max_seq_len (reserve 1 slot for CLS)
    max_genes = max_seq_len - (1 if include_cls else 0)
    sorted_items = sorted_items[:max_genes]

    symbols = [s for s, _ in sorted_items]
    expr_vals = np.array([e for _, e in sorted_items])

    # Step 4: Bin expression
    binned = bin_expression(expr_vals, n_bins=n_bins)

    # Step 5: Build token IDs
    gene_ids = [vocab[s] for s in symbols]
    values = binned.tolist()

    # Step 6: Prepend CLS
    if include_cls and cls_token in vocab:
        gene_ids = [vocab[cls_token]] + gene_ids
        values = [n_bins - 1] + values  # CLS gets max bin value

    return gene_ids, values


# ─── Fold Loading & Tokenization ───────────────────────────────────────────────

def tokenize_fold(
    fold_dir: Path,
    ortholog_map: dict[str, str],
    vocab: dict[str, int],
    config: dict,
    dry_run: bool = False,
) -> dict:
    """
    Tokenize a single LOMO fold (train + test).

    Returns summary stats dict.
    """
    out_dir = fold_dir / "scgpt_tokens"
    out_dir.mkdir(parents=True, exist_ok=True)

    info_file = out_dir / "tokenize_info.json"
    # Only skip if actual .pt data files exist (not just dry-run info)
    train_pt = out_dir / "train_data.pt"
    test_pt = out_dir / "test_data.pt"
    if not dry_run and train_pt.exists() and test_pt.exists() and info_file.exists():
        log.info(f"  Already tokenized: {out_dir} — skipping")
        with open(info_file) as f:
            return json.load(f)

    n_bins = config["n_bins"]
    max_seq_len = config["max_seq_len"]
    pad_token = config["pad_token"]
    pad_value = config["pad_value"]

    results = {}
    for split in ["train", "test"]:
        x_file = fold_dir / f"{split}_X.csv"
        y_file = fold_dir / f"{split}_y.csv"

        if not x_file.exists():
            log.warning(f"  Missing {x_file} — skipping {split}")
            continue

        log.info(f"  Loading {split}_X.csv...")
        df = pd.read_csv(x_file, index_col=0)
        y = pd.read_csv(y_file, index_col=0).iloc[:, 0].values
        n_samples, n_genes = df.shape
        log.info(f"    {n_samples} samples × {n_genes} genes")

        if dry_run:
            # Just check coverage without full tokenization
            gene_ids_example = df.columns.tolist()
            n_mapped = sum(
                1 for g in gene_ids_example
                if g in ortholog_map and ortholog_map[g] in vocab
            )
            log.info(f"    Dry run: {n_mapped}/{n_genes} genes would map to scGPT vocab ({100*n_mapped/n_genes:.1f}%)")
            results[split] = {"n_samples": n_samples, "n_genes": n_genes, "n_mapped": n_mapped}
            continue

        all_gene_ids, all_values, all_labels = [], [], []
        seq_lengths = []

        for i, (sample_id, row) in enumerate(df.iterrows()):
            expression = {gene: float(val) for gene, val in row.items()}
            gene_ids, values = tokenize_sample(
                expression, ortholog_map, vocab,
                n_bins=n_bins, max_seq_len=max_seq_len,
            )
            all_gene_ids.append(gene_ids)
            all_values.append(values)
            all_labels.append(int(y[i]))
            seq_lengths.append(len(gene_ids))

            if (i + 1) % 20 == 0:
                log.info(f"    Tokenized {i+1}/{n_samples} samples...")

        # Pad to max seq len in this split
        max_len = max(len(ids) for ids in all_gene_ids) if all_gene_ids else 1
        pad_id = vocab[pad_token]

        gene_ids_padded = torch.zeros(n_samples, max_len, dtype=torch.long)
        values_padded = torch.full((n_samples, max_len), pad_value, dtype=torch.long)
        labels_tensor = torch.tensor(all_labels, dtype=torch.long)

        for i, (gids, vals) in enumerate(zip(all_gene_ids, all_values)):
            seq_len = len(gids)
            if seq_len > 0:
                gene_ids_padded[i, :seq_len] = torch.tensor(gids, dtype=torch.long)
                values_padded[i, :seq_len] = torch.tensor(vals, dtype=torch.long)

        # Save
        out_file = out_dir / f"{split}_data.pt"
        torch.save({
            "gene_ids": gene_ids_padded,
            "values": values_padded,
            "labels": labels_tensor,
        }, out_file)

        log.info(f"    Saved: {out_file} ({gene_ids_padded.shape})")

        results[split] = {
            "n_samples": n_samples,
            "n_genes": n_genes,
            "max_seq_len_actual": max_len,
            "mean_seq_len": float(np.mean(seq_lengths)),
            "median_seq_len": float(np.median(seq_lengths)),
            "gene_ids_shape": list(gene_ids_padded.shape),
        }

    # Save info
    info = {
        "model_version": "whole_human",
        "n_bins": n_bins,
        "max_seq_len": max_seq_len,
        "fold_dir": str(fold_dir),
        "splits": results,
    }
    with open(info_file, "w") as f:
        json.dump(info, f, indent=2)

    return info


# ─── Main ───────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="scGPT tokenization for mouse bulk RNA-seq")
    p.add_argument("--task", choices=list(TISSUE_TASKS.keys()),
                   help="Task ID (A1-A6)")
    p.add_argument("--task-dir", help="Override task directory name")
    p.add_argument("--fold", default="lomo",
                   help="Fold (test mission name) or 'lomo' for all folds")
    p.add_argument("--model-version", default="whole_human",
                   choices=list(MODEL_CONFIGS.keys()))
    p.add_argument("--all", action="store_true",
                   help="Tokenize all 6 tissues (all LOMO folds)")
    p.add_argument("--dry-run", action="store_true",
                   help="Check coverage without tokenizing")
    # Path overrides (for local testing)
    p.add_argument("--model-base", default=None, help="Override model base path")
    p.add_argument("--data-base", default=None, help="Override data base path")
    p.add_argument("--tasks-base", default=None, help="Override tasks base path")
    return p.parse_args()


def main():
    args = parse_args()

    # Path overrides
    model_base = Path(args.model_base) if args.model_base else MODEL_BASE
    data_base = Path(args.data_base) if args.data_base else DATA_BASE
    tasks_base = Path(args.tasks_base) if args.tasks_base else HPC_TASKS_BASE

    config = MODEL_CONFIGS[args.model_version]
    model_dir = Path(config["model_dir"].replace(str(MODEL_BASE), str(model_base)))

    log.info(f"Model version: {args.model_version}")
    log.info(f"Model dir: {model_dir}")
    log.info(f"Data base: {data_base}")

    # Load ortholog map and vocab
    ortholog_map = load_ortholog_map(data_base)
    vocab = load_vocab(model_dir, config["vocab_file"])

    # Determine which tasks/folds to tokenize
    if args.all:
        task_list = list(TISSUE_TASKS.keys())
    elif args.task:
        task_list = [args.task]
    else:
        raise ValueError("Specify --task <A1-A6> or --all")

    total_folds = 0
    for task_id in task_list:
        task_dir_name, all_folds = TISSUE_TASKS[task_id]
        task_dir = tasks_base / task_dir_name

        if not task_dir.exists():
            log.warning(f"Task dir not found: {task_dir} — skipping")
            continue

        # Which folds to tokenize
        if args.fold == "lomo" or args.all:
            folds_to_run = all_folds
        else:
            folds_to_run = [args.fold]

        log.info(f"\n{'='*60}")
        log.info(f"Task: {task_id} ({task_dir_name})")
        log.info(f"Folds: {folds_to_run}")
        log.info(f"{'='*60}")

        for fold_mission in folds_to_run:
            # Try _test and _holdout suffixes
            fold_dir = None
            for suffix in ["_test", "_holdout"]:
                candidate = task_dir / f"fold_{fold_mission}{suffix}"
                if candidate.exists():
                    fold_dir = candidate
                    break

            if fold_dir is None:
                log.warning(f"  Fold dir not found for {fold_mission} — skipping")
                continue

            log.info(f"\n  Fold: {fold_mission} → {fold_dir.name}")
            t0 = time.time()
            info = tokenize_fold(
                fold_dir, ortholog_map, vocab, config,
                dry_run=args.dry_run,
            )
            elapsed = time.time() - t0
            log.info(f"  Done in {elapsed:.1f}s")
            total_folds += 1

    log.info(f"\nTotal folds processed: {total_folds}")


if __name__ == "__main__":
    main()

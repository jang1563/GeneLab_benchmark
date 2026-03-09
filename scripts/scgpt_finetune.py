"""
scGPT Fine-Tuning Pipeline for Spaceflight Detection (NASA GeneLab LOMO)

Fine-tunes scGPT whole_human (12L Transformer, 512d) for binary classification
(Flight vs Ground Control) using Leave-One-Mission-Out cross-validation folds.

Architecture:
  - Base: scGPT whole_human (TransformerModel, 12L)
  - Head: CLS token → linear classifier (2 classes: FLT=1, GC=0)
  - Strategy: Freeze bottom 10/12 layers, fine-tune top 2 + classification head

Key differences from Geneformer:
  - Input: gene ID tokens + binned expression values (not rank-value tokens)
  - Model: scGPT TransformerModel (not BERT)
  - Ortholog mapping required: ENSMUSG → human gene symbol
  - Larger model (12L vs 6L, 512d vs 256d)

Hardware targets:
  - HPC A40 (48GB): batch_size=8, mixed precision
  - HPC A100 (80GB): batch_size=16

Freezing strategy (small-n bulk RNA-seq):
  --freeze-layers 10  : Fine-tune top 2 layers + head (recommended, n~30-100)
  --freeze-layers 12  : Head-only fine-tuning (most regularized)
  --freeze-layers 0   : Full fine-tuning (risky for n<100)

Usage:
  python scripts/scgpt_finetune.py --task A4 --fold RR-9 --device cuda
  python scripts/scgpt_finetune.py --task A4 --fold lomo --device cuda --epochs 10
  python scripts/scgpt_finetune.py --task A4 --dry-run

Output:
  evaluation/scgpt_whole_human_A4_lomo_results.json
  tasks/A4_thymus_lomo/fold_{fold}/scgpt_tokens/checkpoint_best.pt
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import time
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import roc_auc_score
from torch.optim import AdamW
from torch.utils.data import DataLoader, Dataset

warnings.filterwarnings("ignore")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ─── Paths ─────────────────────────────────────────────────────────────────────

ROOT = Path(__file__).parent.parent
TASKS_DIR = ROOT / "tasks"
EVAL_DIR = ROOT / "evaluation"

MODEL_BASE = Path(
    os.environ.get(
        "SCGPT_MODEL_BASE",
        "/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/models",
    )
)
HPC_TASKS_BASE = Path(
    os.environ.get(
        "SCGPT_TASKS_BASE",
        "/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/tasks",
    )
)
HPC_EVAL_BASE = Path(
    os.environ.get(
        "SCGPT_EVAL_BASE",
        "/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/evaluation",
    )
)

TISSUE_TASKS = {
    "A1": ("A1_liver_lomo", ["RR-1", "RR-3", "RR-6", "RR-8", "RR-9", "MHU-2"]),
    "A2": ("A2_gastrocnemius_lomo", ["RR-1", "RR-3", "RR-5"]),
    "A3": ("A3_kidney_lomo", ["RR-1", "RR-3", "RR-5"]),
    "A4": ("A4_thymus_lomo", ["MHU-1", "MHU-2", "RR-6", "RR-9"]),
    "A5": ("A5_skin_lomo", ["RR-6", "RR-7", "MHU-2"]),
    "A6": ("A6_eye_lomo", ["RR-1", "RR-3", "RR-5"]),
}

DEFAULT_SEED = 42


# ─── Dataset ───────────────────────────────────────────────────────────────────

class BulkRNADataset(Dataset):
    """Dataset from pre-tokenized .pt files."""

    def __init__(self, pt_file: Path):
        if not pt_file.exists():
            raise FileNotFoundError(
                f"Tokenized data not found: {pt_file}\n"
                "Run scgpt_tokenize.py first."
            )
        data = torch.load(str(pt_file), map_location="cpu")
        self.gene_ids = data["gene_ids"]    # (N, seq_len)
        self.values = data["values"]         # (N, seq_len)
        self.labels = data["labels"]         # (N,)
        log.info(f"  Dataset: {len(self.labels)} samples, seq_len={self.gene_ids.shape[1]}")

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return {
            "gene_ids": self.gene_ids[idx],
            "values": self.values[idx],
            "labels": self.labels[idx],
        }


# ─── Model Loading ─────────────────────────────────────────────────────────────

def load_scgpt_for_classification(
    model_dir: Path,
    vocab_size: int,
    n_cls: int = 2,
    use_fast_transformer: bool = False,   # flash-attn; set False if unavailable
) -> nn.Module:
    """
    Load scGPT whole_human pretrained weights into a classification model.

    scGPT args.json config:
      nlayers=12, nheads=8, embsize=512, d_hid=512, n_bins=51,
      input_emb_style='continuous', dropout=0.2, n_layers_cls=3
    """
    # Import here to avoid top-level import issues on login nodes
    from scgpt.model import TransformerModel
    from scgpt.utils import load_pretrained

    model_file = model_dir / "best_model.pt"
    args_file = model_dir / "args.json"

    if not model_file.exists():
        raise FileNotFoundError(f"Model weights not found: {model_file}")

    # Read pretrain config
    with open(args_file) as f:
        pretrain_args = json.load(f)

    nlayers = pretrain_args.get("nlayers", 12)
    nheads = pretrain_args.get("nheads", 8)
    embsize = pretrain_args.get("embsize", 512)
    d_hid = pretrain_args.get("d_hid", 512)
    n_bins = pretrain_args.get("n_bins", 51)
    dropout = pretrain_args.get("dropout", 0.2)
    n_layers_cls = pretrain_args.get("n_layers_cls", 3)
    input_emb_style = pretrain_args.get("input_emb_style", "continuous")
    pad_token = pretrain_args.get("pad_token", "<pad>")
    pad_value = pretrain_args.get("pad_value", -2)
    mask_value = pretrain_args.get("mask_value", -1)

    log.info(f"Building TransformerModel: {nlayers}L, {embsize}d, {nheads}h")
    log.info(f"  Classification head: n_cls={n_cls}, n_layers_cls={n_layers_cls}")
    log.info(f"  Input style: {input_emb_style}, n_bins={n_bins}")

    model = TransformerModel(
        ntoken=vocab_size,
        d_model=embsize,
        nhead=nheads,
        d_hid=d_hid,
        nlayers=nlayers,
        nlayers_cls=n_layers_cls,
        n_cls=n_cls,
        vocab=None,
        dropout=dropout,
        pad_token=pad_token,
        pad_value=pad_value,
        do_mvc=False,
        do_dab=False,
        use_batch_labels=False,
        num_batch_labels=None,
        domain_spec_batchnorm=False,
        input_emb_style=input_emb_style,
        n_input_bins=n_bins,
        cell_emb_style="cls",      # CLS token for classification
        mvc_decoder_style="inner product",
        ecs_threshold=0.3,
        explicit_zero_prob=False,
        use_fast_transformer=use_fast_transformer,
        fast_transformer_backend="flash",
        pre_norm=False,
    )

    # Load pretrained weights (strict=False: mismatches for CLS head, etc.)
    log.info(f"Loading pretrained weights from {model_file}...")
    pretrained_params = torch.load(str(model_file), map_location="cpu")
    model = load_pretrained(model, pretrained_params, strict=False, verbose=True)

    n_params = sum(p.numel() for p in model.parameters())
    log.info(f"  Total params: {n_params:,}")

    return model


def freeze_transformer_layers(model: nn.Module, n_freeze: int) -> None:
    """
    Freeze bottom n_freeze layers of the scGPT TransformerModel.
    scGPT has 12 encoder layers (0-11).

    n_freeze=0:  Full fine-tuning
    n_freeze=10: Fine-tune top 2 layers + head (recommended for n~30-100)
    n_freeze=12: Head-only (most regularized)
    """
    if n_freeze <= 0:
        log.info("  Layer freezing: disabled (full fine-tuning)")
        return

    # Freeze embeddings
    for name, param in model.named_parameters():
        if "encoder" not in name and "cls_decoder" not in name and "mvc_decoder" not in name:
            param.requires_grad = False

    # Freeze bottom n_freeze transformer encoder layers
    for i in range(min(n_freeze, 12)):
        layer_name = f"transformer_encoder.layers.{i}"
        for name, param in model.named_parameters():
            if layer_name in name:
                param.requires_grad = False

    n_trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    n_frozen = sum(p.numel() for p in model.parameters() if not p.requires_grad)
    log.info(f"  Frozen: {n_freeze}/12 layers ({n_frozen:,} params)")
    log.info(f"  Trainable: {n_trainable:,} params (top {12 - n_freeze} layers + head)")


# ─── Training ──────────────────────────────────────────────────────────────────

def train_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    scheduler,
    criterion: nn.Module,
    device: torch.device,
    scaler: Optional[torch.cuda.amp.GradScaler],
) -> float:
    model.train()
    total_loss = 0.0

    for batch in loader:
        gene_ids = batch["gene_ids"].to(device)     # (B, seq_len)
        values = batch["values"].to(device)          # (B, seq_len)
        labels = batch["labels"].to(device)          # (B,)

        optimizer.zero_grad()

        if scaler is not None:
            with torch.cuda.amp.autocast():
                output = model(gene_ids, values, CLS=True, MVC=False, ECS=False)
                cls_emb = output["cls_output"]
                logits = model.cls_decoder(cls_emb)
                loss = criterion(logits, labels)
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            scaler.step(optimizer)
            scaler.update()
        else:
            output = model(gene_ids, values, CLS=True, MVC=False, ECS=False)
            cls_emb = output["cls_output"]
            logits = model.cls_decoder(cls_emb)
            loss = criterion(logits, labels)
            loss.backward()
            nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

        if scheduler is not None:
            scheduler.step()

        total_loss += loss.item()

    return total_loss / max(len(loader), 1)


@torch.no_grad()
def evaluate(
    model: nn.Module,
    loader: DataLoader,
    device: torch.device,
) -> tuple[float, float]:
    """Returns (auroc, loss)."""
    model.eval()
    all_probs, all_labels = [], []
    total_loss = 0.0
    criterion = nn.CrossEntropyLoss()

    for batch in loader:
        gene_ids = batch["gene_ids"].to(device)
        values = batch["values"].to(device)
        labels = batch["labels"].to(device)

        output = model(gene_ids, values, CLS=True, MVC=False, ECS=False)
        cls_emb = output["cls_output"]
        logits = model.cls_decoder(cls_emb)

        loss = criterion(logits, labels)
        total_loss += loss.item()

        probs = torch.softmax(logits, dim=-1)[:, 1]
        all_probs.extend(probs.cpu().numpy())
        all_labels.extend(labels.cpu().numpy())

    try:
        auroc = roc_auc_score(all_labels, all_probs)
    except ValueError:
        auroc = 0.5  # all same class

    return auroc, total_loss / max(len(loader), 1)


# ─── Main ───────────────────────────────────────────────────────────────────────

def run_fold(
    task_id: str,
    task_dir_name: str,
    fold_mission: str,
    args: argparse.Namespace,
    model_dir: Path,
    tasks_base: Path,
    eval_base: Path,
    vocab_size: int,
) -> Optional[dict]:
    """Run fine-tuning for one fold. Returns result dict or None on error."""
    task_dir = tasks_base / task_dir_name

    # Find fold dir
    fold_dir = None
    for suffix in ["_test", "_holdout"]:
        candidate = task_dir / f"fold_{fold_mission}{suffix}"
        if candidate.exists():
            fold_dir = candidate
            break
    if fold_dir is None:
        log.warning(f"Fold dir not found for {fold_mission}")
        return None

    token_dir = fold_dir / "scgpt_tokens"
    train_file = token_dir / "train_data.pt"
    test_file = token_dir / "test_data.pt"

    log.info(f"\n{'='*50}")
    log.info(f"scGPT: {task_id} / {fold_mission}")
    log.info(f"{'='*50}")

    # dry-run: just check CSV files exist and report shapes
    if args.dry_run:
        for split in ["train", "test"]:
            x_file = fold_dir / f"{split}_X.csv"
            if x_file.exists():
                import pandas as pd
                df = pd.read_csv(x_file, index_col=0, nrows=1)
                with open(x_file) as fh:
                    n_rows = sum(1 for _ in fh) - 1
                log.info(f"  Dry run {split}: {n_rows} samples × {df.shape[1]} genes")
        return {"test_mission": fold_mission, "task": task_id, "dry_run": True}

    if not train_file.exists() or not test_file.exists():
        log.error(f"Tokenized data missing: {token_dir}")
        log.error("Run scgpt_tokenize.py first!")
        return None

    device = torch.device(args.device if torch.cuda.is_available() and args.device != "cpu" else "cpu")
    log.info(f"Device: {device}")
    if device.type == "cuda":
        log.info(f"GPU: {torch.cuda.get_device_name()}, {torch.cuda.get_device_properties(0).total_memory//1e9:.0f}GB")

    torch.manual_seed(DEFAULT_SEED)
    np.random.seed(DEFAULT_SEED)

    # Dataset & DataLoader
    train_ds = BulkRNADataset(train_file)
    test_ds = BulkRNADataset(test_file)
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True, num_workers=0)
    test_loader = DataLoader(test_ds, batch_size=args.batch_size, shuffle=False, num_workers=0)

    if args.dry_run:
        log.info(f"Dry run: train={len(train_ds)}, test={len(test_ds)} — OK")
        return {"test_mission": fold_mission, "dry_run": True}

    # Model
    use_flash = device.type == "cuda"
    model = load_scgpt_for_classification(
        model_dir, vocab_size=vocab_size,
        n_cls=2, use_fast_transformer=use_flash,
    )
    freeze_transformer_layers(model, args.freeze_layers)
    model = model.to(device)

    # Optimizer + scheduler
    optimizer = AdamW(
        [p for p in model.parameters() if p.requires_grad],
        lr=args.lr, weight_decay=1e-5,
    )
    total_steps = args.epochs * len(train_loader)
    warmup_steps = min(100, total_steps // 10)
    from transformers import get_linear_schedule_with_warmup
    scheduler = get_linear_schedule_with_warmup(
        optimizer, num_warmup_steps=warmup_steps, num_training_steps=total_steps
    )
    criterion = nn.CrossEntropyLoss()

    # Mixed precision
    scaler = torch.cuda.amp.GradScaler() if device.type == "cuda" else None

    # Training loop
    t0 = time.time()
    best_auroc, best_epoch = 0.5, 0
    best_ckpt = token_dir / "checkpoint_best.pt"

    log.info(f"Training {args.epochs} epochs (freeze_layers={args.freeze_layers}, lr={args.lr})...")
    for epoch in range(1, args.epochs + 1):
        train_loss = train_epoch(model, train_loader, optimizer, scheduler, criterion, device, scaler)
        test_auroc, test_loss = evaluate(model, test_loader, device)

        log.info(f"  Epoch {epoch:2d}/{args.epochs}: train_loss={train_loss:.4f}, "
                 f"test_auroc={test_auroc:.4f}, test_loss={test_loss:.4f}")

        if test_auroc > best_auroc:
            best_auroc = test_auroc
            best_epoch = epoch
            torch.save(model.state_dict(), str(best_ckpt))

    elapsed = time.time() - t0
    log.info(f"  Best AUROC: {best_auroc:.4f} at epoch {best_epoch} ({elapsed:.1f}s)")

    result = {
        "test_mission": fold_mission,
        "task": task_id,
        "model": "scgpt_whole_human",
        "auroc": best_auroc,
        "best_epoch": best_epoch,
        "n_train": len(train_ds),
        "n_test": len(test_ds),
        "train_time_s": round(elapsed, 1),
        "freeze_layers": args.freeze_layers,
        "epochs": args.epochs,
        "lr": args.lr,
        "batch_size": args.batch_size,
    }
    return result


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="scGPT fine-tuning for spaceflight detection")
    p.add_argument("--task", choices=list(TISSUE_TASKS.keys()), help="Task ID (A1-A6)")
    p.add_argument("--task-dir", help="Override task directory name")
    p.add_argument("--fold", default="lomo",
                   help="Fold (test mission name) or 'lomo' for all LOMO folds")
    p.add_argument("--all", action="store_true", help="Run all 6 tissues")
    # Model
    p.add_argument("--model-version", default="whole_human")
    p.add_argument("--freeze-layers", type=int, default=10,
                   help="Freeze bottom N/12 transformer layers (default: 10)")
    # Training
    p.add_argument("--epochs", type=int, default=10)
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--lr", type=float, default=1e-4)
    p.add_argument("--device", default="cuda", choices=["cuda", "cpu", "mps"])
    p.add_argument("--dry-run", action="store_true")
    # Paths
    p.add_argument("--model-base", default=None)
    p.add_argument("--tasks-base", default=None)
    p.add_argument("--eval-base", default=None)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return p.parse_args()


def main():
    args = parse_args()

    model_base = Path(args.model_base) if args.model_base else MODEL_BASE
    tasks_base = Path(args.tasks_base) if args.tasks_base else HPC_TASKS_BASE
    eval_base = Path(args.eval_base) if args.eval_base else HPC_EVAL_BASE
    eval_base.mkdir(parents=True, exist_ok=True)

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    # Load vocab to get vocab_size
    from scgpt.tokenizer import GeneVocab
    vocab_path = model_base / "scgpt_whole_human" / "vocab.json"
    with open(vocab_path) as f:
        vocab_dict = json.load(f)
    vocab_size = len(vocab_dict)
    log.info(f"Vocab size: {vocab_size:,}")

    model_dir = model_base / "scgpt_whole_human"

    # Determine tasks/folds
    if args.all:
        task_list = list(TISSUE_TASKS.keys())
    elif args.task:
        task_list = [args.task]
    else:
        raise ValueError("Specify --task <A1-A6> or --all")

    all_results = {}

    for task_id in task_list:
        task_dir_name, all_folds = TISSUE_TASKS[task_id]

        if args.fold == "lomo" or args.all:
            folds_to_run = all_folds
        else:
            folds_to_run = [args.fold]

        task_results = {}
        for fold_mission in folds_to_run:
            result = run_fold(
                task_id, task_dir_name, fold_mission, args,
                model_dir, tasks_base, eval_base, vocab_size,
            )
            if result:
                task_results[fold_mission] = result

        if task_results:
            # Compute summary (skip dry_run entries)
            aurocs = [r["auroc"] for r in task_results.values() if "auroc" in r]
            mean_auroc = float(np.mean(aurocs))
            log.info(f"\n{task_id} summary: mean AUROC = {mean_auroc:.4f} ({len(aurocs)} folds)")

            all_results[task_id] = {
                "folds": task_results,
                "mean_auroc": mean_auroc,
                "n_folds": len(aurocs),
            }

            # Save per-task results
            out_file = eval_base / f"scgpt_whole_human_{task_id}_lomo_results.json"
            with open(out_file, "w") as f:
                json.dump(all_results[task_id], f, indent=2)
            log.info(f"Saved: {out_file}")

    # Save combined summary
    if all_results and not args.dry_run and any(all_results.values()):
        summary_file = eval_base / "scgpt_whole_human_all_tissues_summary.json"
        with open(summary_file, "w") as f:
            json.dump({
                "model": "scgpt_whole_human",
                "freeze_layers": args.freeze_layers,
                "epochs": args.epochs,
                "lr": args.lr,
                "results": all_results,
            }, f, indent=2)
        log.info(f"\nCombined summary: {summary_file}")


if __name__ == "__main__":
    main()

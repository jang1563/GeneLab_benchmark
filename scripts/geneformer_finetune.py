"""
Geneformer Fine-Tuning Pipeline for Spaceflight Detection (A4 Thymus LOMO)

Fine-tunes Geneformer V1-10M for binary classification (Flight vs Ground Control)
using Leave-One-Mission-Out (LOMO) cross-validation folds.

Architecture:
  - Base: Geneformer-V1-10M (BertForMaskedLM → weights extracted)
  - Head: BertForSequenceClassification (2 classes: FLT=1, GC=0)
  - Strategy: Full fine-tuning OR head-only (--freeze-layers) for small-n datasets

Hardware targets:
  - Development: MacBook Air Apple Silicon (MPS, float32, batch_size=4)
  - Production: HPC A40/A100 (CUDA, AMP mixed precision, batch_size=16-32)

Freezing strategy (important for small n):
  --freeze-layers 0   : Full fine-tuning (all 6 BERT layers + head) — risky for n<100
  --freeze-layers 4   : Fine-tune top 2 layers + head (recommended for n~50)
  --freeze-layers 6   : Head-only (classifier) fine-tuning (most regularized)

Usage:
  # Fine-tune single fold (MPS, development) — head-only for small n
  python scripts/geneformer_finetune.py --fold RR-9 --device mps --epochs 10 --freeze-layers 6

  # Fine-tune all LOMO folds (dry-run, check shapes)
  python scripts/geneformer_finetune.py --fold lomo --dry-run

  # Fine-tune for HPC (CUDA, mixed precision)
  python scripts/geneformer_finetune.py --fold lomo --device cuda --epochs 10 --batch-size 16 --freeze-layers 4

Output:
  evaluation/geneformer_v1_A4_lomo_results.json
  tasks/A4_thymus_lomo/fold_{fold}/geneformer_tokens/v1/checkpoints/
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
import torch
from datasets import load_from_disk
from sklearn.metrics import roc_auc_score
from torch.optim import AdamW
from torch.utils.data import DataLoader
from transformers import (
    BertConfig,
    BertForSequenceClassification,
    get_linear_schedule_with_warmup,
)
from huggingface_hub import hf_hub_download

# Default random seed for reproducibility
DEFAULT_SEED = 42

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

ROOT = Path(__file__).parent.parent
TASKS_DIR = ROOT / "tasks"
EVAL_DIR = ROOT / "evaluation"
HF_REPO = "ctheodoris/Geneformer"

# Mouse-Geneformer local path (Cayuga)
MOUSE_GF_MODEL_DIR = Path(
    os.environ.get(
        "MOUSE_GF_MODEL_DIR",
        "/athena/masonlab/scratch/users/jak4013/huggingface/benchmark/GeneLab_benchmark/models/mouse_gf_base",
    )
)

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

# ─── Model Loading ─────────────────────────────────────────────────────────────

def load_geneformer_for_classification(
    n_labels: int = 2,
    model_version: str = "v1",
    mouse_gf_model_dir: Optional[Path] = None,
) -> BertForSequenceClassification:
    """
    Load Geneformer pretrained weights into a BertForSequenceClassification model.
    The MLM head is discarded; a new 2-class classification head is added.

    Supports:
      v1/v2: Download from HuggingFace (ctheodoris/Geneformer)
      mouse_gf: Load from local Cayuga path (MOUSE_GF_MODEL_DIR)
    """
    if model_version == "mouse_gf":
        # Load Mouse-Geneformer from local path
        model_dir = mouse_gf_model_dir or MOUSE_GF_MODEL_DIR
        log.info(f"Loading Mouse-Geneformer from {model_dir}")
        config_path = model_dir / "config.json"
        weights_path = model_dir / "pytorch_model.bin"
        missing = [str(p) for p in (config_path, weights_path) if not p.exists()]
        if missing:
            raise FileNotFoundError(
                "Mouse-Geneformer model files not found:\n"
                + "\n".join(f"  - {m}" for m in missing)
                + "\nSet --mouse-gf-model-dir or MOUSE_GF_MODEL_DIR to the directory containing these files."
            )
        with open(config_path) as f:
            config_dict = json.load(f)

        config_dict["architectures"] = ["BertForSequenceClassification"]
        config_dict["num_labels"] = n_labels
        config_dict["problem_type"] = "single_label_classification"
        config = BertConfig(**{k: v for k, v in config_dict.items() if k != "architectures"})
        config.num_labels = n_labels

        state_dict = torch.load(str(weights_path), map_location="cpu")
    else:
        model_dirs = {
            "v1": "Geneformer-V1-10M",
            "v2": "Geneformer-V2-104M",
        }
        model_dir = model_dirs[model_version]

        log.info(f"Loading Geneformer {model_version} config...")
        config_path = hf_hub_download(repo_id=HF_REPO, filename=f"{model_dir}/config.json")

        with open(config_path) as f:
            config_dict = json.load(f)

        config_dict["architectures"] = ["BertForSequenceClassification"]
        config_dict["num_labels"] = n_labels
        config_dict["problem_type"] = "single_label_classification"
        config = BertConfig(**{k: v for k, v in config_dict.items() if k != "architectures"})
        config.num_labels = n_labels

        log.info(f"Loading pretrained weights from {model_dir}...")
        try:
            weights_path = hf_hub_download(
                repo_id=HF_REPO,
                filename=f"{model_dir}/model.safetensors"
            )
            from safetensors.torch import load_file
            state_dict = load_file(weights_path)
        except Exception:
            weights_path = hf_hub_download(
                repo_id=HF_REPO,
                filename=f"{model_dir}/pytorch_model.bin"
            )
            state_dict = torch.load(weights_path, map_location="cpu")

    # Initialize classification model and load compatible weights
    model = BertForSequenceClassification(config)

    # Filter out MLM-specific layers (cls.predictions.*) — not compatible with classification
    compatible_keys = {
        k: v for k, v in state_dict.items()
        if not k.startswith("cls.predictions") and k in model.state_dict()
    }
    missing, unexpected = model.load_state_dict(compatible_keys, strict=False)

    log.info(f"  Loaded {len(compatible_keys)} layers from pretrained weights")
    if missing:
        log.info(f"  Missing (new classification head): {missing[:5]}...")
    if unexpected:
        log.info(f"  Unexpected (MLM head, discarded): {unexpected[:3]}...")

    return model


def freeze_bert_layers(model: BertForSequenceClassification, n_freeze: int) -> None:
    """
    Freeze the bottom n_freeze BERT encoder layers.
    Geneformer-V1-10M has 6 encoder layers (0-5).

    n_freeze=0: Full fine-tuning (all layers trainable)
    n_freeze=4: Fine-tune top 2 layers + pooler + head
    n_freeze=6: Head-only (most regularized, best for n<50)

    Always keeps classifier head trainable.
    """
    if n_freeze <= 0:
        log.info("  Layer freezing: disabled (full fine-tuning)")
        return

    # Freeze embeddings
    for param in model.bert.embeddings.parameters():
        param.requires_grad = False

    # Freeze bottom n_freeze encoder layers
    for i in range(min(n_freeze, len(model.bert.encoder.layer))):
        for param in model.bert.encoder.layer[i].parameters():
            param.requires_grad = False

    n_total = len(model.bert.encoder.layer)
    n_trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
    n_frozen = sum(p.numel() for p in model.parameters() if not p.requires_grad)
    log.info(f"  Frozen: {n_freeze}/{n_total} BERT layers + embeddings ({n_frozen:,} params)")
    log.info(f"  Trainable: {n_trainable:,} params (top {n_total - n_freeze} layers + head)")


# ─── DataLoader ────────────────────────────────────────────────────────────────

def collate_fn(batch: list[dict], pad_token_id: int = 0, max_length: int = 2048):
    """
    Pad sequences to the same length within a batch and create attention masks.
    """
    input_ids_list = [torch.tensor(item["input_ids"], dtype=torch.long) for item in batch]
    labels = torch.tensor([item["label"] for item in batch], dtype=torch.long)

    # Pad to max length in batch (not global max_length for efficiency)
    max_len = min(max(len(x) for x in input_ids_list), max_length)
    padded_input_ids = torch.zeros(len(batch), max_len, dtype=torch.long)
    attention_mask = torch.zeros(len(batch), max_len, dtype=torch.long)

    for i, ids in enumerate(input_ids_list):
        seq_len = min(len(ids), max_len)
        padded_input_ids[i, :seq_len] = ids[:seq_len]
        attention_mask[i, :seq_len] = 1

    return {
        "input_ids": padded_input_ids,
        "attention_mask": attention_mask,
        "labels": labels,
    }


def make_dataloader(
    dataset,
    batch_size: int,
    shuffle: bool = True,
    pad_token_id: int = 0,
    device_type: str = "cpu",
):
    # num_workers=0 required for MPS (macOS multiprocessing issues)
    # num_workers=4 recommended for CUDA (enables prefetching)
    nw = 0 if device_type in ("mps", "cpu") else 4
    pin = device_type == "cuda"
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        collate_fn=lambda batch: collate_fn(batch, pad_token_id=pad_token_id),
        num_workers=nw,
        pin_memory=pin,
    )


# ─── Training ──────────────────────────────────────────────────────────────────

def train_epoch(model, dataloader, optimizer, scheduler, device, scaler=None):
    model.train()
    total_loss = 0
    n_batches = 0
    use_amp = scaler is not None

    for batch in dataloader:
        batch = {k: v.to(device) for k, v in batch.items()}
        optimizer.zero_grad()

        if use_amp:
            with torch.cuda.amp.autocast():
                outputs = model(**batch)
                loss = outputs.loss
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            scaler.step(optimizer)
            scaler.update()
        else:
            outputs = model(**batch)
            loss = outputs.loss
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

        scheduler.step()
        total_loss += loss.item()
        n_batches += 1

    return total_loss / n_batches


def evaluate(model, dataloader, device):
    model.eval()
    all_probs = []
    all_labels = []

    with torch.no_grad():
        for batch in dataloader:
            batch = {k: v.to(device) for k, v in batch.items()}
            labels = batch.pop("labels")
            outputs = model(**batch)
            probs = torch.softmax(outputs.logits, dim=-1)[:, 1]  # P(Flight)
            all_probs.extend(probs.cpu().numpy().tolist())
            all_labels.extend(labels.cpu().numpy().tolist())

    if len(set(all_labels)) < 2:
        log.warning("  Only one class in labels — AUROC undefined")
        return float("nan"), all_probs, all_labels

    auroc = roc_auc_score(all_labels, all_probs)
    return auroc, all_probs, all_labels


# ─── Fine-tune Fold ────────────────────────────────────────────────────────────

def finetune_fold(
    fold_dir: Path,
    device: torch.device,
    model_version: str = "v1",
    mouse_gf_model_dir: Optional[Path] = None,
    n_epochs: int = 5,
    batch_size: int = 4,
    lr: float = 2e-5,
    warmup_ratio: float = 0.1,
    freeze_layers: int = 0,
    seed: int = DEFAULT_SEED,
    overwrite: bool = False,
    dry_run: bool = False,
    save_checkpoint: bool = True,
) -> dict:
    """
    Fine-tune Geneformer on one LOMO fold and evaluate on held-out test set.
    """
    token_dir = fold_dir / "geneformer_tokens" / model_version
    checkpoint_dir = token_dir / "checkpoints"
    result_path = token_dir / "finetune_result.json"

    if result_path.exists() and not overwrite:
        log.info(f"  Skipping {fold_dir.name} (already fine-tuned; use --overwrite)")
        with open(result_path) as f:
            return json.load(f)

    if not (token_dir / "train_dataset").exists():
        log.error(f"  Tokenized data not found: {token_dir}")
        log.error("  Run geneformer_tokenize.py first.")
        return {"status": "missing_tokens", "fold": fold_dir.name}

    # Set random seed for reproducibility
    torch.manual_seed(seed)
    np.random.seed(seed)
    if device.type == "cuda":
        torch.cuda.manual_seed_all(seed)

    log.info(f"\n=== Fine-tuning fold: {fold_dir.name} ===")
    log.info(f"  Device: {device} | epochs: {n_epochs} | batch_size: {batch_size} | lr: {lr} | freeze_layers: {freeze_layers} | seed: {seed}")

    # Load tokenized datasets
    train_ds = load_from_disk(str(token_dir / "train_dataset"))
    test_ds = load_from_disk(str(token_dir / "test_dataset"))

    log.info(f"  Train: {len(train_ds)} samples | Test: {len(test_ds)} samples")
    log.info(f"  Train label dist: {sum(train_ds['label'])} FLT / {len(train_ds) - sum(train_ds['label'])} GC")
    log.info(f"  Test  label dist: {sum(test_ds['label'])} FLT / {len(test_ds) - sum(test_ds['label'])} GC")

    if dry_run:
        log.info("  [DRY RUN] Model shape check only...")
        model = load_geneformer_for_classification(
            n_labels=2,
            model_version=model_version,
            mouse_gf_model_dir=mouse_gf_model_dir,
        )
        n_params = sum(p.numel() for p in model.parameters())
        log.info(f"  Model parameters: {n_params:,}")
        return {"status": "dry_run", "fold": fold_dir.name, "n_params": n_params}

    # Load model + apply layer freezing
    model = load_geneformer_for_classification(
        n_labels=2,
        model_version=model_version,
        mouse_gf_model_dir=mouse_gf_model_dir,
    )
    freeze_bert_layers(model, n_freeze=freeze_layers)
    model = model.to(device)

    # DataLoaders (device-aware num_workers and pin_memory)
    train_loader = make_dataloader(train_ds, batch_size=batch_size, shuffle=True, device_type=device.type)
    test_loader = make_dataloader(test_ds, batch_size=batch_size, shuffle=False, device_type=device.type)

    # Optimizer: exclude bias and LayerNorm from weight decay
    no_decay = {"bias", "LayerNorm.weight"}
    param_groups = [
        {
            "params": [p for n, p in model.named_parameters()
                       if not any(nd in n for nd in no_decay) and p.requires_grad],
            "weight_decay": 0.01,
        },
        {
            "params": [p for n, p in model.named_parameters()
                       if any(nd in n for nd in no_decay) and p.requires_grad],
            "weight_decay": 0.0,
        },
    ]
    optimizer = AdamW(param_groups, lr=lr)

    total_steps = len(train_loader) * n_epochs
    warmup_steps = int(warmup_ratio * total_steps)
    scheduler = get_linear_schedule_with_warmup(
        optimizer, num_warmup_steps=warmup_steps, num_training_steps=total_steps
    )

    # Mixed precision scaler for CUDA (not supported on MPS)
    scaler = torch.cuda.amp.GradScaler() if device.type == "cuda" else None
    if scaler:
        log.info("  Mixed precision (AMP) enabled")

    log.info(f"  Training steps: {total_steps} | Warmup: {warmup_steps}")

    # Training loop
    best_test_auroc = float("-inf")
    best_probs = []
    best_labels = []
    last_probs = []
    last_labels = []
    history = []
    start_time = time.time()

    for epoch in range(1, n_epochs + 1):
        epoch_start = time.time()
        train_loss = train_epoch(model, train_loader, optimizer, scheduler, device, scaler=scaler)
        test_auroc, test_probs, test_labels = evaluate(model, test_loader, device)
        last_probs, last_labels = test_probs, test_labels
        epoch_time = time.time() - epoch_start

        log.info(f"  Epoch {epoch}/{n_epochs} | loss={train_loss:.4f} | test_auroc={test_auroc:.4f} | {epoch_time:.1f}s")
        history.append({"epoch": epoch, "train_loss": train_loss, "test_auroc": test_auroc})

        if np.isfinite(test_auroc) and test_auroc > best_test_auroc:
            best_test_auroc = test_auroc
            best_probs = test_probs
            best_labels = test_labels

            # Save best checkpoint
            if save_checkpoint:
                checkpoint_dir.mkdir(parents=True, exist_ok=True)
                model.save_pretrained(str(checkpoint_dir / "best_model"))
                log.info(f"  Saved best checkpoint (auroc={best_test_auroc:.4f})")

    total_time = time.time() - start_time
    if not best_probs:
        best_probs = last_probs
        best_labels = last_labels
    if not np.isfinite(best_test_auroc):
        if best_labels and len(set(best_labels)) >= 2:
            best_test_auroc = float(roc_auc_score(best_labels, best_probs))
        else:
            best_test_auroc = float("nan")
    best_auroc_str = f"{best_test_auroc:.4f}" if np.isfinite(best_test_auroc) else "nan"
    log.info(f"  Done. Best AUROC: {best_auroc_str} | Total time: {total_time:.0f}s")

    # Bootstrap CI
    n_bootstrap = 1000
    rng = np.random.default_rng(seed)
    boot_aurocs = []
    best_probs_arr = np.array(best_probs)
    best_labels_arr = np.array(best_labels)
    if len(best_labels_arr) > 0 and len(np.unique(best_labels_arr)) >= 2:
        for _ in range(n_bootstrap):
            idx = rng.integers(0, len(best_labels_arr), size=len(best_labels_arr))
            if len(set(best_labels_arr[idx])) < 2:
                continue
            boot_aurocs.append(roc_auc_score(best_labels_arr[idx], best_probs_arr[idx]))

    ci_low = float(np.percentile(boot_aurocs, 2.5)) if boot_aurocs else float("nan")
    ci_high = float(np.percentile(boot_aurocs, 97.5)) if boot_aurocs else float("nan")

    result = {
        "status": "ok",
        "fold": fold_dir.name,
        "model_version": model_version,
        "device": str(device),
        "n_epochs": n_epochs,
        "batch_size": batch_size,
        "lr": lr,
        "best_test_auroc": best_test_auroc,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "n_train": len(train_ds),
        "n_test": len(test_ds),
        "n_flight_test": int(sum(test_ds["label"])),
        "n_ground_test": int(len(test_ds) - sum(test_ds["label"])),
        "training_time_sec": total_time,
        "history": history,
        "best_probs": best_probs,
        "best_labels": best_labels,
    }

    if save_checkpoint:
        with open(result_path, "w") as f:
            json.dump({k: v for k, v in result.items() if k not in ("best_probs", "best_labels")}, f, indent=2)

    return result


# ─── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--task", default="A4", help="Task ID (e.g. A1, A4, A6)")
    parser.add_argument("--fold", default="RR-9",
                        help="Fold name (e.g. 'RR-9', 'RR-6', 'MHU-1', 'all'). Default=RR-9")
    parser.add_argument("--model-version", default="mouse_gf", choices=["v1", "v2", "mouse_gf"])
    parser.add_argument("--device", default="auto",
                        help="Device: 'mps', 'cuda', 'cpu', or 'auto'. Default=auto")
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--lr", type=float, default=2e-5)
    parser.add_argument("--freeze-layers", type=int, default=0,
                        help="Freeze bottom N BERT encoder layers (0=full fine-tune, 4=top-2+head, 6=head-only). "
                             "Recommended: 6 for n<50, 4 for n~50-100, 0 for n>200.")
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="Load model + check shapes only, no training")
    parser.add_argument("--no-save", action="store_true", help="Don't save checkpoint to disk")
    parser.add_argument("--task-dir", default=None,
                        help="Explicit task directory name under tasks/ (recommended for A1 variants)")
    parser.add_argument("--mouse-gf-model-dir", default=str(MOUSE_GF_MODEL_DIR),
                        help="Directory containing Mouse-Geneformer config.json and pytorch_model.bin")
    args = parser.parse_args()

    # Device selection
    if args.device == "auto":
        if torch.backends.mps.is_available():
            device = torch.device("mps")
        elif torch.cuda.is_available():
            device = torch.device("cuda")
        else:
            device = torch.device("cpu")
    else:
        device = torch.device(args.device)

    log.info(f"Device: {device}")
    if device.type == "mps":
        log.info("  MacBook Apple Silicon (MPS) — using float32 (4-bit quant not supported)")

    # Find folds (dynamic task directory lookup)
    try:
        task_dir = resolve_task_dir(args.task, task_dir_name=args.task_dir)
    except (FileNotFoundError, ValueError) as e:
        log.error(str(e))
        return
    log.info(f"Task directory: {task_dir.name}")

    if args.fold in ("all", "lomo"):
        # Exclude held-out fold (fold_RR-23_holdout) from LOMO CV
        fold_dirs = sorted([d for d in task_dir.iterdir()
                            if d.is_dir() and d.name.startswith("fold_")
                            and "holdout" not in d.name])
    else:
        if args.fold.startswith("fold_"):
            candidates = [task_dir / args.fold]
        else:
            candidates = [
                task_dir / f"fold_{args.fold}_test",
                task_dir / f"fold_{args.fold}_holdout",
                task_dir / f"fold_{args.fold}",
            ]
        fold_dirs = [d for d in candidates if d.exists() and d.is_dir()]
        # Guard: warn if the held-out evaluation fold was matched
        if any("holdout" in d.name for d in fold_dirs):
            log.warning("⚠  Matched the held-out evaluation fold (fold_RR-23_holdout).")
            log.warning("   This fold is reserved for final evaluation — do NOT use it for LOMO CV training.")
            log.warning("   Use --fold lomo to run LOMO CV on all non-held-out folds.")

    if not fold_dirs:
        log.error(f"No fold found matching '{args.fold}' in {task_dir}")
        return

    log.info(f"Folds to fine-tune: {[d.name for d in fold_dirs]}")

    # Fine-tune each fold
    all_results = []
    for fold_dir in fold_dirs:
        result = finetune_fold(
            fold_dir=fold_dir,
            device=device,
            model_version=args.model_version,
            mouse_gf_model_dir=Path(args.mouse_gf_model_dir),
            n_epochs=args.epochs,
            batch_size=args.batch_size,
            lr=args.lr,
            freeze_layers=args.freeze_layers,
            seed=args.seed,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
            save_checkpoint=not args.no_save,
        )
        all_results.append(result)

    # Print summary
    print("\n=== Geneformer Fine-Tuning Summary ===")
    print(f"Task: {args.task} | Model: Geneformer-{args.model_version} | Device: {device}")
    print(f"{'Fold':<20} {'Status':<10} {'AUROC':<8} {'95% CI':<20} {'n_test'}")
    print("-" * 70)
    for r in all_results:
        if r.get("status") == "ok":
            ci = f"[{r['ci_low']:.3f}, {r['ci_high']:.3f}]"
            print(f"{r['fold']:<20} {'OK':<10} {r['best_test_auroc']:.4f}   {ci:<20} {r['n_test']}")
        else:
            print(f"{r['fold']:<20} {r.get('status', '?'):<10}")

    # Save full LOMO results if all folds done
    if args.fold in ("all", "lomo") and not args.dry_run:
        all_aurocs = [r["best_test_auroc"] for r in all_results if r.get("status") == "ok"]
        if all_aurocs:
            mean_auroc = float(np.mean(all_aurocs))
            std_auroc = float(np.std(all_aurocs))
            EVAL_DIR.mkdir(exist_ok=True)
            # Strip large per-sample arrays (best_probs/best_labels) from fold summaries
            fold_summaries = [
                {k: v for k, v in r.items() if k not in ("best_probs", "best_labels")}
                for r in all_results if r.get("status") == "ok"
            ]
            lomo_result = {
                "task": args.task,
                "model": f"Geneformer-{args.model_version}",
                "freeze_layers": args.freeze_layers,
                "seed": args.seed,
                "mean_auroc": mean_auroc,
                "std_auroc": std_auroc,
                "n_folds": len(all_aurocs),
                "fold_results": fold_summaries,
            }
            out_path = EVAL_DIR / f"geneformer_{args.model_version}_{args.task}_lomo_results.json"
            log.info(f"  Model: {args.model_version} | Task: {args.task}")
            with open(out_path, "w") as f:
                json.dump(lomo_result, f, indent=2)
            log.info(f"\nLOMO summary saved to {out_path}")
            print(f"\nLOMO Mean AUROC: {mean_auroc:.4f} ± {std_auroc:.4f}")


if __name__ == "__main__":
    main()

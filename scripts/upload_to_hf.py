#!/usr/bin/env python3
"""
upload_to_hf.py — Upload GeneLab benchmark feature matrices to HuggingFace Dataset

Uploads train_X.csv and test_X.csv for GO tasks to a HuggingFace Dataset repo.
Labels, metadata, and fold structure (train_y.csv, test_y.csv, fold_info.json, etc.)
are kept in the GitHub repo. Feature matrices are hosted on HuggingFace due to size.

HuggingFace repo: jang1563/genelab-benchmark

Usage:
    python scripts/upload_to_hf.py --task A5 --dry-run
    python scripts/upload_to_hf.py --task all
    python scripts/upload_to_hf.py --task A2 A4 A5
    python scripts/upload_to_hf.py --card-only   # Upload dataset card (README.md) only

Authentication:
    export HF_TOKEN="hf_..."
    python scripts/upload_to_hf.py --task A5
"""

import os
import sys
import argparse
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
TASKS_DIR = BASE_DIR / "tasks"

HF_REPO_ID = "jang1563/genelab-benchmark"

# GO tasks only (A1/A3 NO-GO excluded from primary benchmark)
GO_TASKS = [
    "A2_gastrocnemius_lomo",
    "A4_thymus_lomo",
    "A5_skin_lomo",
    "A6_eye_lomo",
]

UPLOAD_FILES = ["train_X.csv", "test_X.csv"]
CARD_SRC = BASE_DIR / "docs" / "hf_dataset_card.md"


def upload_card(api, dry_run: bool = False) -> None:
    """Upload docs/hf_dataset_card.md as README.md to the HF dataset repo."""
    if not CARD_SRC.exists():
        print(f"  [SKIP] Dataset card not found: {CARD_SRC}")
        return
    tag = "[DRY-RUN] " if dry_run else ""
    print(f"  {tag}Upload: README.md (dataset card, {CARD_SRC.stat().st_size // 1024} KB)")
    if not dry_run:
        api.upload_file(
            path_or_fileobj=str(CARD_SRC),
            path_in_repo="README.md",
            repo_id=HF_REPO_ID,
            repo_type="dataset",
        )


def upload_task(api, task_name: str, dry_run: bool = False) -> int:
    """Upload all fold feature matrices for one task. Returns file count."""
    task_dir = TASKS_DIR / task_name
    if not task_dir.exists():
        print(f"  [SKIP] {task_name}: directory not found")
        return 0

    n_uploaded = 0
    for fold_dir in sorted(task_dir.glob("fold_*")):
        for fname in UPLOAD_FILES:
            fpath = fold_dir / fname
            if not fpath.exists():
                print(f"  [SKIP] {fpath.relative_to(BASE_DIR)}: not found (gitignored?)")
                continue
            size_mb = fpath.stat().st_size / 1e6
            hf_path = f"{task_name}/{fold_dir.name}/{fname}"
            tag = "[DRY-RUN] " if dry_run else ""
            print(f"  {tag}Upload: {hf_path} ({size_mb:.0f} MB)")
            if not dry_run:
                api.upload_file(
                    path_or_fileobj=str(fpath),
                    path_in_repo=hf_path,
                    repo_id=HF_REPO_ID,
                    repo_type="dataset",
                )
            n_uploaded += 1

    return n_uploaded


def main():
    parser = argparse.ArgumentParser(
        description="Upload GeneLab benchmark feature matrices to HuggingFace Dataset"
    )
    parser.add_argument(
        "--task",
        nargs="+",
        default=["all"],
        help=(
            "Task(s) to upload. Use 'all' for all GO tasks, or specify task names "
            f"(e.g., A2_gastrocnemius_lomo). GO tasks: {', '.join(GO_TASKS)}"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print files that would be uploaded without actually uploading",
    )
    parser.add_argument(
        "--repo",
        default=HF_REPO_ID,
        help=f"HuggingFace dataset repo ID (default: {HF_REPO_ID})",
    )
    parser.add_argument(
        "--card-only",
        action="store_true",
        help="Upload only the dataset card (README.md) to HuggingFace, skip feature matrices",
    )
    parser.add_argument(
        "--skip-card",
        action="store_true",
        help="Skip dataset card upload (default: card is uploaded with --task all)",
    )
    args = parser.parse_args()

    if args.dry_run:
        print("=== DRY-RUN MODE (no files will be uploaded) ===\n")

    # Card-only mode: upload README.md and exit
    if args.card_only:
        try:
            from huggingface_hub import HfApi
        except ImportError:
            print("[ERROR] huggingface_hub not installed. Install: pip install huggingface_hub")
            sys.exit(1)
        token = os.environ.get("HF_TOKEN")
        if not token and not args.dry_run:
            try:
                from huggingface_hub import get_token
                token = get_token()
            except Exception:
                pass
        if not token and not args.dry_run:
            print("[ERROR] No HuggingFace token found. export HF_TOKEN='hf_...' or huggingface-cli login")
            sys.exit(1)
        api = HfApi(token=token) if not args.dry_run else None
        print(f"\nUploading dataset card → {args.repo}")
        print("-" * 50)
        upload_card(api, dry_run=args.dry_run)
        return

    # Resolve task list
    if "all" in args.task:
        tasks = GO_TASKS
    else:
        # Allow short names like "A5" as well as full names
        tasks = []
        for t in args.task:
            if t in GO_TASKS:
                tasks.append(t)
            else:
                # Try prefix match
                matched = [name for name in GO_TASKS if name.startswith(t)]
                if matched:
                    tasks.extend(matched)
                else:
                    print(f"[WARN] Unknown task: {t}. Available: {GO_TASKS}")
        if not tasks:
            print("No valid tasks specified. Exiting.")
            sys.exit(1)

    # HuggingFace API (lazy import so script is usable without huggingface_hub installed)
    try:
        from huggingface_hub import HfApi
    except ImportError:
        print("[ERROR] huggingface_hub not installed.")
        print("  Install: pip install huggingface_hub")
        sys.exit(1)

    token = os.environ.get("HF_TOKEN")
    if not token and not args.dry_run:
        # Fall back to huggingface-cli cached login
        try:
            from huggingface_hub import get_token
            token = get_token()
        except Exception:
            pass
    if not token and not args.dry_run:
        print("[ERROR] No HuggingFace token found.")
        print("  Option 1: export HF_TOKEN='hf_...'")
        print("  Option 2: huggingface-cli login")
        sys.exit(1)

    api = HfApi(token=token) if not args.dry_run else None
    repo_id = args.repo

    # Ensure repo exists (create if needed)
    if not args.dry_run:
        api.create_repo(repo_id=repo_id, repo_type="dataset", private=True, exist_ok=True)

    total = 0
    for task in tasks:
        print(f"\nUploading {task} → {repo_id}")
        print("-" * 50)
        n = upload_task(api, task, dry_run=args.dry_run)
        total += n
        print(f"  → {n} files")

    # Upload dataset card with full uploads (unless skipped)
    if not args.skip_card:
        print(f"\nUploading dataset card → {repo_id}")
        print("-" * 50)
        upload_card(api, dry_run=args.dry_run)

    print(f"\nTotal: {total} feature files {'(dry-run)' if args.dry_run else 'uploaded'}.")


if __name__ == "__main__":
    main()

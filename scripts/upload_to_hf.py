#!/usr/bin/env python3
"""
upload_to_hf.py — Upload GeneLab benchmark fold packages to HuggingFace Dataset

Uploads self-contained GO task folds to a HuggingFace Dataset repo. Each public
fold includes feature matrices, labels, metadata, fold_info.json, and the
training-only selected gene list so the HF repo can be reviewed and used without
cross-referencing the GitHub checkout for labels.

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

FOLD_PACKAGE_FILES = [
    "train_X.csv",
    "test_X.csv",
    "train_y.csv",
    "test_y.csv",
    "train_meta.csv",
    "test_meta.csv",
    "fold_info.json",
    "selected_genes.txt",
]
TASK_ROOT_FILES = ["task_info.json"]
STALE_REMOTE_DIRS = [
    "A6_eye_lomo/fold_TBD_test",
]
CARD_SRC = BASE_DIR / "docs" / "hf_dataset_card.md"
CARD_ASSETS = [
    BASE_DIR / "docs" / "assets" / "hf_benchmark_summary.png",
]


def upload_card(api, repo_id: str, dry_run: bool = False) -> None:
    """Upload docs/hf_dataset_card.md and its rendered assets to the HF dataset repo."""
    if not CARD_SRC.exists():
        print(f"  [SKIP] Dataset card not found: {CARD_SRC}")
        return
    tag = "[DRY-RUN] " if dry_run else ""
    print(f"  {tag}Upload: README.md (dataset card, {CARD_SRC.stat().st_size // 1024} KB)")
    if not dry_run:
        api.upload_file(
            path_or_fileobj=str(CARD_SRC),
            path_in_repo="README.md",
            repo_id=repo_id,
            repo_type="dataset",
        )
    for asset in CARD_ASSETS:
        if not asset.exists():
            print(f"  [SKIP] Card asset not found: {asset.relative_to(BASE_DIR)}")
            continue
        asset_repo_path = Path("assets") / asset.name
        size_kb = asset.stat().st_size / 1024
        print(f"  {tag}Upload: {asset_repo_path} ({size_kb:.1f} KB)")
        if not dry_run:
            api.upload_file(
                path_or_fileobj=str(asset),
                path_in_repo=str(asset_repo_path),
                repo_id=repo_id,
                repo_type="dataset",
            )


def report_package_file(fpath: Path, hf_path: str, dry_run: bool = False) -> bool:
    """Report one package file. Returns True when the file is present locally."""
    if not fpath.exists():
        print(f"  [SKIP] {fpath.relative_to(BASE_DIR)}: not found")
        return False

    size_mb = fpath.stat().st_size / 1e6
    size_label = f"{size_mb:.1f} MB" if size_mb >= 1 else f"{fpath.stat().st_size / 1024:.1f} KB"
    tag = "[DRY-RUN] " if dry_run else ""
    print(f"  {tag}Upload: {hf_path} ({size_label})")
    return True


def upload_task(api, task_name: str, repo_id: str, dry_run: bool = False) -> int:
    """Upload all public fold package files for one task. Returns file count."""
    task_dir = TASKS_DIR / task_name
    if not task_dir.exists():
        print(f"  [SKIP] {task_name}: directory not found")
        return 0

    n_uploaded = 0
    for fname in TASK_ROOT_FILES:
        fpath = task_dir / fname
        hf_path = f"{task_name}/{fname}"
        if report_package_file(fpath, hf_path, dry_run=dry_run):
            n_uploaded += 1

    for fold_dir in sorted(task_dir.glob("fold_*")):
        for fname in FOLD_PACKAGE_FILES:
            fpath = fold_dir / fname
            hf_path = f"{task_name}/{fold_dir.name}/{fname}"
            if report_package_file(fpath, hf_path, dry_run=dry_run):
                n_uploaded += 1

    if n_uploaded and not dry_run:
        allow_patterns = TASK_ROOT_FILES + [
            f"fold_*/{fname}" for fname in FOLD_PACKAGE_FILES
        ]
        api.upload_folder(
            folder_path=str(task_dir),
            path_in_repo=task_name,
            repo_id=repo_id,
            repo_type="dataset",
            allow_patterns=allow_patterns,
            commit_message=f"Upload {task_name} self-contained fold package",
        )

    return n_uploaded


def prune_stale_remote_dirs(api, repo_id: str, dry_run: bool = False) -> None:
    """Remove known stale HF paths that would make the public package ambiguous."""
    for remote_dir in STALE_REMOTE_DIRS:
        tag = "[DRY-RUN] " if dry_run else ""
        print(f"  {tag}Delete stale remote folder: {remote_dir}/")
        if dry_run:
            continue
        try:
            api.delete_folder(
                path_in_repo=remote_dir,
                repo_id=repo_id,
                repo_type="dataset",
                commit_message=f"Remove stale {remote_dir}",
            )
        except Exception as exc:
            print(f"  [WARN] Could not delete {remote_dir}/; it may already be absent: {exc}")


def main():
    parser = argparse.ArgumentParser(
        description="Upload GeneLab benchmark self-contained fold packages to HuggingFace Dataset"
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
        "--private",
        action="store_true",
        help="Create the HF dataset repo as private if it does not already exist",
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
    parser.add_argument(
        "--skip-prune",
        action="store_true",
        help="Skip deletion of known stale remote folders such as A6 fold_TBD_test",
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
        repo_id = args.repo
        print(f"\nUploading dataset card → {repo_id}")
        print("-" * 50)
        upload_card(api, repo_id, dry_run=args.dry_run)
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
        api.create_repo(repo_id=repo_id, repo_type="dataset", private=args.private, exist_ok=True)

    if not args.skip_prune:
        print(f"\nPruning stale remote paths → {repo_id}")
        print("-" * 50)
        prune_stale_remote_dirs(api, repo_id, dry_run=args.dry_run)

    total = 0
    for task in tasks:
        print(f"\nUploading {task} → {repo_id}")
        print("-" * 50)
        n = upload_task(api, task, repo_id, dry_run=args.dry_run)
        total += n
        print(f"  → {n} files")

    # Upload dataset card with full uploads (unless skipped)
    if not args.skip_card:
        print(f"\nUploading dataset card → {repo_id}")
        print("-" * 50)
        upload_card(api, repo_id, dry_run=args.dry_run)

    print(f"\nTotal: {total} package files {'(dry-run)' if args.dry_run else 'uploaded'}.")


if __name__ == "__main__":
    main()

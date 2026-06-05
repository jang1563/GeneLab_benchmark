#!/usr/bin/env python3
"""
download_from_hf.py — Download GeneLab benchmark fold packages from HuggingFace

Downloads self-contained GO task folds from the HuggingFace Dataset repo into
the local tasks/ directory, matching the layout expected by evaluate_submission.py
and run_baselines.py.

HuggingFace repo: jang1563/genelab-benchmark

Usage:
    python scripts/download_from_hf.py --task A5
    python scripts/download_from_hf.py --task all
    python scripts/download_from_hf.py --task A2 A4
    python scripts/download_from_hf.py --task A5 --dry-run

Authentication (optional — public repo can be downloaded without token):
    export HF_TOKEN="hf_..."
"""

import os
import sys
import argparse
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
TASKS_DIR = BASE_DIR / "tasks"

HF_REPO_ID = "jang1563/genelab-benchmark"

GO_TASKS = [
    "A2_gastrocnemius_lomo",
    "A4_thymus_lomo",
    "A5_skin_lomo",
    "A6_eye_lomo",
]

FOLD_PACKAGE_FILES = {
    "train_X.csv",
    "test_X.csv",
    "train_y.csv",
    "test_y.csv",
    "train_meta.csv",
    "test_meta.csv",
    "fold_info.json",
    "selected_genes.txt",
}
TASK_ROOT_FILES = {"task_info.json"}
STALE_REMOTE_DIRS = {"A6_eye_lomo/fold_TBD_test"}


def download_task(task_name: str, repo_id: str, token: str = None, dry_run: bool = False) -> int:
    """Download all public fold package files for one task. Returns file count."""
    from huggingface_hub import hf_hub_download, list_repo_files

    # List available task package files from HF.
    try:
        all_files = list(list_repo_files(repo_id, repo_type="dataset", token=token))
    except Exception as e:
        print(f"  [ERROR] Failed to list HF repo files: {e}")
        return 0

    task_files = []
    for hf_path in all_files:
        if not hf_path.startswith(f"{task_name}/"):
            continue
        if any(hf_path.startswith(f"{stale}/") for stale in STALE_REMOTE_DIRS):
            continue
        parts = hf_path.split("/")
        if len(parts) == 2 and parts[1] in TASK_ROOT_FILES:
            task_files.append(hf_path)
        elif len(parts) == 3 and parts[1].startswith("fold_") and parts[2] in FOLD_PACKAGE_FILES:
            task_files.append(hf_path)

    if not task_files:
        print(f"  [SKIP] No files found in HF repo for {task_name}")
        return 0

    n_downloaded = 0
    for hf_path in sorted(task_files):
        local_path = TASKS_DIR / hf_path

        if local_path.exists():
            print(f"  [SKIP] Already exists: {local_path.relative_to(BASE_DIR)}")
            n_downloaded += 1
            continue

        tag = "[DRY-RUN] " if dry_run else ""
        print(f"  {tag}Download: {hf_path} → {local_path.relative_to(BASE_DIR)}")

        if not dry_run:
            local_path.parent.mkdir(parents=True, exist_ok=True)
            hf_hub_download(
                repo_id=repo_id,
                filename=hf_path,
                repo_type="dataset",
                token=token,
                local_dir=str(TASKS_DIR),
            )

        n_downloaded += 1

    return n_downloaded


def main():
    parser = argparse.ArgumentParser(
        description="Download GeneLab benchmark self-contained fold packages from HuggingFace"
    )
    parser.add_argument(
        "--task",
        nargs="+",
        default=["all"],
        help=(
            "Task(s) to download. Use 'all' for all GO tasks, or specify task names. "
            f"GO tasks: {', '.join(GO_TASKS)}"
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print files that would be downloaded without downloading",
    )
    parser.add_argument(
        "--repo",
        default=HF_REPO_ID,
        help=f"HuggingFace dataset repo ID (default: {HF_REPO_ID})",
    )
    args = parser.parse_args()

    if args.dry_run:
        print("=== DRY-RUN MODE (no files will be downloaded) ===\n")

    # Resolve task list
    if "all" in args.task:
        tasks = GO_TASKS
    else:
        tasks = []
        for t in args.task:
            if t in GO_TASKS:
                tasks.append(t)
            else:
                matched = [name for name in GO_TASKS if name.startswith(t)]
                if matched:
                    tasks.extend(matched)
                else:
                    print(f"[WARN] Unknown task: {t}. Available: {GO_TASKS}")
        if not tasks:
            print("No valid tasks specified. Exiting.")
            sys.exit(1)

    try:
        import huggingface_hub  # noqa: F401
    except ImportError:
        print("[ERROR] huggingface_hub not installed.")
        print("  Install: pip install huggingface_hub")
        sys.exit(1)

    token = os.environ.get("HF_TOKEN")  # Optional for public repo

    total = 0
    repo_id = args.repo
    for task in tasks:
        print(f"\nDownloading {task} ← {repo_id}")
        print("-" * 50)
        n = download_task(task, repo_id, token=token, dry_run=args.dry_run)
        total += n
        print(f"  → {n} files")

    print(f"\nTotal: {total} files {'(dry-run)' if args.dry_run else 'downloaded'}.")


if __name__ == "__main__":
    main()

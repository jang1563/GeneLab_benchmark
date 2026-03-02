#!/usr/bin/env python3
"""
download_from_hf.py — Download GeneLab benchmark feature matrices from HuggingFace

Downloads train_X.csv and test_X.csv for GO tasks from the HuggingFace Dataset
repo into the local tasks/ directory, matching the layout expected by
evaluate_submission.py and run_baselines.py.

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

DOWNLOAD_FILES = ["train_X.csv", "test_X.csv"]


def download_task(task_name: str, token: str = None, dry_run: bool = False) -> int:
    """Download all fold feature matrices for one task. Returns file count."""
    from huggingface_hub import hf_hub_download, list_repo_files

    task_dir = TASKS_DIR / task_name
    if not task_dir.exists():
        print(f"  [WARN] Local task directory not found: {task_dir}")
        print(f"         Clone the GitHub repo first: git clone https://github.com/jak4013/GeneLab_benchmark")
        return 0

    # List available fold directories from HF
    try:
        all_files = list(list_repo_files(HF_REPO_ID, repo_type="dataset", token=token))
    except Exception as e:
        print(f"  [ERROR] Failed to list HF repo files: {e}")
        return 0

    task_files = [f for f in all_files if f.startswith(f"{task_name}/fold_")]

    if not task_files:
        print(f"  [SKIP] No files found in HF repo for {task_name}")
        return 0

    n_downloaded = 0
    for hf_path in sorted(task_files):
        # hf_path: "A5_skin_lomo/fold_MHU-2_test/train_X.csv"
        parts = hf_path.split("/")
        if len(parts) != 3 or parts[2] not in DOWNLOAD_FILES:
            continue

        _, fold_name, fname = parts
        local_path = task_dir / fold_name / fname

        if local_path.exists():
            print(f"  [SKIP] Already exists: {local_path.relative_to(BASE_DIR)}")
            n_downloaded += 1
            continue

        tag = "[DRY-RUN] " if dry_run else ""
        print(f"  {tag}Download: {hf_path} → {local_path.relative_to(BASE_DIR)}")

        if not dry_run:
            local_path.parent.mkdir(parents=True, exist_ok=True)
            downloaded = hf_hub_download(
                repo_id=HF_REPO_ID,
                filename=hf_path,
                repo_type="dataset",
                token=token,
                local_dir=str(BASE_DIR),
            )
            # hf_hub_download saves to local_dir/filename; rename if needed
            expected = BASE_DIR / hf_path
            if expected != local_path and expected.exists():
                expected.rename(local_path)

        n_downloaded += 1

    return n_downloaded


def main():
    parser = argparse.ArgumentParser(
        description="Download GeneLab benchmark feature matrices from HuggingFace"
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
    for task in tasks:
        print(f"\nDownloading {task} ← {HF_REPO_ID}")
        print("-" * 50)
        n = download_task(task, token=token, dry_run=args.dry_run)
        total += n
        print(f"  → {n} files")

    print(f"\nTotal: {total} files {'(dry-run)' if args.dry_run else 'downloaded'}.")


if __name__ == "__main__":
    main()

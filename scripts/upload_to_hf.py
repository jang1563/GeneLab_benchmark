#!/usr/bin/env python3
"""Upload SpaceBio-Bench / GeneLab Benchmark public fold packages to HF.

The script is designed as a release workflow, not only a file copier:

- dry-run upload plans work without `huggingface_hub` installed;
- release manifests are included in the HF package by default;
- optional remote diffs show create/replace actions against the HF repo;
- optional smoke validation downloads small uploaded files and checks SHA-256.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any


BASE_DIR = Path(__file__).resolve().parent.parent
TASKS_DIR = BASE_DIR / "tasks"
HF_REPO_ID = "jang1563/genelab-benchmark"

# GO tasks only (A1/A3 NO-GO excluded from the primary public fold package).
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
RELEASE_METADATA_FILES = [
    (BASE_DIR / "release" / "release_manifest.json", "manifest.json"),
    (
        BASE_DIR / "release" / "release_manifest.schema.json",
        "metadata/release_manifest.schema.json",
    ),
]
SMOKE_SMALL_FILENAMES = {
    "README.md",
    "manifest.json",
    "metadata/release_manifest.schema.json",
    "task_info.json",
    "fold_info.json",
    "selected_genes.txt",
}


@dataclass(frozen=True)
class UploadEntry:
    local_path: Path
    repo_path: str
    kind: str

    @property
    def size_bytes(self) -> int:
        return self.local_path.stat().st_size

    @property
    def sha256(self) -> str:
        return sha256_file(self.local_path)

    def as_plan_row(self, remote_files: set[str] | None = None) -> dict[str, Any]:
        action = "upload"
        if remote_files is not None:
            action = "replace" if self.repo_path in remote_files else "create"
        return {
            "action": action,
            "kind": self.kind,
            "repo_path": self.repo_path,
            "local_path": self.local_path.relative_to(BASE_DIR).as_posix(),
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
        }


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat().replace("+00:00", "Z")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def size_label(size_bytes: int) -> str:
    size_mb = size_bytes / 1e6
    if size_mb >= 1:
        return f"{size_mb:.1f} MB"
    return f"{size_bytes / 1024:.1f} KB"


def resolve_tasks(values: list[str]) -> list[str]:
    if "all" in values:
        return GO_TASKS.copy()

    tasks: list[str] = []
    seen: set[str] = set()
    for value in values:
        if value in GO_TASKS:
            matched = [value]
        else:
            matched = [name for name in GO_TASKS if name.startswith(value)]
        if not matched:
            print(f"[WARN] Unknown task: {value}. Available: {GO_TASKS}")
            continue
        for task in matched:
            if task not in seen:
                tasks.append(task)
                seen.add(task)
    return tasks


def collect_task_entries(task_name: str) -> tuple[list[UploadEntry], list[str]]:
    task_dir = TASKS_DIR / task_name
    warnings: list[str] = []
    entries: list[UploadEntry] = []
    if not task_dir.exists():
        warnings.append(f"{task_name}: directory not found")
        return entries, warnings

    for filename in TASK_ROOT_FILES:
        path = task_dir / filename
        if path.exists():
            entries.append(UploadEntry(path, f"{task_name}/{filename}", "task_file"))
        else:
            warnings.append(f"{path.relative_to(BASE_DIR)}: not found")

    for fold_dir in sorted(task_dir.glob("fold_*")):
        for filename in FOLD_PACKAGE_FILES:
            path = fold_dir / filename
            if path.exists():
                entries.append(
                    UploadEntry(
                        path,
                        f"{task_name}/{fold_dir.name}/{filename}",
                        "task_file",
                    )
                )
            else:
                warnings.append(f"{path.relative_to(BASE_DIR)}: not found")

    return entries, warnings


def collect_card_entries() -> tuple[list[UploadEntry], list[str]]:
    entries: list[UploadEntry] = []
    warnings: list[str] = []
    if CARD_SRC.exists():
        entries.append(UploadEntry(CARD_SRC, "README.md", "dataset_card"))
    else:
        warnings.append(f"{CARD_SRC.relative_to(BASE_DIR)}: not found")

    for asset in CARD_ASSETS:
        if asset.exists():
            entries.append(
                UploadEntry(asset, f"assets/{asset.name}", "dataset_card_asset")
            )
        else:
            warnings.append(f"{asset.relative_to(BASE_DIR)}: not found")
    return entries, warnings


def collect_release_metadata_entries() -> tuple[list[UploadEntry], list[str]]:
    entries: list[UploadEntry] = []
    warnings: list[str] = []
    for local_path, repo_path in RELEASE_METADATA_FILES:
        if local_path.exists():
            entries.append(UploadEntry(local_path, repo_path, "release_metadata"))
        else:
            warnings.append(f"{local_path.relative_to(BASE_DIR)}: not found")
    return entries, warnings


def collect_upload_entries(
    tasks: list[str],
    include_card: bool,
    include_release_metadata: bool,
) -> tuple[list[UploadEntry], list[str]]:
    entries: list[UploadEntry] = []
    warnings: list[str] = []

    for task in tasks:
        task_entries, task_warnings = collect_task_entries(task)
        entries.extend(task_entries)
        warnings.extend(task_warnings)

    if include_card:
        card_entries, card_warnings = collect_card_entries()
        entries.extend(card_entries)
        warnings.extend(card_warnings)

    if include_release_metadata:
        metadata_entries, metadata_warnings = collect_release_metadata_entries()
        entries.extend(metadata_entries)
        warnings.extend(metadata_warnings)

    return entries, warnings


def import_hf_api():
    try:
        from huggingface_hub import HfApi
    except ImportError:
        print("[ERROR] huggingface_hub not installed.")
        print("  Install: pip install huggingface_hub")
        sys.exit(1)
    return HfApi


def get_token(required: bool) -> str | None:
    token = os.environ.get("HF_TOKEN")
    if not token:
        try:
            from huggingface_hub import get_token as hf_get_token

            token = hf_get_token()
        except Exception:
            token = None
    if required and not token:
        print("[ERROR] No Hugging Face token found.")
        print("  Option 1: export HF_TOKEN='hf_...'")
        print("  Option 2: huggingface-cli login")
        sys.exit(1)
    return token


def list_remote_files(repo_id: str, token: str | None) -> set[str]:
    HfApi = import_hf_api()
    api = HfApi(token=token)
    try:
        return set(api.list_repo_files(repo_id=repo_id, repo_type="dataset"))
    except Exception as exc:
        print(f"[ERROR] Failed to list HF repo files for {repo_id}: {exc}")
        sys.exit(1)


def print_plan(
    entries: list[UploadEntry],
    repo_id: str,
    stale_remote_dirs: list[str],
    remote_files: set[str] | None,
    dry_run: bool,
) -> None:
    mode = "DRY-RUN" if dry_run else "UPLOAD"
    total_size = sum(entry.size_bytes for entry in entries)
    print(f"=== {mode} PLAN → {repo_id} ===")
    print(f"Files: {len(entries)} | Total local size: {size_label(total_size)}")
    if stale_remote_dirs:
        print(f"Stale remote dirs to prune: {len(stale_remote_dirs)}")
    print()

    by_kind: dict[str, int] = {}
    for entry in entries:
        by_kind[entry.kind] = by_kind.get(entry.kind, 0) + 1
    for kind, count in sorted(by_kind.items()):
        print(f"  {kind}: {count}")
    print()

    for entry in entries:
        row = entry.as_plan_row(remote_files=remote_files)
        print(
            f"  [{row['action']}] {entry.repo_path} "
            f"({entry.kind}, {size_label(entry.size_bytes)}, sha256={entry.sha256[:12]}...)"
        )
    for remote_dir in stale_remote_dirs:
        print(f"  [delete] {remote_dir}/")


def write_upload_plan(
    path: Path,
    repo_id: str,
    entries: list[UploadEntry],
    stale_remote_dirs: list[str],
    remote_files: set[str] | None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "generated_at": utc_now(),
        "repo_id": repo_id,
        "repo_type": "dataset",
        "files": [entry.as_plan_row(remote_files=remote_files) for entry in entries],
        "stale_remote_dirs": stale_remote_dirs,
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"\nWrote upload plan: {path}")


def upload_task_folder(api: Any, task_name: str, repo_id: str) -> int:
    task_entries, _ = collect_task_entries(task_name)
    if not task_entries:
        print(f"  [SKIP] {task_name}: no package files found")
        return 0

    allow_patterns = TASK_ROOT_FILES + [f"fold_*/{name}" for name in FOLD_PACKAGE_FILES]
    api.upload_folder(
        folder_path=str(TASKS_DIR / task_name),
        path_in_repo=task_name,
        repo_id=repo_id,
        repo_type="dataset",
        allow_patterns=allow_patterns,
        commit_message=f"Upload {task_name} self-contained fold package",
    )
    return len(task_entries)


def upload_file_entries(api: Any, entries: list[UploadEntry], repo_id: str) -> int:
    uploaded = 0
    for entry in entries:
        api.upload_file(
            path_or_fileobj=str(entry.local_path),
            path_in_repo=entry.repo_path,
            repo_id=repo_id,
            repo_type="dataset",
            commit_message=f"Upload {entry.repo_path}",
        )
        uploaded += 1
        print(f"  Uploaded: {entry.repo_path}")
    return uploaded


def prune_stale_remote_dirs(api: Any, repo_id: str, stale_remote_dirs: list[str]) -> None:
    for remote_dir in stale_remote_dirs:
        print(f"  Delete stale remote folder: {remote_dir}/")
        try:
            api.delete_folder(
                path_in_repo=remote_dir,
                repo_id=repo_id,
                repo_type="dataset",
                commit_message=f"Remove stale {remote_dir}",
            )
        except Exception as exc:
            print(f"  [WARN] Could not delete {remote_dir}/; it may already be absent: {exc}")


def choose_smoke_entries(entries: list[UploadEntry], limit: int) -> list[UploadEntry]:
    candidates = [
        entry
        for entry in entries
        if entry.size_bytes <= 1_000_000
        and (
            entry.repo_path in SMOKE_SMALL_FILENAMES
            or Path(entry.repo_path).name in SMOKE_SMALL_FILENAMES
            or entry.kind in {"dataset_card", "dataset_card_asset", "release_metadata"}
        )
    ]
    priority = {
        "release_metadata": 0,
        "dataset_card": 1,
        "dataset_card_asset": 2,
        "task_file": 3,
    }
    candidates.sort(key=lambda entry: (priority.get(entry.kind, 9), entry.repo_path))
    return candidates[:limit]


def validate_remote_files(
    repo_id: str,
    token: str | None,
    entries: list[UploadEntry],
    limit: int,
) -> None:
    try:
        from huggingface_hub import hf_hub_download
    except ImportError:
        print("[ERROR] huggingface_hub not installed.")
        print("  Install: pip install huggingface_hub")
        sys.exit(1)

    smoke_entries = choose_smoke_entries(entries, limit=limit)
    if not smoke_entries:
        print("[WARN] No small files selected for remote smoke validation.")
        return

    print(f"\nRemote smoke validation → {repo_id}")
    print("-" * 50)
    with tempfile.TemporaryDirectory(prefix="spacebiobench_hf_smoke_") as tmpdir:
        for entry in smoke_entries:
            downloaded = Path(
                hf_hub_download(
                    repo_id=repo_id,
                    filename=entry.repo_path,
                    repo_type="dataset",
                    token=token,
                    local_dir=tmpdir,
                    force_download=True,
                )
            )
            actual = sha256_file(downloaded)
            expected = entry.sha256
            if actual != expected:
                print(
                    f"[FAIL] {entry.repo_path}: sha256 mismatch "
                    f"expected {expected}, got {actual}",
                    file=sys.stderr,
                )
                sys.exit(1)
            print(f"  [OK] {entry.repo_path} ({size_label(entry.size_bytes)})")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Upload SpaceBio-Bench public fold packages to Hugging Face"
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
        help="Print the upload plan without changing the remote repo",
    )
    parser.add_argument(
        "--remote-diff",
        action="store_true",
        help="List remote files and label planned uploads as create/replace",
    )
    parser.add_argument(
        "--write-upload-plan",
        type=Path,
        help="Write a JSON upload plan with paths, sizes, and SHA-256 digests",
    )
    parser.add_argument(
        "--repo",
        default=HF_REPO_ID,
        help=f"Hugging Face dataset repo ID (default: {HF_REPO_ID})",
    )
    parser.add_argument(
        "--private",
        action="store_true",
        help="Create the HF dataset repo as private if it does not already exist",
    )
    parser.add_argument(
        "--card-only",
        action="store_true",
        help="Upload only the dataset card and card assets",
    )
    parser.add_argument(
        "--manifest-only",
        action="store_true",
        help="Upload only release manifest metadata files",
    )
    parser.add_argument(
        "--skip-card",
        action="store_true",
        help="Skip dataset card upload during a full package upload",
    )
    parser.add_argument(
        "--skip-manifest",
        action="store_true",
        help="Skip release manifest metadata during a full package upload",
    )
    parser.add_argument(
        "--skip-prune",
        action="store_true",
        help="Skip deletion of known stale remote folders such as A6 fold_TBD_test",
    )
    parser.add_argument(
        "--validate-remote",
        action="store_true",
        help=(
            "Download small planned files from HF and compare SHA-256; intended "
            "for post-upload or already-synced remote checks"
        ),
    )
    parser.add_argument(
        "--smoke-limit",
        type=int,
        default=12,
        help="Maximum number of files to download during --validate-remote",
    )
    args = parser.parse_args(argv)

    if args.card_only and args.manifest_only:
        print(
            "[ERROR] --card-only and --manifest-only are mutually exclusive.",
            file=sys.stderr,
        )
        return 1

    tasks = [] if args.card_only or args.manifest_only else resolve_tasks(args.task)
    if not tasks and not args.card_only and not args.manifest_only:
        print("No valid tasks specified. Exiting.")
        return 1

    include_card = args.card_only or (not args.manifest_only and not args.skip_card)
    include_release_metadata = args.manifest_only or (
        not args.card_only and not args.skip_manifest
    )
    stale_remote_dirs = (
        []
        if args.skip_prune or args.card_only or args.manifest_only
        else STALE_REMOTE_DIRS
    )

    entries, warnings = collect_upload_entries(
        tasks=tasks,
        include_card=include_card,
        include_release_metadata=include_release_metadata,
    )
    for warning in warnings:
        print(f"[WARN] {warning}")
    if not entries:
        print("No files selected. Exiting.")
        return 1

    token = None
    remote_files = None
    if args.remote_diff:
        token = get_token(required=False)
        remote_files = list_remote_files(args.repo, token=token)

    print_plan(
        entries=entries,
        repo_id=args.repo,
        stale_remote_dirs=stale_remote_dirs,
        remote_files=remote_files,
        dry_run=args.dry_run,
    )

    if args.write_upload_plan:
        write_upload_plan(
            path=args.write_upload_plan,
            repo_id=args.repo,
            entries=entries,
            stale_remote_dirs=stale_remote_dirs,
            remote_files=remote_files,
        )

    if not args.dry_run:
        HfApi = import_hf_api()
        token = token or get_token(required=True)
        api = HfApi(token=token)
        api.create_repo(
            repo_id=args.repo,
            repo_type="dataset",
            private=args.private,
            exist_ok=True,
        )

        if stale_remote_dirs:
            print(f"\nPruning stale remote paths → {args.repo}")
            print("-" * 50)
            prune_stale_remote_dirs(api, args.repo, stale_remote_dirs=stale_remote_dirs)

        total_task_files = 0
        for task in tasks:
            print(f"\nUploading {task} → {args.repo}")
            print("-" * 50)
            count = upload_task_folder(api, task, args.repo)
            total_task_files += count
            print(f"  → {count} task files")

        non_task_entries = [entry for entry in entries if entry.kind != "task_file"]
        if non_task_entries:
            print(f"\nUploading metadata/card files → {args.repo}")
            print("-" * 50)
            upload_file_entries(api, non_task_entries, args.repo)

        print(f"\nUploaded {len(entries)} planned files ({total_task_files} task files).")

    if args.validate_remote:
        token = token or get_token(required=False)
        validate_remote_files(
            repo_id=args.repo,
            token=token,
            entries=entries,
            limit=args.smoke_limit,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

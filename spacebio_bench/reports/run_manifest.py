"""Run-manifest and report writers for SpaceBio-Bench evaluations."""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence


SPACEBIO_BENCH_VERSION = "0.1.0-alpha"


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_run_manifest(
    *,
    task_manifest: Mapping[str, Any],
    task_manifest_path: str | Path | None,
    submission_path: str | Path,
    evaluation_result: Mapping[str, Any],
    command: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Build a provenance record for one evaluation run."""

    input_files: list[dict[str, str]] = [
        {
            "role": "submission",
            "path": Path(submission_path).as_posix(),
            "sha256": file_sha256(submission_path),
        }
    ]
    if task_manifest_path is not None:
        input_files.append(
            {
                "role": "task_manifest",
                "path": Path(task_manifest_path).as_posix(),
                "sha256": file_sha256(task_manifest_path),
            }
        )

    return {
        "schema_version": "0.1.0",
        "run_type": "spacebio_bench_evaluation",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "spacebio_bench_version": SPACEBIO_BENCH_VERSION,
        "task_id": task_manifest["task_id"],
        "task_family": task_manifest["task_family"],
        "status": evaluation_result.get("status"),
        "command": list(command if command is not None else sys.argv),
        "input_files": input_files,
        "metrics": evaluation_result.get("metrics", {}),
        "validation": evaluation_result.get("validation", {}),
    }


def write_evaluation_report(
    *,
    evaluation_result: Mapping[str, Any],
    task_manifest: Mapping[str, Any],
    task_manifest_path: str | Path | None,
    submission_path: str | Path,
    output_dir: str | Path,
    command: Sequence[str] | None = None,
) -> dict[str, Path]:
    """Write `metrics.json` and `run_manifest.json` for an evaluation."""

    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    metrics_path = outdir / "metrics.json"
    run_manifest_path = outdir / "run_manifest.json"

    metrics_path.write_text(json.dumps(evaluation_result, indent=2, sort_keys=True) + "\n")
    run_manifest = build_run_manifest(
        task_manifest=task_manifest,
        task_manifest_path=task_manifest_path,
        submission_path=submission_path,
        evaluation_result=evaluation_result,
        command=command,
    )
    run_manifest_path.write_text(json.dumps(run_manifest, indent=2, sort_keys=True) + "\n")

    return {
        "metrics": metrics_path,
        "run_manifest": run_manifest_path,
    }

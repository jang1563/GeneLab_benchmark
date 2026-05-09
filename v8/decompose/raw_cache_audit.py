"""Audit DECOMPOSE raw analog cache availability.

This does not fit models or regenerate DECOMPOSE outputs. It records whether
the raw count and sample-table files required by ``factorial_analog.py`` are
present in the current checkout/HPC staging bundle, plus hashes for any files
that are available.
"""
from __future__ import annotations

import hashlib
import json
import runpy
from pathlib import Path
from typing import Any

FACTORIAL_SCRIPT = Path(__file__).resolve().with_name("factorial_analog.py")
OUT_DIR = Path(__file__).resolve().parent / "evaluation"


def _sha256_file(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _rel(path: Path, repo_root: Path) -> str:
    try:
        return str(path.relative_to(repo_root))
    except ValueError:
        return str(path)


def _file_record(path: Path, repo_root: Path) -> dict[str, Any]:
    exists = path.exists()
    return {
        "path": _rel(path, repo_root),
        "exists": exists,
        "size_bytes": path.stat().st_size if exists else None,
        "sha256": _sha256_file(path) if exists else None,
    }


def main() -> None:
    ns = runpy.run_path(str(FACTORIAL_SCRIPT))
    datasets = ns["DATASETS"]
    repo_root = Path(ns["REPO_ROOT"])

    records = {}
    missing: list[str] = []
    present_files = 0
    expected_files = 0

    for analog, spec in datasets.items():
        analog_records = {}
        for role in ("counts", "samples"):
            expected_files += 1
            rec = _file_record(Path(spec[role]), repo_root)
            analog_records[role] = rec
            if rec["exists"]:
                present_files += 1
            else:
                missing.append(f"{analog}:{role}:{rec['path']}")
        records[analog] = {
            "design": spec.get("design"),
            "required_files": analog_records,
        }

    manifest = {
        "schema_version": "0.1.0",
        "audit_id": "decompose.raw_cache.current",
        "status": "complete" if not missing else "incomplete_raw_cache",
        "full_rerun_ready": not missing,
        "source_script": _rel(FACTORIAL_SCRIPT, repo_root),
        "expected_files": expected_files,
        "present_files": present_files,
        "missing_files": missing,
        "datasets": records,
        "limitations": [
            "This audit checks raw cache availability and file identity only; it does not regenerate factorial, Mars, or bootstrap outputs.",
            "A full DECOMPOSE beta freeze still requires a clean HPC rerun from fetched OSDR raw caches.",
        ],
    }

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUT_DIR / "raw_cache_audit.json"
    out_path.write_text(json.dumps(manifest, indent=2))
    print(json.dumps({
        "audit": str(out_path),
        "status": manifest["status"],
        "present_files": present_files,
        "expected_files": expected_files,
    }, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Audit an external raw-response archive for v8 INTERVENE API calls.

The raw JSON response dumps are intentionally kept out of Git. This script
checks that the archive contains one L1000CDS2 response per tissue and Enrichr
addList/enrich responses for each tissue/signature direction, then records
file-level checksums in a compact tracked audit manifest.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import sys
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO_ROOT / "v8" / "provenance" / "api_raw_archive_audit.json"
TISSUES = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]
DIRECTIONS = ["up", "dn"]
ENRICHR_STAGES = ["addList", "enrich"]
ENRICHR_LIBRARY = "LINCS_L1000_CRISPR_KO_Consensus_Sigs"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_json(path: Path) -> Any:
    return json.loads(path.read_text())


def expected_files() -> list[dict[str, str]]:
    files: list[dict[str, str]] = []
    for tissue in TISSUES:
        files.append({"api": "l1000cds2", "tissue": tissue, "path": f"l1000cds2_{tissue}_response.json"})
    for tissue in TISSUES:
        for direction in DIRECTIONS:
            for stage in ENRICHR_STAGES:
                files.append(
                    {
                        "api": "enrichr",
                        "tissue": tissue,
                        "direction": direction,
                        "stage": stage,
                        "path": f"enrichr_{tissue}_{direction}_{stage}.json",
                    }
                )
    return files


def summarize_record(spec: dict[str, str], path: Path) -> dict[str, Any]:
    item: dict[str, Any] = {
        **spec,
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
        "present": True,
    }
    data = load_json(path)
    item["has_error"] = "error" in data
    item["response_status_code"] = data.get("response_status_code")
    item["endpoint"] = data.get("endpoint")

    if spec["api"] == "l1000cds2":
        cfg = data.get("request_payload", {}).get("config", {})
        response = data.get("response_json", {})
        item["request_db_version"] = cfg.get("db-version")
        item["n_top_meta"] = len(response.get("topMeta", [])) if isinstance(response, dict) else None
    elif spec["api"] == "enrichr":
        response = data.get("response_json", {})
        if spec.get("stage") == "addList":
            item["has_user_list_id"] = isinstance(response, dict) and "userListId" in response
        else:
            item["library"] = ENRICHR_LIBRARY
            item["n_terms"] = len(response.get(ENRICHR_LIBRARY, [])) if isinstance(response, dict) else None
    return item


def build_audit(raw_dir: Path, archive_id: str, generated_at: str, root_hint: str | None) -> dict[str, Any]:
    raw_root = raw_dir.expanduser().resolve()
    records: list[dict[str, Any]] = []
    missing: list[str] = []
    errored: list[str] = []
    bundle = hashlib.sha256()

    for spec in expected_files():
        rel = spec["path"]
        path = raw_root / rel
        if not path.exists():
            missing.append(rel)
            records.append({**spec, "present": False})
            continue
        item = summarize_record(spec, path)
        records.append(item)
        if item.get("has_error"):
            errored.append(rel)
        bundle.update(rel.encode("utf-8"))
        bundle.update(b"\0")
        bundle.update(item["sha256"].encode("ascii"))
        bundle.update(b"\n")

    status = "release_candidate" if not missing and not errored else "incomplete"
    db_versions = sorted(
        {
            str(record.get("request_db_version"))
            for record in records
            if record.get("api") == "l1000cds2" and record.get("request_db_version")
        }
    )
    return {
        "schema_version": "0.1.0",
        "audit_id": "v8.intervene.api_raw_archive",
        "archive_id": archive_id,
        "status": status,
        "generated_at": generated_at,
        "raw_archive_root_hint": root_hint or str(raw_root),
        "bundle": {
            "kind": "raw_api_response_set",
            "sha256": bundle.hexdigest() if records else None,
            "n_expected_files": len(expected_files()),
            "n_present_files": sum(1 for record in records if record.get("present")),
            "digest_algorithm": "sha256(relpath + NUL + file_sha256 + LF, expected file order)",
        },
        "l1000cds2": {
            "n_tissues": len(TISSUES),
            "request_db_versions": db_versions,
            "note": "If this remains 'latest', preserve raw response dumps externally for beta freeze.",
        },
        "enrichr_crispr": {
            "n_tissues": len(TISSUES),
            "directions": DIRECTIONS,
            "library": ENRICHR_LIBRARY,
        },
        "files": records,
        "missing_files": missing,
        "errored_files": errored,
        "limitations": [
            "This audit records checksums for external raw API response dumps; the dumps themselves are not tracked in Git.",
            "The archive should be copied to the final Hugging Face, Zenodo, or object-storage target before marking v8 beta frozen.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", required=True, help="Directory containing raw API response JSON files.")
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Audit JSON output path.")
    parser.add_argument("--archive-id", default="v8-intervene-api-raw-archive-rc1")
    parser.add_argument("--root-hint", default=None, help="Stable external path or object-storage URI to record.")
    parser.add_argument(
        "--generated-at",
        default=dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
    )
    args = parser.parse_args()

    audit = build_audit(Path(args.raw_dir), args.archive_id, args.generated_at, args.root_hint)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(audit, indent=2) + "\n")
    print(f"Wrote {out}")
    print(f"API raw archive status: {audit['status']} ({audit['bundle']['n_present_files']}/{audit['bundle']['n_expected_files']} files)")
    if audit["status"] != "release_candidate":
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

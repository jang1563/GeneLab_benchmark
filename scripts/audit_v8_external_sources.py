#!/usr/bin/env python3
"""Audit v8 external processed inputs without copying raw data into Git.

The beta release needs a stable identity for the SpaceOmicsBench/SOMA files used
by promoted v8 outputs. This script records file-level SHA-256 digests and a
deterministic bundle digest over the required processed files. It intentionally
does not archive the files themselves.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SPACEOMICS_ROOT = REPO_ROOT.parent / "SpaceOmicsBench"
DEFAULT_OUT = REPO_ROOT / "v8" / "provenance" / "external_source_audit.json"

SPACEOMICS_PROCESSED_FILES = [
    {
        "path": "v2_public/data/processed/gt_conserved_pathways_GeneLab.csv",
        "required_for": ["bridge.species_transfer_nes"],
    },
    {
        "path": "v2_public/data/processed/gt_conserved_pathways_i4_pbmc.csv",
        "required_for": ["bridge.species_transfer_nes", "bridge.tissue_nes_spearman"],
    },
    {
        "path": "v2_public/data/processed/gt_conserved_pathways_NASA_Twins.csv",
        "required_for": ["bridge.species_transfer_nes", "bridge.tissue_nes_spearman"],
    },
    {
        "path": "v2_public/data/processed/cross_mission_pathway_features.csv",
        "required_for": ["bridge.supervised_conservation", "bridge.leakage_audit"],
    },
    {
        "path": "v2_public/data/processed/conserved_pbmc_to_skin.csv",
        "required_for": ["future.bridge.pbmc_skin_gene_mapping"],
    },
    {
        "path": "v2_public/data/processed/proteomics_plasma_matrix.csv",
        "required_for": ["future.multiomics.propagation"],
    },
    {
        "path": "v2_public/data/processed/metabolomics_anppos_matrix.csv",
        "required_for": ["future.multiomics.propagation"],
    },
    {
        "path": "v2_public/data/processed/clinical_merged_serum.csv",
        "required_for": ["future.multiomics.propagation"],
    },
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def git_metadata(root: Path) -> dict[str, Any]:
    proc = subprocess.run(
        ["git", "rev-parse", "--is-inside-work-tree"],
        cwd=root,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
        check=False,
    )
    if proc.returncode != 0 or proc.stdout.strip() != "true":
        return {
            "available": False,
            "notes": "No .git metadata found; source identity is the processed-file bundle checksum.",
        }

    def git(*args: str) -> str | None:
        out = subprocess.run(
            ["git", *args],
            cwd=root,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            check=False,
        )
        return out.stdout.strip() if out.returncode == 0 else None

    status = git("status", "--short", "--untracked-files=no")
    return {
        "available": True,
        "commit": git("rev-parse", "HEAD"),
        "tags_at_head": [tag for tag in (git("tag", "--points-at", "HEAD") or "").splitlines() if tag],
        "tracked_changes": bool(status),
    }


def load_json(path: Path) -> Any | None:
    try:
        return json.loads(path.read_text())
    except FileNotFoundError:
        return None
    except json.JSONDecodeError:
        return None


def evidence_checks(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    checks: list[dict[str, Any]] = []
    by_path = {record["path"]: record for record in records}
    leakage = load_json(REPO_ROOT / "v8" / "bridge" / "evaluation" / "bridge_leakage_audit.json")
    if leakage:
        expected = (
            leakage.get("inputs", {})
            .get("i2_features", {})
            .get("sha256")
        )
        rel = "v2_public/data/processed/cross_mission_pathway_features.csv"
        actual = by_path.get(rel, {}).get("sha256")
        checks.append(
            {
                "name": "bridge_leakage_audit_i2_features_sha256",
                "path": rel,
                "expected_sha256": expected,
                "actual_sha256": actual,
                "status": "pass" if expected and actual == expected else "fail",
            }
        )
    return checks


def build_audit(spaceomics_root: Path, generated_at: str) -> dict[str, Any]:
    root = spaceomics_root.expanduser().resolve()
    if not root.exists():
        raise FileNotFoundError(f"SpaceOmicsBench root does not exist: {root}")

    records: list[dict[str, Any]] = []
    bundle = hashlib.sha256()
    for spec in SPACEOMICS_PROCESSED_FILES:
        rel = spec["path"]
        path = root / rel
        if not path.exists():
            raise FileNotFoundError(f"Required SpaceOmicsBench input is missing: {path}")
        digest = sha256_file(path)
        record = {
            "path": rel,
            "sha256": digest,
            "bytes": path.stat().st_size,
            "required_for": spec["required_for"],
        }
        records.append(record)
        bundle.update(rel.encode("utf-8"))
        bundle.update(b"\0")
        bundle.update(digest.encode("ascii"))
        bundle.update(b"\n")

    checks = evidence_checks(records)
    return {
        "schema_version": "0.1.0",
        "audit_id": "v8.external_sources.spaceomicsbench_processed",
        "status": "release_candidate",
        "generated_at": generated_at,
        "source_id": "spaceomicsbench.v2_public",
        "root_hint": "$SPACEOMICS_ROOT or ../SpaceOmicsBench",
        "git": git_metadata(root),
        "bundle": {
            "kind": "processed_file_set",
            "sha256": bundle.hexdigest(),
            "n_files": len(records),
            "digest_algorithm": "sha256(relpath + NUL + file_sha256 + LF, declared file order)",
        },
        "files": records,
        "evidence_checks": checks,
        "limitations": [
            "This audit identifies the processed SpaceOmicsBench inputs used by v8; it does not archive raw SOMA/OSDR/API payloads.",
            "If a future SpaceOmicsBench checkout has Git metadata, record its commit or release tag in addition to this bundle checksum.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--spaceomics-root",
        default=os.environ.get("SPACEOMICS_ROOT", str(DEFAULT_SPACEOMICS_ROOT)),
        help="SpaceOmicsBench checkout or processed-source bundle root.",
    )
    parser.add_argument("--out", default=str(DEFAULT_OUT), help="Audit JSON output path.")
    parser.add_argument(
        "--generated-at",
        default=dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z"),
        help="ISO timestamp to record in the audit manifest.",
    )
    args = parser.parse_args()

    audit = build_audit(Path(args.spaceomics_root), args.generated_at)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(audit, indent=2) + "\n")
    print(f"Wrote {out}")
    print(f"SpaceOmicsBench processed bundle sha256: {audit['bundle']['sha256']}")
    failed = [check for check in audit["evidence_checks"] if check["status"] != "pass"]
    if failed:
        print("Evidence check failures:", file=sys.stderr)
        for check in failed:
            print(f"  - {check['name']}: expected {check['expected_sha256']} got {check['actual_sha256']}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

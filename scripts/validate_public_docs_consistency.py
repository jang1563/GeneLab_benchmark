#!/usr/bin/env python3
"""Validate public README, Hugging Face card, and citation consistency."""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
README = REPO_ROOT / "README.md"
HF_CARD = REPO_ROOT / "docs" / "hf_dataset_card.md"
CITATION = REPO_ROOT / "CITATION.cff"
ZENODO = REPO_ROOT / ".zenodo.json"
RELEASE_MANIFEST = REPO_ROOT / "release" / "release_manifest.json"


def load_json(path: Path) -> Any:
    return json.loads(path.read_text())


def parse_cff_scalar(text: str, key: str) -> str | None:
    pattern = re.compile(rf"^{re.escape(key)}:\s*\"?([^\"\n]+)\"?\s*$", re.MULTILINE)
    match = pattern.search(text)
    return match.group(1).strip() if match else None


def parse_front_matter(text: str) -> dict[str, str]:
    if not text.startswith("---\n"):
        return {}
    end = text.find("\n---\n", 4)
    if end == -1:
        return {}
    values: dict[str, str] = {}
    for line in text[4:end].splitlines():
        if not line or line.startswith(" ") or line.lstrip().startswith("-"):
            continue
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        values[key.strip()] = value.strip().strip('"')
    return values


def require_contains(errors: list[str], label: str, text: str, expected: str) -> None:
    if expected not in text:
        errors.append(f"{label}: missing expected text {expected!r}")


def require_absent(errors: list[str], label: str, text: str, forbidden: str) -> None:
    if forbidden in text:
        errors.append(f"{label}: stale or forbidden text present {forbidden!r}")


def validate_public_docs() -> list[str]:
    errors: list[str] = []
    manifest = load_json(RELEASE_MANIFEST)
    readme = README.read_text()
    hf_card = HF_CARD.read_text()
    citation = CITATION.read_text()
    zenodo = load_json(ZENODO)
    hf_front_matter = parse_front_matter(hf_card)

    lanes = {lane["lane_id"]: lane for lane in manifest["release_lanes"]}
    v7 = lanes.get("v7.1", {})
    v9 = lanes.get("v9_public_bulk", {})

    citation_version = v7.get("citation", {}).get("citation_version")
    if citation_version:
        actual_version = parse_cff_scalar(citation, "version")
        if actual_version != citation_version:
            errors.append(
                f"CITATION.cff: version {actual_version!r} does not match "
                f"release manifest citation_version {citation_version!r}"
            )

    require_contains(errors, "README.md", readme, "# SpaceBio-Bench")
    require_contains(errors, "README.md", readme, "Former public name: **GeneLab Benchmark**")
    require_contains(errors, "README.md", readme, "v7.1 GeneLab Benchmark")
    require_contains(errors, "README.md", readme, "release/release_manifest.json")
    require_contains(errors, "README.md", readme, "https://huggingface.co/datasets/jang1563/genelab-benchmark")
    require_contains(errors, "README.md", readme, "docs/assets/hf_benchmark_summary.png")
    require_absent(errors, "README.md", readme, "Version: v7.0 (2026-04-12)")
    require_absent(errors, "README.md", readme, "Status: **v1–v7 Complete**")

    expected_pretty_name = "SpaceBio-Bench / GeneLab Benchmark v7.1.2 Public Fold Package"
    if hf_front_matter.get("pretty_name") != expected_pretty_name:
        errors.append(
            "docs/hf_dataset_card.md: pretty_name does not match expected "
            f"{expected_pretty_name!r}"
        )
    if hf_front_matter.get("license") != "cc-by-4.0":
        errors.append("docs/hf_dataset_card.md: expected license cc-by-4.0")
    if hf_front_matter.get("viewer") != "false":
        errors.append("docs/hf_dataset_card.md: expected viewer: false")

    require_contains(errors, "README.md", readme, "v7.1.2 public-card/metadata patch")
    require_contains(errors, "docs/hf_dataset_card.md", hf_card, "Public status: **v7.1.2 public-card/metadata patch")
    require_contains(errors, "docs/hf_dataset_card.md", hf_card, "Dataset freeze: **2026-03-01**")
    require_contains(errors, "docs/hf_dataset_card.md", hf_card, "repo_id = \"jang1563/genelab-benchmark\"")

    zenodo_text = json.dumps(zenodo, sort_keys=True)
    require_contains(errors, ".zenodo.json", zenodo.get("title", ""), "SpaceBio-Bench")
    require_contains(errors, ".zenodo.json", zenodo.get("version", ""), "7.1.2")
    require_contains(errors, ".zenodo.json", zenodo.get("description", ""), "public-card and metadata patch")
    for forbidden in ("evidence-visibility", "release-candidate surface", "hiring-manager"):
        require_absent(errors, ".zenodo.json", zenodo_text, forbidden)
    if v7.get("public_label"):
        require_contains(errors, "README.md", readme.lower(), v7["public_label"].lower())
        require_contains(errors, "docs/hf_dataset_card.md", hf_card.lower(), v7["public_label"].lower())
    if v9.get("public_label"):
        require_contains(errors, "README.md", readme.lower(), v9["public_label"].lower())

    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate public README, HF card, and citation consistency"
    )
    parser.parse_args(argv)

    errors = validate_public_docs()
    if errors:
        print(f"[FAIL] public docs consistency: {len(errors)} issue(s)", file=sys.stderr)
        for error in errors:
            print(f"  - {error}", file=sys.stderr)
        return 1

    print("[OK] public docs consistency")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

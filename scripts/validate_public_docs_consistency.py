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
V9_HF_CARD = REPO_ROOT / "docs" / "v9_hf_dataset_card.md"
DOC_MAP = REPO_ROOT / "docs" / "SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md"
PORTFOLIO_BRIEF = REPO_ROOT / "docs" / "SPACEBIOBENCH_PORTFOLIO_BRIEF.md"
SYSTEM_CARD = REPO_ROOT / "docs" / "SPACEBIOBENCH_SYSTEM_CARD.md"
EVALUATION_CARD = REPO_ROOT / "docs" / "SPACEBIOBENCH_EVALUATION_CARD.md"
RELEASE_STATUS_CARD = REPO_ROOT / "docs" / "SPACEBIOBENCH_RELEASE_READINESS_CARD.md"
STATEMENT_GUIDE = REPO_ROOT / "docs" / "SPACEBIOBENCH_CLAIM_REGISTER.md"
ARCHIVE_STATUS = REPO_ROOT / "docs" / "RELEASE_ARCHIVE_CARD.md"
ARCHIVE_MANIFEST = REPO_ROOT / "docs" / "RELEASE_ARCHIVE_MANIFEST.md"
ARCHIVE_CHECKLIST = REPO_ROOT / "docs" / "RELEASE_ARCHIVE_CHECKLIST.md"
V9_README = REPO_ROOT / "v9" / "README.md"
V9_REPORTS_README = REPO_ROOT / "v9" / "reports" / "README.md"
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
    v9_hf_card = V9_HF_CARD.read_text()
    doc_map = DOC_MAP.read_text()
    portfolio_brief = PORTFOLIO_BRIEF.read_text()
    system_card = SYSTEM_CARD.read_text()
    evaluation_card = EVALUATION_CARD.read_text()
    release_status_card = RELEASE_STATUS_CARD.read_text()
    statement_guide = STATEMENT_GUIDE.read_text()
    archive_status = ARCHIVE_STATUS.read_text()
    archive_manifest = ARCHIVE_MANIFEST.read_text()
    archive_checklist = ARCHIVE_CHECKLIST.read_text()
    v9_readme = V9_README.read_text()
    v9_reports_readme = V9_REPORTS_README.read_text()
    citation = CITATION.read_text()
    zenodo = load_json(ZENODO)
    hf_front_matter = parse_front_matter(hf_card)
    v9_front_matter = parse_front_matter(v9_hf_card)

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
    citation_title = parse_cff_scalar(citation, "title")
    if citation_title and "SpaceBio-Bench / GeneLab Benchmark" not in citation_title:
        errors.append("CITATION.cff: title should include SpaceBio-Bench / GeneLab Benchmark")

    require_contains(errors, "README.md", readme, "# SpaceBio-Bench")
    require_contains(errors, "README.md", readme, "Former public name: **GeneLab Benchmark**")
    require_contains(errors, "README.md", readme, "v7.1 GeneLab Benchmark")
    require_contains(errors, "README.md", readme, "release/release_manifest.json")
    require_contains(errors, "README.md", readme, "https://huggingface.co/datasets/jang1563/genelab-benchmark")
    require_contains(errors, "README.md", readme, "docs/assets/hf_benchmark_summary.png")
    require_contains(errors, "README.md", readme, "CONTRIBUTING.md")
    require_contains(errors, "README.md", readme, "docs/submission_format.md")
    require_contains(errors, "README.md", readme, "SpaceBio-Bench / GeneLab Benchmark: Mission-Held-Out")
    require_contains(errors, "README.md", readme, "public documentation map")
    require_contains(errors, "README.md", readme, "Release status")
    require_contains(errors, "README.md", readme, "Public statement guide")
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
    require_contains(errors, "docs/hf_dataset_card.md", hf_card, "SpaceBio-Bench / GeneLab Benchmark: Mission-Held-Out")
    require_absent(
        errors,
        "docs/hf_dataset_card.md",
        hf_card,
        "Held-out mission, train missions, sample counts, and provenance",
    )

    expected_v9_pretty_name = "SpaceBio-Bench v9 Public Bulk Metadata Catalog"
    if v9_front_matter.get("pretty_name") != expected_v9_pretty_name:
        errors.append(
            "docs/v9_hf_dataset_card.md: pretty_name does not match expected "
            f"{expected_v9_pretty_name!r}"
        )
    require_contains(
        errors,
        "docs/v9_hf_dataset_card.md",
        v9_hf_card,
        "# SpaceBio-Bench v9 Public Bulk Metadata Catalog",
    )
    require_contains(
        errors,
        "docs/v9_hf_dataset_card.md",
        v9_hf_card,
        "Current catalog scope:",
    )
    require_contains(
        errors,
        "v9/README.md",
        v9_readme,
        "# SpaceBio-Bench v9 Public Bulk Metadata Catalog",
    )
    require_contains(
        errors,
        "v9/reports/README.md",
        v9_reports_readme,
        "# v9 Public Bulk Reports",
    )
    v9_public_text = "\n".join([v9_hf_card, v9_readme, v9_reports_readme])
    for forbidden in (
        "Public Bulk Metadata Alpha",
        "metadata-alpha",
        "alpha snapshot",
        "alpha-boundary",
        "claim boundary",
        "snapshot_decision",
        "allowed language",
        "blocked language",
        "Maintenance Notes",
        "Provenance And Integrity",
        "not frozen",
        "should not be used",
        "release-readiness blockers",
    ):
        require_absent(errors, "v9 public docs", v9_public_text, forbidden)

    public_card_text = "\n".join(
        [
            doc_map,
            portfolio_brief,
            system_card,
            evaluation_card,
            release_status_card,
            statement_guide,
        ]
    )
    for label, text, expected in (
        (
            "docs/SPACEBIOBENCH_TRANSPARENCY_CARD_PACK.md",
            doc_map,
            "# SpaceBio-Bench Public Documentation Map",
        ),
        (
            "docs/SPACEBIOBENCH_PORTFOLIO_BRIEF.md",
            portfolio_brief,
            "# SpaceBio-Bench Portfolio Brief",
        ),
        ("docs/SPACEBIOBENCH_SYSTEM_CARD.md", system_card, "# SpaceBio-Bench System Card"),
        (
            "docs/SPACEBIOBENCH_EVALUATION_CARD.md",
            evaluation_card,
            "# SpaceBio-Bench Evaluation Card",
        ),
        (
            "docs/SPACEBIOBENCH_RELEASE_READINESS_CARD.md",
            release_status_card,
            "# SpaceBio-Bench Release Status Card",
        ),
        (
            "docs/SPACEBIOBENCH_CLAIM_REGISTER.md",
            statement_guide,
            "# SpaceBio-Bench Public Statement Guide",
        ),
    ):
        require_contains(errors, label, text, expected)
    for forbidden in (
        "metadata-alpha",
        "metadata alpha",
        "Public Bulk Metadata Alpha",
        "claim boundary",
        "blocked wording",
        "future-only",
        "provenance boundary",
        "Provenance And Integrity",
        "release-readiness blockers",
        "should not be used",
        "not a frozen",
    ):
        require_absent(errors, "linked public cards", public_card_text, forbidden)

    archive_text = "\n".join([archive_status, archive_manifest, archive_checklist])
    for label, text, expected in (
        (
            "docs/RELEASE_ARCHIVE_CARD.md",
            archive_status,
            "# SpaceBio-Bench Release Archive Status",
        ),
        (
            "docs/RELEASE_ARCHIVE_MANIFEST.md",
            archive_manifest,
            "# SpaceBio-Bench Release Archive Manifest",
        ),
        (
            "docs/RELEASE_ARCHIVE_CHECKLIST.md",
            archive_checklist,
            "# SpaceBio-Bench Release Archive Checklist",
        ),
    ):
        require_contains(errors, label, text, expected)
    for forbidden in (
        "metadata-alpha",
        "metadata alpha",
        "claim-boundary",
        "claim boundary",
        "release_candidate",
        "release-candidate",
        "public_review_ready",
        "blocked",
        "not a frozen",
        "should not",
        "/Users/",
        "~/.claude",
    ):
        require_absent(errors, "release archive docs", archive_text, forbidden)

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

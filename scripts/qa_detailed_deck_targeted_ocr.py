#!/usr/bin/env python3
"""Targeted OCR QA for selected detailed-deck slides.

Use this after editing a slide batch so visible-copy regressions can be checked
without rebuilding the full 60-slide walkthrough contact sheets.
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
import sys
from pathlib import Path

from PIL import Image


ROOT = Path(__file__).resolve().parent.parent
WORKSPACE = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
)
PLAN = WORKSPACE / "spacebiobench_detailed_deck_assembly_plan_v2.tsv"
OUT_DIR = WORKSPACE / "assets" / "full_walkthrough_qa" / "targeted_ocr"
OUT_DIR.mkdir(parents=True, exist_ok=True)

EXPECTED_SIZE = (3840, 2160)
FLAG_PATTERNS = [
    ("source/provenance", re.compile(r"\b(source|sources|provenance|source files|source notes)\b", re.I)),
    ("defensive-claim", re.compile(r"\b(claim|claimed|claiming|claim strength|claim status|claim scope)\b", re.I)),
    ("boundary/scope", re.compile(r"\b(boundary|boundaries|scope|caveat|guardrail|guardrails)\b", re.I)),
    ("reader-instruction", re.compile(r"\b(reader rule|how to read|read the stack|teaching rule)\b", re.I)),
    ("internal-gate", re.compile(r"\b(evidence gates?|diagnostic-only|payload|blocked|blocker|decision rule|internal)\b", re.I)),
    ("proof-language", re.compile(r"\b(proof object|proof layer|proof slides?)\b", re.I)),
]


def parse_slides(args: list[str]) -> list[int]:
    slides: set[int] = set()
    for arg in args:
        for part in arg.split(","):
            part = part.strip()
            if not part:
                continue
            if "-" in part:
                start, end = part.split("-", 1)
                slides.update(range(int(start), int(end) + 1))
            else:
                slides.add(int(part))
    if not slides:
        raise SystemExit("usage: qa_detailed_deck_targeted_ocr.py SLIDE[,SLIDE] or START-END")
    return sorted(slides)


def load_rows(slides: list[int]) -> list[dict[str, str]]:
    wanted = {str(slide) for slide in slides}
    rows = []
    with PLAN.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["slide"] not in wanted:
                continue
            asset = ROOT / row["asset_or_action"]
            if not asset.exists():
                raise FileNotFoundError(asset)
            with Image.open(asset) as image:
                if image.size != EXPECTED_SIZE or image.mode != "RGB":
                    raise ValueError(f"slide {row['slide']} has {image.size} {image.mode}: {asset}")
            row["asset_abs"] = str(asset)
            rows.append(row)
    missing = sorted(wanted - {row["slide"] for row in rows}, key=int)
    if missing:
        raise ValueError(f"slides not found in assembly plan: {', '.join(missing)}")
    return sorted(rows, key=lambda row: int(row["slide"]))


def ocr_slide(row: dict[str, str], batch_dir: Path) -> str:
    src = Path(row["asset_abs"])
    downscaled = batch_dir / f"slide_{int(row['slide']):02d}_ocr_input.png"
    with Image.open(src) as image:
        image = image.convert("RGB")
        resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
        image.thumbnail((2200, 1240), resample)
        image.save(downscaled, quality=95)
    try:
        result = subprocess.run(
            ["tesseract", str(downscaled), "stdout", "--psm", "6"],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=25,
        )
    except subprocess.TimeoutExpired:
        return "[OCR timed out]"
    if result.returncode != 0:
        return f"[OCR failed: {result.stderr.strip()}]"
    return result.stdout


def main() -> None:
    slides = parse_slides(sys.argv[1:])
    batch_name = "_".join(f"{slide:02d}" for slide in slides)
    batch_dir = OUT_DIR / batch_name
    batch_dir.mkdir(parents=True, exist_ok=True)
    rows = load_rows(slides)

    text_blocks: list[str] = []
    flags: list[dict[str, str]] = []
    for row in rows:
        text_out = ocr_slide(row, batch_dir)
        text_blocks.append(f"--- slide {int(row['slide']):02d}: {row['title']} ---\n{text_out.strip()}\n")
        for label, pattern in FLAG_PATTERNS:
            matches = sorted(set(match.group(0) for match in pattern.finditer(text_out)))
            if matches:
                flags.append(
                    {
                        "slide": f"{int(row['slide']):02d}",
                        "title": row["title"],
                        "category": label,
                        "matches": ", ".join(matches),
                    }
                )

    text_path = batch_dir / "targeted_ocr_text.txt"
    flags_path = batch_dir / "targeted_ocr_flags.tsv"
    text_path.write_text("\n".join(text_blocks), encoding="utf-8")
    with flags_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["slide", "title", "category", "matches"], delimiter="\t")
        writer.writeheader()
        writer.writerows(flags)

    print(
        json.dumps(
            {
                "slides": slides,
                "flag_count": len(flags),
                "flags": str(flags_path.relative_to(ROOT)),
                "text": str(text_path.relative_to(ROOT)),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Full-walkthrough QA for the 60-slide detailed SpaceBio-Bench deck.

The script reads the assembly plan, validates slide assets, renders overview
and section contact sheets, and optionally OCRs each rendered slide to catch
defensive or internal-facing visible copy before PPTX assembly.
"""

from __future__ import annotations

import csv
import json
import re
import subprocess
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
WORKSPACE = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
)
PLAN = WORKSPACE / "spacebiobench_detailed_deck_assembly_plan_v2.tsv"
OUT_DIR = WORKSPACE / "assets" / "full_walkthrough_qa"
OUT_DIR.mkdir(parents=True, exist_ok=True)

OVERVIEW = OUT_DIR / "spacebiobench_detailed_deck_60_slide_overview_contact_sheet.png"
OVERVIEW_GRAY = OUT_DIR / "spacebiobench_detailed_deck_60_slide_overview_contact_sheet_grayscale.png"
MANIFEST = OUT_DIR / "spacebiobench_detailed_deck_full_walkthrough_manifest.json"
REPORT = OUT_DIR / "spacebiobench_detailed_deck_full_walkthrough_QA.md"
OCR_TEXT = OUT_DIR / "spacebiobench_detailed_deck_ocr_text.txt"
OCR_FLAGS = OUT_DIR / "spacebiobench_detailed_deck_ocr_flags.tsv"
OCR_DIR = OUT_DIR / "ocr_downscaled"
OCR_DIR.mkdir(parents=True, exist_ok=True)

EXPECTED_SIZE = (3840, 2160)

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "line": "#2A394D",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#6D7A8A",
    "open": "#70A7FF",
    "method": "#5FD3C4",
    "core": "#8BD17C",
    "models": "#B39DFF",
    "robust": "#F4C26B",
    "biology": "#F27D8E",
    "translation": "#73D6FF",
    "v8": "#D9B36A",
    "platform": "#A6D672",
    "extension": "#F59E6C",
    "close": "#E7EEF8",
}

ACT_COLORS = {
    "Open": COLORS["open"],
    "Method": COLORS["method"],
    "Core Result": COLORS["core"],
    "Models": COLORS["models"],
    "Robustness": COLORS["robust"],
    "Biology": COLORS["biology"],
    "Translation": COLORS["translation"],
    "v8": COLORS["v8"],
    "Platform": COLORS["platform"],
    "Extension": COLORS["extension"],
    "Close": COLORS["close"],
}

FLAG_PATTERNS = [
    ("source/provenance", re.compile(r"\b(source|sources|provenance|source files|source notes)\b", re.I)),
    ("defensive-claim", re.compile(r"\b(claim|claimed|claiming|claim strength|claim status|claim scope)\b", re.I)),
    ("boundary/scope", re.compile(r"\b(boundary|boundaries|scope|caveat|guardrail|guardrails)\b", re.I)),
    ("reader-instruction", re.compile(r"\b(reader rule|how to read|read the stack|teaching rule)\b", re.I)),
    ("internal-gate", re.compile(r"\b(evidence gates?|diagnostic-only|payload|blocked|blocker|decision rule|internal)\b", re.I)),
    ("proof-language", re.compile(r"\b(proof object|proof layer|proof slides?)\b", re.I)),
]


def font(size: int, *, bold: bool = False) -> ImageFont.ImageFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Helvetica.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            continue
    return ImageFont.load_default()


F_TITLE = font(54, bold=True)
F_SUB = font(24)
F_LABEL = font(22, bold=True)
F_SMALL = font(18)
F_TINY = font(15)
F_NUM = font(28, bold=True)


def wrap(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.ImageFont, width: int, max_lines: int) -> list[str]:
    words = text.split()
    lines: list[str] = []
    current: list[str] = []
    for word in words:
        trial = " ".join(current + [word])
        if draw.textlength(trial, font=fnt) <= width:
            current.append(word)
        else:
            if current:
                lines.append(" ".join(current))
            current = [word]
    if current:
        lines.append(" ".join(current))
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        if lines[-1] and not lines[-1].endswith("..."):
            lines[-1] = lines[-1].rstrip(".") + "..."
    return lines


def fit(path: Path, size: tuple[int, int]) -> Image.Image:
    image = Image.open(path).convert("RGB")
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    image.thumbnail(size, resample)
    canvas = Image.new("RGB", size, COLORS["bg"])
    canvas.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return canvas


def load_rows() -> list[dict[str, str]]:
    rows = list(csv.DictReader(PLAN.open(), delimiter="\t"))
    if len(rows) != 60:
        raise ValueError(f"Expected 60 assembly rows, found {len(rows)}")
    for row in rows:
        asset = ROOT / row["asset_or_action"]
        row["asset_abs"] = str(asset)
        if row["status"] != "ready":
            raise ValueError(f"slide {row['slide']} is not ready: {row['status']}")
        if not asset.exists():
            raise FileNotFoundError(asset)
        with Image.open(asset) as image:
            if image.size != EXPECTED_SIZE or image.mode != "RGB":
                raise ValueError(f"slide {row['slide']} has {image.size} {image.mode}: {asset}")
    return rows


def draw_card(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    row: dict[str, str],
    x: int,
    y: int,
    thumb: tuple[int, int],
    label_h: int,
) -> None:
    slide = int(row["slide"])
    act = row["act"]
    color = ACT_COLORS.get(act, COLORS["muted"])
    image = fit(Path(row["asset_abs"]), thumb)
    draw.rounded_rectangle((x - 5, y - 5, x + thumb[0] + 5, y + thumb[1] + 5), radius=16, fill=COLORS["panel"], outline=COLORS["line"], width=2)
    canvas.paste(image, (x, y))
    draw.rounded_rectangle((x + 14, y + 14, x + 82, y + 52), radius=15, fill=color)
    draw.text((x + 30, y + 19), f"{slide:02d}", font=F_NUM, fill="#081018")
    label_y = y + thumb[1] + 14
    draw.text((x, label_y), act.upper(), font=F_TINY, fill=color)
    for i, line in enumerate(wrap(draw, row["title"], F_LABEL, thumb[0] - 8, 2)):
        draw.text((x, label_y + 26 + i * 26), line, font=F_LABEL, fill=COLORS["text"])
    for i, line in enumerate(wrap(draw, row["main_question"], F_TINY, thumb[0] - 8, 2)):
        draw.text((x, label_y + 84 + i * 20), line, font=F_TINY, fill=COLORS["muted"])
    if label_h > 120:
        draw.text((x, label_y + 126), row["proof_object"], font=F_TINY, fill=COLORS["dim"])


def render_overview(rows: list[dict[str, str]]) -> None:
    thumb = (520, 293)
    label_h = 122
    cols = 4
    row_gap = 46
    col_gap = 48
    pad_x = 70
    header_h = 218
    footer_h = 95
    rows_n = (len(rows) + cols - 1) // cols
    width = pad_x * 2 + cols * thumb[0] + (cols - 1) * col_gap
    height = header_h + rows_n * (thumb[1] + label_h) + (rows_n - 1) * row_gap + footer_h

    canvas = Image.new("RGB", (width, height), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)
    draw.text((pad_x, 38), "SpaceBio-Bench detailed deck walkthrough QA", font=F_TITLE, fill=COLORS["text"])
    draw.text((pad_x, 104), "60-slide overview: rhythm, density, repeated grammar, and section handoffs before PPTX assembly.", font=F_SUB, fill=COLORS["muted"])

    legend_x = pad_x
    legend_y = 154
    for act, color in ACT_COLORS.items():
        w = max(78, int(draw.textlength(act, font=F_TINY)) + 34)
        draw.rounded_rectangle((legend_x, legend_y, legend_x + w, legend_y + 28), radius=12, fill="#121B28", outline=color, width=1)
        draw.text((legend_x + 16, legend_y + 7), act, font=F_TINY, fill=COLORS["text"])
        legend_x += w + 10

    for i, row in enumerate(rows):
        r, c = divmod(i, cols)
        x = pad_x + c * (thumb[0] + col_gap)
        y = header_h + r * (thumb[1] + label_h + row_gap)
        draw_card(canvas, draw, row, x, y, thumb, label_h)

    footer_y = height - footer_h + 24
    draw.line((pad_x, footer_y - 18, width - pad_x, footer_y - 18), fill=COLORS["line"], width=2)
    draw.text((pad_x, footer_y), "QA gate: no missing assets; all slides are 3840x2160 RGB; OCR flags written separately for visible-copy review.", font=F_SMALL, fill=COLORS["muted"])
    draw.text((width - pad_x, footer_y), "60 slides", font=F_LABEL, fill=COLORS["text"], anchor="ra")

    canvas.save(OVERVIEW, quality=95)
    canvas.convert("L").convert("RGB").save(OVERVIEW_GRAY, quality=95)


def render_section_sheets(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    specs = [(1, 20), (21, 40), (41, 60)]
    outputs: list[dict[str, str]] = []
    for start, end in specs:
        subset = [row for row in rows if start <= int(row["slide"]) <= end]
        thumb = (720, 405)
        label_h = 138
        cols = 2
        row_gap = 54
        col_gap = 60
        pad_x = 70
        header_h = 180
        footer_h = 90
        rows_n = (len(subset) + cols - 1) // cols
        width = pad_x * 2 + cols * thumb[0] + (cols - 1) * col_gap
        height = header_h + rows_n * (thumb[1] + label_h) + (rows_n - 1) * row_gap + footer_h
        canvas = Image.new("RGB", (width, height), COLORS["bg"])
        draw = ImageDraw.Draw(canvas)
        draw.text((pad_x, 36), f"Detailed deck walkthrough: slides {start}-{end}", font=F_TITLE, fill=COLORS["text"])
        draw.text((pad_x, 100), "Larger thumbnails for margin, overlap, text-size, and section-flow QA.", font=F_SUB, fill=COLORS["muted"])
        for i, row in enumerate(subset):
            r, c = divmod(i, cols)
            x = pad_x + c * (thumb[0] + col_gap)
            y = header_h + r * (thumb[1] + label_h + row_gap)
            draw_card(canvas, draw, row, x, y, thumb, label_h)
        footer_y = height - footer_h + 24
        draw.line((pad_x, footer_y - 18, width - pad_x, footer_y - 18), fill=COLORS["line"], width=2)
        draw.text((pad_x, footer_y), "Review these at full size for clipped labels, line-text collisions, and contrast loss.", font=F_SMALL, fill=COLORS["muted"])
        out = OUT_DIR / f"spacebiobench_detailed_deck_slides_{start:02d}_{end:02d}_walkthrough_contact_sheet.png"
        gray = OUT_DIR / f"spacebiobench_detailed_deck_slides_{start:02d}_{end:02d}_walkthrough_contact_sheet_grayscale.png"
        canvas.save(out, quality=95)
        canvas.convert("L").convert("RGB").save(gray, quality=95)
        outputs.append({"range": f"{start}-{end}", "contact_sheet": str(out.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT))})
    return outputs


def ocr_slide(row: dict[str, str]) -> str:
    src = Path(row["asset_abs"])
    downscaled = OCR_DIR / f"slide_{int(row['slide']):02d}_ocr_input.png"
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


def run_ocr(rows: list[dict[str, str]]) -> tuple[list[dict[str, object]], list[dict[str, str]]]:
    all_text: list[str] = []
    flags: list[dict[str, str]] = []
    slide_texts: list[dict[str, object]] = []
    for row in rows:
        slide = int(row["slide"])
        text_out = ocr_slide(row)
        slide_texts.append({"slide": slide, "title": row["title"], "text": text_out})
        all_text.append(f"--- slide {slide:02d}: {row['title']} ---\n{text_out.strip()}\n")
        for label, pattern in FLAG_PATTERNS:
            matches = sorted(set(match.group(0) for match in pattern.finditer(text_out)))
            if matches:
                flags.append(
                    {
                        "slide": f"{slide:02d}",
                        "title": row["title"],
                        "category": label,
                        "matches": ", ".join(matches),
                    }
                )
    OCR_TEXT.write_text("\n".join(all_text), encoding="utf-8")
    with OCR_FLAGS.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["slide", "title", "category", "matches"], delimiter="\t")
        writer.writeheader()
        writer.writerows(flags)
    return slide_texts, flags


def write_report(rows: list[dict[str, str]], section_outputs: list[dict[str, str]], flags: list[dict[str, str]]) -> None:
    flag_counts: dict[str, int] = {}
    for flag in flags:
        flag_counts[flag["category"]] = flag_counts.get(flag["category"], 0) + 1
    section_counts: dict[str, int] = {}
    for row in rows:
        section_counts[row["act"]] = section_counts.get(row["act"], 0) + 1
    lines = [
        "# Detailed Deck Full-Walkthrough QA",
        "",
        "## Asset Gate",
        "",
        f"- Slides checked: {len(rows)}",
        "- Status: all ready",
        "- Dimensions: all 3840 x 2160 RGB",
        "- Missing assets: 0",
        "",
        "## Contact Sheets",
        "",
        f"- 60-slide overview: `{OVERVIEW.relative_to(ROOT)}`",
        f"- 60-slide grayscale overview: `{OVERVIEW_GRAY.relative_to(ROOT)}`",
    ]
    for item in section_outputs:
        lines.append(f"- Slides {item['range']}: `{item['contact_sheet']}`")
        lines.append(f"- Slides {item['range']} grayscale: `{item['grayscale']}`")
    lines.extend(["", "## Section Rhythm", ""])
    for act, count in section_counts.items():
        lines.append(f"- {act}: {count} slide(s)")
    lines.extend(["", "## OCR Visible-Copy Flags", ""])
    if not flags:
        lines.append("- No flagged visible-copy terms found by OCR.")
    else:
        for category, count in sorted(flag_counts.items()):
            lines.append(f"- {category}: {count}")
        lines.append("")
        lines.append("Flagged slides:")
        for flag in flags:
            lines.append(f"- Slide {flag['slide']}: {flag['category']} -> {flag['matches']} | {flag['title']}")
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    rows = load_rows()
    render_overview(rows)
    section_outputs = render_section_sheets(rows)
    _slide_texts, flags = run_ocr(rows)
    write_report(rows, section_outputs, flags)
    MANIFEST.write_text(
        json.dumps(
            {
                "overview_contact_sheet": str(OVERVIEW.relative_to(ROOT)),
                "overview_grayscale": str(OVERVIEW_GRAY.relative_to(ROOT)),
                "section_contact_sheets": section_outputs,
                "ocr_text": str(OCR_TEXT.relative_to(ROOT)),
                "ocr_flags": str(OCR_FLAGS.relative_to(ROOT)),
                "qa_report": str(REPORT.relative_to(ROOT)),
                "slides": [
                    {
                        "slide": int(row["slide"]),
                        "act": row["act"],
                        "title": row["title"],
                        "asset": row["asset_or_action"],
                    }
                    for row in rows
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "overview": str(OVERVIEW.relative_to(ROOT)),
                "sections": [item["contact_sheet"] for item in section_outputs],
                "ocr_flags": str(OCR_FLAGS.relative_to(ROOT)),
                "report": str(REPORT.relative_to(ROOT)),
                "flag_count": len(flags),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

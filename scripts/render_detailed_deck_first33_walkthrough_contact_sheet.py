#!/usr/bin/env python3
"""Render the canonical first-33 detailed-deck walkthrough contact sheet."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, ImageOps


ROOT = Path(__file__).resolve().parent.parent
WORKSPACE = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
)
PLAN = WORKSPACE / "spacebiobench_detailed_deck_assembly_plan_v2.tsv"
OUT_DIR = WORKSPACE / "assets" / "first33_walkthrough_assembly"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CONTACT = OUT_DIR / "slides_1_33_walkthrough_contact_sheet.png"
GRAY = OUT_DIR / "slides_1_33_walkthrough_contact_sheet_grayscale.png"
MANIFEST = OUT_DIR / "slides_1_33_walkthrough_contact_sheet_manifest.json"

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "line": "#2A394D",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#667487",
    "teal": "#5FD3C4",
    "blue": "#66A6E8",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "violet": "#B39DFF",
    "rose": "#E17882",
}

ACT_COLORS = {
    "Open": COLORS["blue"],
    "Method": COLORS["teal"],
    "Core Result": COLORS["green"],
    "Models": COLORS["violet"],
    "Robustness": COLORS["amber"],
}

SECTION_BREAKS = {
    1: "OPENING CONTRACT",
    4: "METHOD SETUP",
    13: "CORE RESULT",
    22: "MODEL COMPARISON",
    30: "ROBUSTNESS + TAKEAWAY",
}


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


F_TITLE = font(52, bold=True)
F_SUB = font(23)
F_LABEL = font(21, bold=True)
F_SMALL = font(17)
F_TINY = font(15)
F_ACT = font(17, bold=True)
F_SECTION = font(18, bold=True)


def wrap(draw: ImageDraw.ImageDraw, value: str, fnt: ImageFont.ImageFont, width: int) -> list[str]:
    words = value.split()
    lines: list[str] = []
    current: list[str] = []
    for word in words:
        trial = " ".join(current + [word])
        if draw.textlength(trial, font=fnt) <= width:
            current.append(word)
            continue
        if current:
            lines.append(" ".join(current))
        current = [word]
    if current:
        lines.append(" ".join(current))
    return lines


def fit(path: Path, size: tuple[int, int]) -> Image.Image:
    image = Image.open(path).convert("RGB")
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    image.thumbnail(size, resample)
    canvas = Image.new("RGB", size, COLORS["bg"])
    canvas.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return canvas


def load_rows() -> list[dict[str, str]]:
    rows = list(csv.DictReader(PLAN.open(encoding="utf-8"), delimiter="\t"))
    first33 = [row for row in rows if 1 <= int(row["slide"]) <= 33]
    if len(first33) != 33:
        raise ValueError(f"expected 33 rows, found {len(first33)}")
    for row in first33:
        if row["status"] != "ready":
            raise ValueError(f"slide {row['slide']} is not ready: {row['status']}")
        asset = ROOT / row["asset_or_action"]
        if not asset.exists():
            raise FileNotFoundError(asset)
        with Image.open(asset) as image:
            if image.size != (3840, 2160):
                raise ValueError(f"slide {row['slide']} has unexpected size {image.size}: {asset}")
    return first33


def draw_card(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    row: dict[str, str],
    x: int,
    y: int,
    thumb: tuple[int, int],
) -> dict[str, object]:
    slide = int(row["slide"])
    act = row["act"]
    color = ACT_COLORS.get(act, COLORS["teal"])
    image = fit(ROOT / row["asset_or_action"], thumb)

    draw.rounded_rectangle(
        (x - 5, y - 5, x + thumb[0] + 5, y + thumb[1] + 5),
        radius=18,
        fill=COLORS["panel"],
        outline=COLORS["line"],
        width=2,
    )
    canvas.paste(image, (x, y))
    label_y = y + thumb[1] + 16
    draw.rounded_rectangle((x, label_y, x + 82, label_y + 34), radius=12, fill="#132132", outline=color, width=2)
    draw.text((x + 14, label_y + 8), f"{slide:02d}", font=F_LABEL, fill=COLORS["text"])
    draw.text((x + 100, label_y + 8), act.upper(), font=F_ACT, fill=color)

    title_y = label_y + 42
    for i, line in enumerate(wrap(draw, row["title"], F_LABEL, thumb[0] - 8)[:2]):
        draw.text((x, title_y + i * 26), line, font=F_LABEL, fill=COLORS["text"])
    question_y = title_y + 58
    for i, line in enumerate(wrap(draw, row["main_question"], F_TINY, thumb[0] - 8)[:2]):
        draw.text((x, question_y + i * 20), line, font=F_TINY, fill=COLORS["muted"])

    return {
        "slide": slide,
        "act": row["act"],
        "title": row["title"],
        "main_question": row["main_question"],
        "path": row["asset_or_action"],
    }


def main() -> None:
    rows = load_rows()
    thumb = (590, 332)
    cols = 3
    card_h = thumb[1] + 128
    pad_x = 72
    gap_x = 54
    gap_y = 50
    header_h = 255
    footer_h = 128
    rows_n = (len(rows) + cols - 1) // cols
    width = pad_x * 2 + cols * thumb[0] + (cols - 1) * gap_x
    height = header_h + rows_n * card_h + (rows_n - 1) * gap_y + footer_h

    canvas = Image.new("RGB", (width, height), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)
    draw.text((pad_x, 36), "SpaceBio-Bench detailed deck: first 33-slide walkthrough", font=F_TITLE, fill=COLORS["text"])
    draw.text(
        (pad_x, 102),
        "Canonical contact-sheet QA after adding slides 1-6: thesis -> method -> result -> models -> robustness.",
        font=F_SUB,
        fill=COLORS["muted"],
    )

    legend_x = pad_x
    for act, color in ACT_COLORS.items():
        draw.rounded_rectangle((legend_x, 150, legend_x + 160, 184), radius=12, fill="#132132", outline=color, width=2)
        draw.text((legend_x + 14, 159), act, font=F_TINY, fill=COLORS["text"])
        legend_x += 176

    draw.text(
        (pad_x, 210),
        "QA focus: flow, rhythm, visual density, recurring reading rules, section transitions, grayscale hierarchy, and slide-to-slide cognitive load.",
        font=F_SMALL,
        fill=COLORS["dim"],
    )

    entries = []
    for i, row in enumerate(rows):
        slide = int(row["slide"])
        r, c = divmod(i, cols)
        x = pad_x + c * (thumb[0] + gap_x)
        y = header_h + r * (card_h + gap_y)
        if c == 0 and slide in SECTION_BREAKS:
            section_color = ACT_COLORS.get(row["act"], COLORS["teal"])
            draw.rounded_rectangle((x, y - 32, x + 430, y - 6), radius=10, fill="#132132", outline=section_color, width=1)
            draw.text((x + 14, y - 27), SECTION_BREAKS[slide], font=F_SECTION, fill=section_color)
        entries.append(draw_card(canvas, draw, row, x, y, thumb))

    footer_y = height - footer_h + 26
    draw.line((pad_x, footer_y - 20, width - pad_x, footer_y - 20), fill=COLORS["line"], width=2)
    draw.text((pad_x, footer_y), "Ready assets: 33 / 33", font=F_LABEL, fill=COLORS["teal"])
    draw.text(
        (pad_x + 235, footer_y + 3),
        "Next gate: visual QA notes for the first-33 sequence, then PPT assembly / export when the sequence is approved.",
        font=F_SMALL,
        fill=COLORS["muted"],
    )
    draw.text((width - pad_x, footer_y), "slides 1-33", font=F_LABEL, fill=COLORS["teal"], anchor="ra")

    canvas.save(CONTACT, quality=95)
    ImageOps.grayscale(canvas).convert("RGB").save(GRAY, quality=95)
    MANIFEST.write_text(
        json.dumps(
            {
                "contact_sheet": str(CONTACT.relative_to(ROOT)),
                "grayscale": str(GRAY.relative_to(ROOT)),
                "scope": "canonical first 33-slide detailed-deck walkthrough; all assets ready at 3840x2160",
                "slides": entries,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"contact_sheet": str(CONTACT.relative_to(ROOT)), "grayscale": str(GRAY.relative_to(ROOT)), "manifest": str(MANIFEST.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

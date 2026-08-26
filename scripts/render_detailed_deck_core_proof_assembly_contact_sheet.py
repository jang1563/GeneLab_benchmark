#!/usr/bin/env python3
"""Render ordered assembly QA contact sheets for ready detailed-deck proof assets."""

from __future__ import annotations

import csv
import json
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
OUT_DIR = WORKSPACE / "assets" / "core_proof_assembly"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CONTACT = OUT_DIR / "slides_7_33_core_proof_assembly_contact_sheet.png"
GRAY = OUT_DIR / "slides_7_33_core_proof_assembly_contact_sheet_grayscale.png"
MANIFEST = OUT_DIR / "slides_7_33_core_proof_assembly_contact_sheet_manifest.json"


COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "grid": "#2A3546",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#667487",
    "teal": "#5FD3C4",
    "blue": "#66A6E8",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "violet": "#B39DFF",
}

ACT_COLORS = {
    "Method": COLORS["teal"],
    "Core Result": COLORS["green"],
    "Models": COLORS["violet"],
    "Robustness": COLORS["amber"],
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


F_TITLE = font(48, bold=True)
F_SUB = font(21)
F_LABEL = font(21, bold=True)
F_SMALL = font(17)
F_TINY = font(15)
F_ACT = font(17, bold=True)


def wrap(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.ImageFont, width: int) -> list[str]:
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
    return lines


def fit(path: Path, size: tuple[int, int]) -> Image.Image:
    image = Image.open(path).convert("RGB")
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    image.thumbnail(size, resample)
    canvas = Image.new("RGB", size, COLORS["bg"])
    canvas.paste(image, ((size[0] - image.width) // 2, (size[1] - image.height) // 2))
    return canvas


def load_rows() -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    rows = list(csv.DictReader(PLAN.open(), delimiter="\t"))
    pending = [row for row in rows if int(row["slide"]) <= 6]
    ready = [row for row in rows if 7 <= int(row["slide"]) <= 33]
    for row in ready:
        if row["status"] != "ready":
            raise ValueError(f"slide {row['slide']} is not ready: {row['status']}")
        asset = ROOT / row["asset_or_action"]
        if not asset.exists():
            raise FileNotFoundError(asset)
    return pending, ready


def draw_card(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    row: dict[str, str],
    x: int,
    y: int,
    thumb: tuple[int, int],
) -> None:
    slide = int(row["slide"])
    act = row["act"]
    color = ACT_COLORS.get(act, COLORS["teal"])
    image = fit(ROOT / row["asset_or_action"], thumb)
    draw.rounded_rectangle((x - 5, y - 5, x + thumb[0] + 5, y + thumb[1] + 5), radius=18, fill=COLORS["panel"], outline=COLORS["grid"], width=2)
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


def main() -> None:
    pending, ready = load_rows()
    thumb = (620, 349)
    cols = 3
    card_h = thumb[1] + 126
    pad_x = 70
    gap_x = 58
    gap_y = 48
    header_h = 210
    footer_h = 110
    rows_n = (len(ready) + cols - 1) // cols
    width = pad_x * 2 + cols * thumb[0] + (cols - 1) * gap_x
    height = header_h + rows_n * card_h + (rows_n - 1) * gap_y + footer_h

    canvas = Image.new("RGB", (width, height), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)
    draw.text((pad_x, 36), "Detailed deck core proof-asset assembly QA", font=F_TITLE, fill=COLORS["text"])
    draw.text((pad_x, 96), "Ordered sequence review for slides 7-33. Slides 1-6 remain opening/method build items and are outside this proof-asset gate.", font=F_SUB, fill=COLORS["muted"])
    ready_note = f"Ready proof assets: {len(ready)} / 27 | Pending opening/method setup: " + ", ".join(f"{r['slide']} {r['status']}" for r in pending)
    draw.text((pad_x, 132), ready_note, font=F_SMALL, fill=COLORS["dim"])

    legend_x = pad_x
    for act, color in ACT_COLORS.items():
        draw.rounded_rectangle((legend_x, 164, legend_x + 148, 196), radius=12, fill="#132132", outline=color, width=2)
        draw.text((legend_x + 14, 173), act, font=F_TINY, fill=COLORS["text"])
        legend_x += 164

    entries = []
    for i, row in enumerate(ready):
        r, c = divmod(i, cols)
        x = pad_x + c * (thumb[0] + gap_x)
        y = header_h + r * (card_h + gap_y)
        draw_card(canvas, draw, row, x, y, thumb)
        entries.append(
            {
                "slide": int(row["slide"]),
                "act": row["act"],
                "title": row["title"],
                "path": row["asset_or_action"],
            }
        )

    footer_y = height - footer_h + 24
    draw.line((pad_x, footer_y - 20, width - pad_x, footer_y - 20), fill=COLORS["grid"], width=2)
    draw.text((pad_x, footer_y), "QA focus: rhythm, slide-to-slide logic, thumbnail distinctness, repeated visual grammar, and readiness for first 33-slide assembly.", font=F_SMALL, fill=COLORS["muted"])
    draw.text((width - pad_x, footer_y), "slides 7-33", font=F_LABEL, fill=COLORS["teal"], anchor="ra")

    canvas.save(CONTACT, quality=95)
    canvas.convert("L").convert("RGB").save(GRAY, quality=95)
    MANIFEST.write_text(
        json.dumps(
            {
                "contact_sheet": str(CONTACT.relative_to(ROOT)),
                "grayscale": str(GRAY.relative_to(ROOT)),
                "scope": "slides 7-33 ready proof assets; slides 1-6 pending opening/method setup",
                "pending_slides_1_6": [
                    {"slide": int(row["slide"]), "title": row["title"], "status": row["status"], "asset_or_action": row["asset_or_action"]}
                    for row in pending
                ],
                "slides": entries,
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps({"contact_sheet": str(CONTACT.relative_to(ROOT)), "grayscale": str(GRAY.relative_to(ROOT)), "manifest": str(MANIFEST.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

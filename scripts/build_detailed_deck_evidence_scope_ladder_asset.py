#!/usr/bin/env python3
"""Build the detailed-deck evidence scope ladder teaching asset."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "evidence_scope_ladder"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "panel2": "#151F2D",
    "grid": "#2A3546",
    "axis": "#98A7BA",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#667487",
    "teal": "#5FD3C4",
    "sky": "#73A7FF",
    "green": "#84D278",
    "amber": "#F4C26B",
    "rose": "#E17882",
    "violet": "#B39DFF",
    "white": "#FFFFFF",
}


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(value: str, alpha: int) -> tuple[int, int, int, int]:
    return (*hex_to_rgb(value), alpha)


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Helvetica.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            continue
    return ImageFont.load_default()


F = {
    "kicker": load_font(34, True),
    "title": load_font(82, True),
    "subtitle": load_font(36, False),
    "h2": load_font(44, True),
    "h3": load_font(34, True),
    "body": load_font(30, False),
    "body_bold": load_font(30, True),
    "small": load_font(25, False),
    "small_bold": load_font(25, True),
    "tiny": load_font(21, False),
    "tiny_bold": load_font(21, True),
    "number": load_font(58, True),
}


def rounded(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], radius: int, fill: str, outline: str | None = None, width: int = 1) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def text(draw: ImageDraw.ImageDraw, xy: tuple[int | float, int | float], value: str, font: ImageFont.ImageFont, fill: str = COLORS["text"], anchor: str | None = None) -> None:
    draw.text(xy, value, font=font, fill=fill, anchor=anchor)


def wrap_lines(draw: ImageDraw.ImageDraw, body: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
    words = body.split()
    lines: list[str] = []
    cur: list[str] = []
    for word in words:
        trial = " ".join(cur + [word])
        if draw.textlength(trial, font=font) <= max_width:
            cur.append(word)
        else:
            if cur:
                lines.append(" ".join(cur))
            cur = [word]
    if cur:
        lines.append(" ".join(cur))
    return lines


def multiline(draw: ImageDraw.ImageDraw, xy: tuple[int, int], lines: Iterable[str], font: ImageFont.ImageFont, fill: str, leading: int = 8) -> int:
    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=font, fill=fill)
        y += font.size + leading
    return y


def paragraph(draw: ImageDraw.ImageDraw, xy: tuple[int, int], body: str, font: ImageFont.ImageFont, max_width: int, fill: str, leading: int = 8) -> int:
    return multiline(draw, xy, wrap_lines(draw, body, font, max_width), font, fill, leading)


def draw_badges(draw: ImageDraw.ImageDraw) -> None:
    badges = [
        ("AUDIENCE GUIDE", "ML + biology readers", 430),
        ("LAYERS", "task to translation", 365),
        ("READOUT", "scope before score", 390),
    ]
    bx = 2250
    for kicker, body, badge_w in badges:
        rounded(draw, (bx, 72, bx + badge_w, 176), 26, "#121B28", "#2A394D", 2)
        text(draw, (bx + 28, 96), kicker, F["tiny"], COLORS["sky"])
        text(draw, (bx + 28, 126), body, F["small_bold"], COLORS["text"])
        bx += badge_w + 30


def draw_chip(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, color: str) -> int:
    pad_x = 20
    w = int(draw.textlength(label, font=F["tiny_bold"]) + pad_x * 2)
    rounded(draw, (x, y, x + w, y + 42), 18, "#172231", color, 2)
    text(draw, (x + pad_x, y + 10), label, F["tiny_bold"], COLORS["text"])
    return w


def draw_reader_contract(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "A. Reader contract", F["h2"], COLORS["text"])
    paragraph(draw, (x0 + 50, y0 + 100), "Before reading a result, identify the task, the hidden unit, and the evidence layer.", F["small"], x1 - x0 - 100, COLORS["muted"], 8)

    rows = [
        ("1", "Task fixed", "input view, label, split, and metric are defined", COLORS["teal"]),
        ("2", "Hidden unit", "the evaluation unit is one unseen mission", COLORS["sky"]),
        ("3", "Readout layer", "score, robustness, biology, and follow-up stay distinct", COLORS["amber"]),
    ]
    y = y0 + 255
    for num, title, body, color in rows:
        rounded(draw, (x0 + 50, y, x1 - 50, y + 145), 24, "#151F2D", "#2A394D", 2)
        draw.ellipse((x0 + 82, y + 40, x0 + 132, y + 90), fill=rgba(color, 230))
        text(draw, (x0 + 107, y + 48), num, F["tiny_bold"], "#081018", anchor="ma")
        text(draw, (x0 + 160, y + 30), title, F["h3"], COLORS["text"])
        paragraph(draw, (x0 + 160, y + 78), body, F["small"], x1 - x0 - 230, COLORS["muted"], 6)
        y += 175

    rounded(draw, (x0 + 50, y1 - 190, x1 - 50, y1 - 48), 24, "#211E17", "#69532B", 2)
    text(draw, (x0 + 82, y1 - 155), "Use this slide", F["small_bold"], COLORS["amber"])
    paragraph(draw, (x0 + 82, y1 - 112), "It is the legend for the next result slides: read the layer first, then read the number.", F["small"], x1 - x0 - 165, COLORS["muted"], 6)


def draw_layer_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    number: str,
    title: str,
    question: str,
    chips: list[str],
    readout: str,
    color: str,
) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 28, COLORS["panel"], "#29374A", 2)
    draw.ellipse((x0 + 36, y0 + 42, x0 + 126, y0 + 132), fill=rgba(color, 235))
    text(draw, (x0 + 81, y0 + 56), number, F["number"], "#081018", anchor="ma")
    text(draw, (x0 + 160, y0 + 36), title, F["h2"], COLORS["text"])
    paragraph(draw, (x0 + 160, y0 + 91), question, F["small"], x1 - x0 - 230, COLORS["muted"], 6)

    chip_x = x0 + 160
    chip_y = y0 + 175
    for chip in chips:
        chip_w = draw_chip(draw, chip_x, chip_y, chip, color)
        chip_x += chip_w + 18
        if chip_x > x1 - 275:
            chip_x = x0 + 160
            chip_y += 56

    rounded(draw, (x0 + 160, y1 - 70, x1 - 42, y1 - 24), 18, "#151F2D", "#2A394D", 1)
    text(draw, (x0 + 185, y1 - 59), readout, F["small_bold"], color)


def draw_ladder(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    text(draw, (x0, y0 - 58), "B. Evidence ladder", F["h2"], COLORS["text"])
    text(draw, (x0, y0 - 15), "Each later result belongs to one primary evidence layer.", F["small"], COLORS["muted"])
    layers = [
        (
            "1",
            "Benchmark task evidence",
            "Was an unseen mission scored under a fixed task contract?",
            ["hidden mission", "train-only pipeline", "AUROC + CI"],
            "transfer readout",
            COLORS["teal"],
        ),
        (
            "2",
            "Robustness evidence",
            "Does the signal persist across controls and method choices?",
            ["negative controls", "reserved missions", "DGE callers"],
            "reliability readout",
            COLORS["sky"],
        ),
        (
            "3",
            "Biological support",
            "Does known spaceflight biology point in the same direction?",
            ["NES conservation", "Cell 2020 pathways", "SHAP genes"],
            "plausibility readout",
            COLORS["green"],
        ),
        (
            "4",
            "Follow-up layer",
            "Which findings guide biology, translation, or release work next?",
            ["temporal context", "mouse-human bridge", "payload readiness"],
            "follow-up map",
            COLORS["amber"],
        ),
    ]
    card_h = 305
    gap = 38
    for i, layer in enumerate(layers):
        y = y0 + i * (card_h + gap)
        draw_layer_card(draw, (x0, y, x1, y + card_h), *layer)
        if i < len(layers) - 1:
            cx = x0 + 80
            draw.line((cx, y + card_h + 6, cx, y + card_h + gap - 6), fill=rgba(layer[-1], 170), width=6)


def draw_unlock_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "C. What this unlocks", F["h2"], COLORS["text"])
    paragraph(draw, (x0 + 50, y0 + 100), "The deck can be dense because every section declares its readout layer before the key visual.", F["small"], x1 - x0 - 100, COLORS["muted"], 8)

    rows = [
        ("Methods", "task, split, metric, model inputs", COLORS["teal"]),
        ("Core result", "tissue transfer and pathway conservation", COLORS["green"]),
        ("Robustness", "holdout, DGE callers, external biology", COLORS["sky"]),
        ("Extension", "biology, translation, platform readiness", COLORS["amber"]),
    ]
    y = y0 + 245
    for label, body, color in rows:
        rounded(draw, (x0 + 50, y, x1 - 50, y + 138), 24, "#151F2D", "#2A394D", 2)
        draw.rectangle((x0 + 50, y, x0 + 62, y + 138), fill=color)
        text(draw, (x0 + 88, y + 28), label, F["h3"], color)
        paragraph(draw, (x0 + 88, y + 75), body, F["small"], x1 - x0 - 165, COLORS["text"], 6)
        y += 166

    rounded(draw, (x0 + 50, y1 - 245, x1 - 50, y1 - 48), 28, "#211E17", "#69532B", 2)
    text(draw, (x0 + 82, y1 - 205), "Reading move", F["h3"], COLORS["amber"])
    paragraph(draw, (x0 + 82, y1 - 154), "Ask which layer owns the result. The deck keeps score, robustness, biology, and translation in separate lanes.", F["small"], x1 - x0 - 165, COLORS["muted"], 8)


def main() -> None:
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 52), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 42), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "READOUT PRIMER", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "Evidence layers keep the benchmark readable", F["title"], COLORS["text"])
    text(draw, (150, 216), "Scores, robustness checks, biological support, and follow-up layers answer different questions.", F["subtitle"], COLORS["muted"])
    draw_badges(draw)

    draw_reader_contract(draw, (150, 350, 1005, 1800))
    draw_ladder(draw, (1080, 410, 2460, 1800))
    draw_unlock_panel(draw, (2535, 350, 3690, 1800))

    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    text(draw, (205, 1925), "Takeaway", F["small_bold"], COLORS["sky"])
    source = "Scores, robustness, biology, and follow-up layers answer different questions and should not be collapsed."
    paragraph(draw, (390, 1925), source, F["small"], 3220, COLORS["muted"], 6)
    text(draw, (205, 1995), "Next", F["small_bold"], COLORS["amber"])
    scope = "The core result section starts by reading one tissue score, then scales up to tissue ranking and pathway conservation."
    paragraph(draw, (390, 1995), scope, F["small"], 3220, COLORS["muted"], 6)

    png = OUT_DIR / "evidence_scope_ladder_premium.png"
    canvas.save(png, quality=95)
    gray = OUT_DIR / "evidence_scope_ladder_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    manifest = {
        "title": "Evidence layers keep the benchmark readable",
        "outputs": {"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT))},
        "slide_role": "Readout primer before result slides",
        "layers": [
            "Benchmark task evidence",
            "Robustness evidence",
            "Biological support",
            "Follow-up layer",
        ],
        "scope": "Guide for interpreting later result slides by evidence layer.",
    }
    manifest_path = OUT_DIR / "evidence_scope_ladder_premium_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

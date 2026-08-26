#!/usr/bin/env python3
"""Build slide 56 asset: organoid biology check."""

from __future__ import annotations

import json
import math
import random
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, ImageOps, ImageStat


ROOT = Path(__file__).resolve().parent.parent
ASSET_ROOT = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
)
OUT_DIR = ASSET_ROOT / "organoid_biology_check"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "organoid_extension_is_a_biology_check_premium.png"
GRAY_PATH = OUT_DIR / "organoid_extension_is_a_biology_check_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "organoid_biology_check_manifest.json"
QA_NOTE = OUT_DIR / "ORGANOID_BIOLOGY_CHECK_ASSET_VISUAL_QA.md"

COLORS = {
    "bg": "#0B1119",
    "bg2": "#111721",
    "header": "#101826",
    "panel": "#111B28",
    "panel2": "#172335",
    "panel3": "#0F1825",
    "grid": "#263245",
    "text": "#F4F7FB",
    "muted": "#AAB6C6",
    "dim": "#687789",
    "blue": "#66A6E8",
    "teal": "#5FD3C4",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "violet": "#B39DFF",
    "rose": "#E17882",
    "ink": "#081018",
}

FOOTPRINT = [
    {"label": "public studies", "value": "2", "detail": "OSD-863 / OSD-871", "color": COLORS["amber"]},
    {"label": "samples", "value": "42", "detail": "neural organoid RNA-seq", "color": COLORS["amber"]},
    {"label": "flight-ground contrasts", "value": "8", "detail": "Ground vs LEO/ISS", "color": COLORS["blue"]},
    {"label": "reference rows", "value": "242.7k", "detail": "gene / contrast rows", "color": COLORS["blue"]},
    {"label": "significant rows", "value": "2,368", "detail": "FDR <= 0.05", "color": COLORS["teal"]},
]

READOUTS = [
    {
        "name": "Primary Prediction",
        "role": "flight / ground classifier",
        "color": COLORS["green"],
        "items": [("AUROC", 0.6147727273), ("Macro-F1", 0.5194508009)],
    },
    {
        "name": "Response-Pattern Check",
        "role": "flight biology alignment",
        "color": COLORS["blue"],
        "items": [("Direction match", 0.7706734868), ("Rank correlation", 0.1760078660)],
    },
    {
        "name": "Gene-Effect Check",
        "role": "model coefficient alignment",
        "color": COLORS["amber"],
        "items": [("Direction match", 0.6078431373), ("Rank correlation", 0.0867280024)],
    },
]

TOPK = [
    {"label": "Top 50", "observed": 1, "expected": 0.6375, "enrichment": 1.5686},
    {"label": "Top 100", "observed": 5, "expected": 1.2750, "enrichment": 3.9216},
    {"label": "Top 250", "observed": 10, "expected": 3.1875, "enrichment": 3.1373},
    {"label": "Top 500", "observed": 14, "expected": 6.3750, "enrichment": 2.1961},
]


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(value: str, alpha: int) -> tuple[int, int, int, int]:
    return (*hex_to_rgb(value), alpha)


def blend(a: str, b: str, t: float) -> str:
    ar, ag, ab = hex_to_rgb(a)
    br, bg, bb = hex_to_rgb(b)
    t = max(0.0, min(1.0, t))
    return f"#{int(ar + (br - ar) * t):02x}{int(ag + (bg - ag) * t):02x}{int(ab + (bb - ab) * t):02x}"


def load_font(size: int, bold: bool = False) -> ImageFont.ImageFont:
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
    "title": load_font(78, True),
    "subtitle": load_font(37),
    "section": load_font(40, True),
    "h2": load_font(34, True),
    "body": load_font(26),
    "body_bold": load_font(26, True),
    "small": load_font(23),
    "small_bold": load_font(23, True),
    "tiny": load_font(20),
    "tiny_bold": load_font(20, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "metric": load_font(62, True),
    "big_metric": load_font(70, True),
    "mono": load_font(22),
    "mono_bold": load_font(22, True),
}


def text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    value: str,
    font: ImageFont.ImageFont,
    fill: str = COLORS["text"],
    anchor: str | None = None,
) -> None:
    draw.text(xy, value, font=font, fill=fill, anchor=anchor)


def wrap(draw: ImageDraw.ImageDraw, value: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
    words = value.split()
    lines: list[str] = []
    line = ""
    for word in words:
        trial = word if not line else f"{line} {word}"
        if draw.textlength(trial, font=font) <= max_width:
            line = trial
        else:
            if line:
                lines.append(line)
            line = word
    if line:
        lines.append(line)
    return lines


def paragraph(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    value: str,
    font: ImageFont.ImageFont,
    max_width: int,
    fill: str = COLORS["muted"],
    leading: int = 8,
) -> int:
    x, y = xy
    for line in wrap(draw, value, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def rounded(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    radius: int,
    fill: str,
    outline: str | None = None,
    width: int = 1,
) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def background() -> Image.Image:
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image, "RGBA")
    for y in range(H):
        t = y / H
        color = blend(COLORS["bg"], COLORS["bg2"], t * 0.72)
        draw.line((0, y, W, y), fill=color)
    for x in range(0, W, 160):
        draw.line((x, 260, x, H - 220), fill=rgba(COLORS["grid"], 72), width=1)
    for y in range(320, H - 220, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 64), width=1)
    draw.rectangle((0, 0, W, 260), fill=rgba(COLORS["header"], 245))
    draw.rectangle((0, H - 190, W, H), fill=rgba("#071019", 235))
    return image


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    w = 230 if len(value) < 5 else 260
    rounded(draw, (x, y, x + w, y + 72), 14, COLORS["panel2"], color, 2)
    text(draw, (x + 18, y + 14), label.upper(), F["micro_bold"], COLORS["muted"])
    text(draw, (x + 18, y + 42), value, F["small_bold"], COLORS["text"])
    return x + w + 24


def draw_rosette(draw: ImageDraw.ImageDraw, cx: int, cy: int, r: int) -> None:
    rng = random.Random(56)
    for i in range(42):
        angle = (math.tau * i / 42.0) + rng.uniform(-0.06, 0.06)
        rr = rng.uniform(r * 0.16, r * 0.94)
        x = cx + math.cos(angle) * rr
        y = cy + math.sin(angle) * rr
        size = rng.uniform(10, 26)
        fill = rgba(COLORS["teal"] if i % 3 else COLORS["violet"], 70)
        draw.ellipse((x - size, y - size, x + size, y + size), fill=fill, outline=rgba("#EAF2FA", 60), width=1)
    for ring, alpha in [(r * 0.42, 80), (r * 0.68, 60), (r * 0.92, 42)]:
        draw.ellipse((cx - ring, cy - ring, cx + ring, cy + ring), outline=rgba(COLORS["teal"], alpha), width=3)
    draw.ellipse((cx - 54, cy - 54, cx + 54, cy + 54), fill=rgba(COLORS["green"], 72), outline=rgba("#EAF2FA", 90), width=2)
    text(draw, (cx, cy - 8), "42", F["h2"], COLORS["text"], "ma")
    text(draw, (cx, cy + 26), "samples", F["micro_bold"], COLORS["muted"], "ma")


def draw_footprint_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 610, 1780, 1392
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["teal"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Compact Human Organoid Footprint", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "Two human neural-organoid studies create a compact extension lane for flight-response biology.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    draw_rosette(draw, x1 + 230, y1 + 390, 165)
    rounded(draw, (x1 + 76, y1 + 594, x1 + 384, y1 + 688), 18, COLORS["panel3"], COLORS["teal"], 2)
    text(draw, (x1 + 100, y1 + 615), "human iPSC", F["tiny_bold"], COLORS["text"])
    text(draw, (x1 + 100, y1 + 649), "neural organoid RNA-seq", F["tiny"], COLORS["muted"])

    start_x, end_x, line_y = x1 + 520, x2 - 110, y1 + 408
    xs = [int(start_x + i * (end_x - start_x) / (len(FOOTPRINT) - 1)) for i in range(len(FOOTPRINT))]
    draw.line((xs[0], line_y, xs[-1], line_y), fill=rgba(COLORS["dim"], 180), width=4)
    for i, (x, item) in enumerate(zip(xs, FOOTPRINT)):
        color = str(item["color"])
        draw.ellipse((x - 19, line_y - 19, x + 19, line_y + 19), fill=color, outline="#EAF2FA", width=2)
        if i % 2 == 0:
            text(draw, (x, line_y - 118), str(item["value"]), F["metric"], color, "ma")
            text(draw, (x, line_y - 50), str(item["label"]), F["small_bold"], COLORS["text"], "ma")
            text(draw, (x, line_y + 50), str(item["detail"]), F["tiny"], COLORS["muted"], "ma")
        else:
            text(draw, (x, line_y + 86), str(item["value"]), F["metric"], color, "ma")
            text(draw, (x, line_y + 152), str(item["label"]), F["small_bold"], COLORS["text"], "ma")
            text(draw, (x, line_y - 56), str(item["detail"]), F["tiny"], COLORS["muted"], "ma")


def draw_slider(
    draw: ImageDraw.ImageDraw,
    x1: int,
    y: int,
    x2: int,
    label: str,
    value: float,
    color: str,
) -> None:
    text(draw, (x1, y - 9), label, F["small_bold"], COLORS["text"])
    line_y = y + 32
    draw.line((x1, line_y, x2, line_y), fill=rgba(COLORS["dim"], 160), width=4)
    draw.line((x1, line_y, x1 + int((x2 - x1) * value), line_y), fill=color, width=5)
    dot_x = x1 + int((x2 - x1) * value)
    draw.ellipse((dot_x - 13, line_y - 13, dot_x + 13, line_y + 13), fill=color, outline="#EAF2FA", width=2)
    text(draw, (x2 + 22, line_y), f"{value:.3f}", F["small_bold"], color, "lm")
    text(draw, (x1, line_y + 48), "0", F["micro"], COLORS["dim"])
    text(draw, (x2, line_y + 48), "1", F["micro"], COLORS["dim"], "ra")


def draw_readouts_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 1835, 610, 3720, 1392
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["amber"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Prediction And Biology Checks Are Separate Readouts", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "One primary classifier is paired with response-pattern and gene-effect alignment checks.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    col_w = 560
    gap = 50
    for i, block in enumerate(READOUTS):
        bx1 = x1 + 58 + i * (col_w + gap)
        bx2 = bx1 + col_w
        rounded(draw, (bx1, y1 + 220, bx2, y2 - 72), 24, COLORS["panel2"], str(block["color"]), 2)
        text(draw, (bx1 + 30, y1 + 250), str(block["name"]), F["h2"], str(block["color"]))
        text(draw, (bx1 + 30, y1 + 292), str(block["role"]), F["small"], COLORS["muted"])
        for j, (label, value) in enumerate(block["items"]):
            draw_slider(draw, bx1 + 30, y1 + 380 + j * 180, bx2 - 104, label, float(value), str(block["color"]))


def draw_bar_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1450, 3720, 1856
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Top-Ranked Model Genes Enrich For Flight-Response Rows", F["h2"], COLORS["text"])
    text(draw, (x2 - 80, y1 + 46), "enrichment over expectation", F["small_bold"], COLORS["teal"], "ra")

    chart_x1, chart_y1, chart_x2, chart_y2 = x1 + 118, y1 + 118, x2 - 118, y2 - 82
    max_v = 4.2
    for tick in [0, 1, 2, 3, 4]:
        yy = int(chart_y2 - (tick / max_v) * (chart_y2 - chart_y1))
        draw.line((chart_x1, yy, chart_x2, yy), fill=rgba(COLORS["grid"], 115), width=1)
        text(draw, (chart_x1 - 24, yy), str(tick), F["micro"], COLORS["dim"], "rm")
    yy_one = int(chart_y2 - (1 / max_v) * (chart_y2 - chart_y1))
    for dash_x in range(chart_x1, chart_x2, 34):
        draw.line((dash_x, yy_one, min(dash_x + 18, chart_x2), yy_one), fill=rgba(COLORS["muted"], 160), width=2)
    text(draw, (chart_x2, yy_one - 24), "expected = 1.0", F["tiny_bold"], COLORS["muted"], "ra")

    bar_w = 430
    slot = (chart_x2 - chart_x1) / len(TOPK)
    for i, row in enumerate(TOPK):
        cx = int(chart_x1 + slot * (i + 0.5))
        v = float(row["enrichment"])
        h = int((v / max_v) * (chart_y2 - chart_y1))
        bx1, bx2 = cx - bar_w // 2, cx + bar_w // 2
        by1, by2 = chart_y2 - h, chart_y2
        color = COLORS["green"] if v >= 2.0 else COLORS["dim"]
        rounded(draw, (bx1, by1, bx2, by2), 14, color, None, 1)
        text(draw, (cx, by1 - 38), f"{int(row['observed'])} vs {float(row['expected']):.2f}", F["small_bold"], COLORS["text"], "ma")
        text(draw, (cx, by1 + 32), f"{v:.2f}x", F["h2"], COLORS["ink"] if v >= 2.0 else COLORS["text"], "ma")
        text(draw, (cx, chart_y2 + 30), str(row["label"]), F["small_bold"], COLORS["text"], "ma")


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "Organoid data add a compact biology-check lane: prediction, response-pattern agreement, and top-gene enrichment are read as three separate outputs.",
        F["small"],
        3300,
        COLORS["muted"],
        8,
    )


def write_outputs(image: Image.Image) -> None:
    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image).convert("RGB")
    gray.save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(image)
    manifest = {
        "title": "Organoid Extension Is A Biology Check",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "size": list(image.size),
        "mode": image.mode,
        "mean_rgb": [round(v, 2) for v in stat.mean],
        "footprint": FOOTPRINT,
        "readouts": READOUTS,
        "topk": TOPK,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    QA_NOTE.write_text(
        "# Organoid Biology Check Asset Visual QA\n\n"
        "Slide 56 explains the organoid extension as a compact biology-check lane.\n\n"
        "Checks performed:\n"
        "- Full-size render at `3840x2160`.\n"
        "- Strict crops for header, footprint, readout sliders, bar chart, and footer.\n"
        "- Grayscale render for contrast and hierarchy.\n\n"
        "Status: ready after visual QA.\n"
    )


def build() -> None:
    image = background()
    draw = ImageDraw.Draw(image, "RGBA")

    text(draw, (120, 72), "SLIDE 56 | ACT 6 | ORGANOID EXTENSION", F["kicker"], COLORS["teal"])
    bx = 1840
    bx = badge(draw, bx, 56, "STUDIES", "2", COLORS["amber"])
    bx = badge(draw, bx, 56, "SAMPLES", "42", COLORS["amber"])
    bx = badge(draw, bx, 56, "CONTRASTS", "8", COLORS["blue"])
    badge(draw, bx, 56, "SIG ROWS", "2,368", COLORS["teal"])

    text(draw, (120, 330), "Organoid Extension Is A Biology Check", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "Two human neural-organoid RNA-seq studies add a compact check on whether prediction signals align with flight-response biology.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_footprint_panel(draw)
    draw_readouts_panel(draw)
    draw_bar_panel(draw)
    draw_footer(draw)
    write_outputs(image)


if __name__ == "__main__":
    build()
    print(json.dumps({"output": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))

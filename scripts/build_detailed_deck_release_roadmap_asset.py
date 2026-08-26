#!/usr/bin/env python3
"""Build slide 59 asset: roadmap to release."""

from __future__ import annotations

import json
import math
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
OUT_DIR = ASSET_ROOT / "release_roadmap"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "roadmap_to_release_premium.png"
GRAY_PATH = OUT_DIR / "roadmap_to_release_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "release_roadmap_manifest.json"
QA_NOTE = OUT_DIR / "RELEASE_ROADMAP_ASSET_VISUAL_QA.md"

COLORS = {
    "bg": "#0B1119",
    "bg2": "#101722",
    "header": "#111B28",
    "panel": "#111B28",
    "panel2": "#162236",
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

BADGES = [
    ("READY", "58 assets", COLORS["teal"]),
    ("BULK TASKS", "8", COLORS["green"]),
    ("CHECKSUMS", "22 / 22", COLORS["blue"]),
    ("NEXT", "2 close slides", COLORS["violet"]),
]

LANES = [
    {
        "lane": "Data Package",
        "tag": "tables",
        "color": COLORS["teal"],
        "in_hand": ["descriptor + 21 entries", "18 CSV tables", "8 task files"],
        "next": ["local file mirror", "22 hash checks", "package version tag"],
        "release": ["versioned package", "portable tables", "readable data card"],
    },
    {
        "lane": "Evaluator Runs",
        "tag": "scores",
        "color": COLORS["green"],
        "in_hand": ["task registry", "metric profiles", "run record shape"],
        "next": ["canonical reruns", "metrics + run manifests", "report cross-check"],
        "release": ["reproducible score bundle", "comparable reports", "archived run rows"],
    },
    {
        "lane": "Extension Inputs",
        "tag": "tracks",
        "color": COLORS["amber"],
        "in_hand": ["organoid check", "OSD-120 interaction check", "single-cell metric contract"],
        "next": ["fixed single-cell matrix", "extension run records", "appendix table set"],
        "release": ["extension appendix", "input checklist", "scoring examples"],
    },
    {
        "lane": "Paper And Deck",
        "tag": "figures",
        "color": COLORS["violet"],
        "in_hand": ["58 detailed assets", "copy-tone QA", "grayscale checks"],
        "next": ["final close slide", "figure/table crosswalk", "methods appendix"],
        "release": ["paper-aligned package", "presentation appendix", "review-ready figures"],
    },
]

MILESTONES = [
    ("Metadata alpha", "indexed now", COLORS["teal"]),
    ("File mirror", "next object", COLORS["amber"]),
    ("Hash checks", "22 rows", COLORS["blue"]),
    ("Evaluator runs", "metrics.json", COLORS["green"]),
    ("Paper package", "figures + tables", COLORS["violet"]),
    ("Release", "versioned bundle", COLORS["rose"]),
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
    "title": load_font(84, True),
    "subtitle": load_font(37),
    "section": load_font(40, True),
    "h2": load_font(35, True),
    "h3": load_font(30, True),
    "body": load_font(26),
    "body_bold": load_font(26, True),
    "small": load_font(23),
    "small_bold": load_font(23, True),
    "tiny": load_font(20),
    "tiny_bold": load_font(20, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "metric": load_font(62, True),
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
        color = blend(COLORS["bg"], COLORS["bg2"], (y / H) * 0.72)
        draw.line((0, y, W, y), fill=color)
    for x in range(0, W, 160):
        draw.line((x, 260, x, H - 220), fill=rgba(COLORS["grid"], 72), width=1)
    for y in range(320, H - 220, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 64), width=1)
    draw.rectangle((0, 0, W, 260), fill=rgba(COLORS["header"], 245))
    draw.rectangle((0, H - 190, W, H), fill=rgba("#071019", 235))
    return image


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    w = max(230, int(draw.textlength(value, font=F["small_bold"]) + 74))
    rounded(draw, (x, y, x + w, y + 72), 14, COLORS["panel2"], color, 2)
    text(draw, (x + 18, y + 14), label.upper(), F["micro_bold"], COLORS["muted"])
    text(draw, (x + 18, y + 42), value, F["small_bold"], COLORS["text"])
    return x + w + 24


def arrow(draw: ImageDraw.ImageDraw, x1: int, y1: int, x2: int, y2: int, color: str, width: int = 4) -> None:
    draw.line((x1, y1, x2, y2), fill=color, width=width)
    ang = math.atan2(y2 - y1, x2 - x1)
    size = 18
    p1 = (x2, y2)
    p2 = (x2 - size * math.cos(ang - 0.45), y2 - size * math.sin(ang - 0.45))
    p3 = (x2 - size * math.cos(ang + 0.45), y2 - size * math.sin(ang + 0.45))
    draw.polygon([p1, p2, p3], fill=color)


def draw_lane_cell(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], items: list[str], color: str) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 18, COLORS["panel2"], color, 2)
    for i, item in enumerate(items):
        y = y1 + 22 + i * 40
        draw.ellipse((x1 + 22, y + 4, x1 + 34, y + 16), fill=color)
        text(draw, (x1 + 48, y), item, F["tiny"], COLORS["muted"])


def draw_roadmap_table(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 610, 3720, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["teal"], width=5)

    col_x = [x1 + 360, x1 + 1355, x1 + 2350]
    col_w = 910
    headers = [
        ("Already in hand", COLORS["teal"]),
        ("Next fixed object", COLORS["amber"]),
        ("Release-ready output", COLORS["green"]),
    ]
    for (label, color), cx in zip(headers, col_x):
        rounded(draw, (cx, y1 + 58, cx + col_w, y1 + 112), 16, COLORS["panel3"], color, 2)
        text(draw, (cx + 28, y1 + 74), label, F["small_bold"], color)

    lane_top = y1 + 150
    lane_h = 154
    for i, lane in enumerate(LANES):
        y = lane_top + i * 168
        color = str(lane["color"])
        rounded(draw, (x1 + 44, y, x1 + 318, y + lane_h), 22, COLORS["panel3"], color, 2)
        text(draw, (x1 + 72, y + 34), str(lane["lane"]), F["h3"], color)
        text(draw, (x1 + 72, y + 78), str(lane["tag"]), F["tiny_bold"], COLORS["muted"])
        draw_lane_cell(draw, (col_x[0], y, col_x[0] + col_w, y + lane_h), list(lane["in_hand"]), color)
        draw_lane_cell(draw, (col_x[1], y, col_x[1] + col_w, y + lane_h), list(lane["next"]), color)
        draw_lane_cell(draw, (col_x[2], y, col_x[2] + col_w, y + lane_h), list(lane["release"]), color)
        arrow(draw, col_x[0] + col_w + 30, y + lane_h // 2, col_x[1] - 40, y + lane_h // 2, COLORS["dim"], 3)
        arrow(draw, col_x[1] + col_w + 30, y + lane_h // 2, col_x[2] - 40, y + lane_h // 2, COLORS["dim"], 3)


def draw_milestone_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1530, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Release Runway", F["h2"], COLORS["text"])

    start_x, end_x, line_y = x1 + 650, x2 - 250, y1 + 172
    draw.line((start_x, line_y, end_x, line_y), fill=rgba(COLORS["dim"], 170), width=6)
    gap = (end_x - start_x) / (len(MILESTONES) - 1)
    for i, (label, detail, color) in enumerate(MILESTONES):
        cx = int(start_x + i * gap)
        draw.ellipse((cx - 28, line_y - 28, cx + 28, line_y + 28), fill=color, outline="#EAF2FA", width=2)
        rounded(draw, (cx - 165, line_y - 118, cx + 165, line_y - 58), 16, COLORS["panel2"], color, 2)
        text(draw, (cx, line_y - 102), label, F["tiny_bold"], COLORS["text"], "ma")
        rounded(draw, (cx - 150, line_y + 56, cx + 150, line_y + 108), 14, COLORS["panel3"], color, 2)
        text(draw, (cx, line_y + 72), detail, F["micro_bold"], COLORS["muted"], "ma")


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "The release path is now a packaging sequence: fix files, verify hashes, run evaluators, and align paper figures to the same objects.",
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
        "title": "Roadmap To Release",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "size": list(image.size),
        "mode": image.mode,
        "mean_rgb": [round(v, 2) for v in stat.mean],
        "badges": BADGES,
        "lanes": LANES,
        "milestones": MILESTONES,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    QA_NOTE.write_text(
        "# Release Roadmap Asset Visual QA\n\n"
        "Slide 59 converts the platform and extension work into four release lanes.\n\n"
        "Checks performed:\n"
        "- Full-size render at `3840x2160`.\n"
        "- Strict crops for title, roadmap lanes, arrows, milestones, and footer.\n"
        "- Grayscale render for contrast and hierarchy.\n\n"
        "Status: ready after visual QA.\n"
    )


def build() -> None:
    image = background()
    draw = ImageDraw.Draw(image, "RGBA")

    text(draw, (120, 72), "SLIDE 59 | ACT 6 | RELEASE ROADMAP", F["kicker"], COLORS["teal"])
    bx = 1760
    for label, value, color in BADGES:
        bx = badge(draw, bx, 56, label, value, color)

    text(draw, (120, 330), "Roadmap To Release", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "Task contracts, metric profiles, metadata alpha, and extension checks are in place; release work now becomes a sequence of fixed files, evaluator runs, and paper alignment.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_roadmap_table(draw)
    draw_milestone_panel(draw)
    draw_footer(draw)
    write_outputs(image)


if __name__ == "__main__":
    build()
    print(json.dumps({"output": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))

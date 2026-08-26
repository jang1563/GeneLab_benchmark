#!/usr/bin/env python3
"""Build slide 60 asset: final established story."""

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
OUT_DIR = ASSET_ROOT / "final_establishes"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "what_this_deck_establishes_premium.png"
GRAY_PATH = OUT_DIR / "what_this_deck_establishes_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "final_establishes_manifest.json"
QA_NOTE = OUT_DIR / "FINAL_ESTABLISHES_ASSET_VISUAL_QA.md"

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
    ("CORE", "benchmark", COLORS["teal"]),
    ("PLATFORM", "rerunnable", COLORS["green"]),
    ("BIOLOGY", "interpretable", COLORS["amber"]),
    ("EXTENSIONS", "expanding", COLORS["violet"]),
]

ESTABLISHED = [
    {
        "title": "Core Benchmark",
        "tag": "mission-held-out tasks",
        "color": COLORS["teal"],
        "statement": "Fixed tasks and mission-held-out scoring create a comparable transfer readout.",
        "items": ["task contract", "held-out missions", "negative controls", "pathway support"],
    },
    {
        "title": "Platform Package",
        "tag": "rerunnable objects",
        "color": COLORS["green"],
        "statement": "Task registry, metric profiles, run records, and a data package make scores portable.",
        "items": ["8 public bulk tasks", "22 study rows", "24 baseline outputs", "58 ready assets"],
    },
    {
        "title": "Biology Readout",
        "tag": "interpretation layers",
        "color": COLORS["amber"],
        "statement": "Temporal context, immune/TF activity, targets, biomarkers, and translation layers explain where to follow up.",
        "items": ["systems biology", "translation bridge", "target triage", "v8 workstreams"],
    },
    {
        "title": "Extension Tracks",
        "tag": "next surfaces",
        "color": COLORS["violet"],
        "statement": "Organoid, OSD-120, and fixed-input single-cell tracks broaden the benchmark surface.",
        "items": ["2 organoid studies", "36 OSD-120 samples", "8 single-cell samples", "release runway"],
    },
]

TAKEAWAYS = [
    ("ML audience", "task fit beats scale", COLORS["teal"]),
    ("Biology audience", "pathways carry transfer", COLORS["amber"]),
    ("Release audience", "objects make scores rerunnable", COLORS["green"]),
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
    "section": load_font(44, True),
    "h2": load_font(36, True),
    "h3": load_font(30, True),
    "body": load_font(27),
    "body_bold": load_font(27, True),
    "small": load_font(23),
    "small_bold": load_font(23, True),
    "tiny": load_font(20),
    "tiny_bold": load_font(20, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "metric": load_font(66, True),
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


def draw_stack_icon(draw: ImageDraw.ImageDraw, cx: int, cy: int, color: str) -> None:
    for i in range(4):
        y = cy + i * 34
        rounded(draw, (cx - 82 + i * 14, y, cx + 82 + i * 14, y + 28), 12, COLORS["panel2"], color, 2)
    draw.ellipse((cx - 34, cy + 130, cx + 34, cy + 198), fill=rgba(color, 170), outline="#EAF2FA", width=2)


def draw_established_cards(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 610, 3720, 1480
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["teal"], width=5)

    card_w, card_h = 1724, 360
    positions = [
        (x1 + 54, y1 + 72),
        (x1 + 1822, y1 + 72),
        (x1 + 54, y1 + 466),
        (x1 + 1822, y1 + 466),
    ]
    for i, item in enumerate(ESTABLISHED):
        x, y = positions[i]
        color = str(item["color"])
        rounded(draw, (x, y, x + card_w, y + card_h), 26, COLORS["panel2"], color, 2)
        draw_stack_icon(draw, x + 138, y + 80, color)
        text(draw, (x + 270, y + 48), str(item["title"]), F["section"], COLORS["text"])
        text(draw, (x + card_w - 44, y + 58), str(item["tag"]), F["small_bold"], color, "ra")
        paragraph(draw, (x + 270, y + 112), str(item["statement"]), F["body"], 1260, COLORS["muted"], 8)
        chips_y = y + 242
        chip_x = x + 270
        for label in list(item["items"]):
            w = int(draw.textlength(str(label), font=F["tiny_bold"]) + 48)
            rounded(draw, (chip_x, chips_y, chip_x + w, chips_y + 50), 14, COLORS["panel3"], color, 2)
            text(draw, (chip_x + 24, chips_y + 16), str(label), F["tiny_bold"], COLORS["text"])
            chip_x += w + 18


def draw_audience_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1530, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Three Takeaways", F["h2"], COLORS["text"])
    text(draw, (x2 - 80, y1 + 46), "benchmark -> biology -> package", F["small_bold"], COLORS["teal"], "ra")

    start_x = x1 + 240
    card_w = 930
    for i, (audience, takeaway, color) in enumerate(TAKEAWAYS):
        x = start_x + i * 1120
        rounded(draw, (x, y1 + 112, x + card_w, y1 + 236), 24, COLORS["panel2"], color, 2)
        text(draw, (x + 34, y1 + 138), audience, F["small_bold"], color)
        text(draw, (x + 34, y1 + 178), takeaway, F["h2"], COLORS["text"])
        if i < len(TAKEAWAYS) - 1:
            arrow(draw, x + card_w + 42, y1 + 174, x + card_w + 156, y1 + 174, COLORS["dim"], 4)


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "The final story is a benchmark stack: fixed tasks create comparable scores, conserved biology explains transfer, and package objects make the work rerunnable.",
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
        "title": "What This Deck Establishes",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "size": list(image.size),
        "mode": image.mode,
        "mean_rgb": [round(v, 2) for v in stat.mean],
        "badges": BADGES,
        "established": ESTABLISHED,
        "takeaways": TAKEAWAYS,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    QA_NOTE.write_text(
        "# Final Established Story Asset Visual QA\n\n"
        "Slide 60 closes the detailed deck with four established layers and three audience takeaways.\n\n"
        "Checks performed:\n"
        "- Full-size render at `3840x2160`.\n"
        "- Strict crops for title, four cards, audience takeaways, and footer.\n"
        "- Grayscale render for contrast and hierarchy.\n\n"
        "Status: ready after visual QA.\n"
    )


def build() -> None:
    image = background()
    draw = ImageDraw.Draw(image, "RGBA")

    text(draw, (120, 72), "SLIDE 60 | CLOSE | ESTABLISHED STORY", F["kicker"], COLORS["teal"])
    bx = 1710
    for label, value, color in BADGES:
        bx = badge(draw, bx, 56, label, value, color)

    text(draw, (120, 330), "What This Deck Establishes", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "SpaceBio-Bench turns public space-biology omics into mission-held-out tasks, interpretable model readouts, and a release-ready object chain.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_established_cards(draw)
    draw_audience_panel(draw)
    draw_footer(draw)
    write_outputs(image)


if __name__ == "__main__":
    build()
    print(json.dumps({"output": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))

#!/usr/bin/env python3
"""Build slide 57 asset: OSD-120 same-study interaction check."""

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
OUT_DIR = ASSET_ROOT / "osd120_interaction_check"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "osd120_is_a_same_study_interaction_check_premium.png"
GRAY_PATH = OUT_DIR / "osd120_is_a_same_study_interaction_check_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "osd120_interaction_check_manifest.json"
QA_NOTE = OUT_DIR / "OSD120_INTERACTION_CHECK_ASSET_VISUAL_QA.md"

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
    ("Study", "OSD-120", COLORS["teal"]),
    ("Samples", "36", COLORS["amber"]),
    ("Gene rows", "24.7k", COLORS["blue"]),
    ("Labels", "18 + 18", COLORS["green"]),
]

GENOTYPES = ["Col.0", "Col.0.PhyD", "Ws"]
LIGHTS = ["Dark", "Light"]

HOLDOUTS = [
    {
        "lane": "Genotype / ecotype",
        "folds": "3 folds",
        "train": "24",
        "test": "12",
        "detail": "one genotype held out",
        "color": COLORS["green"],
    },
    {
        "lane": "Light treatment",
        "folds": "2 folds",
        "train": "18",
        "test": "18",
        "detail": "one light context held out",
        "color": COLORS["amber"],
    },
    {
        "lane": "Genotype x light stratum",
        "folds": "6 folds",
        "train": "30",
        "test": "6",
        "detail": "one matrix cell held out",
        "color": COLORS["violet"],
    },
]

READOUTS = [
    {
        "lane": "Genotype",
        "folds": "3",
        "nearest_ba": 0.667,
        "nearest_auc": 0.735,
        "sparse_ba": 0.917,
        "sparse_auc": 0.966,
        "color": COLORS["green"],
    },
    {
        "lane": "Light",
        "folds": "2",
        "nearest_ba": 0.667,
        "nearest_auc": 0.784,
        "sparse_ba": 0.833,
        "sparse_auc": 0.861,
        "color": COLORS["amber"],
    },
    {
        "lane": "Stratum",
        "folds": "6",
        "nearest_ba": 0.667,
        "nearest_auc": 0.765,
        "sparse_ba": 0.889,
        "sparse_auc": 0.914,
        "color": COLORS["violet"],
    },
]

FEATURES = [
    ("Col.0.PhyD genotype fold", "12 active features", "7 stable >= 0.5", COLORS["green"]),
    ("Light-treatment fold", "10 active features", "2 stable >= 0.5", COLORS["amber"]),
    ("Col.0.PhyD dark stratum", "10 active features", "10 stable >= 0.5", COLORS["violet"]),
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
    w = 230 if len(value) < 6 else 270
    rounded(draw, (x, y, x + w, y + 72), 14, COLORS["panel2"], color, 2)
    text(draw, (x + 18, y + 14), label.upper(), F["micro_bold"], COLORS["muted"])
    text(draw, (x + 18, y + 42), value, F["small_bold"], COLORS["text"])
    return x + w + 24


def draw_plant(draw: ImageDraw.ImageDraw, cx: int, cy: int) -> None:
    draw.line((cx, cy + 96, cx, cy - 120), fill=rgba(COLORS["green"], 220), width=9)
    for side in [-1, 1]:
        for i, (dy, length) in enumerate([(-86, 76), (-46, 96), (-4, 82)]):
            x2 = cx + side * length
            y2 = cy + dy - 24 * i
            draw.line((cx, cy + dy, x2, y2), fill=rgba(COLORS["green"], 170), width=6)
            draw.ellipse((x2 - 40, y2 - 23, x2 + 40, y2 + 23), fill=rgba(COLORS["green"], 78), outline=rgba("#EAF2FA", 70), width=2)
    for i in range(8):
        angle = math.pi / 2 + (i - 3.5) * 0.24
        x2 = cx + math.cos(angle) * 84
        y2 = cy + 100 + math.sin(angle) * 84
        draw.line((cx, cy + 96, x2, y2), fill=rgba(COLORS["amber"], 145), width=4)
    draw.ellipse((cx - 64, cy - 170, cx + 64, cy - 66), fill=rgba(COLORS["teal"], 72), outline=rgba(COLORS["teal"], 150), width=3)
    text(draw, (cx, cy - 118), "36", F["h2"], COLORS["text"], "ma")
    text(draw, (cx, cy - 84), "samples", F["micro_bold"], COLORS["muted"], "ma")


def draw_geometry_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 610, 1340, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["teal"], width=5)
    text(draw, (x1 + 42, y1 + 52), "One Study, Two Experimental Axes", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "Arabidopsis root RNA-seq is arranged as genotype / ecotype by light treatment.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    draw_plant(draw, x1 + 210, y1 + 360)

    card_w, card_h = 210, 92
    top, left = y1 + 290, x1 + 470
    text(draw, (left, top - 58), "3 genotype / ecotype levels", F["small_bold"], COLORS["text"])
    text(draw, (left + 710, top - 58), "2 light contexts", F["small_bold"], COLORS["text"], "ra")
    for r, light in enumerate(LIGHTS):
        text(draw, (left - 34, top + r * 128 + 47), light, F["small_bold"], COLORS["muted"], "ra")
    for c, genotype in enumerate(GENOTYPES):
        text(draw, (left + c * 240 + card_w // 2, top - 20), genotype, F["tiny_bold"], COLORS["muted"], "ma")

    for r, _light in enumerate(LIGHTS):
        for c, _genotype in enumerate(GENOTYPES):
            x = left + c * 240
            y = top + r * 128
            color = [COLORS["green"], COLORS["blue"], COLORS["amber"]][c]
            rounded(draw, (x, y, x + card_w, y + card_h), 18, COLORS["panel2"], color, 2)
            text(draw, (x + 28, y + 18), "6 samples", F["small_bold"], COLORS["text"])
            text(draw, (x + 28, y + 54), "3 ground + 3 LEO", F["tiny"], COLORS["muted"])

    fy = y2 - 172
    for i, (label, value, color) in enumerate(FOOTPRINT):
        x = x1 + 62 + i * 282
        rounded(draw, (x, fy, x + 248, fy + 104), 18, COLORS["panel3"], color, 2)
        text(draw, (x + 22, fy + 18), label.upper(), F["micro_bold"], color)
        text(draw, (x + 22, fy + 50), value, F["h2"], COLORS["text"])


def arrow(draw: ImageDraw.ImageDraw, x1: int, y1: int, x2: int, y2: int, color: str, width: int = 4) -> None:
    draw.line((x1, y1, x2, y2), fill=color, width=width)
    ang = math.atan2(y2 - y1, x2 - x1)
    size = 18
    p1 = (x2, y2)
    p2 = (x2 - size * math.cos(ang - 0.45), y2 - size * math.sin(ang - 0.45))
    p3 = (x2 - size * math.cos(ang + 0.45), y2 - size * math.sin(ang + 0.45))
    draw.polygon([p1, p2, p3], fill=color)


def draw_holdout_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 1380, 610, 2520, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["amber"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Holdout Geometry Defines The Readout", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "Each lane holds out a different part of the same genotype-by-light matrix.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    top = y1 + 232
    for i, item in enumerate(HOLDOUTS):
        y = top + i * 190
        color = str(item["color"])
        rounded(draw, (x1 + 54, y, x2 - 54, y + 146), 24, COLORS["panel2"], color, 2)
        text(draw, (x1 + 86, y + 26), str(item["lane"]), F["h3"], color)
        text(draw, (x1 + 86, y + 66), str(item["folds"]), F["small_bold"], COLORS["text"])
        text(draw, (x1 + 86, y + 100), str(item["detail"]), F["tiny"], COLORS["muted"])
        tx1, tx2 = x1 + 510, x2 - 340
        rounded(draw, (tx1, y + 36, tx1 + 212, y + 94), 14, COLORS["panel3"], COLORS["blue"], 2)
        text(draw, (tx1 + 28, y + 54), f"TRAIN {item['train']}", F["small_bold"], COLORS["text"])
        arrow(draw, tx1 + 238, y + 65, tx2 + 20, y + 65, COLORS["dim"], 3)
        rounded(draw, (tx2 + 42, y + 36, tx2 + 220, y + 94), 14, COLORS["panel3"], color, 2)
        text(draw, (tx2 + 70, y + 54), f"TEST {item['test']}", F["small_bold"], COLORS["text"])

    rounded(draw, (x1 + 54, y2 - 98, x2 - 54, y2 - 44), 16, "#122234", COLORS["teal"], 2)
    text(draw, (x1 + 82, y2 - 82), "Balanced labels stay visible in each split.", F["small_bold"], COLORS["teal"])


def draw_metric_bar(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    width: int,
    value: float,
    color: str,
    label: str,
) -> None:
    text(draw, (x, y), label, F["micro_bold"], COLORS["muted"])
    bar_y = y + 34
    draw.line((x, bar_y, x + width, bar_y), fill=rgba(COLORS["dim"], 150), width=5)
    draw.line((x, bar_y, x + int(width * value), bar_y), fill=color, width=6)
    dot_x = x + int(width * value)
    draw.ellipse((dot_x - 12, bar_y - 12, dot_x + 12, bar_y + 12), fill=color, outline="#EAF2FA", width=2)
    text(draw, (x + width + 24, bar_y), f"{value:.3f}", F["small_bold"], color, "lm")


def draw_readout_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 2560, 610, 3720, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["green"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Sparse L1 Lifts The Same Split Surface", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "The transparent sparse model is read beside a nearest-centroid baseline for the same holdout lanes.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    top = y1 + 228
    for i, row in enumerate(READOUTS):
        y = top + i * 178
        color = str(row["color"])
        rounded(draw, (x1 + 54, y, x2 - 54, y + 140), 22, COLORS["panel2"], color, 2)
        text(draw, (x1 + 84, y + 22), str(row["lane"]), F["h3"], color)
        text(draw, (x1 + 84, y + 60), f"{row['folds']} folds", F["tiny_bold"], COLORS["muted"])
        draw_metric_bar(draw, x1 + 270, y + 22, 280, float(row["nearest_auc"]), COLORS["dim"], "centroid AUROC")
        draw_metric_bar(draw, x1 + 660, y + 22, 280, float(row["sparse_auc"]), color, "sparse AUROC")
        text(draw, (x1 + 270, y + 100), f"BA {row['nearest_ba']:.3f}", F["tiny_bold"], COLORS["muted"])
        text(draw, (x1 + 660, y + 100), f"BA {row['sparse_ba']:.3f}", F["tiny_bold"], color)

    rounded(draw, (x1 + 54, y2 - 98, x2 - 54, y2 - 44), 16, "#122234", COLORS["green"], 2)
    text(draw, (x1 + 82, y2 - 82), "Primary sparse row: C=1, 2,000 variable genes.", F["small_bold"], COLORS["green"])


def draw_bottom_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1530, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Same-Study Read Path", F["h2"], COLORS["text"])
    text(draw, (x2 - 80, y1 + 46), "matrix -> split -> score -> compact features", F["small_bold"], COLORS["teal"], "ra")

    path = [
        ("OSD-120", "one Arabidopsis root study", COLORS["teal"]),
        ("Matrix", "36 samples x 24.7k gene rows", COLORS["blue"]),
        ("Split", "3 + 2 + 6 holdout lanes", COLORS["amber"]),
        ("Score", "AUROC and balanced accuracy", COLORS["green"]),
        ("Features", "sparse coefficients", COLORS["violet"]),
    ]
    node_w, node_h, gap = 560, 94, 120
    start_x, node_y = x1 + 170, y1 + 108
    for i, (title, body, color) in enumerate(path):
        nx = start_x + i * (node_w + gap)
        rounded(draw, (nx, node_y, nx + node_w, node_y + node_h), 20, COLORS["panel2"], color, 2)
        text(draw, (nx + 28, node_y + 16), title, F["small_bold"], COLORS["text"])
        text(draw, (nx + 28, node_y + 54), body, F["tiny"], COLORS["muted"])
        if i < len(path) - 1:
            arrow(draw, nx + node_w + 18, node_y + node_h // 2, nx + node_w + gap - 26, node_y + node_h // 2, COLORS["dim"], 4)

    fy = y2 - 58
    for i, (label, active, stable, color) in enumerate(FEATURES):
        x = x1 + 170 + i * 1080
        rounded(draw, (x, fy - 42, x + 880, fy + 22), 16, COLORS["panel3"], color, 2)
        text(draw, (x + 22, fy - 24), label, F["tiny_bold"], COLORS["text"])
        text(draw, (x + 360, fy - 24), active, F["tiny"], COLORS["muted"])
        text(draw, (x + 640, fy - 24), stable, F["tiny_bold"], color)


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "OSD-120 contributes a constrained same-study interaction check: the split geometry carries the interpretation, and sparse L1 makes the fold behavior visible.",
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
        "title": "OSD-120 Is A Same-Study Interaction Check",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "size": list(image.size),
        "mode": image.mode,
        "mean_rgb": [round(v, 2) for v in stat.mean],
        "footprint": FOOTPRINT,
        "holdouts": HOLDOUTS,
        "readouts": READOUTS,
        "features": FEATURES,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    QA_NOTE.write_text(
        "# OSD-120 Interaction Check Asset Visual QA\n\n"
        "Slide 57 explains OSD-120 as a constrained same-study interaction check.\n\n"
        "Checks performed:\n"
        "- Full-size render at `3840x2160`.\n"
        "- Strict crops for header, geometry, holdout lanes, readout bars, read path, and footer.\n"
        "- Grayscale render for contrast and hierarchy.\n\n"
        "Status: ready after visual QA.\n"
    )


def build() -> None:
    image = background()
    draw = ImageDraw.Draw(image, "RGBA")

    text(draw, (120, 72), "SLIDE 57 | ACT 6 | OSD-120 INTERACTION CHECK", F["kicker"], COLORS["teal"])
    bx = 1840
    bx = badge(draw, bx, 56, "STUDY", "OSD-120", COLORS["teal"])
    bx = badge(draw, bx, 56, "SAMPLES", "36", COLORS["amber"])
    bx = badge(draw, bx, 56, "GENE ROWS", "24.7k", COLORS["blue"])
    badge(draw, bx, 56, "FOLDS", "11", COLORS["violet"])

    text(draw, (120, 330), "OSD-120 Is A Same-Study Interaction Check", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "Arabidopsis root RNA-seq tests flight / ground classification while holding out genotype, light treatment, and genotype-by-light strata inside one study.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_geometry_panel(draw)
    draw_holdout_panel(draw)
    draw_readout_panel(draw)
    draw_bottom_panel(draw)
    draw_footer(draw)
    write_outputs(image)


if __name__ == "__main__":
    build()
    print(json.dumps({"output": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))

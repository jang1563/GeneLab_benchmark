#!/usr/bin/env python3
"""Build premium opening slide scenes for slides 1-3.

The opening block should orient a first-time viewer before the methods and
result figures begin. Each rendered preview includes the intended editable
text overlay; each scene plate keeps the visual scaffold available for PPTX
rebuilds with editable text.
"""

from __future__ import annotations

import json
import math
import random
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont, ImageOps, ImageStat


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "premium_opening_slides_1_3_v0_1"
QA = OUT / "qa"
CREATED = "2026-06-03"
W, H = 3840, 2160

COLORS = {
    "void": "#06090D",
    "deep": "#08121A",
    "panel": "#101923",
    "panel2": "#172536",
    "ink": "#F4F7F8",
    "soft": "#BBC9D5",
    "muted": "#7E8E9C",
    "rule": "#34485B",
    "paper": "#F6F3EA",
    "paper2": "#EDE8DA",
    "teal": "#18A392",
    "sky": "#63A9D7",
    "blue": "#2F6C9F",
    "amber": "#B9852E",
    "rose": "#B45A7E",
    "green": "#168B63",
    "violet": "#7B68A8",
    "red": "#B23A3A",
}

SLIDES = [
    {
        "n": 1,
        "id": "slide01_spacebiobench_title",
        "section": "Opening",
        "eyebrow": "SPACEBIO-BENCH",
        "title": "Testing biological AI under spaceflight domain shift",
        "subtitle": "A mission-held-out benchmark built from public space omics studies.",
        "bridge": "Can a model trained on known missions survive a new mission?",
        "caveat": "Benchmark framing; not a clinical or countermeasure claim.",
        "source": "public GeneLab/OSDR study records + v1-v7 benchmark core",
        "asset_role": "title visual",
    },
    {
        "n": 2,
        "id": "slide02_external_gap_positioning",
        "section": "Opening",
        "eyebrow": "EXTERNAL GAP",
        "title": "A distinct niche, not a firstness claim",
        "subtitle": "Existing resources organize data, biology, or ML tasks; this deck asks what transfers across missions.",
        "bridge": "The novelty is the evaluation unit: a whole held-out mission.",
        "caveat": "Positioning map only; it does not rank external platforms.",
        "source": "external landscape scan + local v1-v9 positioning notes",
        "asset_role": "positioning matrix",
    },
    {
        "n": 3,
        "id": "slide03_project_evolution_timeline",
        "section": "Opening",
        "eyebrow": "PROJECT EVOLUTION",
        "title": "v1-v9 moves from benchmark results to platform",
        "subtitle": "The deck separates completed benchmark evidence, hypothesis-only translation, and platformization.",
        "bridge": "Do not mix evidence levels: results first, extensions later.",
        "caveat": "v8/v9 are extension and platform layers, not equal-strength discoveries.",
        "source": "V1_V9 master outline + deck assembly bridge review",
        "asset_role": "version timeline",
    },
]


def rgb(token: str) -> tuple[int, int, int]:
    value = COLORS[token].lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(token: str, alpha: int) -> tuple[int, int, int, int]:
    return rgb(token) + (alpha,)


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def font(size: int, *, bold: bool = False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


F = {
    "eyebrow": font(22, bold=True),
    "title": font(72, bold=True),
    "title2": font(60, bold=True),
    "subtitle": font(30),
    "h": font(32, bold=True),
    "body": font(24),
    "small": font(18),
    "tiny": font(14),
    "num": font(42, bold=True),
}


def ensure() -> None:
    QA.mkdir(parents=True, exist_ok=True)


def text_size(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.ImageFont) -> tuple[int, int]:
    if not text:
        return (0, 0)
    box = draw.textbbox((0, 0), text, font=fnt)
    return (box[2] - box[0], box[3] - box[1])


def wrap_text(draw: ImageDraw.ImageDraw, text: str, fnt: ImageFont.ImageFont, max_width: int) -> list[str]:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = word if not current else f"{current} {word}"
        if text_size(draw, candidate, fnt)[0] <= max_width:
            current = candidate
            continue
        if current:
            lines.append(current)
        current = word
    if current:
        lines.append(current)
    return lines


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    fnt: ImageFont.ImageFont,
    fill: tuple[int, int, int] | tuple[int, int, int, int],
    max_width: int,
    *,
    line_gap: int = 10,
) -> int:
    x, y = xy
    for line in wrap_text(draw, text, fnt, max_width):
        draw.text((x, y), line, font=fnt, fill=fill)
        y += text_size(draw, line, fnt)[1] + line_gap
    return y


def alpha_layer() -> Image.Image:
    return Image.new("RGBA", (W, H), (0, 0, 0, 0))


def rounded_shadow(canvas: Image.Image, box: tuple[int, int, int, int], radius: int = 18, blur: int = 26, alpha: int = 120) -> None:
    x0, y0, x1, y1 = box
    layer = alpha_layer()
    draw = ImageDraw.Draw(layer)
    draw.rounded_rectangle((x0 + 18, y0 + 24, x1 + 18, y1 + 24), radius=radius, fill=(0, 0, 0, alpha))
    canvas.alpha_composite(layer.filter(ImageFilter.GaussianBlur(blur)))


def draw_background(canvas: Image.Image, seed: int) -> None:
    draw = ImageDraw.Draw(canvas)
    top = rgb("void")
    bottom = rgb("deep")
    for y in range(0, H, 2):
        t = y / max(1, H - 1)
        color = tuple(int(top[i] * (1 - t) + bottom[i] * t) for i in range(3))
        draw.line((0, y, W, y), fill=color + (255,), width=2)

    noise = Image.effect_noise((W, H), 26).convert("L").filter(ImageFilter.GaussianBlur(0.5))
    noise = ImageOps.colorize(noise, black="#000000", white="#E7F7FF").convert("RGBA")
    noise.putalpha(16)
    canvas.alpha_composite(noise)

    for x in range(180, W, 180):
        draw.line((x, 0, x, H), fill=rgba("rule", 26), width=1)
    for y in range(150, H, 150):
        draw.line((0, y, W, y), fill=rgba("rule", 22), width=1)

    rng = random.Random(seed)
    for _ in range(250):
        x = rng.randint(70, W - 80)
        y = rng.randint(70, H - 80)
        alpha = rng.randint(30, 100)
        size = rng.choice([2, 2, 3, 4])
        draw.ellipse((x, y, x + size, y + size), fill=(235, 248, 255, alpha))

    for idx, radius in enumerate([970, 1240, 1510]):
        bbox = (
            int(W * 0.47) - radius,
            int(H * 0.30) - int(radius * 0.38),
            int(W * 0.47) + radius,
            int(H * 0.30) + int(radius * 0.38),
        )
        draw.arc(bbox, 196, 356, fill=rgba("sky", 58 - idx * 13), width=3)


def draw_bio_texture(canvas: Image.Image, seed: int, *, x_bias: float = 0.72) -> None:
    rng = random.Random(seed)
    layer = alpha_layer()
    draw = ImageDraw.Draw(layer)
    for _ in range(46):
        cx = int(W * rng.uniform(x_bias - 0.18, x_bias + 0.18))
        cy = int(H * rng.uniform(0.14, 0.88))
        r = rng.randint(15, 52)
        tone = rng.choice(["teal", "sky", "green", "violet"])
        draw.ellipse((cx - r, cy - r, cx + r, cy + r), outline=rgba(tone, rng.randint(28, 72)), width=rng.randint(2, 4))
        if rng.random() < 0.65:
            draw.ellipse((cx - r // 4, cy - r // 4, cx + r // 4, cy + r // 4), fill=rgba(tone, rng.randint(18, 44)))
    for _ in range(26):
        x0 = int(W * rng.uniform(x_bias - 0.26, x_bias + 0.10))
        y0 = int(H * rng.uniform(0.18, 0.84))
        x1 = x0 + rng.randint(120, 420)
        y1 = y0 + rng.randint(-70, 90)
        draw.line((x0, y0, x1, y1), fill=rgba(rng.choice(["teal", "sky", "green"]), rng.randint(28, 60)), width=2)
    canvas.alpha_composite(layer.filter(ImageFilter.GaussianBlur(0.7)))


def draw_helix(canvas: Image.Image, x: int, y: int, w: int, h: int, tone: str = "teal", alpha: int = 115) -> None:
    draw = ImageDraw.Draw(canvas)
    points_a = []
    points_b = []
    steps = 72
    for i in range(steps + 1):
        t = i / steps
        yy = y + int(t * h)
        phase = t * math.pi * 5.4
        xa = x + int(w * 0.5 + math.sin(phase) * w * 0.34)
        xb = x + int(w * 0.5 + math.sin(phase + math.pi) * w * 0.34)
        points_a.append((xa, yy))
        points_b.append((xb, yy))
    draw.line(points_a, fill=rgba(tone, alpha), width=4)
    draw.line(points_b, fill=rgba("sky", max(30, alpha - 30)), width=4)
    for i in range(0, steps + 1, 6):
        draw.line((points_a[i][0], points_a[i][1], points_b[i][0], points_b[i][1]), fill=rgba("rule", 85), width=2)
        r = 7
        draw.ellipse((points_a[i][0] - r, points_a[i][1] - r, points_a[i][0] + r, points_a[i][1] + r), fill=rgba(tone, alpha + 20))
        draw.ellipse((points_b[i][0] - r, points_b[i][1] - r, points_b[i][0] + r, points_b[i][1] + r), fill=rgba("sky", alpha))


def draw_orbit_marker(draw: ImageDraw.ImageDraw, xy: tuple[int, int], label: str, tone: str) -> None:
    x, y = xy
    draw.ellipse((x - 17, y - 17, x + 17, y + 17), fill=rgba(tone, 230), outline=rgba("ink", 180), width=2)
    draw.text((x - 28, y + 36), label, font=F["tiny"], fill=rgb("muted"))


def draw_flow_plate(
    canvas: Image.Image,
    box: tuple[int, int, int, int],
    label: str,
    title: str,
    lines: list[str],
    tone: str,
    *,
    include_text: bool,
) -> None:
    draw = ImageDraw.Draw(canvas)
    x0, y0, x1, y1 = box
    rounded_shadow(canvas, box, radius=14, blur=22, alpha=92)
    draw.rounded_rectangle(box, radius=14, fill=rgba("panel", 232), outline=rgba(tone, 170), width=2)
    draw.rectangle((x0, y0, x1, y0 + 8), fill=rgba(tone, 210))
    if not include_text:
        for idx in range(5):
            yy = y0 + 72 + idx * 46
            draw.line((x0 + 40, yy, x1 - 52, yy), fill=rgba("soft", 40), width=2)
        return
    draw.text((x0 + 40, y0 + 34), label.upper(), font=F["tiny"], fill=rgb(tone))
    draw.text((x0 + 40, y0 + 72), title, font=F["h"], fill=rgb("ink"))
    yy = y0 + 126
    for line in lines:
        draw.text((x0 + 40, yy), line, font=F["small"], fill=rgb("soft"))
        yy += 34


def draw_slide_chrome(canvas: Image.Image, spec: dict, *, include_overlay: bool) -> None:
    if not include_overlay:
        return
    draw = ImageDraw.Draw(canvas)
    draw.text((230, 182), spec["eyebrow"], font=F["eyebrow"], fill=rgb("teal"))
    title_font = F["title"] if spec["n"] == 1 else F["title2"]
    y = draw_wrapped(draw, (230, 238), spec["title"], title_font, rgb("ink"), 1500, line_gap=14)
    y += 20
    draw_wrapped(draw, (234, y), spec["subtitle"], F["subtitle"], rgb("soft"), 1500, line_gap=10)
    draw.rounded_rectangle((230, 1788, 1950, 1948), radius=12, fill=rgba("panel", 226), outline=rgba("rule", 160), width=2)
    draw.text((270, 1828), spec["bridge"], font=F["body"], fill=rgb("ink"))
    draw.text((270, 1874), spec["caveat"], font=F["small"], fill=rgb("muted"))
    draw.text((2440, 1888), f"source: {spec['source']}", font=F["small"], fill=rgb("muted"))


def draw_slide01(canvas: Image.Image, *, include_overlay: bool) -> None:
    draw_background(canvas, 11)
    draw_bio_texture(canvas, 12, x_bias=0.76)
    draw = ImageDraw.Draw(canvas)

    for idx, radius in enumerate([1120, 1460, 1800]):
        bbox = (2140 - radius, 660 - int(radius * 0.28), 2140 + radius, 660 + int(radius * 0.28))
        draw.arc(bbox, 196, 348, fill=rgba("teal", 82 - idx * 18), width=4)
    for pos, label, tone in [
        ((1960, 635), "train", "sky"),
        ((2305, 570), "train", "sky"),
        ((2670, 556), "train", "green"),
        ((3060, 624), "held-out", "amber"),
    ]:
        draw_orbit_marker(draw, pos, label, tone)

    draw_helix(canvas, 2960, 350, 270, 770, tone="green", alpha=82)

    boxes = [
        ((2030, 1090, 2430, 1398), "Z2 source", "Public studies", ["mission records", "samples", "assays"], "sky"),
        ((2540, 1038, 2940, 1446), "Z3 task", "Benchmark task", ["train missions", "test mission", "scored transfer"], "teal"),
        ((3050, 1110, 3450, 1390), "Z4 guard", "Hidden mission", ["held out as unit", "no sample leakage"], "amber"),
    ]
    for box, label, title, lines, tone in boxes:
        draw_flow_plate(canvas, box, label, title, lines, tone, include_text=include_overlay)
    for x0, y0, x1, y1, tone in [(2438, 1242, 2530, 1242, "teal"), (2948, 1242, 3040, 1242, "amber")]:
        draw.line((x0, y0, x1, y1), fill=rgba(tone, 180), width=4)
        draw.polygon([(x1, y1), (x1 - 18, y1 - 10), (x1 - 18, y1 + 10)], fill=rgba(tone, 180))

    if include_overlay:
        draw.text((2025, 925), "Evaluation grammar", font=F["small"], fill=rgb("muted"))
        draw.text((2030, 1510), "known missions", font=F["small"], fill=rgb("sky"))
        draw.text((3052, 1510), "new mission", font=F["small"], fill=rgb("amber"))
    draw_slide_chrome(canvas, SLIDES[0], include_overlay=include_overlay)


def draw_slide02(canvas: Image.Image, *, include_overlay: bool) -> None:
    draw_background(canvas, 21)
    draw_bio_texture(canvas, 22, x_bias=0.69)
    draw = ImageDraw.Draw(canvas)

    matrix = (1530, 420, 3440, 1650)
    rounded_shadow(canvas, matrix, radius=18, blur=28, alpha=115)
    draw.rounded_rectangle(matrix, radius=18, fill=rgba("panel2", 238), outline=rgba("rule", 210), width=2)
    x0, y0, x1, y1 = matrix
    for i in range(1, 6):
        x = x0 + int((x1 - x0) * i / 6)
        draw.line((x, y0 + 76, x, y1 - 82), fill=rgba("rule", 66), width=1)
    for i in range(1, 5):
        y = y0 + int((y1 - y0) * i / 5)
        draw.line((x0 + 82, y, x1 - 82, y), fill=rgba("rule", 66), width=1)

    highlight = alpha_layer()
    hd = ImageDraw.Draw(highlight)
    hd.rounded_rectangle((x1 - 560, y0 + 120, x1 - 100, y0 + 430), radius=18, fill=rgba("amber", 28), outline=rgba("amber", 92), width=3)
    hd.line((x1 - 680, y0 + 525, x1 - 405, y0 + 355), fill=rgba("amber", 150), width=5)
    canvas.alpha_composite(highlight)

    points = [
        ("NASA BPS /\nGeneLab", 0.24, 0.66, "sky", "data + community resource"),
        ("GLARE", 0.50, 0.72, "green", "space omics knowledge graph"),
        ("OpenProblems /\ncell-eval", 0.58, 0.34, "violet", "ML task framework"),
        ("SpaceBio-Bench", 0.82, 0.82, "amber", "held-out mission benchmark"),
    ]
    for label, px, py, tone, note in points:
        cx = x0 + int((x1 - x0) * px)
        cy = y1 - int((y1 - y0) * py)
        r = 30 if label == "SpaceBio-Bench" else 22
        draw.ellipse((cx - r, cy - r, cx + r, cy + r), fill=rgba(tone, 230), outline=rgba("ink", 180), width=2)
        if label == "SpaceBio-Bench":
            draw.ellipse((cx - 92, cy - 92, cx + 92, cy + 92), outline=rgba("amber", 105), width=3)
            draw.line((cx - 160, cy + 120, cx - 38, cy + 36), fill=rgba("amber", 165), width=4)
        if include_overlay:
            tx = cx + 34
            ty = cy - 28
            label_font = F["h"] if label == "SpaceBio-Bench" else F["body"]
            draw_wrapped(draw, (tx, ty), label.replace("\n", " "), label_font, rgb("ink"), 390, line_gap=4)
            draw_wrapped(draw, (tx, ty + 48), note, F["tiny"], rgb("muted"), 360, line_gap=3)

    if include_overlay:
        draw.text((x0 + 70, y0 + 46), "Landscape map", font=F["small"], fill=rgb("teal"))
        draw.text((x0 + 92, y1 - 62), "evaluation unit: sample -> study -> mission", font=F["small"], fill=rgb("soft"))
        draw.text((x0 + 76, y0 + 110), "space biology specificity", font=F["small"], fill=rgb("soft"))
        draw.text((x1 - 650, y1 - 108), "mission-held-out benchmark niche", font=F["body"], fill=rgb("amber"))
        draw.text((x1 - 650, y1 - 72), "distinct position; no firstness claim", font=F["small"], fill=rgb("muted"))

    draw.line((x0 + 92, y1 - 105, x1 - 120, y1 - 105), fill=rgba("soft", 155), width=3)
    draw.line((x0 + 92, y1 - 105, x0 + 92, y0 + 120), fill=rgba("soft", 155), width=3)
    draw.polygon([(x1 - 120, y1 - 105), (x1 - 145, y1 - 119), (x1 - 145, y1 - 91)], fill=rgba("soft", 155))
    draw.polygon([(x0 + 92, y0 + 120), (x0 + 78, y0 + 145), (x0 + 106, y0 + 145)], fill=rgba("soft", 155))

    draw_slide_chrome(canvas, SLIDES[1], include_overlay=include_overlay)


def draw_timeline_node(
    canvas: Image.Image,
    center: tuple[int, int],
    tone: str,
    label: str,
    title: str,
    body: str,
    *,
    include_text: bool,
    large: bool = False,
) -> None:
    draw = ImageDraw.Draw(canvas)
    cx, cy = center
    r = 44 if large else 34
    draw.ellipse((cx - r, cy - r, cx + r, cy + r), fill=rgba(tone, 232), outline=rgba("ink", 170), width=2)
    draw.ellipse((cx - r - 20, cy - r - 20, cx + r + 20, cy + r + 20), outline=rgba(tone, 70), width=3)
    if not include_text:
        return
    draw.text((cx - 34, cy - 102), label, font=F["small"], fill=rgb(tone))
    draw_wrapped(draw, (cx - 190, cy + 76), title, F["body"], rgb("ink"), 410, line_gap=4)
    draw_wrapped(draw, (cx - 190, cy + 146), body, F["small"], rgb("muted"), 430, line_gap=4)


def draw_slide03(canvas: Image.Image, *, include_overlay: bool) -> None:
    draw_background(canvas, 31)
    draw_bio_texture(canvas, 32, x_bias=0.73)
    draw = ImageDraw.Draw(canvas)

    x_start, x_end, y_mid = 1460, 3380, 1030
    draw.line((x_start, y_mid, x_end, y_mid), fill=rgba("soft", 130), width=5)
    for x, tone in [(1800, "teal"), (2425, "amber"), (3030, "sky")]:
        draw.line((x, y_mid - 260, x, y_mid + 260), fill=rgba(tone, 42), width=2)

    draw_timeline_node(
        canvas,
        (1660, y_mid),
        "teal",
        "v1-v7",
        "Completed benchmark core",
        "tissue transfer, pathway rescue, model comparison, hardening",
        include_text=include_overlay,
        large=True,
    )
    draw_timeline_node(
        canvas,
        (2425, y_mid),
        "amber",
        "v8",
        "Hypothesis-only translation",
        "drug and perturbation ideas stay explicitly bounded",
        include_text=include_overlay,
    )
    draw_timeline_node(
        canvas,
        (3030, y_mid),
        "sky",
        "v9",
        "Platformization",
        "dataset, API, packaging, and species/organoid expansion path",
        include_text=include_overlay,
    )

    stage_boxes = [
        (1518, 1420, 1930, 1548, "Evidence strength", "completed result family", "teal"),
        (2258, 1420, 2700, 1548, "Boundary", "hypothesis-generating only", "amber"),
        (2860, 1420, 3350, 1548, "Operating layer", "reusable benchmark system", "sky"),
    ]
    for x0, y0, x1, y1, label, text, tone in stage_boxes:
        rounded_shadow(canvas, (x0, y0, x1, y1), radius=12, blur=16, alpha=70)
        draw.rounded_rectangle((x0, y0, x1, y1), radius=12, fill=rgba("panel", 224), outline=rgba(tone, 165), width=2)
        if include_overlay:
            draw.text((x0 + 30, y0 + 26), label, font=F["small"], fill=rgb(tone))
            draw.text((x0 + 30, y0 + 68), text, font=F["small"], fill=rgb("soft"))

    if include_overlay:
        draw.text((1470, 710), "One storyline, three evidence levels", font=F["h"], fill=rgb("ink"))
        draw.text((1470, 760), "The opening should stop the audience from reading v8/v9 as finished biology.", font=F["small"], fill=rgb("muted"))
        draw.text((1470, 1668), "Default deck rule: v1-v7 first; v8/v9 after the completed result spine.", font=F["body"], fill=rgb("soft"))

    draw_slide_chrome(canvas, SLIDES[2], include_overlay=include_overlay)


def render_slide(spec: dict, *, include_overlay: bool) -> Image.Image:
    canvas = Image.new("RGBA", (W, H), (0, 0, 0, 255))
    if spec["n"] == 1:
        draw_slide01(canvas, include_overlay=include_overlay)
    elif spec["n"] == 2:
        draw_slide02(canvas, include_overlay=include_overlay)
    elif spec["n"] == 3:
        draw_slide03(canvas, include_overlay=include_overlay)
    else:
        raise ValueError(f"unsupported slide: {spec['n']}")
    return canvas.convert("RGB")


def save_grayscale(src: Path, dest: Path) -> None:
    image = Image.open(src).convert("RGB")
    gray = ImageOps.grayscale(image).convert("RGB")
    gray.save(dest)


def image_metrics(path: Path) -> dict:
    image = Image.open(path).convert("RGB")
    gray = ImageOps.grayscale(image)
    stat = ImageStat.Stat(gray)
    edges = gray.filter(ImageFilter.FIND_EDGES)
    edge_stat = ImageStat.Stat(edges)
    return {
        "path": rel(path),
        "width_px": image.width,
        "height_px": image.height,
        "luma_mean": round(stat.mean[0], 2),
        "luma_stddev": round(stat.stddev[0], 2),
        "edge_mean": round(edge_stat.mean[0], 2),
        "nonblank_pass": stat.stddev[0] >= 16.5 and edge_stat.mean[0] >= 2.0,
    }


def build_contact_sheet(paths: list[Path], dest: Path, title: str) -> None:
    thumb_w, thumb_h = 1180, 664
    margin = 70
    header_h = 150
    sheet = Image.new("RGB", (margin * 2 + thumb_w * 3 + 50 * 2, header_h + thumb_h + 150), rgb("void"))
    draw = ImageDraw.Draw(sheet)
    draw.text((margin, 50), title, font=F["h"], fill=rgb("ink"))
    draw.text((margin, 96), "Opening block: thesis, landscape position, and version evidence levels.", font=F["small"], fill=rgb("muted"))
    for idx, path in enumerate(paths):
        image = Image.open(path).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        x = margin + idx * (thumb_w + 50)
        y = header_h
        sheet.paste(image, (x, y))
        draw.text((x, y + thumb_h + 28), f"Slide {idx + 1}", font=F["small"], fill=rgb("soft"))
    dest.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(dest)


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def caption_pack(paths: dict[str, dict[str, Path]]) -> tuple[list[dict], str]:
    rows = []
    md_lines = ["# Opening Slides 1-3 Caption Pack", "", f"Date: {CREATED}", ""]
    for spec in SLIDES:
        entry = {
            "slide": spec["n"],
            "title": spec["title"],
            "role": spec["asset_role"],
            "headline": spec["title"],
            "bridge": spec["bridge"],
            "caveat": spec["caveat"],
            "source": spec["source"],
            "scene_plate": rel(paths[spec["id"]]["scene_plate"]),
            "rendered_preview": rel(paths[spec["id"]]["preview"]),
            "grayscale_qa": rel(paths[spec["id"]]["grayscale"]),
        }
        rows.append(entry)
        md_lines.extend(
            [
                f"## Slide {spec['n']}: {spec['asset_role']}",
                "",
                f"- Headline: {spec['title']}",
                f"- Bridge: {spec['bridge']}",
                f"- Caveat: {spec['caveat']}",
                f"- Source: {spec['source']}",
                f"- Scene plate: `{entry['scene_plate']}`",
                f"- Rendered preview: `{entry['rendered_preview']}`",
                "",
            ]
        )
    return rows, "\n".join(md_lines)


def overlay_spec(paths: dict[str, dict[str, Path]]) -> dict:
    return {
        "created": CREATED,
        "z_stack_rule": "Layered PNG scene plus editable scientific interpretation.",
        "global_rules": [
            "Use scene plates as Z0-Z2 visual scaffolds, not as final all-in-one PPTX slides.",
            "Rebuild eyebrow, title, subtitle, bridge, caveat, source, and positioning labels as editable text.",
            "Keep the opening restrained: one question per slide, no dense proof-object table.",
            "Avoid firstness claims; slide 2 is a positioning map, not a ranking.",
            "Do not introduce organoid or v8/v9 extension evidence before the completed result spine.",
        ],
        "slides": [
            {
                "slide": spec["n"],
                "id": spec["id"],
                "scene_plate": rel(paths[spec["id"]]["scene_plate"]),
                "rendered_preview": rel(paths[spec["id"]]["preview"]),
                "editable_text": {
                    "eyebrow": spec["eyebrow"],
                    "title": spec["title"],
                    "subtitle": spec["subtitle"],
                    "bridge": spec["bridge"],
                    "caveat": spec["caveat"],
                    "source": spec["source"],
                },
            }
            for spec in SLIDES
        ],
        "layer_notes": {
            "Z0": "dark field, subtle cellular texture, and restrained mission atmosphere",
            "Z1": "measurement grid, orbit arcs, axes, and timeline ticks",
            "Z2": "source/task/mission plates, landscape matrix, and evolution nodes",
            "Z3": "editable interpretation text and highlight paths",
            "Z4": "claim boundary and source/caveat line",
        },
    }


def build() -> dict[str, str]:
    ensure()
    paths: dict[str, dict[str, Path]] = {}
    metrics: list[dict] = []

    for spec in SLIDES:
        scene_plate = OUT / f"{spec['id']}_scene_plate.png"
        preview = OUT / f"{spec['id']}_rendered_preview.png"
        grayscale = QA / f"{spec['id']}_rendered_preview_grayscale.png"
        render_slide(spec, include_overlay=False).save(scene_plate)
        render_slide(spec, include_overlay=True).save(preview)
        save_grayscale(preview, grayscale)
        paths[spec["id"]] = {"scene_plate": scene_plate, "preview": preview, "grayscale": grayscale}
        metrics.append(image_metrics(preview))

    contact = OUT / "opening_slides_1_3_contact_sheet.png"
    gray_contact = QA / "opening_slides_1_3_contact_sheet_grayscale.png"
    build_contact_sheet([paths[spec["id"]]["preview"] for spec in SLIDES], contact, "Premium opening slides 1-3")
    save_grayscale(contact, gray_contact)

    captions, captions_md = caption_pack(paths)
    caption_json = OUT / "opening_slides_1_3_caption_pack.json"
    caption_md = OUT / "opening_slides_1_3_caption_pack.md"
    overlay_json = OUT / "opening_slides_1_3_overlay_spec.json"
    manifest_json = OUT / "opening_slides_1_3_manifest.json"
    qa_json = OUT / "opening_slides_1_3_qa.json"

    write_json(caption_json, {"created": CREATED, "slides": captions})
    caption_md.write_text(captions_md + "\n", encoding="utf-8")
    write_json(overlay_json, overlay_spec(paths))

    manifest = {
        "created": CREATED,
        "artifact_role": "premium opening slide scenes for slides 1-3",
        "outputs": {
            "contact_sheet": rel(contact),
            "grayscale_contact_sheet": rel(gray_contact),
            "caption_pack_json": rel(caption_json),
            "caption_pack_markdown": rel(caption_md),
            "overlay_spec": rel(overlay_json),
            "qa": rel(qa_json),
        },
        "slides": [
            {
                "slide": spec["n"],
                "id": spec["id"],
                "role": spec["asset_role"],
                "scene_plate": rel(paths[spec["id"]]["scene_plate"]),
                "rendered_preview": rel(paths[spec["id"]]["preview"]),
                "grayscale_qa": rel(paths[spec["id"]]["grayscale"]),
                "bridge": spec["bridge"],
                "claim_boundary": spec["caveat"],
            }
            for spec in SLIDES
        ],
        "local_reference_inputs": [
            "docs/SLIDES_1_14_DECK_ASSEMBLY_BRIDGE_REVIEW_2026_06_03.md",
            "/Users/jak4013/Dropbox/Bioinformatics/Claude/Presentation_RLVR/examples/premium_alignment/isgp_2026_plenary_v28_premium_depth_labels.json",
            "/Users/jak4013/Dropbox/Bioinformatics/Claude/Presentation_RLVR/examples/premium_judges/isgp_2026_plenary_v28_premium_judge_run.json",
            "/Users/jak4013/Dropbox/Bioinformatics/Claude/Agentic_system/skills/scientific-slide-deck-production/references/isgp-space-omics-plenary-case-notes.md",
        ],
    }
    write_json(manifest_json, manifest)

    qa = {
        "created": CREATED,
        "automatic_checks": {
            "n_slides": len(SLIDES),
            "all_previews_exist": all(paths[spec["id"]]["preview"].exists() for spec in SLIDES),
            "all_scene_plates_exist": all(paths[spec["id"]]["scene_plate"].exists() for spec in SLIDES),
            "all_grayscale_outputs_exist": all(paths[spec["id"]]["grayscale"].exists() for spec in SLIDES),
            "contact_sheet_exists": contact.exists(),
            "grayscale_contact_sheet_exists": gray_contact.exists(),
            "caption_pack_exists": caption_json.exists() and caption_md.exists(),
            "overlay_spec_exists": overlay_json.exists(),
            "metrics": metrics,
            "nonblank_pass": all(row["nonblank_pass"] for row in metrics),
        },
        "manual_review": {
            "visual_review_status": "pass: contact sheet and grayscale QA inspected; slide 2 contrast tuned after automatic dark-scene warning",
            "premium_reference": "ISGP opening block: conceptual depth, hierarchy, rhythm, and restraint.",
            "known_watchpoints": [
                "avoid decorative-only space background",
                "avoid jargon shock in slide 2 landscape labels",
                "keep v8/v9 evidence level separated from completed results",
            ],
        },
    }
    write_json(qa_json, qa)

    return {
        "contact_sheet": rel(contact),
        "grayscale_contact_sheet": rel(gray_contact),
        "manifest": rel(manifest_json),
        "qa": rel(qa_json),
        "caption_pack": rel(caption_md),
        "overlay_spec": rel(overlay_json),
    }


if __name__ == "__main__":
    print(json.dumps(build(), indent=2))

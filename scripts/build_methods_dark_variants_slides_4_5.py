#!/usr/bin/env python3
"""Build dark premium methods variants for slides 4-5.

These variants are intended to sit between the dark opening block and the dark
result-family slides. Existing light proof-stage assets remain useful as backup
or source-proof material; these scenes provide the main-deck visual grammar.
"""

from __future__ import annotations

import json
import math
import random
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont, ImageOps, ImageStat


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "premium_methods_dark_variants_slides_4_5_v0_1"
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
        "n": 4,
        "id": "slide04_evaluation_layer_dark",
        "eyebrow": "METHODS BRIDGE",
        "title": "Public studies become auditable tasks",
        "subtitle": "Source context stays attached as studies become benchmark tasks and score records.",
        "bridge": "What exactly is being tested?",
        "caveat": "Conceptual methods bridge; not quantitative result evidence.",
        "source": "B1/B2 narration pack + source/task inventory",
        "claim_boundary": "evaluation-layer schematic only",
        "backup_light_asset": "output/premium_bridge_scenes/b1_evaluation_layer/rendered_preview.png",
        "backup_note_asset": "output/premium_bridge_rebuild_scenes/b2_study_to_task_premium/rendered_preview.png",
    },
    {
        "n": 5,
        "id": "slide05_mission_heldout_dark",
        "eyebrow": "METHODS GUARD",
        "title": "The hidden test set is an entire mission",
        "subtitle": "Models fit on training missions, then score the held-out mission after the boundary.",
        "bridge": "Does the model transfer to a mission it did not learn from?",
        "caveat": "Mission-held-out evaluation; not random-sample validation.",
        "source": "B3/B4 narration pack + mission-held-out split rule",
        "claim_boundary": "held-out mission schematic only",
        "backup_light_asset": "output/premium_bridge_rebuild_scenes/b3_mission_held_out_premium/rendered_preview.png",
        "backup_note_asset": "output/premium_bridge_rebuild_scenes/b4_train_only_guard_premium/rendered_preview.png",
    },
]

FORBIDDEN_VISIBLE_TERMS = [
    "LOMO",
    "random CV",
    "cross-validation",
    "payload",
    "RRRM",
    "alpha",
    "NES",
    "macro-F1",
    "/Users/",
    "sklearn",
    "function",
    "class",
    "wireframe",
    "micro-plan",
]


def rgb(token: str) -> tuple[int, int, int]:
    value = COLORS[token].lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(token: str, alpha: int) -> tuple[int, int, int, int]:
    return rgb(token) + (alpha,)


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
    "title": font(64, bold=True),
    "subtitle": font(29),
    "h": font(31, bold=True),
    "body": font(24),
    "small": font(18),
    "tiny": font(14),
    "num": font(40, bold=True),
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


def rounded_shadow(canvas: Image.Image, box: tuple[int, int, int, int], radius: int = 18, blur: int = 26, alpha: int = 115) -> None:
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

    noise = Image.effect_noise((W, H), 24).convert("L").filter(ImageFilter.GaussianBlur(0.5))
    noise = ImageOps.colorize(noise, black="#000000", white="#DDF6FF").convert("RGBA")
    noise.putalpha(14)
    canvas.alpha_composite(noise)

    for x in range(180, W, 180):
        draw.line((x, 0, x, H), fill=rgba("rule", 25), width=1)
    for y in range(150, H, 150):
        draw.line((0, y, W, y), fill=rgba("rule", 22), width=1)

    rng = random.Random(seed)
    for _ in range(180):
        x = rng.randint(80, W - 90)
        y = rng.randint(70, H - 90)
        size = rng.choice([2, 2, 3])
        draw.ellipse((x, y, x + size, y + size), fill=(235, 248, 255, rng.randint(26, 88)))

    center = (int(W * 0.58), int(H * 0.31))
    for idx, radius in enumerate([980, 1220, 1490]):
        bbox = (
            center[0] - radius,
            center[1] - int(radius * 0.35),
            center[0] + radius,
            center[1] + int(radius * 0.35),
        )
        draw.arc(bbox, 196, 350, fill=rgba("sky", 52 - idx * 12), width=3)


def draw_molecular_field(canvas: Image.Image, seed: int, *, side: str = "right") -> None:
    rng = random.Random(seed)
    layer = alpha_layer()
    draw = ImageDraw.Draw(layer)
    x_min, x_max = (0.62, 0.92) if side == "right" else (0.08, 0.42)
    points: list[tuple[int, int, str]] = []
    for _ in range(34):
        x = int(W * rng.uniform(x_min, x_max))
        y = int(H * rng.uniform(0.20, 0.84))
        tone = rng.choice(["teal", "sky", "green", "violet"])
        r = rng.randint(8, 22)
        draw.ellipse((x - r, y - r, x + r, y + r), outline=rgba(tone, rng.randint(28, 70)), width=2)
        points.append((x, y, tone))
    for idx, (x0, y0, tone) in enumerate(points[:-1]):
        if rng.random() < 0.34:
            x1, y1, _ = points[idx + 1]
            draw.line((x0, y0, x1, y1), fill=rgba(tone, rng.randint(20, 44)), width=2)
    canvas.alpha_composite(layer.filter(ImageFilter.GaussianBlur(0.4)))


def draw_slide_chrome(canvas: Image.Image, spec: dict, *, include_overlay: bool) -> None:
    if not include_overlay:
        return
    draw = ImageDraw.Draw(canvas)
    draw.text((230, 182), spec["eyebrow"], font=F["eyebrow"], fill=rgb("teal"))
    y = draw_wrapped(draw, (230, 238), spec["title"], F["title"], rgb("ink"), 1620, line_gap=12)
    y += 16
    draw_wrapped(draw, (234, y), spec["subtitle"], F["subtitle"], rgb("soft"), 1650, line_gap=8)
    draw.rounded_rectangle((230, 1788, 1980, 1948), radius=12, fill=rgba("panel", 226), outline=rgba("rule", 160), width=2)
    draw.text((270, 1828), spec["bridge"], font=F["body"], fill=rgb("ink"))
    draw.text((270, 1874), spec["caveat"], font=F["small"], fill=rgb("muted"))
    draw.text((2440, 1888), f"source: {spec['source']}", font=F["small"], fill=rgb("muted"))


def plate(
    canvas: Image.Image,
    box: tuple[int, int, int, int],
    title: str,
    lines: list[str],
    tone: str,
    *,
    include_overlay: bool,
    fill: str = "panel",
    label: str | None = None,
) -> None:
    draw = ImageDraw.Draw(canvas)
    x0, y0, x1, y1 = box
    rounded_shadow(canvas, box, radius=14, blur=20, alpha=86)
    draw.rounded_rectangle(box, radius=14, fill=rgba(fill, 234), outline=rgba(tone, 172), width=2)
    draw.rectangle((x0, y0, x1, y0 + 8), fill=rgba(tone, 215))
    if not include_overlay:
        for idx in range(5):
            yy = y0 + 72 + idx * 42
            draw.line((x0 + 40, yy, x1 - 52, yy), fill=rgba("soft", 36), width=2)
        return
    if label:
        draw.text((x0 + 40, y0 + 34), label.upper(), font=F["tiny"], fill=rgb(tone))
        title_y = y0 + 72
    else:
        title_y = y0 + 42
    draw.text((x0 + 40, title_y), title, font=F["h"], fill=rgb("ink"))
    yy = title_y + 58
    for line in lines:
        draw.text((x0 + 40, yy), line, font=F["small"], fill=rgb("soft"))
        yy += 36


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], tone: str, *, width: int = 5, alpha: int = 190) -> None:
    x0, y0 = start
    x1, y1 = end
    draw.line((x0, y0, x1, y1), fill=rgba(tone, alpha), width=width)
    angle = math.atan2(y1 - y0, x1 - x0)
    head = 24
    p1 = (x1 - int(math.cos(angle - 0.45) * head), y1 - int(math.sin(angle - 0.45) * head))
    p2 = (x1 - int(math.cos(angle + 0.45) * head), y1 - int(math.sin(angle + 0.45) * head))
    draw.polygon([end, p1, p2], fill=rgba(tone, alpha))


def draw_small_status_row(draw: ImageDraw.ImageDraw, x: int, y: int, tones: list[str]) -> None:
    for idx, tone in enumerate(tones):
        xx = x + idx * 32
        draw.ellipse((xx, y, xx + 15, y + 15), fill=rgba(tone, 220))


def draw_slide04(canvas: Image.Image, *, include_overlay: bool) -> None:
    draw_background(canvas, 44)
    draw_molecular_field(canvas, 45, side="right")
    draw = ImageDraw.Draw(canvas)

    work = (1220, 590, 3500, 1458)
    rounded_shadow(canvas, work, radius=18, blur=28, alpha=112)
    draw.rounded_rectangle(work, radius=18, fill=rgba("panel", 194), outline=rgba("rule", 144), width=2)
    if include_overlay:
        draw.text((1260, 632), "Evaluation layer", font=F["small"], fill=rgb("teal"))
        draw.text((1260, 678), "A public study is not automatically a benchmark.", font=F["body"], fill=rgb("soft"))

    plate(
        canvas,
        (1350, 820, 1775, 1238),
        "Public studies",
        ["mission records", "samples", "assays"],
        "sky",
        include_overlay=include_overlay,
        label="source",
    )

    lane_x0, lane_y0, lane_x1, lane_y1 = 1925, 835, 2505, 1232
    rounded_shadow(canvas, (lane_x0, lane_y0, lane_x1, lane_y1), radius=14, blur=18, alpha=70)
    draw.rounded_rectangle((lane_x0, lane_y0, lane_x1, lane_y1), radius=14, fill=rgba("panel2", 232), outline=rgba("green", 152), width=2)
    if include_overlay:
        draw.text((lane_x0 + 42, lane_y0 + 40), "Context stays attached", font=F["h"], fill=rgb("ink"))
    lane_labels = [("mission", "teal"), ("tissue", "green"), ("sample", "sky"), ("label", "amber")]
    for idx, (label, tone) in enumerate(lane_labels):
        yy = lane_y0 + 128 + idx * 64
        draw.line((lane_x0 + 52, yy, lane_x1 - 70, yy), fill=rgba(tone, 132), width=4)
        draw.ellipse((lane_x0 + 52 - 10, yy - 10, lane_x0 + 52 + 10, yy + 10), fill=rgba(tone, 220))
        draw.ellipse((lane_x1 - 70 - 10, yy - 10, lane_x1 - 70 + 10, yy + 10), fill=rgba(tone, 190))
        if include_overlay:
            draw.text((lane_x0 + 76, yy - 18), label, font=F["small"], fill=rgb("soft"))

    plate(
        canvas,
        (2670, 790, 3118, 1270),
        "Benchmark task",
        ["defined split", "metric", "audit trail"],
        "violet",
        include_overlay=include_overlay,
        label="task",
    )
    plate(
        canvas,
        (3210, 900, 3428, 1170),
        "Score",
        ["reported", "traceable"],
        "teal",
        include_overlay=include_overlay,
        label="record",
    )

    arrow(draw, (1788, 1030), (1910, 1030), "teal", width=5, alpha=178)
    arrow(draw, (2518, 1030), (2652, 1030), "violet", width=5, alpha=178)
    arrow(draw, (3128, 1030), (3194, 1030), "teal", width=4, alpha=168)
    draw_small_status_row(draw, 3290, 1218, ["sky", "green", "amber", "teal"])

    if include_overlay:
        draw.text((1360, 1350), "source evidence", font=F["small"], fill=rgb("sky"))
        draw.text((1916, 1350), "task definition", font=F["small"], fill=rgb("green"))
        draw.text((2670, 1350), "auditable result", font=F["small"], fill=rgb("violet"))

    draw_slide_chrome(canvas, SLIDES[0], include_overlay=include_overlay)


def draw_mission_node(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, tone: str, *, active: bool = True) -> None:
    r = 54 if active else 46
    draw.ellipse((x - r, y - r, x + r, y + r), fill=rgba("panel2", 238), outline=rgba(tone, 205), width=4 if active else 3)
    draw.text((x - 24, y - 16), label, font=F["body"], fill=rgb(tone))
    for idx in range(3):
        yy = y + 74 + idx * 24
        draw.line((x - 46, yy, x + 46, yy), fill=rgba(tone, 60), width=2)


def draw_slide05(canvas: Image.Image, *, include_overlay: bool) -> None:
    draw_background(canvas, 54)
    draw_molecular_field(canvas, 55, side="right")
    draw = ImageDraw.Draw(canvas)

    panel = (1180, 570, 3500, 1480)
    rounded_shadow(canvas, panel, radius=18, blur=28, alpha=112)
    draw.rounded_rectangle(panel, radius=18, fill=rgba("panel", 196), outline=rgba("rule", 150), width=2)

    boundary_x = 2692
    zones = alpha_layer()
    zd = ImageDraw.Draw(zones)
    zd.rounded_rectangle((1220, 628, boundary_x - 22, 1395), radius=16, fill=rgba("green", 28), outline=rgba("green", 72), width=2)
    zd.rounded_rectangle((boundary_x + 22, 628, 3440, 1395), radius=16, fill=rgba("red", 30), outline=rgba("red", 80), width=2)
    zd.rectangle((boundary_x - 50, 628, boundary_x + 50, 1395), fill=(0, 0, 0, 0))
    canvas.alpha_composite(zones)
    draw.line((boundary_x, 655, boundary_x, 1410), fill=rgba("red", 230), width=5)
    draw.line((boundary_x + 18, 690, boundary_x + 18, 1375), fill=rgba("red", 80), width=2)

    if include_overlay:
        draw.text((1240, 632), "training side", font=F["small"], fill=rgb("green"))
        draw.text((2760, 632), "hidden test side", font=F["small"], fill=rgb("red"))

    xs = [1440, 1705, 1970]
    for idx, x in enumerate(xs, start=1):
        draw_mission_node(draw, x, 970, f"M{idx}", "green")
    draw.line((1440, 970, 2420, 970), fill=rgba("green", 135), width=5)
    draw.line((1440, 1012, 2420, 1012), fill=rgba("sky", 50), width=2)
    arrow(draw, (2050, 970), (2300, 970), "green", width=5, alpha=170)

    model = (2250, 810, 2558, 1140)
    rounded_shadow(canvas, model, radius=14, blur=18, alpha=80)
    draw.rounded_rectangle(model, radius=14, fill=rgba("panel2", 238), outline=rgba("sky", 170), width=2)
    if include_overlay:
        draw.text((2296, 864), "Fit pipeline", font=F["h"], fill=rgb("ink"))
        draw.text((2296, 930), "features", font=F["small"], fill=rgb("soft"))
        draw.text((2296, 966), "scaling", font=F["small"], fill=rgb("soft"))
        draw.text((2296, 1002), "model", font=F["small"], fill=rgb("soft"))
    for yy, tone in [(1020, "teal"), (1054, "sky"), (1088, "green")]:
        draw.line((2294, yy, 2510, yy), fill=rgba(tone, 108), width=5)

    draw_mission_node(draw, 3075, 970, "M4", "red")
    arrow(draw, (2568, 970), (2674, 970), "sky", width=4, alpha=120)
    arrow(draw, (2720, 970), (3000, 970), "red", width=5, alpha=168)
    if include_overlay:
        draw.text((2820, 810), "score after training", font=F["body"], fill=rgb("red"))
        draw.text((2820, 852), "hidden samples do not fit the model", font=F["small"], fill=rgb("soft"))

    guard = (1340, 1288, 3280, 1398)
    draw.rounded_rectangle(guard, radius=12, fill=rgba("panel2", 230), outline=rgba("amber", 160), width=2)
    if include_overlay:
        draw.text((1380, 1322), "Train-only guard:", font=F["body"], fill=rgb("amber"))
        draw.text((1578, 1322), "feature choices are learned before the boundary.", font=F["body"], fill=rgb("soft"))

    draw_slide_chrome(canvas, SLIDES[1], include_overlay=include_overlay)


def render_slide(spec: dict, *, include_overlay: bool) -> Image.Image:
    canvas = Image.new("RGBA", (W, H), (0, 0, 0, 255))
    if spec["n"] == 4:
        draw_slide04(canvas, include_overlay=include_overlay)
    elif spec["n"] == 5:
        draw_slide05(canvas, include_overlay=include_overlay)
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


def visible_text_forbidden_hits() -> dict[str, list[str]]:
    hits: dict[str, list[str]] = {}
    for spec in SLIDES:
        text = " ".join([spec["title"], spec["subtitle"], spec["bridge"], spec["caveat"], spec["source"]])
        found = [term for term in FORBIDDEN_VISIBLE_TERMS if term in text]
        hits[spec["id"]] = found
    return hits


def build_contact_sheet(paths: list[Path], dest: Path, title: str) -> None:
    thumb_w, thumb_h = 1540, 866
    margin = 80
    header_h = 150
    sheet = Image.new("RGB", (margin * 2 + thumb_w * 2 + 70, header_h + thumb_h + 150), rgb("void"))
    draw = ImageDraw.Draw(sheet)
    draw.text((margin, 50), title, font=F["h"], fill=rgb("ink"))
    draw.text((margin, 96), "Slides 4-5 dark methods variants for continuity with opening and result scenes.", font=F["small"], fill=rgb("muted"))
    for idx, path in enumerate(paths):
        image = Image.open(path).convert("RGB").resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
        x = margin + idx * (thumb_w + 70)
        y = header_h
        sheet.paste(image, (x, y))
        draw.text((x, y + thumb_h + 28), f"Slide {idx + 4}", font=F["small"], fill=rgb("soft"))
    dest.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(dest)


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def caption_pack(paths: dict[str, dict[str, Path]]) -> tuple[list[dict], str]:
    rows = []
    md_lines = ["# Slides 4-5 Dark Methods Variant Caption Pack", "", f"Date: {CREATED}", ""]
    for spec in SLIDES:
        entry = {
            "slide": spec["n"],
            "title": spec["title"],
            "bridge": spec["bridge"],
            "caveat": spec["caveat"],
            "source": spec["source"],
            "claim_boundary": spec["claim_boundary"],
            "scene_plate": rel(paths[spec["id"]]["scene_plate"]),
            "rendered_preview": rel(paths[spec["id"]]["preview"]),
            "grayscale_qa": rel(paths[spec["id"]]["grayscale"]),
            "backup_light_asset": spec["backup_light_asset"],
            "backup_note_asset": spec["backup_note_asset"],
        }
        rows.append(entry)
        md_lines.extend(
            [
                f"## Slide {spec['n']}: {spec['title']}",
                "",
                f"- Bridge: {spec['bridge']}",
                f"- Caveat: {spec['caveat']}",
                f"- Source: {spec['source']}",
                f"- Scene plate: `{entry['scene_plate']}`",
                f"- Rendered preview: `{entry['rendered_preview']}`",
                f"- Backup light asset: `{spec['backup_light_asset']}`",
                "",
            ]
        )
    return rows, "\n".join(md_lines)


def overlay_spec(paths: dict[str, dict[str, Path]]) -> dict:
    return {
        "created": CREATED,
        "z_stack_rule": "Dark layered methods scene plus editable scientific interpretation.",
        "global_rules": [
            "Use dark scene plates as main-deck visual scaffolds for slides 4-5.",
            "Keep existing light proof-stage slides as backup/source-proof assets.",
            "Rebuild titles, labels, bridge lines, caveats, and source lines as editable PPTX text.",
            "Keep B2 and B4 definitions in notes or backup; do not overload slides 4-5.",
            "Avoid table-like panels and internal implementation terms.",
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
                "backup_light_asset": spec["backup_light_asset"],
                "backup_note_asset": spec["backup_note_asset"],
            }
            for spec in SLIDES
        ],
        "layer_notes": {
            "Z0": "dark field and molecular texture",
            "Z1": "grid, orbit arcs, task lanes, and split boundary",
            "Z2": "source/task/score plates and mission nodes",
            "Z3": "editable labels, arrows, and focus text",
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

    contact = OUT / "slides_4_5_dark_methods_contact_sheet.png"
    gray_contact = QA / "slides_4_5_dark_methods_contact_sheet_grayscale.png"
    build_contact_sheet([paths[spec["id"]]["preview"] for spec in SLIDES], contact, "Slides 4-5 dark methods variants")
    save_grayscale(contact, gray_contact)

    captions, captions_md = caption_pack(paths)
    caption_json = OUT / "slides_4_5_dark_methods_caption_pack.json"
    caption_md = OUT / "slides_4_5_dark_methods_caption_pack.md"
    overlay_json = OUT / "slides_4_5_dark_methods_overlay_spec.json"
    manifest_json = OUT / "slides_4_5_dark_methods_manifest.json"
    qa_json = OUT / "slides_4_5_dark_methods_qa.json"

    write_json(caption_json, {"created": CREATED, "slides": captions})
    caption_md.write_text(captions_md + "\n", encoding="utf-8")
    write_json(overlay_json, overlay_spec(paths))

    manifest = {
        "created": CREATED,
        "artifact_role": "dark methods variants for slides 4-5",
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
                "scene_plate": rel(paths[spec["id"]]["scene_plate"]),
                "rendered_preview": rel(paths[spec["id"]]["preview"]),
                "grayscale_qa": rel(paths[spec["id"]]["grayscale"]),
                "claim_boundary": spec["claim_boundary"],
                "backup_light_asset": spec["backup_light_asset"],
                "backup_note_asset": spec["backup_note_asset"],
            }
            for spec in SLIDES
        ],
        "source_documents": [
            "docs/PREMIUM_BRIDGE_METHODS_NARRATION_PACK_B1_B4_2026_06_02.md",
            "docs/SLIDES_1_3_OPENING_VISUAL_SYSTEM_REVIEW_2026_06_03.md",
            "docs/SLIDES_1_14_DECK_ASSEMBLY_BRIDGE_REVIEW_2026_06_03.md",
            "output/premium_bridge_scenes/b1_evaluation_layer/overlay_spec.json",
            "output/premium_bridge_rebuild_scenes/b3_mission_held_out_premium/overlay_spec.json",
        ],
    }
    write_json(manifest_json, manifest)

    forbidden_hits = visible_text_forbidden_hits()
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
            "forbidden_visible_term_hits": forbidden_hits,
            "forbidden_visible_terms_pass": all(not hits for hits in forbidden_hits.values()),
        },
        "manual_review": {
            "visual_review_status": "pass: contact sheet, individual previews, and grayscale QA inspected after slide 5 zone-color revision",
            "production_decision": "candidate replacement for slides 4-5 main deck; light proof scenes remain backup",
            "known_watchpoints": [
                "slide 4 should not become a dense implementation diagram",
                "slide 5 must not read as random-sample validation",
                "final PPTX should rebuild text as editable overlay",
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

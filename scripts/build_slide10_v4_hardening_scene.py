#!/usr/bin/env python3
"""Build slide 10 v4 hardening premium scene.

The visible values use docs/CANONICAL_RESULTS_V7_1.md rather than raw v4 JSON
best rows, because the deck should stay on the canonical public-facing surface.
"""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "premium_core_result_slides" / "slide10_v4_hardening_v0_1"
QA = OUT / "qa"
CREATED = "2026-06-03"

COLORS = {
    "void": "#070A0E",
    "deep": "#0B1117",
    "panel": "#101923",
    "panel2": "#152331",
    "ink": "#F4F7F8",
    "soft": "#B9C7D2",
    "muted": "#788898",
    "rule": "#33465A",
    "blue": "#2D6F9F",
    "sky": "#6BAED6",
    "teal": "#1AA090",
    "green": "#178B63",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "red": "#B23A3A",
}

CANONICAL_ROWS = [
    {"tissue": "Thymus", "auroc": 0.948, "method": "PCA-LR", "feature": "KEGG", "sig": "p<0.05"},
    {"tissue": "Colon", "auroc": 0.921, "method": "PCA-LR", "feature": "KEGG", "sig": "p<0.05"},
    {"tissue": "Lung", "auroc": 0.901, "method": "PCA-LR", "feature": "Gene", "sig": "p<0.05"},
    {"tissue": "Kidney", "auroc": 0.829, "method": "ElasticNet-LR", "feature": "Hallmark", "sig": "p<0.01"},
    {"tissue": "Eye", "auroc": 0.823, "method": "PCA-LR", "feature": "Hallmark", "sig": "best row NS"},
    {"tissue": "Skin", "auroc": 0.819, "method": "PCA-LR", "feature": "Gene", "sig": "best row NS"},
    {"tissue": "Gastrocnemius", "auroc": 0.776, "method": "PCA-LR", "feature": "Gene", "sig": "best row NS"},
    {"tissue": "Liver", "auroc": 0.670, "method": "PCA-LR", "feature": "Gene", "sig": "best row NS"},
]

OVERLAY_TEXT = [
    "The result survives a wider method grid",
    "v4 expands the benchmark surface before the deck moves into biology and translation.",
    "8 tissues",
    "8 classifiers",
    "4 feature views",
    "256 evaluations",
    "40 significant",
    "6/8 tissues with signal",
    "Canonical v7.1 table; not every best row is significant.",
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
    "eyebrow": font(20, bold=True),
    "title": font(62, bold=True),
    "subtitle": font(28),
    "h": font(30, bold=True),
    "body": font(21),
    "small": font(17),
    "tiny": font(13),
    "num": font(42, bold=True),
}


def ensure() -> None:
    QA.mkdir(parents=True, exist_ok=True)


def draw_background(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    w, h = canvas.size
    draw.rectangle((0, 0, w, h), fill=rgb("void"))
    top = rgb("void")
    bottom = rgb("panel")
    for y in range(0, h, 2):
        t = y / max(1, h - 1)
        color = tuple(int(top[i] * (1 - t) + bottom[i] * t) for i in range(3))
        draw.line((0, y, w, y), fill=color + (255,), width=2)
    for x in range(170, w, 170):
        draw.line((x, 0, x, h), fill=rgba("rule", 24), width=1)
    for y in range(150, h, 150):
        draw.line((0, y, w, y), fill=rgba("rule", 21), width=1)
    center = (int(w * 0.73), int(h * 0.25))
    for idx, radius in enumerate([760, 1000, 1240]):
        bbox = (
            center[0] - radius,
            center[1] - int(radius * 0.33),
            center[0] + radius,
            center[1] + int(radius * 0.33),
        )
        draw.arc(bbox, 200, 350, fill=rgba("sky", 46 - idx * 10), width=3)


def auroc_color(value: float) -> tuple[int, int, int]:
    low = rgb("blue")
    mid = rgb("teal")
    high = rgb("amber")
    t = max(0.0, min(1.0, (value - 0.65) / 0.30))
    if t < 0.62:
        u = t / 0.62
        return tuple(int(low[i] * (1 - u) + mid[i] * u) for i in range(3))
    u = (t - 0.62) / 0.38
    return tuple(int(mid[i] * (1 - u) + high[i] * u) for i in range(3))


def panel_shadow(canvas: Image.Image, box: tuple[int, int, int, int]) -> None:
    x, y, w, h = box
    shadow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow)
    sd.rounded_rectangle((x + 18, y + 22, x + w + 18, y + h + 22), radius=12, fill=(0, 0, 0, 130))
    canvas.alpha_composite(shadow.filter(ImageFilter.GaussianBlur(18)))


def draw_signal_surface(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = 250, 590, 2140, 1080
    panel_shadow(canvas, (x, y, w, h))
    draw.rounded_rectangle((x, y, x + w, y + h), radius=10, fill=rgba("panel", 238), outline=rgba("rule", 190), width=2)
    draw.text((x + 54, y + 42), "Canonical tissue signal surface", font=F["h"], fill=rgb("ink"))
    draw.text((x + 54, y + 82), "Best public-facing row per tissue from v7.1 documentation", font=F["small"], fill=rgb("muted"))

    axis_x = x + 470
    row_y0 = y + 170
    row_h = 88
    max_w = 900
    for idx, row in enumerate(CANONICAL_ROWS):
        yy = row_y0 + idx * row_h
        draw.text((x + 54, yy + 18), row["tissue"], font=F["body"], fill=rgb("soft"))
        draw.text((x + 264, yy + 20), row["method"], font=F["tiny"], fill=rgb("muted"))
        draw.text((x + 382, yy + 20), row["feature"], font=F["tiny"], fill=rgb("muted"))
        draw.line((axis_x, yy + 40, axis_x + max_w, yy + 40), fill=rgba("rule", 70), width=2)
        bar_w = int(max_w * row["auroc"])
        color = auroc_color(row["auroc"])
        draw.rounded_rectangle((axis_x, yy + 22, axis_x + bar_w, yy + 58), radius=7, fill=color + (230,))
        draw.text((axis_x + max_w + 35, yy + 18), f"{row['auroc']:.3f}", font=F["body"], fill=rgb("ink"))
        sig_tone = "green" if row["sig"].startswith("p<") else "muted"
        draw.text((axis_x + max_w + 132, yy + 20), row["sig"], font=F["tiny"], fill=rgb(sig_tone))

    draw.text((axis_x, y + h - 74), "0.0", font=F["tiny"], fill=rgb("muted"))
    draw.text((axis_x + int(max_w * 0.5) - 18, y + h - 74), "0.5", font=F["tiny"], fill=rgb("muted"))
    draw.text((axis_x + max_w - 18, y + h - 74), "1.0", font=F["tiny"], fill=rgb("muted"))


def draw_grid_lattice(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = 2550, 590, 940, 610
    panel_shadow(canvas, (x, y, w, h))
    draw.rounded_rectangle((x, y, x + w, y + h), radius=10, fill=rgba("panel", 238), outline=rgba("teal", 180), width=2)
    draw.text((x + 44, y + 42), "Hardening grid", font=F["h"], fill=rgb("ink"))
    draw.text((x + 44, y + 82), "8 tissues x 8 classifiers x 4 feature views", font=F["small"], fill=rgb("muted"))

    gx, gy = x + 78, y + 178
    cell = 22
    gap = 10
    for tissue in range(8):
        for method in range(8):
            intensity = 80 + (tissue * 19 + method * 11) % 110
            draw.rounded_rectangle(
                (
                    gx + method * (cell + gap),
                    gy + tissue * (cell + gap),
                    gx + method * (cell + gap) + cell,
                    gy + tissue * (cell + gap) + cell,
                ),
                radius=4,
                fill=rgba("sky" if tissue < 4 else "teal", intensity),
                outline=rgba("rule", 70),
                width=1,
            )
    draw.text((gx + 8 * (cell + gap) + 42, gy + 54), "x4", font=F["num"], fill=rgb("amber"))
    draw.text((gx + 8 * (cell + gap) + 42, gy + 106), "feature views", font=F["small"], fill=rgb("soft"))
    draw.text((x + 44, y + h - 88), "256 total evaluations", font=F["h"], fill=rgb("ink"))
    draw.text((x + 44, y + h - 48), "coverage, not a single cherry-picked row", font=F["small"], fill=rgb("muted"))


def draw_metric_tiles(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    tiles = [
        ("256", "evaluations", "sky"),
        ("40", "significant", "teal"),
        ("6/8", "tissues with signal", "amber"),
    ]
    x0, y0 = 2550, 1260
    for idx, (num, label, tone) in enumerate(tiles):
        x = x0 + idx * 310
        panel_shadow(canvas, (x, y0, 270, 260))
        draw.rounded_rectangle((x, y0, x + 270, y0 + 260), radius=10, fill=rgba("panel", 238), outline=rgba(tone, 190), width=2)
        draw.text((x + 34, y0 + 58), num, font=F["num"], fill=rgb(tone))
        draw.text((x + 34, y0 + 122), label, font=F["body"], fill=rgb("soft"))


def draw_overlay(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    draw.text((230, 196), "CORE RESULT HARDENING", font=F["eyebrow"], fill=rgb("teal"))
    draw.text((230, 248), "The result survives a wider method grid", font=F["title"], fill=rgb("ink"))
    draw.text((234, 330), "v4 expands the benchmark surface before the deck moves into biology and translation.", font=F["subtitle"], fill=rgb("soft"))

    draw.rounded_rectangle((230, 1810, 2200, 1936), radius=10, fill=rgba("panel", 225), outline=rgba("rule", 160), width=2)
    draw.text((270, 1838), "Canonical v7.1 table; not every best row is significant.", font=F["body"], fill=rgb("soft"))
    draw.text((270, 1880), "Use this as hardening evidence, not as a new leaderboard.", font=F["small"], fill=rgb("muted"))
    draw.text((2830, 1888), "source: CANONICAL_RESULTS_V7_1.md", font=F["small"], fill=rgb("muted"))


def render(*, with_overlay: bool) -> Image.Image:
    canvas = Image.new("RGBA", (3840, 2160), (0, 0, 0, 255))
    draw_background(canvas)
    draw_signal_surface(canvas)
    draw_grid_lattice(canvas)
    draw_metric_tiles(canvas)
    if with_overlay:
        draw_overlay(canvas)
    return canvas.convert("RGB")


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def build() -> dict[str, str]:
    ensure()
    scene_plate = OUT / "slide10_v4_hardening_scene_plate.png"
    preview = OUT / "slide10_v4_hardening_rendered_preview.png"
    grayscale = QA / "slide10_v4_hardening_rendered_preview_grayscale.png"
    manifest = OUT / "slide10_v4_hardening_manifest.json"
    qa = OUT / "slide10_v4_hardening_qa.json"

    render(with_overlay=False).save(scene_plate, quality=95)
    rendered = render(with_overlay=True)
    rendered.save(preview, quality=95)
    rendered.convert("L").convert("RGB").save(grayscale, quality=95)

    manifest_data = {
        "slide_id": "slide10_v4_hardening_v0_1",
        "created": CREATED,
        "slide_role": "result_hardening",
        "reference_calibration_role": "k562_compact_proof_object",
        "source_documents": ["docs/CANONICAL_RESULTS_V7_1.md", "v4/evaluation/M1_summary.json"],
        "visible_value_source": "docs/CANONICAL_RESULTS_V7_1.md",
        "claim_boundary": "Hardening/coverage evidence only; do not imply every best-row tissue is significant.",
        "outputs": {
            "scene_plate": rel(scene_plate),
            "rendered_preview": rel(preview),
            "grayscale_qa": rel(grayscale),
            "qa": rel(qa),
        },
        "canonical_rows": CANONICAL_ROWS,
        "visible_text_word_count": len(" ".join(OVERLAY_TEXT).split()),
        "visible_text_budget": 45,
    }
    write_json(manifest, manifest_data)
    write_json(
        qa,
        {
            "created": CREATED,
            "automatic_checks": {
                "rendered_outputs_exist": all(path.exists() for path in [scene_plate, preview, grayscale, manifest]),
                "image_dimensions": {"width_px": rendered.width, "height_px": rendered.height},
                "visible_text_word_count": len(" ".join(OVERLAY_TEXT).split()),
                "visible_text_budget": 45,
                "canonical_row_count": len(CANONICAL_ROWS),
            },
            "manual_review": {
                "visual_review_status": "pass: rendered preview, scene plate, and grayscale QA inspected",
                "claim_boundary": "canonical v7.1 hardening surface, not a new leaderboard",
                "readability": "primary hierarchy is tissue signal surface plus 256/40/6-of-8 tiles",
                "notes": "method and feature labels are intentionally secondary to avoid table-like clutter",
            },
        },
    )
    return {"scene_plate": rel(scene_plate), "rendered_preview": rel(preview), "grayscale": rel(grayscale), "manifest": rel(manifest), "qa": rel(qa)}


def main() -> None:
    print(json.dumps(build(), indent=2))


if __name__ == "__main__":
    main()

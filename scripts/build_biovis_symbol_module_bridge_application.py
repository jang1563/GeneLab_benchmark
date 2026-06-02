#!/usr/bin/env python3
"""Apply the v0.3 biological symbol modules to real bridge-style layouts.

The output is a visual QA pack, not a final deck. It tests whether the reusable
symbols/modules help explain SpaceBio-Bench methods on light and dark canvases.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "output" / ".matplotlib"))

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "biovis_symbol_module_bridge_application"
ASSET = ROOT / "assets" / "biovis_symbol_module_pack_v0_3"
MOTIF = ROOT / "assets" / "biovis_motif_pack_v0_3" / "png"

CREATED = "2026-06-02"
W, H = 3840, 2160

COLORS = {
    "paper": "#F7F4EC",
    "paper2": "#FBFAF7",
    "ink": "#17212B",
    "muted": "#5D6978",
    "rule": "#AEB8C5",
    "blue": "#2D6F9F",
    "sky": "#6BAED6",
    "green": "#178B63",
    "teal": "#1AA090",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "purple": "#6750A4",
    "red": "#B23A3A",
    "deep": "#0D1720",
    "deep2": "#111D27",
    "deep3": "#172636",
    "white": "#F2F6F8",
    "soft": "#B8C4CF",
}


def rgb(token: str) -> tuple[int, int, int]:
    value = COLORS[token].lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(token: str, alpha: int) -> tuple[int, int, int, int]:
    return rgb(token) + (alpha,)


def font(size: int, *, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


FONTS = {
    "eyebrow": font(28, bold=True),
    "title": font(74, bold=True),
    "subtitle": font(34),
    "h2": font(38, bold=True),
    "h3": font(28, bold=True),
    "body": font(26),
    "small": font(20),
    "tiny": font(16),
    "badge": font(24, bold=True),
}


def ensure() -> None:
    OUT.mkdir(parents=True, exist_ok=True)


def text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], value: str, style: str, color: str, *, anchor: str | None = None) -> None:
    draw.text(xy, value, font=FONTS[style], fill=rgb(color), anchor=anchor)


def line(draw: ImageDraw.ImageDraw, xy: tuple[int, int, int, int], color: str, *, width: int = 3, alpha: int = 255) -> None:
    draw.line(xy, fill=rgba(color, alpha), width=width)


def paste_fit(canvas: Image.Image, path: Path, box: tuple[int, int, int, int], *, opacity: float = 1.0) -> None:
    image = Image.open(path).convert("RGBA")
    x, y, w, h = box
    scale = min(w / image.width, h / image.height)
    size = (max(1, int(image.width * scale)), max(1, int(image.height * scale)))
    image = image.resize(size, Image.Resampling.LANCZOS)
    if opacity < 1.0:
        alpha = image.getchannel("A")
        alpha = alpha.point(lambda p: int(p * opacity))
        image.putalpha(alpha)
    px = x + (w - image.width) // 2
    py = y + (h - image.height) // 2
    canvas.alpha_composite(image, (px, py))


def paste_cover(canvas: Image.Image, path: Path, box: tuple[int, int, int, int], *, opacity: float = 1.0) -> None:
    image = Image.open(path).convert("RGBA")
    x, y, w, h = box
    scale = max(w / image.width, h / image.height)
    size = (max(1, int(image.width * scale)), max(1, int(image.height * scale)))
    image = image.resize(size, Image.Resampling.LANCZOS)
    left = (image.width - w) // 2
    top = (image.height - h) // 2
    image = image.crop((left, top, left + w, top + h))
    if opacity < 1.0:
        alpha = image.getchannel("A")
        alpha = alpha.point(lambda p: int(p * opacity))
        image.putalpha(alpha)
    canvas.alpha_composite(image, (x, y))


def draw_badge(draw: ImageDraw.ImageDraw, xy: tuple[int, int], label: str, tone: str, *, dark: bool = False) -> None:
    x, y = xy
    label_w = max(136, int(draw.textlength(label, font=FONTS["badge"])) + 54)
    fill = "deep2" if dark else "paper2"
    text_color = "white" if dark else "ink"
    draw.rounded_rectangle((x, y, x + label_w, y + 52), radius=22, fill=rgba(fill, 235), outline=rgb(tone), width=3)
    draw.ellipse((x + 18, y + 18, x + 32, y + 32), fill=rgb(tone))
    text(draw, (x + 44, y + 14), label, "badge", text_color)


def draw_proof_slot(draw: ImageDraw.ImageDraw, canvas: Image.Image, xywh: tuple[int, int, int, int], motif: str, title: str, note: str, *, dark: bool = False) -> None:
    x, y, w, h = xywh
    fill = "deep2" if dark else "paper2"
    outline = "soft" if dark else "rule"
    label = "white" if dark else "ink"
    body = "soft" if dark else "muted"
    draw.rounded_rectangle((x, y, x + w, y + h), radius=18, fill=rgba(fill, 232), outline=rgba(outline, 160), width=2)
    paste_cover(canvas, MOTIF / motif, (x + 18, y + 18, w - 36, h - 118), opacity=0.88)
    draw.rectangle((x + 18, y + h - 96, x + w - 18, y + h - 18), fill=rgba(fill, 230))
    text(draw, (x + 36, y + h - 84), title, "h3", label)
    text(draw, (x + 36, y + h - 48), note, "small", body)


def draw_focus_ring(draw: ImageDraw.ImageDraw, xywh: tuple[int, int, int, int], tone: str) -> None:
    x, y, w, h = xywh
    draw.rounded_rectangle((x, y, x + w, y + h), radius=28, outline=rgb(tone), width=5)
    line(draw, (x + 22, y + h - 22, x + 140, y + h - 22), tone, width=5)


def light_slide() -> Image.Image:
    canvas = Image.new("RGBA", (W, H), rgb("paper") + (255,))
    draw = ImageDraw.Draw(canvas)
    for y in [248, 1070, 1910]:
        line(draw, (160, y, W - 160, y), "rule", width=2, alpha=120)

    text(draw, (160, 112), "BIOVIS MODULE APPLICATION TEST", "eyebrow", "blue")
    text(draw, (160, 164), "Can reusable biological modules explain the method?", "title", "ink")
    text(draw, (160, 258), "Light-surface bridge using v0.3 symbols, status badges, and method modules", "subtitle", "muted")

    paste_fit(canvas, ASSET / "modules" / "png" / "sample_to_feature_stack.png", (170, 332, 2120, 620))
    paste_fit(canvas, ASSET / "modules" / "png" / "claim_boundary_bar.png", (170, 1725, 1900, 390))

    draw_proof_slot(
        draw,
        canvas,
        (2380, 350, 650, 515),
        "omics_matrix_texture_field.png",
        "Reserved proof plate",
        "Replace schematic context with source-derived evidence",
    )
    draw_proof_slot(
        draw,
        canvas,
        (3100, 350, 560, 515),
        "pathway_network_summary_field.png",
        "Feature layer context",
        "Generated schematic, not measured values",
    )

    draw_focus_ring(draw, (170, 332, 2120, 620), "teal")
    draw_badge(draw, (2380, 900), "schematic", "muted")
    draw_badge(draw, (2570, 900), "source slot", "blue")
    draw_badge(draw, (2790, 900), "caveat", "red")

    text(draw, (2380, 1060), "What improves vs plain boxes", "h2", "ink")
    bullets = [
        ("Biological scale is visible", "sample, assay, matrix, program, task are no longer abstract rectangles"),
        ("Claim boundary is explicit", "schematic context, source proof, processed result, and validation are separated"),
        ("Future proof slots are clear", "the module reserves room for real dataset screenshots or result panels"),
    ]
    y = 1130
    for idx, (head, body) in enumerate(bullets, start=1):
        draw.ellipse((2380, y + 7, 2418, y + 45), fill=rgb(["blue", "green", "amber"][idx - 1]))
        text(draw, (2448, y), head, "h3", "ink")
        text(draw, (2448, y + 42), body, "body", "muted")
        y += 150

    paste_fit(canvas, ASSET / "modules" / "png" / "species_coverage_strip.png", (2180, 1580, 1500, 345))
    text(draw, (2180, 1995), "Use species symbols as compact scope markers, not hero graphics.", "small", "muted")

    return canvas.convert("RGB")


def dark_slide() -> Image.Image:
    canvas = Image.new("RGBA", (W, H), rgb("deep") + (255,))
    draw = ImageDraw.Draw(canvas)
    for x in [360, 1920, 3360]:
        line(draw, (x, 180, x, H - 190), "deep3", width=2, alpha=165)
    for y in [300, 1180, 1890]:
        line(draw, (160, y, W - 160, y), "deep3", width=2, alpha=200)

    text(draw, (160, 112), "BIOVIS MODULE APPLICATION TEST", "eyebrow", "sky")
    text(draw, (160, 164), "Dark-field modules work when the canvas carries depth", "title", "white")
    text(draw, (160, 258), "Dark variants preserve text contrast and keep the evidence hierarchy readable", "subtitle", "soft")

    paste_fit(canvas, ASSET / "modules_dark" / "png" / "space_bio_assay_lane.png", (170, 378, 2240, 610))
    paste_fit(canvas, ASSET / "modules_dark" / "png" / "trust_status_ribbon.png", (170, 1245, 2040, 500))
    paste_fit(canvas, ASSET / "modules_dark" / "png" / "claim_boundary_bar.png", (1860, 1740, 1560, 350))

    draw_proof_slot(
        draw,
        canvas,
        (2530, 380, 620, 520),
        "cell_micrograph_texture_field.png",
        "Biology context plate",
        "Generated schematic; source proof still needed",
        dark=True,
    )
    draw_proof_slot(
        draw,
        canvas,
        (3210, 380, 520, 520),
        "single_cell_embedding_field.png",
        "Modality context",
        "Use real embedding for result claims",
        dark=True,
    )

    draw_focus_ring(draw, (170, 378, 2240, 610), "sky")
    draw_badge(draw, (2530, 932), "schematic", "soft", dark=True)
    draw_badge(draw, (2730, 932), "processed", "teal", dark=True)
    draw_badge(draw, (2935, 932), "held-out", "sky", dark=True)

    text(draw, (2530, 1110), "Dark-canvas verdict", "h2", "white")
    y = 1180
    for head, body, tone in [
        ("Pass: contrast", "module labels and lane structure remain readable on deep canvas", "green"),
        ("Pass: hierarchy", "source/proof/status layers do not collapse into decoration", "blue"),
        ("Caution: density", "use one module plus one proof plate; avoid stacking too many controls", "amber"),
    ]:
        draw.rounded_rectangle((2530, y, 3650, y + 98), radius=18, fill=rgba("deep2", 235), outline=rgba(tone, 180), width=2)
        draw.ellipse((2560, y + 34, 2588, y + 62), fill=rgb(tone))
        text(draw, (2618, y + 17), head, "h3", "white")
        text(draw, (2618, y + 54), body, "small", "soft")
        y += 126

    text(draw, (170, 1888), "Design rule", "h2", "white")
    text(draw, (170, 1944), "Layered PNG scene + editable scientific interpretation. Symbols explain status; proof plates carry claims.", "body", "soft")

    return canvas.convert("RGB")


def contact_sheet(light: Image.Image, dark: Image.Image) -> Image.Image:
    width, height = 2400, 1500
    canvas = Image.new("RGB", (width, height), rgb("paper"))
    draw = ImageDraw.Draw(canvas)
    text(draw, (72, 46), "v0.3 symbol/module bridge application A/B", "h2", "ink")
    text(draw, (72, 94), "Light and dark bridge stress tests using reusable biological modules.", "body", "muted")
    thumb_w, thumb_h = 1060, 596
    light_thumb = light.resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
    dark_thumb = dark.resize((thumb_w, thumb_h), Image.Resampling.LANCZOS)
    canvas.paste(light_thumb, (72, 190))
    canvas.paste(dark_thumb, (72 + thumb_w + 110, 190))
    text(draw, (72, 835), "01 light bridge application", "h3", "ink")
    text(draw, (72 + thumb_w + 110, 835), "02 dark bridge application", "h3", "ink")

    notes = [
        ("Strongest use", "sample-to-feature and space-biology assay lane clarify the method faster than generic boxes"),
        ("Safe use", "claim-boundary bar should travel with proof/result slides to prevent overclaiming"),
        ("Caution", "species strip is scope grammar; use as compact marker, not as a main visual thesis"),
        ("v0.4 target", "build source-proof placeholder modules for real screenshots, plots, and source-derived plates"),
    ]
    y = 930
    for idx, (head, body) in enumerate(notes):
        x = 72 + (idx % 2) * 1130
        if idx == 2:
            y += 180
        draw.rounded_rectangle((x, y, x + 1010, y + 130), radius=16, fill=rgb("paper2"), outline=rgb("rule"), width=1)
        text(draw, (x + 28, y + 24), head, "h3", "ink")
        text(draw, (x + 28, y + 68), body, "small", "muted")
    return canvas


def write_manifest() -> None:
    manifest = {
        "created": CREATED,
        "purpose": "Stress-test the biological symbol module v0.3 pack in real bridge-style layouts.",
        "inputs": {
            "symbol_module_pack": "assets/biovis_symbol_module_pack_v0_3",
            "motif_pack": "assets/biovis_motif_pack_v0_3",
        },
        "outputs": {
            "light_bridge": "output/biovis_symbol_module_bridge_application/01_light_bridge_application.png",
            "dark_bridge": "output/biovis_symbol_module_bridge_application/02_dark_bridge_application.png",
            "contact_sheet": "output/biovis_symbol_module_bridge_application/03_light_dark_ab_contact_sheet.png",
        },
        "claim_boundary": "Generated motif plates and symbols are context/status assets, not source proof.",
    }
    (OUT / "biovis_symbol_module_bridge_application_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    qa = {
        "created": CREATED,
        "automatic_checks": {
            "light_png_exists": (OUT / "01_light_bridge_application.png").exists(),
            "dark_png_exists": (OUT / "02_dark_bridge_application.png").exists(),
            "contact_sheet_exists": (OUT / "03_light_dark_ab_contact_sheet.png").exists(),
        },
        "manual_review": {
            "initial_verdict": "pending visual inspection",
            "review_focus": "method comprehension, dark/light contrast, badge density, and claim-boundary clarity",
        },
    }
    (OUT / "biovis_symbol_module_bridge_application_qa.json").write_text(json.dumps(qa, indent=2) + "\n")


def main() -> None:
    ensure()
    light = light_slide()
    dark = dark_slide()
    sheet = contact_sheet(light, dark)
    light.save(OUT / "01_light_bridge_application.png")
    dark.save(OUT / "02_dark_bridge_application.png")
    sheet.save(OUT / "03_light_dark_ab_contact_sheet.png")
    write_manifest()
    print(json.dumps({"output": str(OUT.relative_to(ROOT)), "slides": 2, "contact_sheet": True}, indent=2))


if __name__ == "__main__":
    main()

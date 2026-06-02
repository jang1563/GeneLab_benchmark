#!/usr/bin/env python3
"""Build OSDR source-record proof panels from captured official pages.

The panels combine official OSDR page crops with local proof objects. They are
slide-scale PNG scenes, not final manuscript result figures.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont, ImageOps


ROOT = Path(__file__).resolve().parents[1]
CAPTURE_DIR = ROOT / "output" / "biovis_osdr_source_record_captures_v0_1"
LOCAL_PROOF_DIR = ROOT / "output" / "biovis_source_object_proof_crops_v0_1" / "proof_crops"
OUT = ROOT / "output" / "biovis_osdr_source_record_proof_panel_v0_1"
CROPS = OUT / "crops"
PANELS = OUT / "panels"
QA = OUT / "qa"
CREATED = "2026-06-02"

COLORS = {
    "void": "#080B0F",
    "deep": "#0C1218",
    "deep2": "#101922",
    "deep3": "#172231",
    "ink": "#F4F7F8",
    "soft": "#B7C5D1",
    "muted": "#718293",
    "rule": "#33465A",
    "paper": "#F6F4EF",
    "blue": "#2D6F9F",
    "sky": "#6BAED6",
    "teal": "#1AA090",
    "green": "#178B63",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "red": "#B23A3A",
    "violet": "#6750A4",
}


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
    "eyebrow": font(23, bold=True),
    "title": font(66, bold=True),
    "subtitle": font(32),
    "h1": font(40, bold=True),
    "h2": font(30, bold=True),
    "body": font(25),
    "small": font(19),
    "tiny": font(15),
}


def ensure() -> None:
    for directory in [OUT, CROPS, PANELS, QA]:
        directory.mkdir(parents=True, exist_ok=True)


def load_manifest() -> dict:
    with (CAPTURE_DIR / "osdr_source_record_capture_manifest.json").open() as handle:
        return json.load(handle)


def record_by_id(manifest: dict) -> dict[str, dict]:
    return {record["source_id"]: record for record in manifest["records"]}


def fit_contain(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    copy = image.convert("RGB")
    copy.thumbnail(size, Image.Resampling.LANCZOS)
    return copy


def fit_cover(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    return ImageOps.fit(image.convert("RGB"), size, method=Image.Resampling.LANCZOS, centering=(0.5, 0.5))


def fit_cover_top(image: Image.Image, size: tuple[int, int]) -> Image.Image:
    return ImageOps.fit(image.convert("RGB"), size, method=Image.Resampling.LANCZOS, centering=(0.5, 0.0))


def crop_source_header(source_id: str) -> Path:
    source = Image.open(CAPTURE_DIR / f"{source_id.lower().replace('-', '_')}_viewport.png").convert("RGB")
    # Remove NASA top chrome and the left navigation. Keep the study title,
    # accession fields, description heading, organism, and assay context.
    crop = source.crop((292, 145, 1580, 1120))
    out = CROPS / f"{source_id.lower().replace('-', '_')}_source_record_content_crop.png"
    crop.save(out)
    return out


def make_source_crops(source_ids: list[str]) -> dict[str, Path]:
    return {source_id: crop_source_header(source_id) for source_id in source_ids}


def shadowed_surface(
    canvas: Image.Image,
    image: Image.Image,
    xy: tuple[int, int],
    *,
    width: int,
    height: int,
    radius: int = 30,
    shadow: int = 34,
    outline: tuple[int, int, int, int] | None = None,
) -> None:
    x, y = xy
    surface = fit_cover_top(image, (width, height)).convert("RGBA")
    mask = Image.new("L", (width, height), 0)
    md = ImageDraw.Draw(mask)
    md.rounded_rectangle((0, 0, width, height), radius=radius, fill=255)

    shadow_layer = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow_layer)
    sd.rounded_rectangle((x + 18, y + 26, x + width + 18, y + height + 26), radius=radius, fill=(0, 0, 0, 150))
    shadow_layer = shadow_layer.filter(ImageFilter.GaussianBlur(shadow))
    canvas.alpha_composite(shadow_layer)
    canvas.paste(surface, (x, y), mask)
    if outline is not None:
        draw = ImageDraw.Draw(canvas)
        draw.rounded_rectangle((x, y, x + width, y + height), radius=radius, outline=outline, width=3)


def contained_surface(
    canvas: Image.Image,
    image: Image.Image,
    xy: tuple[int, int],
    *,
    width: int,
    height: int,
    radius: int = 30,
    fill: tuple[int, int, int, int] = rgba("deep2", 240),
    outline: tuple[int, int, int, int] = rgba("rule", 255),
) -> None:
    x, y = xy
    draw = ImageDraw.Draw(canvas)
    shadow_layer = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow_layer)
    sd.rounded_rectangle((x + 18, y + 24, x + width + 18, y + height + 24), radius=radius, fill=(0, 0, 0, 150))
    shadow_layer = shadow_layer.filter(ImageFilter.GaussianBlur(32))
    canvas.alpha_composite(shadow_layer)
    draw.rounded_rectangle((x, y, x + width, y + height), radius=radius, fill=fill, outline=outline, width=2)
    fitted = fit_contain(image, (width - 70, height - 70)).convert("RGBA")
    px = x + (width - fitted.width) // 2
    py = y + (height - fitted.height) // 2
    canvas.alpha_composite(fitted, (px, py))


def draw_background(canvas: Image.Image, *, accent: str = "teal") -> None:
    draw = ImageDraw.Draw(canvas)
    width, height = canvas.size
    draw.rectangle((0, 0, width, height), fill=rgb("void"))
    for y in range(0, height, 2):
        t = y / max(1, height - 1)
        base = tuple(int(rgb("void")[i] * (1 - t) + rgb("deep2")[i] * t) for i in range(3))
        draw.line((0, y, width, y), fill=base + (255,), width=2)
    for x in range(160, width, 160):
        draw.line((x, 0, x, height), fill=rgba("rule", 35), width=1)
    for y in range(140, height, 140):
        draw.line((0, y, width, y), fill=rgba("rule", 28), width=1)

    center = (int(width * 0.68), int(height * 0.36))
    for idx, radius in enumerate([820, 1040, 1260]):
        bbox = (center[0] - radius, center[1] - int(radius * 0.38), center[0] + radius, center[1] + int(radius * 0.38))
        draw.arc(bbox, 196, 348, fill=rgba(accent, 48 - idx * 10), width=3)
    for i in range(42):
        x = int((i * 233) % width)
        y = int((i * 379) % height)
        r = 2 + (i % 3)
        draw.ellipse((x - r, y - r, x + r, y + r), fill=rgba("soft", 34))


def draw_cell_cluster(draw: ImageDraw.ImageDraw, xy: tuple[int, int], *, scale: float = 1.0, tone: str = "teal") -> None:
    x, y = xy
    nodes = [(0, 0, 32), (54, -10, 24), (94, 26, 30), (34, 54, 22), (86, 82, 18), (132, 66, 20)]
    for i, (dx, dy, r) in enumerate(nodes):
        cx = x + int(dx * scale)
        cy = y + int(dy * scale)
        rr = int(r * scale)
        draw.ellipse((cx - rr, cy - rr, cx + rr, cy + rr), outline=rgba(tone, 120), width=max(2, int(3 * scale)))
        draw.ellipse((cx - rr // 3, cy - rr // 3, cx + rr // 3, cy + rr // 3), fill=rgba(tone, 70))
        if i > 0:
            px = x + int(nodes[i - 1][0] * scale)
            py = y + int(nodes[i - 1][1] * scale)
            draw.line((px, py, cx, cy), fill=rgba(tone, 72), width=max(1, int(2 * scale)))


def draw_root_motif(draw: ImageDraw.ImageDraw, xy: tuple[int, int], *, scale: float = 1.0) -> None:
    x, y = xy
    stroke = rgba("green", 118)
    draw.line((x, y, x, y + int(110 * scale)), fill=stroke, width=max(2, int(4 * scale)))
    for idx, (dy, direction) in enumerate([(28, -1), (48, 1), (72, -1), (94, 1)]):
        x2 = x + direction * int((38 + idx * 10) * scale)
        y2 = y + int((dy + 24) * scale)
        draw.line((x, y + int(dy * scale), x2, y2), fill=stroke, width=max(1, int(3 * scale)))
    draw.ellipse((x - int(28 * scale), y - int(34 * scale), x + int(10 * scale), y + int(2 * scale)), outline=rgba("green", 100), width=max(1, int(3 * scale)))
    draw.ellipse((x - int(2 * scale), y - int(38 * scale), x + int(42 * scale), y - int(2 * scale)), outline=rgba("green", 100), width=max(1, int(3 * scale)))


def draw_text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    font_obj,
    fill,
    *,
    width: int | None = None,
    line_height: int = 34,
    max_lines: int = 3,
) -> int:
    x, y = xy
    if width is None:
        draw.text((x, y), text, font=font_obj, fill=fill)
        return y + line_height
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if draw.textlength(candidate, font=font_obj) <= width:
            current = candidate
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1].rstrip(".") + "..."
    for idx, line in enumerate(lines):
        draw.text((x, y + idx * line_height), line, font=font_obj, fill=fill)
    return y + len(lines) * line_height


def pill(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, tone: str, *, fill: str = "deep2") -> None:
    x, y = xy
    width = max(130, int(draw.textlength(text, font=F["small"])) + 54)
    draw.rounded_rectangle((x, y, x + width, y + 40), radius=18, fill=rgba(fill, 232), outline=rgb(tone), width=2)
    draw.ellipse((x + 18, y + 15, x + 28, y + 25), fill=rgb(tone))
    draw.text((x + 40, y + 11), text, font=F["small"], fill=rgb("ink"))


def caption_block(draw: ImageDraw.ImageDraw, xywh: tuple[int, int, int, int], title: str, body: str, *, tone: str = "teal") -> None:
    x, y, w, h = xywh
    draw.rounded_rectangle((x, y, x + w, y + h), radius=24, fill=rgba("deep2", 228), outline=rgba(tone, 150), width=2)
    draw.text((x + 30, y + 28), title, font=F["h2"], fill=rgb("ink"))
    draw_text(draw, (x + 30, y + 76), body, F["small"], rgb("soft"), width=w - 60, line_height=28, max_lines=4)


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], *, tone: str = "teal") -> None:
    x1, y1 = start
    x2, y2 = end
    mx = (x1 + x2) // 2
    points = []
    for i in range(44):
        t = i / 43
        x = int((1 - t) ** 2 * x1 + 2 * (1 - t) * t * mx + t**2 * x2)
        y = int((1 - t) ** 2 * y1 + 2 * (1 - t) * t * (y1 - 80) + t**2 * y2)
        points.append((x, y))
    draw.line(points, fill=rgba(tone, 190), width=5)
    angle = math.atan2(y2 - points[-3][1], x2 - points[-3][0])
    size = 22
    p1 = (x2, y2)
    p2 = (int(x2 - size * math.cos(angle - 0.45)), int(y2 - size * math.sin(angle - 0.45)))
    p3 = (int(x2 - size * math.cos(angle + 0.45)), int(y2 - size * math.sin(angle + 0.45)))
    draw.polygon([p1, p2, p3], fill=rgba(tone, 210))


def footer(draw: ImageDraw.ImageDraw, text: str) -> None:
    draw.rounded_rectangle((120, 1948, 3720, 2040), radius=22, fill=rgba("deep2", 230), outline=rgba("rule", 220), width=2)
    draw.text((158, 1981), "Claim boundary", font=F["h2"], fill=rgb("ink"))
    draw_text(draw, (410, 1984), text, F["small"], rgb("soft"), width=3190, line_height=27, max_lines=2)


def render_contact_sheet(records: dict[str, dict], crops: dict[str, Path]) -> Path:
    canvas = Image.new("RGBA", (3200, 2050), rgba("void", 255))
    draw_background(canvas, accent="sky")
    draw = ImageDraw.Draw(canvas)
    draw.text((120, 96), "Official OSDR source-record captures", font=F["title"], fill=rgb("ink"))
    draw.text((122, 178), "Captured page crops used as source-proof objects before claim slides are composed.", font=F["subtitle"], fill=rgb("soft"))
    placements = {
        "OSD-863": (120, 320, "Human cortical organoid", "teal"),
        "OSD-871": (1640, 320, "Human dopaminergic organoid", "sky"),
        "OSD-120": (120, 1120, "Arabidopsis spaceflight transcriptome", "green"),
        "OSD-918": (1640, 1120, "Mouse single-cell blood data", "amber"),
    }
    for source_id, (x, y, label, tone) in placements.items():
        crop = Image.open(crops[source_id])
        shadowed_surface(canvas, crop, (x, y), width=1380, height=610, radius=30, outline=rgba(tone, 230))
        draw.rounded_rectangle((x + 28, y + 28, x + 282, y + 86), radius=24, fill=rgba("void", 220), outline=rgba(tone, 230), width=2)
        draw.text((x + 54, y + 45), source_id, font=F["h2"], fill=rgb("ink"))
        draw.text((x + 28, y + 644), label, font=F["h2"], fill=rgb("ink"))
        hit_count = sum(1 for value in records[source_id]["expected_term_hits"].values() if value)
        draw.text((x + 28, y + 684), f"{hit_count}/{len(records[source_id]['expected_terms'])} accession terms visible", font=F["small"], fill=rgb("soft"))
        pill(draw, (x + 28, y + 728), "official page", tone)
        pill(draw, (x + 220, y + 728), "reviewed crop", "blue")
    footer(draw, "These screenshots prove official record availability and accession context only; they do not validate local payload completeness or benchmark performance.")
    out = PANELS / "01_osdr_source_record_contact_sheet.png"
    canvas.convert("RGB").save(out)
    return out


def render_organoid_panel(records: dict[str, dict], crops: dict[str, Path]) -> Path:
    canvas = Image.new("RGBA", (3840, 2160), rgba("void", 255))
    draw_background(canvas, accent="teal")
    draw = ImageDraw.Draw(canvas)
    draw_cell_cluster(draw, (3200, 260), scale=2.3, tone="teal")
    draw.text((120, 92), "Human organoid source records", font=F["title"], fill=rgb("ink"))
    draw_text(
        draw,
        (122, 176),
        "Official OSDR pages anchor the cortical and dopaminergic organoid matrices.",
        F["subtitle"],
        rgb("soft"),
        width=1900,
        line_height=42,
        max_lines=2,
    )
    pill(draw, (122, 260), "OSD-863", "teal")
    pill(draw, (304, 260), "OSD-871", "sky")
    pill(draw, (486, 260), "human neural organoids", "violet")

    source_863 = Image.open(crops["OSD-863"])
    source_871 = Image.open(crops["OSD-871"])
    matrix = Image.open(LOCAL_PROOF_DIR / "organoid_matrix_audit_proof.png")
    shadowed_surface(canvas, source_863, (120, 370), width=1480, height=560, radius=34, outline=rgba("teal", 230))
    shadowed_surface(canvas, source_871, (245, 945), width=1480, height=560, radius=34, outline=rgba("sky", 230))
    contained_surface(canvas, matrix, (1810, 425), width=1830, height=1040, radius=38, outline=rgba("teal", 230))

    draw.text((120, 328), "Official source pages", font=F["h2"], fill=rgb("ink"))
    draw.text((1810, 378), "Downloaded matrix proof", font=F["h2"], fill=rgb("ink"))
    arrow(draw, (1600, 648), (1810, 716), tone="teal")
    arrow(draw, (1725, 1224), (1810, 1184), tone="sky")
    caption_block(
        draw,
        (1960, 1534, 1320, 240),
        "What this proves",
        "The studies are official human neural organoid OSDR records with GeneLab IDs, GEO accession, organism, and RNA-seq context. The local matrix crop shows downloaded analysis inputs, not a finalized benchmark claim.",
        tone="teal",
    )
    caption_block(
        draw,
        (120, 1588, 1480, 190),
        "Readable bridge",
        "First-time viewers can follow the path from source page to local matrix without parsing internal file names.",
        tone="sky",
    )
    footer(draw, "Use as a methods or source-evidence slide. Keep the source/date line visible and avoid presenting the organoid pilot as a finished benchmark result.")
    out = PANELS / "02_dark_organoid_source_record_proof.png"
    canvas.convert("RGB").save(out)
    return out


def render_osd120_panel(records: dict[str, dict], crops: dict[str, Path]) -> Path:
    canvas = Image.new("RGBA", (3840, 2160), rgba("void", 255))
    draw_background(canvas, accent="green")
    draw = ImageDraw.Draw(canvas)
    draw_root_motif(draw, (3300, 270), scale=3.0)
    draw.text((120, 92), "OSD-120 source to split proof", font=F["title"], fill=rgb("ink"))
    draw_text(
        draw,
        (122, 176),
        "A plant spaceflight record becomes an auditable local train/test task.",
        F["subtitle"],
        rgb("soft"),
        width=1880,
        line_height=42,
        max_lines=2,
    )
    pill(draw, (122, 260), "Arabidopsis", "green")
    pill(draw, (344, 260), "RNA-seq", "blue")
    pill(draw, (518, 260), "same-study diagnostic", "amber")

    source = Image.open(crops["OSD-120"])
    split = Image.open(LOCAL_PROOF_DIR / "osd120_interaction_split_proof.png")
    shadowed_surface(canvas, source, (120, 400), width=1460, height=980, radius=36, outline=rgba("green", 230))
    contained_surface(canvas, split, (1770, 400), width=1900, height=1080, radius=38, outline=rgba("amber", 230))
    draw.text((120, 350), "Official source page", font=F["h2"], fill=rgb("ink"))
    draw.text((1770, 350), "Local split record", font=F["h2"], fill=rgb("ink"))
    arrow(draw, (1580, 890), (1770, 910), tone="green")
    caption_block(
        draw,
        (120, 1518, 1460, 250),
        "Why this matters",
        "The slide can explain data collection and task construction visually: official OSDR record first, then a local held-out split object with sample counts.",
        tone="green",
    )
    caption_block(
        draw,
        (1770, 1546, 1440, 222),
        "Plain-language caution",
        "This is a same-study task check. It should not be described as leave-one-mission-out or cross-species generalization.",
        tone="amber",
    )
    footer(draw, "This panel supports method explanation and provenance. It does not establish a final generalization benchmark by itself.")
    out = PANELS / "03_dark_osd120_source_record_proof.png"
    canvas.convert("RGB").save(out)
    return out


def render_extension_board(records: dict[str, dict], crops: dict[str, Path]) -> Path:
    canvas = Image.new("RGBA", (3840, 2160), rgba("void", 255))
    draw_background(canvas, accent="amber")
    draw = ImageDraw.Draw(canvas)
    draw_cell_cluster(draw, (3040, 260), scale=1.7, tone="violet")
    draw_root_motif(draw, (3440, 320), scale=2.1)
    draw.text((120, 92), "Extension lanes beyond mouse tissue", font=F["title"], fill=rgb("ink"))
    draw_text(
        draw,
        (122, 176),
        "Real source pages show where organoid, plant, and single-cell expansions can be anchored.",
        F["subtitle"],
        rgb("soft"),
        width=2100,
        line_height=42,
        max_lines=2,
    )

    layout = [
        ("OSD-863", "Human cortical organoid", "teal", 120, 390),
        ("OSD-871", "Human dopaminergic organoid", "sky", 1960, 390),
        ("OSD-120", "Plant transcriptome task", "green", 120, 1185),
        ("OSD-918", "Mouse single-cell expansion", "amber", 1960, 1185),
    ]
    for source_id, label, tone, x, y in layout:
        crop = Image.open(crops[source_id])
        shadowed_surface(canvas, crop, (x, y), width=1640, height=610, radius=34, outline=rgba(tone, 225))
        draw.rounded_rectangle((x + 32, y + 32, x + 286, y + 90), radius=24, fill=rgba("void", 220), outline=rgba(tone, 230), width=2)
        draw.text((x + 58, y + 49), source_id, font=F["h2"], fill=rgb("ink"))
        draw.text((x + 34, y + 646), label, font=F["h2"], fill=rgb("ink"))
        terms = ", ".join(records[source_id]["expected_terms"])
        draw_text(draw, (x + 34, y + 686), f"Visible anchors: {terms}", F["small"], rgb("soft"), width=1420, line_height=28, max_lines=2)

    caption_block(
        draw,
        (122, 282, 1260, 84),
        "Design use",
        "Use this as an extension-map or appendix proof board, not as a scored result figure.",
        tone="amber",
    )
    footer(draw, "Source availability is real; benchmark readiness still requires payload freeze, leakage checks, task definitions, and reviewed metrics.")
    out = PANELS / "04_dark_extension_source_lanes.png"
    canvas.convert("RGB").save(out)
    return out


def grayscale_exports(paths: list[Path]) -> list[Path]:
    exports = []
    for path in paths:
        image = Image.open(path).convert("L")
        out = QA / f"{path.stem}_grayscale.png"
        image.save(out)
        exports.append(out)
    return exports


def luminance(rgb_value: tuple[int, int, int]) -> float:
    srgb = [value / 255 for value in rgb_value]
    linear = [v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4 for v in srgb]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def contrast(c1: tuple[int, int, int], c2: tuple[int, int, int]) -> float:
    l1 = luminance(c1)
    l2 = luminance(c2)
    lighter = max(l1, l2)
    darker = min(l1, l2)
    return (lighter + 0.05) / (darker + 0.05)


def write_metadata(manifest: dict, crops: dict[str, Path], panels: list[Path], grayscale: list[Path]) -> None:
    records = record_by_id(manifest)
    panel_records = [
        {
            "path": rel(path),
            "resolution": Image.open(path).size,
            "role": {
                "01_osdr_source_record_contact_sheet": "review-board/contact-sheet",
                "02_dark_organoid_source_record_proof": "main-deck candidate for organoid source-to-matrix bridge",
                "03_dark_osd120_source_record_proof": "main-deck candidate for OSD-120 source-to-split bridge",
                "04_dark_extension_source_lanes": "appendix or extension-map candidate",
            }.get(path.stem, "proof-panel"),
        }
        for path in panels
    ]
    data = {
        "created": CREATED,
        "input_manifest": rel(CAPTURE_DIR / "osdr_source_record_capture_manifest.json"),
        "source_urls": {source_id: records[source_id]["url"] for source_id in records},
        "source_crops": {source_id: rel(path) for source_id, path in crops.items()},
        "panels": panel_records,
        "grayscale_qa": [rel(path) for path in grayscale],
        "claim_boundary": "Panels prove source-record availability and local proof-object continuity only; they do not validate final benchmark claims.",
    }
    (OUT / "osdr_source_record_proof_panel_manifest.json").write_text(json.dumps(data, indent=2) + "\n")
    qa = {
        "created": CREATED,
        "contrast_checks": {
            "ink_on_deep": round(contrast(rgb("ink"), rgb("deep")), 2),
            "soft_on_deep": round(contrast(rgb("soft"), rgb("deep")), 2),
            "ink_on_void": round(contrast(rgb("ink"), rgb("void")), 2),
        },
        "visual_checks": [
            "Raw OSDR browser chrome is cropped away from main panels.",
            "Source page, local proof object, and claim boundary occupy separate z-layers.",
            "Color is reinforced by text labels; green/amber/teal are not sole carriers of meaning.",
            "Panels are slide-scale PNG scenes with editable interpretation expected in the final deck.",
        ],
        "review_required": [
            "Before manuscript use, replace or accompany screenshots with permanent source citations.",
            "Do not use the extension-lane board as evidence that all lanes are benchmark-ready.",
        ],
    }
    (OUT / "osdr_source_record_proof_panel_qa.json").write_text(json.dumps(qa, indent=2) + "\n")


def main() -> None:
    ensure()
    manifest = load_manifest()
    records = record_by_id(manifest)
    source_ids = ["OSD-863", "OSD-871", "OSD-120", "OSD-918"]
    missing = [source_id for source_id in source_ids if source_id not in records]
    if missing:
        raise SystemExit(f"Missing capture records: {missing}")
    crops = make_source_crops(source_ids)
    panels = [
        render_contact_sheet(records, crops),
        render_organoid_panel(records, crops),
        render_osd120_panel(records, crops),
        render_extension_board(records, crops),
    ]
    grayscale = grayscale_exports(panels)
    write_metadata(manifest, crops, panels, grayscale)
    print(
        json.dumps(
            {
                "output": rel(OUT),
                "panels": [rel(path) for path in panels],
                "grayscale_qa": [rel(path) for path in grayscale],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

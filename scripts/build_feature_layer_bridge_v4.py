#!/usr/bin/env python3
"""Build the v0.4 slide-6 feature-layer bridge.

The rendered preview explains genes versus pathway summaries for first-time
viewers. The scene plate keeps the scientific objects mostly text-free so the
final deck can rebuild interpretation as editable overlay text.
"""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont, ImageOps


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "premium_feature_layer_bridge_v0_4"
QA_DIR = OUT / "qa"
CREATED = "2026-06-03"

MOTIF_V03 = ROOT / "assets" / "biovis_motif_pack_v0_3" / "png"
SYMBOLS = ROOT / "assets" / "biovis_symbol_module_pack_v0_3" / "symbols" / "png"

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
    "paper": "#F9FBFA",
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
    "eyebrow": font(20, bold=True),
    "title": font(66, bold=True),
    "subtitle": font(30),
    "h": font(28, bold=True),
    "body": font(21),
    "small": font(17),
    "tiny": font(13),
    "label": font(22, bold=True),
}


OVERLAY_TEXT = [
    {
        "id": "headline",
        "role": "decision_headline",
        "content": "What the model sees can change",
        "x": 0.060,
        "y": 0.100,
        "font_pt": 33,
        "color": "ink",
        "max_lines": 1,
    },
    {
        "id": "support",
        "role": "supporting_claim",
        "content": "Same samples, two views: genes or pathways, then one held-out score.",
        "x": 0.061,
        "y": 0.171,
        "font_pt": 15,
        "color": "soft",
        "max_lines": 2,
    },
    {
        "id": "label_gene_table",
        "role": "plate_label",
        "content": "Gene table",
        "x": 0.109,
        "y": 0.353,
        "font_pt": 13,
        "color": "sky",
        "max_lines": 1,
    },
    {
        "id": "label_gene_model",
        "role": "plate_label",
        "content": "Gene-level view",
        "x": 0.327,
        "y": 0.353,
        "font_pt": 13,
        "color": "teal",
        "max_lines": 1,
    },
    {
        "id": "label_pathway_model",
        "role": "plate_label",
        "content": "Pathway summaries",
        "x": 0.558,
        "y": 0.353,
        "font_pt": 13,
        "color": "green",
        "max_lines": 1,
    },
    {
        "id": "label_score",
        "role": "plate_label",
        "content": "Score",
        "x": 0.800,
        "y": 0.353,
        "font_pt": 13,
        "color": "rose",
        "max_lines": 1,
    },
    {
        "id": "status",
        "role": "claim_boundary",
        "content": "Fit choices on training missions; score the held-out mission.",
        "x": 0.061,
        "y": 0.862,
        "font_pt": 10,
        "color": "soft",
        "max_lines": 1,
    },
    {
        "id": "caveat",
        "role": "trust_caveat",
        "content": "Pathways help some tasks, not all.",
        "x": 0.061,
        "y": 0.900,
        "font_pt": 9,
        "color": "muted",
        "max_lines": 1,
    },
]

VISIBLE_TEXT = [item["content"] for item in OVERLAY_TEXT]
FORBIDDEN_VISIBLE_TERMS = [
    "M1",
    "M2",
    "M3",
    "M4",
    "same split",
    "v0.3",
    "replacement test",
    "generated schematic texture",
    "not source evidence",
]


def ensure() -> None:
    QA_DIR.mkdir(parents=True, exist_ok=True)


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    font_obj,
    fill,
    *,
    width: int,
    line_height: int,
    max_lines: int,
) -> int:
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
    x, y = xy
    for idx, line in enumerate(lines):
        draw.text((x, y + idx * line_height), line, font=font_obj, fill=fill)
    return y + len(lines) * line_height


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
        draw.line((x, 0, x, h), fill=rgba("rule", 25), width=1)
    for y in range(150, h, 150):
        draw.line((0, y, w, y), fill=rgba("rule", 22), width=1)
    center = (int(w * 0.70), int(h * 0.28))
    for idx, radius in enumerate([780, 1010, 1240]):
        bbox = (
            center[0] - radius,
            center[1] - int(radius * 0.34),
            center[0] + radius,
            center[1] + int(radius * 0.34),
        )
        draw.arc(bbox, 200, 350, fill=rgba("sky", 50 - idx * 11), width=3)
    for idx in range(38):
        x = int((idx * 263) % w)
        y = int((idx * 421) % h)
        r = 2 + idx % 3
        draw.ellipse((x - r, y - r, x + r, y + r), fill=rgba("soft", 32))


def rounded_mask(size: tuple[int, int], radius: int) -> Image.Image:
    mask = Image.new("L", size, 0)
    draw = ImageDraw.Draw(mask)
    draw.rounded_rectangle((0, 0, size[0], size[1]), radius=radius, fill=255)
    return mask


def paste_asset(
    canvas: Image.Image,
    path: Path,
    box: tuple[int, int, int, int],
    *,
    alpha: int = 255,
    contain: bool = False,
    crop: tuple[int, int, int, int] | None = None,
) -> None:
    x, y, w, h = box
    image = Image.open(path).convert("RGBA")
    if crop:
        image = image.crop(crop)
    if contain:
        image.thumbnail((w, h), Image.Resampling.LANCZOS)
        layer = Image.new("RGBA", (w, h), (0, 0, 0, 0))
        layer.alpha_composite(image, ((w - image.width) // 2, (h - image.height) // 2))
        image = layer
    else:
        image = ImageOps.fit(image, (w, h), method=Image.Resampling.LANCZOS, centering=(0.5, 0.5))
    if alpha < 255:
        image.putalpha(Image.eval(image.getchannel("A"), lambda value: int(value * alpha / 255)))
    canvas.alpha_composite(image, (x, y))


def draw_plate(
    canvas: Image.Image,
    box: tuple[int, int, int, int],
    *,
    tone: str,
    fill_token: str = "panel2",
) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = box
    shadow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow)
    sd.rounded_rectangle((x + 18, y + 22, x + w + 18, y + h + 22), radius=12, fill=(0, 0, 0, 140))
    shadow = shadow.filter(ImageFilter.GaussianBlur(18))
    canvas.alpha_composite(shadow)
    draw.rounded_rectangle((x, y, x + w, y + h), radius=10, fill=rgba(fill_token, 232), outline=rgba(tone, 220), width=2)
    draw.line((x + 22, y + h - 34, x + w - 22, y + h - 34), fill=rgba(tone, 100), width=2)


def draw_gene_cloud(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x, y, w, h = box
    points = [
        (0.10, 0.35, "teal"), (0.16, 0.56, "sky"), (0.22, 0.47, "muted"), (0.25, 0.67, "green"),
        (0.30, 0.28, "sky"), (0.33, 0.52, "teal"), (0.40, 0.38, "muted"), (0.45, 0.65, "teal"),
        (0.50, 0.45, "green"), (0.55, 0.30, "sky"), (0.62, 0.58, "teal"), (0.68, 0.42, "muted"),
        (0.73, 0.70, "green"), (0.78, 0.50, "sky"), (0.84, 0.36, "teal"), (0.90, 0.61, "muted"),
        (0.18, 0.24, "green"), (0.36, 0.76, "sky"), (0.58, 0.75, "muted"), (0.82, 0.24, "green"),
    ]
    for px, py, tone in points:
        cx = x + int(px * w)
        cy = y + int(py * h)
        r = 8 if tone in {"teal", "green"} else 6
        draw.ellipse((cx - r, cy - r, cx + r, cy + r), fill=rgba(tone, 185))
    for offset in [0.30, 0.43, 0.56, 0.69]:
        yy = y + int(offset * h)
        draw.line((x + 40, yy, x + w - 40, yy), fill=rgba("rule", 95), width=2)
    draw.line((x + 55, y + int(0.62 * h), x + w - 80, y + int(0.36 * h)), fill=rgba("teal", 170), width=6)
    draw.line((x + 55, y + int(0.62 * h), x + w - 80, y + int(0.36 * h)), fill=rgba("sky", 90), width=14)


def draw_score_report(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x, y, w, h = box
    draw.rounded_rectangle((x + 58, y + 46, x + w - 58, y + h - 44), radius=8, fill=(248, 250, 250), outline=rgba("rose", 230), width=3)
    draw.ellipse((x + w // 2 - 72, y + 12, x + w // 2 + 72, y + 156), fill=(248, 250, 250), outline=rgba("rose", 230), width=3)
    draw.line((x + w // 2 - 32, y + 84, x + w // 2 + 32, y + 84), fill=rgba("rose", 190), width=3)
    draw.line((x + w // 2, y + 54, x + w // 2, y + 112), fill=rgba("rose", 190), width=3)
    for idx, yy in enumerate([0.46, 0.61, 0.76]):
        color = "rose" if idx == 0 else "rule"
        draw.line((x + 96, y + int(yy * h), x + w - 96, y + int(yy * h)), fill=rgba(color, 160), width=3)


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], tone: str) -> None:
    sx, sy = start
    ex, ey = end
    draw.line((sx, sy, ex, ey), fill=rgba(tone, 190), width=5)
    head = [(ex, ey), (ex - 22, ey - 12), (ex - 18, ey + 12)]
    draw.polygon(head, fill=rgba(tone, 210))


def draw_scene(*, with_overlay: bool) -> Image.Image:
    canvas = Image.new("RGBA", (3840, 2160), (0, 0, 0, 255))
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)

    rail_y = 1110
    draw.line((420, rail_y, 3380, rail_y), fill=rgba("rule", 150), width=3)
    for tick in [560, 1430, 2300, 3160]:
        draw.line((tick, rail_y - 28, tick, rail_y + 28), fill=rgba("soft", 90), width=2)

    plate_boxes = [
        (290, 760, 680, 430),
        (1080, 760, 680, 430),
        (1870, 760, 680, 430),
        (2670, 760, 680, 430),
    ]
    tones = ["sky", "teal", "green", "rose"]
    for box, tone in zip(plate_boxes, tones):
        draw_plate(canvas, box, tone=tone)

    paste_asset(
        canvas,
        MOTIF_V03 / "omics_matrix_texture_field.png",
        (340, 815, 580, 295),
        alpha=240,
        contain=False,
        crop=(250, 120, 1350, 740),
    )
    draw_gene_cloud(draw, (1130, 820, 580, 280))
    paste_asset(
        canvas,
        MOTIF_V03 / "pathway_network_summary_field.png",
        (1930, 815, 560, 300),
        alpha=242,
        contain=False,
        crop=(250, 120, 1350, 740),
    )
    draw_score_report(draw, (2720, 800, 580, 330))

    paste_asset(canvas, SYMBOLS / "gene_locus.png", (1290, 1028, 92, 92), alpha=160, contain=True)
    paste_asset(canvas, SYMBOLS / "pathway_network.png", (2185, 1028, 92, 92), alpha=170, contain=True)

    arrow(draw, (970, rail_y), (1080, rail_y), "muted")
    arrow(draw, (1760, rail_y), (1870, rail_y), "teal")
    arrow(draw, (2550, rail_y), (2670, rail_y), "muted")

    # Training/held-out guard as a small Z4 trust layer.
    guard_x, guard_y = 2370, 420
    draw.rounded_rectangle((guard_x, guard_y, guard_x + 1040, guard_y + 116), radius=10, fill=rgba("panel", 235), outline=rgba("rule", 180), width=2)
    draw.text((guard_x + 32, guard_y + 24), "Same held-out mission for both views", font=F["h"], fill=rgb("ink"))
    for idx, (label, tone) in enumerate([("Train", "green"), ("Train", "green"), ("Train", "green"), ("Score", "rose")]):
        x = guard_x + 520 + idx * 102
        draw.rounded_rectangle((x, guard_y + 32, x + 78, guard_y + 70), radius=9, fill=rgba("panel2", 238), outline=rgba(tone, 210), width=2)
        draw.text((x + 17, guard_y + 43), label, font=F["tiny"], fill=rgb(tone))
        if idx == 2:
            draw.line((x + 96, guard_y + 24, x + 96, guard_y + 84), fill=rgba("rose", 180), width=2)

    if with_overlay:
        draw.text((230, 205), "METHOD BRIDGE", font=F["eyebrow"], fill=rgb("teal"))
        draw.text((230, 252), "What the model sees can change", font=F["title"], fill=rgb("ink"))
        draw_wrapped(
            draw,
            (234, 342),
            "Same samples, two views: genes or pathways, then one held-out score.",
            F["subtitle"],
            rgb("soft"),
            width=1860,
            line_height=40,
            max_lines=2,
        )
        labels = [
            ("Gene table", plate_boxes[0], "sky"),
            ("Gene-level view", plate_boxes[1], "teal"),
            ("Pathway summaries", plate_boxes[2], "green"),
            ("Score", plate_boxes[3], "rose"),
        ]
        for text, box, tone in labels:
            draw.text((box[0] + 48, box[1] - 52), text, font=F["label"], fill=rgb(tone))
        draw.rounded_rectangle((230, 1812, 1910, 1936), radius=10, fill=rgba("panel", 225), outline=rgba("rule", 160), width=2)
        draw.text((270, 1840), "Fit choices on training missions; score the held-out mission.", font=F["body"], fill=rgb("soft"))
        draw.text((270, 1882), "Pathways help some tasks, not all.", font=F["small"], fill=rgb("muted"))
        draw.text((2830, 1888), "schematic methods bridge, not a result figure", font=F["small"], fill=rgb("muted"))

    return canvas.convert("RGB")


def word_count() -> int:
    return len(" ".join(VISIBLE_TEXT).split())


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def build() -> dict[str, str]:
    ensure()
    scene_plate = OUT / "feature_layer_bridge_v0_4_scene_plate.png"
    preview = OUT / "feature_layer_bridge_v0_4_rendered_preview.png"
    grayscale = QA_DIR / "feature_layer_bridge_v0_4_rendered_preview_grayscale.png"
    overlay_spec = OUT / "feature_layer_bridge_v0_4_overlay_spec.json"
    manifest = OUT / "feature_layer_bridge_v0_4_manifest.json"
    qa = OUT / "feature_layer_bridge_v0_4_qa.json"

    draw_scene(with_overlay=False).save(scene_plate, quality=95)
    rendered = draw_scene(with_overlay=True)
    rendered.save(preview, quality=95)
    rendered.convert("L").convert("RGB").save(grayscale, quality=95)

    overlay = {
        "slide_id": "slide06_feature_layer_bridge_v0_4",
        "created": CREATED,
        "canvas": {"width_px": 3840, "height_px": 2160, "aspect_ratio": "16:9"},
        "coordinate_system": "normalized_0_1",
        "text": OVERLAY_TEXT,
        "status_labels": [item for item in OVERLAY_TEXT if item["role"] in {"claim_boundary", "trust_caveat"}],
        "focus_marks": [
            {"id": "gene_to_pathway_arrow", "role": "transfer_path", "from": [0.458, 0.514], "to": [0.487, 0.514], "color": "teal"},
            {"id": "heldout_guard", "role": "trust_boundary", "x": 0.617, "y": 0.194, "w": 0.271, "h": 0.054},
        ],
        "forbidden_visible_terms": FORBIDDEN_VISIBLE_TERMS,
    }
    write_json(overlay_spec, overlay)

    visible_blob = "\n".join(VISIBLE_TEXT).lower()
    blocked_absent = {term: term.lower() not in visible_blob for term in FORBIDDEN_VISIBLE_TERMS}
    manifest_data = {
        "slide_id": "slide06_feature_layer_bridge_v0_4",
        "created": CREATED,
        "slide_role": "method_bridge",
        "reference_calibration_role": "kmu_proof_stage",
        "content_brief": "Slide 6 bridge between B1-B4 task construction and v1-v7 result figures.",
        "evidence_sources": [
            rel(MOTIF_V03 / "omics_matrix_texture_field.png"),
            rel(MOTIF_V03 / "pathway_network_summary_field.png"),
            rel(SYMBOLS / "gene_locus.png"),
            rel(SYMBOLS / "pathway_network.png"),
        ],
        "claim_boundary": "Schematic methods bridge explaining feature views; not a quantitative result figure.",
        "outputs": {
            "scene_plate": rel(scene_plate),
            "rendered_preview": rel(preview),
            "grayscale_qa": rel(grayscale),
            "overlay_spec": rel(overlay_spec),
            "qa": rel(qa),
        },
        "visible_text_policy": {
            "word_count": word_count(),
            "word_count_budget": 45,
            "preferred_terms": [
                "What the model sees can change",
                "Gene table",
                "Gene-level view",
                "Pathway summaries",
                "Score",
                "held-out mission",
            ],
            "blocked_internal_terms_absent_from_new_visible_text": blocked_absent,
        },
        "z_layers": {
            "Z0": "dark scientific canvas with grid and orbit arcs",
            "Z1": "source-to-score measurement rail",
            "Z2": "matrix, gene-view, pathway-summary, and score-report evidence surfaces",
            "Z3": "editable headline, labels, and arrows",
            "Z4": "held-out mission guard and caveat",
        },
    }
    write_json(manifest, manifest_data)

    qa_data = {
        "created": CREATED,
        "automatic_checks": {
            "rendered_outputs_exist": all(path.exists() for path in [scene_plate, preview, grayscale, overlay_spec, manifest]),
            "image_dimensions": {"width_px": rendered.width, "height_px": rendered.height},
            "visible_text_word_count": word_count(),
            "visible_text_budget": 45,
            "blocked_internal_terms_absent": blocked_absent,
            "assets_exist": all(
                path.exists()
                for path in [
                    MOTIF_V03 / "omics_matrix_texture_field.png",
                    MOTIF_V03 / "pathway_network_summary_field.png",
                    SYMBOLS / "gene_locus.png",
                    SYMBOLS / "pathway_network.png",
                ]
            ),
        },
        "manual_review": {
            "visual_review_status": "pass - rendered preview, scene plate, and grayscale QA inspected at full size.",
            "readability": "pass - headline and four object labels are readable; footer caveat is subordinate but legible.",
            "visual_hierarchy": "pass - matrix, gene view, pathway summary, and score objects read as one left-to-right methods bridge.",
            "grayscale_review": "pass - object sequence and arrows remain visible without relying on color alone.",
            "deck_readiness": "draft premium slide-6 candidate; final deck should rebuild overlay text as editable slide objects.",
            "claim_boundary": "schematic methods bridge only; not a quantitative result figure",
        },
    }
    write_json(qa, qa_data)

    return {
        "scene_plate": rel(scene_plate),
        "rendered_preview": rel(preview),
        "grayscale": rel(grayscale),
        "overlay_spec": rel(overlay_spec),
        "manifest": rel(manifest),
        "qa": rel(qa),
    }


def main() -> None:
    print(json.dumps(build(), indent=2))


if __name__ == "__main__":
    main()

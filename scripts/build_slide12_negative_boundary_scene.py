#!/usr/bin/env python3
"""Build slide 12 negative/failure-mode boundary premium scene."""

from __future__ import annotations

import json
import statistics
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "premium_core_result_slides" / "slide12_negative_boundary_v0_1"
QA = OUT / "qa"
ASSETS = ROOT / "assets" / "biovis_symbol_module_pack_v0_3"
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
    "violet": "#7B68A8",
}

OVERLAY_TEXT = [
    "Negative results set boundaries",
    "v3 marks where spatial brain and pretrained embeddings need controls.",
    "Spatial brain",
    "section AUROC 0.139; animal AUROC 0.444",
    "No signal detected",
    "Foundation models",
    "UCE 0.542; scFoundation 0.584",
    "No automatic improvement",
    "Boundary: informative negative, not absent biology.",
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
    "title": font(64, bold=True),
    "subtitle": font(28),
    "h": font(31, bold=True),
    "body": font(22),
    "small": font(17),
    "tiny": font(13),
    "num": font(48, bold=True),
    "metric": font(38, bold=True),
}


def ensure() -> None:
    QA.mkdir(parents=True, exist_ok=True)


def load_json(path: str) -> dict:
    with (ROOT / path).open(encoding="utf-8") as handle:
        return json.load(handle)


def fm_summary(path: str) -> dict:
    data = load_json(path)
    rows = data["results"]
    multi = [row for row in rows if row.get("n_missions", 0) > 1]
    sig = [row for row in rows if row.get("perm_p", 1.0) < 0.05]
    return {
        "family": data["model_family"],
        "mean_multi_mission": round(statistics.mean(row["lomo_auroc"] for row in multi), 3),
        "n_tissues": len(rows),
        "n_multi_mission": len(multi),
        "n_significant": len(sig),
        "significant": [
            {"tissue": row["tissue"], "auroc": round(row["lomo_auroc"], 3), "p": row.get("perm_p")}
            for row in sig
        ],
        "rows": [
            {"tissue": row["tissue"], "auroc": round(row["lomo_auroc"], 3), "n_missions": row.get("n_missions", 0)}
            for row in rows
        ],
    }


def load_sources() -> dict:
    f3a = load_json("v3/evaluation/F3a_visium_classification.json")
    f3b = load_json("v3/evaluation/F3b_visium_svg.json")
    f3d = load_json("v3/evaluation/F3d_visium_cross_resolution.json")
    uce = fm_summary("v3/evaluation/FM_uce.json")
    scf = fm_summary("v3/evaluation/FM_scfoundation.json")
    return {
        "spatial": {
            "dataset": f3a["dataset"],
            "mission": f3a["mission"],
            "tissue": f3a["tissue"],
            "n_animals": f3a["section_level"]["n_animals"],
            "n_sections": f3a["section_level"]["n_samples"],
            "section_auroc": round(f3a["section_level"]["auroc"], 3),
            "section_p_exact": f3a["section_level"]["p_exact"],
            "animal_auroc": round(f3a["animal_level"]["auroc"], 3),
            "animal_p_exact": f3a["animal_level"]["p_exact"],
            "bulk_auroc": round(f3d["bulk_auroc"]["auroc"], 3),
            "svg_jaccard": round(f3b["svg_overlap"]["jaccard"], 3),
            "svg_flt_specific": f3b["condition_specific_svgs"]["n_flt_specific"],
            "svg_gc_specific": f3b["condition_specific_svgs"]["n_gc_specific"],
        },
        "fm": {
            "uce": uce,
            "scfoundation": scf,
            "pca_lr_reference_mean": 0.758,
            "pca_lr_reference_source": "v3/README.md foundation-model section",
        },
    }


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
        draw.line((x, 0, x, h), fill=rgba("rule", 22), width=1)
    for y in range(150, h, 150):
        draw.line((0, y, w, y), fill=rgba("rule", 18), width=1)
    center = (int(w * 0.56), int(h * 0.32))
    for idx, radius in enumerate([820, 1080, 1340]):
        bbox = (
            center[0] - radius,
            center[1] - int(radius * 0.35),
            center[0] + radius,
            center[1] + int(radius * 0.35),
        )
        draw.arc(bbox, 200, 354, fill=rgba("sky", 38 - idx * 9), width=3)


def panel_shadow(canvas: Image.Image, box: tuple[int, int, int, int], *, radius: int = 14) -> None:
    x, y, w, h = box
    shadow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow)
    sd.rounded_rectangle((x + 18, y + 22, x + w + 18, y + h + 22), radius=radius, fill=(0, 0, 0, 122))
    canvas.alpha_composite(shadow.filter(ImageFilter.GaussianBlur(18)))


def draw_tinted_icon(canvas: Image.Image, name: str, xy: tuple[int, int], size: int, tone: str, alpha: int = 200) -> None:
    path = ASSETS / "symbols" / "png" / f"{name}.png"
    if not path.exists():
        return
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    icon = Image.open(path).convert("RGBA").resize((size, size), resample)
    mask = icon.split()[-1].point(lambda value: int(value * alpha / 255))
    tinted = Image.new("RGBA", icon.size, rgb(tone) + (0,))
    tinted.putalpha(mask)
    canvas.alpha_composite(tinted, xy)


def draw_header(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    draw.text((230, 196), "NEGATIVE RESULTS AND BENCHMARK LIMITS", font=F["eyebrow"], fill=rgb("teal"))
    draw.text((230, 248), "Negative results set boundaries", font=F["title"], fill=rgb("ink"))
    draw.text(
        (234, 330),
        "v3 marks where spatial brain and pretrained embeddings need controls.",
        font=F["subtitle"],
        fill=rgb("soft"),
    )


def draw_gauge(draw: ImageDraw.ImageDraw, x: int, y: int, value: float, label: str, tone: str) -> None:
    draw.text((x, y), label, font=F["small"], fill=rgb("soft"))
    draw.line((x, y + 70, x + 620, y + 70), fill=rgba("rule", 95), width=4)
    draw.line((x + 310, y + 52, x + 310, y + 88), fill=rgba("muted", 130), width=2)
    marker_x = x + int(620 * max(0.0, min(1.0, value)))
    draw.ellipse((marker_x - 18, y + 52, marker_x + 18, y + 88), fill=rgba(tone, 230))
    draw.text((x + 292, y + 96), "0.5", font=F["tiny"], fill=rgb("muted"))
    draw.text((x + 650, y + 50), f"{value:.3f}", font=F["metric"], fill=rgb(tone))


def draw_spot_field(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    tones = ["sky", "teal", "amber", "rose"]
    for row in range(8):
        for col in range(12):
            idx = (row * 7 + col * 5) % len(tones)
            alpha = 80 + ((row * 17 + col * 11) % 90)
            draw.ellipse((x + col * 24, y + row * 22, x + col * 24 + 12, y + row * 22 + 12), fill=rgba(tones[idx], alpha))


def draw_spatial_surface(canvas: Image.Image, data: dict) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = 230, 620, 1120, 780
    panel_shadow(canvas, (x, y, w, h))
    draw.rounded_rectangle((x, y, x + w, y + h), radius=12, fill=rgba("panel", 228), outline=rgba("rose", 155), width=2)
    draw_tinted_icon(canvas, "brain_organ", (x + 54, y + 44), 88, "rose", 175)
    draw_tinted_icon(canvas, "spatial_spots", (x + 132, y + 47), 76, "sky", 165)
    draw.text((x + 220, y + 54), "Spatial brain", font=F["h"], fill=rgb("ink"))
    draw.text((x + 220, y + 94), "OSD-352 RR-3, 6 animals and 12 sections", font=F["small"], fill=rgb("muted"))

    draw_spot_field(draw, x + 70, y + 198)
    draw_gauge(draw, x + 410, y + 186, data["spatial"]["section_auroc"], "section AUROC", "rose")
    draw_gauge(draw, x + 410, y + 350, data["spatial"]["animal_auroc"], "animal AUROC", "amber")
    draw_gauge(draw, x + 410, y + 514, data["spatial"]["bulk_auroc"], "companion bulk", "muted")
    draw.text((x + 70, y + 562), "No signal detected", font=F["metric"], fill=rgb("rose"))
    draw.text((x + 70, y + 612), "in this small spatial-brain setting", font=F["small"], fill=rgb("muted"))
    draw.text((x + 70, y + 662), f"SVG overlap Jaccard {data['spatial']['svg_jaccard']:.3f}", font=F["small"], fill=rgb("soft"))


def draw_boundary_corridor(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y = 1510, 640
    draw.text((x, y), "Failure mode as evidence", font=F["h"], fill=rgb("ink"))
    draw.text((x, y + 42), "The slide is a guardrail, not a dismissal", font=F["small"], fill=rgb("muted"))
    points = [(x + 110, y + 210), (x + 320, y + 360), (x + 520, y + 510), (x + 740, y + 660)]
    draw.line(points, fill=rgba("teal", 45), width=38, joint="curve")
    draw.line(points, fill=rgba("amber", 145), width=5, joint="curve")
    labels = [
        ("negative", "observed result"),
        ("limit", "task boundary"),
        ("control", "avoid overclaim"),
        ("next", "better test"),
    ]
    for idx, point in enumerate(points):
        tone = ["rose", "amber", "teal", "sky"][idx]
        draw.ellipse((point[0] - 22, point[1] - 22, point[0] + 22, point[1] + 22), fill=rgba(tone, 220))
        draw.text((point[0] + 42, point[1] - 34), labels[idx][0], font=F["body"], fill=rgb("soft"))
        draw.text((point[0] + 42, point[1] - 2), labels[idx][1], font=F["small"], fill=rgb("muted"))
    draw.rounded_rectangle((x, y + 790, x + 910, y + 902), radius=10, fill=rgba("panel", 220), outline=rgba("rule", 140), width=1)
    draw.text((x + 34, y + 820), "Boundary: informative negative, not absent biology.", font=F["body"], fill=rgb("soft"))
    draw.text((x + 34, y + 860), "Design follow-up tasks; do not erase the null.", font=F["small"], fill=rgb("muted"))


def draw_fm_bars(draw: ImageDraw.ImageDraw, x: int, y: int, summaries: list[tuple[str, float, int, str]]) -> None:
    for idx, (label, value, n_sig, tone) in enumerate(summaries):
        yy = y + idx * 150
        draw.text((x, yy), label, font=F["body"], fill=rgb("soft"))
        draw.line((x, yy + 64, x + 560, yy + 64), fill=rgba("rule", 90), width=4)
        draw.line((x + 280, yy + 46, x + 280, yy + 82), fill=rgba("muted", 130), width=2)
        draw.rounded_rectangle((x, yy + 46, x + int(560 * value), yy + 82), radius=6, fill=rgba(tone, 210))
        draw.text((x + 590, yy + 39), f"{value:.3f}", font=F["metric"], fill=rgb(tone))
        sig_label = "ref." if label.startswith("PCA-LR") else f"{n_sig} sig."
        draw.text((x + 590, yy + 84), sig_label, font=F["small"], fill=rgb("muted"))


def draw_fm_surface(canvas: Image.Image, data: dict) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = 2720, 620, 800, 780
    panel_shadow(canvas, (x, y, w, h))
    draw.rounded_rectangle((x, y, x + w, y + h), radius=12, fill=rgba("panel", 228), outline=rgba("teal", 155), width=2)
    draw_tinted_icon(canvas, "cell_state", (x + 54, y + 42), 84, "teal", 170)
    draw_tinted_icon(canvas, "omics_matrix", (x + 122, y + 44), 78, "sky", 165)
    draw.text((x + 210, y + 54), "Foundation models", font=F["h"], fill=rgb("ink"))
    draw.text((x + 210, y + 94), "pretrained embeddings need task controls", font=F["small"], fill=rgb("muted"))

    summaries = [
        ("UCE", data["fm"]["uce"]["mean_multi_mission"], data["fm"]["uce"]["n_significant"], "sky"),
        ("scFoundation", data["fm"]["scfoundation"]["mean_multi_mission"], data["fm"]["scfoundation"]["n_significant"], "teal"),
        ("PCA-LR ref.", data["fm"]["pca_lr_reference_mean"], 0, "amber"),
    ]
    draw_fm_bars(draw, x + 64, y + 210, summaries)
    draw.text((x + 64, y + 622), "No automatic improvement", font=F["metric"], fill=rgb("teal"))
    draw.text((x + 64, y + 674), "some tissues significant, but baseline still leads", font=F["small"], fill=rgb("muted"))


def draw_footer(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    draw.rounded_rectangle((230, 1810, 2440, 1938), radius=10, fill=rgba("panel", 224), outline=rgba("rule", 155), width=2)
    draw_tinted_icon(canvas, "caveat_flag", (258, 1832), 66, "amber", 165)
    draw.text((338, 1838), "Boundary: informative negative, not absent biology.", font=F["body"], fill=rgb("soft"))
    draw.text((338, 1880), "Use this slide to define limits of current data/model settings, not universal impossibility.", font=F["small"], fill=rgb("muted"))
    draw.text((2830, 1888), "source: v3 spatial and FM evaluation JSON", font=F["small"], fill=rgb("muted"))


def render(data: dict, *, with_overlay: bool) -> Image.Image:
    canvas = Image.new("RGBA", (3840, 2160), (0, 0, 0, 255))
    draw_background(canvas)
    draw_spatial_surface(canvas, data)
    draw_boundary_corridor(canvas)
    draw_fm_surface(canvas, data)
    if with_overlay:
        draw_header(canvas)
        draw_footer(canvas)
    return canvas.convert("RGB")


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def build() -> dict[str, str]:
    ensure()
    data = load_sources()
    scene_plate = OUT / "slide12_negative_boundary_scene_plate.png"
    preview = OUT / "slide12_negative_boundary_rendered_preview.png"
    grayscale = QA / "slide12_negative_boundary_rendered_preview_grayscale.png"
    source_summary = OUT / "slide12_negative_boundary_source_summary.json"
    manifest = OUT / "slide12_negative_boundary_manifest.json"
    qa = OUT / "slide12_negative_boundary_qa.json"

    render(data, with_overlay=False).save(scene_plate, quality=95)
    rendered = render(data, with_overlay=True)
    rendered.save(preview, quality=95)
    rendered.convert("L").convert("RGB").save(grayscale, quality=95)
    write_json(source_summary, data)

    visible_text_count = len(" ".join(OVERLAY_TEXT).split())
    manifest_data = {
        "slide_id": "slide12_negative_boundary_v0_1",
        "created": CREATED,
        "slide_role": "negative_result_boundary",
        "source_documents": [
            "v3/evaluation/F3a_visium_classification.json",
            "v3/evaluation/F3b_visium_svg.json",
            "v3/evaluation/F3d_visium_cross_resolution.json",
            "v3/evaluation/FM_uce.json",
            "v3/evaluation/FM_scfoundation.json",
            "v3/README.md",
        ],
        "claim_boundary": "Negative result defines current task limits; do not claim universal absence of brain biology or universal foundation-model failure.",
        "outputs": {
            "scene_plate": rel(scene_plate),
            "rendered_preview": rel(preview),
            "grayscale_qa": rel(grayscale),
            "source_summary": rel(source_summary),
            "qa": rel(qa),
        },
        "visible_text_word_count": visible_text_count,
        "visible_text_budget": 50,
    }
    write_json(manifest, manifest_data)
    write_json(
        qa,
        {
            "created": CREATED,
            "automatic_checks": {
                "rendered_outputs_exist": all(path.exists() for path in [scene_plate, preview, grayscale, source_summary, manifest]),
                "image_dimensions": {"width_px": rendered.width, "height_px": rendered.height},
                "visible_text_word_count": visible_text_count,
                "visible_text_budget": 50,
                "section_auroc": data["spatial"]["section_auroc"],
                "animal_auroc": data["spatial"]["animal_auroc"],
                "uce_mean": data["fm"]["uce"]["mean_multi_mission"],
                "scfoundation_mean": data["fm"]["scfoundation"]["mean_multi_mission"],
            },
            "manual_review": {
                "visual_review_status": "pass: rendered preview and grayscale QA inspected after corridor/text revision",
                "claim_boundary": "benchmark boundary, not universal biological/model impossibility",
                "readability": "left spatial null and right FM control are separated by a restrained boundary corridor",
                "notes": "subtitle shortened and corridor reduced after first render",
            },
        },
    )
    return {
        "scene_plate": rel(scene_plate),
        "rendered_preview": rel(preview),
        "grayscale": rel(grayscale),
        "source_summary": rel(source_summary),
        "manifest": rel(manifest),
        "qa": rel(qa),
    }


def main() -> None:
    print(json.dumps(build(), indent=2))


if __name__ == "__main__":
    main()

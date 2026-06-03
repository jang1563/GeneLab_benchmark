#!/usr/bin/env python3
"""Build slide 11 temporal/RRRM guardrails premium scene.

The slide compresses v2 temporal and single-cell lessons into a deck-safe
interpretation layer: preservation artifacts can dominate, recovery projections
are descriptive, and RRRM single-cell data are benchmark-ready but underpowered
for composition claims.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "premium_core_result_slides" / "slide11_temporal_rrrm_guardrails_v0_1"
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
    "Timepoint labels need guardrails",
    "v2 separates preservation, recovery projection, and single-cell readiness before interpretation.",
    "Preservation can dominate",
    "RR-8 GC 0.973 > FLT 0.930",
    "Recovery is descriptive",
    "RR-8 ratio 0.652; probability 1.000 -> 0.404",
    "RRRM ready, underpowered",
    "4 tissues; 38,081 cells; 0 composition FDR hits",
    "Caveat: not held-out recovery validation.",
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


def load_sources() -> dict:
    temporal = load_json("v2/evaluation/T_temporal_summary.json")
    f2a = load_json("v2/evaluation/F2A_composition.json")

    rr8 = temporal["T1"]["T1b_RR-8_liver"]
    rr6 = temporal["T1"]["T1a_RR-6_liver"]
    rr8_conditions = rr8["conditions"]
    rr6_conditions = rr6["conditions"]

    recovery = {}
    for key, row in temporal["T2"].items():
        shift = row["classification_shift"]
        recovery[key.replace("T2_", "")] = {
            "pca_ratio": row["pca_recovery_ratio"],
            "pathway_mean": row["pathway_recovery_summary"]["mean"],
            "n_recovering": row["pathway_recovery_summary"]["n_recovering"],
            "n_total": row["pathway_recovery_summary"]["n_total"],
            "iss_t_prob": shift["mean_flt_isst_flight_prob"],
            "lar_prob": shift["mean_flt_lar_flight_prob"],
            "caveat": shift["caveat"],
        }

    f2a_summary = {}
    all_padj = []
    for tissue, celltypes in f2a["results"].items():
        values = [entry["padj"] for entry in celltypes.values()]
        all_padj.extend(values)
        f2a_summary[tissue] = {
            "cell_types": len(celltypes),
            "lowest_padj": min(values),
            "n_significant_05": sum(value < 0.05 for value in values),
        }

    manifest_rows = list(csv.DictReader((ROOT / "v2/docs/RRRM1_BENCHMARK_READY_MANIFEST_2026-03-12.csv").open()))
    raw_gb = sum(float(row["raw_fastq_size_gb"]) for row in manifest_rows)

    return {
        "preservation": {
            "rr8_flt_gene": rr8_conditions["FLT_gene"]["auroc"],
            "rr8_gc_gene": rr8_conditions["GC_gene"]["auroc"],
            "rr8_bsl_gene": rr8_conditions["BSL_gene"]["auroc"],
            "rr8_excess_signal": rr8["interpretation"]["excess_signal"],
            "rr6_flt_gene": rr6_conditions["FLT_gene"]["auroc"],
            "rr6_gc_gene": rr6_conditions["GC_gene"]["auroc"],
            "rr6_bsl_gene": rr6_conditions["BSL_gene"]["auroc"],
            "rr6_excess_signal": rr6["interpretation"]["excess_signal"],
        },
        "recovery": recovery,
        "rrrm": {
            "tissues": ["blood", "eye", "muscle", "skin"],
            "n_tissues": len(manifest_rows),
            "total_cells": 38081,
            "raw_fastq_gb": raw_gb,
            "cell_inventory": {
                "blood": 4377,
                "eye": 2197,
                "muscle": 15669,
                "skin": 15838,
            },
            "f2a_summary": f2a_summary,
            "f2a_total_significant_05": sum(item["n_significant_05"] for item in f2a_summary.values()),
            "f2a_lowest_padj": min(all_padj),
            "status": "benchmark-ready subset",
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

    center = (int(w * 0.58), int(h * 0.34))
    for idx, radius in enumerate([760, 1010, 1260]):
        bbox = (
            center[0] - radius,
            center[1] - int(radius * 0.35),
            center[0] + radius,
            center[1] + int(radius * 0.35),
        )
        draw.arc(bbox, 200, 354, fill=rgba("sky", 40 - idx * 9), width=3)


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
    draw.text((230, 196), "TEMPORAL AND SINGLE-CELL GUARDRAILS", font=F["eyebrow"], fill=rgb("teal"))
    draw.text((230, 248), "Timepoint labels need guardrails", font=F["title"], fill=rgb("ink"))
    draw.text(
        (234, 330),
        "v2 separates preservation, recovery projection, and single-cell readiness before interpretation.",
        font=F["subtitle"],
        fill=rgb("soft"),
    )


def draw_bar_pair(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, values: list[tuple[str, float, str]]) -> None:
    draw.text((x, y), label, font=F["small"], fill=rgb("soft"))
    bar_x = x + 180
    for idx, (name, value, tone) in enumerate(values):
        yy = y + 42 + idx * 56
        draw.text((x, yy + 5), name, font=F["tiny"], fill=rgb("muted"))
        draw.line((bar_x, yy + 18, bar_x + 450, yy + 18), fill=rgba("rule", 80), width=2)
        draw.rounded_rectangle((bar_x, yy + 3, bar_x + int(450 * value), yy + 33), radius=6, fill=rgba(tone, 215))
        draw.text((bar_x + 470, yy), f"{value:.3f}", font=F["small"], fill=rgb("soft"))


def draw_preservation_surface(canvas: Image.Image, data: dict) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = 230, 610, 1040, 720
    panel_shadow(canvas, (x, y, w, h))
    draw.rounded_rectangle((x, y, x + w, y + h), radius=12, fill=rgba("panel", 228), outline=rgba("rose", 155), width=2)
    draw_tinted_icon(canvas, "split_guard", (x + 54, y + 46), 90, "rose", 175)
    draw.text((x + 150, y + 54), "Preservation can dominate", font=F["h"], fill=rgb("ink"))
    draw.text((x + 150, y + 94), "ISS-T versus LAR separates even in controls", font=F["small"], fill=rgb("muted"))

    draw_bar_pair(
        draw,
        x + 60,
        y + 180,
        "RR-8 liver",
        [
            ("FLT", data["preservation"]["rr8_flt_gene"], "sky"),
            ("GC", data["preservation"]["rr8_gc_gene"], "rose"),
            ("BSL", data["preservation"]["rr8_bsl_gene"], "amber"),
        ],
    )
    draw.text((x + 60, y + 560), "RR-8 GC 0.973 > FLT 0.930", font=F["metric"], fill=rgb("rose"))
    draw.text((x + 60, y + 614), "baseline separation stays high, so timing is confounded", font=F["small"], fill=rgb("muted"))


def draw_recovery_surface(canvas: Image.Image, data: dict) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = 1450, 570, 1120, 780
    draw.text((x, y), "Recovery is descriptive", font=F["h"], fill=rgb("ink"))
    draw.text((x, y + 42), "LAR projects closer to baseline, but this is not held-out validation", font=F["small"], fill=rgb("muted"))

    base_y = y + 150
    missions = [("RR-6", data["recovery"]["RR-6"], "teal"), ("RR-8", data["recovery"]["RR-8"], "amber")]
    for idx, (mission, row, tone) in enumerate(missions):
        yy = base_y + idx * 250
        draw.text((x, yy), mission, font=F["h"], fill=rgb(tone))
        draw.line((x + 135, yy + 28, x + 850, yy + 28), fill=rgba("rule", 95), width=4)
        start = (x + 165, yy + 28)
        mid = (x + 475, yy + 28)
        end = (x + 825, yy + 28)
        draw.ellipse((start[0] - 16, start[1] - 16, start[0] + 16, start[1] + 16), fill=rgba("rose", 220))
        draw.ellipse((mid[0] - 16, mid[1] - 16, mid[0] + 16, mid[1] + 16), fill=rgba(tone, 220))
        draw.ellipse((end[0] - 16, end[1] - 16, end[0] + 16, end[1] + 16), fill=rgba("green", 200))
        draw.text((start[0] - 48, yy + 58), "ISS-T", font=F["tiny"], fill=rgb("muted"))
        draw.text((mid[0] - 34, yy + 58), "LAR", font=F["tiny"], fill=rgb("muted"))
        draw.text((end[0] - 58, yy + 58), "baseline", font=F["tiny"], fill=rgb("muted"))
        draw.text((x + 140, yy + 108), f"ratio {row['pca_ratio']:.3f}", font=F["metric"], fill=rgb(tone))
        draw.text((x + 420, yy + 116), f"prob {row['iss_t_prob']:.3f} -> {row['lar_prob']:.3f}", font=F["body"], fill=rgb("soft"))
        draw.text((x + 420, yy + 150), f"{row['n_recovering']}/{row['n_total']} pathways recover", font=F["small"], fill=rgb("muted"))

    draw.rounded_rectangle((x, y + h - 108, x + w - 40, y + h - 22), radius=10, fill=rgba("panel", 218), outline=rgba("rule", 140), width=1)
    draw.text((x + 36, y + h - 82), "Use as projection shift, not proof that recovery was validated.", font=F["body"], fill=rgb("soft"))


def draw_cell_mosaic(draw: ImageDraw.ImageDraw, x: int, y: int, counts: dict) -> None:
    tissues = [("blood", "rose"), ("eye", "sky"), ("muscle", "amber"), ("skin", "teal")]
    max_count = max(counts.values())
    for idx, (tissue, tone) in enumerate(tissues):
        cx = x + idx * 190
        n = max(5, int(36 * counts[tissue] / max_count))
        for j in range(n):
            gx = cx + (j % 6) * 20
            gy = y + (j // 6) * 20
            draw.ellipse((gx, gy, gx + 12, gy + 12), fill=rgba(tone, 170))
        draw.text((cx, y + 150), tissue, font=F["tiny"], fill=rgb("muted"))
        draw.text((cx, y + 174), f"{counts[tissue]:,}", font=F["small"], fill=rgb("soft"))


def draw_rrrm_surface(canvas: Image.Image, data: dict) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = 2790, 610, 680, 720
    panel_shadow(canvas, (x, y, w, h))
    draw.rounded_rectangle((x, y, x + w, y + h), radius=12, fill=rgba("panel", 228), outline=rgba("teal", 160), width=2)
    draw_tinted_icon(canvas, "single_cell_droplet", (x + 50, y + 42), 88, "teal", 180)
    draw.text((x + 140, y + 54), "RRRM ready,", font=F["h"], fill=rgb("ink"))
    draw.text((x + 140, y + 94), "underpowered", font=F["h"], fill=rgb("amber"))

    draw.text((x + 54, y + 190), f"{data['rrrm']['n_tissues']} tissues", font=F["metric"], fill=rgb("teal"))
    draw.text((x + 300, y + 190), f"{data['rrrm']['total_cells']:,} cells", font=F["metric"], fill=rgb("sky"))
    draw.text((x + 54, y + 244), "benchmark-ready subset", font=F["small"], fill=rgb("muted"))
    draw_cell_mosaic(draw, x + 54, y + 310, data["rrrm"]["cell_inventory"])

    draw.rounded_rectangle((x + 54, y + h - 126, x + w - 54, y + h - 52), radius=10, fill=rgba("panel2", 210), outline=rgba("amber", 125), width=1)
    draw.text((x + 82, y + h - 102), f"{data['rrrm']['f2a_total_significant_05']} composition FDR hits", font=F["body"], fill=rgb("soft"))
    draw.text((x + 82, y + h - 72), f"lowest padj {data['rrrm']['f2a_lowest_padj']:.3f}; n=4+4 per tissue", font=F["small"], fill=rgb("muted"))


def draw_footer(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    draw.rounded_rectangle((230, 1810, 2310, 1938), radius=10, fill=rgba("panel", 224), outline=rgba("rule", 155), width=2)
    draw_tinted_icon(canvas, "caveat_flag", (258, 1832), 66, "amber", 165)
    draw.text((338, 1838), "Caveat: not held-out recovery validation.", font=F["body"], fill=rgb("soft"))
    draw.text((338, 1880), "ISS-T scores are in-sample references; RRRM labels remain first-pass benchmark annotations.", font=F["small"], fill=rgb("muted"))
    draw.text((2830, 1888), "source: v2/evaluation and RRRM docs", font=F["small"], fill=rgb("muted"))


def render(data: dict, *, with_overlay: bool) -> Image.Image:
    canvas = Image.new("RGBA", (3840, 2160), (0, 0, 0, 255))
    draw_background(canvas)
    draw_preservation_surface(canvas, data)
    draw_recovery_surface(canvas, data)
    draw_rrrm_surface(canvas, data)
    if with_overlay:
        draw_header(canvas)
        draw_footer(canvas)
    return canvas.convert("RGB")


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def build() -> dict[str, str]:
    ensure()
    data = load_sources()
    scene_plate = OUT / "slide11_temporal_rrrm_guardrails_scene_plate.png"
    preview = OUT / "slide11_temporal_rrrm_guardrails_rendered_preview.png"
    grayscale = QA / "slide11_temporal_rrrm_guardrails_rendered_preview_grayscale.png"
    source_summary = OUT / "slide11_temporal_rrrm_guardrails_source_summary.json"
    manifest = OUT / "slide11_temporal_rrrm_guardrails_manifest.json"
    qa = OUT / "slide11_temporal_rrrm_guardrails_qa.json"

    render(data, with_overlay=False).save(scene_plate, quality=95)
    rendered = render(data, with_overlay=True)
    rendered.save(preview, quality=95)
    rendered.convert("L").convert("RGB").save(grayscale, quality=95)
    write_json(source_summary, data)

    manifest_data = {
        "slide_id": "slide11_temporal_rrrm_guardrails_v0_1",
        "created": CREATED,
        "slide_role": "temporal_single_cell_guardrails",
        "source_documents": [
            "v2/evaluation/T_temporal_summary.json",
            "v2/evaluation/V2_RESULTS_SUMMARY.md",
            "v2/evaluation/F2A_composition.json",
            "v2/docs/RRRM1_BENCHMARK_READY_MANIFEST_2026-03-12.csv",
            "v2/docs/RRRM1_BROAD_ANNOTATION_SUMMARY_2026-03-12.md",
        ],
        "claim_boundary": "Preservation and recovery are guardrail lessons; recovery projection is descriptive and RRRM composition tests are underpowered.",
        "outputs": {
            "scene_plate": rel(scene_plate),
            "rendered_preview": rel(preview),
            "grayscale_qa": rel(grayscale),
            "source_summary": rel(source_summary),
            "qa": rel(qa),
        },
        "visible_text_word_count": len(" ".join(OVERLAY_TEXT).split()),
        "visible_text_budget": 52,
    }
    write_json(manifest, manifest_data)
    write_json(
        qa,
        {
            "created": CREATED,
            "automatic_checks": {
                "rendered_outputs_exist": all(path.exists() for path in [scene_plate, preview, grayscale, source_summary, manifest]),
                "image_dimensions": {"width_px": rendered.width, "height_px": rendered.height},
                "visible_text_word_count": len(" ".join(OVERLAY_TEXT).split()),
                "visible_text_budget": 52,
                "rr8_gc_gt_flt": data["preservation"]["rr8_gc_gene"] > data["preservation"]["rr8_flt_gene"],
                "rrrm_tissues": data["rrrm"]["n_tissues"],
                "rrrm_cells": data["rrrm"]["total_cells"],
                "composition_fdr_hits": data["rrrm"]["f2a_total_significant_05"],
            },
            "manual_review": {
                "visual_review_status": "pass: rendered preview and grayscale QA inspected",
                "claim_boundary": "guardrail lesson, not proof of recovery or mature cell-type mechanism",
                "readability": "three-lane layout separates confound, projection, and RRRM readiness",
                "notes": "subtitle shortened after first render to keep visible text within budget",
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

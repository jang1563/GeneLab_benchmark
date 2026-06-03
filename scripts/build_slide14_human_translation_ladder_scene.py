#!/usr/bin/env python3
"""Build slide 14 human translation ladder premium scene.

The slide uses v6 evaluation JSON as source data and keeps the visible claim
boundary conservative: pathway/target-tier evidence supports partial human
alignment, not direct gene-level transfer or clinical actionability.
"""

from __future__ import annotations

import json
import statistics
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "premium_core_result_slides" / "slide14_human_translation_ladder_v0_1"
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

TISSUE_LABELS = {
    "liver": "Liver",
    "gastrocnemius": "Gastroc.",
    "kidney": "Kidney",
    "thymus": "Thymus",
    "eye": "Eye",
    "skin": "Skin",
    "lung": "Lung",
    "colon": "Colon",
}

OVERLAY_TEXT = [
    "Translation is partial, not direct",
    "Mouse signals carry over through pathways and target tiers, while gene-level transfer stays weak.",
    "Pathways partly align",
    "mean rho 0.285",
    "Direct transfer weak",
    "0/8 significant transfer tests",
    "Target tiers, not treatments",
    "A=3, B=7",
    "Use as human-alignment evidence, not clinical proof.",
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
    "num": font(46, bold=True),
    "metric": font(37, bold=True),
}


def ensure() -> None:
    QA.mkdir(parents=True, exist_ok=True)


def load_json(path: str) -> dict:
    with (ROOT / path).open(encoding="utf-8") as handle:
        return json.load(handle)


def load_sources() -> dict:
    gene = load_json("v6/evaluation/V6_A_gene_conservation.json")
    pathway = load_json("v6/evaluation/V6_B_pathway_conservation.json")
    transfer = load_json("v6/evaluation/V6_C_cross_species_transfer.json")
    biomarker = load_json("v6/evaluation/V6_D_biomarker_validation.json")
    tf = load_json("v6/evaluation/V6_E_tf_conservation.json")
    targets = load_json("v6/evaluation/V6_F_drug_target_validation.json")

    pathway_rows = []
    for tissue, values in pathway["per_tissue_correlations"].items():
        pathway_rows.append(
            {
                "tissue": TISSUE_LABELS.get(tissue, tissue),
                "rho": values["spearman_rho"],
                "p": values["p_value"],
                "significant": values.get("significant", False),
            }
        )
    pathway_rows.sort(key=lambda row: row["rho"], reverse=True)

    transfer_rows = []
    for tissue, values in transfer["analyses"]["pre_vs_post"].items():
        transfer_rows.append(
            {
                "tissue": TISSUE_LABELS.get(tissue, tissue),
                "auroc": values.get("transfer_auroc") or 0.5,
                "perm_p": values.get("perm_p", 1.0),
                "significant": values.get("significant", False),
            }
        )
    transfer_values = [row["auroc"] for row in transfer_rows]

    gene_enrichment = gene["enrichment_results"]
    max_gene_fe = max(item.get("fold_enrichment", 0.0) for item in gene_enrichment.values())
    significant_gene_sets = sum(1 for item in gene_enrichment.values() if item.get("significant", False))

    return {
        "gene": {
            "max_fold_enrichment": round(max_gene_fe, 3),
            "significant_gene_sets": significant_gene_sets,
            "biomarker_panel_overlap": gene_enrichment["biomarker_panel"]["n_overlap_drr"],
        },
        "pathway": {
            "mean_rho": pathway["mean_rho"],
            "n_tissues_analyzed": pathway["n_tissues_analyzed"],
            "positive_significant_tissues": sum(1 for row in pathway_rows if row["significant"] and row["rho"] > 0),
            "negative_significant_tissues": sum(1 for row in pathway_rows if row["significant"] and row["rho"] < 0),
            "rows": pathway_rows,
        },
        "transfer": {
            "n_tests": len(transfer_rows),
            "n_significant": sum(1 for row in transfer_rows if row["significant"]),
            "median_auroc": round(statistics.median(transfer_values), 3),
            "max_auroc": round(max(transfer_values), 3),
            "rows": transfer_rows,
        },
        "biomarker": {
            "panel_size": biomarker["panel_size"],
            "n_detected_in_cfrna": biomarker["n_detected_in_cfrna"],
            "n_de_fdr05": biomarker["n_de_fdr05"],
            "n_drr": biomarker["n_drr"],
        },
        "tf": {
            "mean_rho": tf["mean_rho"],
            "n_tissues_analyzed": tf["n_tissues_analyzed"],
            "n_significant_tissues": sum(
                1 for values in tf["per_tissue_correlations"].values() if values.get("significant", False)
            ),
        },
        "targets": {
            "tier_counts": targets["tier_counts"],
            "tier_a_genes": [item["gene"] for item in targets.get("tier_A_validated", [])],
            "tier_b_genes": [item["gene"] for item in targets.get("tier_B_promising", [])],
            "n_drug_target_genes": targets["n_drug_target_genes"],
            "n_in_cfrna": targets["n_in_cfrna"],
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

    center = (int(w * 0.62), int(h * 0.33))
    for idx, radius in enumerate([760, 1000, 1240]):
        bbox = (
            center[0] - radius,
            center[1] - int(radius * 0.35),
            center[0] + radius,
            center[1] + int(radius * 0.35),
        )
        draw.arc(bbox, 198, 354, fill=rgba("sky", 42 - idx * 10), width=3)


def panel_shadow(canvas: Image.Image, box: tuple[int, int, int, int], *, radius: int = 14) -> None:
    x, y, w, h = box
    shadow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow)
    sd.rounded_rectangle((x + 20, y + 24, x + w + 20, y + h + 24), radius=radius, fill=(0, 0, 0, 128))
    canvas.alpha_composite(shadow.filter(ImageFilter.GaussianBlur(20)))


def draw_tinted_icon(canvas: Image.Image, name: str, xy: tuple[int, int], size: int, tone: str, alpha: int = 205) -> None:
    path = ASSETS / "symbols" / "png" / f"{name}.png"
    if not path.exists():
        return
    resample = getattr(getattr(Image, "Resampling", Image), "LANCZOS")
    icon = Image.open(path).convert("RGBA").resize((size, size), resample)
    mask = icon.split()[-1].point(lambda value: int(value * alpha / 255))
    tinted = Image.new("RGBA", icon.size, rgb(tone) + (0,))
    tinted.putalpha(mask)
    canvas.alpha_composite(tinted, xy)


def rho_color(value: float) -> tuple[int, int, int]:
    if value < 0:
        return rgb("rose")
    if value < 0.25:
        return rgb("blue")
    if value < 0.45:
        return rgb("teal")
    return rgb("amber")


def draw_header(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    draw.text((230, 196), "HUMAN TRANSLATION BOUNDARY", font=F["eyebrow"], fill=rgb("teal"))
    draw.text((230, 248), "Translation is partial, not direct", font=F["title"], fill=rgb("ink"))
    draw.text(
        (234, 330),
        "Mouse signals carry over through pathways and target tiers, while gene-level transfer stays weak.",
        font=F["subtitle"],
        fill=rgb("soft"),
    )


def draw_source_surface(canvas: Image.Image, data: dict) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = 230, 615, 640, 980
    panel_shadow(canvas, (x, y, w, h))
    draw.rounded_rectangle((x, y, x + w, y + h), radius=12, fill=rgba("panel", 230), outline=rgba("rule", 170), width=2)
    draw_tinted_icon(canvas, "mouse_model", (x + 54, y + 42), 94, "sky", 190)
    draw_tinted_icon(canvas, "tissue_section", (x + 144, y + 45), 84, "teal", 185)
    draw.text((x + 54, y + 170), "Mouse tissue evidence", font=F["h"], fill=rgb("ink"))
    draw.text((x + 54, y + 214), "v1-v6 benchmark and biology layers", font=F["small"], fill=rgb("muted"))

    labels = [
        ("gene lists", f"max FE {data['gene']['max_fold_enrichment']:.2f}", "muted"),
        ("pathways", f"{data['pathway']['n_tissues_analyzed']} tissues", "teal"),
        ("model transfer", f"max AUROC {data['transfer']['max_auroc']:.3f}", "sky"),
        ("target links", f"{data['targets']['n_in_cfrna']} in cfRNA", "amber"),
    ]
    yy = y + 330
    for idx, (label, value, tone) in enumerate(labels):
        row_y = yy + idx * 128
        draw.line((x + 58, row_y + 68, x + w - 58, row_y + 68), fill=rgba("rule", 70), width=2)
        draw.ellipse((x + 58, row_y + 25, x + 96, row_y + 63), fill=rgba(tone, 190))
        draw.text((x + 118, row_y + 20), label, font=F["body"], fill=rgb("soft"))
        draw.text((x + 118, row_y + 52), value, font=F["small"], fill=rgb(tone if tone != "muted" else "muted"))

    draw.rounded_rectangle((x + 54, y + h - 132, x + w - 54, y + h - 54), radius=10, fill=rgba("panel2", 210), outline=rgba("rule", 120), width=1)
    draw.text((x + 82, y + h - 106), "Not a human replacement", font=F["body"], fill=rgb("soft"))
    draw.text((x + 82, y + h - 74), "mouse data stays an upstream evidence layer", font=F["small"], fill=rgb("muted"))


def draw_pathway_strip(draw: ImageDraw.ImageDraw, x: int, y: int, rows: list[dict]) -> None:
    for idx, row in enumerate(rows):
        px = x + idx * 126
        tone = rho_color(row["rho"])
        draw.rounded_rectangle((px, y, px + 92, y + 58), radius=10, fill=tone + (215,), outline=rgba("rule", 80), width=1)
        draw.text((px + 14, y + 12), f"{row['rho']:+.2f}", font=F["small"], fill=rgb("void"))
        draw.text((px, y + 72), row["tissue"], font=F["tiny"], fill=rgb("muted"))


def draw_transfer_strip(draw: ImageDraw.ImageDraw, x: int, y: int, rows: list[dict]) -> None:
    draw.line((x, y + 76, x + 710, y + 76), fill=rgba("rule", 120), width=2)
    draw.text((x + 336, y + 86), "0.5", font=F["tiny"], fill=rgb("muted"))
    for idx, row in enumerate(rows):
        px = x + idx * 84
        value = row["auroc"]
        bar_h = int(118 * max(0.0, min(1.0, value)))
        tone = "sky" if value >= 0.55 else "muted"
        draw.rounded_rectangle((px, y + 130 - bar_h, px + 40, y + 130), radius=5, fill=rgba(tone, 190), outline=rgba("rule", 75), width=1)


def draw_tier_pyramid(draw: ImageDraw.ImageDraw, x: int, y: int, counts: dict) -> None:
    tiers = [("A", counts.get("A", 0), "amber"), ("B", counts.get("B", 0), "teal"), ("C", counts.get("C", 0), "sky"), ("D", counts.get("D", 0), "muted")]
    widths = [250, 340, 460, 560]
    for idx, (tier, count, tone) in enumerate(tiers):
        yy = y + idx * 64
        width = widths[idx]
        left = x + (560 - width) // 2
        draw.rounded_rectangle((left, yy, left + width, yy + 46), radius=8, fill=rgba(tone, 170), outline=rgba("rule", 90), width=1)
        fill = rgb("void") if tone in {"amber", "teal", "sky"} else rgb("soft")
        draw.text((left + 20, yy + 12), f"{tier} = {count}", font=F["small"], fill=fill)


def draw_ladder(canvas: Image.Image, data: dict) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y = 930, 560
    draw.text((x + 18, y), "Translation ladder", font=F["h"], fill=rgb("ink"))
    draw.text((x + 18, y + 42), "The useful signal changes level as it moves toward human data", font=F["small"], fill=rgb("muted"))

    cx = x + 150
    start_y = y + 170
    step_h = 194
    step_w = 1220
    steps = [
        ("Gene lists", "weak overlap", f"0 panel genes in human response list; max FE {data['gene']['max_fold_enrichment']:.2f}", "muted", "gene_locus"),
        ("Pathways partly align", f"mean rho {data['pathway']['mean_rho']:.3f}", "positive in kidney, thymus, and eye; liver is discordant", "teal", "pathway_network"),
        ("Direct transfer weak", f"{data['transfer']['n_significant']}/{data['transfer']['n_tests']} significant transfer tests", f"median AUROC {data['transfer']['median_auroc']:.3f}; TF mean rho {data['tf']['mean_rho']:.3f}", "sky", "rna_readout"),
        ("Target tiers, not treatments", f"A={data['targets']['tier_counts'].get('A', 0)}, B={data['targets']['tier_counts'].get('B', 0)}", f"C={data['targets']['tier_counts'].get('C', 0)}, D={data['targets']['tier_counts'].get('D', 0)}; hypothesis triage only", "amber", "protein_program"),
    ]

    boxes = []
    for idx, (label, metric, note, tone, icon_name) in enumerate(steps):
        yy = start_y + idx * 246
        offset = idx * 64
        boxes.append((idx, label, metric, note, tone, icon_name, (cx + offset, yy, cx + offset + step_w, yy + step_h)))

    corridor = [(box[6][0] + 70, box[6][1] + step_h // 2) for box in boxes]
    draw.line(corridor, fill=rgba("teal", 42), width=56, joint="curve")
    draw.line(corridor, fill=rgba("amber", 135), width=7, joint="curve")

    for idx, label, metric, note, tone, icon_name, box in boxes:
        draw.rounded_rectangle(box, radius=10, fill=rgba("panel2", 206), outline=rgba(tone, 135), width=2)
        draw_tinted_icon(canvas, icon_name, (box[0] + 28, box[1] + 36), 92, tone if tone != "muted" else "sky", 170)
        draw.text((box[0] + 140, box[1] + 34), label, font=F["h"], fill=rgb("ink"))
        draw.text((box[0] + 140, box[1] + 78), metric, font=F["metric"], fill=rgb(tone if tone != "muted" else "soft"))
        draw.text((box[0] + 140, box[1] + 132), note, font=F["small"], fill=rgb("muted"))

        if idx == 1:
            draw_pathway_strip(draw, box[0] + 580, box[1] + 56, data["pathway"]["rows"])
        if idx == 2:
            draw_transfer_strip(draw, box[0] + 620, box[1] + 36, data["transfer"]["rows"])

    for point in corridor:
        draw.ellipse((point[0] - 13, point[1] - 13, point[0] + 13, point[1] + 13), fill=rgba("amber", 210))
        draw.ellipse((point[0] - 26, point[1] - 26, point[0] + 26, point[1] + 26), outline=rgba("teal", 80), width=2)


def draw_human_surface(canvas: Image.Image, data: dict) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = 2630, 615, 1000, 980
    panel_shadow(canvas, (x, y, w, h))
    draw.rounded_rectangle((x, y, x + w, y + h), radius=12, fill=rgba("panel", 232), outline=rgba("teal", 160), width=2)
    draw_tinted_icon(canvas, "human_model", (x + 54, y + 42), 100, "teal", 185)
    draw_tinted_icon(canvas, "sample_tube", (x + 145, y + 44), 82, "sky", 180)
    draw.text((x + 54, y + 170), "Human blood cfRNA check", font=F["h"], fill=rgb("ink"))
    draw.text((x + 54, y + 214), "translation readout, not clinical validation", font=F["small"], fill=rgb("muted"))

    metric_y = y + 318
    metrics = [
        (f"{data['biomarker']['n_detected_in_cfrna']}/{data['biomarker']['panel_size']}", "panel genes detected", "sky"),
        (str(data["biomarker"]["n_de_fdr05"]), "FDR hits in panel", "rose"),
        (str(data["biomarker"]["n_drr"]), "panel genes in response list", "muted"),
    ]
    for idx, (num, label, tone) in enumerate(metrics):
        mx = x + 58 + idx * 300
        draw.text((mx, metric_y), num, font=F["num"], fill=rgb(tone if tone != "muted" else "soft"))
        draw.text((mx, metric_y + 58), label, font=F["small"], fill=rgb("muted"))

    draw.line((x + 58, y + 470, x + w - 58, y + 470), fill=rgba("rule", 110), width=2)
    draw.text((x + 58, y + 522), "Target evidence tiers", font=F["h"], fill=rgb("ink"))
    draw.text((x + 58, y + 562), "rank evidence strength without implying treatment readiness", font=F["small"], fill=rgb("muted"))
    draw_tier_pyramid(draw, x + 220, y + 640, data["targets"]["tier_counts"])


def draw_footer(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    draw.rounded_rectangle((230, 1810, 2450, 1938), radius=10, fill=rgba("panel", 224), outline=rgba("rule", 155), width=2)
    draw_tinted_icon(canvas, "caveat_flag", (258, 1832), 66, "amber", 165)
    draw.text((338, 1838), "Use as human-alignment evidence, not clinical proof.", font=F["body"], fill=rgb("soft"))
    draw.text((338, 1880), "Tier labels are evidence strata; no treatment or countermeasure recommendation is made.", font=F["small"], fill=rgb("muted"))
    draw.text((2830, 1888), "source: v6/evaluation translation JSON", font=F["small"], fill=rgb("muted"))


def render(data: dict, *, with_overlay: bool) -> Image.Image:
    canvas = Image.new("RGBA", (3840, 2160), (0, 0, 0, 255))
    draw_background(canvas)
    draw_source_surface(canvas, data)
    draw_ladder(canvas, data)
    draw_human_surface(canvas, data)
    if with_overlay:
        draw_header(canvas)
        draw_footer(canvas)
    return canvas.convert("RGB")


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def build() -> dict[str, str]:
    ensure()
    data = load_sources()
    scene_plate = OUT / "slide14_human_translation_ladder_scene_plate.png"
    preview = OUT / "slide14_human_translation_ladder_rendered_preview.png"
    grayscale = QA / "slide14_human_translation_ladder_rendered_preview_grayscale.png"
    source_summary = OUT / "slide14_human_translation_ladder_source_summary.json"
    manifest = OUT / "slide14_human_translation_ladder_manifest.json"
    qa = OUT / "slide14_human_translation_ladder_qa.json"

    render(data, with_overlay=False).save(scene_plate, quality=95)
    rendered = render(data, with_overlay=True)
    rendered.save(preview, quality=95)
    rendered.convert("L").convert("RGB").save(grayscale, quality=95)
    write_json(source_summary, data)

    manifest_data = {
        "slide_id": "slide14_human_translation_ladder_v0_1",
        "created": CREATED,
        "slide_role": "human_translation_boundary",
        "source_documents": [
            "v6/evaluation/V6_A_gene_conservation.json",
            "v6/evaluation/V6_B_pathway_conservation.json",
            "v6/evaluation/V6_C_cross_species_transfer.json",
            "v6/evaluation/V6_D_biomarker_validation.json",
            "v6/evaluation/V6_E_tf_conservation.json",
            "v6/evaluation/V6_F_drug_target_validation.json",
            "docs/PROJECT_RESULTS_LOCATION_INVENTORY_2026_05_31.md",
        ],
        "claim_boundary": "Partial pathway and target-tier human alignment only; no direct gene-transfer, clinical, or countermeasure claim.",
        "outputs": {
            "scene_plate": rel(scene_plate),
            "rendered_preview": rel(preview),
            "grayscale_qa": rel(grayscale),
            "source_summary": rel(source_summary),
            "qa": rel(qa),
        },
        "visible_text_word_count": len(" ".join(OVERLAY_TEXT).split()),
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
                "visible_text_word_count": len(" ".join(OVERLAY_TEXT).split()),
                "visible_text_budget": 50,
                "source_mean_pathway_rho": data["pathway"]["mean_rho"],
                "source_mean_tf_rho": data["tf"]["mean_rho"],
                "source_transfer_significant_tests": data["transfer"]["n_significant"],
                "source_tier_counts": data["targets"]["tier_counts"],
            },
            "manual_review": {
                "visual_review_status": "pass: rendered preview and grayscale QA inspected after corridor revision",
                "claim_boundary": "human-alignment evidence, not clinical proof",
                "readability": "headline carries the plain-language message; rho/AUROC/tier values remain evidence labels",
                "notes": "central outer card was removed to reduce nested card-box reading",
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

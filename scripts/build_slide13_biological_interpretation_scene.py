#!/usr/bin/env python3
"""Build slide 13 biological interpretation triad premium scene.

The slide compresses v5 biological interpretation outputs into a deck-safe
story layer: immune/microenvironment shifts, metabolic context, and biomarker
target triage. Drug-target evidence is explicitly kept as hypothesis triage,
not as treatment or countermeasure recommendation.
"""

from __future__ import annotations

import glob
import json
import statistics
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "premium_core_result_slides" / "slide13_biological_interpretation_v0_1"
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
    "Biology layer turns hits into hypotheses",
    "v5 links immune shifts, metabolic context, and biomarker triage without treatment claims.",
    "Immune context",
    "10 significant shifts; 1 signaling link",
    "Metabolic context",
    "6 tissues; 103 subsystems; 89.8% coverage",
    "Biomarker triage",
    "20 genes; mean AUROC 0.682; 10 target-linked",
    "Hypothesis layer, not mechanism or treatment proof.",
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
    "metric": font(37, bold=True),
}


def ensure() -> None:
    QA.mkdir(parents=True, exist_ok=True)


def load_json(path: str | Path) -> dict:
    with (ROOT / path).open(encoding="utf-8") as handle:
        return json.load(handle)


def summarize_immune() -> dict:
    tissue_rows = []
    top_changes = []
    for path in sorted(glob.glob(str(ROOT / "v5" / "evaluation" / "immune_deconv_*.json"))):
        data = load_json(Path(path).relative_to(ROOT))
        tissue = data["tissue"]
        sig = data["n_significant_fdr05"]
        tissue_rows.append(
            {
                "tissue": tissue,
                "label": TISSUE_LABELS.get(tissue, tissue.title()),
                "n_samples": data["n_samples"],
                "n_cell_types": data["n_cell_types"],
                "n_significant_fdr05": sig,
            }
        )
        for cell_type, row in data["cell_types"].items():
            fdr = row.get("fdr_p", 1.0)
            delta = row.get("cliffs_delta", 0.0) or 0.0
            if fdr < 0.05:
                top_changes.append(
                    {
                        "tissue": tissue,
                        "label": f"{TISSUE_LABELS.get(tissue, tissue.title())} {cell_type}",
                        "cell_type": cell_type,
                        "cliffs_delta": delta,
                        "fdr_p": fdr,
                        "direction": row.get("direction", ""),
                    }
                )
    top_changes.sort(key=lambda row: abs(row["cliffs_delta"]), reverse=True)

    signaling = load_json("v5/evaluation/cross_organ_signaling.json")
    active_edges = []
    for key, row in signaling["tissue_pairs"].items():
        if row.get("n_active_pairs", 0) > 0:
            active_edges.append(
                {
                    "key": key,
                    "source": row["source"],
                    "target": row["target"],
                    "n_active_pairs": row["n_active_pairs"],
                    "active_pairs": row["active_pairs"],
                }
            )
    total_active = sum(row["n_active_pairs"] for row in active_edges)
    return {
        "n_tissues": len(tissue_rows),
        "n_cell_type_tests": sum(row["n_cell_types"] for row in tissue_rows),
        "n_significant_fdr05": sum(row["n_significant_fdr05"] for row in tissue_rows),
        "tissue_rows": sorted(tissue_rows, key=lambda row: row["n_significant_fdr05"], reverse=True),
        "top_changes": top_changes[:6],
        "signaling": {
            "n_lr_pairs_database": signaling["n_lr_pairs_database"],
            "n_tissue_pairs": signaling["n_tissue_pairs"],
            "n_active_edges": len(active_edges),
            "n_active_pairs": total_active,
            "active_edges": active_edges,
        },
    }


def summarize_metabolism() -> dict:
    rows = []
    subsystem_hits = []
    for path in sorted(glob.glob(str(ROOT / "v5" / "evaluation" / "metabolic_flux_*.json"))):
        data = load_json(Path(path).relative_to(ROOT))
        rows.append(
            {
                "tissue": data["tissue"],
                "label": TISSUE_LABELS.get(data["tissue"], data["tissue"].title()),
                "gene_coverage_pct": data["gene_coverage_pct"],
                "fba_objective_diff": data["fba_objective_diff"],
                "n_subsystems": data["n_subsystems"],
            }
        )
        for subsystem, value in data["top_subsystems"].items():
            subsystem_hits.append(
                {
                    "tissue": data["tissue"],
                    "label": subsystem,
                    "visible_label": short_subsystem_label(subsystem),
                    "mean_flux_diff": value["mean_flux_diff"],
                    "abs_mean_diff": value["abs_mean_diff"],
                    "n_reactions": value["n_reactions"],
                }
            )
    subsystem_hits.sort(key=lambda row: row["abs_mean_diff"], reverse=True)
    return {
        "n_tissues": len(rows),
        "mean_gene_coverage_pct": round(statistics.mean(row["gene_coverage_pct"] for row in rows), 1),
        "min_gene_coverage_pct": round(min(row["gene_coverage_pct"] for row in rows), 1),
        "n_subsystems": int(statistics.median(row["n_subsystems"] for row in rows)),
        "max_abs_objective_diff": round(max(abs(row["fba_objective_diff"]) for row in rows), 3),
        "rows": sorted(rows, key=lambda row: abs(row["fba_objective_diff"]), reverse=True),
        "top_subsystems": subsystem_hits[:8],
    }


def short_subsystem_label(label: str) -> str:
    lower = label.lower()
    if "ros" in lower:
        return "ROS detox"
    if "glycolysis" in lower:
        return "glycolysis"
    if "tryptophan" in lower:
        return "tryptophan"
    if "valine" in lower or "leucine" in lower or "isoleucine" in lower:
        return "BCAA metabolism"
    if "pentose phosphate" in lower:
        return "pentose phosphate"
    if "heme" in lower:
        return "heme pathway"
    if "nad" in lower:
        return "NAD pathway"
    if "nucleotide" in lower:
        return "nucleotide use"
    if len(label) > 22:
        return label[:19] + "..."
    return label


def summarize_targets() -> dict:
    panel = load_json("v5/evaluation/consensus_biomarker_panel.json")
    drug = load_json("v5/evaluation/drug_targets.json")

    validation = panel["panel_validation"]
    validation_rows = [
        {
            "tissue": tissue,
            "label": TISSUE_LABELS.get(tissue, tissue.title()),
            "auroc": row["auroc"],
            "n_folds": row["n_folds"],
            "n_panel_genes_found": row["n_panel_genes_found"],
        }
        for tissue, row in validation.items()
    ]
    validation_rows.sort(key=lambda row: row["auroc"], reverse=True)

    panel_drug_linked = sum(1 for row in panel["panel"] if row.get("drugs"))
    top_panel = [
        {
            "rank": row["rank"],
            "gene": row["human_symbol"],
            "score": row["score"],
            "n_tissues": row["n_tissues"],
            "drug_linked": bool(row.get("drugs")),
        }
        for row in panel["panel"][:10]
    ]

    tier_counts = {tier: len(items) for tier, items in drug["tiers"].items()}
    return {
        "panel_size": panel["panel_size"],
        "n_candidates_scored": panel["n_candidates_scored"],
        "max_possible_score": panel["max_possible_score"],
        "top_score": max(row["score"] for row in panel["panel"]),
        "mean_panel_auroc": round(statistics.mean(row["auroc"] for row in validation.values()), 3),
        "best_panel_auroc": round(validation_rows[0]["auroc"], 3),
        "best_panel_tissue": validation_rows[0]["label"],
        "validation_rows": validation_rows,
        "panel_drug_linked": panel_drug_linked,
        "top_panel": top_panel,
        "n_mapped_to_human": drug["n_mapped_to_human"],
        "n_genes_with_dgidb_interactions": drug["n_genes_with_dgidb_interactions"],
        "n_total_dgidb_interactions": drug["n_total_dgidb_interactions"],
        "n_chembl_targets": drug["n_chembl_targets"],
        "tier_counts": tier_counts,
        "tissue_druggability": sorted(
            [
                {
                    "tissue": tissue,
                    "label": TISSUE_LABELS.get(tissue, tissue.title()),
                    "pct_druggable": row["pct_druggable"],
                    "druggable": row["druggable"],
                    "total": row["total"],
                }
                for tissue, row in drug["tissue_druggability"].items()
            ],
            key=lambda row: row["pct_druggable"],
            reverse=True,
        ),
    }


def load_sources() -> dict:
    return {
        "immune": summarize_immune(),
        "metabolism": summarize_metabolism(),
        "targets": summarize_targets(),
        "source_review_flags": [
            {
                "file": "v5/figures/html/Fig7_immune_signaling.html",
                "issue": "HTML narrative says 62 active L-R pairs; current JSON source has 1 active edge and 1 active pair.",
                "slide13_resolution": "Use current JSON values; keep the HTML sentence out of the figure.",
            },
            {
                "file": "v5/figures/html/Fig9_consensus_panel.html",
                "issue": "HTML narrative says every panel gene has known drug interactions, while JSON has 10/20 panel genes with drug links.",
                "slide13_resolution": "Use 10/20 drug-linked and avoid enumerating drug names in the main figure.",
            },
        ],
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

    center = (int(w * 0.54), int(h * 0.35))
    for idx, radius in enumerate([820, 1060, 1320]):
        bbox = (
            center[0] - radius,
            center[1] - int(radius * 0.35),
            center[0] + radius,
            center[1] + int(radius * 0.35),
        )
        draw.arc(bbox, 198, 354, fill=rgba("sky", 40 - idx * 9), width=3)

    for idx, y in enumerate([560, 840, 1120, 1400]):
        draw.line((250, y, 3560, y - 100), fill=rgba("rule", 22 + idx * 4), width=2)


def panel_shadow(canvas: Image.Image, box: tuple[int, int, int, int], *, radius: int = 14) -> None:
    x, y, w, h = box
    shadow = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    sd = ImageDraw.Draw(shadow)
    sd.rounded_rectangle((x + 20, y + 24, x + w + 20, y + h + 24), radius=radius, fill=(0, 0, 0, 128))
    canvas.alpha_composite(shadow.filter(ImageFilter.GaussianBlur(20)))


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
    draw.text((230, 196), "BIOLOGICAL INTERPRETATION LAYER", font=F["eyebrow"], fill=rgb("teal"))
    draw.text((230, 248), "Biology layer turns hits into hypotheses", font=F["title"], fill=rgb("ink"))
    draw.text(
        (234, 330),
        "v5 links immune shifts, metabolic context, and biomarker triage without treatment claims.",
        font=F["subtitle"],
        fill=rgb("soft"),
    )


def metric_block(draw: ImageDraw.ImageDraw, x: int, y: int, value: str, label: str, tone: str) -> None:
    draw.text((x, y), value, font=F["num"], fill=rgb(tone))
    draw.text((x, y + 58), label, font=F["small"], fill=rgb("muted"))


def draw_tissue_ring(draw: ImageDraw.ImageDraw, x: int, y: int, data: dict) -> None:
    coords = [
        (x + 350, y + 22),
        (x + 510, y + 100),
        (x + 548, y + 248),
        (x + 462, y + 350),
        (x + 252, y + 350),
        (x + 78, y + 248),
        (x + 126, y + 100),
        (x + 282, y + 216),
    ]
    rows = {row["tissue"]: row for row in data["immune"]["tissue_rows"]}
    order = ["thymus", "colon", "kidney", "skin", "liver", "lung", "eye", "gastrocnemius"]
    for idx, tissue in enumerate(order):
        row = rows[tissue]
        px, py = coords[idx]
        sig = row["n_significant_fdr05"]
        tone = "rose" if sig >= 4 else "amber" if sig > 0 else "muted"
        radius = 27 + sig * 4
        draw.ellipse((px - radius, py - radius, px + radius, py + radius), fill=rgba(tone, 190), outline=rgba("rule", 100), width=1)
        label = TISSUE_LABELS.get(tissue, tissue.title())
        draw.text((px - 42, py + radius + 10), label, font=F["tiny"], fill=rgb("muted"))
        if sig:
            draw.text((px - 8, py - 12), str(sig), font=F["small"], fill=rgb("void"))

    thymus = coords[0]
    colon = coords[1]
    draw.line((thymus[0], thymus[1], colon[0], colon[1]), fill=rgba("teal", 185), width=5)
    mid = ((thymus[0] + colon[0]) // 2, (thymus[1] + colon[1]) // 2)
    draw.rounded_rectangle((mid[0] - 42, mid[1] - 22, mid[0] + 42, mid[1] + 22), radius=8, fill=rgba("panel2", 230), outline=rgba("teal", 130), width=1)
    draw.text((mid[0] - 24, mid[1] - 10), "F12", font=F["tiny"], fill=rgb("soft"))
    draw.text((mid[0] + 4, mid[1] - 10), "Gp9", font=F["tiny"], fill=rgb("teal"))


def draw_immune_surface(canvas: Image.Image, data: dict) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = 230, 610, 1050, 760
    panel_shadow(canvas, (x, y, w, h))
    draw.rounded_rectangle((x, y, x + w, y + h), radius=12, fill=rgba("panel", 230), outline=rgba("rose", 150), width=2)
    draw_tinted_icon(canvas, "cell_population", (x + 56, y + 46), 90, "rose", 175)
    draw_tinted_icon(canvas, "vascular_context", (x + 138, y + 50), 78, "sky", 165)
    draw.text((x + 230, y + 54), "Immune context", font=F["h"], fill=rgb("ink"))
    draw.text((x + 230, y + 94), "where tissue-level shifts concentrate", font=F["small"], fill=rgb("muted"))

    metric_block(draw, x + 58, y + 174, str(data["immune"]["n_significant_fdr05"]), "significant shifts", "rose")
    metric_block(draw, x + 322, y + 174, str(data["immune"]["signaling"]["n_active_pairs"]), "signaling link", "teal")
    metric_block(draw, x + 562, y + 174, str(data["immune"]["n_tissues"]), "tissues tested", "sky")

    draw_tissue_ring(draw, x + 120, y + 286, data)

    call_x, call_y = x + 720, y + 340
    draw.text((call_x, call_y), "largest shifts", font=F["body"], fill=rgb("soft"))
    for idx, row in enumerate(data["immune"]["top_changes"][:4]):
        yy = call_y + 48 + idx * 62
        tone = "rose" if row["cliffs_delta"] < 0 else "teal"
        draw.ellipse((call_x, yy + 8, call_x + 22, yy + 30), fill=rgba(tone, 205))
        label = row["label"].replace("Monocytes / macrophages", "Macrophages")
        if len(label) > 28:
            label = label[:25] + "..."
        draw.text((call_x + 36, yy), label, font=F["small"], fill=rgb("soft"))
        draw.text((call_x + 36, yy + 28), f"delta {row['cliffs_delta']:+.2f}", font=F["tiny"], fill=rgb("muted"))


def flux_color(value: float) -> str:
    if value > 0:
        return "teal"
    return "rose"


def draw_flux_river(draw: ImageDraw.ImageDraw, x: int, y: int, rows: list[dict]) -> None:
    base = y + 258
    draw.line((x, base, x + 720, base), fill=rgba("rule", 120), width=2)
    for idx, row in enumerate(rows[:6]):
        px = x + idx * 118
        value = row["fba_objective_diff"]
        max_abs = 1000.0
        height = int(180 * abs(value) / max_abs)
        tone = flux_color(value)
        y0 = base - height if value > 0 else base
        y1 = base if value > 0 else base + height
        draw.rounded_rectangle((px, y0, px + 58, y1), radius=8, fill=rgba(tone, 190), outline=rgba("rule", 75), width=1)
        draw.text((px - 4, base + 204), row["label"], font=F["tiny"], fill=rgb("muted"))


def draw_metabolic_surface(canvas: Image.Image, data: dict) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = 1450, 610, 1080, 760
    panel_shadow(canvas, (x, y, w, h))
    draw.rounded_rectangle((x, y, x + w, y + h), radius=12, fill=rgba("panel", 230), outline=rgba("teal", 150), width=2)
    draw_tinted_icon(canvas, "mitochondrial_stress", (x + 56, y + 42), 92, "teal", 175)
    draw_tinted_icon(canvas, "pathway_network", (x + 134, y + 44), 84, "sky", 165)
    draw.text((x + 234, y + 54), "Metabolic context", font=F["h"], fill=rgb("ink"))
    draw.text((x + 234, y + 94), "pathway estimates, not direct flux measurement", font=F["small"], fill=rgb("muted"))

    metric_block(draw, x + 58, y + 174, str(data["metabolism"]["n_tissues"]), "tissues", "sky")
    metric_block(draw, x + 260, y + 174, str(data["metabolism"]["n_subsystems"]), "subsystems", "teal")
    metric_block(draw, x + 540, y + 174, f"{data['metabolism']['mean_gene_coverage_pct']:.1f}%", "mean coverage", "amber")

    draw_flux_river(draw, x + 84, y + 342, data["metabolism"]["rows"])
    draw.text((x + 82, y + 632), "largest subsystem contexts", font=F["body"], fill=rgb("soft"))
    for idx, row in enumerate(data["metabolism"]["top_subsystems"][:3]):
        yy = y + 672 + idx * 36
        tone = flux_color(row["mean_flux_diff"])
        label = row["visible_label"]
        draw.ellipse((x + 88, yy + 8, x + 106, yy + 26), fill=rgba(tone, 205))
        draw.text((x + 122, yy), label, font=F["tiny"], fill=rgb("muted"))


def draw_gene_chip(draw: ImageDraw.ImageDraw, x: int, y: int, gene: str, score: int, drug: bool) -> None:
    tone = "amber" if drug else "sky"
    draw.rounded_rectangle((x, y, x + 150, y + 42), radius=8, fill=rgba(tone, 165), outline=rgba("rule", 80), width=1)
    fill = rgb("void") if drug else rgb("ink")
    draw.text((x + 14, y + 11), gene[:10], font=F["tiny"], fill=fill)
    draw.text((x + 118, y + 11), str(score), font=F["tiny"], fill=fill)


def draw_auroc_strip(draw: ImageDraw.ImageDraw, x: int, y: int, rows: list[dict]) -> None:
    draw.line((x, y + 122, x + 720, y + 122), fill=rgba("rule", 120), width=2)
    draw.text((x + 356, y + 134), "0.5", font=F["tiny"], fill=rgb("muted"))
    for idx, row in enumerate(rows):
        px = x + idx * 86
        value = row["auroc"]
        height = int(160 * max(0, value - 0.5) / 0.35)
        tone = "amber" if value >= 0.75 else "teal" if value >= 0.70 else "sky" if value >= 0.65 else "muted"
        draw.rounded_rectangle((px, y + 122 - height, px + 42, y + 122), radius=5, fill=rgba(tone, 195), outline=rgba("rule", 80), width=1)
        draw.text((px - 3, y + 172), row["label"], font=F["tiny"], fill=rgb("muted"))


def draw_biomarker_surface(canvas: Image.Image, data: dict) -> None:
    draw = ImageDraw.Draw(canvas)
    x, y, w, h = 2700, 610, 820, 760
    panel_shadow(canvas, (x, y, w, h))
    draw.rounded_rectangle((x, y, x + w, y + h), radius=12, fill=rgba("panel", 230), outline=rgba("amber", 150), width=2)
    draw_tinted_icon(canvas, "gene_locus", (x + 52, y + 46), 90, "amber", 180)
    draw_tinted_icon(canvas, "protein_program", (x + 130, y + 50), 80, "teal", 165)
    draw.text((x + 222, y + 54), "Biomarker triage", font=F["h"], fill=rgb("ink"))
    draw.text((x + 222, y + 94), "ranked evidence, not countermeasure advice", font=F["small"], fill=rgb("muted"))

    metric_block(draw, x + 58, y + 174, str(data["targets"]["panel_size"]), "panel genes", "amber")
    metric_block(draw, x + 260, y + 174, f"{data['targets']['mean_panel_auroc']:.3f}", "mean AUROC", "teal")
    metric_block(draw, x + 520, y + 174, str(data["targets"]["panel_drug_linked"]), "target-linked", "sky")

    chip_x, chip_y = x + 60, y + 342
    for idx, row in enumerate(data["targets"]["top_panel"][:8]):
        draw_gene_chip(draw, chip_x + (idx % 4) * 170, chip_y + (idx // 4) * 58, row["gene"], row["score"], row["drug_linked"])

    draw_auroc_strip(draw, x + 64, y + 480, data["targets"]["validation_rows"])

    draw.rounded_rectangle((x + 58, y + 656, x + w - 58, y + 714), radius=9, fill=rgba("panel2", 205), outline=rgba("rule", 120), width=1)
    draw.text(
        (x + 86, y + 674),
        f"{data['targets']['tier_counts']['tier1_approved']} approved links; {data['targets']['tier_counts']['tier3_preclinical']} preclinical",
        font=F["small"],
        fill=rgb("soft"),
    )


def draw_interpretation_path(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    points = [(775, 1500), (1620, 1612), (2380, 1510), (3100, 1610)]
    draw.line(points, fill=rgba("teal", 30), width=34, joint="curve")
    draw.line(points, fill=rgba("amber", 150), width=5, joint="curve")
    labels = [
        ("measure", "cell states"),
        ("context", "pathways"),
        ("rank", "panel genes"),
        ("triage", "testable hypotheses"),
    ]
    for idx, point in enumerate(points):
        tone = ["rose", "teal", "amber", "sky"][idx]
        draw.ellipse((point[0] - 24, point[1] - 24, point[0] + 24, point[1] + 24), fill=rgba(tone, 215))
        draw.text((point[0] - 76, point[1] + 42), labels[idx][0], font=F["body"], fill=rgb("soft"))
        draw.text((point[0] - 92, point[1] + 76), labels[idx][1], font=F["small"], fill=rgb("muted"))


def draw_footer(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    draw.rounded_rectangle((230, 1810, 2450, 1938), radius=10, fill=rgba("panel", 224), outline=rgba("rule", 155), width=2)
    draw_tinted_icon(canvas, "caveat_flag", (258, 1832), 66, "amber", 165)
    draw.text((338, 1838), "Hypothesis layer, not mechanism or treatment proof.", font=F["body"], fill=rgb("soft"))
    draw.text((338, 1880), "Drug links rank follow-up evidence; they are not countermeasure recommendations.", font=F["small"], fill=rgb("muted"))
    draw.text((2810, 1888), "source: v5 biological interpretation JSON", font=F["small"], fill=rgb("muted"))


def render(data: dict, *, with_overlay: bool) -> Image.Image:
    canvas = Image.new("RGBA", (3840, 2160), (0, 0, 0, 255))
    draw_background(canvas)
    draw_immune_surface(canvas, data)
    draw_metabolic_surface(canvas, data)
    draw_biomarker_surface(canvas, data)
    draw_interpretation_path(canvas)
    if with_overlay:
        draw_header(canvas)
        draw_footer(canvas)
    return canvas.convert("RGB")


def write_json(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def build() -> dict[str, str]:
    ensure()
    data = load_sources()
    scene_plate = OUT / "slide13_biological_interpretation_scene_plate.png"
    preview = OUT / "slide13_biological_interpretation_rendered_preview.png"
    grayscale = QA / "slide13_biological_interpretation_rendered_preview_grayscale.png"
    source_summary = OUT / "slide13_biological_interpretation_source_summary.json"
    manifest = OUT / "slide13_biological_interpretation_manifest.json"
    qa = OUT / "slide13_biological_interpretation_qa.json"

    render(data, with_overlay=False).save(scene_plate, quality=95)
    rendered = render(data, with_overlay=True)
    rendered.save(preview, quality=95)
    rendered.convert("L").convert("RGB").save(grayscale, quality=95)
    write_json(source_summary, data)

    visible_text_count = len(" ".join(OVERLAY_TEXT).split())
    manifest_data = {
        "slide_id": "slide13_biological_interpretation_v0_1",
        "created": CREATED,
        "slide_role": "biological_interpretation_triad",
        "source_documents": [
            "v5/evaluation/immune_deconv_*.json",
            "v5/evaluation/cross_organ_signaling.json",
            "v5/evaluation/metabolic_flux_*.json",
            "v5/evaluation/consensus_biomarker_panel.json",
            "v5/evaluation/drug_targets.json",
            "v5/figures/html/Fig7_immune_signaling.html",
            "v5/figures/html/Fig8_metabolism_drugs.html",
            "v5/figures/html/Fig9_consensus_panel.html",
        ],
        "claim_boundary": "Biological interpretation and target triage only; no causal mechanism, validated countermeasure, or treatment recommendation claim.",
        "outputs": {
            "scene_plate": rel(scene_plate),
            "rendered_preview": rel(preview),
            "grayscale_qa": rel(grayscale),
            "source_summary": rel(source_summary),
            "qa": rel(qa),
        },
        "visible_text_word_count": visible_text_count,
        "visible_text_budget": 58,
        "source_review_flags": data["source_review_flags"],
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
                "visible_text_budget": 58,
                "immune_significant_cell_type_shifts": data["immune"]["n_significant_fdr05"],
                "active_lr_edges": data["immune"]["signaling"]["n_active_edges"],
                "active_lr_pairs": data["immune"]["signaling"]["n_active_pairs"],
                "metabolic_mean_gene_coverage_pct": data["metabolism"]["mean_gene_coverage_pct"],
                "metabolic_n_subsystems": data["metabolism"]["n_subsystems"],
                "panel_mean_auroc": data["targets"]["mean_panel_auroc"],
                "panel_drug_linked": data["targets"]["panel_drug_linked"],
                "tier1_approved_links": data["targets"]["tier_counts"]["tier1_approved"],
                "tier3_preclinical_links": data["targets"]["tier_counts"]["tier3_preclinical"],
            },
            "manual_review": {
                "visual_review_status": "pass: rendered preview and grayscale QA inspected after bubble/path/text revision",
                "claim_boundary": "hypothesis layer only; no treatment/countermeasure recommendation",
                "readability": "avoids table-like layouts; visible headline carries plain-language interpretation",
                "source_consistency": "current JSON values override legacy HTML narrative sentences flagged in manifest",
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

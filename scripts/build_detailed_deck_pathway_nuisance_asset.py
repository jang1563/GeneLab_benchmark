#!/usr/bin/env python3
"""Build the detailed-deck pathway nuisance-signal explainer asset."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
D_SUMMARY = ROOT / "evaluation" / "D_condition_summary.json"
J5_SOURCE = ROOT / "evaluation" / "J5_gene_vs_pathway.json"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "pathway_nuisance"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "panel2": "#151F2D",
    "grid": "#2A3546",
    "axis": "#98A7BA",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#667487",
    "blue": "#66A6E8",
    "sky": "#73A7FF",
    "teal": "#5FD3C4",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "rose": "#E17882",
    "violet": "#A887FF",
    "white": "#FFFFFF",
}

LABEL_TASKS = [
    ("D3", "Mission ID", "liver", "6-way mission label"),
    ("D5_liver", "Hardware", "liver", "RR vs MHU"),
    ("D5_thymus", "Hardware", "thymus", "RR vs MHU"),
    ("D4", "Mouse strain", "thymus", "strain label"),
    ("D6_liver", "Gravity", "liver", "AG / GC / uG"),
    ("D6_thymus", "Gravity", "thymus", "AG / GC / uG"),
]

TISSUE_LABELS = {
    "A_kidney": "Kidney",
    "A_eye": "Eye",
    "A_thymus": "Thymus",
    "A_liver": "Liver",
    "A_gastrocnemius": "Gastrocnemius",
}


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(value: str, alpha: int) -> tuple[int, int, int, int]:
    return (*hex_to_rgb(value), alpha)


def mix(a: str, b: str, t: float) -> tuple[int, int, int]:
    ar, ag, ab = hex_to_rgb(a)
    br, bg, bb = hex_to_rgb(b)
    return (int(ar + (br - ar) * t), int(ag + (bg - ag) * t), int(ab + (bb - ab) * t))


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Helvetica.ttf",
        "/Library/Fonts/Arial Bold.ttf" if bold else "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            continue
    return ImageFont.load_default()


F = {
    "kicker": load_font(34, True),
    "title": load_font(76, True),
    "subtitle": load_font(36, False),
    "h2": load_font(42, True),
    "h3": load_font(33, True),
    "body": load_font(29, False),
    "body_bold": load_font(29, True),
    "small": load_font(25, False),
    "small_bold": load_font(25, True),
    "tiny": load_font(21, False),
    "tiny_bold": load_font(21, True),
    "micro": load_font(17, False),
    "micro_bold": load_font(17, True),
    "stat": load_font(88, True),
}


def rounded(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    radius: int,
    fill: str | tuple[int, int, int, int],
    outline: str | tuple[int, int, int, int] | None = None,
    width: int = 1,
) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int | float, int | float],
    value: str,
    font: ImageFont.ImageFont,
    fill: str | tuple[int, int, int] | tuple[int, int, int, int] = COLORS["text"],
    anchor: str | None = None,
) -> None:
    draw.text(xy, value, font=font, fill=fill, anchor=anchor)


def wrap_lines(draw: ImageDraw.ImageDraw, label: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
    words = label.split()
    lines: list[str] = []
    current: list[str] = []
    for word in words:
        trial = " ".join(current + [word])
        if draw.textlength(trial, font=font) <= max_width:
            current.append(word)
        else:
            if current:
                lines.append(" ".join(current))
            current = [word]
    if current:
        lines.append(" ".join(current))
    return lines


def paragraph(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    body: str,
    font: ImageFont.ImageFont,
    max_width: int,
    fill: str | tuple[int, int, int, int],
    leading: int = 8,
) -> int:
    x, y = xy
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def dashed_vertical(draw: ImageDraw.ImageDraw, x: int, y0: int, y1: int, color: str, width: int = 3, dash: int = 18) -> None:
    y = y0
    while y < y1:
        draw.line((x, y, x, min(y + dash, y1)), fill=color, width=width)
        y += dash * 2


def load_data() -> dict[str, object]:
    d_summary = json.loads(D_SUMMARY.read_text())
    tasks = d_summary["tasks"] if "tasks" in d_summary else d_summary
    label_rows = []
    for task_id, label, tissue, note in LABEL_TASKS:
        entry = tasks[task_id]
        gene = entry["feature_modes"]["gene"]
        pathway = entry["feature_modes"]["pathway"]
        label_rows.append(
            {
                "task": task_id,
                "label": label,
                "tissue": tissue,
                "note": note,
                "gene": float(gene["macro_f1"]),
                "pathway": float(pathway["macro_f1"]),
                "delta": float(pathway["macro_f1"]) - float(gene["macro_f1"]),
                "gene_features": int(gene["n_features"]),
                "pathway_features": int(pathway["n_features"]),
                "pathway_p": float(pathway["perm_p"]),
            }
        )

    j5 = json.loads(J5_SOURCE.read_text())["j5_table"]
    detection_rows = []
    for rec in j5:
        if rec["category"] != "A" or rec["task"] not in TISSUE_LABELS:
            continue
        detection_rows.append(
            {
                "task": rec["task"],
                "label": TISSUE_LABELS[rec["task"]],
                "gene": float(rec["gene"]),
                "pathway": float(rec["pathway"]),
                "delta": float(rec["diff"]),
                "winner": rec["winner"],
            }
        )
    detection_rows.sort(key=lambda row: float(row["delta"]), reverse=True)
    return {"label_rows": label_rows, "detection_rows": detection_rows}


def stat_badge(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, label: str, value: str, accent: str) -> None:
    rounded(draw, (x, y, x + w, y + 104), 26, "#121B28", "#2A394D", 2)
    text(draw, (x + 28, y + 24), label, F["tiny"], accent)
    text(draw, (x + 28, y + 56), value, F["small_bold"], COLORS["text"])


def draw_header(draw: ImageDraw.ImageDraw) -> None:
    badges = [
        ("GENE VIEW", "23k-28k genes", 330, COLORS["blue"]),
        ("PATHWAY VIEW", "50 Hallmark", 315, COLORS["teal"]),
        ("LABEL CHECKS", "6", 220, COLORS["amber"]),
        ("DETECTION TASKS", "5", 300, COLORS["green"]),
    ]
    bx = 2180
    for label, value, width, accent in badges:
        stat_badge(draw, bx, 72, width, label, value, accent)
        bx += width + 28


def draw_gene_matrix(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, h: int) -> None:
    rounded(draw, (x, y, x + w, y + h), 22, "#121B28", "#2A394D", 2)
    text(draw, (x + 28, y + 24), "Gene-level view", F["small_bold"], COLORS["blue"])
    text(draw, (x + 28, y + 60), "samples x genes", F["tiny"], COLORS["muted"])
    grid_x, grid_y = x + 34, y + 116
    cell = 18
    for r in range(8):
        for c in range(16):
            value = (r * 7 + c * 11) % 23
            color = COLORS["rose"] if value < 6 else (COLORS["amber"] if value < 12 else (COLORS["blue"] if value < 18 else COLORS["teal"]))
            alpha = 130 + (value % 5) * 20
            rounded(draw, (grid_x + c * (cell + 6), grid_y + r * (cell + 6), grid_x + c * (cell + 6) + cell, grid_y + r * (cell + 6) + cell), 5, rgba(color, alpha))
    text(draw, (x + 28, y + h - 70), "High-dimensional signal surface", F["tiny_bold"], COLORS["text"])
    text(draw, (x + 28, y + h - 38), "many genes can encode mission context", F["tiny"], COLORS["muted"])


def draw_pathway_vector(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, h: int) -> None:
    rounded(draw, (x, y, x + w, y + h), 22, "#121B28", "#2A394D", 2)
    text(draw, (x + 28, y + 24), "Pathway-summary view", F["small_bold"], COLORS["teal"])
    text(draw, (x + 28, y + 60), "samples x 50 pathways", F["tiny"], COLORS["muted"])
    bar_x, bar_y = x + 44, y + 130
    for i in range(6):
        bar_w = int((0.28 + ((i * 17) % 61) / 100) * (w - 88))
        color = [COLORS["teal"], COLORS["blue"], COLORS["green"], COLORS["amber"]][i % 4]
        y_i = bar_y + i * 29
        rounded(draw, (bar_x, y_i, bar_x + bar_w, y_i + 18), 8, rgba(color, 210))
        draw.line((bar_x, y_i + 28, bar_x + w - 88, y_i + 28), fill=rgba(COLORS["grid"], 120), width=1)
    text(draw, (x + 28, y + h - 70), "Biology-aligned compression", F["tiny_bold"], COLORS["text"])
    text(draw, (x + 28, y + h - 38), "less room for mission fingerprinting", F["tiny"], COLORS["muted"])


def draw_feature_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "A. Same samples, different feature view", F["h2"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "The benchmark can score the same train/test split from raw genes or from compact pathway activity summaries.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    draw_gene_matrix(draw, x0 + 58, y0 + 250, x1 - x0 - 116, 420)
    arrow_y = y0 + 750
    draw.line((x0 + 180, arrow_y, x1 - 180, arrow_y), fill=COLORS["axis"], width=6)
    draw.polygon([(x1 - 180, arrow_y), (x1 - 224, arrow_y - 22), (x1 - 224, arrow_y + 22)], fill=COLORS["axis"])
    text(draw, ((x0 + x1) / 2, arrow_y - 48), "summarize into pathway scores", F["small_bold"], COLORS["amber"], anchor="mm")
    draw_pathway_vector(draw, x0 + 58, y0 + 790, x1 - x0 - 116, 360)

    readout = (x0 + 58, y1 - 280, x1 - 58, y1 - 56)
    rounded(draw, readout, 28, "#211E17", "#69532B", 2)
    text(draw, (readout[0] + 34, readout[1] + 34), "Takeaway", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (readout[0] + 34, readout[1] + 88),
        "If a diagnostic label is easy to recover from genes but less recoverable from pathways, the pathway view is carrying less of that coupled label.",
        F["small"],
        readout[2] - readout[0] - 68,
        COLORS["muted"],
        8,
    )


def draw_label_score_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], rows: list[dict[str, object]]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "B. Coupled-label signal falls in pathway view", F["h2"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "Label-recovery score uses macro-F1: higher means the model can more easily identify mission, hardware, strain, or gravity labels.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    plot_x0, plot_x1 = x0 + 430, x1 - 170
    plot_y0, plot_y1 = y0 + 300, y1 - 230
    sx = lambda v: int(plot_x0 + max(0.0, min(1.0, v)) * (plot_x1 - plot_x0))

    for value, label in [(0.0, "0.0"), (0.25, "0.25"), (0.50, "0.50"), (0.75, "0.75"), (1.0, "1.0")]:
        x = sx(value)
        draw.line((x, plot_y0 - 18, x, plot_y1), fill=rgba(COLORS["grid"], 145), width=2)
        text(draw, (x, plot_y0 - 44), label, F["micro_bold"], COLORS["axis"], anchor="mm")
    text(draw, ((plot_x0 + plot_x1) / 2, plot_y1 + 42), "label-recovery score (macro-F1)", F["tiny_bold"], COLORS["muted"], anchor="mm")

    row_gap = (plot_y1 - plot_y0) / len(rows)
    for idx, row in enumerate(rows):
        cy = int(plot_y0 + row_gap * idx + row_gap * 0.52)
        label_x = x0 + 58
        accent = COLORS["violet"] if row["task"] == "D4" else COLORS["teal"]
        text(draw, (label_x, cy - 34), str(row["label"]), F["small_bold"], accent)
        text(draw, (label_x, cy - 2), str(row["tissue"]), F["tiny_bold"], COLORS["text"])
        text(draw, (label_x, cy + 29), str(row["note"]), F["micro"], COLORS["muted"])

        gx = sx(float(row["gene"]))
        px = sx(float(row["pathway"]))
        rounded(draw, (plot_x0, cy - 34, plot_x1, cy + 34), 18, "#121B28", rgba(accent, 62), 1)
        draw.line((min(gx, px), cy, max(gx, px), cy), fill=rgba(COLORS["axis"], 185), width=7)
        draw.ellipse((gx - 15, cy - 15, gx + 15, cy + 15), fill=COLORS["blue"], outline=COLORS["bg"], width=3)
        draw.ellipse((px - 15, cy - 15, px + 15, cy + 15), fill=COLORS["teal"], outline=COLORS["bg"], width=3)
        delta = float(row["delta"])
        d_text = f"{delta:+.3f}"
        d_color = COLORS["teal"] if delta > 0 else COLORS["amber"]
        text(draw, (x1 - 58, cy - 10), d_text, F["tiny_bold"], d_color, anchor="ra")

    legend_y = y1 - 165
    draw.ellipse((x0 + 58, legend_y, x0 + 84, legend_y + 26), fill=COLORS["blue"])
    text(draw, (x0 + 98, legend_y + 2), "gene view", F["tiny"], COLORS["muted"])
    draw.ellipse((x0 + 250, legend_y, x0 + 276, legend_y + 26), fill=COLORS["teal"])
    text(draw, (x0 + 290, legend_y + 2), "pathway view", F["tiny"], COLORS["muted"])
    text(draw, (x1 - 58, legend_y + 2), "delta = pathway - gene", F["tiny"], COLORS["muted"], anchor="ra")


def draw_detection_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], rows: list[dict[str, object]]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 44, y0 + 45), "C. Detection changes are task-specific", F["h2"])
    paragraph(
        draw,
        (x0 + 44, y0 + 98),
        "A positive delta means the pathway view improves the held-out detection score for that tissue.",
        F["small"],
        x1 - x0 - 88,
        COLORS["muted"],
        7,
    )

    plot_x0, plot_x1 = x0 + 260, x1 - 88
    plot_y0, plot_y1 = y0 + 320, y0 + 1035
    min_v, max_v = -0.16, 0.34
    sx = lambda v: int(plot_x0 + (v - min_v) / (max_v - min_v) * (plot_x1 - plot_x0))
    zero = sx(0.0)
    dashed_vertical(draw, zero, plot_y0 - 50, plot_y1 + 40, COLORS["axis"], 3, 16)
    for value in [-0.15, 0.0, 0.15, 0.30]:
        x = sx(value)
        draw.line((x, plot_y0 - 24, x, plot_y1 + 24), fill=rgba(COLORS["grid"], 145), width=2)
        text(draw, (x, plot_y0 - 48), f"{value:+.2f}", F["micro_bold"], COLORS["axis"], anchor="mm")

    row_gap = (plot_y1 - plot_y0) / len(rows)
    for idx, row in enumerate(rows):
        cy = int(plot_y0 + row_gap * idx + row_gap * 0.50)
        label = str(row["label"])
        delta = float(row["delta"])
        bar_color = COLORS["teal"] if delta >= 0 else COLORS["dim"]
        text(draw, (x0 + 54, cy - 15), label, F["small_bold"], COLORS["text"])
        text(draw, (x0 + 54, cy + 18), f"gene {float(row['gene']):.3f} -> path {float(row['pathway']):.3f}", F["micro"], COLORS["muted"])
        x_start, x_end = (zero, sx(delta)) if delta >= 0 else (sx(delta), zero)
        rounded(draw, (x_start, cy - 18, x_end, cy + 18), 14, rgba(bar_color, 210))
        draw.ellipse((sx(delta) - 13, cy - 13, sx(delta) + 13, cy + 13), fill=COLORS["teal"] if delta >= 0 else COLORS["rose"], outline=COLORS["bg"], width=3)
        text(draw, (plot_x1, cy - 15), f"{delta:+.3f}", F["tiny_bold"], COLORS["teal"] if delta >= 0 else COLORS["muted"], anchor="ra")

    text(draw, ((plot_x0 + plot_x1) / 2, plot_y1 + 85), "pathway - gene AUROC", F["tiny_bold"], COLORS["muted"], anchor="mm")

    scope = (x0 + 48, y1 - 410, x1 - 48, y1 - 56)
    rounded(draw, scope, 30, "#211E17", "#69532B", 2)
    text(draw, (scope[0] + 34, scope[1] + 38), "Readout frame", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (scope[0] + 34, scope[1] + 94),
        "Use pathway summaries as a selected feature view: they weaken several coupled labels and lift kidney/eye detection, while some tissues remain stronger in the gene view.",
        F["small"],
        scope[2] - scope[0] - 68,
        COLORS["muted"],
        8,
    )


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    footer = (130, 1870, 3710, 2050)
    rounded(draw, footer, 30, "#0F1722", "#27364A", 2)
    text(draw, (footer[0] + 42, footer[1] + 36), "Takeaway", F["h3"], COLORS["teal"])
    paragraph(
        draw,
        (footer[0] + 42, footer[1] + 92),
        "Pathway summaries can reduce selected coupled-label signals while preserving or improving task-specific detection in some tissues.",
        F["small"],
        footer[2] - footer[0] - 84,
        COLORS["muted"],
        6,
    )


def build_manifest(data: dict[str, object], image_path: Path) -> None:
    manifest = {
        "asset": str(image_path.relative_to(ROOT)),
        "sources": [str(D_SUMMARY.relative_to(ROOT)), str(J5_SOURCE.relative_to(ROOT))],
        "label_rows": [
            {
                "task": row["task"],
                "label": row["label"],
                "tissue": row["tissue"],
                "gene_macro_f1": round(float(row["gene"]), 4),
                "pathway_macro_f1": round(float(row["pathway"]), 4),
                "delta": round(float(row["delta"]), 4),
            }
            for row in data["label_rows"]  # type: ignore[index]
        ],
        "detection_rows": [
            {
                "task": row["task"],
                "label": row["label"],
                "gene_auroc": round(float(row["gene"]), 4),
                "pathway_auroc": round(float(row["pathway"]), 4),
                "delta": round(float(row["delta"]), 4),
            }
            for row in data["detection_rows"]  # type: ignore[index]
        ],
        "visible_readout": "Pathway features reduce selected coupled-label signals and provide task-specific detection gains.",
    }
    (OUT_DIR / "pathway_nuisance_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    data = load_data()
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image, "RGBA")

    draw.rectangle((0, 0, W, 250), fill="#0B1018")
    text(draw, (140, 72), "PATHWAY FEATURE CHECK", F["kicker"], COLORS["teal"])
    text(draw, (140, 118), "Pathway Features Reduce Selected Nuisance Signals", F["title"], COLORS["text"])
    text(
        draw,
        (140, 212),
        "The same samples can be scored through genes or compact pathway activity; the feature view changes which labels remain recoverable.",
        F["subtitle"],
        COLORS["muted"],
    )
    draw_header(draw)

    draw_feature_panel(draw, (130, 330, 1035, 1820))
    draw_label_score_panel(draw, (1070, 330, 2545, 1820), data["label_rows"])  # type: ignore[arg-type]
    draw_detection_panel(draw, (2580, 330, 3710, 1820), data["detection_rows"])  # type: ignore[arg-type]
    draw_footer(draw)

    out = OUT_DIR / "pathway_features_reduce_selected_nuisance_premium.png"
    image.save(out, quality=96)
    build_manifest(data, out)
    print(json.dumps({"asset": str(out.relative_to(ROOT)), "manifest": str((OUT_DIR / "pathway_nuisance_manifest.json").relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

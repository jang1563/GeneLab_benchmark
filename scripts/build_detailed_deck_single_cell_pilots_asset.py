#!/usr/bin/env python3
"""Build slide 35 asset: single-cell pilots provide context."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, ImageOps, ImageStat


ROOT = Path(__file__).resolve().parent.parent
ASSET_ROOT = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
)
OUT_DIR = ASSET_ROOT / "single_cell_pilots"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "single_cell_pilots_provide_context_premium.png"
GRAY_PATH = OUT_DIR / "single_cell_pilots_provide_context_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "single_cell_pilots_manifest.json"

RRRM1_QC = ROOT / "v2" / "docs" / "RRRM1_SAMPLE_QC_SUMMARY_2026-03-12.csv"
RRRM1_READY = ROOT / "v2" / "docs" / "RRRM1_BENCHMARK_READY_MANIFEST_2026-03-12.csv"
RRRM1_COMPOSITION = ROOT / "v2" / "evaluation" / "F2A_composition.json"
RRRM1_LOAO = ROOT / "v2" / "evaluation" / "F2C_loao_classifier.json"
RRRM2_LOAO = ROOT / "v3" / "evaluation" / "F5C_rrrm2_loao.json"
RRRM2_BM = ROOT / "v3" / "evaluation" / "F5E_rrrm2_bone_marrow.json"

COLORS = {
    "bg": "#0C111A",
    "bg2": "#091019",
    "header": "#0F1824",
    "footer": "#080D14",
    "panel": "#101823",
    "panel2": "#151F2D",
    "panel3": "#0F1A26",
    "grid": "#263244",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#667487",
    "blue": "#66A6E8",
    "teal": "#5FD3C4",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "orange": "#E69F00",
    "violet": "#B39DFF",
    "rose": "#E17882",
    "red": "#F17C88",
    "ink": "#080D14",
}

TISSUE_COLORS = {
    "blood": "#F17C88",
    "eye": "#66A6E8",
    "muscle": "#F4C26B",
    "skin": "#8BD17C",
}


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(value: str, alpha: int) -> tuple[int, int, int, int]:
    return (*hex_to_rgb(value), alpha)


def blend(a: str, b: str, t: float) -> str:
    ar, ag, ab = hex_to_rgb(a)
    br, bg, bb = hex_to_rgb(b)
    t = max(0.0, min(1.0, t))
    return f"#{int(ar + (br - ar) * t):02x}{int(ag + (bg - ag) * t):02x}{int(ab + (bb - ab) * t):02x}"


def load_font(size: int, bold: bool = False) -> ImageFont.ImageFont:
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
    "title": load_font(84, True),
    "subtitle": load_font(37),
    "h1": load_font(50, True),
    "h2": load_font(40, True),
    "h3": load_font(32, True),
    "body": load_font(29),
    "body_bold": load_font(29, True),
    "small": load_font(25),
    "small_bold": load_font(25, True),
    "tiny": load_font(21),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "stat": load_font(72, True),
    "stat2": load_font(62, True),
    "axis": load_font(20),
}


def rounded(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    radius: int,
    fill: str | tuple[int, int, int, int],
    outline: str | None = None,
    width: int = 1,
) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int | float, int | float],
    value: str,
    font: ImageFont.ImageFont,
    fill: str = COLORS["text"],
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
            continue
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
    fill: str,
    leading: int = 8,
) -> int:
    x, y = xy
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def pill(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, color: str, font: ImageFont.ImageFont = F["tiny_bold"]) -> int:
    width = int(draw.textlength(label, font=font) + 42)
    rounded(draw, (x, y, x + width, y + 42), 18, "#172335", color, 2)
    text(draw, (x + 21, y + 11), label, font, COLORS["text"])
    return x + width + 14


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(230, int(draw.textlength(value, font=F["tiny_bold"]) + 76))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def load_data() -> dict[str, object]:
    with RRRM1_QC.open(newline="", encoding="utf-8") as handle:
        qc = list(csv.DictReader(handle))
    with RRRM1_READY.open(newline="", encoding="utf-8") as handle:
        ready = list(csv.DictReader(handle))

    tissue_cells = {row["tissue"]: int(row["n_cells"]) for row in qc}
    total_cells = sum(tissue_cells.values())
    ready_rows = sum(1 for row in ready if row["benchmark_readiness"] == "ready_for_benchmark")

    composition = json.loads(RRRM1_COMPOSITION.read_text())
    comp_entries = []
    for tissue, cell_types in composition["results"].items():
        for cell_type, result in cell_types.items():
            comp_entries.append(
                {
                    "tissue": tissue,
                    "cell_type": cell_type,
                    "delta": float(result["delta"]),
                    "padj": float(result["padj"]),
                    "p_raw": float(result["p_raw"]),
                }
            )
    comp_entries.sort(key=lambda row: (row["padj"], -abs(row["delta"])))
    fdr_calls = sum(1 for row in comp_entries if row["padj"] < 0.05)

    rrrm1_loao = json.loads(RRRM1_LOAO.read_text())
    rrrm1_top = []
    for tissue, cell_types in rrrm1_loao["results"].items():
        for cell_type, result in cell_types.items():
            rrrm1_top.append(
                {
                    "tissue": tissue,
                    "cell_type": cell_type,
                    "auroc": float(result["auroc"]),
                    "n_cells_total": int(result["n_cells_total"]),
                }
            )
    rrrm1_top.sort(key=lambda row: row["auroc"], reverse=True)

    rrrm2_loao = json.loads(RRRM2_LOAO.read_text())
    pbmc = []
    for cell_type, result in rrrm2_loao["results"]["pbmc"].items():
        pbmc.append(
            {
                "cell_type": cell_type,
                "auroc": float(result["auroc"]),
                "n_cells_total": int(result["n_cells_total"]),
                "ci_low": float(result["ci_low"]),
                "ci_high": float(result["ci_high"]),
            }
        )
    pbmc.sort(key=lambda row: row["auroc"], reverse=True)

    bone = json.loads(RRRM2_BM.read_text())["results"]
    bone_summary = []
    for site, result in bone.items():
        if "cell_type_results" not in result:
            continue
        vals = [float(row["auroc"]) for row in result["cell_type_results"].values()]
        vals_sorted = sorted(vals)
        median = (vals_sorted[len(vals_sorted) // 2 - 1] + vals_sorted[len(vals_sorted) // 2]) / 2
        bone_summary.append(
            {
                "site": site,
                "n_cells_total": int(result["n_cells_total"]),
                "n_cell_types": int(result["n_cell_types"]),
                "median": median,
                "min": min(vals),
                "max": max(vals),
            }
        )

    return {
        "tissue_cells": tissue_cells,
        "total_cells": total_cells,
        "ready_rows": ready_rows,
        "comp_entries": comp_entries,
        "fdr_calls": fdr_calls,
        "rrrm1_top": rrrm1_top,
        "pbmc": pbmc,
        "bone": bone_summary,
    }


def draw_background(draw: ImageDraw.ImageDraw) -> None:
    draw.rectangle((0, 0, W, H), fill=COLORS["bg"])
    for y in range(H):
        draw.line((0, y, W, y), fill=blend(COLORS["bg"], COLORS["bg2"], y / H))
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=COLORS["grid"], width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill="#172234", width=1)
    draw.rectangle((0, 0, W, 315), fill=COLORS["header"])
    draw.rectangle((0, 1840, W, H), fill=COLORS["footer"])
    draw.line((0, 315, W, 315), fill="#1E2B3D", width=2)
    draw.line((0, 1840, W, 1840), fill="#1E2B3D", width=2)


def draw_header(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 76), "BIOLOGY / SINGLE-CELL CONTEXT", F["kicker"], COLORS["teal"])
    x = 2145
    x = badge(draw, x, 66, "RRRM-1", "38,495 cells", COLORS["blue"])
    x = badge(draw, x, 66, "surface", "4 tissues", COLORS["green"])
    x = badge(draw, x, 66, "RRRM-2", "PBMC + marrow", COLORS["amber"])
    badge(draw, x, 66, "readout", "cell-type AUROC", COLORS["violet"])


def draw_title(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 382), "Single-Cell Pilots Provide Context", F["title"], COLORS["text"])
    paragraph(
        draw,
        (155, 493),
        "RRRM pilots add cell-type resolution around bulk signals, then separate signal-rich surfaces from biology that stays closer to baseline.",
        F["subtitle"],
        2250,
        COLORS["muted"],
        10,
    )


def draw_panel_header(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, step: str, title: str, color: str) -> None:
    rounded(draw, (x, y, x + w, y + 1010), 34, COLORS["panel"], "#29374A", 2)
    rounded(draw, (x + 34, y + 34, x + 98, y + 98), 20, "#172335", color, 2)
    text(draw, (x + 66, y + 51), step, F["h3"], COLORS["text"], "ma")
    text(draw, (x + 120, y + 44), title, F["h2"], COLORS["text"])


def draw_inventory_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, data: dict[str, object]) -> None:
    draw_panel_header(draw, x, y, w, "1", "Curated pilot inventory", COLORS["teal"])
    tissue_cells = data["tissue_cells"]
    total = int(data["total_cells"])
    ready_rows = int(data["ready_rows"])
    rrrm1_top = data["rrrm1_top"]

    text(draw, (x + 54, y + 142), f"{total:,}", F["stat"], COLORS["text"])
    text(draw, (x + 56, y + 220), "single cells", F["h3"], COLORS["muted"])
    pill(draw, x + 56, y + 278, f"{ready_rows} ready RRRM-1 tissues", COLORS["green"])
    pill(draw, x + 330, y + 278, "10x 3' scRNA", COLORS["blue"])

    paragraph(
        draw,
        (x + 54, y + 352),
        "The pilot layer gives the talk a biological vocabulary: blood, eye, muscle, and skin each carry distinct cell-type structure.",
        F["body"],
        w - 108,
        COLORS["muted"],
        9,
    )

    bar_x = x + 270
    bar_y = y + 540
    bar_w = w - 395
    max_cells = max(tissue_cells.values())
    order = ["blood", "eye", "muscle", "skin"]
    for i, tissue in enumerate(order):
        yy = bar_y + i * 91
        count = int(tissue_cells[tissue])
        pct = count / total
        color = TISSUE_COLORS[tissue]
        text(draw, (x + 56, yy + 12), tissue.upper(), F["small_bold"], color)
        rounded(draw, (bar_x, yy + 8, bar_x + bar_w, yy + 43), 16, "#0B111A", "#263244", 1)
        fill_w = int(bar_w * count / max_cells)
        rounded(draw, (bar_x, yy + 8, bar_x + fill_w, yy + 43), 16, color)
        text(draw, (bar_x + bar_w + 28, yy + 6), f"{count:,}", F["small_bold"], COLORS["text"])
        text(draw, (bar_x + bar_w + 28, yy + 36), f"{pct:.0%}", F["tiny"], COLORS["muted"])

    examples = [
        next(row for row in rrrm1_top if row["tissue"] == "blood" and row["cell_type"] == "b_cell"),
        next(row for row in rrrm1_top if row["tissue"] == "muscle" and row["cell_type"] == "endothelial"),
        next(row for row in rrrm1_top if row["tissue"] == "muscle" and row["cell_type"] == "t_nk_lymphocyte"),
    ]
    rounded(draw, (x + 54, y + 888, x + w - 54, y + 968), 22, "#172335", "#2A394D", 1)
    text(draw, (x + 82, y + 908), "Cell-type AUROC examples", F["tiny_bold"], COLORS["teal"])
    xx = x + 82
    yy = y + 936
    labels = [
        f"blood B {examples[0]['auroc']:.3f}",
        f"muscle endothelial {examples[1]['auroc']:.3f}",
        f"muscle T/NK {examples[2]['auroc']:.3f}",
    ]
    for label in labels:
        xx = pill(draw, xx, yy, label, COLORS["blue"], F["tiny_bold"])


def draw_composition_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, data: dict[str, object]) -> None:
    draw_panel_header(draw, x, y, w, "2", "Composition screen", COLORS["amber"])
    comp_entries = data["comp_entries"]
    fdr_calls = int(data["fdr_calls"])

    text(draw, (x + 54, y + 145), "34", F["stat2"], COLORS["text"])
    text(draw, (x + 160, y + 155), "cell-type tests", F["h3"], COLORS["muted"])
    text(draw, (x + 590, y + 145), str(fdr_calls), F["stat2"], COLORS["amber"])
    text(draw, (x + 655, y + 155), "FDR < 0.05 calls", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "Composition shifts are visible; the FDR screen marks them as context cues that guide cell-type and pathway reading.",
        F["body"],
        w - 108,
        COLORS["muted"],
        9,
    )

    chart_x0 = x + 410
    chart_x1 = x + w - 165
    zero_x = chart_x0 + int((0 - (-1.05)) / 2.1 * (chart_x1 - chart_x0))
    chart_y0 = y + 438
    row_h = 70
    draw.line((zero_x, chart_y0 - 34, zero_x, chart_y0 + row_h * 7 + 14), fill="#78879A", width=2)
    text(draw, (zero_x, chart_y0 - 60), "delta 0", F["axis"], COLORS["dim"], "mm")
    text(draw, (chart_x1 - 170, chart_y0 - 60), "top examples q=0.088", F["axis"], COLORS["dim"])
    top = comp_entries[:7]
    for i, row in enumerate(top):
        yy = chart_y0 + i * row_h
        tissue = row["tissue"]
        cell_type = row["cell_type"].replace("_", " ")
        delta = float(row["delta"])
        color = COLORS["green"] if delta > 0 else COLORS["rose"]
        tx = chart_x0 + int((delta - (-1.05)) / 2.1 * (chart_x1 - chart_x0))
        text(draw, (x + 56, yy - 7), tissue.upper(), F["micro_bold"], TISSUE_COLORS[tissue])
        label_lines = wrap_lines(draw, cell_type, F["tiny"], 265)
        for j, line in enumerate(label_lines[:2]):
            text(draw, (x + 56, yy + 18 + j * 24), line, F["tiny"], COLORS["muted"])
        draw.line((zero_x, yy + 22, tx, yy + 22), fill=color, width=6)
        draw.ellipse((tx - 13, yy + 9, tx + 13, yy + 35), fill=color)
        text(draw, (chart_x1 + 22, yy + 14), f"{delta:+.2f}", F["small_bold"], COLORS["text"])

    rounded(draw, (x + 54, y + 905, x + w - 54, y + 964), 20, "#172335", "#2A394D", 1)
    text(draw, (x + 82, y + 922), "Reader cue", F["tiny_bold"], COLORS["amber"])
    text(draw, (x + 218, y + 922), "large deltas point to compartments worth following in pathway space", F["tiny"], COLORS["muted"])


def draw_auroc_bar(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    label: str,
    value: float,
    color: str,
    axis_x0: int,
    axis_x1: int,
) -> None:
    text(draw, (x, y + 8), label, F["small_bold"], COLORS["text"])
    draw.line((axis_x0, y + 28, axis_x1, y + 28), fill="#263244", width=15)
    marker = axis_x0 + int((0.5 - 0.3) / 0.7 * (axis_x1 - axis_x0))
    draw.line((marker, y + 4, marker, y + 52), fill="#78879A", width=2)
    fill_x = axis_x0 + int((value - 0.3) / 0.7 * (axis_x1 - axis_x0))
    draw.line((axis_x0, y + 28, fill_x, y + 28), fill=color, width=15)
    draw.ellipse((fill_x - 17, y + 11, fill_x + 17, y + 45), fill=color)
    text(draw, (axis_x1 + 25, y + 7), f"{value:.3f}", F["small_bold"], COLORS["text"])


def draw_extension_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, data: dict[str, object]) -> None:
    draw_panel_header(draw, x, y, w, "3", "Extension check", COLORS["violet"])
    pbmc = data["pbmc"]
    bone = data["bone"]
    nk = pbmc[0]

    text(draw, (x + 54, y + 145), f"{nk['auroc']:.3f}", F["stat2"], COLORS["violet"])
    text(draw, (x + 248, y + 155), "PBMC NK AUROC", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "RRRM-2 shows where cell-type signal separates cleanly and where marrow readouts stay closer to AUROC 0.5.",
        F["body"],
        w - 108,
        COLORS["muted"],
        9,
    )

    axis_x0 = x + 360
    axis_x1 = x + w - 145
    text(draw, (axis_x0, y + 388), "AUROC guide: 0.5 overlap, 1.0 separation", F["tiny"], COLORS["dim"])
    for i, row in enumerate(pbmc[:5]):
        label = row["cell_type"].replace("_", " ")
        color = COLORS["violet"] if i == 0 else COLORS["blue"] if i < 3 else COLORS["teal"]
        draw_auroc_bar(draw, x + 64, y + 432 + i * 80, label, float(row["auroc"]), color, axis_x0, axis_x1)

    y2 = y + 870
    for i, row in enumerate(bone):
        card_x = x + 56 + i * 514
        rounded(draw, (card_x, y2, card_x + 470, y2 + 105), 24, "#172335", "#2A394D", 2)
        label = "Femur BM" if row["site"] == "femur_bm" else "Humerus BM"
        text(draw, (card_x + 25, y2 + 20), label, F["small_bold"], COLORS["text"])
        text(draw, (card_x + 25, y2 + 58), f"median {row['median']:.3f}", F["small"], COLORS["muted"])
        text(draw, (card_x + 275, y2 + 20), f"{row['n_cells_total']:,}", F["small_bold"], COLORS["green"])
        text(draw, (card_x + 275, y2 + 58), f"{row['n_cell_types']} cell types", F["tiny"], COLORS["muted"])


def draw_reader_rule(draw: ImageDraw.ImageDraw) -> None:
    box = (150, 1695, 3690, 1828)
    rounded(draw, box, 30, "#111A27", "#2A394D", 2)
    text(draw, (196, 1736), "Layer decoder", F["h3"], COLORS["teal"])
    steps = [
        ("Inventory", "what cell types exist?"),
        ("Context screen", "which compartments move?"),
        ("Extension check", "where does signal travel?"),
        ("Scoring readiness", "platform section handles status"),
    ]
    x = 690
    for i, (head, body) in enumerate(steps):
        color = [COLORS["teal"], COLORS["amber"], COLORS["violet"], COLORS["green"]][i]
        rounded(draw, (x, 1718, x + 630, 1808), 24, "#172335", color, 2)
        text(draw, (x + 28, 1733), head, F["small_bold"], COLORS["text"])
        text(draw, (x + 28, 1768), body, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            ax = x + 650
            draw.line((ax, 1763, ax + 82, 1763), fill="#6F7E90", width=4)
            draw.polygon([(ax + 82, 1763), (ax + 62, 1752), (ax + 62, 1774)], fill="#6F7E90")
        x += 735


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    source = "Takeaway: RRRM pilots add cell-type context around bulk pathway and tissue-transfer signals."
    scope = "Next: public scoring readiness appears later in the platform extension section."
    paragraph(draw, (150, 1888), source, F["small"], 3440, COLORS["muted"], 7)
    paragraph(draw, (150, 1993), scope, F["small_bold"], 3440, COLORS["text"], 7)


def main() -> None:
    data = load_data()
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)
    draw_background(draw)
    draw_header(draw)
    draw_title(draw)

    panel_y = 625
    panel_w = 1160
    gap = 55
    xs = [150, 150 + panel_w + gap, 150 + (panel_w + gap) * 2]
    draw_inventory_panel(draw, xs[0], panel_y, panel_w, data)
    draw_composition_panel(draw, xs[1], panel_y, panel_w, data)
    draw_extension_panel(draw, xs[2], panel_y, panel_w, data)
    draw_reader_rule(draw)
    draw_footer(draw)

    canvas.save(OUT_PATH, quality=95)
    gray = ImageOps.grayscale(canvas).convert("RGB")
    gray.save(GRAY_PATH, quality=95)

    stat = ImageStat.Stat(gray)
    manifest = {
        "slide": 35,
        "title": "Single-Cell Pilots Provide Context",
        "readout": "RRRM pilots add cell-type context around bulk signals.",
        "asset": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "sources": [
            str(RRRM1_QC.relative_to(ROOT)),
            str(RRRM1_READY.relative_to(ROOT)),
            str(RRRM1_COMPOSITION.relative_to(ROOT)),
            str(RRRM1_LOAO.relative_to(ROOT)),
            str(RRRM2_LOAO.relative_to(ROOT)),
            str(RRRM2_BM.relative_to(ROOT)),
        ],
        "data": {
            "rrrm1_cells": data["total_cells"],
            "rrrm1_tissues": len(data["tissue_cells"]),
            "composition_tests": len(data["comp_entries"]),
            "composition_fdr_calls_lt_0_05": data["fdr_calls"],
            "rrrm2_pbmc_top": data["pbmc"][:5],
            "rrrm2_bone_marrow": data["bone"],
        },
        "visual_qa": {
            "dimensions": [W, H],
            "mean_grayscale": round(stat.mean[0], 2),
            "stddev_grayscale": round(stat.stddev[0], 2),
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"asset": str(OUT_PATH.relative_to(ROOT)), "grayscale": str(GRAY_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

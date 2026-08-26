#!/usr/bin/env python3
"""Build the detailed-deck bulk RNA-seq adaptation-surface asset."""

from __future__ import annotations

import json
import math
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
GENEFORMER_SUMMARY = ROOT / "evaluation" / "geneformer_mouse_gf_all_tissues_summary.json"
SCGPT_SUMMARY = ROOT / "evaluation" / "scgpt_whole_human_all_tissues_summary.json"
V1_CONTENT = ROOT / "docs" / "V1_PAPER_CONTENT.md"
RESULTS_SUMMARY = ROOT / "evaluation" / "RESULTS_SUMMARY.md"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "bulk_surface_hard"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "bulk_rnaseq_hard_surface_for_cell_fms_premium.png"
GRAY_PATH = OUT_DIR / "bulk_rnaseq_hard_surface_for_cell_fms_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "bulk_surface_hard_manifest.json"

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "panel2": "#151F2D",
    "panel3": "#121B28",
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
    "violet": "#B39DFF",
    "magenta": "#D56DFF",
    "white": "#FFFFFF",
}


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(value: str, alpha: int) -> tuple[int, int, int, int]:
    return (*hex_to_rgb(value), alpha)


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
    "subtitle": load_font(35, False),
    "h2": load_font(42, True),
    "h3": load_font(32, True),
    "body": load_font(29, False),
    "body_bold": load_font(29, True),
    "small": load_font(25, False),
    "small_bold": load_font(25, True),
    "tiny": load_font(21, False),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18, False),
    "micro_bold": load_font(18, True),
    "stat": load_font(66, True),
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
    leading: int = 7,
) -> int:
    x, y = xy
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def draw_background(canvas: Image.Image) -> None:
    overlay = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba("#1E2A3A", 28), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba("#1E2A3A", 22), width=1)
    draw.rectangle((0, 0, W, 310), fill=rgba("#121B28", 170))
    draw.rectangle((0, 1855, W, H), fill=rgba("#090E15", 145))
    canvas.alpha_composite(overlay)


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: str, width: int = 4) -> None:
    x0, y0 = start
    x1, y1 = end
    draw.line((x0, y0, x1, y1), fill=rgba(color, 180), width=width)
    direction = 1 if x1 >= x0 else -1
    draw.polygon([(x1, y1), (x1 - direction * 18, y1 - 11), (x1 - direction * 18, y1 + 11)], fill=rgba(color, 210))


def load_data() -> dict[str, object]:
    gf = json.loads(GENEFORMER_SUMMARY.read_text())
    scgpt = json.loads(SCGPT_SUMMARY.read_text())
    deltas = [float(row["delta_vs_baseline"]) for row in scgpt["tasks"].values()]
    return {
        "n_tasks": int(scgpt["n_tasks"]),
        "scgpt_folds": int(scgpt["n_folds_total"]),
        "geneformer_folds": sum(int(row["geneformer_n_folds"]) for row in gf["tissues"].values()),
        "scgpt_mean": float(scgpt["overall_mean_auroc"]),
        "baseline_mean": float(scgpt["baseline_mean_auroc"]),
        "geneformer_mean": float(gf["overall"]["geneformer_mean"]),
        "scgpt_local_gains": sum(1 for delta in deltas if delta > 0),
        "pretrain_cells": "30M / 33M",
    }


def stat_badge(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, label: str, value: str, accent: str) -> None:
    rounded(draw, (x, y, x + w, y + 104), 24, "#121B28", "#2A394D", 2)
    text(draw, (x + 26, y + 23), label, F["tiny_bold"], accent)
    text(draw, (x + 26, y + 57), value, F["small_bold"], COLORS["text"])


def fmt3(value: float) -> str:
    return str(Decimal(str(value)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))


def draw_header(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    text(draw, (136, 76), "MODELS / ADAPTATION SURFACE", F["kicker"], COLORS["sky"])
    text(draw, (136, 122), "Bulk RNA-seq is a hard adaptation surface", F["title"], COLORS["text"])
    text(
        draw,
        (138, 222),
        "Single-cell foundation models begin with cell-level pretraining; SpaceBio-Bench asks for sample-level mission transfer.",
        F["subtitle"],
        COLORS["muted"],
    )
    badges = [
        ("PRETRAIN CELLS", str(data["pretrain_cells"]), 250, COLORS["teal"]),
        ("TASKS", f"{data['n_tasks']} tissues", 210, COLORS["sky"]),
        ("FM FOLDS", f"{data['scgpt_folds']} + {data['geneformer_folds']}", 210, COLORS["violet"]),
        ("LOCAL scGPT", f"{data['scgpt_local_gains']}/6 gains", 218, COLORS["green"]),
        ("SURFACE", "bulk samples", 224, COLORS["amber"]),
    ]
    bx = 2396
    for label, value, width, accent in badges:
        stat_badge(draw, bx, 72, width, label, value, accent)
        bx += width + 18


def draw_cell_cloud(draw: ImageDraw.ImageDraw, cx: int, cy: int, w: int, h: int) -> None:
    colors = [COLORS["teal"], COLORS["violet"], COLORS["sky"], COLORS["green"], COLORS["amber"]]
    for i in range(58):
        angle = i * 2.39996323
        radius = math.sqrt(i / 58) * min(w, h) * 0.46
        x = int(cx + math.cos(angle) * radius * 1.22)
        y = int(cy + math.sin(angle) * radius * 0.80)
        r = 11 + (i % 4) * 3
        color = colors[i % len(colors)]
        draw.ellipse((x - r, y - r, x + r, y + r), fill=rgba(color, 215), outline=rgba("#FFFFFF", 70))
    draw.rounded_rectangle((cx - w // 2, cy - h // 2, cx + w // 2, cy + h // 2), radius=32, outline="#2A394D", width=2)


def draw_token_strip(draw: ImageDraw.ImageDraw, x: int, y: int, labels: list[str], color: str) -> None:
    cx = x
    for label in labels:
        tw = int(draw.textlength(label, font=F["micro_bold"])) + 36
        rounded(draw, (cx, y, cx + tw, y + 42), 17, "#121B28", rgba(color, 150), 1)
        text(draw, (cx + tw / 2, y + 22), label, F["micro_bold"], COLORS["text"], "mm")
        cx += tw + 12


def draw_pretrain_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 46), "A. Pretraining surface", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "The prior is learned from cell-level expression examples before the spaceflight task appears.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        leading=6,
    )
    draw_cell_cloud(draw, x0 + 460, y0 + 430, 720, 360)
    text(draw, (x0 + 94, y0 + 658), "Cell-level atlas", F["h3"], COLORS["teal"])
    paragraph(
        draw,
        (x0 + 94, y0 + 706),
        "Each training example is a cell or cell-like expression profile with gene-context structure.",
        F["tiny"],
        x1 - x0 - 188,
        COLORS["muted"],
        leading=4,
    )
    rounded(draw, (x0 + 72, y0 + 835, x1 - 72, y0 + 1036), 28, COLORS["panel2"], "#2A394D", 2)
    text(draw, (x0 + 112, y0 + 875), "Model inputs", F["h3"], COLORS["text"])
    draw_token_strip(draw, x0 + 112, y0 + 934, ["rank tokens", "expression bins", "gene IDs"], COLORS["violet"])

    rows = [
        ("Geneformer route", "~30M mouse scRNA cells", COLORS["violet"]),
        ("scGPT route", "33M human CellXGene cells", COLORS["teal"]),
    ]
    y = y0 + 1090
    for title, value, accent in rows:
        rounded(draw, (x0 + 72, y, x1 - 72, y + 148), 26, COLORS["panel2"], "#2A394D", 2)
        text(draw, (x0 + 112, y + 34), title, F["h3"], COLORS["text"])
        text(draw, (x0 + 112, y + 82), value, F["small_bold"], accent)
        y += 180


def draw_pinch_card(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    title: str,
    body: str,
    accent: str,
    idx: str,
) -> None:
    rounded(draw, (x, y, x + w, y + 140), 26, COLORS["panel2"], "#2A394D", 2)
    rounded(draw, (x + 24, y + 36, x + 86, y + 98), 17, rgba(accent, 48), rgba(accent, 175), 2)
    text(draw, (x + 55, y + 58), idx, F["tiny_bold"], accent, "mm")
    text(draw, (x + 112, y + 30), title, F["h3"], COLORS["text"])
    paragraph(draw, (x + 112, y + 76), body, F["tiny"], w - 152, COLORS["muted"], leading=4)


def draw_bottleneck_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 46), "B. Adaptation pinch points", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "The task is a sequence of translations from pretraining language into benchmark language.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        leading=6,
    )

    ribbon_y = y0 + 228
    rounded(draw, (x0 + 82, ribbon_y, x1 - 82, ribbon_y + 108), 34, "#121B28", rgba("#73A7FF", 130), 2)
    text(draw, (x0 + 126, ribbon_y + 38), "cell-level prior", F["small_bold"], COLORS["sky"])
    text(draw, (x1 - 126, ribbon_y + 38), "sample-level score", F["small_bold"], COLORS["green"], "ra")
    arrow(draw, (x0 + 360, ribbon_y + 54), (x1 - 360, ribbon_y + 54), COLORS["sky"], width=5)

    cards = [
        ("Aggregation", "Single cells are mixed into one tissue-level bulk profile.", COLORS["teal"]),
        ("Token translation", "Ranks, bins, gene IDs, and ortholog maps define what the encoder receives.", COLORS["violet"]),
        ("Species and corpus", "Mouse bulk missions meet mouse or human single-cell pretraining corpora.", COLORS["amber"]),
        ("Small adaptation set", "The task adapter learns from tens to low hundreds of samples per fold.", COLORS["rose"]),
        ("Mission shift", "The hidden test unit is an entire mission rather than a random sample split.", COLORS["green"]),
    ]
    y = y0 + 410
    for i, (title, body, accent) in enumerate(cards, start=1):
        draw_pinch_card(draw, x0 + 78, y, x1 - x0 - 156, title, body, accent, str(i))
        y += 166

    rounded(draw, (x0 + 78, y1 - 190, x1 - 78, y1 - 58), 28, "#121B28", rgba("#F4C26B", 145), 2)
    text(draw, (x0 + 116, y1 - 150), "Takeaway", F["small_bold"], COLORS["amber"])
    paragraph(
        draw,
        (x0 + 116, y1 - 112),
        "A model-tier result localizes which translation step carried signal.",
        F["body_bold"],
        x1 - x0 - 232,
        COLORS["text"],
        leading=5,
    )


def draw_bulk_vector(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, h: int) -> None:
    rounded(draw, (x, y, x + w, y + h), 28, "#121B28", "#2A394D", 2)
    colors = [COLORS["teal"], COLORS["violet"], COLORS["amber"], COLORS["green"], COLORS["rose"]]
    base_y = y + 46
    for i in range(5):
        row_y = base_y + i * 46
        text(draw, (x + 32, row_y - 4), f"cell mix {i + 1}", F["micro_bold"], rgba(colors[i], 230))
        start = x + 190
        for j in range(20):
            value = 22 + ((i * 7 + j * 5) % 64)
            bar_h = int(value * 0.55)
            bx = start + j * 24
            draw.rounded_rectangle((bx, row_y + 18 - bar_h, bx + 14, row_y + 18), radius=5, fill=rgba(colors[(i + j) % len(colors)], 205))
    draw.line((x + 50, y + h - 82, x + w - 50, y + h - 82), fill=rgba("#98A7BA", 80), width=2)
    text(draw, (x + 52, y + h - 50), "one bulk sample vector", F["tiny_bold"], COLORS["text"])
    text(draw, (x + w - 52, y + h - 50), "genes", F["tiny_bold"], COLORS["axis"], "ra")


def draw_benchmark_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "C. Benchmark surface", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 50, y0 + 98),
        "The score is assigned to bulk RNA-seq samples held out by mission.",
        F["small"],
        x1 - x0 - 100,
        COLORS["muted"],
        leading=6,
    )
    draw_bulk_vector(draw, x0 + 72, y0 + 210, x1 - x0 - 144, 420)

    split_y = y0 + 700
    rounded(draw, (x0 + 72, split_y, x1 - 72, split_y + 276), 30, COLORS["panel2"], "#2A394D", 2)
    text(draw, (x0 + 112, split_y + 38), "Mission-held-out task", F["h3"], COLORS["text"])
    train_box = (x0 + 112, split_y + 106, x0 + 472, split_y + 210)
    test_box = (x1 - 472, split_y + 106, x1 - 112, split_y + 210)
    rounded(draw, train_box, 24, "#121B28", rgba("#73A7FF", 150), 2)
    rounded(draw, test_box, 24, "#121B28", rgba("#8BD17C", 150), 2)
    text(draw, ((train_box[0] + train_box[2]) / 2, train_box[1] + 38), "train missions", F["small_bold"], COLORS["sky"], "mm")
    text(draw, ((test_box[0] + test_box[2]) / 2, test_box[1] + 38), "hidden mission", F["small_bold"], COLORS["green"], "mm")
    arrow(draw, (train_box[2] + 32, split_y + 158), (test_box[0] - 32, split_y + 158), COLORS["axis"], width=4)

    metrics_y = y0 + 1030
    rows = [
        ("shared result", f"scGPT {fmt3(float(data['scgpt_mean']))} vs ref {fmt3(float(data['baseline_mean']))}", COLORS["teal"]),
        ("local gains", f"{data['scgpt_local_gains']}/6 tissues for scGPT", COLORS["green"]),
        ("Geneformer mean", fmt3(float(data["geneformer_mean"])), COLORS["rose"]),
    ]
    metric_gap = 18
    metric_w = (x1 - x0 - 144 - metric_gap * 2) // 3
    for i, (title, value, accent) in enumerate(rows):
        mx = x0 + 72 + i * (metric_w + metric_gap)
        rounded(draw, (mx, metrics_y, mx + metric_w, metrics_y + 150), 24, COLORS["panel2"], "#2A394D", 2)
        text(draw, (mx + 24, metrics_y + 32), title, F["micro_bold"], accent)
        paragraph(draw, (mx + 24, metrics_y + 70), value, F["small_bold"], metric_w - 48, COLORS["text"], leading=4)

    rounded(draw, (x0 + 72, y1 - 198, x1 - 72, y1 - 58), 28, "#121B28", rgba("#66A6E8", 145), 2)
    text(draw, (x0 + 110, y1 - 158), "Use this slide as", F["small_bold"], COLORS["sky"])
    paragraph(
        draw,
        (x0 + 110, y1 - 118),
        "A guide for adaptation work before adding model complexity.",
        F["body_bold"],
        x1 - x0 - 220,
        COLORS["text"],
        leading=5,
    )


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    y = 1872
    rounded(draw, (136, y, 3704, 2042), 32, "#101823", "#2A394D", 2)
    text(draw, (180, y + 38), "Slide 27 readout", F["small_bold"], COLORS["sky"])
    paragraph(
        draw,
        (180, y + 76),
        "The model-tier result is a surface-fit question: cell-level priors must survive aggregation, token translation, species mapping, and mission shift.",
        F["body_bold"],
        2400,
        COLORS["text"],
        leading=6,
    )
    paragraph(
        draw,
        (2730, y + 46),
        "Next: return from model interpretation to method hardening and v4 result-surface expansion.",
        F["tiny_bold"],
        820,
        COLORS["muted"],
        leading=5,
    )
    text(
        draw,
        (140, 2102),
        "Takeaway: bulk RNA-seq creates a hard adaptation surface for models pretrained on single-cell expression.",
        F["micro"],
        COLORS["dim"],
    )
    text(draw, (3704, 2102), "SPACEBIO-BENCH DETAILED DECK / MODELS", F["micro_bold"], COLORS["dim"], "ra")


def main() -> None:
    data = load_data()
    canvas = Image.new("RGBA", (W, H), COLORS["bg"])
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)

    draw_header(draw, data)
    draw_pretrain_panel(draw, (136, 345, 1088, 1815))
    draw_bottleneck_panel(draw, (1128, 345, 2478, 1815))
    draw_benchmark_panel(draw, (2518, 345, 3704, 1815), data)
    draw_footer(draw)

    rgb = canvas.convert("RGB")
    rgb.save(OUT_PATH, quality=95)
    rgb.convert("L").convert("RGB").save(GRAY_PATH, quality=95)
    MANIFEST_PATH.write_text(
        json.dumps(
            {
                "asset": str(OUT_PATH.relative_to(ROOT)),
                "grayscale": str(GRAY_PATH.relative_to(ROOT)),
                "source_files": [
                    str(V1_CONTENT.relative_to(ROOT)),
                    str(RESULTS_SUMMARY.relative_to(ROOT)),
                    str(GENEFORMER_SUMMARY.relative_to(ROOT)),
                    str(SCGPT_SUMMARY.relative_to(ROOT)),
                ],
                "metrics": data,
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps({"asset": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

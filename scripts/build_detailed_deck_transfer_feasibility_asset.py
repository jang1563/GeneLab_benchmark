#!/usr/bin/env python3
"""Build the detailed-deck transfer-feasibility preflight asset.

The slide translates the NES-conservation result into a practical screening
workflow: DGE signatures can be converted to pathway activity, compared across
missions, and used to prioritize full mission-held-out modeling before fitting a
classifier.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "evaluation" / "NES_conservation_vs_transfer.json"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "transfer_feasibility"
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
    "violet": "#B39DFF",
    "white": "#FFFFFF",
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
    "micro": load_font(18, False),
    "micro_bold": load_font(18, True),
    "stat": load_font(82, True),
}

TISSUE_ORDER = ["thymus", "eye", "skin", "liver", "kidney", "gastrocnemius"]

TISSUE_COLORS = {
    "thymus": COLORS["teal"],
    "eye": COLORS["sky"],
    "skin": COLORS["amber"],
    "liver": COLORS["dim"],
    "kidney": COLORS["dim"],
    "gastrocnemius": COLORS["rose"],
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


def rank(vals: list[float]) -> list[float]:
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j + 2) / 2
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def pearson(xs: list[float], ys: list[float]) -> float:
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy)


def spearman(xs: list[float], ys: list[float]) -> float:
    return pearson(rank(xs), rank(ys))


def load_data() -> tuple[list[dict[str, object]], dict[str, float]]:
    source = json.loads(SOURCE.read_text())
    rows = []
    for tissue in TISSUE_ORDER:
        entry = source["data"][tissue]
        label = "Gastrocnemius" if tissue == "gastrocnemius" else tissue.capitalize()
        n_nes = int(entry["n_missions_nes"])
        nes = float(entry["nes_mean_r"])
        auroc = float(entry["transfer_auroc"])
        if tissue in {"thymus", "eye"}:
            lane = "Run first"
            lane_note = "high agreement"
            lane_color = COLORS["teal"] if tissue == "thymus" else COLORS["sky"]
        elif tissue == "skin":
            lane = "Audit mixed signal"
            lane_note = "moderate screen"
            lane_color = COLORS["amber"]
        elif tissue == "gastrocnemius":
            lane = "Add DGE coverage"
            lane_note = "2 NES missions"
            lane_color = COLORS["rose"]
        else:
            lane = "Expect heterogeneity"
            lane_note = "low agreement"
            lane_color = COLORS["dim"]
        rows.append(
            {
                "tissue": tissue,
                "label": label,
                "nes": nes,
                "auroc": auroc,
                "n_nes": n_nes,
                "n_transfer": int(entry["n_missions_transfer"]),
                "lane": lane,
                "lane_note": lane_note,
                "lane_color": lane_color,
            }
        )
    primary = [r for r in rows if r["tissue"] != "gastrocnemius"]
    core = [r for r in rows if r["tissue"] not in {"gastrocnemius", "skin"}]
    stats = {
        "spearman_primary_5": spearman([float(r["nes"]) for r in primary], [float(r["auroc"]) for r in primary]),
        "spearman_core_4": spearman([float(r["nes"]) for r in core], [float(r["auroc"]) for r in core]),
    }
    return rows, stats


def draw_badge(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, label: str, value: str, accent: str) -> None:
    rounded(draw, (x, y, x + w, y + 104), 26, "#121B28", "#2A394D", 2)
    text(draw, (x + 28, y + 24), label, F["tiny_bold"], accent)
    text(draw, (x + 28, y + 56), value, F["small_bold"], COLORS["text"])


def draw_header(draw: ImageDraw.ImageDraw) -> None:
    badges = [
        ("INPUT", "2+ missions", 255, COLORS["sky"]),
        ("SCREEN", "fGSEA NES", 275, COLORS["teal"]),
        ("OUTPUT", "triage lane", 270, COLORS["amber"]),
        ("TIMING", "before model fit", 350, COLORS["green"]),
    ]
    bx = 2390
    for label, value, width, accent in badges:
        draw_badge(draw, bx, 72, width, label, value, accent)
        bx += width + 28


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: str, width: int = 7) -> None:
    x0, y0 = start
    x1, y1 = end
    draw.line((x0, y0, x1, y1), fill=color, width=width)
    if y1 > y0:
        pts = [(x1 - 18, y1 - 24), (x1 + 18, y1 - 24), (x1, y1 + 10)]
    else:
        pts = [(x1 - 18, y1 + 24), (x1 + 18, y1 + 24), (x1, y1 - 10)]
    draw.polygon(pts, fill=color)


def draw_workflow_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "A. Preflight screen", F["h2"], COLORS["text"])
    text(draw, (x0 + 50, y0 + 96), "Pathway agreement is computed before fitting a classifier.", F["small"], COLORS["muted"])

    steps = [
        ("1", "Mission contrasts", "Run task-specific DGE per mission.", COLORS["blue"]),
        ("2", "Pathway activity", "Convert ranked genes to Hallmark NES.", COLORS["teal"]),
        ("3", "Agreement score", "Compare NES vectors across missions.", COLORS["amber"]),
        ("4", "Triage lane", "Choose run-first, audit, or coverage actions.", COLORS["green"]),
    ]
    card_x0, card_x1 = x0 + 56, x1 - 56
    y = y0 + 188
    for idx, title, body, color in steps:
        rounded(draw, (card_x0, y, card_x1, y + 160), 26, "#131D2A", "#2A394D", 2)
        draw.ellipse((card_x0 + 34, y + 38, card_x0 + 98, y + 102), fill=rgba(color, 55), outline=color, width=4)
        text(draw, (card_x0 + 66, y + 56), idx, F["small_bold"], COLORS["text"], anchor="mm")
        text(draw, (card_x0 + 126, y + 34), title, F["h3"], COLORS["text"])
        paragraph(draw, (card_x0 + 126, y + 82), body, F["small"], card_x1 - card_x0 - 170, COLORS["muted"], 6)
        if idx != "4":
            arrow(draw, ((card_x0 + card_x1) // 2, y + 168), ((card_x0 + card_x1) // 2, y + 207), COLORS["dim"], 5)
        y += 220

    readout = (x0 + 56, y1 - 208, x1 - 56, y1 - 52)
    rounded(draw, readout, 24, "#14241F", "#315E51", 2)
    text(draw, (readout[0] + 32, readout[1] + 30), "Takeaway", F["small_bold"], COLORS["teal"])
    body = "Use this screen before model training: decide where full mission-held-out training is worth running first."
    paragraph(draw, (readout[0] + 32, readout[1] + 70), body, F["small"], readout[2] - readout[0] - 70, COLORS["text"], 7)


def bar(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, h: int, value: float, vmin: float, vmax: float, color: str) -> int:
    rounded(draw, (x, y, x + w, y + h), h // 2, "#0D141F", "#26354A", 1)
    frac = max(0.0, min(1.0, (value - vmin) / (vmax - vmin)))
    fill_w = int(w * frac)
    if fill_w > 0:
        rounded(draw, (x, y, x + fill_w, y + h), h // 2, color, None, 0)
    return fill_w


def draw_calibration_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], rows: list[dict[str, object]], stats: dict[str, float]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "B. Benchmark calibration", F["h2"], COLORS["text"])
    text(draw, (x0 + 50, y0 + 96), "Observed transfer calibrates the screen retrospectively.", F["small"], COLORS["muted"])

    score_box = (x0 + 50, y0 + 148, x1 - 50, y0 + 292)
    rounded(draw, score_box, 26, "#151F2D", "#2A394D", 2)
    text(draw, (score_box[0] + 34, score_box[1] + 28), "Primary rank check", F["small_bold"], COLORS["teal"])
    text(draw, (score_box[0] + 34, score_box[1] + 62), f"rho = {stats['spearman_primary_5']:.2f}", F["stat"], COLORS["text"])
    draw.line((score_box[0] + 400, score_box[1] + 30, score_box[0] + 400, score_box[3] - 30), fill="#2A394D", width=2)
    paragraph(
        draw,
        (score_box[0] + 440, score_box[1] + 38),
        "Five-tissue rank set links pathway agreement to observed transfer. Gastrocnemius is shown separately because NES coverage is limited.",
        F["small"],
        score_box[2] - score_box[0] - 485,
        COLORS["muted"],
        7,
    )

    table_x = x0 + 50
    table_y = y0 + 355
    text(draw, (table_x, table_y), "Tissue", F["tiny_bold"], COLORS["axis"])
    text(draw, (table_x + 265, table_y), "Pathway agreement", F["tiny_bold"], COLORS["axis"])
    text(draw, (table_x + 625, table_y), "Observed AUROC", F["tiny_bold"], COLORS["axis"])
    text(draw, (table_x + 875, table_y), "Preflight lane", F["tiny_bold"], COLORS["axis"])

    row_y = table_y + 44
    row_h = 128
    for row in rows:
        color = str(row["lane_color"])
        tissue = str(row["tissue"])
        rounded(draw, (table_x - 10, row_y - 14, x1 - 54, row_y + row_h - 12), 22, "#121B28", "#223048", 1)
        draw.ellipse((table_x + 8, row_y + 24, table_x + 42, row_y + 58), fill=TISSUE_COLORS[tissue])
        text(draw, (table_x + 60, row_y + 10), str(row["label"]), F["small_bold"], COLORS["text"])
        text(draw, (table_x + 60, row_y + 49), f"{int(row['n_nes'])} NES missions", F["tiny"], COLORS["muted"])

        bar(draw, table_x + 265, row_y + 20, 270, 34, float(row["nes"]), 0.0, 0.65, color)
        text(draw, (table_x + 552, row_y + 17), f"{float(row['nes']):.3f}", F["small_bold"], COLORS["text"])
        text(draw, (table_x + 265, row_y + 66), "mean pairwise NES r", F["micro"], COLORS["dim"])

        bar(draw, table_x + 625, row_y + 20, 180, 34, float(row["auroc"]), 0.50, 0.90, TISSUE_COLORS[tissue])
        text(draw, (table_x + 858, row_y + 17), f"{float(row['auroc']):.3f}", F["small_bold"], COLORS["text"], anchor="ra")
        text(draw, (table_x + 625, row_y + 66), "held-out transfer", F["micro"], COLORS["dim"])

        lane_box = (table_x + 875, row_y + 14, x1 - 86, row_y + 59)
        rounded(draw, lane_box, 20, rgba(color, 38), color, 2)
        text(draw, (lane_box[0] + 22, lane_box[1] + 11), str(row["lane"]), F["tiny_bold"], COLORS["text"])
        text(draw, (table_x + 897, row_y + 68), str(row["lane_note"]), F["micro"], COLORS["muted"])
        row_y += row_h + 12

    note_box = (x0 + 50, y1 - 160, x1 - 50, y1 - 50)
    rounded(draw, note_box, 24, "#211E17", "#69532B", 2)
    text(draw, (note_box[0] + 30, note_box[1] + 25), "Use", F["small_bold"], COLORS["amber"])
    paragraph(
        draw,
        (note_box[0] + 105, note_box[1] + 25),
        "Rank tissues or tasks before expensive model sweeps; treat discordant or low-coverage cases as audit targets.",
        F["small"],
        note_box[2] - note_box[0] - 145,
        COLORS["text"],
        6,
    )


def draw_preflight_board(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    rows: list[dict[str, object]],
) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "C. Triage board", F["h2"], COLORS["text"])
    text(draw, (x0 + 50, y0 + 96), "Two pre-training signals: coverage and pathway agreement.", F["small"], COLORS["muted"])

    px0, py0, px1, py1 = x0 + 90, y0 + 210, x1 - 100, y0 + 780
    rounded(draw, (px0, py0, px1, py1), 22, "#0E1520", "#243347", 2)

    def sx(v: float) -> int:
        return int(px0 + (v - 2) / (6 - 2) * (px1 - px0))

    def sy(v: float) -> int:
        return int(py1 - (v - 0.0) / (0.65 - 0.0) * (py1 - py0))

    for n in [2, 3, 4, 5, 6]:
        x = sx(n)
        draw.line((x, py0, x, py1), fill="#243247", width=2)
        text(draw, (x, py1 + 28), str(n), F["tiny"], COLORS["axis"], anchor="ma")
    for val in [0.0, 0.3, 0.6]:
        y = sy(val)
        draw.line((px0, y, px1, y), fill="#243247", width=2)
        text(draw, (px0 - 26, y - 10), f"{val:.1f}", F["tiny"], COLORS["axis"], anchor="rm")

    # Decision zones.
    draw.rectangle((sx(3), py0, px1, sy(0.3)), fill=rgba(COLORS["teal"], 22))
    draw.rectangle((px0, sy(0.3), sx(3), py1), fill=rgba(COLORS["rose"], 18))
    draw.line((sx(3), py0, sx(3), py1), fill=rgba(COLORS["axis"], 160), width=3)
    draw.line((px0, sy(0.3), px1, sy(0.3)), fill=rgba(COLORS["axis"], 160), width=3)
    text(draw, (sx(4.55), sy(0.53)), "Run first", F["tiny_bold"], COLORS["teal"], anchor="mm")
    text(draw, (sx(2.45), sy(0.12)), "Add coverage", F["tiny_bold"], COLORS["rose"], anchor="mm")
    text(draw, (sx(4.8), sy(0.12)), "Audit biology", F["tiny_bold"], COLORS["amber"], anchor="mm")

    offsets = {
        "thymus": (-58, -48),
        "eye": (24, -44),
        "skin": (22, -38),
        "liver": (-78, 22),
        "kidney": (-82, 20),
        "gastrocnemius": (24, 22),
    }
    for row in rows:
        tissue = str(row["tissue"])
        x, y = sx(float(row["n_nes"])), sy(float(row["nes"]))
        color = TISSUE_COLORS[tissue]
        draw.ellipse((x - 22, y - 22, x + 22, y + 22), fill=color, outline=COLORS["white"], width=3)
        label = "Gastro" if tissue == "gastrocnemius" else str(row["label"])
        dx, dy = offsets[tissue]
        text(draw, (x + dx, y + dy), label, F["micro_bold"], COLORS["text"])

    text(draw, ((px0 + px1) // 2, py1 + 70), "NES coverage: number of missions with fGSEA", F["tiny"], COLORS["muted"], anchor="ma")
    y_label = Image.new("RGBA", (510, 54), (0, 0, 0, 0))
    yd = ImageDraw.Draw(y_label)
    yd.text((255, 8), "Pathway activity agreement", font=F["tiny"], fill=COLORS["muted"], anchor="ma")
    y_label = y_label.rotate(90, expand=True)
    canvas_x = x0 + 22
    canvas_y = y0 + 335
    canvas.paste(y_label, (canvas_x, canvas_y), y_label)

    cards = [
        ("High agreement", "Prioritize full mission-held-out training and inspect the learned signal.", COLORS["teal"]),
        ("Low agreement", "Pair the model run with a tissue-specific heterogeneity audit.", COLORS["amber"]),
        ("Low coverage", "Add DGE/fGSEA inputs before assigning feasibility priority.", COLORS["rose"]),
    ]
    cy = y0 + 930
    for title, body, color in cards:
        rounded(draw, (x0 + 60, cy, x1 - 60, cy + 150), 24, "#131D2A", "#2A394D", 2)
        draw.rectangle((x0 + 60, cy, x0 + 74, cy + 150), fill=color)
        text(draw, (x0 + 98, cy + 26), title, F["small_bold"], color)
        paragraph(draw, (x0 + 98, cy + 66), body, F["small"], x1 - x0 - 205, COLORS["muted"], 6)
        cy += 178


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    label_x = 205
    body_x = 390
    text(draw, (label_x, 1925), "Takeaway", F["small_bold"], COLORS["sky"])
    source = (
        "Pathway agreement can prioritize which tissues deserve full mission-held-out training first."
    )
    paragraph(draw, (body_x, 1925), source, F["small"], 3180, COLORS["muted"], 6)
    text(draw, (label_x, 1995), "Readout", F["small_bold"], COLORS["teal"])
    readout = (
        "The screen prioritizes where to spend full modeling effort; it flags feasibility and coverage gaps before classifier training."
    )
    paragraph(draw, (body_x, 1995), readout, F["small"], 3180, COLORS["muted"], 6)


def main() -> None:
    rows, stats = load_data()

    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)

    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 48), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 38), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "CORE BENCHMARK IMPLICATION", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "Screen transfer feasibility before training", F["title"])
    text(
        draw,
        (150, 214),
        "Use fGSEA pathway activity agreement as a preflight signal before expensive mission-held-out model sweeps.",
        F["subtitle"],
        COLORS["muted"],
    )
    draw_header(draw)

    draw_workflow_panel(draw, (150, 350, 1260, 1800))
    draw_calibration_panel(draw, (1320, 350, 2600, 1800), rows, stats)
    draw_preflight_board(canvas, draw, (2660, 350, 3690, 1800), rows)
    draw_footer(draw)

    png = OUT_DIR / "screen_transfer_feasibility_before_training_premium.png"
    canvas.save(png, quality=95)

    gray = OUT_DIR / "screen_transfer_feasibility_before_training_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    manifest = {
        "title": "Screen transfer feasibility before training",
        "source": str(SOURCE.relative_to(ROOT)),
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "screening_workflow": [
            "task-specific DGE per mission",
            "fGSEA Hallmark NES per mission",
            "pairwise pathway activity agreement",
            "preflight triage lane before classifier training",
        ],
        "data": [
            {
                "tissue": row["tissue"],
                "label": row["label"],
                "nes_mean_r": row["nes"],
                "transfer_auroc": row["auroc"],
                "n_missions_nes": row["n_nes"],
                "n_missions_transfer": row["n_transfer"],
                "preflight_lane": row["lane"],
                "lane_note": row["lane_note"],
            }
            for row in rows
        ],
        "statistics": {
            "spearman_primary_5_rank_set": round(stats["spearman_primary_5"], 3),
            "spearman_core_4_rank_set": round(stats["spearman_core_4"], 3),
        },
        "scope": (
            "Transfer AUROC is used as retrospective benchmark calibration. The deployable pre-training "
            "screen uses only fGSEA NES coverage and pathway activity agreement."
        ),
    }
    manifest_path = OUT_DIR / "transfer_feasibility_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

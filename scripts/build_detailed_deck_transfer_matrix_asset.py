#!/usr/bin/env python3
"""Build the detailed-deck transfer-matrix explainer asset."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
TRANSFER_DIR = ROOT / "processed" / "B_cross_mission"
SUMMARY = ROOT / "evaluation" / "B_cross_mission_summary.json"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "transfer_matrix"
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

TISSUES = [
    {"id": "thymus", "label": "Thymus", "accent": COLORS["teal"]},
    {"id": "gastrocnemius", "label": "Gastrocnemius", "accent": COLORS["green"]},
    {"id": "skin", "label": "Skin", "accent": COLORS["blue"]},
    {"id": "eye", "label": "Eye", "accent": COLORS["sky"]},
    {"id": "liver", "label": "Liver", "accent": COLORS["amber"]},
    {"id": "kidney", "label": "Kidney", "accent": COLORS["rose"]},
]

MISSION_ALIAS = {"TBD": "OSD-397"}


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
    "title": load_font(78, True),
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


def dashed_horizontal(draw: ImageDraw.ImageDraw, x0: int, x1: int, y: int, color: str, width: int = 2, dash: int = 18) -> None:
    x = x0
    while x < x1:
        draw.line((x, y, min(x + dash, x1), y), fill=color, width=width)
        x += dash * 2


def read_matrix(tissue: str) -> dict[str, object]:
    path = TRANSFER_DIR / tissue / "B_transfer_matrix_pca_lr.csv"
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle))
    missions = [MISSION_ALIAS.get(m, m) for m in rows[0][1:]]
    pairs: list[dict[str, object]] = []
    matrix: dict[tuple[str, str], float | None] = {}
    for raw_row in rows[1:]:
        train = MISSION_ALIAS.get(raw_row[0], raw_row[0])
        for test, raw in zip(missions, raw_row[1:]):
            value = float(raw) if raw else None
            matrix[(train, test)] = value
            if value is not None:
                pairs.append({"train": train, "test": test, "auroc": value})
    return {"path": path, "missions": missions, "pairs": pairs, "matrix": matrix}


def load_data() -> list[dict[str, object]]:
    summary = json.loads(SUMMARY.read_text())
    rows: list[dict[str, object]] = []
    for item in TISSUES:
        tissue = item["id"]
        matrix_data = read_matrix(tissue)
        pairs = matrix_data["pairs"]
        values = [float(rec["auroc"]) for rec in pairs]  # type: ignore[index]
        stats = summary[tissue]["methods"]["pca_lr"]
        rows.append(
            {
                **item,
                "missions": matrix_data["missions"],
                "pairs": pairs,
                "matrix": matrix_data["matrix"],
                "matrix_path": matrix_data["path"],
                "mean": float(stats["mean_auroc"]),
                "ci_low": float(stats["ci_low"]),
                "ci_high": float(stats["ci_high"]),
                "n_pairs": int(stats["n_valid_pairs"]),
                "n_missions": int(summary[tissue]["n_missions"]),
                "min": min(values),
                "max": max(values),
            }
        )
    return rows


def score_color(value: float) -> str:
    if value >= 0.85:
        return COLORS["teal"]
    if value >= 0.70:
        return COLORS["blue"]
    if value >= 0.50:
        return COLORS["amber"]
    return COLORS["rose"]


def heat_color(value: float | None) -> str:
    if value is None:
        return "#222B38"
    if value < 0.50:
        return "#6B3540"
    if value < 0.70:
        t = (value - 0.50) / 0.20
        return "#%02x%02x%02x" % mix("#6B5130", COLORS["amber"], max(0.0, min(1.0, t)))
    t = (value - 0.70) / 0.30
    return "#%02x%02x%02x" % mix(COLORS["blue"], COLORS["teal"], max(0.0, min(1.0, t)))


def stat_badge(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, label: str, value: str, accent: str) -> None:
    rounded(draw, (x, y, x + w, y + 104), 26, "#121B28", "#2A394D", 2)
    text(draw, (x + 28, y + 24), label, F["tiny"], accent)
    text(draw, (x + 28, y + 56), value, F["small_bold"], COLORS["text"])


def draw_header(draw: ImageDraw.ImageDraw, rows: list[dict[str, object]]) -> None:
    total_pairs = sum(int(row["n_pairs"]) for row in rows)
    badges = [
        ("TISSUES", str(len(rows)), 215, COLORS["teal"]),
        ("DIRECTED PAIRS", str(total_pairs), 340, COLORS["blue"]),
        ("SCORE UNIT", "AUROC", 260, COLORS["amber"]),
        ("REFERENCE", "0.70", 245, COLORS["green"]),
    ]
    bx = 2370
    for label, value, width, accent in badges:
        stat_badge(draw, bx, 72, width, label, value, accent)
        bx += width + 28


def draw_heatmap(
    draw: ImageDraw.ImageDraw,
    row: dict[str, object],
    x: int,
    y: int,
    cell: int,
    label: str,
) -> None:
    missions = list(row["missions"])  # type: ignore[arg-type]
    matrix = row["matrix"]  # type: ignore[assignment]
    accent = str(row["accent"])
    text(draw, (x, y - 46), label, F["small_bold"], accent)
    for i, mission in enumerate(missions):
        short = mission.replace("OSD-397", "OSD")
        text(draw, (x + i * cell + cell / 2, y - 12), short, F["micro"], COLORS["muted"], anchor="mm")
        text(draw, (x - 12, y + i * cell + cell / 2), short, F["micro"], COLORS["muted"], anchor="rm")
    for i, train in enumerate(missions):
        for j, test_mission in enumerate(missions):
            x0, y0 = x + j * cell, y + i * cell
            value = matrix[(train, test_mission)]  # type: ignore[index]
            if i == j:
                fill = "#243041"
            else:
                fill = heat_color(value)
            rounded(draw, (x0 + 2, y0 + 2, x0 + cell - 2, y0 + cell - 2), 9, fill, None)
            if value is not None and i != j:
                fg = COLORS["bg"] if value >= 0.70 else COLORS["text"]
                font = F["micro_bold"] if cell >= 54 else F["micro"]
                text(draw, (x0 + cell / 2, y0 + cell / 2), f"{value:.2f}", font, fg, anchor="mm")


def draw_matrix_grammar(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], rows: list[dict[str, object]]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "A. Directed pairs create the mean", F["h2"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "Rows train the model; columns test it on a different mission. The tissue score is the average of those directed cells.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    thymus = next(row for row in rows if row["id"] == "thymus")
    liver = next(row for row in rows if row["id"] == "liver")
    draw_heatmap(draw, thymus, x0 + 210, y0 + 330, 92, "Thymus matrix")
    draw_heatmap(draw, liver, x0 + 178, y0 + 870, 62, "Liver matrix")

    legend = (x0 + 62, y1 - 210, x1 - 62, y1 - 56)
    rounded(draw, legend, 26, "#151F2D", "#2A394D", 2)
    text(draw, (legend[0] + 34, legend[1] + 32), "Matrix reading move", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (legend[0] + 34, legend[1] + 82),
        "A high mean is strongest when many cells cluster on the right side of the AUROC scale.",
        F["tiny"],
        legend[2] - legend[0] - 68,
        COLORS["muted"],
        7,
    )


def draw_axis_labels(draw: ImageDraw.ImageDraw, sx, x0: int, x1: int, y_top: int, y_bottom: int) -> None:
    for value, label in [(0.0, "0.0"), (0.25, "0.25"), (0.50, "0.50"), (0.70, "0.70"), (0.85, "0.85"), (1.0, "1.0")]:
        x = sx(value)
        color = COLORS["blue"] if value == 0.70 else (COLORS["amber"] if value == 0.50 else COLORS["axis"])
        if value in {0.50, 0.70}:
            dashed_vertical(draw, x, y_top, y_bottom, color, 3, 16)
        else:
            draw.line((x, y_top, x, y_bottom), fill=rgba(COLORS["grid"], 130), width=2)
        text(draw, (x, y_top - 22), label, F["micro_bold"], color, anchor="mm")
    text(draw, (sx(0.50) + 10, y_bottom + 26), "chance", F["micro_bold"], COLORS["amber"])
    text(draw, (sx(0.70) + 10, y_bottom + 26), "transfer reference", F["micro_bold"], COLORS["blue"])


def draw_pair_spread(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], rows: list[dict[str, object]]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "B. Ranking as pair-level evidence", F["h2"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "Each dot is one train-to-test mission pair. The diamond is the tissue mean and the line is its bootstrap interval.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    plot_x0, plot_x1 = x0 + 430, x1 - 185
    plot_y0, plot_y1 = y0 + 275, y1 - 215
    sx = lambda v: int(plot_x0 + max(0.0, min(1.0, v)) * (plot_x1 - plot_x0))
    draw_axis_labels(draw, sx, plot_x0, plot_x1, plot_y0, plot_y1)

    row_gap = (plot_y1 - plot_y0) / len(rows)
    lanes = [-38, -22, -8, 8, 22, 38]
    for idx, row in enumerate(rows):
        cy = int(plot_y0 + row_gap * idx + row_gap * 0.52)
        accent = str(row["accent"])
        row_top = int(cy - row_gap * 0.37)
        row_bottom = int(cy + row_gap * 0.37)
        rounded(draw, (plot_x0, row_top, plot_x1, row_bottom), 24, "#121B28", rgba(accent, 70), 2)
        dashed_horizontal(draw, plot_x0, plot_x1, cy, "#233044", 2, 22)

        label_x = x0 + 54
        text(draw, (label_x, cy - 58), str(row["label"]), F["h3"], accent)
        text(draw, (label_x, cy - 15), f"mean {float(row['mean']):.3f}", F["small_bold"], COLORS["text"])
        text(draw, (label_x, cy + 22), f"{int(row['n_missions'])} missions / {int(row['n_pairs'])} pairs", F["tiny"], COLORS["muted"])

        ci_low = float(row["ci_low"])
        ci_high = float(row["ci_high"])
        mean = float(row["mean"])
        draw.line((sx(ci_low), cy + 52, sx(ci_high), cy + 52), fill=COLORS["axis"], width=6)
        draw.line((sx(ci_low), cy + 42, sx(ci_low), cy + 62), fill=COLORS["axis"], width=5)
        draw.line((sx(ci_high), cy + 42, sx(ci_high), cy + 62), fill=COLORS["axis"], width=5)

        pair_values = sorted([float(rec["auroc"]) for rec in row["pairs"]])  # type: ignore[index]
        for j, value in enumerate(pair_values):
            lane = lanes[j % len(lanes)]
            px = sx(value)
            py = cy + lane
            fill = score_color(value)
            draw.ellipse((px - 14, py - 14, px + 14, py + 14), fill=rgba(fill, 218), outline="#0B1018", width=3)

        mx = sx(mean)
        diamond = [(mx, cy + 24), (mx + 22, cy + 52), (mx, cy + 80), (mx - 22, cy + 52)]
        draw.polygon(diamond, fill=accent, outline=COLORS["bg"])
        text(draw, (x1 - 52, cy + 43), f"{float(row['min']):.2f}-{float(row['max']):.2f}", F["micro_bold"], COLORS["muted"], anchor="ra")

    text(draw, (x1 - 52, plot_y0 - 22), "pair range", F["micro"], COLORS["dim"], anchor="ra")
    legend_y = y1 - 150
    legend_items = [
        ("<0.50", COLORS["rose"]),
        ("0.50-0.69", COLORS["amber"]),
        ("0.70-0.84", COLORS["blue"]),
        (">=0.85", COLORS["teal"]),
    ]
    lx = x0 + 52
    for label, color in legend_items:
        draw.ellipse((lx, legend_y, lx + 26, legend_y + 26), fill=color)
        text(draw, (lx + 40, legend_y + 2), label, F["tiny"], COLORS["muted"])
        lx += 210
    text(draw, (x1 - 600, legend_y + 2), "diamond = mean; bar = 95% CI", F["tiny"], COLORS["muted"])


def draw_insight_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], rows: list[dict[str, object]]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 44, y0 + 45), "C. What the spread tells you", F["h2"])
    paragraph(
        draw,
        (x0 + 44, y0 + 98),
        "The mean is the ranking number; the pair cloud shows whether the signal is repeated or heterogeneous.",
        F["small"],
        x1 - x0 - 88,
        COLORS["muted"],
        7,
    )

    cards = [
        ("High-transfer rows", "Pair mass sits to the right of 0.70.", COLORS["teal"]),
        ("Middle rows", "Small mission inventories still show supportive pair clusters.", COLORS["blue"]),
        ("Heterogeneous rows", "Wide spread pulls the row mean toward chance.", COLORS["amber"]),
    ]
    cy = y0 + 300
    for idx, (title, body, accent) in enumerate(cards):
        rounded(draw, (x0 + 48, cy, x1 - 48, cy + 210), 28, "#151F2D", "#2A394D", 2)
        text(draw, (x0 + 84, cy + 34), title, F["h3"], accent)
        paragraph(draw, (x0 + 84, cy + 84), body, F["small"], x1 - x0 - 168, COLORS["muted"], 8)
        strip_x0, strip_x1 = x0 + 84, x1 - 84
        strip_y = cy + 164
        draw.line((strip_x0, strip_y, strip_x1, strip_y), fill="#2A3546", width=7)
        for value in [0.55, 0.68, 0.74, 0.88, 0.96] if idx == 0 else ([0.51, 0.66, 0.73, 0.81, 0.93] if idx == 1 else [0.19, 0.44, 0.61, 0.75, 1.0]):
            px = int(strip_x0 + value * (strip_x1 - strip_x0))
            draw.ellipse((px - 11, strip_y - 11, px + 11, strip_y + 11), fill=score_color(value), outline=COLORS["bg"], width=2)
        cy += 252

    callout = (x0 + 48, y1 - 345, x1 - 48, y1 - 56)
    rounded(draw, callout, 30, "#211E17", "#69532B", 2)
    text(draw, (callout[0] + 34, callout[1] + 40), "Slide role", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (callout[0] + 34, callout[1] + 96),
        "This slide turns the ranked tissue list into an auditable transfer matrix story before moving into pathway conservation.",
        F["small"],
        callout[2] - callout[0] - 68,
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
        "The matrix makes directionality visible: rows are train missions, columns are held-out test missions, and each cell is AUROC.",
        F["small"],
        footer[2] - footer[0] - 84,
        COLORS["muted"],
        6,
    )


def build_manifest(rows: list[dict[str, object]], image_path: Path) -> None:
    manifest_rows = []
    for row in rows:
        manifest_rows.append(
            {
                "tissue": row["id"],
                "label": row["label"],
                "mean": round(float(row["mean"]), 4),
                "ci": [round(float(row["ci_low"]), 4), round(float(row["ci_high"]), 4)],
                "n_missions": int(row["n_missions"]),
                "n_pairs": int(row["n_pairs"]),
                "pair_range": [round(float(row["min"]), 4), round(float(row["max"]), 4)],
                "matrix_source": str(Path(row["matrix_path"]).relative_to(ROOT)),
            }
        )
    manifest = {
        "asset": str(image_path.relative_to(ROOT)),
        "summary_source": str(SUMMARY.relative_to(ROOT)),
        "score_unit": "AUROC",
        "reference_line": 0.70,
        "rows": manifest_rows,
    }
    (OUT_DIR / "transfer_matrix_behind_ranking_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    rows = load_data()
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image, "RGBA")

    draw.rectangle((0, 0, W, 250), fill="#0B1018")
    text(draw, (140, 72), "TRANSFER MATRIX CHECK", F["kicker"], COLORS["teal"])
    text(draw, (140, 118), "The Transfer Matrix Behind The Ranking", F["title"], COLORS["text"])
    text(
        draw,
        (140, 212),
        "A tissue score is an average over directed mission-pair tests; the spread shows stable versus heterogeneous transfer.",
        F["subtitle"],
        COLORS["muted"],
    )
    draw_header(draw, rows)

    draw_matrix_grammar(draw, (130, 330, 1050, 1820), rows)
    draw_pair_spread(draw, (1085, 330, 2840, 1820), rows)
    draw_insight_panel(draw, (2875, 330, 3710, 1820), rows)
    draw_footer(draw)

    out = OUT_DIR / "transfer_matrix_behind_ranking_premium.png"
    image.save(out, quality=96)
    build_manifest(rows, out)
    print(json.dumps({"asset": str(out.relative_to(ROOT)), "manifest": str((OUT_DIR / "transfer_matrix_behind_ranking_manifest.json").relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

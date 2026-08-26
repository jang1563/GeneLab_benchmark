#!/usr/bin/env python3
"""Build the detailed-deck liver heterogeneity explainer asset."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
TRANSFER_SUMMARY = ROOT / "evaluation" / "B_cross_mission_summary.json"
NES_TRANSFER = ROOT / "evaluation" / "NES_conservation_vs_transfer.json"
MISSION_ID = ROOT / "v4" / "evaluation" / "D3_v4_liver_pca_lr_gene.json"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "liver_heterogeneity"
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

MISSION_COLORS = {
    "MHU-2": COLORS["violet"],
    "RR-1": COLORS["blue"],
    "RR-3": COLORS["sky"],
    "RR-6": COLORS["teal"],
    "RR-8": COLORS["green"],
    "RR-9": COLORS["amber"],
}

PASS_COUNTS = {
    "liver": {"pass": 13, "pairs": 30},
    "thymus": {"pass": 9, "pairs": 12},
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
    "h2": load_font(44, True),
    "h3": load_font(34, True),
    "body": load_font(30, False),
    "body_bold": load_font(30, True),
    "small": load_font(26, False),
    "small_bold": load_font(26, True),
    "tiny": load_font(22, False),
    "tiny_bold": load_font(22, True),
    "micro": load_font(18, False),
    "micro_bold": load_font(18, True),
    "stat": load_font(90, True),
    "mega": load_font(118, True),
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
    fill: str,
    leading: int = 8,
) -> int:
    x, y = xy
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def stat_badge(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, label: str, value: str, accent: str) -> None:
    rounded(draw, (x, y, x + w, y + 104), 26, "#121B28", "#2A394D", 2)
    text(draw, (x + 28, y + 24), label, F["tiny"], accent)
    text(draw, (x + 28, y + 56), value, F["small_bold"], COLORS["text"])


def load_data() -> dict[str, object]:
    transfer = json.loads(TRANSFER_SUMMARY.read_text())
    nes = json.loads(NES_TRANSFER.read_text())["data"]
    mission_id = json.loads(MISSION_ID.read_text())
    return {"transfer": transfer, "nes": nes, "mission_id": mission_id}


def draw_header(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    liver = data["transfer"]["liver"]  # type: ignore[index]
    liver_pca = liver["methods"]["pca_lr"]  # type: ignore[index]
    nes_liver = data["nes"]["liver"]  # type: ignore[index]
    badges = [
        ("MISSIONS", str(liver["n_missions"]), 230, COLORS["teal"]),
        ("DIRECTED PAIRS", str(liver_pca["n_valid_pairs"]), 330, COLORS["blue"]),
        ("MEAN AUROC", f"{liver_pca['mean_auroc']:.3f}", 300, COLORS["amber"]),
        ("NES MEAN r", f"{nes_liver['nes_mean_r']:.3f}", 300, COLORS["green"]),
    ]
    bx = 2320
    for label, value, width, accent in badges:
        stat_badge(draw, bx, 72, width, label, value, accent)
        bx += width + 28


def draw_inventory_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    liver = data["transfer"]["liver"]  # type: ignore[index]
    missions = list(liver["missions"])  # type: ignore[index]
    sizes = liver["mission_sizes"]  # type: ignore[index]
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "A. More missions create many tests", F["h2"])
    paragraph(
        draw,
        (x0 + 48, y0 + 100),
        "Each ordered train-to-test mission pair becomes a separate transfer question.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    cx, cy = (x0 + x1) // 2, y0 + 520
    radius = 260
    coords: dict[str, tuple[int, int]] = {}
    for i, mission in enumerate(missions):
        angle = -math.pi / 2 + i * 2 * math.pi / len(missions)
        coords[mission] = (int(cx + radius * math.cos(angle)), int(cy + radius * math.sin(angle)))

    for i, a in enumerate(missions):
        for b in missions[i + 1 :]:
            ax, ay = coords[a]
            bx, by = coords[b]
            draw.line((ax, ay, bx, by), fill=rgba("#4D5B70", 70), width=3)

    for mission in missions:
        mx, my = coords[mission]
        total = int(sizes[mission]["total"])
        color = MISSION_COLORS[mission]
        node_r = 34 + min(44, int(total / 3))
        draw.ellipse((mx - node_r, my - node_r, mx + node_r, my + node_r), fill=rgba(color, 210), outline=COLORS["white"], width=3)
        text(draw, (mx, my - 10), mission, F["tiny_bold"], "#091019", anchor="ma")
        text(draw, (mx, my + 18), f"n={total}", F["micro_bold"], "#091019", anchor="ma")

    formula = (x0 + 78, y1 - 300, x1 - 78, y1 - 72)
    rounded(draw, formula, 28, "#151F2D", "#2A394D", 2)
    text(draw, (formula[0] + 36, formula[1] + 34), "Directed-pair grammar", F["h3"], COLORS["blue"])
    text(draw, (formula[0] + 36, formula[1] + 92), "6 missions", F["body_bold"], COLORS["text"])
    text(draw, (formula[0] + 206, formula[1] + 92), "x", F["body_bold"], COLORS["muted"])
    text(draw, (formula[0] + 246, formula[1] + 92), "5 targets", F["body_bold"], COLORS["text"])
    text(draw, (formula[0] + 420, formula[1] + 92), "=", F["body_bold"], COLORS["muted"])
    text(draw, (formula[0] + 465, formula[1] + 82), "30", F["stat"], COLORS["amber"])
    text(draw, (formula[0] + 590, formula[1] + 98), "directed tests", F["small_bold"], COLORS["text"])
    paragraph(
        draw,
        (formula[0] + 36, formula[1] + 168),
        "The liver row samples a broader mission landscape than three- or four-mission tissues.",
        F["small"],
        formula[2] - formula[0] - 72,
        COLORS["muted"],
        7,
    )


def draw_pair_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    liver = data["transfer"]["liver"]  # type: ignore[index]
    thymus = data["transfer"]["thymus"]  # type: ignore[index]
    liver_pca = liver["methods"]["pca_lr"]  # type: ignore[index]
    thymus_pca = thymus["methods"]["pca_lr"]  # type: ignore[index]
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "B. The average compresses the pair mix", F["h2"])
    paragraph(
        draw,
        (x0 + 48, y0 + 100),
        "The Category B score is a mean over directed pairs, so inconsistent pairs pull the row toward chance.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    score_box = (x0 + 58, y0 + 220, x0 + 560, y0 + 500)
    rounded(draw, score_box, 26, "#151F2D", "#2A394D", 2)
    text(draw, (score_box[0] + 34, score_box[1] + 28), "Liver transfer", F["h3"], COLORS["amber"])
    text(draw, (score_box[0] + 34, score_box[1] + 86), f"{liver_pca['mean_auroc']:.3f}", F["mega"], COLORS["text"])
    text(
        draw,
        (score_box[0] + 34, score_box[1] + 208),
        f"95% CI {liver_pca['ci_low']:.3f}-{liver_pca['ci_high']:.3f}",
        F["small_bold"],
        COLORS["muted"],
    )

    compare_box = (x0 + 590, y0 + 220, x1 - 58, y0 + 500)
    rounded(draw, compare_box, 26, "#151F2D", "#2A394D", 2)
    text(draw, (compare_box[0] + 34, compare_box[1] + 28), "Thymus reference", F["h3"], COLORS["teal"])
    text(draw, (compare_box[0] + 34, compare_box[1] + 86), f"{thymus_pca['mean_auroc']:.3f}", F["mega"], COLORS["text"])
    text(draw, (compare_box[0] + 34, compare_box[1] + 208), "9/12 pairs at or above 0.70", F["small_bold"], COLORS["muted"])

    grid_box = (x0 + 58, y0 + 565, x1 - 58, y0 + 1015)
    rounded(draw, grid_box, 28, "#151F2D", "#2A394D", 2)
    text(draw, (grid_box[0] + 34, grid_box[1] + 30), "30 liver directed pairs", F["h3"], COLORS["blue"])
    text(draw, (grid_box[0] + 34, grid_box[1] + 78), "Each tile is one train-mission -> test-mission run.", F["small"], COLORS["muted"])
    tile = 56
    gap = 12
    gx = grid_box[0] + 36
    gy = grid_box[1] + 142
    for i in range(PASS_COUNTS["liver"]["pairs"]):
        row, col = divmod(i, 10)
        color = COLORS["teal"] if i < PASS_COUNTS["liver"]["pass"] else COLORS["amber"]
        alpha = 230 if i < PASS_COUNTS["liver"]["pass"] else 150
        rounded(draw, (gx + col * (tile + gap), gy + row * (tile + gap), gx + col * (tile + gap) + tile, gy + row * (tile + gap) + tile), 12, rgba(color, alpha), None, 0)
    legend_x = grid_box[0] + 760
    text(draw, (legend_x, gy + 4), "Pair-count readout", F["small_bold"], COLORS["text"])
    rounded(draw, (legend_x, gy + 58, legend_x + 220, gy + 100), 16, rgba(COLORS["teal"], 230), None, 0)
    text(draw, (legend_x + 250, gy + 64), "13/30 at >= 0.70", F["small_bold"], COLORS["text"])
    rounded(draw, (legend_x, gy + 122, legend_x + 220, gy + 164), 16, rgba(COLORS["amber"], 150), None, 0)
    text(draw, (legend_x + 250, gy + 128), "17/30 shape mean", F["small_bold"], COLORS["text"])

    callout = (x0 + 58, y1 - 300, x1 - 58, y1 - 58)
    rounded(draw, callout, 28, "#211E17", "#69532B", 2)
    text(draw, (callout[0] + 34, callout[1] + 34), "Reading move", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (callout[0] + 34, callout[1] + 92),
        "For a broad inventory tissue, the mean transfer score is a consistency test across many train-test directions.",
        F["small_bold"],
        callout[2] - callout[0] - 68,
        COLORS["text"],
        8,
    )
    paragraph(
        draw,
        (callout[0] + 34, callout[1] + 174),
        "The same AUROC grammar now asks whether liver response is shared across mission contexts.",
        F["small"],
        callout[2] - callout[0] - 68,
        COLORS["muted"],
        7,
    )


def matrix_color(value: float) -> tuple[int, int, int]:
    if value < 0:
        t = min(1.0, abs(value) / 0.45)
        return mix("#151F2D", COLORS["rose"], 0.35 + 0.65 * t)
    t = min(1.0, value / 0.45)
    return mix("#151F2D", COLORS["teal"], 0.25 + 0.75 * t)


def draw_nes_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    liver = data["transfer"]["liver"]  # type: ignore[index]
    missions = list(liver["missions"])  # type: ignore[index]
    nes_liver = data["nes"]["liver"]  # type: ignore[index]
    pairs = nes_liver["nes_pairs"]  # type: ignore[index]
    mission_id = data["mission_id"]  # type: ignore[assignment]
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "C. Liver carries visible mission context", F["h2"])
    paragraph(
        draw,
        (x0 + 48, y0 + 100),
        "Pathway-rank agreement varies across liver mission pairs, and mission identity is easy to recover.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    matrix_box = (x0 + 58, y0 + 230, x1 - 58, y0 + 875)
    rounded(draw, matrix_box, 28, "#151F2D", "#2A394D", 2)
    text(draw, (matrix_box[0] + 32, matrix_box[1] + 28), "NES conservation map", F["h3"], COLORS["green"])
    text(draw, (matrix_box[0] + 32, matrix_box[1] + 76), "Spearman r across 50 Hallmark pathways", F["small"], COLORS["muted"])

    cell = 72
    gap = 10
    mx = matrix_box[0] + 190
    my = matrix_box[1] + 150
    for i, mission in enumerate(missions):
        text(draw, (mx + i * (cell + gap) + cell / 2, my - 38), mission, F["micro_bold"], COLORS["axis"], anchor="ma")
        text(draw, (mx - 18, my + i * (cell + gap) + cell / 2 - 10), mission, F["micro_bold"], COLORS["axis"], anchor="ra")
    values: list[float] = []
    for i, a in enumerate(missions):
        for j, b in enumerate(missions):
            bx = mx + j * (cell + gap)
            by = my + i * (cell + gap)
            if i == j:
                rounded(draw, (bx, by, bx + cell, by + cell), 14, "#101823", "#2A394D", 1)
                text(draw, (bx + cell / 2, by + 25), "-", F["tiny_bold"], COLORS["dim"], anchor="ma")
                continue
            key = f"{a}_{b}" if f"{a}_{b}" in pairs else f"{b}_{a}"
            value = float(pairs[key])
            values.append(value)
            rounded(draw, (bx, by, bx + cell, by + cell), 14, matrix_color(value), None, 0)
            text(draw, (bx + cell / 2, by + 24), f"{value:+.2f}", F["micro_bold"], COLORS["text"], anchor="ma")

    unique_values = [float(v) for v in pairs.values()]
    neg = sum(1 for v in unique_values if v < 0)
    pos = sum(1 for v in unique_values if v > 0)
    stat_x = matrix_box[0] + 720
    stats = [
        ("mean r", f"{nes_liver['nes_mean_r']:.3f}", COLORS["green"]),
        ("range", f"{min(unique_values):+.3f} to {max(unique_values):+.3f}", COLORS["blue"]),
        ("pair signs", f"{pos} positive / {neg} negative", COLORS["amber"]),
    ]
    sy = matrix_box[1] + 150
    for label, value, color in stats:
        rounded(draw, (stat_x, sy, matrix_box[2] - 34, sy + 105), 20, "#101823", "#2A394D", 2)
        text(draw, (stat_x + 28, sy + 22), label, F["tiny"], color)
        text(draw, (stat_x + 28, sy + 54), value, F["small_bold"], COLORS["text"])
        sy += 126

    clue = (x0 + 58, y0 + 930, x1 - 58, y0 + 1165)
    rounded(draw, clue, 28, "#151F2D", "#2A394D", 2)
    text(draw, (clue[0] + 32, clue[1] + 32), "Mission-identity clue", F["h3"], COLORS["violet"])
    text(draw, (clue[0] + 36, clue[1] + 92), "macro-F1", F["small_bold"], COLORS["muted"])
    text(draw, (clue[0] + 36, clue[1] + 124), f"{mission_id['macro_f1']:.3f}", F["stat"], COLORS["text"])
    paragraph(
        draw,
        (clue[0] + 360, clue[1] + 58),
        "In the expanded liver panel, gene expression cleanly separates six mission labels.",
        F["small_bold"],
        clue[2] - clue[0] - 400,
        COLORS["text"],
        8,
    )
    paragraph(
        draw,
        (clue[0] + 360, clue[1] + 150),
        "That gives models a strong mission-context signal alongside the flight-vs-ground signal.",
        F["small"],
        clue[2] - clue[0] - 400,
        COLORS["muted"],
        7,
    )

    take = (x0 + 58, y1 - 300, x1 - 58, y1 - 58)
    rounded(draw, take, 28, "#211E17", "#69532B", 2)
    text(draw, (take[0] + 34, take[1] + 34), "Carry forward", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (take[0] + 34, take[1] + 92),
        "Liver becomes the heterogeneity example; the next slide generalizes this with NES conservation across tissues.",
        F["small_bold"],
        take[2] - take[0] - 68,
        COLORS["text"],
        8,
    )


def main() -> None:
    data = load_data()
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)

    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 48), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 42), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "LIVER HETEROGENEITY EXPLAINER", F["kicker"], COLORS["amber"])
    text(draw, (150, 122), "More liver missions expose lower transfer consistency", F["title"])
    text(
        draw,
        (150, 216),
        "The same Category B grammar averages many train-test directions, so mission-to-mission disagreement becomes visible.",
        F["subtitle"],
        COLORS["muted"],
    )
    draw_header(draw, data)

    draw_inventory_panel(draw, (150, 350, 1070, 1800), data)
    draw_pair_panel(draw, (1110, 350, 2450, 1800), data)
    draw_nes_panel(draw, (2490, 350, 3690, 1800), data)

    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    text(draw, (205, 1925), "Takeaway", F["small_bold"], COLORS["blue"])
    source = "More mission contexts expose why liver has lower transfer consistency than the top-ranked tissues."
    paragraph(draw, (390, 1925), source, F["small"], 3140, COLORS["muted"], 7)
    text(draw, (205, 1995), "Readout", F["small_bold"], COLORS["amber"])
    readout = "This slide explains how a broad mission inventory can reveal tissue-specific transfer heterogeneity."
    paragraph(draw, (390, 1995), readout, F["small"], 3140, COLORS["muted"], 7)

    png = OUT_DIR / "liver_mission_heterogeneity_premium.png"
    canvas.save(png, quality=95)
    gray = OUT_DIR / "liver_mission_heterogeneity_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    liver = data["transfer"]["liver"]  # type: ignore[index]
    liver_pca = liver["methods"]["pca_lr"]  # type: ignore[index]
    nes_liver = data["nes"]["liver"]  # type: ignore[index]
    mission_id = data["mission_id"]  # type: ignore[assignment]
    manifest = {
        "title": "More liver missions expose lower transfer consistency",
        "sources": [
            "evaluation/B_cross_mission_summary.json",
            "README.md",
            "evaluation/NES_conservation_vs_transfer.json",
            "v4/evaluation/D3_v4_liver_pca_lr_gene.json",
        ],
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "liver_transfer": {
            "missions": int(liver["n_missions"]),
            "directed_pairs": int(liver_pca["n_valid_pairs"]),
            "mean_auroc": round(float(liver_pca["mean_auroc"]), 3),
            "ci": [round(float(liver_pca["ci_low"]), 3), round(float(liver_pca["ci_high"]), 3)],
            "std_auroc": round(float(liver_pca["std_auroc"]), 3),
            "pairs_at_or_above_0_70": PASS_COUNTS["liver"]["pass"],
        },
        "liver_nes": {
            "mean_r": round(float(nes_liver["nes_mean_r"]), 3),
            "pair_count": len(nes_liver["nes_pairs"]),
            "min_r": round(min(float(v) for v in nes_liver["nes_pairs"].values()), 3),
            "max_r": round(max(float(v) for v in nes_liver["nes_pairs"].values()), 3),
        },
        "mission_identity": {
            "source": "v4/evaluation/D3_v4_liver_pca_lr_gene.json",
            "macro_f1": round(float(mission_id["macro_f1"]), 3),
            "n_samples": int(mission_id["n_samples"]),
        },
        "readout": "This slide explains how a broad mission inventory can reveal tissue-specific transfer heterogeneity.",
    }
    manifest_path = OUT_DIR / "liver_mission_heterogeneity_premium_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build the detailed-deck tissue transfer hierarchy asset."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
SOURCE_TRANSFER = ROOT / "evaluation" / "NES_conservation_vs_transfer.json"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "tissue_hierarchy"
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
    "violet": "#8C63F7",
    "white": "#FFFFFF",
}


CANONICAL_ROWS = [
    {"tissue": "thymus", "label": "Thymus", "pairs": 12, "missions": 4, "pass": 9, "tier": 1, "color": COLORS["teal"]},
    {"tissue": "gastrocnemius", "label": "Gastrocnemius", "pairs": 6, "missions": 3, "pass": 4, "tier": 1, "color": COLORS["green"]},
    {"tissue": "skin", "label": "Skin", "pairs": 6, "missions": 3, "pass": 5, "tier": 2, "color": COLORS["blue"]},
    {"tissue": "eye", "label": "Eye", "pairs": 6, "missions": 3, "pass": 5, "tier": 2, "color": COLORS["sky"]},
    {"tissue": "liver", "label": "Liver", "pairs": 30, "missions": 6, "pass": 13, "tier": 3, "color": COLORS["amber"]},
    {"tissue": "kidney", "label": "Kidney", "pairs": 6, "missions": 3, "pass": 2, "tier": 3, "color": COLORS["rose"]},
]


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
    "title": load_font(80, True),
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
    "stat": load_font(92, True),
    "mega": load_font(118, True),
}


def rounded(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    radius: int,
    fill: str,
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
        else:
            if current:
                lines.append(" ".join(current))
            current = [word]
    if current:
        lines.append(" ".join(current))
    return lines


def multiline(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    lines: Iterable[str],
    font: ImageFont.ImageFont,
    fill: str,
    leading: int = 8,
) -> int:
    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=font, fill=fill)
        y += font.size + leading
    return y


def paragraph(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    body: str,
    font: ImageFont.ImageFont,
    max_width: int,
    fill: str,
    leading: int = 8,
) -> int:
    return multiline(draw, xy, wrap_lines(draw, body, font, max_width), font, fill, leading)


def dashed_vertical(draw: ImageDraw.ImageDraw, x: int, y0: int, y1: int, color: str, width: int = 3, dash: int = 18) -> None:
    y = y0
    while y < y1:
        draw.line((x, y, x, min(y + dash, y1)), fill=color, width=width)
        y += dash * 2


def load_data() -> list[dict[str, object]]:
    transfer = json.loads(SOURCE_TRANSFER.read_text())["data"]
    rows: list[dict[str, object]] = []
    for row in CANONICAL_ROWS:
        entry = transfer[row["tissue"]]
        rows.append(
            {
                **row,
                "auroc": float(entry["transfer_auroc"]),
                "ci_lower": float(entry["transfer_ci"][0]),
                "ci_upper": float(entry["transfer_ci"][1]),
                "nes_mean_r": float(entry["nes_mean_r"]),
            }
        )
    return rows


def draw_header_badges(draw: ImageDraw.ImageDraw, rows: list[dict[str, object]]) -> None:
    total_pairs = sum(int(row["pairs"]) for row in rows)
    badges = [
        ("CATEGORY", "B transfer", 320),
        ("TISSUES", str(len(rows)), 240),
        ("DIRECTED PAIRS", str(total_pairs), 350),
        ("LEAD ROW", "thymus 0.860", 390),
    ]
    bx = 2225
    for kicker, body, badge_w in badges:
        rounded(draw, (bx, 72, bx + badge_w, 176), 26, "#121B28", "#2A394D", 2)
        text(draw, (bx + 28, 96), kicker, F["tiny"], COLORS["teal"])
        text(draw, (bx + 28, 126), body, F["small_bold"], COLORS["text"])
        bx += badge_w + 30


def draw_reading_guide(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "A. Read every row the same way", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "The worked thymus example becomes a ranking when the same score grammar is applied to all tissues.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    cards = [
        ("dot", "mean transfer AUROC", COLORS["teal"]),
        ("line", "bootstrap interval", COLORS["axis"]),
        ("0.50", "chance reference", COLORS["amber"]),
        ("0.70", "transfer reference", COLORS["blue"]),
        ("pass", "pairs at or above 0.70", COLORS["green"]),
    ]
    cy = y0 + 240
    for title, body, color in cards:
        rounded(draw, (x0 + 58, cy, x1 - 58, cy + 130), 24, "#151F2D", "#2A394D", 2)
        rounded(draw, (x0 + 92, cy + 32, x0 + 188, cy + 88), 18, "#101823", color, 2)
        text(draw, (x0 + 140, cy + 47), title, F["small_bold"], color, anchor="ma")
        text(draw, (x0 + 220, cy + 36), body, F["small_bold"], COLORS["text"])
        cy += 158

    callout = (x0 + 58, y1 - 310, x1 - 58, y1 - 56)
    rounded(draw, callout, 28, "#211E17", "#69532B", 2)
    text(draw, (callout[0] + 34, callout[1] + 36), "Reading move", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (callout[0] + 34, callout[1] + 92),
        "Rows near the top combine high score, supportive interval, and repeated pair-level transfer.",
        F["small_bold"],
        callout[2] - callout[0] - 68,
        COLORS["text"],
        8,
    )
    paragraph(
        draw,
        (callout[0] + 34, callout[1] + 190),
        "That makes the hierarchy a tissue-specific generalization map.",
        F["small"],
        callout[2] - callout[0] - 68,
        COLORS["muted"],
        7,
    )


def draw_forest_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], rows: list[dict[str, object]]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "B. Cross-mission transfer hierarchy", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "Tissues separate into high, intermediate, and low-transfer groups under the same directed-pair protocol.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    axis_x0, axis_x1 = x0 + 460, x1 - 335
    row_y0 = y0 + 355
    row_gap = 185
    min_v, max_v = 0.35, 1.00

    def sx(value: float) -> int:
        return int(axis_x0 + (value - min_v) / (max_v - min_v) * (axis_x1 - axis_x0))

    for tick in [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        tx = sx(tick)
        draw.line((tx, row_y0 - 120, tx, row_y0 + row_gap * (len(rows) - 1) + 94), fill=rgba(COLORS["grid"], 135), width=2)
        text(draw, (tx, row_y0 + row_gap * len(rows) - 30), f"{tick:.1f}", F["tiny"], COLORS["axis"], anchor="ma")
    dashed_vertical(draw, sx(0.5), row_y0 - 132, row_y0 + row_gap * (len(rows) - 1) + 96, COLORS["amber"], 4, 18)
    draw.line((sx(0.7), row_y0 - 132, sx(0.7), row_y0 + row_gap * (len(rows) - 1) + 96), fill=COLORS["blue"], width=4)
    text(draw, (sx(0.5), row_y0 - 170), "chance", F["tiny_bold"], COLORS["amber"], anchor="ma")
    text(draw, (sx(0.7), row_y0 - 170), "0.70", F["tiny_bold"], COLORS["blue"], anchor="ma")

    tier_bg = {
        1: rgba(COLORS["teal"], 28),
        2: rgba(COLORS["blue"], 24),
        3: rgba(COLORS["amber"], 22),
    }
    for i, row in enumerate(rows):
        y = row_y0 + i * row_gap
        color = str(row["color"])
        rounded(draw, (x0 + 44, y - 76, x1 - 44, y + 82), 22, tier_bg[int(row["tier"])], "#26364A", 1)
        text(draw, (x0 + 72, y - 52), str(row["label"]), F["h3"], COLORS["text"])
        text(draw, (x0 + 72, y - 6), f"{row['missions']} missions | {row['pairs']} pairs", F["small"], COLORS["text"])
        rounded(draw, (x0 + 270, y + 26, x0 + 405, y + 68), 15, "#101823", color, 2)
        text(draw, (x0 + 338, y + 37), f"TIER {row['tier']}", F["micro"], color, anchor="ma")

        ci0, ci1 = sx(float(row["ci_lower"])), sx(float(row["ci_upper"]))
        dot = sx(float(row["auroc"]))
        draw.line((axis_x0, y, axis_x1, y), fill="#263344", width=10)
        draw.line((ci0, y, ci1, y), fill=COLORS["axis"], width=10)
        draw.line((ci0, y - 22, ci0, y + 22), fill=COLORS["axis"], width=5)
        draw.line((ci1, y - 22, ci1, y + 22), fill=COLORS["axis"], width=5)
        draw.ellipse((dot - 25, y - 25, dot + 25, y + 25), fill=color, outline=COLORS["white"], width=3)

        stat_box = (x1 - 300, y - 58, x1 - 58, y + 58)
        rounded(draw, stat_box, 18, "#0F1722", "#2A394D", 1)
        text(draw, (stat_box[0] + 22, y - 43), "AUROC", F["micro"], COLORS["muted"])
        text(draw, (stat_box[0] + 22, y - 16), f"{float(row['auroc']):.3f}", F["body_bold"], COLORS["text"])
        ci_label = f"CI {float(row['ci_lower']):.3f}-{float(row['ci_upper']):.3f}"
        text(draw, (stat_box[0] + 22, y + 22), ci_label, F["tiny"], COLORS["axis"])

    text(draw, ((axis_x0 + axis_x1) / 2, y1 - 72), "Category B mean transfer AUROC", F["small_bold"], COLORS["muted"], anchor="ma")


def draw_summary_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], rows: list[dict[str, object]]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "C. What the hierarchy says", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "The ranking separates reactive tissues from tissues with broader mission heterogeneity.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    tier_cards = [
        ("Tier 1", "high transfer", "Thymus + gastrocnemius", COLORS["teal"]),
        ("Tier 2", "intermediate transfer", "Skin + eye", COLORS["blue"]),
        ("Tier 3", "heterogeneous transfer", "Liver + kidney", COLORS["amber"]),
    ]
    cy = y0 + 225
    for tier, label, tissues, color in tier_cards:
        rounded(draw, (x0 + 58, cy, x1 - 58, cy + 150), 24, "#151F2D", "#2A394D", 2)
        text(draw, (x0 + 92, cy + 26), tier, F["h3"], color)
        text(draw, (x0 + 92, cy + 72), label, F["small_bold"], COLORS["text"])
        text(draw, (x0 + 92, cy + 108), tissues, F["small"], COLORS["muted"])
        cy += 180

    pass_box = (x0 + 58, y0 + 805, x1 - 58, y0 + 1138)
    rounded(draw, pass_box, 28, "#151F2D", "#2A394D", 2)
    text(draw, (pass_box[0] + 32, pass_box[1] + 32), "Pair consistency", F["h3"], COLORS["green"])
    bar_x0, bar_x1 = pass_box[0] + 32, pass_box[2] - 32
    bar_y = pass_box[1] + 102
    for row in rows:
        frac = int(row["pass"]) / int(row["pairs"])
        color = str(row["color"])
        text(draw, (bar_x0, bar_y - 5), str(row["label"])[:12], F["tiny"], COLORS["text"])
        rounded(draw, (bar_x0 + 190, bar_y, bar_x1 - 90, bar_y + 26), 13, "#101823", "#2A394D", 1)
        rounded(draw, (bar_x0 + 190, bar_y, bar_x0 + 190 + int((bar_x1 - bar_x0 - 280) * frac), bar_y + 26), 13, color, None, 0)
        text(draw, (bar_x1, bar_y - 5), f"{row['pass']}/{row['pairs']}", F["tiny"], COLORS["muted"], anchor="ra")
        bar_y += 38

    take = (x0 + 58, y1 - 310, x1 - 58, y1 - 56)
    rounded(draw, take, 28, "#211E17", "#69532B", 2)
    text(draw, (take[0] + 34, take[1] + 36), "Carry forward", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (take[0] + 34, take[1] + 92),
        "Thymus becomes the lead worked example; liver becomes the heterogeneity contrast.",
        F["small_bold"],
        take[2] - take[0] - 68,
        COLORS["text"],
        8,
    )
    paragraph(
        draw,
        (take[0] + 34, take[1] + 188),
        "Next: explain why the largest mission inventory yields a lower transfer row.",
        F["small"],
        take[2] - take[0] - 68,
        COLORS["muted"],
        7,
    )


def main() -> None:
    rows = load_data()
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)

    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 48), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 42), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "TISSUE TRANSFER HIERARCHY", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "Thymus leads the cross-mission transfer hierarchy", F["title"])
    text(
        draw,
        (150, 216),
        "Six tissues are ranked with the same Category B grammar: mean AUROC, uncertainty, and mission-pair consistency.",
        F["subtitle"],
        COLORS["muted"],
    )
    draw_header_badges(draw, rows)

    draw_reading_guide(draw, (150, 350, 910, 1800))
    draw_forest_panel(draw, (955, 350, 2600, 1800), rows)
    draw_summary_panel(draw, (2645, 350, 3690, 1800), rows)

    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    text(draw, (205, 1925), "Takeaway", F["small_bold"], COLORS["blue"])
    source = "Ranking tissues with one score grammar shows where cross-mission transfer is strong and where it becomes fragile."
    paragraph(draw, (390, 1925), source, F["small"], 3140, COLORS["muted"], 7)
    text(draw, (205, 1995), "Next", F["small_bold"], COLORS["amber"])
    scope = "Liver heterogeneity and pathway conservation explain why the hierarchy is not just a leaderboard."
    paragraph(draw, (390, 1995), scope, F["small"], 3140, COLORS["muted"], 7)

    png = OUT_DIR / "tissue_transfer_hierarchy_premium.png"
    canvas.save(png, quality=95)
    gray = OUT_DIR / "tissue_transfer_hierarchy_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    manifest = {
        "title": "Thymus leads the cross-mission transfer hierarchy",
        "sources": [
            "evaluation/RESULTS_SUMMARY.md",
            "README.md",
            "evaluation/NES_conservation_vs_transfer.json",
        ],
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "rows": [
            {
                "tissue": str(row["label"]),
                "auroc": round(float(row["auroc"]), 3),
                "ci": [round(float(row["ci_lower"]), 3), round(float(row["ci_upper"]), 3)],
                "missions": int(row["missions"]),
                "pairs": int(row["pairs"]),
                "pairs_at_or_above_0_70": int(row["pass"]),
                "tier": int(row["tier"]),
            }
            for row in rows
        ],
        "comparisons": {
            "thymus_minus_liver": round(float(rows[0]["auroc"]) - float(rows[4]["auroc"]), 3),
            "thymus_vs_liver_p": 0.001,
            "gastrocnemius_vs_liver_p": 0.048,
            "skin_vs_liver_p": 0.032,
        },
        "scope": "This slide ranks tissue-level cross-mission transfer; companion slides explain liver heterogeneity and pathway conservation.",
    }
    manifest_path = OUT_DIR / "tissue_transfer_hierarchy_premium_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

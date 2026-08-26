#!/usr/bin/env python3
"""Build the detailed-deck worked tissue-score example asset."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
SOURCE_B4 = ROOT / "evaluation" / "submission_PCA-LR_baseline_v1_B4_eval.json"
SOURCE_TRANSFER = ROOT / "evaluation" / "NES_conservation_vs_transfer.json"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "worked_tissue_score"
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


def load_data() -> dict:
    b4 = json.loads(SOURCE_B4.read_text())
    transfer = json.loads(SOURCE_TRANSFER.read_text())["data"]
    folds = []
    for rec in b4["fold_results"]:
        _, train, test = rec["fold"].split("_", 2)
        folds.append(
            {
                "label": f"{train} -> {test}",
                "train": train,
                "test": test,
                "auroc": float(rec["auroc"]),
                "ci_lower": float(rec["ci_lower"]),
                "ci_upper": float(rec["ci_upper"]),
                "perm_p": float(rec["perm_p"]),
                "n_test": int(rec["n_test"]),
                "n_flight": int(rec["n_flight"]),
                "n_ground": int(rec["n_ground"]),
            }
        )
    thymus = transfer["thymus"]
    liver = transfer["liver"]
    return {
        "folds": folds,
        "missions": ["MHU-1", "MHU-2", "RR-6", "RR-9"],
        "mean_auroc": float(thymus["transfer_auroc"]),
        "ci": [float(v) for v in thymus["transfer_ci"]],
        "n_pairs": len(folds),
        "n_pass": sum(1 for rec in folds if rec["auroc"] >= 0.70),
        "liver_auroc": float(liver["transfer_auroc"]),
        "delta_vs_liver": float(thymus["transfer_auroc"]) - float(liver["transfer_auroc"]),
        "thymus_vs_liver_p": 0.001,
    }


def score_color(value: float) -> str:
    if value >= 0.90:
        return COLORS["teal"]
    if value >= 0.70:
        return COLORS["blue"]
    if value >= 0.50:
        return COLORS["amber"]
    return COLORS["rose"]


def draw_header_badges(draw: ImageDraw.ImageDraw, data: dict) -> None:
    badges = [
        ("TISSUE", "thymus", 280),
        ("MISSIONS", "4", 220),
        ("DIRECTED PAIRS", str(data["n_pairs"]), 340),
        ("TRANSFER AUROC", f"{data['mean_auroc']:.3f}", 370),
    ]
    bx = 2190
    for kicker, body, badge_w in badges:
        rounded(draw, (bx, 72, bx + badge_w, 176), 26, "#121B28", "#2A394D", 2)
        text(draw, (bx + 28, 96), kicker, F["tiny"], COLORS["teal"])
        text(draw, (bx + 28, 126), body, F["small_bold"], COLORS["text"])
        bx += badge_w + 30


def draw_matrix_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "A. One tissue row expands into pairs", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "Train on one thymus mission, test on a different thymus mission, then summarize the directed pair scores.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    missions = data["missions"]
    by_pair = {(rec["train"], rec["test"]): rec for rec in data["folds"]}
    cell = 145
    gap = 12
    mx0 = x0 + 250
    my0 = y0 + 355
    text(draw, (mx0 + cell * 2 + gap * 2, my0 - 126), "test mission", F["small_bold"], COLORS["sky"], anchor="ma")
    y_label = Image.new("RGBA", (300, 62), (0, 0, 0, 0))
    yd = ImageDraw.Draw(y_label)
    yd.text((150, 10), "train mission", font=F["small_bold"], fill=COLORS["teal"], anchor="ma")
    y_label = y_label.rotate(270, expand=True)
    draw.bitmap((x0 + 54, my0 + 188), y_label, fill=COLORS["teal"])

    for i, mission in enumerate(missions):
        text(draw, (mx0 + i * (cell + gap) + cell / 2, my0 - 58), mission, F["tiny_bold"], COLORS["text"], anchor="ma")
        text(draw, (mx0 - 34, my0 + i * (cell + gap) + cell / 2 - 12), mission, F["tiny_bold"], COLORS["text"], anchor="ra")

    for r, train in enumerate(missions):
        for c, test_mission in enumerate(missions):
            x = mx0 + c * (cell + gap)
            y = my0 + r * (cell + gap)
            if train == test_mission:
                rounded(draw, (x, y, x + cell, y + cell), 22, "#0F1723", "#26364A", 1)
                text(draw, (x + cell / 2, y + 54), "same", F["micro"], COLORS["dim"], anchor="ma")
                text(draw, (x + cell / 2, y + 82), "mission", F["micro"], COLORS["dim"], anchor="ma")
                continue
            rec = by_pair[(train, test_mission)]
            color = score_color(rec["auroc"])
            rounded(draw, (x, y, x + cell, y + cell), 22, rgba(color, 210), None, 0)
            rounded(draw, (x, y, x + cell, y + cell), 22, "#D8E9F8" if rec["auroc"] >= 0.70 else "#2A394D", 2)
            label_fill = COLORS["bg"] if rec["auroc"] >= 0.70 else COLORS["text"]
            text(draw, (x + cell / 2, y + 42), f"{rec['auroc']:.2f}", F["body_bold"], label_fill, anchor="ma")
            text(draw, (x + cell / 2, y + 88), f"n={rec['n_test']}", F["tiny_bold"], label_fill, anchor="ma")

    legend_y = y0 + 1085
    legend = [
        ("0.90+", COLORS["teal"]),
        ("0.70-0.89", COLORS["blue"]),
        ("0.50-0.69", COLORS["amber"]),
        ("<0.50", COLORS["rose"]),
    ]
    lx = x0 + 90
    for label, color in legend:
        rounded(draw, (lx, legend_y, lx + 170, legend_y + 52), 16, rgba(color, 220), None, 0)
        text(draw, (lx + 85, legend_y + 14), label, F["tiny_bold"], COLORS["bg"], anchor="ma")
        lx += 190

    callout = (x0 + 54, y1 - 260, x1 - 54, y1 - 56)
    rounded(draw, callout, 28, "#211E17", "#69532B", 2)
    text(draw, (callout[0] + 34, callout[1] + 36), "Reading move", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (callout[0] + 34, callout[1] + 90),
        "The matrix reveals the row behind the average: most thymus mission pairs clear the 0.70 transfer reference.",
        F["small_bold"],
        callout[2] - callout[0] - 68,
        COLORS["text"],
        8,
    )


def draw_pair_evidence_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "B. Pair evidence explains the row", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "Each dot is one train mission -> test mission transfer score; the line is the bootstrap interval.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    chart_x0, chart_x1 = x0 + 345, x1 - 88
    row_y0 = y0 + 315
    row_gap = 82
    min_v, max_v = 0.0, 1.0

    def sx(value: float) -> int:
        return int(chart_x0 + (value - min_v) / (max_v - min_v) * (chart_x1 - chart_x0))

    for tick in [0.0, 0.5, 0.7, 1.0]:
        tx = sx(tick)
        draw.line((tx, row_y0 - 95, tx, row_y0 + row_gap * 12 - 18), fill=rgba(COLORS["grid"], 140), width=2)
        text(draw, (tx, row_y0 + row_gap * 12 + 12), f"{tick:.1f}", F["tiny"], COLORS["axis"], anchor="ma")
    dashed_vertical(draw, sx(0.5), row_y0 - 104, row_y0 + row_gap * 12 - 20, COLORS["amber"], 4, 18)
    draw.line((sx(0.7), row_y0 - 104, sx(0.7), row_y0 + row_gap * 12 - 20), fill=COLORS["teal"], width=4)
    text(draw, (sx(0.5), row_y0 - 138), "chance", F["tiny_bold"], COLORS["amber"], anchor="ma")
    text(draw, (sx(0.7), row_y0 - 138), "0.70", F["tiny_bold"], COLORS["teal"], anchor="ma")

    for i, rec in enumerate(data["folds"]):
        y = row_y0 + i * row_gap
        color = score_color(rec["auroc"])
        text(draw, (x0 + 58, y - 21), rec["label"], F["tiny_bold"], COLORS["text"])
        text(draw, (x0 + 58, y + 7), f"{rec['n_flight']} FLT / {rec['n_ground']} GC", F["micro"], COLORS["muted"])
        draw.line((chart_x0, y, chart_x1, y), fill="#263344", width=6)
        draw.line((sx(rec["ci_lower"]), y, sx(rec["ci_upper"]), y), fill=COLORS["axis"], width=8)
        draw.line((sx(rec["ci_lower"]), y - 15, sx(rec["ci_lower"]), y + 15), fill=COLORS["axis"], width=4)
        draw.line((sx(rec["ci_upper"]), y - 15, sx(rec["ci_upper"]), y + 15), fill=COLORS["axis"], width=4)
        px = sx(rec["auroc"])
        draw.ellipse((px - 20, y - 20, px + 20, y + 20), fill=color, outline=COLORS["white"], width=3)
        text(draw, (x1 - 46, y - 18), f"{rec['auroc']:.3f}", F["tiny_bold"], COLORS["text"], anchor="ra")

    strip = (x0 + 54, y1 - 150, x1 - 54, y1 - 56)
    rounded(draw, strip, 24, "#151F2D", "#2A394D", 2)
    text(draw, (strip[0] + 34, strip[1] + 30), "Pair count", F["small_bold"], COLORS["teal"])
    text(draw, (strip[0] + 205, strip[1] + 30), f"{data['n_pass']}/{data['n_pairs']} directed pairs are at or above AUROC 0.70.", F["small_bold"], COLORS["text"])


def draw_readout_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "C. Tissue-level readout", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "The tissue row becomes a compact score after the mission-pair evidence is summarized.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    score_box = (x0 + 60, y0 + 225, x1 - 60, y0 + 555)
    rounded(draw, score_box, 32, "#15231F", COLORS["teal"], 3)
    text(draw, (score_box[0] + 46, score_box[1] + 38), "Transfer AUROC", F["h3"], COLORS["teal"])
    text(draw, (score_box[0] + 46, score_box[1] + 90), f"{data['mean_auroc']:.3f}", F["mega"], COLORS["text"])
    text(draw, (score_box[0] + 50, score_box[1] + 230), f"95% CI [{data['ci'][0]:.3f}, {data['ci'][1]:.3f}]", F["body_bold"], COLORS["muted"])

    cards = [
        (f"{data['n_pass']}/{data['n_pairs']}", "pairs at or above 0.70", COLORS["blue"]),
        (f"+{data['delta_vs_liver']:.3f}", "vs liver transfer AUROC", COLORS["amber"]),
        ("p=0.001", "thymus vs liver test", COLORS["green"]),
    ]
    cy = y0 + 620
    for value, label, color in cards:
        rounded(draw, (x0 + 60, cy, x1 - 60, cy + 140), 24, "#151F2D", "#2A394D", 2)
        text(draw, (x0 + 92, cy + 8), value, F["stat"], color)
        text(draw, (x0 + 92, cy + 96), label, F["small_bold"], COLORS["text"])
        cy += 165

    take = (x0 + 60, y1 - 310, x1 - 60, y1 - 56)
    rounded(draw, take, 28, "#211E17", "#69532B", 2)
    text(draw, (take[0] + 34, take[1] + 36), "Carry forward", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (take[0] + 34, take[1] + 92),
        "Thymus enters the hierarchy as a high-transfer tissue. Pair spread still shows which mission contexts soften transfer.",
        F["small_bold"],
        take[2] - take[0] - 68,
        COLORS["text"],
        8,
    )
    paragraph(
        draw,
        (take[0] + 34, take[1] + 188),
        "Next: rank every tissue with the same score grammar.",
        F["small"],
        take[2] - take[0] - 68,
        COLORS["muted"],
        7,
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

    text(draw, (150, 74), "WORKED RESULT EXAMPLE", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "Read one tissue score as mission-pair evidence", F["title"])
    text(
        draw,
        (150, 216),
        "Before ranking tissues, unpack the thymus row into directed train-to-test mission transfers.",
        F["subtitle"],
        COLORS["muted"],
    )
    draw_header_badges(draw, data)

    draw_matrix_panel(draw, (150, 350, 1135, 1800), data)
    draw_pair_evidence_panel(draw, (1185, 350, 2630, 1800), data)
    draw_readout_panel(draw, (2680, 350, 3690, 1800), data)

    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    text(draw, (205, 1925), "Takeaway", F["small_bold"], COLORS["blue"])
    source = "A tissue score is easiest to trust when the audience can see the train-to-test mission pairs behind the mean."
    paragraph(draw, (390, 1925), source, F["small"], 3140, COLORS["muted"], 7)
    text(draw, (205, 1995), "Next", F["small_bold"], COLORS["amber"])
    scope = "The next slides rank all tissues, explain liver heterogeneity, and connect transfer to pathway biology."
    paragraph(draw, (390, 1995), scope, F["small"], 3140, COLORS["muted"], 7)

    png = OUT_DIR / "worked_tissue_score_thymus_transfer_premium.png"
    canvas.save(png, quality=95)
    gray = OUT_DIR / "worked_tissue_score_thymus_transfer_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    manifest = {
        "title": "Read one tissue score as mission-pair evidence",
        "sources": [
            "evaluation/submission_PCA-LR_baseline_v1_B4_eval.json",
            "evaluation/NES_conservation_vs_transfer.json",
            "evaluation/RESULTS_SUMMARY.md",
        ],
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "summary": {
            "tissue": "thymus",
            "transfer_auroc": round(data["mean_auroc"], 3),
            "transfer_ci": [round(v, 3) for v in data["ci"]],
            "n_missions": len(data["missions"]),
            "n_pairs": data["n_pairs"],
            "pairs_at_or_above_0_70": data["n_pass"],
            "liver_transfer_auroc": round(data["liver_auroc"], 3),
            "delta_vs_liver": round(data["delta_vs_liver"], 3),
            "thymus_vs_liver_p": data["thymus_vs_liver_p"],
        },
        "readout_frame": "This slide teaches a tissue transfer row; companion slides rank all tissues and connect the result to pathway biology.",
    }
    manifest_path = OUT_DIR / "worked_tissue_score_thymus_transfer_premium_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

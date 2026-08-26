#!/usr/bin/env python3
"""Build the detailed-deck AUROC / uncertainty metric-primer asset."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "metric_primer"
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
    "violet": "#9C7CFF",
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
    "micro": load_font(19, False),
    "metric": load_font(88, True),
    "score": load_font(54, True),
    "gate": load_font(40, True),
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


def dashed_horizontal(draw: ImageDraw.ImageDraw, x0: int, y: int, x1: int, color: str, width: int = 3, dash: int = 18) -> None:
    x = x0
    while x < x1:
        draw.line((x, y, min(x + dash, x1), y), fill=color, width=width)
        x += dash * 2


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: str, width: int = 5) -> None:
    x0, y0 = start
    x1, y1 = end
    draw.line((x0, y0, x1, y1), fill=color, width=width)
    draw.polygon([(x1, y1), (x1 - 22, y1 - 13), (x1 - 22, y1 + 13)], fill=color)


def load_examples() -> list[dict[str, object]]:
    a4_holdout = json.loads((ROOT / "evaluation" / "A4_holdout_results.json").read_text())["lr"]
    a5_holdout = json.loads((ROOT / "evaluation" / "A5_holdout_results.json").read_text())
    a3 = json.loads((ROOT / "evaluation" / "A3_baseline_results.json").read_text())
    a3_rf_rr7 = next(f for f in a3["rf"]["folds"] if f["test_mission"] == "RR-7")
    return [
        {
            "label": "A5 skin RR-7 LR",
            "auroc": a5_holdout["lr"]["auroc"],
            "ci_lower": a5_holdout["lr"]["ci_lower"],
            "ci_upper": a5_holdout["lr"]["ci_upper"],
            "p": a5_holdout["lr"]["perm_pvalue"],
            "tag": "passes all gates",
            "color": COLORS["teal"],
        },
        {
            "label": "A4 thymus RR-23 LR",
            "auroc": a4_holdout["auroc"],
            "ci_lower": a4_holdout["ci_lower"],
            "ci_upper": a4_holdout["ci_upper"],
            "p": a4_holdout["perm_pvalue"],
            "tag": "passes all gates",
            "color": COLORS["green"],
        },
        {
            "label": "A5 skin RR-7 RF",
            "auroc": a5_holdout["rf"]["auroc"],
            "ci_lower": a5_holdout["rf"]["ci_lower"],
            "ci_upper": a5_holdout["rf"]["ci_upper"],
            "p": a5_holdout["rf"]["perm_pvalue"],
            "tag": "CI check active",
            "color": COLORS["amber"],
        },
        {
            "label": "A3 kidney RR-7 RF",
            "auroc": a3_rf_rr7["auroc"],
            "ci_lower": a3_rf_rr7["ci_lower"],
            "ci_upper": a3_rf_rr7["ci_upper"],
            "p": a3_rf_rr7["perm_pvalue"],
            "tag": "review zone",
            "color": COLORS["rose"],
        },
    ]


def fmt_p(value: float) -> str:
    if value < 0.001:
        return "<0.001"
    return f"{value:.3f}"


def draw_header_badges(draw: ImageDraw.ImageDraw) -> None:
    badges = [
        ("CHANCE", "AUROC 0.50", 330),
        ("UNCERTAINTY", "bootstrap CI", 370),
        ("DECISION", "3-gate GO", 330),
    ]
    bx = 2540
    for kicker, body, badge_w in badges:
        rounded(draw, (bx, 72, bx + badge_w, 176), 26, "#121B28", "#2A394D", 2)
        text(draw, (bx + 28, 96), kicker, F["tiny"], COLORS["sky"])
        text(draw, (bx + 28, 126), body, F["small_bold"], COLORS["text"])
        bx += badge_w + 30


def draw_rank_strip(
    draw: ImageDraw.ImageDraw,
    x0: int,
    y: int,
    w: int,
    title: str,
    ground_scores: list[float],
    flight_scores: list[float],
    auroc: str,
    color: str,
) -> None:
    rounded(draw, (x0, y, x0 + w, y + 260), 24, "#151F2D", "#2A394D", 2)
    text(draw, (x0 + 28, y + 28), title, F["small_bold"], COLORS["text"])
    text(draw, (x0 + w - 28, y + 28), auroc, F["small_bold"], color, anchor="ra")
    axis_x0 = x0 + 55
    axis_x1 = x0 + w - 55
    axis_y = y + 164
    draw.line((axis_x0, axis_y, axis_x1, axis_y), fill=COLORS["grid"], width=8)
    text(draw, (axis_x0 + 42, axis_y + 68), "lower score", F["tiny"], COLORS["dim"], anchor="la")
    text(draw, (axis_x1 - 42, axis_y + 68), "higher score", F["tiny"], COLORS["dim"], anchor="ra")
    for score in ground_scores:
        x = int(axis_x0 + score * (axis_x1 - axis_x0))
        draw.ellipse((x - 15, axis_y - 43, x + 15, axis_y - 13), fill=COLORS["dim"], outline="#B8C3D0", width=2)
    for score in flight_scores:
        x = int(axis_x0 + score * (axis_x1 - axis_x0))
        draw.ellipse((x - 17, axis_y + 14, x + 17, axis_y + 48), fill=color, outline="#DDF7F4", width=2)
    text(draw, (x0 + 54, y + 83), "Ground", F["tiny_bold"], COLORS["dim"])
    text(draw, (x0 + 54, y + 186), "Flight", F["tiny_bold"], color)


def draw_auroc_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "A. What AUROC asks", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "AUROC reads how well the model ranks Flight samples above Ground samples.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    formula = (x0 + 54, y0 + 205, x1 - 54, y0 + 380)
    rounded(draw, formula, 26, "#171D25", "#2A394D", 2)
    text(draw, (formula[0] + 35, formula[1] + 35), "AUROC", F["h3"], COLORS["teal"])
    text(draw, (formula[0] + 208, formula[1] + 34), "=", F["h3"], COLORS["muted"])
    text(draw, (formula[0] + 260, formula[1] + 28), "P(Flight score > Ground score)", F["body_bold"], COLORS["text"])
    paragraph(
        draw,
        (formula[0] + 35, formula[1] + 93),
        "A random Flight/Ground pair is compared by model score; more correct rank order pushes AUROC upward.",
        F["tiny"],
        formula[2] - formula[0] - 70,
        COLORS["muted"],
        5,
    )

    draw_rank_strip(
        draw,
        x0 + 54,
        y0 + 455,
        x1 - x0 - 108,
        "Mixed ranks",
        [0.16, 0.32, 0.50, 0.63, 0.80],
        [0.20, 0.38, 0.57, 0.70, 0.84],
        "AUROC approx 0.50",
        COLORS["amber"],
    )
    draw_rank_strip(
        draw,
        x0 + 54,
        y0 + 770,
        x1 - x0 - 108,
        "Separated ranks",
        [0.11, 0.23, 0.34, 0.46, 0.55],
        [0.58, 0.68, 0.77, 0.86, 0.94],
        "AUROC approx 0.90",
        COLORS["teal"],
    )

    insight = (x0 + 54, y1 - 275, x1 - 54, y1 - 56)
    rounded(draw, insight, 26, "#211E17", "#69532B", 2)
    text(draw, (insight[0] + 32, insight[1] + 34), "Reader intuition", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (insight[0] + 32, insight[1] + 88),
        "AUROC ignores the absolute probability scale and focuses on rank separation between classes.",
        F["small"],
        insight[2] - insight[0] - 64,
        COLORS["text"],
        8,
    )


def draw_result_axis(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], examples: list[dict[str, object]]) -> None:
    x0, y0, x1, y1 = box
    axis_x0 = x0 + 315
    axis_x1 = x1 - 385
    row_y0 = y0 + 325
    row_gap = 205
    min_v, max_v = 0.45, 1.00

    def sx(value: float) -> int:
        return int(axis_x0 + (value - min_v) / (max_v - min_v) * (axis_x1 - axis_x0))

    # Shared x-axis and reference gates.
    for v in [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
        x = sx(v)
        draw.line((x, row_y0 - 100, x, row_y0 + row_gap * (len(examples) - 1) + 112), fill="#223044", width=2)
        text(draw, (x, row_y0 + row_gap * len(examples) - 30), f"{v:.1f}", F["tiny"], COLORS["axis"], anchor="ma")
    dashed_vertical(draw, sx(0.50), row_y0 - 118, row_y0 + row_gap * (len(examples) - 1) + 115, COLORS["amber"], 4, 18)
    draw.line((sx(0.60), row_y0 - 118, sx(0.60), row_y0 + row_gap * (len(examples) - 1) + 115), fill=rgba(COLORS["blue"], 170), width=4)
    draw.line((sx(0.70), row_y0 - 118, sx(0.70), row_y0 + row_gap * (len(examples) - 1) + 115), fill=rgba(COLORS["teal"], 190), width=5)
    text(draw, (sx(0.50), row_y0 - 160), "chance", F["tiny_bold"], COLORS["amber"], anchor="ma")
    text(draw, (sx(0.60), row_y0 - 160), "CI gate", F["tiny_bold"], COLORS["blue"], anchor="ma")
    text(draw, (sx(0.70), row_y0 - 160), "score gate", F["tiny_bold"], COLORS["teal"], anchor="ma")

    legend_y = y0 + 187
    rounded(draw, (x0 + 48, legend_y, x1 - 48, legend_y + 82), 20, "#151F2D", "#2A394D", 2)
    text(draw, (x0 + 78, legend_y + 26), "dot", F["small_bold"], COLORS["teal"])
    text(draw, (x0 + 130, legend_y + 26), "= AUROC", F["small"], COLORS["muted"])
    dashed_horizontal(draw, x0 + 300, legend_y + 42, x0 + 420, COLORS["axis"], 5, 16)
    text(draw, (x0 + 445, legend_y + 26), "= bootstrap 95% CI", F["small"], COLORS["muted"])
    text(draw, (x1 - 78, legend_y + 26), "higher = stronger separation", F["small_bold"], COLORS["text"], anchor="ra")

    for i, rec in enumerate(examples):
        y = row_y0 + i * row_gap
        color = str(rec["color"])
        ci0 = sx(float(rec["ci_lower"]))
        ci1 = sx(float(rec["ci_upper"]))
        dot = sx(float(rec["auroc"]))

        rounded(draw, (x0 + 48, y - 70, x1 - 48, y + 90), 20, "#111A26", "#26364A", 1)
        text(draw, (x0 + 78, y - 42), str(rec["label"]), F["small_bold"], COLORS["text"])
        text(draw, (x0 + 78, y - 5), f"p={fmt_p(float(rec['p']))}", F["tiny_bold"], COLORS["muted"])
        rounded(draw, (x0 + 78, y + 28, x0 + 245, y + 68), 14, "#151F2D", color, 2)
        text(draw, (x0 + 162, y + 38), str(rec["tag"]), F["micro"], color, anchor="ma")

        draw.line((ci0, y, ci1, y), fill=COLORS["axis"], width=10)
        draw.line((ci0, y - 20, ci0, y + 20), fill=COLORS["axis"], width=5)
        draw.line((ci1, y - 20, ci1, y + 20), fill=COLORS["axis"], width=5)
        draw.ellipse((dot - 25, y - 25, dot + 25, y + 25), fill=color, outline=COLORS["white"], width=3)
        value = f"{float(rec['auroc']):.3f} [{float(rec['ci_lower']):.3f}, {float(rec['ci_upper']):.3f}]"
        value_x1 = x1 - 70
        value_w = int(draw.textlength(value, font=F["tiny_bold"])) + 24
        rounded(draw, (value_x1 - value_w, y - 34, value_x1, y + 6), 12, "#0F1722", "#2A394D", 1)
        text(draw, (value_x1 - 12, y - 28), value, F["tiny_bold"], COLORS["text"], anchor="ra")


def draw_reading_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], examples: list[dict[str, object]]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "B. Result row anatomy", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "Read the score, the interval, and the reference gates together before moving to biology.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )
    draw_result_axis(draw, box, examples)

    callout = (x0 + 48, y1 - 212, x1 - 48, y1 - 56)
    rounded(draw, callout, 24, "#151F2D", "#2A394D", 2)
    text(draw, (callout[0] + 32, callout[1] + 30), "Readout", F["h3"], COLORS["sky"])
    paragraph(
        draw,
        (callout[0] + 32, callout[1] + 78),
        "Rows above chance with tight intervals make the task result ready for biological interpretation.",
        F["small"],
        callout[2] - callout[0] - 64,
        COLORS["text"],
        7,
    )


def draw_gate(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    number: str,
    title: str,
    body: str,
    color: str,
) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 24, "#151F2D", "#2A394D", 2)
    rounded(draw, (x0 + 28, y0 + 28, x0 + 88, y0 + 88), 18, "#101823", color, 2)
    text(draw, (x0 + 58, y0 + 43), number, F["small_bold"], color, anchor="ma")
    text(draw, (x0 + 116, y0 + 27), title, F["gate"], COLORS["text"])
    paragraph(draw, (x0 + 116, y0 + 83), body, F["small"], x1 - x0 - 150, COLORS["muted"], 6)


def draw_decision_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "C. Benchmark score gate", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "A GO result passes score, uncertainty, and permutation checks together.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    gates = [
        ("1", "Mean AUROC > 0.700", "Average separation across folds clears the benchmark score gate.", COLORS["teal"]),
        ("2", "CI lower > 0.600", "Bootstrap uncertainty keeps the lower bound above the stability gate.", COLORS["blue"]),
        ("3", "perm p < 0.050", "Label permutation places the result beyond the chance reference.", COLORS["green"]),
    ]
    gate_y = y0 + 228
    for idx, gate in enumerate(gates):
        draw_gate(draw, (x0 + 54, gate_y, x1 - 54, gate_y + 198), *gate)
        gate_y += 238
        if idx < len(gates) - 1:
            text(draw, ((x0 + x1) // 2, gate_y - 30), "AND", F["small_bold"], COLORS["amber"], anchor="ma")

    final_box = (x0 + 54, y0 + 930, x1 - 54, y0 + 1088)
    rounded(draw, final_box, 28, "#18231E", COLORS["green"], 3)
    text(draw, (final_box[0] + 46, final_box[1] + 42), "GO", F["score"], COLORS["green"])
    paragraph(
        draw,
        (final_box[0] + 205, final_box[1] + 45),
        "The task result moves forward as reproducible held-out separation.",
        F["small_bold"],
        final_box[2] - final_box[0] - 250,
        COLORS["text"],
        8,
    )

    scope_box = (x0 + 54, y0 + 1134, x1 - 54, y1 - 56)
    rounded(draw, scope_box, 28, "#211E17", "#69532B", 2)
    text(draw, (scope_box[0] + 34, scope_box[1] + 36), "Readout frame", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (scope_box[0] + 34, scope_box[1] + 90),
        "AUROC is the held-out discrimination score. Companion slides connect this readout to mechanism and external biology.",
        F["small"],
        scope_box[2] - scope_box[0] - 68,
        COLORS["text"],
        8,
    )


def main() -> None:
    examples = load_examples()
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)

    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 48), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 42), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "METRIC PRIMER", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "Read AUROC as separation plus uncertainty", F["title"])
    text(
        draw,
        (150, 216),
        "Each score is interpreted against chance, interval width, and the benchmark score gate.",
        F["subtitle"],
        COLORS["muted"],
    )
    draw_header_badges(draw)

    draw_auroc_panel(draw, (150, 350, 1210, 1800))
    draw_reading_panel(draw, (1250, 350, 2605, 1800), examples)
    draw_decision_panel(draw, (2645, 350, 3690, 1800))

    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    text(draw, (205, 1925), "Takeaway", F["small_bold"], COLORS["sky"])
    source = "AUROC is useful only when the audience sees chance, uncertainty, and permutation context together."
    paragraph(draw, (390, 1925), source, F["small"], 3140, COLORS["muted"], 7)
    text(draw, (205, 1995), "Next", F["small_bold"], COLORS["amber"])
    scope = "The next method slides use this score grammar to compare classical ML, foundation models, and text LLM checks."
    paragraph(draw, (390, 1995), scope, F["small"], 3140, COLORS["muted"], 7)

    png = OUT_DIR / "metric_primer_auroc_uncertainty_premium.png"
    canvas.save(png, quality=95)
    gray = OUT_DIR / "metric_primer_auroc_uncertainty_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    manifest = {
        "title": "Read AUROC as separation plus uncertainty",
        "sources": [
            "README.md",
            "docs/submission_format.md",
            "evaluation/A4_holdout_results.json",
            "evaluation/A5_holdout_results.json",
            "evaluation/A3_baseline_results.json",
        ],
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "examples": [
            {
                "label": str(item["label"]),
                "auroc": round(float(item["auroc"]), 3),
                "ci_lower": round(float(item["ci_lower"]), 3),
                "ci_upper": round(float(item["ci_upper"]), 3),
                "perm_p": float(item["p"]),
                "tag": str(item["tag"]),
            }
            for item in examples
        ],
        "decision_rule": {
            "mean_auroc": "> 0.700",
            "ci_lower": "> 0.600",
            "perm_p": "< 0.050",
            "logic": "AND",
        },
        "scope": (
            "AUROC is the benchmark readout for held-out separation; CI and permutation p-value "
            "make the decision rule reproducible."
        ),
    }
    manifest_path = OUT_DIR / "metric_primer_auroc_uncertainty_premium_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

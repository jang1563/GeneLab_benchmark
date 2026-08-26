#!/usr/bin/env python3
"""Build the detailed-deck model-family bridge teaching asset."""

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
    / "model_family_bridge"
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
    "teal": "#5FD3C4",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "violet": "#8C63F7",
    "violet2": "#B69CFF",
    "magenta": "#C447D8",
    "rose": "#E17882",
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
    "stat": load_font(84, True),
    "lane": load_font(40, True),
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


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: str, width: int = 5) -> None:
    x0, y0 = start
    x1, y1 = end
    draw.line((x0, y0, x1, y1), fill=color, width=width)
    if x1 >= x0:
        draw.polygon([(x1, y1), (x1 - 22, y1 - 13), (x1 - 22, y1 + 13)], fill=color)
    else:
        draw.polygon([(x1, y1), (x1 + 22, y1 - 13), (x1 + 22, y1 + 13)], fill=color)


def dashed_vertical(draw: ImageDraw.ImageDraw, x: int, y0: int, y1: int, color: str, width: int = 3, dash: int = 18) -> None:
    y = y0
    while y < y1:
        draw.line((x, y, x, min(y + dash, y1)), fill=color, width=width)
        y += dash * 2


def draw_header_badges(draw: ImageDraw.ImageDraw) -> None:
    badges = [
        ("FIXED TASK", "split + view + metric", 405),
        ("FAMILY GUIDE", "classical + FM + LLM", 430),
        ("READOUT", "one score surface", 350),
    ]
    bx = 2300
    for kicker, body, badge_w in badges:
        rounded(draw, (bx, 72, bx + badge_w, 176), 26, "#121B28", "#2A394D", 2)
        text(draw, (bx + 28, 96), kicker, F["tiny"], COLORS["teal"])
        text(draw, (bx + 28, 126), body, F["small_bold"], COLORS["text"])
        bx += badge_w + 30


def draw_matrix(draw: ImageDraw.ImageDraw, x: int, y: int, cell: int = 24, rows: int = 6, cols: int = 9) -> None:
    palette = [COLORS["blue"], COLORS["teal"], COLORS["amber"], COLORS["violet"], "#3E5066"]
    for r in range(rows):
        for c in range(cols):
            fill = rgba(palette[(r * 2 + c * 3) % len(palette)], 210 if (r + c) % 3 else 130)
            draw.rounded_rectangle(
                (x + c * cell, y + r * cell, x + c * cell + cell - 5, y + r * cell + cell - 5),
                radius=5,
                fill=fill,
            )


def draw_token_stack(draw: ImageDraw.ImageDraw, x: int, y: int, color: str, width: int = 280) -> None:
    labels = ["Ifng", "Myc", "Dbp", "E2f1", "Apoe"]
    for i, label in enumerate(labels):
        yy = y + i * 43
        bar_w = width - i * 24
        rounded(draw, (x, yy, x + bar_w, yy + 27), 12, rgba(color, 215), None, 0)
        text(draw, (x + bar_w + 16, yy), label, F["tiny"], COLORS["muted"])


def draw_prompt_card(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, h: int) -> None:
    rounded(draw, (x, y, x + w, y + h), 22, "#0F1723", "#2A394D", 2)
    text(draw, (x + 26, y + 26), "Prompt summary", F["small_bold"], COLORS["text"])
    lines = [
        "Top variable genes:",
        "Ifng, Myc, Dbp...",
        "",
        "Predict Flight",
        "or Ground",
    ]
    multiline(draw, (x + 26, y + 68), lines, F["tiny"], COLORS["muted"], 6)


def draw_contract_step(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    number: str,
    title: str,
    body: str,
    color: str,
    icon: str,
) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 24, "#151F2D", "#2A394D", 2)
    rounded(draw, (x0 + 28, y0 + 28, x0 + 92, y0 + 92), 18, "#101823", color, 2)
    text(draw, (x0 + 60, y0 + 43), number, F["small_bold"], color, anchor="ma")
    text(draw, (x0 + 122, y0 + 28), title, F["h3"], COLORS["text"])
    paragraph(draw, (x0 + 122, y0 + 76), body, F["small"], x1 - x0 - 372, COLORS["muted"], 6)

    ix = x1 - 140
    iy = y0 + 60
    if icon == "split":
        rounded(draw, (ix - 90, iy - 38, ix - 20, iy + 38), 18, "#101823", COLORS["teal"], 2)
        rounded(draw, (ix + 20, iy - 38, ix + 90, iy + 38), 18, "#241719", COLORS["rose"], 2)
        text(draw, (ix - 55, iy - 12), "train", F["micro"], COLORS["teal"], anchor="ma")
        text(draw, (ix + 55, iy - 12), "test", F["micro"], COLORS["rose"], anchor="ma")
    elif icon == "view":
        draw_matrix(draw, ix - 92, iy - 56, cell=18, rows=5, cols=8)
    else:
        dashed_vertical(draw, ix - 30, iy - 55, iy + 55, COLORS["amber"], 4, 16)
        draw.ellipse((ix + 18, iy - 24, ix + 66, iy + 24), fill=color, outline=COLORS["white"], width=3)
        text(draw, (ix - 30, iy + 67), "0.5", F["micro"], COLORS["amber"], anchor="ma")


def draw_contract_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "A. Fix the task contract", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "Before model names matter, the benchmark fixes the split, input view, and score grammar.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    step_y = y0 + 210
    steps = [
        ("1", "Held-out split", "One hidden mission becomes the test unit.", COLORS["teal"], "split"),
        ("2", "Input view", "The model receives a defined feature representation.", COLORS["blue"], "view"),
        ("3", "Metric rule", "AUROC, CI, and permutation p produce one readout.", COLORS["amber"], "metric"),
    ]
    for number, title, body, color, icon in steps:
        draw_contract_step(draw, (x0 + 54, step_y, x1 - 54, step_y + 245), number, title, body, color, icon)
        step_y += 300

    rule_box = (x0 + 54, y1 - 270, x1 - 54, y1 - 56)
    rounded(draw, rule_box, 28, "#211E17", "#69532B", 2)
    text(draw, (rule_box[0] + 34, rule_box[1] + 36), "Comparison path", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (rule_box[0] + 34, rule_box[1] + 92),
        "Compare model families only after these three pieces stay fixed.",
        F["small_bold"],
        rule_box[2] - rule_box[0] - 68,
        COLORS["text"],
        8,
    )


def draw_lane(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    tier: str,
    title: str,
    definition: str,
    sees: str,
    examples: str,
    color: str,
    visual: str,
) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 28, "#151F2D", "#2A394D", 2)
    rounded(draw, (x0 + 32, y0 + 30, x0 + 218, y0 + 78), 17, "#101823", color, 2)
    text(draw, (x0 + 55, y0 + 43), tier, F["tiny_bold"], color)
    text(draw, (x0 + 32, y0 + 100), title, F["lane"], COLORS["text"])
    paragraph(draw, (x0 + 32, y0 + 155), definition, F["small"], 440, COLORS["muted"], 6)
    text(draw, (x0 + 32, y1 - 86), "Examples:", F["tiny_bold"], color)
    paragraph(draw, (x0 + 166, y1 - 88), examples, F["tiny"], 420, COLORS["text"], 4)

    vx = x0 + 560
    vy = y0 + 90
    if visual == "matrix":
        text(draw, (vx, vy - 34), "numeric table", F["small_bold"], COLORS["text"])
        draw_matrix(draw, vx, vy, cell=27, rows=6, cols=9)
        arrow(draw, (vx + 275, vy + 82), (vx + 390, vy + 82), color, 4)
        rounded(draw, (vx + 420, vy + 18, vx + 650, vy + 148), 22, "#0F1723", color, 2)
        text(draw, (vx + 454, vy + 55), "train", F["small_bold"], COLORS["text"])
        text(draw, (vx + 454, vy + 92), "on task", F["small"], COLORS["muted"])
    elif visual == "tokens":
        text(draw, (vx, vy - 34), "rank / embedding", F["small_bold"], COLORS["text"])
        draw_token_stack(draw, vx, vy, color, width=300)
        arrow(draw, (vx + 365, vy + 95), (vx + 482, vy + 95), color, 4)
        rounded(draw, (vx + 515, vy + 28, vx + 710, vy + 158), 22, "#0F1723", color, 2)
        text(draw, (vx + 548, vy + 64), "adapt", F["small_bold"], COLORS["text"])
        text(draw, (vx + 548, vy + 101), "to bulk", F["small"], COLORS["muted"])
    else:
        text(draw, (vx, vy - 34), "language prompt", F["small_bold"], COLORS["text"])
        draw_prompt_card(draw, vx, vy, 345, 210)
        arrow(draw, (vx + 375, vy + 105), (vx + 480, vy + 105), color, 4)
        rounded(draw, (vx + 515, vy + 45, vx + 705, vy + 165), 22, "#0F1723", color, 2)
        text(draw, (vx + 550, vy + 77), "A/B", F["h3"], COLORS["text"])
        text(draw, (vx + 548, vy + 116), "confidence", F["tiny"], COLORS["muted"])

    rounded(draw, (x1 - 330, y1 - 92, x1 - 44, y1 - 38), 18, "#101823", color, 2)
    text(draw, (x1 - 306, y1 - 77), sees, F["tiny_bold"], COLORS["text"])


def draw_family_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "B. Three input languages", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "The same Flight-vs-Ground question is translated into the input each model family can use.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    lanes = [
        (
            "TIER 1",
            "Classical ML",
            "Task-trained numeric learner.",
            "sees matrix",
            "PCA-LR / RF / LR",
            COLORS["blue"],
            "matrix",
        ),
        (
            "TIER 2",
            "Foundation model",
            "Pretrained expression model, then adapted.",
            "sees tokens",
            "scGPT / Geneformer / UCE",
            COLORS["violet"],
            "tokens",
        ),
        (
            "TIER 3",
            "Text LLM check",
            "Prompt-only language model.",
            "sees prompt",
            "DeepSeek / Gemini / Llama",
            COLORS["magenta"],
            "prompt",
        ),
    ]
    lane_y = y0 + 190
    for lane in lanes:
        draw_lane(draw, (x0 + 48, lane_y, x1 - 48, lane_y + 370), *lane)
        lane_y += 410


def draw_score_surface(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "C. Read one result surface", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "After the task contract is fixed, each family returns a comparable Flight probability.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        7,
    )

    lock_y = y0 + 210
    locks = [
        ("SPLIT", "held-out mission", COLORS["teal"]),
        ("VIEW", "feature protocol", COLORS["blue"]),
        ("METRIC", "AUROC rule", COLORS["amber"]),
    ]
    lx = x0 + 58
    for title, body, color in locks:
        rounded(draw, (lx, lock_y, lx + 292, lock_y + 130), 24, "#151F2D", "#2A394D", 2)
        text(draw, (lx + 28, lock_y + 28), title, F["tiny_bold"], color)
        text(draw, (lx + 28, lock_y + 64), body, F["small_bold"], COLORS["text"])
        lx += 322

    chart = (x0 + 68, y0 + 450, x1 - 68, y0 + 1045)
    rounded(draw, chart, 28, "#0F1723", "#2A394D", 2)
    text(draw, (chart[0] + 34, chart[1] + 34), "Example score grammar", F["h3"], COLORS["text"])
    axis_x0, axis_x1 = chart[0] + 260, chart[2] - 80
    row_y0, row_gap = chart[1] + 180, 122
    min_v, max_v = 0.40, 0.80

    def sx(value: float) -> int:
        return int(axis_x0 + (value - min_v) / (max_v - min_v) * (axis_x1 - axis_x0))

    for tick in [0.4, 0.5, 0.6, 0.7, 0.8]:
        tx = sx(tick)
        draw.line((tx, row_y0 - 88, tx, row_y0 + row_gap * 3 - 35), fill=rgba(COLORS["grid"], 150), width=2)
        text(draw, (tx, row_y0 + row_gap * 3 - 4), f"{tick:.1f}", F["tiny"], COLORS["axis"], anchor="ma")
    dashed_vertical(draw, sx(0.5), row_y0 - 92, row_y0 + row_gap * 3 - 35, COLORS["amber"], 4, 18)
    text(draw, (sx(0.5), row_y0 - 78), "chance", F["tiny_bold"], COLORS["amber"], anchor="ma")

    rows = [
        ("Classical ML", 0.758, COLORS["blue"], "shared 6-tissue mean"),
        ("Foundation model", 0.666, COLORS["violet"], "scGPT shared mean"),
        ("Text LLM check", (0.47, 0.51), COLORS["magenta"], "zero-shot range"),
    ]
    for i, (label, value, color, note) in enumerate(rows):
        y = row_y0 + i * row_gap
        text(draw, (chart[0] + 34, y - 25), label, F["small_bold"], COLORS["text"])
        text(draw, (chart[0] + 34, y + 8), note, F["micro"], COLORS["muted"])
        draw.line((axis_x0, y, axis_x1, y), fill="#263344", width=9)
        if isinstance(value, tuple):
            x_low, x_high = sx(value[0]), sx(value[1])
            rounded(draw, (x_low, y - 17, x_high, y + 17), 16, color, None, 0)
            draw.ellipse((x_low - 13, y - 13, x_low + 13, y + 13), fill=color, outline=COLORS["white"], width=2)
            draw.ellipse((x_high - 13, y - 13, x_high + 13, y + 13), fill=color, outline=COLORS["white"], width=2)
            text(draw, (axis_x1, y - 33), "0.47-0.51", F["tiny_bold"], COLORS["text"], anchor="ra")
        else:
            px = sx(value)
            draw.line((axis_x0, y, px, y), fill=color, width=13)
            draw.ellipse((px - 18, y - 18, px + 18, y + 18), fill=color, outline=COLORS["white"], width=3)
            text(draw, (axis_x1, y - 33), f"{value:.3f}", F["tiny_bold"], COLORS["text"], anchor="ra")

    takeaway = (x0 + 68, y0 + 1120, x1 - 68, y1 - 56)
    rounded(draw, takeaway, 28, "#211E17", "#69532B", 2)
    text(draw, (takeaway[0] + 34, takeaway[1] + 36), "Takeaway", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (takeaway[0] + 34, takeaway[1] + 92),
        "This is a controlled model-family stress test: same task contract, different input languages, one benchmark readout.",
        F["small_bold"],
        takeaway[2] - takeaway[0] - 68,
        COLORS["text"],
        9,
    )
    paragraph(
        draw,
        (takeaway[0] + 34, takeaway[1] + 185),
        "The next result slide can focus on performance because the audience already knows what each family received.",
        F["small"],
        takeaway[2] - takeaway[0] - 68,
        COLORS["muted"],
        7,
    )


def main() -> None:
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)

    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 48), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 42), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "MODEL-FAMILY BRIDGE", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "Model families share one task contract", F["title"])
    text(
        draw,
        (150, 216),
        "Classical ML, foundation models, and text LLM checks are compared after split, input view, and AUROC rule are fixed.",
        F["subtitle"],
        COLORS["muted"],
    )
    draw_header_badges(draw)

    draw_contract_panel(draw, (150, 350, 1030, 1800))
    draw_family_panel(draw, (1080, 350, 2585, 1800))
    draw_score_surface(draw, (2635, 350, 3690, 1800))

    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    text(draw, (205, 1925), "Takeaway", F["small_bold"], COLORS["blue"])
    source = "Model names are comparable only after split, input view, and score grammar are fixed."
    paragraph(draw, (390, 1925), source, F["small"], 3140, COLORS["muted"], 7)
    text(draw, (205, 1995), "Next", F["small_bold"], COLORS["amber"])
    scope = "The next slide shows the three input surfaces side by side before the result slides begin."
    paragraph(draw, (390, 1995), scope, F["small"], 3140, COLORS["muted"], 7)

    png = OUT_DIR / "model_family_bridge_premium.png"
    canvas.save(png, quality=95)
    gray = OUT_DIR / "model_family_bridge_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    manifest = {
        "title": "Model families share one task contract",
        "sources": [
            "README.md",
            "docs/text_llm_format.md",
            "evaluation/RESULTS_SUMMARY.md",
            "docs/SPACEBIOBENCH_CONFERENCE_DECK_METHOD_FOUNDATION_MODEL_FLOW_REVIEW_2026_06_09.md",
        ],
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "families": [
            {"tier": "Tier 1", "label": "Classical ML", "input": "numeric gene or pathway table"},
            {"tier": "Tier 2", "label": "Foundation model", "input": "ranked genes or embeddings adapted from expression pretraining"},
            {"tier": "Tier 3", "label": "Text LLM check", "input": "natural-language gene-list prompt"},
        ],
        "example_scores": {
            "pca_lr_shared_mean": 0.758,
            "scgpt_shared_mean": 0.666,
            "text_llm_zero_shot_range": [0.47, 0.51],
        },
        "scope": "This bridge defines model-family inputs; the companion model-family result slide reports the full comparison.",
    }
    manifest_path = OUT_DIR / "model_family_bridge_premium_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

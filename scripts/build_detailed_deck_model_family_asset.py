#!/usr/bin/env python3
"""Build the detailed-deck model-family methodology proof asset.

This asset is a teaching slide for mixed ML / biology audiences. It explains
that classical ML, gene-expression foundation models, and text LLM checks enter
the benchmark through different input representations on related spaceflight
detection tasks.
"""

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
    / "model_family"
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
    "magenta": "#C447D8",
    "violet": "#8C63F7",
    "violet2": "#B69CFF",
    "rose": "#E17882",
    "white": "#FFFFFF",
}


MODEL_ROWS = [
    {
        "name": "PCA-LR baseline",
        "family": "Classical ML",
        "score": 0.758,
        "kind": "point",
        "color": COLORS["blue"],
        "note": "shared 6-tissue mean",
    },
    {
        "name": "scGPT",
        "family": "Foundation model",
        "score": 0.666,
        "kind": "point",
        "color": COLORS["violet"],
        "note": "shared 6-tissue mean",
    },
    {
        "name": "Mouse-Geneformer",
        "family": "Foundation model",
        "score": 0.476,
        "kind": "point",
        "color": COLORS["violet2"],
        "note": "shared 6-tissue mean",
    },
    {
        "name": "Text LLMs",
        "family": "Prompt diagnostic",
        "score": (0.47, 0.51),
        "kind": "range",
        "color": COLORS["magenta"],
        "note": "DeepSeek/Gemini/Llama range",
    },
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
    "title": load_font(82, True),
    "subtitle": load_font(36, False),
    "h2": load_font(44, True),
    "h3": load_font(34, True),
    "body": load_font(30, False),
    "body_bold": load_font(30, True),
    "small": load_font(27, False),
    "small_bold": load_font(27, True),
    "tiny": load_font(22, False),
    "stat": load_font(82, True),
    "lane": load_font(38, True),
}


def rounded(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], radius: int, fill: str, outline: str | None = None, width: int = 1) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], value: str, font: ImageFont.ImageFont, fill: str = COLORS["text"], anchor: str | None = None) -> None:
    draw.text(xy, value, font=font, fill=fill, anchor=anchor)


def wrap_lines(draw: ImageDraw.ImageDraw, label: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
    words = label.split()
    lines: list[str] = []
    cur: list[str] = []
    for word in words:
        trial = " ".join(cur + [word])
        if draw.textlength(trial, font=font) <= max_width:
            cur.append(word)
        else:
            if cur:
                lines.append(" ".join(cur))
            cur = [word]
    if cur:
        lines.append(" ".join(cur))
    return lines


def multiline(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    lines: Iterable[str],
    font: ImageFont.ImageFont,
    fill: str = COLORS["muted"],
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
    fill: str = COLORS["muted"],
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


def draw_badges(draw: ImageDraw.ImageDraw) -> None:
    badges = [
        ("PRIMER", "model inputs", 310),
        ("INPUT", "shared rows", 300),
        ("INPUT", "different modalities", 380),
    ]
    bx = 2380
    for kicker, body, badge_w in badges:
        rounded(draw, (bx, 72, bx + badge_w, 176), 26, "#121B28", "#2A394D", 2)
        text(draw, (bx + 28, 96), kicker, F["tiny"], COLORS["teal"])
        text(draw, (bx + 28, 126), body, F["small_bold"], COLORS["text"])
        bx += badge_w + 30


def draw_matrix(draw: ImageDraw.ImageDraw, x: int, y: int, cell: int = 28, rows: int = 7, cols: int = 10) -> None:
    colors = [COLORS["blue"], COLORS["teal"], COLORS["amber"], COLORS["violet"], "#3E5066"]
    for r in range(rows):
        for c in range(cols):
            idx = (r * 3 + c * 2) % len(colors)
            fill = rgba(colors[idx], 210 if (r + c) % 3 else 125)
            draw.rounded_rectangle(
                (x + c * cell, y + r * cell, x + c * cell + cell - 5, y + r * cell + cell - 5),
                radius=5,
                fill=fill,
            )
    text(draw, (x, y - 38), "sample x gene matrix", F["small_bold"], COLORS["text"])
    text(draw, (x, y + rows * cell + 8), "numeric expression values", F["tiny"], COLORS["muted"])


def draw_token_stack(draw: ImageDraw.ImageDraw, x: int, y: int, width: int = 330) -> None:
    genes = ["Ifng", "Myc", "Dbp", "E2f1", "Apoe", "Npas2"]
    for i, gene in enumerate(genes):
        yy = y + i * 38
        bar_w = width - i * 28
        rounded(draw, (x, yy, x + bar_w, yy + 24), 11, rgba(COLORS["violet"], 210), None, 0)
        text(draw, (x + bar_w + 15, yy - 2), gene, F["tiny"], COLORS["muted"])
    text(draw, (x, y - 38), "ranked genes / embeddings", F["small_bold"], COLORS["text"])


def draw_prompt(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    rounded(draw, (x, y, x + 430, y + 245), 24, "#0F1723", "#2A394D", 2)
    text(draw, (x + 28, y + 28), "Prompt summary", F["small_bold"], COLORS["text"])
    lines = [
        "Top variable genes:",
        "Ifng, Myc, Dbp, Apoe...",
        "",
        "Predict:",
        "Flight or Ground",
    ]
    multiline(draw, (x + 28, y + 72), lines, F["tiny"], COLORS["muted"], 7)


def draw_task_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 45), "A. One benchmark task", F["h2"], COLORS["text"])
    paragraph(draw, (x0 + 48, y0 + 98), "A task fixes the samples, labels, hidden mission split, feature view, and score before a model is read.", F["small"], x1 - x0 - 96, COLORS["muted"], 7)

    draw_matrix(draw, x0 + 78, y0 + 235, cell=34, rows=7, cols=10)

    train_box = (x0 + 64, y0 + 595, x1 - 64, y0 + 748)
    test_box = (x0 + 64, y0 + 790, x1 - 64, y0 + 943)
    rounded(draw, train_box, 24, "#151F2D", "#2A394D", 2)
    rounded(draw, test_box, 24, "#211E17", "#69532B", 2)
    text(draw, (x0 + 94, y0 + 622), "Training side", F["h3"], COLORS["teal"])
    paragraph(draw, (x0 + 94, y0 + 670), "Seen missions define scaling, features, and model fit.", F["small"], x1 - x0 - 185, COLORS["text"], 6)
    text(draw, (x0 + 94, y0 + 817), "Hidden mission", F["h3"], COLORS["amber"])
    paragraph(draw, (x0 + 94, y0 + 865), "The test mission stays unseen until scoring.", F["small"], x1 - x0 - 185, COLORS["text"], 6)

    rounded(draw, (x0 + 64, y0 + 1040, x1 - 64, y0 + 1238), 24, "#151F2D", "#2A394D", 2)
    text(draw, (x0 + 94, y0 + 1068), "Metric grammar", F["h3"], COLORS["text"])
    text(draw, (x0 + 94, y0 + 1120), "AUROC", F["stat"], COLORS["white"])
    text(draw, (x0 + 430, y0 + 1138), "0.5 = chance", F["small_bold"], COLORS["amber"])
    text(draw, (x0 + 430, y0 + 1176), "higher = better separation", F["small"], COLORS["muted"])
    rounded(draw, (x0 + 64, y1 - 128, x1 - 64, y1 - 48), 20, "#101823", "#2A394D", 2)
    text(draw, (x0 + 94, y1 - 112), "The task is the contract.", F["small_bold"], COLORS["teal"])
    text(draw, (x0 + 94, y1 - 82), "Models only compete cleanly inside it.", F["small"], COLORS["muted"])


def draw_lane_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    kicker: str,
    title: str,
    input_text: str,
    question: str,
    color: str,
    visual: str,
    result: str,
) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 30, COLORS["panel2"], "#2A394D", 2)
    rounded(draw, (x0 + 34, y0 + 32, x0 + 280, y0 + 82), 18, "#101823", color, 2)
    text(draw, (x0 + 58, y0 + 45), kicker, F["tiny"], color)
    text(draw, (x0 + 34, y0 + 104), title, F["lane"], COLORS["text"])

    text(draw, (x0 + 34, y0 + 154), "Input view", F["small_bold"], COLORS["teal"])
    paragraph(draw, (x0 + 34, y0 + 190), input_text, F["small"], 510, COLORS["text"], 6)
    text(draw, (x0 + 34, y0 + 255), "Question", F["small_bold"], COLORS["amber"])
    paragraph(draw, (x0 + 34, y0 + 286), question, F["small"], 560, COLORS["muted"], 4)

    visual_x = x0 + 690
    if visual == "matrix":
        draw_matrix(draw, visual_x, y0 + 130, cell=24, rows=6, cols=9)
        arrow(draw, (visual_x + 275, y0 + 220), (visual_x + 405, y0 + 220), color, 4)
        rounded(draw, (visual_x + 430, y0 + 155, visual_x + 700, y0 + 285), 22, "#0F1723", color, 2)
        text(draw, (visual_x + 462, y0 + 190), "fit model", F["small_bold"], COLORS["text"])
        text(draw, (visual_x + 462, y0 + 230), "score mission", F["small"], COLORS["muted"])
    elif visual == "tokens":
        draw_token_stack(draw, visual_x, y0 + 104, width=315)
        arrow(draw, (visual_x + 355, y0 + 205), (visual_x + 470, y0 + 205), color, 4)
        rounded(draw, (visual_x + 500, y0 + 132, visual_x + 700, y0 + 262), 22, "#0F1723", color, 2)
        text(draw, (visual_x + 530, y0 + 172), "adapt", F["small_bold"], COLORS["text"])
        text(draw, (visual_x + 530, y0 + 212), "to bulk task", F["small"], COLORS["muted"])
    else:
        draw_prompt(draw, visual_x, y0 + 112)
        arrow(draw, (visual_x + 455, y0 + 230), (visual_x + 550, y0 + 230), color, 4)
        rounded(draw, (visual_x + 580, y0 + 165, visual_x + 725, y0 + 292), 22, "#0F1723", color, 2)
        text(draw, (visual_x + 615, y0 + 198), "A/B", F["h3"], COLORS["text"])
        text(draw, (visual_x + 610, y0 + 238), "label", F["small"], COLORS["muted"])

    rounded(draw, (x1 - 355, y0 + 34, x1 - 40, y0 + 94), 18, "#101823", color, 2)
    text(draw, (x1 - 330, y0 + 50), result, F["small_bold"], COLORS["text"])


def draw_family_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 45, y0 + 45), "B. What each model family actually sees", F["h2"], COLORS["text"])
    paragraph(draw, (x0 + 45, y0 + 98), "The same biological question becomes three different input representations.", F["small"], x1 - x0 - 90, COLORS["muted"], 6)

    draw_lane_card(
        draw,
        (x0 + 38, y0 + 162, x1 - 38, y0 + 510),
        "TIER 1",
        "Classical ML",
        "A numeric sample-by-gene or pathway table.",
        "Can expression separate flight from ground?",
        COLORS["blue"],
        "matrix",
        "mean 0.758",
    )
    draw_lane_card(
        draw,
        (x0 + 38, y0 + 555, x1 - 38, y0 + 903),
        "TIER 2",
        "Gene-expression foundation model",
        "Ranked genes or embeddings adapted from pretrained single-cell models.",
        "Does pretraining transfer to small-n bulk data?",
        COLORS["violet"],
        "tokens",
        "best shared 0.666",
    )
    draw_lane_card(
        draw,
        (x0 + 38, y0 + 948, x1 - 38, y0 + 1296),
        "TIER 3",
        "Text LLM check",
        "A natural-language prompt summary of top genes.",
        "Can gene-list priors infer the label?",
        COLORS["magenta"],
        "prompt",
        "range 0.47-0.51",
    )

    rounded(draw, (x0 + 38, y1 - 118, x1 - 38, y1 - 42), 20, "#101823", "#2A394D", 2)
    text(draw, (x0 + 66, y1 - 92), "Comparison path", F["small_bold"], COLORS["amber"])
    text(draw, (x0 + 235, y1 - 92), "Compare results only after checking input modality and task rows.", F["small"], COLORS["muted"])


def draw_score_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 45, y0 + 45), "C. Reading the score surface", F["h2"], COLORS["text"])
    paragraph(draw, (x0 + 45, y0 + 98), "Shared rows are direct comparisons. Prompt-only and single-tissue rows are contextual.", F["small"], x1 - x0 - 90, COLORS["muted"], 6)

    axis_x0, axis_x1 = x0 + 360, x1 - 105
    row_y0 = y0 + 255
    row_gap = 138
    axis_min, axis_max = 0.40, 0.80

    def sx(value: float) -> int:
        return int(axis_x0 + (value - axis_min) / (axis_max - axis_min) * (axis_x1 - axis_x0))

    for tick in [0.40, 0.50, 0.60, 0.70, 0.80]:
        tx = sx(tick)
        draw.line((tx, row_y0 - 60, tx, row_y0 + row_gap * 4 - 35), fill=rgba(COLORS["grid"], 140), width=2)
        text(draw, (tx, row_y0 + row_gap * 4 - 8), f"{tick:.2f}", F["tiny"], COLORS["axis"], anchor="ma")
    chance_x = sx(0.50)
    draw.line((chance_x, row_y0 - 70, chance_x, row_y0 + row_gap * 4 - 28), fill=COLORS["amber"], width=3)
    text(draw, (chance_x, row_y0 - 104), "chance", F["tiny"], COLORS["amber"], anchor="ma")

    for i, row in enumerate(MODEL_ROWS):
        y = row_y0 + i * row_gap
        color = row["color"]
        text(draw, (x0 + 45, y - 28), row["name"], F["small_bold"], COLORS["text"])
        text(draw, (x0 + 45, y + 6), row["family"], F["tiny"], COLORS["muted"])
        draw.line((axis_x0, y, axis_x1, y), fill="#263344", width=10)
        if row["kind"] == "point":
            score = float(row["score"])
            px = sx(score)
            draw.line((axis_x0, y, px, y), fill=color, width=14)
            draw.ellipse((px - 17, y - 17, px + 17, y + 17), fill=color, outline=COLORS["white"], width=3)
            text(draw, (x1 - 42, y - 18), f"{score:.3f}", F["small_bold"], COLORS["text"], anchor="ra")
        else:
            low, high = row["score"]
            lx, hx = sx(low), sx(high)
            rounded(draw, (lx, y - 18, hx, y + 18), 16, color, None, 0)
            draw.ellipse((lx - 13, y - 13, lx + 13, y + 13), fill=color, outline=COLORS["white"], width=2)
            draw.ellipse((hx - 13, y - 13, hx + 13, y + 13), fill=color, outline=COLORS["white"], width=2)
            text(draw, (x1 - 42, y - 18), "0.47-0.51", F["small_bold"], COLORS["text"], anchor="ra")
        text(draw, (axis_x0, y + 32), row["note"], F["tiny"], COLORS["muted"])

    y = y0 + 900
    rounded(draw, (x0 + 45, y, x1 - 45, y + 202), 28, "#151F2D", "#2A394D", 2)
    text(draw, (x0 + 80, y + 32), "Context rows", F["h3"], COLORS["violet2"])
    paragraph(draw, (x0 + 80, y + 82), "UCE 0.632 and scFoundation 0.635 are best single-tissue extension rows; useful context beside the shared 6-task surface.", F["small"], x1 - x0 - 160, COLORS["text"], 7)

    y2 = y + 255
    rounded(draw, (x0 + 45, y2, x1 - 45, y2 + 300), 28, "#211E17", "#69532B", 2)
    text(draw, (x0 + 80, y2 + 36), "Readout frame", F["h3"], COLORS["amber"])
    lines = [
        "Result: classical baselines remain strong on this surface.",
        "Interpretation: pretraining needs task-specific adaptation.",
        "General FM capability remains outside this comparison.",
        "Text LLMs saw prompt summaries rather than matrices.",
    ]
    multiline(draw, (x0 + 80, y2 + 90), lines, F["small"], COLORS["text"], 10)


def main() -> None:
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)

    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 48), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 42), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "MODEL-FAMILY PRIMER", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "Same benchmark, three different model inputs", F["title"])
    text(
        draw,
        (150, 216),
        "Classical ML, gene-expression foundation models, and text LLMs answer related questions through different input views.",
        F["subtitle"],
        COLORS["muted"],
    )
    draw_badges(draw)

    draw_task_panel(draw, (150, 350, 895, 1800))
    draw_family_panel(draw, (930, 350, 2550, 1800))
    draw_score_panel(draw, (2585, 350, 3690, 1800))

    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    text(draw, (205, 1925), "Takeaway", F["small_bold"], COLORS["blue"])
    source = "The three model families answer related questions through different input views, so the input surface stays visible."
    paragraph(draw, (390, 1925), source, F["small"], 3140, COLORS["muted"], 7)
    text(draw, (205, 1995), "Next", F["small_bold"], COLORS["amber"])
    scope = "The following model result slides keep shared rows and context rows visually separated."
    paragraph(draw, (390, 1995), scope, F["small"], 3140, COLORS["muted"], 7)

    png = OUT_DIR / "model_family_input_surface_premium.png"
    canvas.save(png, quality=95)
    gray = OUT_DIR / "model_family_input_surface_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    manifest = {
        "title": "Same benchmark, three different model inputs",
        "source_files": [
            "README.md",
            "evaluation/RESULTS_SUMMARY.md",
            "docs/text_llm_format.md",
        ],
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "scores": {
            "pca_lr_shared_mean": 0.758,
            "scgpt_shared_mean": 0.666,
            "mouse_geneformer_shared_mean": 0.476,
            "text_llm_range": [0.47, 0.51],
            "uce_best_single_tissue": 0.632,
            "scfoundation_best_single_tissue": 0.635,
        },
        "scope": (
            "Classical ML, gene-expression foundation models, and text LLMs use different input representations. "
            "The current benchmark shows that pretrained and prompt-based approaches need task-specific adaptation to exceed tuned classical baselines under small-n bulk mission shift."
        ),
    }
    manifest_path = OUT_DIR / "model_family_input_surface_premium_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

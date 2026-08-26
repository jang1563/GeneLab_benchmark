#!/usr/bin/env python3
"""Build the detailed-deck foundation-model concept proof asset."""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
GENEFORMER_SUMMARY = ROOT / "evaluation" / "geneformer_mouse_gf_all_tissues_summary.json"
SCGPT_SUMMARY = ROOT / "evaluation" / "scgpt_whole_human_all_tissues_summary.json"
V1_CONTENT = ROOT / "docs" / "V1_PAPER_CONTENT.md"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "foundation_model_adds"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "what_a_foundation_model_adds_premium.png"
GRAY_PATH = OUT_DIR / "what_a_foundation_model_adds_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "foundation_model_adds_manifest.json"

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
    "stat": load_font(72, True),
    "huge": load_font(96, True),
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


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: str, width: int = 4) -> None:
    x0, y0 = start
    x1, y1 = end
    draw.line((x0, y0, x1, y1), fill=rgba(color, 180), width=width)
    if abs(x1 - x0) >= abs(y1 - y0):
        direction = 1 if x1 >= x0 else -1
        points = [(x1, y1), (x1 - direction * 18, y1 - 11), (x1 - direction * 18, y1 + 11)]
    else:
        direction = 1 if y1 >= y0 else -1
        points = [(x1, y1), (x1 - 11, y1 - direction * 18), (x1 + 11, y1 - direction * 18)]
    draw.polygon(points, fill=rgba(color, 210))


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


def load_data() -> dict[str, object]:
    geneformer = json.loads(GENEFORMER_SUMMARY.read_text())
    scgpt = json.loads(SCGPT_SUMMARY.read_text())
    gf_folds = sum(int(row["geneformer_n_folds"]) for row in geneformer["tissues"].values())
    return {
        "geneformer_pretrain_cells_m": 30,
        "scgpt_pretrain_cells_m": 33,
        "geneformer_folds": gf_folds,
        "scgpt_folds": int(scgpt["n_folds_total"]),
        "n_tasks": int(scgpt["n_tasks"]),
        "scgpt_mean": float(scgpt["overall_mean_auroc"]),
        "geneformer_mean": float(geneformer["overall"]["geneformer_mean"]),
        "baseline_mean": float(geneformer["overall"]["baseline_mean"]),
    }


def stat_badge(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, label: str, value: str, accent: str) -> None:
    rounded(draw, (x, y, x + w, y + 104), 24, "#121B28", "#2A394D", 2)
    text(draw, (x + 26, y + 23), label, F["tiny_bold"], accent)
    text(draw, (x + 26, y + 57), value, F["small_bold"], COLORS["text"])


def draw_header(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    text(draw, (136, 76), "MODELS / FOUNDATION MODEL PRIMER", F["kicker"], COLORS["sky"])
    text(draw, (136, 122), "A foundation model adds a learned representation", F["title"], COLORS["text"])
    text(
        draw,
        (138, 222),
        "Pretraining changes the input language and encoder; the task split, labels, and AUROC gate stay fixed.",
        F["subtitle"],
        COLORS["muted"],
    )
    badges = [
        ("FM ROWS", "2 models", 235, COLORS["violet"]),
        ("PRETRAIN", "30M / 33M cells", 330, COLORS["teal"]),
        ("TASKS", f"{int(data['n_tasks'])} tissues", 235, COLORS["sky"]),
        ("FOLDS", f"{int(data['geneformer_folds'])}+{int(data['scgpt_folds'])}", 215, COLORS["amber"]),
        ("GATE", "AUROC", 198, COLORS["green"]),
    ]
    bx = 2380
    for label, value, width, accent in badges:
        stat_badge(draw, bx, 72, width, label, value, accent)
        bx += width + 18


def draw_contract_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 46), "A. What stays fixed", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "A foundation model enters the same benchmark contract as the classical surface.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        leading=6,
    )

    rows = [
        ("1", "Same samples", "Bulk RNA-seq samples remain the evaluation unit.", COLORS["sky"]),
        ("2", "Same split", "Train and test missions are defined before model fitting.", COLORS["teal"]),
        ("3", "Same labels", "Spaceflight-vs-control labels define the task.", COLORS["amber"]),
        ("4", "Same score", "AUROC and uncertainty carry the readout.", COLORS["green"]),
    ]
    y = y0 + 240
    for idx, title, body, accent in rows:
        rounded(draw, (x0 + 60, y, x1 - 60, y + 162), 28, "#121B28", "#2A394D", 2)
        rounded(draw, (x0 + 92, y + 41, x0 + 160, y + 109), 18, rgba(accent, 52), rgba(accent, 190), 2)
        text(draw, (x0 + 126, y + 64), idx, F["small_bold"], accent, "mm")
        text(draw, (x0 + 190, y + 38), title, F["h3"], COLORS["text"])
        paragraph(draw, (x0 + 190, y + 84), body, F["tiny"], x1 - x0 - 286, COLORS["muted"], leading=4)
        y += 202

    readout_y = y1 - 178
    rounded(draw, (x0 + 60, readout_y, x1 - 60, y1 - 54), 30, rgba("#66A6E8", 34), rgba("#66A6E8", 150), 2)
    text(draw, (x0 + 96, readout_y + 35), "Readout", F["small_bold"], COLORS["sky"])
    paragraph(
        draw,
        (x0 + 96, readout_y + 72),
        "The comparison changes the model language, then keeps the benchmark gate constant.",
        F["body_bold"],
        x1 - x0 - 192,
        COLORS["text"],
        leading=5,
    )


def draw_token_strip(draw: ImageDraw.ImageDraw, x: int, y: int, labels: list[str], color: str) -> None:
    cx = x
    for label in labels:
        w = int(draw.textlength(label, font=F["micro_bold"])) + 28
        rounded(draw, (cx, y, cx + w, y + 38), 16, rgba(color, 48), rgba(color, 150), 1)
        text(draw, (cx + w / 2, y + 20), label, F["micro_bold"], COLORS["text"], "mm")
        cx += w + 10


def draw_encoder_icon(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], accent: str) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 28, "#121B28", rgba(accent, 170), 2)
    cols = 5
    rows = 4
    dx = (x1 - x0 - 120) / (cols - 1)
    dy = (y1 - y0 - 105) / (rows - 1)
    points: list[tuple[int, int]] = []
    for r in range(rows):
        for c in range(cols):
            px = int(x0 + 60 + c * dx)
            py = int(y0 + 54 + r * dy)
            points.append((px, py))
    for i, (px, py) in enumerate(points):
        for j in [i + 1, i + cols]:
            if j < len(points):
                qx, qy = points[j]
                if abs(qx - px) < 260 and abs(qy - py) < 160:
                    draw.line((px, py, qx, qy), fill=rgba(accent, 55), width=2)
    for i, (px, py) in enumerate(points):
        fill = accent if i % 3 == 0 else COLORS["blue"]
        draw.ellipse((px - 10, py - 10, px + 10, py + 10), fill=rgba(fill, 210), outline=rgba("#FFFFFF", 80))


def draw_added_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 46), "B. What gets added", F["h2"], COLORS["text"])
    text(draw, (x0 + 48, y0 + 98), "Pretraining supplies an encoder before task-specific adaptation.", F["small"], COLORS["muted"])

    lane_w = (x1 - x0 - 150) // 2
    left_x = x0 + 56
    right_x = left_x + lane_w + 38
    lane_y = y0 + 176
    for lx, title, corpus, tokens, accent in [
        (left_x, "Geneformer route", "~30M mouse scRNA cells", ["rank", "ENSMUSG", "top 2048"], COLORS["violet"]),
        (right_x, "scGPT route", "33M human CellXGene cells", ["bins", "ortholog", "human ckpt"], COLORS["teal"]),
    ]:
        rounded(draw, (lx, lane_y, lx + lane_w, lane_y + 242), 28, COLORS["panel2"], "#2A394D", 2)
        text(draw, (lx + 30, lane_y + 34), title, F["h3"], COLORS["text"])
        text(draw, (lx + 30, lane_y + 78), corpus, F["small_bold"], accent)
        paragraph(draw, (lx + 30, lane_y + 116), "A large expression corpus teaches gene-context patterns before this task.", F["tiny"], lane_w - 60, COLORS["muted"], leading=4)
        draw_token_strip(draw, lx + 30, lane_y + 184, tokens, accent)

    arrow(draw, (left_x + lane_w // 2, lane_y + 265), (x0 + 560, y0 + 540), COLORS["violet"])
    arrow(draw, (right_x + lane_w // 2, lane_y + 265), (x0 + 745, y0 + 540), COLORS["teal"])

    encoder_box = (x0 + 310, y0 + 530, x1 - 310, y0 + 840)
    draw_encoder_icon(draw, encoder_box, COLORS["sky"])
    text(draw, ((encoder_box[0] + encoder_box[2]) / 2, encoder_box[1] + 36), "pretrained encoder", F["h3"], COLORS["text"], "mm")
    text(draw, ((encoder_box[0] + encoder_box[2]) / 2, encoder_box[3] - 38), "maps sample input into a reusable representation", F["tiny_bold"], COLORS["muted"], "mm")

    arrow(draw, ((encoder_box[0] + encoder_box[2]) // 2, encoder_box[3] + 18), ((encoder_box[0] + encoder_box[2]) // 2, y0 + 1018), COLORS["sky"])
    head_y = y0 + 1016
    rounded(draw, (x0 + 260, head_y, x1 - 260, head_y + 160), 30, COLORS["panel2"], COLORS["amber"], 2)
    text(draw, (x0 + 306, head_y + 38), "task adapter / classifier head", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (x0 + 306, head_y + 84),
        "The spaceflight labels train the final readout layer under the same split.",
        F["tiny"],
        x1 - x0 - 612,
        COLORS["muted"],
        leading=4,
    )

    arrow(draw, ((x0 + x1) // 2, head_y + 176), ((x0 + x1) // 2, y1 - 286), COLORS["amber"])
    rounded(draw, (x0 + 260, y1 - 286, x1 - 260, y1 - 166), 30, COLORS["panel2"], COLORS["green"], 2)
    text(draw, ((x0 + x1) / 2, y1 - 246), "same AUROC gate", F["h3"], COLORS["green"], "mm")
    text(draw, ((x0 + x1) / 2, y1 - 204), "score the held-out task with the shared metric grammar", F["tiny_bold"], COLORS["text"], "mm")

    rounded(draw, (x0 + 60, y1 - 124, x1 - 60, y1 - 52), 24, "#121B28", "#2A394D", 1)
    text(draw, (x0 + 92, y1 - 84), "Concept: pretrain -> translate input -> encode sample -> adapt readout -> score", F["small_bold"], COLORS["text"])


def draw_lever_card(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    title: str,
    body: str,
    accent: str,
) -> None:
    rounded(draw, (x, y, x + w, y + 162), 28, "#121B28", "#2A394D", 2)
    rounded(draw, (x + 28, y + 42, x + 96, y + 110), 18, rgba(accent, 52), rgba(accent, 190), 2)
    text(draw, (x + 126, y + 36), title, F["h3"], COLORS["text"])
    paragraph(draw, (x + 126, y + 82), body, F["tiny"], w - 166, COLORS["muted"], leading=4)


def draw_levers_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "C. What pretraining can change", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 50, y0 + 98),
        "Read foundation-model rows as a test of which learned prior transfers into this bulk expression task.",
        F["small"],
        x1 - x0 - 100,
        COLORS["muted"],
        leading=6,
    )
    cards = [
        ("Pretraining prior", "Scale, cell-type diversity, and species coverage shape the starting representation.", COLORS["teal"]),
        ("Input translation", "Rank tokens, expression bins, and ortholog maps define what the encoder receives.", COLORS["violet"]),
        ("Domain surface", "Single-cell pretraining meets small-n bulk profiles with mixed cell populations.", COLORS["amber"]),
        ("Adapter budget", "Frozen layers and a classifier head determine how much task signal can move the model.", COLORS["green"]),
    ]
    y = y0 + 248
    for title, body, accent in cards:
        draw_lever_card(draw, x0 + 66, y, x1 - x0 - 132, title, body, accent)
        y += 194

    bridge_y = y1 - 376
    rounded(draw, (x0 + 66, bridge_y, x1 - 66, y1 - 64), 32, COLORS["panel2"], "#2A394D", 2)
    text(draw, (x0 + 104, bridge_y + 42), "Takeaway", F["h3"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 104, bridge_y + 92),
        "A foundation-model gain is read as tested-setting evidence over the matched classical surface.",
        F["body_bold"],
        x1 - x0 - 208,
        COLORS["text"],
        leading=7,
    )
    rounded(draw, (x0 + 104, y1 - 176, x1 - 104, y1 - 96), 22, "#101823", "#2A394D", 1)
    text(draw, (x0 + 130, y1 - 150), "Next", F["micro_bold"], COLORS["axis"])
    text(draw, (x0 + 130, y1 - 121), "Slide 25 separates prompt-only text LLM checks from expression-matrix adapters.", F["tiny_bold"], COLORS["sky"])


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    y = 1872
    rounded(draw, (136, y, 3704, 2042), 32, "#101823", "#2A394D", 2)
    text(draw, (180, y + 38), "Slide 24 readout", F["small_bold"], COLORS["sky"])
    paragraph(
        draw,
        (180, y + 76),
        "Foundation models add a pretrained expression representation; the benchmark reads that representation through the same task readout.",
        F["body_bold"],
        2220,
        COLORS["text"],
        leading=6,
    )
    paragraph(
        draw,
        (2610, y + 46),
        "Result slides should compare deltas over the classical surface and localize where gains appear.",
        F["tiny_bold"],
        920,
        COLORS["muted"],
        leading=5,
    )
    text(
        draw,
        (140, 2102),
        "Takeaway: pretraining adds reusable representation, but the benchmark score still comes from held-out mission adaptation.",
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
    draw_contract_panel(draw, (136, 345, 1068, 1815))
    draw_added_panel(draw, (1108, 345, 2608, 1815), data)
    draw_levers_panel(draw, (2648, 345, 3704, 1815), data)
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

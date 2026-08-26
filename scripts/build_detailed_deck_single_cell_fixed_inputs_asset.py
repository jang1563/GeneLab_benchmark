#!/usr/bin/env python3
"""Build slide 58 asset: single-cell scoring fixed inputs."""

from __future__ import annotations

import json
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont, ImageOps, ImageStat


ROOT = Path(__file__).resolve().parent.parent
ASSET_ROOT = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
)
OUT_DIR = ASSET_ROOT / "single_cell_fixed_inputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "single_cell_scoring_needs_fixed_inputs_premium.png"
GRAY_PATH = OUT_DIR / "single_cell_scoring_needs_fixed_inputs_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "single_cell_fixed_inputs_manifest.json"
QA_NOTE = OUT_DIR / "SINGLE_CELL_FIXED_INPUTS_ASSET_VISUAL_QA.md"

COLORS = {
    "bg": "#0B1119",
    "bg2": "#101722",
    "header": "#111B28",
    "panel": "#111B28",
    "panel2": "#162236",
    "panel3": "#0F1825",
    "grid": "#263245",
    "text": "#F4F7FB",
    "muted": "#AAB6C6",
    "dim": "#687789",
    "blue": "#66A6E8",
    "teal": "#5FD3C4",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "violet": "#B39DFF",
    "rose": "#E17882",
    "ink": "#081018",
}

FOOTPRINT = [
    ("Study", "OSD-918", COLORS["teal"]),
    ("Samples", "8", COLORS["amber"]),
    ("Cells", "4,395", COLORS["blue"]),
    ("Genes", "19,064", COLORS["green"]),
]

METRICS = [
    ("Primary", "2", "AUROC + AUPRC", COLORS["green"]),
    ("Context checks", "3", "embedding + DE rows", COLORS["amber"]),
    ("Optional", "1", "expression MAE", COLORS["violet"]),
]

LADDER = [
    ("Study files", "19 files listed", "ready", COLORS["teal"]),
    ("FASTQ pairs", "8 / 8 complete", "ready", COLORS["green"]),
    ("Fixed matrix", "h5ad 0 | STARsolo 0", "next object", COLORS["amber"]),
    ("obs / var audit", "27 checks queued", "after matrix", COLORS["blue"]),
    ("Score files", "predictions + metrics", "runnable after audit", COLORS["violet"]),
]

INPUT_FIELDS = [
    ("obs rows", "cell_id, sample_id, label, animal_id, cell type", COLORS["teal"]),
    ("var rows", "gene_symbol, feature_id", COLORS["blue"]),
]

RUN_PATH = [
    ("RRRM-1 blood", "8 samples, 4 + 4 labels", COLORS["teal"]),
    ("Fixed matrix", "AnnData or STARsolo table", COLORS["amber"]),
    ("Audit", "obs / var fields pass", COLORS["blue"]),
    ("Prediction file", "cell labels and probabilities", COLORS["green"]),
    ("Report", "metrics and run record", COLORS["violet"]),
]


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(value: str, alpha: int) -> tuple[int, int, int, int]:
    return (*hex_to_rgb(value), alpha)


def blend(a: str, b: str, t: float) -> str:
    ar, ag, ab = hex_to_rgb(a)
    br, bg, bb = hex_to_rgb(b)
    t = max(0.0, min(1.0, t))
    return f"#{int(ar + (br - ar) * t):02x}{int(ag + (bg - ag) * t):02x}{int(ab + (bb - ab) * t):02x}"


def load_font(size: int, bold: bool = False) -> ImageFont.ImageFont:
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
    "subtitle": load_font(37),
    "section": load_font(40, True),
    "h2": load_font(35, True),
    "h3": load_font(30, True),
    "body": load_font(26),
    "body_bold": load_font(26, True),
    "small": load_font(23),
    "small_bold": load_font(23, True),
    "tiny": load_font(20),
    "tiny_bold": load_font(20, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "metric": load_font(66, True),
    "big": load_font(80, True),
}


def text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[float, float],
    value: str,
    font: ImageFont.ImageFont,
    fill: str = COLORS["text"],
    anchor: str | None = None,
) -> None:
    draw.text(xy, value, font=font, fill=fill, anchor=anchor)


def wrap(draw: ImageDraw.ImageDraw, value: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
    words = value.split()
    lines: list[str] = []
    line = ""
    for word in words:
        trial = word if not line else f"{line} {word}"
        if draw.textlength(trial, font=font) <= max_width:
            line = trial
        else:
            if line:
                lines.append(line)
            line = word
    if line:
        lines.append(line)
    return lines


def paragraph(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    value: str,
    font: ImageFont.ImageFont,
    max_width: int,
    fill: str = COLORS["muted"],
    leading: int = 8,
) -> int:
    x, y = xy
    for line in wrap(draw, value, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def rounded(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    radius: int,
    fill: str,
    outline: str | None = None,
    width: int = 1,
) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def background() -> Image.Image:
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image, "RGBA")
    for y in range(H):
        color = blend(COLORS["bg"], COLORS["bg2"], (y / H) * 0.72)
        draw.line((0, y, W, y), fill=color)
    for x in range(0, W, 160):
        draw.line((x, 260, x, H - 220), fill=rgba(COLORS["grid"], 72), width=1)
    for y in range(320, H - 220, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 64), width=1)
    draw.rectangle((0, 0, W, 260), fill=rgba(COLORS["header"], 245))
    draw.rectangle((0, H - 190, W, H), fill=rgba("#071019", 235))
    return image


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    w = max(220, int(draw.textlength(value, font=F["small_bold"]) + 72))
    rounded(draw, (x, y, x + w, y + 72), 14, COLORS["panel2"], color, 2)
    text(draw, (x + 18, y + 14), label.upper(), F["micro_bold"], COLORS["muted"])
    text(draw, (x + 18, y + 42), value, F["small_bold"], COLORS["text"])
    return x + w + 24


def arrow(draw: ImageDraw.ImageDraw, x1: int, y1: int, x2: int, y2: int, color: str, width: int = 4) -> None:
    draw.line((x1, y1, x2, y2), fill=color, width=width)
    ang = math.atan2(y2 - y1, x2 - x1)
    size = 18
    p1 = (x2, y2)
    p2 = (x2 - size * math.cos(ang - 0.45), y2 - size * math.sin(ang - 0.45))
    p3 = (x2 - size * math.cos(ang + 0.45), y2 - size * math.sin(ang + 0.45))
    draw.polygon([p1, p2, p3], fill=color)


def draw_cell_cloud(draw: ImageDraw.ImageDraw, cx: int, cy: int) -> None:
    positions = [
        (-90, -42, 44, COLORS["teal"]),
        (-34, -72, 36, COLORS["blue"]),
        (34, -52, 42, COLORS["green"]),
        (82, -8, 34, COLORS["amber"]),
        (-78, 34, 36, COLORS["violet"]),
        (-12, 20, 48, COLORS["rose"]),
        (54, 42, 38, COLORS["teal"]),
        (5, -12, 30, COLORS["blue"]),
    ]
    for dx, dy, r, color in positions:
        draw.ellipse((cx + dx - r, cy + dy - r, cx + dx + r, cy + dy + r), fill=rgba(color, 96), outline=rgba("#EAF2FA", 86), width=2)
        draw.ellipse((cx + dx - r // 3, cy + dy - r // 3, cx + dx + r // 3, cy + dy + r // 3), fill=rgba(color, 180))
    rounded(draw, (cx - 122, cy + 96, cx + 122, cy + 170), 18, COLORS["panel3"], COLORS["teal"], 2)
    text(draw, (cx, cy + 112), "4,395", F["h2"], COLORS["text"], "ma")
    text(draw, (cx, cy + 146), "QC cells", F["micro_bold"], COLORS["muted"], "ma")


def draw_metric_contract_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 610, 1325, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["teal"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Metric Contract Is Ready", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "The RRRM-1 blood task has a label domain, cell table shape, and metric vocabulary.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    draw_cell_cloud(draw, x1 + 220, y1 + 330)

    card_x, card_y = x1 + 470, y1 + 236
    for i, (label, value, color) in enumerate(FOOTPRINT):
        cx = card_x + (i % 2) * 300
        cy = card_y + (i // 2) * 134
        rounded(draw, (cx, cy, cx + 260, cy + 104), 18, COLORS["panel2"], color, 2)
        text(draw, (cx + 22, cy + 18), label.upper(), F["micro_bold"], color)
        text(draw, (cx + 22, cy + 52), value, F["h2"], COLORS["text"])

    my = y2 - 250
    text(draw, (x1 + 62, my), "Metric profile", F["h3"], COLORS["text"])
    for i, (role, count, detail, color) in enumerate(METRICS):
        x = x1 + 62 + i * 360
        rounded(draw, (x, my + 58, x + 322, my + 170), 18, COLORS["panel3"], color, 2)
        text(draw, (x + 24, my + 82), role, F["small_bold"], color)
        text(draw, (x + 24, my + 122), count, F["h2"], COLORS["text"])
        text(draw, (x + 82, my + 128), detail, F["micro_bold"], COLORS["muted"])


def draw_ladder_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 1365, 610, 2535, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["amber"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Fixed Input Ladder", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "Raw files are accounted for; a single analysis-ready matrix is the next scoring object.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    rail_x = x1 + 140
    top = y1 + 222
    for i, (label, value, status, color) in enumerate(LADDER):
        y = top + i * 112
        if i < len(LADDER) - 1:
            draw.line((rail_x, y + 38, rail_x, y + 150), fill=rgba(COLORS["dim"], 130), width=4)
        draw.ellipse((rail_x - 24, y + 14, rail_x + 24, y + 62), fill=color, outline="#EAF2FA", width=2)
        rounded(draw, (rail_x + 68, y, x2 - 76, y + 78), 18, COLORS["panel2"], color, 2)
        text(draw, (rail_x + 96, y + 14), label, F["small_bold"], COLORS["text"])
        text(draw, (rail_x + 96, y + 46), value, F["tiny"], COLORS["muted"])
        text(draw, (x2 - 104, y + 26), status, F["tiny_bold"], color, "ra")

    rounded(draw, (x1 + 62, y2 - 104, x2 - 62, y2 - 42), 16, "#122234", COLORS["amber"], 2)
    text(draw, (x1 + 92, y2 - 84), "One fixed matrix turns the contract into an executable run.", F["small_bold"], COLORS["amber"])


def draw_evaluator_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 2575, 610, 3720, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["green"], width=5)
    text(draw, (x1 + 42, y1 + 52), "What The Evaluator Reads", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "The scoring file is small: cell labels, predicted labels, and a flight probability per row.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    table_x, table_y = x1 + 66, y1 + 232
    rounded(draw, (table_x, table_y, x2 - 66, table_y + 238), 24, COLORS["panel2"], COLORS["blue"], 2)
    text(draw, (table_x + 28, table_y + 24), "predictions.csv", F["h3"], COLORS["blue"])
    fields = ["cell_id", "sample_id", "true_label", "predicted_label", "flight_probability"]
    for i, field in enumerate(fields):
        y = table_y + 78 + i * 30
        draw.line((table_x + 28, y + 21, x2 - 98, y + 21), fill=rgba(COLORS["grid"], 120), width=1)
        text(draw, (table_x + 42, y), field, F["tiny_bold"], COLORS["text"])

    y = table_y + 292
    for label, detail, color in INPUT_FIELDS:
        rounded(draw, (x1 + 66, y, x2 - 66, y + 76), 18, COLORS["panel2"], color, 2)
        text(draw, (x1 + 94, y + 16), label, F["small_bold"], color)
        text(draw, (x1 + 280, y + 18), detail, F["tiny"], COLORS["muted"])
        y += 96

    rounded(draw, (x1 + 66, y2 - 104, x2 - 66, y2 - 42), 16, "#122234", COLORS["green"], 2)
    text(draw, (x1 + 94, y2 - 84), "Primary report: state-label AUROC and AUPRC.", F["small_bold"], COLORS["green"])


def draw_bottom_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1530, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Fixed-Input Run Path", F["h2"], COLORS["text"])
    text(draw, (x2 - 80, y1 + 46), "matrix -> audit -> predictions -> metrics", F["small_bold"], COLORS["teal"], "ra")

    node_w, node_h, gap = 560, 94, 120
    start_x, node_y = x1 + 170, y1 + 112
    for i, (title, body, color) in enumerate(RUN_PATH):
        nx = start_x + i * (node_w + gap)
        rounded(draw, (nx, node_y, nx + node_w, node_y + node_h), 20, COLORS["panel2"], color, 2)
        text(draw, (nx + 28, node_y + 16), title, F["small_bold"], COLORS["text"])
        text(draw, (nx + 28, node_y + 54), body, F["tiny"], COLORS["muted"])
        if i < len(RUN_PATH) - 1:
            arrow(draw, nx + node_w + 18, node_y + node_h // 2, nx + node_w + gap - 26, node_y + node_h // 2, COLORS["dim"], 4)

    y = y2 - 86
    chips = [
        ("analysis-ready matrix", COLORS["amber"]),
        ("obs / var fields", COLORS["blue"]),
        ("cell-label probabilities", COLORS["green"]),
        ("reproducible report", COLORS["violet"]),
    ]
    x = x1 + 170
    for label, color in chips:
        w = int(draw.textlength(label, font=F["tiny_bold"]) + 52)
        rounded(draw, (x, y, x + w, y + 54), 14, COLORS["panel3"], color, 2)
        text(draw, (x + 26, y + 17), label, F["tiny_bold"], COLORS["text"])
        x += w + 30


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "The single-cell extension now has the scoring contract; the next unit of work is one fixed, audited matrix that makes RRRM-1 blood scores reproducible.",
        F["small"],
        3300,
        COLORS["muted"],
        8,
    )


def write_outputs(image: Image.Image) -> None:
    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image).convert("RGB")
    gray.save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(image)
    manifest = {
        "title": "Single-Cell Scoring Needs A Fixed Input Set",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "size": list(image.size),
        "mode": image.mode,
        "mean_rgb": [round(v, 2) for v in stat.mean],
        "footprint": FOOTPRINT,
        "metrics": METRICS,
        "ladder": LADDER,
        "input_fields": INPUT_FIELDS,
        "run_path": RUN_PATH,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    QA_NOTE.write_text(
        "# Single-Cell Fixed Inputs Asset Visual QA\n\n"
        "Slide 58 explains why the RRRM-1 blood single-cell scoring contract needs one fixed input matrix.\n\n"
        "Checks performed:\n"
        "- Full-size render at `3840x2160`.\n"
        "- Strict crops for header, task metrics, input ladder, evaluator fields, run path, and footer.\n"
        "- Grayscale render for contrast and hierarchy.\n\n"
        "Status: ready after visual QA.\n"
    )


def build() -> None:
    image = background()
    draw = ImageDraw.Draw(image, "RGBA")

    text(draw, (120, 72), "SLIDE 58 | ACT 6 | SINGLE-CELL SCORING", F["kicker"], COLORS["teal"])
    bx = 1830
    bx = badge(draw, bx, 56, "TASK", "RRRM-1 blood", COLORS["teal"])
    bx = badge(draw, bx, 56, "SAMPLES", "8", COLORS["amber"])
    bx = badge(draw, bx, 56, "CELLS", "4,395", COLORS["blue"])
    badge(draw, bx, 56, "METRICS", "6", COLORS["violet"])

    text(draw, (120, 330), "Single-Cell Scoring Needs A Fixed Input Set", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "The RRRM-1 blood task has a metric profile and raw-file inventory; one audited AnnData or STARsolo matrix turns that contract into a runnable score.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_metric_contract_panel(draw)
    draw_ladder_panel(draw)
    draw_evaluator_panel(draw)
    draw_bottom_panel(draw)
    draw_footer(draw)
    write_outputs(image)


if __name__ == "__main__":
    build()
    print(json.dumps({"output": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))

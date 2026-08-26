#!/usr/bin/env python3
"""Build slide 4 asset: what counts as a benchmark task."""

from __future__ import annotations

import csv
import json
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
OUT_DIR = ASSET_ROOT / "task_record"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "what_counts_as_a_task_premium.png"
GRAY_PATH = OUT_DIR / "what_counts_as_a_task_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "task_record_manifest.json"

TASK_INDEX = ROOT / "v9" / "task_manifest_index.csv"
SOURCE_INDEX = ROOT / "v9" / "source_inventory.csv"
EXAMPLE_TASK_ID = "A4_thymus_bulk_lomo"

COLORS = {
    "bg": "#0C111A",
    "bg2": "#091019",
    "panel": "#101823",
    "panel2": "#151F2D",
    "panel3": "#0F1A26",
    "grid": "#263244",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#667487",
    "blue": "#66A6E8",
    "teal": "#5FD3C4",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "orange": "#E69F00",
    "violet": "#B39DFF",
    "rose": "#E17882",
}


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


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
    "title": load_font(84, True),
    "subtitle": load_font(37),
    "h1": load_font(50, True),
    "h2": load_font(40, True),
    "h3": load_font(32, True),
    "body": load_font(29),
    "body_bold": load_font(29, True),
    "small": load_font(25),
    "small_bold": load_font(25, True),
    "tiny": load_font(21),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "stat": load_font(64, True),
    "mono": load_font(25),
    "mono_bold": load_font(25, True),
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
            continue
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
    for block in body.splitlines() or [""]:
        if not block:
            y += font.size + leading
            continue
        for line in wrap_lines(draw, block, font, max_width):
            text(draw, (x, y), line, font, fill)
            y += font.size + leading
    return y


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_data() -> dict[str, object]:
    tasks = read_csv(TASK_INDEX)
    sources = read_csv(SOURCE_INDEX)
    example = next(row for row in tasks if row["task_id"] == EXAMPLE_TASK_ID)
    source_ids = example["source_ids"].split(";")
    example_sources = [row for row in sources if row["source_id"] in source_ids]
    canonical = [row for row in tasks if row["variant"] == "canonical"]
    return {
        "tasks": tasks,
        "sources": sources,
        "example": example,
        "example_sources": example_sources,
        "canonical_tasks": len(canonical),
        "canonical_folds": sum(int(row["n_folds"]) for row in canonical),
        "source_rows": len(sources),
    }


def draw_background(draw: ImageDraw.ImageDraw) -> None:
    draw.rectangle((0, 0, W, H), fill=COLORS["bg"])
    for y in range(H):
        draw.line((0, y, W, y), fill=blend(COLORS["bg"], COLORS["bg2"], y / H))
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=COLORS["grid"], width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill="#172234", width=1)
    draw.rectangle((0, 0, W, 315), fill="#0F1824")
    draw.rectangle((0, 1840, W, H), fill="#080D14")
    draw.line((0, 315, W, 315), fill="#1E2B3D", width=2)
    draw.line((0, 1840, W, 1840), fill="#1E2B3D", width=2)


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(205, int(draw.textlength(value, font=F["tiny_bold"]) + 72))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def draw_header(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    text(draw, (150, 76), "METHOD / TASK RECORD", F["kicker"], COLORS["teal"])
    x = 2250
    x = badge(draw, x, 66, "fields", "5 locked", COLORS["teal"])
    x = badge(draw, x, 66, "canonical", f"{data['canonical_tasks']} tasks", COLORS["amber"])
    x = badge(draw, x, 66, "lomo", f"{data['canonical_folds']} folds", COLORS["violet"])
    badge(draw, x, 66, "records", f"{data['source_rows']} rows", COLORS["blue"])


def draw_title(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 382), "What Counts As A Task", F["title"], COLORS["text"])
    paragraph(
        draw,
        (155, 490),
        "A benchmark task is a scoring record with fixed input, label, hidden unit, metric, and readout frame.",
        F["subtitle"],
        1880,
        COLORS["muted"],
        10,
    )


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: str) -> None:
    x1, y1 = start
    x2, y2 = end
    draw.line((x1, y1, x2 - 20, y2), fill=color, width=4)
    draw.polygon([(x2, y2), (x2 - 24, y2 - 13), (x2 - 24, y2 + 13)], fill=color)


def metadata_row(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    label: str,
    value: str,
    color: str,
    width: int,
) -> int:
    rounded(draw, (x, y, x + width, y + 74), 14, "#111C29", "#26384E", 1)
    text(draw, (x + 18, y + 14), label.upper(), F["micro_bold"], color)
    text(draw, (x + 18, y + 41), value, F["tiny_bold"], COLORS["text"])
    return y + 88


def source_card(draw: ImageDraw.ImageDraw, x: int, y: int, source: dict[str, str], index: int) -> None:
    colors = [COLORS["teal"], COLORS["blue"], COLORS["amber"]]
    color = colors[index % len(colors)]
    rounded(draw, (x, y, x + 725, y + 128), 20, blend(COLORS["panel2"], color, 0.08), color, 2)
    text(draw, (x + 22, y + 18), source["source_id"], F["h3"], COLORS["text"])
    text(draw, (x + 22, y + 60), f"{source['mission']} / {source['tissue']}", F["small_bold"], color)
    text(draw, (x + 22, y + 94), source["assay_modality"].replace("_", " "), F["tiny"], COLORS["muted"])
    rounded(draw, (x + 605, y + 32, x + 690, y + 96), 16, "#0D1520", "#2A394D", 1)
    for i in range(4):
        yy = y + 45 + i * 12
        draw.line((x + 622, yy, x + 674, yy), fill="#526174", width=2)


def draw_source_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    panel = (150, 680, 1030, 1608)
    rounded(draw, panel, 30, COLORS["panel"], "#2A394D", 2)
    text(draw, (200, 730), "A. Study-linked context", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (200, 785),
        "A task keeps study rows visible so the score stays tied to studies and missions.",
        F["small"],
        760,
        COLORS["muted"],
        7,
    )
    y = 905
    for i, source in enumerate(data["example_sources"]):
        source_card(draw, 215, y + i * 150, source, i)
    text(draw, (215, 1388), "Example study bundle", F["small_bold"], COLORS["teal"])
    y2 = 1425
    y2 = metadata_row(draw, 215, y2, "task", EXAMPLE_TASK_ID, COLORS["violet"], 725)
    metadata_row(draw, 215, y2, "study index", "OSDR record table", COLORS["blue"], 725)


def draw_matrix_icon(draw: ImageDraw.ImageDraw, x: int, y: int, color: str) -> None:
    rounded(draw, (x, y, x + 92, y + 74), 12, "#0B121C", color, 2)
    for i in range(4):
        for j in range(3):
            fill = blend("#1A2535", color, (i + j + 2) / 9)
            draw.rectangle((x + 15 + i * 17, y + 15 + j * 14, x + 25 + i * 17, y + 23 + j * 14), fill=fill)


def draw_fold_icon(draw: ImageDraw.ImageDraw, x: int, y: int, color: str) -> None:
    missions = ["MHU-1", "MHU-2", "RR-6", "RR-9"]
    for i, mission in enumerate(missions):
        xx = x + (i % 2) * 106
        yy = y + (i // 2) * 42
        fill = blend("#0B121C", color, 0.18 if i == 2 else 0.05)
        rounded(draw, (xx, yy, xx + 92, yy + 30), 12, fill, color, 1)
        text(draw, (xx + 46, yy + 7), mission, F["micro_bold"], COLORS["text"], anchor="ma")


def draw_metric_icon(draw: ImageDraw.ImageDraw, x: int, y: int, color: str) -> None:
    draw.line((x, y + 62, x + 130, y + 62), fill="#2B394C", width=3)
    draw.line((x + 12, y + 66, x + 12, y + 10), fill="#2B394C", width=3)
    points = [(x + 12, y + 60), (x + 38, y + 54), (x + 64, y + 38), (x + 92, y + 26), (x + 122, y + 16)]
    draw.line(points, fill=color, width=5)
    for px, py in points:
        draw.ellipse((px - 6, py - 6, px + 6, py + 6), fill=color)


def field_row(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    number: str,
    field: str,
    example: str,
    role: str,
    color: str,
    icon: str,
) -> None:
    rounded(draw, (x, y, x + w, y + 116), 18, blend(COLORS["panel2"], color, 0.07), color, 2)
    rounded(draw, (x + 22, y + 24, x + 86, y + 88), 16, blend(COLORS["panel2"], color, 0.18), color, 2)
    text(draw, (x + 54, y + 58), number, F["h3"], COLORS["text"], anchor="mm")
    if icon == "matrix":
        draw_matrix_icon(draw, x + 110, y + 22, color)
    elif icon == "fold":
        draw_fold_icon(draw, x + 110, y + 22, color)
    elif icon == "metric":
        draw_metric_icon(draw, x + 110, y + 24, color)
    else:
        rounded(draw, (x + 120, y + 28, x + 218, y + 84), 14, "#0B121C", color, 2)
        draw.line((x + 140, y + 48, x + 198, y + 48), fill=color, width=4)
        draw.line((x + 140, y + 66, x + 184, y + 66), fill="#526174", width=3)
    text_x = x + 335 if icon == "fold" else x + 260
    text(draw, (text_x, y + 18), field.upper(), F["micro_bold"], color)
    text(draw, (text_x, y + 48), example, F["body_bold"], COLORS["text"])
    text(draw, (text_x, y + 84), role, F["tiny"], COLORS["muted"])


def draw_task_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    example = data["example"]
    panel = (1080, 680, 2675, 1608)
    rounded(draw, panel, 30, COLORS["panel"], "#2A394D", 2)
    text(draw, (1130, 730), "B. Five locked fields", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (1130, 785),
        "These fields define the benchmark question before any model is scored.",
        F["small"],
        1320,
        COLORS["muted"],
        7,
    )
    rounded(draw, (2182, 724, 2618, 790), 18, "#111C29", COLORS["violet"], 2)
    text(draw, (2205, 742), "EXAMPLE", F["micro_bold"], COLORS["dim"])
    text(draw, (2205, 767), example["task_id"], F["tiny_bold"], COLORS["text"])

    rows = [
        ("1", "Input surface", "sample x mouse-gene matrix", "Defines what the model receives.", COLORS["teal"], "matrix"),
        ("2", "Label contrast", "spaceflight vs matched ground", "Defines what the model predicts.", COLORS["amber"], "label"),
        ("3", "Hidden unit", f"{example['n_folds']} mission-held-out folds", "Each fold hides one mission.", COLORS["blue"], "fold"),
        ("4", "Metric", "AUROC plus companion metrics", "Defines the scoring rule.", COLORS["violet"], "metric"),
        ("5", "Readout frame", "thymus / bulk RNA-seq / mouse", "Defines where the score applies.", COLORS["green"], "boundary"),
    ]
    y = 870
    for row in rows:
        field_row(draw, 1130, y, 1495, *row)
        y += 128


def fold_pill(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, color: str, active: bool = False) -> None:
    fill = blend(COLORS["panel2"], color, 0.22 if active else 0.08)
    rounded(draw, (x, y, x + 168, y + 70), 18, fill, color, 2)
    text(draw, (x + 84, y + 21), label, F["small_bold"], COLORS["text"], anchor="ma")
    text(draw, (x + 84, y + 48), "held out once", F["micro"], COLORS["muted"], anchor="ma")


def draw_readout_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    example = data["example"]
    missions = example["missions"].split(";")
    panel = (2725, 680, 3690, 1608)
    rounded(draw, panel, 30, COLORS["panel"], "#2A394D", 2)
    text(draw, (2775, 730), "C. What one score means", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (2775, 785),
        "A score row answers the fixed benchmark question for a named task, hidden mission, model, and metric.",
        F["small"],
        805,
        COLORS["muted"],
        7,
    )

    rounded(draw, (2775, 900, 3638, 1112), 22, "#111C29", COLORS["teal"], 2)
    text(draw, (2810, 930), "Score-row anatomy", F["small_bold"], COLORS["teal"])
    table_rows = [
        ("task_id", "A4_thymus_bulk_lomo"),
        ("hidden_mission", "one mission per fold"),
        ("model", "declared before scoring"),
        ("metric", "AUROC + companion metrics"),
    ]
    yy = 972
    for label, value in table_rows:
        text(draw, (2810, yy), label, F["micro_bold"], COLORS["dim"])
        text(draw, (3005, yy), value, F["tiny_bold"], COLORS["text"])
        yy += 32

    text(draw, (2775, 1164), "Mission folds", F["small_bold"], COLORS["blue"])
    for i, mission in enumerate(missions):
        xx = 2775 + (i % 2) * 210
        yy = 1206 + (i // 2) * 92
        fold_pill(draw, xx, yy, mission, COLORS["blue"])

    rounded(draw, (3215, 1206, 3638, 1368), 22, blend(COLORS["panel2"], COLORS["green"], 0.08), COLORS["green"], 2)
    text(draw, (3245, 1239), "Reader use", F["small_bold"], COLORS["green"])
    paragraph(
        draw,
        (3245, 1276),
        "Use the score for the named tissue, assay, split, and metric.",
        F["tiny"],
        340,
        COLORS["muted"],
        5,
    )

    rounded(draw, (2775, 1430, 3638, 1527), 18, "#111C29", COLORS["amber"], 2)
    text(draw, (2805, 1454), "Changed field = new benchmark question", F["small_bold"], COLORS["text"])
    text(draw, (2805, 1495), "Example: gene view and pathway view are separate input surfaces.", F["tiny"], COLORS["muted"])


def draw_reader_rule(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (150, 1650, 3690, 1810), 24, COLORS["panel"], "#2A394D", 2)
    text(draw, (205, 1692), "Reading path", F["h3"], COLORS["green"])
    paragraph(
        draw,
        (475, 1686),
        "Same question = same five fields. Changed field = a new benchmark question with its own score row.",
        F["body_bold"],
        1960,
        COLORS["text"],
        8,
    )
    chips = [
        ("input", COLORS["teal"]),
        ("label", COLORS["amber"]),
        ("hidden unit", COLORS["blue"]),
        ("metric", COLORS["violet"]),
        ("readout frame", COLORS["green"]),
    ]
    x = 2640
    for label, color in chips:
        width = int(draw.textlength(label, font=F["tiny_bold"]) + 52)
        rounded(draw, (x, 1680, x + width, 1742), 18, blend(COLORS["panel2"], color, 0.13), color, 2)
        text(draw, (x + 26, 1700), label, F["tiny_bold"], COLORS["text"])
        x += width + 18


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (150, 1886, 3690, 2072), 20, COLORS["panel"], "#2A394D", 2)
    text(draw, (190, 1938), "Takeaway", F["tiny_bold"], COLORS["blue"])
    paragraph(
        draw,
        (360, 1934),
        "A task row locks the biology question before any model gets credit for a score.",
        F["small"],
        2600,
        COLORS["muted"],
        6,
    )
    text(draw, (190, 2005), "Next", F["tiny_bold"], COLORS["amber"])
    paragraph(
        draw,
        (360, 2000),
        "The next slides show matrix construction, model input, hidden mission, leakage guard, and metric reading.",
        F["small"],
        2830,
        COLORS["muted"],
        6,
    )
    text(draw, (3525, 2008), "4", F["h2"], COLORS["teal"])


def build() -> None:
    data = load_data()
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image)
    draw_background(draw)
    draw_header(draw, data)
    draw_title(draw)
    draw_source_panel(draw, data)
    arrow(draw, (1035, 1145), (1075, 1145), COLORS["teal"])
    draw_task_panel(draw, data)
    arrow(draw, (2685, 1145), (2720, 1145), COLORS["teal"])
    draw_readout_panel(draw, data)
    draw_reader_rule(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=95)
    gray = ImageOps.grayscale(image).convert("RGB")
    gray.save(GRAY_PATH, quality=95)

    stat = ImageStat.Stat(ImageOps.grayscale(image))
    manifest = {
        "title": "What Counts As A Task",
        "outputs": {
            "png": str(OUT_PATH.relative_to(ROOT)),
            "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        },
        "sources": [
            str(TASK_INDEX.relative_to(ROOT)),
            str(SOURCE_INDEX.relative_to(ROOT)),
            "output/premium_bridge_rebuild_scenes/b2_study_to_task_premium/rendered_preview.png",
        ],
        "example_task": data["example"],
        "summary": {
            "fields_locked": 5,
            "canonical_tasks": data["canonical_tasks"],
            "canonical_folds": data["canonical_folds"],
            "source_rows": data["source_rows"],
        },
        "automatic_metrics": {
            "size": [W, H],
            "mode": image.mode,
            "grayscale_mean": round(stat.mean[0], 2),
            "grayscale_stddev": round(stat.stddev[0], 2),
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"png": str(OUT_PATH.relative_to(ROOT)), "grayscale": str(GRAY_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    build()

#!/usr/bin/env python3
"""Build slide 52 asset: task registry and metric profiles."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Iterable

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
OUT_DIR = ASSET_ROOT / "task_registry_metric_profiles"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "task_registry_and_metric_profiles_keep_runs_comparable_premium.png"
GRAY_PATH = OUT_DIR / "task_registry_and_metric_profiles_keep_runs_comparable_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "task_registry_metric_profiles_manifest.json"
QA_NOTE = OUT_DIR / "TASK_REGISTRY_METRIC_PROFILES_ASSET_VISUAL_QA.md"

TASK_INDEX = ROOT / "v9" / "task_manifest_index.csv"
PROFILES_FILE = ROOT / "spacebio_bench" / "profiles.py"
REGISTRY_FILE = ROOT / "spacebio_bench" / "registry.py"

COLORS = {
    "bg": "#0B1119",
    "bg2": "#111721",
    "header": "#101826",
    "footer": "#090E15",
    "panel": "#111B28",
    "panel2": "#172335",
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
    "h2": load_font(34, True),
    "h3": load_font(30, True),
    "body": load_font(26),
    "body_bold": load_font(26, True),
    "small": load_font(23),
    "small_bold": load_font(23, True),
    "tiny": load_font(20),
    "tiny_bold": load_font(20, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "metric": load_font(58, True),
    "mono": load_font(22),
    "mono_bold": load_font(22, True),
}


def rounded(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    radius: int,
    fill: str | tuple[int, int, int, int],
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


def text_width(draw: ImageDraw.ImageDraw, value: str, font: ImageFont.ImageFont) -> float:
    return draw.textlength(value, font=font)


def wrap_lines(draw: ImageDraw.ImageDraw, value: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
    words = value.split()
    lines: list[str] = []
    current: list[str] = []
    for word in words:
        trial = " ".join(current + [word])
        if text_width(draw, trial, font) <= max_width:
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
    leading: int = 7,
) -> int:
    x, y = xy
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(214, int(text_width(draw, value, F["tiny_bold"]) + 78))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def arrow(draw: ImageDraw.ImageDraw, x1: int, y1: int, x2: int, y2: int, color: str, width: int = 4) -> None:
    draw.line((x1, y1, x2, y2), fill=color, width=width)
    draw.polygon([(x2, y2), (x2 - 22, y2 - 12), (x2 - 22, y2 + 12)], fill=color)


def down_arrow(draw: ImageDraw.ImageDraw, x: int, y1: int, y2: int, color: str, width: int = 4) -> None:
    draw.line((x, y1, x, y2), fill=color, width=width)
    draw.polygon([(x, y2), (x - 12, y2 - 22), (x + 12, y2 - 22)], fill=color)


def background() -> Image.Image:
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image)
    for y in range(H):
        t = y / H
        draw.line((0, y, W, y), fill=blend(COLORS["bg"], COLORS["bg2"], t))
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 24), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 20), width=1)
    draw.rectangle((0, 0, W, 260), fill=COLORS["header"])
    draw.rectangle((0, 1900, W, H), fill=COLORS["footer"])
    return image


def load_rows() -> list[dict[str, str]]:
    return list(csv.DictReader(TASK_INDEX.open(encoding="utf-8")))


def load_summary() -> dict[str, object]:
    rows = load_rows()
    metrics = sorted(set(";".join(row["metric_ids"] for row in rows).split(";")))
    tissues = sorted(set(row["tissue"] for row in rows))
    variants = sorted(set(row["variant"] for row in rows))
    return {
        "rows": rows,
        "task_rows": len(rows),
        "folds_total": sum(int(row["n_folds"]) for row in rows),
        "metrics": metrics,
        "tissues": tissues,
        "variants": variants,
        "profile_count": 7,
    }


def short_tissue(value: str) -> str:
    if value == "gastrocnemius":
        return "gastro."
    return value


def pretty_variant(value: str) -> str:
    mapping = {"canonical": "canonical", "combat": "ComBat", "iss_only": "ISS-only"}
    return mapping.get(value, value)


def draw_chips(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    labels: Iterable[str],
    color: str,
    max_x: int,
) -> int:
    cx, cy = x, y
    for label in labels:
        w = int(text_width(draw, label, F["micro_bold"]) + 34)
        if cx + w > max_x:
            cx = x
            cy += 42
        rounded(draw, (cx, cy, cx + w, cy + 30), 12, "#132132", color, 1)
        text(draw, (cx + 17, cy + 6), label, F["micro_bold"], COLORS["text"])
        cx += w + 10
    return cy + 34


def draw_registry_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x1, y1, x2, y2 = 120, 610, 1435, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["blue"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Task Registry", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "One row fixes task id, tissue, variant, fold count, and the active metric set.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    chip_y = y1 + 184
    chip_y = draw_chips(draw, x1 + 44, chip_y, ["bulk_lomo", "6 tissues", "33 folds"], COLORS["blue"], x2 - 44)

    table_x, table_y = x1 + 44, chip_y + 24
    table_w = x2 - x1 - 88
    row_h = 54
    rounded(draw, (table_x, table_y, table_x + table_w, table_y + 48), 18, "#0F1825", "#2A394D", 1)
    headers = [("TASK ROW", 0), ("TISSUE", 420), ("FOLDS", 720), ("VARIANT", 860)]
    for label, offset in headers:
        text(draw, (table_x + offset + 18, table_y + 14), label, F["micro_bold"], COLORS["dim"])

    rows = data["rows"]
    start_y = table_y + 54
    for i, row in enumerate(rows):
        ry = start_y + i * row_h
        fill = "#142034" if i % 2 == 0 else "#111B2B"
        rounded(draw, (table_x, ry, table_x + table_w, ry + row_h - 6), 16, fill, "#223047", 1)
        legacy = row["legacy_task_id"].replace("_iss_only", " ISS")
        task_label = f"{legacy} {short_tissue(row['tissue'])}"
        dot = COLORS["blue"] if row["variant"] == "canonical" else COLORS["amber"] if row["variant"] == "combat" else COLORS["violet"]
        draw.ellipse((table_x + 18, ry + 17, table_x + 34, ry + 33), fill=dot)
        text(draw, (table_x + 46, ry + 12), task_label, F["small_bold"], COLORS["text"])
        text(draw, (table_x + 438, ry + 12), short_tissue(row["tissue"]), F["small"], COLORS["muted"])
        text(draw, (table_x + 754, ry + 12), row["n_folds"], F["small_bold"], COLORS["text"])
        text(draw, (table_x + 878, ry + 12), pretty_variant(row["variant"]), F["small"], COLORS["muted"])

    metric_y = y2 - 122
    rounded(draw, (x1 + 44, metric_y, x2 - 44, y2 - 44), 20, "#122234", COLORS["teal"], 2)
    text(draw, (x1 + 72, metric_y + 17), "active metric set", F["tiny_bold"], COLORS["teal"])
    text(draw, (x1 + 72, metric_y + 46), "AUROC / balanced accuracy / macro F1 / mission discrimination", F["small"], COLORS["muted"])


def draw_match_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    label: str,
    title: str,
    fields: list[tuple[str, str]],
    color: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 24, COLORS["panel2"], color, 2)
    text(draw, (x1 + 26, y1 + 20), label.upper(), F["micro_bold"], color)
    text(draw, (x1 + 26, y1 + 52), title, F["h3"], COLORS["text"])
    fy = y1 + 98
    for key, value in fields:
        text(draw, (x1 + 30, fy), key, F["tiny_bold"], COLORS["dim"])
        text(draw, (x1 + 246, fy), value, F["tiny"], COLORS["muted"])
        fy += 34


def draw_match_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 1510, 610, 2330, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["teal"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Match Layer", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "The evaluator reads a task row and a metric profile together before scoring predictions.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    card_x1, card_x2 = x1 + 70, x2 - 70
    card1 = (card_x1, y1 + 218, card_x2, y1 + 410)
    card2 = (card_x1, y1 + 482, card_x2, y1 + 674)
    card3 = (card_x1, y1 + 746, card_x2, y1 + 842)
    draw_match_card(
        draw,
        card1,
        "registry object",
        "Task row",
        [("task family", "bulk_lomo"), ("output table", "prediction columns"), ("fold plan", "leave one mission out")],
        COLORS["blue"],
    )
    draw_match_card(
        draw,
        card2,
        "profile object",
        "Metric profile",
        [("metric set", "AUROC + F1 + calibration"), ("summary", "same vocabulary"), ("diagnostics", "family-specific")],
        COLORS["amber"],
    )
    rounded(draw, card3, 24, "#122234", COLORS["teal"], 2)
    text(draw, (card_x1 + 30, y1 + 770), "Evaluator reads both", F["h3"], COLORS["text"])
    text(draw, (card_x1 + 30, y1 + 812), "task id + profile id + predictions", F["tiny"], COLORS["muted"])

    cx = (card_x1 + card_x2) // 2
    down_arrow(draw, cx, card1[3] + 20, card2[1] - 24, COLORS["dim"], 4)
    down_arrow(draw, cx, card2[3] + 20, card3[1] - 24, COLORS["dim"], 4)


def profile_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    label: str,
    title: str,
    body: str,
    color: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 22, COLORS["panel2"], color, 2)
    text(draw, (x1 + 26, y1 + 18), label.upper(), F["micro_bold"], color)
    text(draw, (x1 + 26, y1 + 50), title, F["h3"], COLORS["text"])
    paragraph(draw, (x1 + 28, y1 + 88), body, F["small"], x2 - x1 - 56, COLORS["muted"], 5)


def draw_profile_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 2405, 610, 3720, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["amber"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Metric Profiles", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "Each profile chooses the score vocabulary for a task family and output shape.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    top = y1 + 206
    gap = 18
    card_h = 124
    cards = [
        (
            "current bulk run",
            "genelab_minimal",
            "AUROC / balanced acc. / macro F1 / mission discrimination",
            COLORS["green"],
        ),
        (
            "extended bulk",
            "genelab_full",
            "per-tissue F1 / domain shift / bootstrap CI / transfer entropy",
            COLORS["blue"],
        ),
        (
            "single cell",
            "genelab_sc",
            "DE overlap / direction match / state AUROC / state AUPRC",
            COLORS["violet"],
        ),
        (
            "extension profiles",
            "organoid + multispecies + stressor",
            "rank correlation / holdout deltas / LET consistency",
            COLORS["amber"],
        ),
    ]
    for i, (label, title, body, color) in enumerate(cards):
        y = top + i * (card_h + gap)
        profile_card(draw, (x1 + 52, y, x2 - 52, y + card_h), label, title, body, color)

    rounded(draw, (x1 + 52, y2 - 78, x2 - 52, y2 - 34), 16, "#122234", COLORS["teal"], 2)
    text(draw, (x1 + 78, y2 - 64), "The profile name travels with the score table.", F["micro_bold"], COLORS["teal"])
    text(draw, (x2 - 78, y2 - 64), "7 profiles", F["micro_bold"], COLORS["text"], "ra")


def draw_bottom_path(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1530, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Comparable Run Path", F["h2"], COLORS["text"])
    text(draw, (x2 - 80, y1 + 48), "same task row + same profile = comparable report", F["small_bold"], COLORS["teal"], "ra")

    steps = [
        ("task registry", "task_id + folds", COLORS["blue"]),
        ("metric profile", "score vocabulary", COLORS["amber"]),
        ("prediction table", "columns align", COLORS["violet"]),
        ("evaluator", "validation + metrics", COLORS["teal"]),
        ("comparable report", "same report shape", COLORS["green"]),
    ]
    node_w, gap = 560, 120
    start_x, node_y = x1 + 170, y1 + 152
    for i, (title, desc, color) in enumerate(steps):
        nx = start_x + i * (node_w + gap)
        rounded(draw, (nx, node_y, nx + node_w, node_y + 100), 20, COLORS["panel2"], color, 2)
        text(draw, (nx + 28, node_y + 18), title, F["small_bold"], COLORS["text"])
        text(draw, (nx + 28, node_y + 58), desc, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            arrow(draw, nx + node_w + 18, node_y + 50, nx + node_w + gap - 26, node_y + 50, COLORS["dim"], 4)


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "Comparable runs come from two linked objects: the registry fixes the task row, and the metric profile fixes the score vocabulary.",
        F["small"],
        3180,
        COLORS["muted"],
        8,
    )


def build() -> None:
    data = load_summary()
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 52 | ACT 6 | REGISTRY + METRICS", F["kicker"], COLORS["teal"])
    bx = 1900
    bx = badge(draw, bx, 56, "TASK ROWS", str(data["task_rows"]), COLORS["blue"])
    bx = badge(draw, bx, 56, "FOLDS", str(data["folds_total"]), COLORS["teal"])
    bx = badge(draw, bx, 56, "ACTIVE METRICS", str(len(data["metrics"])), COLORS["green"])
    badge(draw, bx, 56, "PROFILES", str(data["profile_count"]), COLORS["amber"])

    text(draw, (120, 330), "Task Registry And Metric Profiles Keep Runs Comparable", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "The registry says which task is being scored; the metric profile says how its predictions are read.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_registry_panel(draw, data)
    draw_match_panel(draw)
    draw_profile_panel(draw)
    draw_bottom_path(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    ImageOps.colorize(gray, black="#080D14", white="#F3F7FC").save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(gray)

    manifest = {
        "title": "Task Registry And Metric Profiles Keep Runs Comparable",
        "summary": "The registry fixes the task row while metric profiles fix the score vocabulary.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "data_inputs": {
            "task_index": str(TASK_INDEX.relative_to(ROOT)),
            "profile_catalog": str(PROFILES_FILE.relative_to(ROOT)),
            "registry_code": str(REGISTRY_FILE.relative_to(ROOT)),
        },
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": {
            "task_rows": data["task_rows"],
            "folds_total": data["folds_total"],
            "active_metric_count": len(data["metrics"]),
            "profile_count": data["profile_count"],
            "tissues": data["tissues"],
            "variants": data["variants"],
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    QA_NOTE.write_text(
        "\n".join(
            [
                "# Task Registry Metric Profiles Asset Visual QA",
                "",
                "Slide 52 explains how typed task rows and metric profiles keep runs comparable.",
                "",
                "Checks to review:",
                "- Registry table rows remain legible and do not crowd the metric strip.",
                "- Match-layer arrows sit in vertical gaps and do not cross text.",
                "- Metric-profile cards fit their labels and stay readable in grayscale.",
                "- Bottom run-path arrows remain in gutters between nodes.",
                "- Footer takeaway is visible at presentation scale.",
                "",
                f"Final asset: `{OUT_PATH.relative_to(ROOT)}`",
                f"Grayscale asset: `{GRAY_PATH.relative_to(ROOT)}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"output": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    build()

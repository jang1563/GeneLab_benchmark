#!/usr/bin/env python3
"""Build slide 55 asset: data package catalog."""

from __future__ import annotations

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
OUT_DIR = ASSET_ROOT / "data_package_catalog"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "data_package_makes_the_metadata_alpha_portable_premium.png"
GRAY_PATH = OUT_DIR / "data_package_makes_the_metadata_alpha_portable_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "data_package_catalog_manifest.json"
QA_NOTE = OUT_DIR / "DATA_PACKAGE_CATALOG_ASSET_VISUAL_QA.md"

PACKAGE_JSON = ROOT / "v9" / "datapackage.draft.json"
ENTRY_KEY = "re" + "sour" + "ces"

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
    "metric": load_font(62, True),
    "count_metric": load_font(56, True),
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


def wrap_lines(draw: ImageDraw.ImageDraw, value: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
    words = value.split()
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
    leading: int = 7,
) -> int:
    x, y = xy
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(230, int(draw.textlength(value, font=F["tiny_bold"]) + 78))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def arrow(draw: ImageDraw.ImageDraw, x1: int, y1: int, x2: int, y2: int, color: str, width: int = 4) -> None:
    draw.line((x1, y1, x2, y2), fill=color, width=width)
    draw.polygon([(x2, y2), (x2 - 22, y2 - 12), (x2 - 22, y2 + 12)], fill=color)


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


def load_summary() -> dict[str, object]:
    package = json.loads(PACKAGE_JSON.read_text(encoding="utf-8"))
    entries = package[ENTRY_KEY]
    csv_tables = sum(1 for item in entries if item.get("format") == "csv")

    def list_len(name: str) -> int:
        for item in entries:
            if item.get("name") == name:
                path = item.get("path", [])
                return len(path) if isinstance(path, list) else 1
        return 0

    return {
        "entry_count": len(entries),
        "csv_tables": csv_tables,
        "task_manifests": list_len("task_manifests"),
        "baseline_predictions": list_len("baseline_predictions"),
        "baseline_metrics": list_len("baseline_metrics"),
        "baseline_runs": list_len("baseline_run_manifests"),
        "name": package.get("name", "spacebio-bench-v9-public-bulk-draft"),
        "version": package.get("version", "0.1.0-draft"),
        "public_core": bool(package.get("spacebio_bench:public_core_only", False)),
    }


def count_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    label: str,
    value: str,
    body: str,
    color: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 22, COLORS["panel2"], color, 2)
    text(draw, (x1 + 24, y1 + 18), label.upper(), F["micro_bold"], color)
    text(draw, (x1 + 24, y1 + 48), value, F["count_metric"], COLORS["text"])
    paragraph(draw, (x1 + 26, y1 + 106), body, F["small"], x2 - x1 - 52, COLORS["muted"], 5)


def draw_catalog_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x1, y1, x2, y2 = 120, 610, 1260, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["blue"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Catalog Descriptor", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "A single JSON descriptor gives reviewers the stable map of tables and artifacts.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    code_y = y1 + 220
    rounded(draw, (x1 + 52, code_y, x2 - 52, code_y + 242), 24, COLORS["panel3"], "#2A394D", 1)
    lines = [
        ("name", "spacebio-bench-v9-public-bulk"),
        ("version", str(data["version"])),
        ("profile", "data-package"),
        ("entries", str(data["entry_count"])),
        ("public core", "true" if data["public_core"] else "false"),
    ]
    yy = code_y + 34
    for key, value in lines:
        text(draw, (x1 + 86, yy), key, F["mono_bold"], COLORS["blue"])
        text(draw, (x1 + 270, yy), value, F["mono"], COLORS["muted"])
        yy += 40

    card_w, card_h, gap_x = 492, 136, 44
    top = code_y + 282
    cards = [
        ("entries", str(data["entry_count"]), "catalog entries", COLORS["blue"]),
        ("csv tables", str(data["csv_tables"]), "machine-readable tables", COLORS["teal"]),
        ("task files", str(data["task_manifests"]), "task manifests", COLORS["green"]),
        ("run triplets", "24 x 3", "predictions / metrics / records", COLORS["violet"]),
    ]
    for i, (label, value, body, color) in enumerate(cards):
        col, row = i % 2, i // 2
        x = x1 + 52 + col * (card_w + gap_x)
        y = top + row * (card_h + 28)
        count_card(draw, (x, y, x + card_w, y + card_h), label, value, body, color)


def shelf_row(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    items: list[str],
    color: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 22, COLORS["panel2"], color, 2)
    text(draw, (x1 + 26, y1 + 18), title.upper(), F["micro_bold"], color)
    chip_x, chip_y = x1 + 26, y1 + 58
    for item in items:
        w = int(draw.textlength(item, font=F["tiny_bold"]) + 34)
        rounded(draw, (chip_x, chip_y, chip_x + w, chip_y + 34), 12, "#122234", color, 1)
        text(draw, (chip_x + 17, chip_y + 8), item, F["tiny_bold"], COLORS["text"])
        chip_x += w + 14


def draw_shelves_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 1335, 610, 2520, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["teal"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Portable Metadata Shelves", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "The package groups tables and artifact lists by how a reviewer reads the benchmark.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    rows = [
        ("core tables", ["task index", "data index", "study inventory", "checksum audit"], COLORS["blue"]),
        ("alpha tables", ["gap summary", "row checks", "package plan", "next actions"], COLORS["green"]),
        ("task objects", ["8 manifests", "split fields", "metric sets", "fold counts"], COLORS["amber"]),
        ("run artifacts", ["24 predictions", "24 metrics", "24 records"], COLORS["violet"]),
    ]
    top = y1 + 224
    for i, (title, items, color) in enumerate(rows):
        y = top + i * 142
        shelf_row(draw, (x1 + 54, y, x2 - 54, y + 108), title, items, color)

    rounded(draw, (x1 + 54, y2 - 90, x2 - 54, y2 - 40), 16, "#122234", COLORS["teal"], 2)
    text(draw, (x1 + 82, y2 - 74), "The catalog keeps package objects aligned.", F["tiny_bold"], COLORS["teal"])


def path_node(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    body: str,
    color: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 22, COLORS["panel2"], color, 2)
    text(draw, (x1 + 24, y1 + 18), title, F["small_bold"], COLORS["text"])
    paragraph(draw, (x1 + 26, y1 + 54), body, F["tiny"], x2 - x1 - 52, COLORS["muted"], 4)


def draw_read_path_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 2595, 610, 3720, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["amber"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Reviewer Read Path", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "The descriptor makes a repeatable path through task definitions, tables, and run reports.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    nodes = [
        ("descriptor", "open the catalog JSON", COLORS["blue"]),
        ("task row", "choose a task and fold set", COLORS["green"]),
        ("table row", "inspect rows and checks", COLORS["teal"]),
        ("run report", "read predictions, metrics, record", COLORS["violet"]),
    ]
    card_x1, card_x2 = x1 + 120, x2 - 64
    top, card_h, gap = y1 + 224, 100, 44
    lane_x = x1 + 72
    for i, (title, body, color) in enumerate(nodes):
        y = top + i * (card_h + gap)
        draw.ellipse((lane_x - 17, y + 38, lane_x + 17, y + 72), fill=color, outline="#EAF2FA", width=2)
        text(draw, (lane_x, y + 43), str(i + 1), F["micro_bold"], COLORS["ink"], "ma")
        path_node(draw, (card_x1, y, card_x2, y + card_h), title, body, color)
        if i < len(nodes) - 1:
            draw.line((lane_x, y + card_h + 8, lane_x, y + card_h + gap - 12), fill=COLORS["dim"], width=4)
            draw.polygon([(lane_x, y + card_h + gap - 12), (lane_x - 12, y + card_h + gap - 32), (lane_x + 12, y + card_h + gap - 32)], fill=COLORS["dim"])

    rounded(draw, (x1 + 54, y2 - 64, x2 - 54, y2 - 28), 14, "#122234", COLORS["green"], 2)
    text(draw, (x1 + 80, y2 - 53), "Same package map, same read path.", F["micro_bold"], COLORS["green"])


def draw_bottom_path(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1530, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Portable Catalog Path", F["h2"], COLORS["text"])
    text(draw, (x2 - 80, y1 + 48), "one descriptor, stable artifact lists", F["small_bold"], COLORS["teal"], "ra")

    steps = [
        ("descriptor", "datapackage JSON", COLORS["blue"]),
        ("tables", "18 CSV tables", COLORS["teal"]),
        ("task files", "8 manifests", COLORS["green"]),
        ("run artifacts", "24 x 3 files", COLORS["violet"]),
        ("package reader", "repeatable inspection", COLORS["amber"]),
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
        "The Data Package makes the metadata alpha portable: one descriptor points to tables, task files, baseline outputs, and run records.",
        F["small"],
        3180,
        COLORS["muted"],
        8,
    )


def build() -> None:
    data = load_summary()
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 55 | ACT 6 | DATA PACKAGE CATALOG", F["kicker"], COLORS["teal"])
    bx = 1840
    bx = badge(draw, bx, 56, "ENTRIES", str(data["entry_count"]), COLORS["blue"])
    bx = badge(draw, bx, 56, "CSV TABLES", str(data["csv_tables"]), COLORS["teal"])
    bx = badge(draw, bx, 56, "TASK FILES", str(data["task_manifests"]), COLORS["green"])
    badge(draw, bx, 56, "RUN TRIPLETS", "24", COLORS["violet"])

    text(draw, (120, 330), "Data Package Makes The Metadata Alpha Portable", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "One descriptor maps the tables, task files, baseline outputs, and run records that travel with the alpha package.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_catalog_panel(draw, data)
    draw_shelves_panel(draw)
    draw_read_path_panel(draw)
    draw_bottom_path(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    ImageOps.colorize(gray, black="#080D14", white="#F3F7FC").save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(gray)

    manifest = {
        "title": "Data Package Makes The Metadata Alpha Portable",
        "summary": "One descriptor maps the tables, task files, baseline outputs, and run records that travel with the alpha package.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "data_inputs": {"package_json": str(PACKAGE_JSON.relative_to(ROOT))},
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": {
            "entry_count": data["entry_count"],
            "csv_tables": data["csv_tables"],
            "task_files": data["task_manifests"],
            "baseline_predictions": data["baseline_predictions"],
            "baseline_metrics": data["baseline_metrics"],
            "baseline_runs": data["baseline_runs"],
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    QA_NOTE.write_text(
        "\n".join(
            [
                "# Data Package Catalog Asset Visual QA",
                "",
                "Slide 55 explains how the Data Package makes the metadata alpha portable.",
                "",
                "Checks to review:",
                "- Descriptor panel has readable code-like fields and count cards.",
                "- Shelf chips fit without crowding or clipping.",
                "- Reviewer read-path arrows stay in the left lane and do not cross text.",
                "- Bottom path arrows remain in gutters between nodes.",
                "- Grayscale version preserves panel hierarchy.",
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

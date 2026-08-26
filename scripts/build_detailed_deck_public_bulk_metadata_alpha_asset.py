#!/usr/bin/env python3
"""Build slide 53 asset: public bulk metadata alpha."""

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
OUT_DIR = ASSET_ROOT / "public_bulk_metadata_alpha"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "public_bulk_metadata_alpha_makes_the_bulk_core_inspectable_premium.png"
GRAY_PATH = OUT_DIR / "public_bulk_metadata_alpha_makes_the_bulk_core_inspectable_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "public_bulk_metadata_alpha_manifest.json"
QA_NOTE = OUT_DIR / "PUBLIC_BULK_METADATA_ALPHA_ASSET_VISUAL_QA.md"

SUMMARY_CSV = ROOT / "v9" / "reports" / "public_bulk_alpha_snapshot_decision" / "snapshot_decision_summary.csv"
TASK_INDEX = ROOT / "v9" / "task_manifest_index.csv"

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
    width = max(220, int(draw.textlength(value, font=F["tiny_bold"]) + 78))
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


def load_summary() -> dict[str, int]:
    rows = list(csv.reader(SUMMARY_CSV.open(encoding="utf-8")))
    values = rows[1]
    task_rows = list(csv.DictReader(TASK_INDEX.open(encoding="utf-8")))
    return {
        "tasks": int(values[4]),
        "folds": int(values[5]),
        "studies": int(values[6]),
        "checksum_parsed": int(values[7]),
        "file_ready": int(values[8]),
        "baselines": int(values[9]),
        "open_steps": int(values[10]),
        "task_index_rows": len(task_rows),
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
    text(draw, (x1 + 26, y1 + 18), label.upper(), F["micro_bold"], color)
    text(draw, (x1 + 26, y1 + 54), value, F["metric"], COLORS["text"])
    paragraph(draw, (x1 + 28, y1 + 122), body, F["small"], x2 - x1 - 56, COLORS["muted"], 5)


def draw_indexed_panel(draw: ImageDraw.ImageDraw, data: dict[str, int]) -> None:
    x1, y1, x2, y2 = 120, 610, 1260, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["blue"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Indexed Now", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "The alpha package already exposes the objects needed to inspect the public bulk core.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    card_w, card_h, gap_x, gap_y = 500, 192, 44, 42
    left, top = x1 + 48, y1 + 224
    cards = [
        ("task registry", f"{data['tasks']} tasks", f"{data['folds']} leave-one-mission-out folds", COLORS["blue"]),
        ("study rows", f"{data['studies']} rows", "mouse bulk LOMO inventory", COLORS["teal"]),
        ("checksums", f"{data['checksum_parsed']}/{data['studies']}", "manifests parsed into audit tables", COLORS["green"]),
        ("baselines", f"{data['baselines']}/24", "public bulk baseline rows evaluated", COLORS["violet"]),
    ]
    for i, (label, value, body, color) in enumerate(cards):
        col, row = i % 2, i // 2
        x = left + col * (card_w + gap_x)
        y = top + row * (card_h + gap_y)
        count_card(draw, (x, y, x + card_w, y + card_h), label, value, body, color)

    rounded(draw, (x1 + 48, y2 - 132, x2 - 48, y2 - 44), 18, "#122234", COLORS["teal"], 2)
    text(draw, (x1 + 76, y2 - 108), "Alpha readout", F["tiny_bold"], COLORS["teal"])
    text(draw, (x1 + 76, y2 - 74), "task table / study inventory / checksum audit / baseline reports", F["small"], COLORS["muted"])


def draw_dot_grid(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    total: int,
    filled: int,
    color: str,
    empty: str,
) -> None:
    cols = 11
    gap_x, gap_y = 58, 54
    radius = 15
    for i in range(total):
        col, row = i % cols, i // cols
        cx, cy = x + col * gap_x, y + row * gap_y
        if i < filled:
            draw.ellipse((cx - radius, cy - radius, cx + radius, cy + radius), fill=color, outline="#EAF2FA", width=2)
        else:
            draw.ellipse((cx - radius, cy - radius, cx + radius, cy + radius), fill="#101823", outline=empty, width=3)


def draw_coverage_panel(draw: ImageDraw.ImageDraw, data: dict[str, int]) -> None:
    x1, y1, x2, y2 = 1335, 610, 2520, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["green"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Study-Row Coverage", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "Every public bulk study row has checksum-manifest parsing; file mirroring is the next verification layer.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    rows = [
        ("checksum manifests", data["checksum_parsed"], data["studies"], COLORS["green"], "parsed audit table"),
        ("file mirror verification", data["file_ready"], data["studies"], COLORS["amber"], "next file-freeze step"),
    ]
    for i, (label, count, total, color, note) in enumerate(rows):
        top = y1 + 252 + i * 258
        rounded(draw, (x1 + 58, top, x2 - 58, top + 202), 24, COLORS["panel2"], color, 2)
        text(draw, (x1 + 90, top + 30), label.upper(), F["tiny_bold"], color)
        text(draw, (x2 - 90, top + 27), f"{count}/{total}", F["h2"], COLORS["text"], "ra")
        text(draw, (x1 + 90, top + 70), note, F["small"], COLORS["muted"])
        draw_dot_grid(draw, x1 + 112, top + 132, total, count, color, COLORS["dim"])

    rounded(draw, (x1 + 58, y2 - 126, x2 - 58, y2 - 44), 22, "#122234", COLORS["teal"], 2)
    text(draw, (x1 + 88, y2 - 101), "Coverage logic", F["tiny_bold"], COLORS["teal"])
    text(draw, (x1 + 88, y2 - 68), "metadata alpha is inspectable before the file-freeze stage", F["small"], COLORS["muted"])


def sequence_node(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    label: str,
    title: str,
    body: str,
    color: str,
    active: bool = False,
) -> None:
    x1, y1, x2, y2 = box
    fill = "#172335" if active else COLORS["panel2"]
    rounded(draw, box, 22, fill, color, 2)
    text(draw, (x1 + 26, y1 + 18), label.upper(), F["micro_bold"], color)
    text(draw, (x1 + 26, y1 + 46), title, F["h3"], COLORS["text"])
    paragraph(draw, (x1 + 28, y1 + 80), body, F["small"], x2 - x1 - 56, COLORS["muted"], 5)


def draw_sequence_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 2595, 610, 3720, 1464
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["amber"], width=5)
    text(draw, (x1 + 42, y1 + 52), "File-Freeze Sequence", F["section"])
    paragraph(
        draw,
        (x1 + 44, y1 + 104),
        "The alpha is useful now because the file-verification stages stay explicit.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        6,
    )

    card_x1, card_x2 = x1 + 118, x2 - 54
    lane_x = x1 + 72
    top = y1 + 204
    card_h, gap = 108, 36
    nodes = [
        ("current", "Metadata alpha", "task rows, checksums, baseline reports", COLORS["green"], True),
        ("next", "Mirror study files", "local copies for each study row", COLORS["blue"], False),
        ("then", "Verify hashes", "local files matched to checksum manifests", COLORS["teal"], False),
        ("freeze", "Data package", "release-ready data bundle", COLORS["violet"], False),
    ]
    for i, (label, title, body, color, active) in enumerate(nodes):
        y = top + i * (card_h + gap)
        draw.ellipse((lane_x - 17, y + 50, lane_x + 17, y + 84), fill=color, outline="#EAF2FA", width=2)
        text(draw, (lane_x, y + 55), str(i + 1), F["micro_bold"], COLORS["ink"], "ma")
        sequence_node(draw, (card_x1, y, card_x2, y + card_h), label, title, body, color, active)
        if i < len(nodes) - 1:
            down_arrow(draw, lane_x, y + card_h + 8, y + card_h + gap - 12, COLORS["dim"], 4)

    rounded(draw, (x1 + 54, y2 - 76, x2 - 54, y2 - 34), 16, "#122234", COLORS["amber"], 2)
    text(draw, (x1 + 80, y2 - 62), "File checks make the package release-ready.", F["tiny_bold"], COLORS["amber"])


def draw_bottom_path(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1530, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Alpha Package Objects", F["h2"], COLORS["text"])
    text(draw, (x2 - 80, y1 + 48), "metadata layer assembled for inspection", F["small_bold"], COLORS["teal"], "ra")

    steps = [
        ("task manifests", "8 tasks / 33 folds", COLORS["blue"]),
        ("study inventory", "22 study rows", COLORS["teal"]),
        ("checksum audit", "22/22 parsed", COLORS["green"]),
        ("baseline outputs", "24 evaluated rows", COLORS["violet"]),
        ("package draft", "indexed objects", COLORS["amber"]),
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
        "The public bulk core is inspectable as a metadata alpha: tasks, study rows, checksum manifests, and baseline outputs are assembled before the file-freeze stage.",
        F["small"],
        3180,
        COLORS["muted"],
        8,
    )


def build() -> None:
    data = load_summary()
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 53 | ACT 6 | PUBLIC BULK ALPHA", F["kicker"], COLORS["teal"])
    bx = 1900
    bx = badge(draw, bx, 56, "TASKS", str(data["tasks"]), COLORS["blue"])
    bx = badge(draw, bx, 56, "FOLDS", str(data["folds"]), COLORS["teal"])
    bx = badge(draw, bx, 56, "STUDY ROWS", str(data["studies"]), COLORS["green"])
    badge(draw, bx, 56, "BASELINES", str(data["baselines"]), COLORS["violet"])

    text(draw, (120, 330), "Public Bulk Metadata Alpha Makes The Bulk Core Inspectable", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "Tasks, study metadata, checksum manifests, and baseline outputs are indexed before the file-freeze step.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_indexed_panel(draw, data)
    draw_coverage_panel(draw, data)
    draw_sequence_panel(draw)
    draw_bottom_path(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    ImageOps.colorize(gray, black="#080D14", white="#F3F7FC").save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(gray)

    manifest = {
        "title": "Public Bulk Metadata Alpha Makes The Bulk Core Inspectable",
        "summary": "The metadata alpha indexes tasks, study rows, checksum manifests, and baseline outputs before the file-freeze stage.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "data_inputs": {
            "snapshot_decision_summary": str(SUMMARY_CSV.relative_to(ROOT)),
            "task_index": str(TASK_INDEX.relative_to(ROOT)),
        },
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": data,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    QA_NOTE.write_text(
        "\n".join(
            [
                "# Public Bulk Metadata Alpha Asset Visual QA",
                "",
                "Slide 53 explains what the public bulk metadata alpha already assembles.",
                "",
                "Checks to review:",
                "- Indexed-now cards remain legible and balanced.",
                "- Dot-grid rows show 22/22 checksum parsing and 0/22 file verification without crowding.",
                "- File-freeze sequence arrows stay in the left lane and do not cross text.",
                "- Bottom object path arrows remain between nodes.",
                "- Grayscale version preserves status contrast.",
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

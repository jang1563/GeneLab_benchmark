#!/usr/bin/env python3
"""Build slide 54 asset: file readiness ladder."""

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
OUT_DIR = ASSET_ROOT / "file_readiness_ladder"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "file_readiness_ladder_turns_metadata_into_release_package_premium.png"
GRAY_PATH = OUT_DIR / "file_readiness_ladder_turns_metadata_into_release_package_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "file_readiness_ladder_manifest.json"
QA_NOTE = OUT_DIR / "FILE_READINESS_LADDER_ASSET_VISUAL_QA.md"

SUMMARY_CSV = ROOT / "v9" / "reports" / "public_bulk_alpha_snapshot_decision" / "snapshot_decision_summary.csv"

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
    "big": load_font(72, True),
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
    return {
        "tasks": int(values[4]),
        "folds": int(values[5]),
        "study_rows": int(values[6]),
        "checksums": int(values[7]),
        "hash_ready": int(values[8]),
        "baselines": int(values[9]),
    }


def draw_dot_grid(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    total: int,
    filled: int,
    color: str,
    empty: str,
    *,
    cols: int = 11,
    radius: int = 13,
) -> None:
    gap_x, gap_y = 49, 48
    for i in range(total):
        col, row = i % cols, i // cols
        cx, cy = x + col * gap_x, y + row * gap_y
        if i < filled:
            draw.ellipse((cx - radius, cy - radius, cx + radius, cy + radius), fill=color, outline="#EAF2FA", width=2)
        else:
            draw.ellipse((cx - radius, cy - radius, cx + radius, cy + radius), fill="#101823", outline=empty, width=3)


def draw_step_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    number: int,
    label: str,
    title: str,
    count: str,
    body: str,
    color: str,
    state: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 28, COLORS["panel2"], color, 2)
    draw.ellipse((x1 + 26, y1 + 24, x1 + 76, y1 + 74), fill=color, outline="#EAF2FA", width=2)
    text(draw, (x1 + 51, y1 + 35), str(number), F["small_bold"], COLORS["ink"], "ma")
    text(draw, (x1 + 96, y1 + 26), label.upper(), F["micro_bold"], color)
    text(draw, (x2 - 28, y1 + 26), state.upper(), F["micro_bold"], COLORS["muted"], "ra")
    text(draw, (x1 + 28, y1 + 92), title, F["h3"], COLORS["text"])
    text(draw, (x1 + 28, y1 + 140), count, F["big"], COLORS["text"])
    paragraph(draw, (x1 + 30, y1 + 226), body, F["small"], x2 - x1 - 60, COLORS["muted"], 6)


def draw_ladder(draw: ImageDraw.ImageDraw, data: dict[str, int]) -> None:
    x1, y1, x2, y2 = 120, 610, 3720, 1316
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["teal"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Five-Step Readiness Ladder", F["section"], COLORS["text"])
    text(draw, (x2 - 80, y1 + 62), "current position: metadata alpha", F["small_bold"], COLORS["teal"], "ra")

    steps = [
        ("metadata", "Metadata indexed", f"{data['tasks']} tasks", f"{data['folds']} folds and {data['study_rows']} study rows assembled", COLORS["blue"], "ready"),
        ("checksums", "Checksums parsed", f"{data['checksums']}/{data['study_rows']}", "checksum manifests parsed into audit tables", COLORS["green"], "ready"),
        ("files", "Files mirrored", "next", "local study files copied into a fixed folder", COLORS["amber"], "next"),
        ("hashes", "Hashes verified", f"{data['hash_ready']}/{data['study_rows']}", "local file hashes matched to checksum manifests", COLORS["teal"], "check"),
        ("package", "Release package", "bundle", "verified files travel with task objects and reports", COLORS["violet"], "assemble"),
    ]
    card_w, card_h, gap = 624, 384, 78
    start_x, card_y = x1 + 56, y1 + 180
    for i, (label, title, count, body, color, state) in enumerate(steps):
        cx = start_x + i * (card_w + gap)
        draw_step_card(draw, (cx, card_y, cx + card_w, card_y + card_h), i + 1, label, title, count, body, color, state)
        if i < len(steps) - 1:
            arrow(draw, cx + card_w + 18, card_y + 192, cx + card_w + gap - 24, card_y + 192, COLORS["dim"], 4)

    track_y = y2 - 88
    rounded(draw, (x1 + 58, track_y, x2 - 58, track_y + 36), 18, "#0F1825", "#2A394D", 1)
    segment_w = (x2 - x1 - 116) // 5
    track_colors = [COLORS["blue"], COLORS["green"], COLORS["amber"], COLORS["teal"], COLORS["violet"]]
    for i, color in enumerate(track_colors):
        sx = x1 + 58 + i * segment_w
        alpha_fill = color if i < 2 else "#172335"
        rounded(draw, (sx + 6, track_y + 6, sx + segment_w - 8, track_y + 30), 12, alpha_fill, color, 1)
    text(draw, (x1 + 76, track_y + 50), "ready now: indexed metadata + parsed checksums", F["tiny_bold"], COLORS["green"])
    text(draw, (x2 - 76, track_y + 50), "next: mirror files, verify hashes, assemble package", F["tiny_bold"], COLORS["amber"], "ra")


def draw_row_checks(draw: ImageDraw.ImageDraw, data: dict[str, int]) -> None:
    x1, y1, x2, y2 = 120, 1380, 2380, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["green"], width=5)
    text(draw, (x1 + 42, y1 + 52), "22 Study-Row Checks", F["h2"], COLORS["text"])

    rows = [
        ("checksum manifests parsed", data["checksums"], COLORS["green"], "metadata layer complete"),
        ("local hashes verified", data["hash_ready"], COLORS["amber"], "file layer to complete"),
    ]
    for i, (label, count, color, note) in enumerate(rows):
        top = y1 + 134 + i * 160
        rounded(draw, (x1 + 54, top, x2 - 54, top + 118), 22, COLORS["panel2"], color, 2)
        text(draw, (x1 + 84, top + 22), label.upper(), F["tiny_bold"], color)
        text(draw, (x2 - 84, top + 20), f"{count}/{data['study_rows']}", F["h2"], COLORS["text"], "ra")
        text(draw, (x1 + 84, top + 58), note, F["small"], COLORS["muted"])
        draw_dot_grid(draw, x1 + 700, top + 54, data["study_rows"], count, color, COLORS["dim"], radius=12)


def draw_release_requirements(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 2460, 1380, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=COLORS["amber"], width=5)
    text(draw, (x1 + 42, y1 + 52), "Release-Ready Bundle", F["h2"], COLORS["text"])

    items = [
        ("task objects", "manifests, metric profiles, evaluators", COLORS["blue"]),
        ("verified files", "local files with matching hashes", COLORS["teal"]),
        ("reports", "baseline outputs and run records", COLORS["violet"]),
    ]
    top = y1 + 126
    for i, (title, body, color) in enumerate(items):
        y = top + i * 84
        rounded(draw, (x1 + 54, y, x2 - 54, y + 62), 18, COLORS["panel2"], color, 2)
        text(draw, (x1 + 82, y + 14), title, F["small_bold"], COLORS["text"])
        text(draw, (x1 + 396, y + 14), body, F["small"], COLORS["muted"])

    rounded(draw, (x1 + 54, y2 - 68, x2 - 54, y2 - 30), 14, "#122234", COLORS["green"], 2)
    text(draw, (x1 + 80, y2 - 56), "Ready after file checks pass.", F["micro_bold"], COLORS["green"])


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "The alpha is already inspectable, but a release-ready data bundle needs local files, matching hashes, and reports packaged together.",
        F["small"],
        3180,
        COLORS["muted"],
        8,
    )


def build() -> None:
    data = load_summary()
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 54 | ACT 6 | FILE READINESS LADDER", F["kicker"], COLORS["teal"])
    bx = 1840
    bx = badge(draw, bx, 56, "TASKS", str(data["tasks"]), COLORS["blue"])
    bx = badge(draw, bx, 56, "STUDY ROWS", str(data["study_rows"]), COLORS["green"])
    bx = badge(draw, bx, 56, "CHECKSUMS", f"{data['checksums']}/{data['study_rows']}", COLORS["teal"])
    badge(draw, bx, 56, "HASH VERIFIED", f"{data['hash_ready']}/{data['study_rows']}", COLORS["amber"])

    text(draw, (120, 330), "File Readiness Ladder Turns Metadata Into A Release Package", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "Five checks separate indexed metadata from a release-ready data bundle.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_ladder(draw, data)
    draw_row_checks(draw, data)
    draw_release_requirements(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    ImageOps.colorize(gray, black="#080D14", white="#F3F7FC").save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(gray)

    manifest = {
        "title": "File Readiness Ladder Turns Metadata Into A Release Package",
        "summary": "Five checks separate indexed metadata from a release-ready data bundle.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "data_inputs": {"snapshot_decision_summary": str(SUMMARY_CSV.relative_to(ROOT))},
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": data,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    QA_NOTE.write_text(
        "\n".join(
            [
                "# File Readiness Ladder Asset Visual QA",
                "",
                "Slide 54 explains the five checks between metadata alpha and release-ready data bundle.",
                "",
                "Checks to review:",
                "- Five ladder cards read left to right without arrow overlap.",
                "- 22/22 checksum parsing and 0/22 hash verification are visible.",
                "- Dot grids fit inside their cards and remain readable in grayscale.",
                "- Release-ready bundle panel has no text collision.",
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

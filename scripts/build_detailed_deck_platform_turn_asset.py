#!/usr/bin/env python3
"""Build slide 50 asset: platform turn packages results into a reproducible benchmark."""

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
OUT_DIR = ASSET_ROOT / "platform_turn"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "platform_turn_packages_results_into_reproducible_benchmark_premium.png"
GRAY_PATH = OUT_DIR / "platform_turn_packages_results_into_reproducible_benchmark_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "platform_turn_manifest.json"
QA_NOTE = OUT_DIR / "PLATFORM_TURN_ASSET_VISUAL_QA.md"

SNAPSHOT = ROOT / "v9" / "reports" / "public_bulk_alpha_snapshot_decision" / "snapshot_decision_summary.json"

COLORS = {
    "bg": "#0B1119",
    "bg2": "#111721",
    "header": "#101826",
    "footer": "#090E15",
    "panel": "#111B28",
    "panel2": "#172233",
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
    "title": load_font(82, True),
    "subtitle": load_font(37),
    "section": load_font(41, True),
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
    "big": load_font(70, True),
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
    width = max(224, int(draw.textlength(value, font=F["tiny_bold"]) + 78))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def arrow(draw: ImageDraw.ImageDraw, x1: int, y1: int, x2: int, y2: int, color: str, width: int = 5) -> None:
    draw.line((x1, y1, x2, y2), fill=color, width=width)
    draw.polygon([(x2, y2), (x2 - 24, y2 - 14), (x2 - 24, y2 + 14)], fill=color)


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


def load_public_bulk_counts() -> dict[str, str]:
    row = json.loads(SNAPSHOT.read_text(encoding="utf-8"))[0]
    return {
        "tasks": row["bulk_task_count"],
        "folds": row["bulk_fold_count"],
        "baselines": row["baseline_row_count"],
        "bulk_rows": row["bulk_source_count"],
    }


def draw_evidence_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 610, 1345, 1462
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Evidence Workstreams", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 94),
        "The result section has already produced reusable analysis outputs.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )

    rows = [
        ("score readouts", "mission-held-out task performance", COLORS["blue"]),
        ("model comparisons", "classical ML, FMs, text LLM checks", COLORS["violet"]),
        ("biology layers", "systems, immune, target, translation", COLORS["green"]),
        ("v8 workstreams", "BRIDGE, DECOMPOSE, INTERVENE, CAUSAL", COLORS["teal"]),
    ]
    yy = y1 + 216
    for title, desc, color in rows:
        rounded(draw, (x1 + 58, yy, x2 - 58, yy + 112), 24, COLORS["panel2"], color, 2)
        rounded(draw, (x1 + 84, yy + 28, x1 + 136, yy + 80), 17, "#122234", color, 2)
        text(draw, (x1 + 162, yy + 22), title, F["small_bold"], COLORS["text"])
        text(draw, (x1 + 162, yy + 63), desc, F["tiny"], COLORS["muted"])
        yy += 140


def draw_platform_turn(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 1415, 610, 2375, 1462
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 48, y1 + 38), "Platform Turn", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 50, y1 + 96),
        "The same evidence becomes a system of objects that can be rerun, compared, and packaged.",
        F["small"],
        x2 - x1 - 100,
        COLORS["muted"],
        7,
    )

    cy = y1 + 472
    rounded(draw, (x1 + 170, cy - 104, x2 - 170, cy + 104), 42, "#122234", COLORS["teal"], 3)
    text(draw, ((x1 + x2) // 2, cy - 30), "RESULT", F["metric"], COLORS["text"], "mm")
    text(draw, ((x1 + x2) // 2, cy + 32), "TO SYSTEM", F["metric"], COLORS["teal"], "mm")
    arrow(draw, x1 + 96, cy, x1 + 146, cy, COLORS["dim"], 5)
    arrow(draw, x2 - 146, cy, x2 - 96, cy, COLORS["dim"], 5)

    chips = [
        ("comparable", COLORS["blue"]),
        ("repeatable", COLORS["green"]),
        ("reviewable", COLORS["amber"]),
    ]
    cx = x1 + 112
    for label, color in chips:
        rounded(draw, (cx, y2 - 142, cx + 224, y2 - 82), 18, COLORS["panel2"], color, 2)
        text(draw, (cx + 112, y2 - 112), label, F["tiny_bold"], COLORS["text"], "mm")
        cx += 254


def draw_system_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 2445, 610, 3720, 1462
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Benchmark System", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 94),
        "The benchmark becomes runnable through a compact object chain.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )

    layers = [
        ("task registry", "what is being predicted", COLORS["blue"]),
        ("metric profiles", "how each task is scored", COLORS["teal"]),
        ("evaluator", "how predictions become scores", COLORS["green"]),
        ("run record", "what was executed", COLORS["violet"]),
        ("release package", "what gets shared and reviewed", COLORS["amber"]),
    ]
    yy = y1 + 216
    for i, (title, desc, color) in enumerate(layers):
        inset = i * 18
        rounded(draw, (x1 + 72 + inset, yy, x2 - 72 - inset, yy + 84), 22, COLORS["panel2"], color, 2)
        text(draw, (x1 + 104 + inset, yy + 16), title, F["small_bold"], COLORS["text"])
        text(draw, (x1 + 104 + inset, yy + 50), desc, F["tiny"], COLORS["muted"])
        if i < len(layers) - 1:
            arrow(draw, (x1 + x2) // 2, yy + 96, (x1 + x2) // 2, yy + 122, COLORS["dim"], 3)
        yy += 120


def draw_object_chain(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1516, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Act 6 Object Chain", F["h2"], COLORS["text"])

    steps = [
        ("task", "typed prediction unit", COLORS["blue"]),
        ("metric", "shared scoring rule", COLORS["teal"]),
        ("evaluator", "score generator", COLORS["green"]),
        ("run record", "execution trace", COLORS["violet"]),
        ("package", "release artifact", COLORS["amber"]),
    ]
    node_w, gap = 575, 70
    start_x, node_y = x1 + 190, y1 + 150
    for i, (title, desc, color) in enumerate(steps):
        nx = start_x + i * (node_w + gap)
        rounded(draw, (nx, node_y, nx + node_w, node_y + 96), 20, COLORS["panel2"], color, 2)
        text(draw, (nx + 26, node_y + 17), title, F["small_bold"], COLORS["text"])
        text(draw, (nx + 26, node_y + 55), desc, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            arrow(draw, nx + node_w + 10, node_y + 48, nx + node_w + gap - 16, node_y + 48, COLORS["dim"], 4)


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "The platform layer turns analysis into a benchmark system: tasks define comparisons, metrics define scoring, evaluator runs produce records, and packages carry the artifacts.",
        F["small"],
        3180,
        COLORS["muted"],
        8,
    )


def build() -> None:
    counts = load_public_bulk_counts()
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 50 | ACT 6 | PLATFORM TURN", F["kicker"], COLORS["teal"])
    bx = 2050
    bx = badge(draw, bx, 56, "TASKS", counts["tasks"], COLORS["blue"])
    bx = badge(draw, bx, 56, "FOLDS", counts["folds"], COLORS["teal"])
    bx = badge(draw, bx, 56, "BASELINES", counts["baselines"], COLORS["green"])
    badge(draw, bx, 56, "BULK ROWS", counts["bulk_rows"], COLORS["amber"])

    text(draw, (120, 330), "Platform Turn Packages Results Into A Reproducible Benchmark", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "Evidence workstreams become an object chain that makes runs comparable, repeatable, and reviewable.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_evidence_panel(draw)
    draw_platform_turn(draw)
    draw_system_panel(draw)
    draw_object_chain(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    ImageOps.colorize(gray, black="#080D14", white="#F3F7FC").save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(gray)

    manifest = {
        "title": "Platform Turn Packages Results Into A Reproducible Benchmark",
        "claim": "Act 6 turns evidence workstreams into a reproducible benchmark object chain.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "count_input": str(SNAPSHOT.relative_to(ROOT)),
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": counts,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    QA_NOTE.write_text(
        "\n".join(
            [
                "# Platform Turn Asset Visual QA",
                "",
                "Slide 50 transitions from Act 5 evidence workstreams to the Act 6 benchmark object chain.",
                "",
                "Checks to review:",
                "- Three main panels separate evidence, platform turn, and benchmark system.",
                "- Arrows stay in gutters and do not touch text.",
                "- Header badges remain readable at presentation scale.",
                "- Bottom object chain previews Act 6 without crowding.",
                "- Grayscale version preserves the panel hierarchy.",
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

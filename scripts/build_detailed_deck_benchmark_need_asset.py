#!/usr/bin/env python3
"""Build slide 2 asset: why the project needs a benchmark layer."""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont, ImageOps, ImageStat


ROOT = Path(__file__).resolve().parent.parent
ASSET_ROOT = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
)
OUT_DIR = ASSET_ROOT / "benchmark_need"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "why_this_needs_a_benchmark_premium.png"
GRAY_PATH = OUT_DIR / "why_this_needs_a_benchmark_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "why_this_needs_a_benchmark_manifest.json"

EVIDENCE_STACK = ASSET_ROOT / "evidence_stack" / "evidence_stack_manifest.json"

COLORS = {
    "bg": "#0C111A",
    "bg2": "#091019",
    "panel": "#101823",
    "panel2": "#151F2D",
    "panel3": "#1A2534",
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
    "red": "#F06E7D",
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
    "title": load_font(86, True),
    "subtitle": load_font(38),
    "h1": load_font(52, True),
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
    "stat": load_font(76, True),
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


def load_metrics() -> dict[str, object]:
    evidence = json.loads(EVIDENCE_STACK.read_text())["metrics"]
    return {
        "osd_studies": 24,
        "missions": 17,
        "samples": "~450",
        "tissues": evidence["mission_tasks"],
        "lomo_folds": evidence["lomo_folds"],
        "mission_labels": evidence["mission_labels"],
    }


def draw_background(draw: ImageDraw.ImageDraw) -> None:
    draw.rectangle((0, 0, W, H), fill=COLORS["bg"])
    for y in range(H):
        t = y / H
        draw.line((0, y, W, y), fill=blend(COLORS["bg"], COLORS["bg2"], t))
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
    text(draw, (150, 76), "OPENING / WHY A BENCHMARK", F["kicker"], COLORS["teal"])
    x = 2260
    x = badge(draw, x, 66, "study records", "NASA OSDR", COLORS["blue"])
    x = badge(draw, x, 66, "curated set", f"{data['osd_studies']} OSD studies", COLORS["teal"])
    badge(draw, x, 66, "folds", f"{data['lomo_folds']} LOMO", COLORS["amber"])


def draw_title(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 382), "Why This Needs A Benchmark", F["title"], COLORS["text"])
    paragraph(
        draw,
        (155, 490),
        "Public study records explain what was measured; benchmark folds test which signals transfer to a new mission.",
        F["subtitle"],
        1580,
        COLORS["muted"],
        10,
    )


def arrow(draw: ImageDraw.ImageDraw, x1: int, y1: int, x2: int, y2: int, color: str) -> None:
    draw.line((x1, y1, x2, y2), fill=color, width=5)
    draw.polygon([(x2, y2), (x2 - 22, y2 - 13), (x2 - 22, y2 + 13)], fill=color)


def mini_card(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], label: str, body: str, color: str) -> None:
    rounded(draw, box, 18, blend("#101823", color, 0.10), color, 2)
    x1, y1, x2, _ = box
    text(draw, (x1 + 24, y1 + 22), label, F["tiny_bold"], color)
    paragraph(draw, (x1 + 24, y1 + 56), body, F["tiny"], x2 - x1 - 48, COLORS["muted"], 4)


def draw_study_record_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    panel = (150, 680, 1160, 1620)
    rounded(draw, panel, 30, COLORS["panel"], "#2A394D", 2)
    text(draw, (200, 730), "A. Public study records", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (200, 785),
        "OSDR gives the study objects: studies, samples, assays, metadata, and processed files.",
        F["small"],
        860,
        COLORS["muted"],
        7,
    )
    stats = [
        (data["osd_studies"], "OSD studies", COLORS["blue"]),
        (data["missions"], "missions", COLORS["teal"]),
        (data["samples"], "binary samples", COLORS["amber"]),
    ]
    x = 205
    for value, label, color in stats:
        rounded(draw, (x, 930, x + 285, 1125), 22, "#142033", color, 2)
        text(draw, (x + 34, 962), str(value), F["stat"], COLORS["text"])
        text(draw, (x + 34, 1054), label, F["tiny_bold"], COLORS["muted"])
        x += 305
    stack = [
        ("Study metadata", "organism, tissue, mission context"),
        ("Sample table", "Flight / Ground labels and group fields"),
        ("Assay files", "bulk expression matrices and derived products"),
        ("DGE records", "per-study differential-expression outputs"),
    ]
    y = 1190
    for idx, (label, body) in enumerate(stack):
        color = [COLORS["blue"], COLORS["teal"], COLORS["violet"], COLORS["amber"]][idx]
        mini_card(draw, (205, y, 1095, y + 82), label, body, color)
        y += 92


def draw_gap_bridge(draw: ImageDraw.ImageDraw) -> None:
    panel = (1240, 680, 2440, 1620)
    rounded(draw, panel, 30, COLORS["panel"], "#2A394D", 2)
    text(draw, (1290, 730), "B. Evidence layer added by benchmarking", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (1290, 785),
        "The benchmark freezes a task contract, then asks whether a signal moves from known missions to a hidden mission.",
        F["small"],
        990,
        COLORS["muted"],
        7,
    )
    y = 940
    steps = [
        ("Study-level readout", "What changed inside one mission or study?", COLORS["blue"]),
        ("Task contract", "Input, label, split, metric, and readout frame are fixed.", COLORS["teal"]),
        ("Hidden-mission fold", "Known missions train; one whole mission waits outside fitting.", COLORS["amber"]),
        ("Transfer score", "AUROC plus uncertainty turns a model result into benchmark evidence.", COLORS["green"]),
    ]
    for idx, (label, body, color) in enumerate(steps):
        box = (1310, y + idx * 145, 2370, y + idx * 145 + 105)
        rounded(draw, box, 20, blend("#101823", color, 0.10), color, 2)
        text(draw, (box[0] + 28, box[1] + 22), label, F["h3"], COLORS["text"])
        paragraph(draw, (box[0] + 28, box[1] + 62), body, F["tiny"], 970, COLORS["muted"], 3)
        if idx < len(steps) - 1:
            next_top = y + (idx + 1) * 145
            tip_y = next_top - 8
            base_y = tip_y - 16
            draw.line((1840, box[3] + 6, 1840, base_y), fill=color, width=4)
            draw.polygon([(1840, tip_y), (1827, base_y), (1853, base_y)], fill=color)
    rounded(draw, (1310, 1510, 2370, 1572), 16, "#111B28", "#26384E", 2)
    text(draw, (1335, 1528), "Takeaway", F["tiny_bold"], COLORS["amber"])
    text(draw, (1505, 1528), "Study evidence becomes transfer evidence only after the split is explicit.", F["tiny"], COLORS["muted"])


def draw_hierarchy_map(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    panel = (2520, 680, 3690, 1620)
    rounded(draw, panel, 30, COLORS["panel"], "#2A394D", 2)
    text(draw, (2570, 730), "C. Hierarchy of questions", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (2570, 785),
        "Each level keeps the study record intact while changing the question being answered.",
        F["small"],
        960,
        COLORS["muted"],
        7,
    )

    levels = [
        ("Record", "What is available?", "OSDR studies, metadata, samples", COLORS["blue"], 0.26),
        ("Study", "What changed here?", "mission/tissue contrasts and DGE", COLORS["teal"], 0.47),
        ("Task", "What is scored?", "fixed input, label, split, metric", COLORS["violet"], 0.68),
        ("Benchmark", "What transfers?", f"{data['tissues']} tissues, {data['lomo_folds']} hidden-mission folds", COLORS["amber"], 0.89),
    ]
    axis_x = 2625
    y_top, y_bottom = 960, 1485
    draw.line((axis_x, y_top, axis_x, y_bottom), fill=COLORS["dim"], width=4)
    draw.polygon([(axis_x, y_top - 18), (axis_x - 13, y_top + 10), (axis_x + 13, y_top + 10)], fill=COLORS["dim"])
    text(draw, (2655, y_top - 30), "evaluation depth", F["micro_bold"], COLORS["dim"])

    for label, question, body, color, t in levels:
        cy = int(y_bottom - (y_bottom - y_top) * t)
        draw.ellipse((axis_x - 18, cy - 18, axis_x + 18, cy + 18), fill=color, outline=COLORS["text"], width=2)
        rounded(draw, (2695, cy - 62, 3615, cy + 62), 18, blend("#101823", color, 0.10), color, 2)
        text(draw, (2724, cy - 43), label, F["tiny_bold"], color)
        text(draw, (2860, cy - 43), question, F["small_bold"], COLORS["text"])
        paragraph(draw, (2724, cy - 5), body, F["tiny"], 820, COLORS["muted"], 3)

    rounded(draw, (2570, 1510, 3615, 1572), 16, blend("#101823", COLORS["green"], 0.10), COLORS["green"], 2)
    text(draw, (2595, 1528), "Positioning", F["tiny_bold"], COLORS["green"])
    text(draw, (2765, 1528), "Benchmarking is an evaluation layer on top of OSDR records.", F["tiny"], COLORS["muted"])


def draw_flow_line(draw: ImageDraw.ImageDraw) -> None:
    y = 1660
    x1, x2 = 320, 3510
    draw.line((x1, y, x2, y), fill="#26384E", width=4)
    draw.polygon([(x2, y), (x2 - 24, y - 14), (x2 - 24, y + 14)], fill="#26384E")
    labels = [
        (520, "record", COLORS["blue"]),
        (1035, "study evidence", COLORS["teal"]),
        (1680, "task contract", COLORS["violet"]),
        (2320, "hidden mission", COLORS["amber"]),
        (3040, "transfer readout", COLORS["green"]),
    ]
    for x, label, color in labels:
        draw.ellipse((x - 12, y - 12, x + 12, y + 12), fill=color)
        text(draw, (x - 54, y + 28), label, F["tiny_bold"], color)


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (150, 1886, 3690, 2068), 22, "#111B28", "#26384E", 2)
    text(draw, (190, 1928), "Takeaway", F["tiny_bold"], COLORS["blue"])
    paragraph(
        draw,
        (360, 1928),
        "Public study records answer what happened; benchmark tasks ask whether a signal transfers to an unseen mission.",
        F["small"],
        2710,
        COLORS["muted"],
        7,
    )
    text(draw, (190, 1992), "Next", F["tiny_bold"], COLORS["amber"])
    paragraph(
        draw,
        (360, 1992),
        "The next method slides define the fixed task record, sample matrix, hidden mission split, leakage guard, and metric.",
        F["small"],
        2850,
        COLORS["muted"],
        7,
    )
    text(draw, (3510, 1998), "2", F["h2"], COLORS["teal"])


def render() -> Image.Image:
    data = load_metrics()
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image)
    draw_background(draw)
    draw_header(draw, data)
    draw_title(draw)
    draw_study_record_panel(draw, data)
    draw_gap_bridge(draw)
    draw_hierarchy_map(draw, data)
    draw_flow_line(draw)
    draw_footer(draw)
    return image


def image_metrics(path: Path) -> dict[str, object]:
    image = Image.open(path).convert("RGB")
    gray = ImageOps.grayscale(image)
    stat = ImageStat.Stat(gray)
    edges = gray.filter(ImageFilter.FIND_EDGES)
    edge_stat = ImageStat.Stat(edges)
    return {
        "width_px": image.width,
        "height_px": image.height,
        "luma_mean": round(stat.mean[0], 2),
        "luma_stddev": round(stat.stddev[0], 2),
        "edge_mean": round(edge_stat.mean[0], 2),
        "nonblank_pass": stat.stddev[0] >= 16.0 and edge_stat.mean[0] >= 2.0,
    }


def main() -> None:
    image = render()
    image.save(OUT_PATH, quality=95)
    ImageOps.grayscale(image).convert("RGB").save(GRAY_PATH, quality=95)
    data = load_metrics()
    manifest = {
        "title": "Why This Needs A Benchmark",
        "slide": 2,
        "status": "ready",
        "role": "OSDR/study/mission benchmark-need map",
        "created": "2026-06-11",
        "outputs": {
            "png": str(OUT_PATH.relative_to(ROOT)),
            "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        },
        "metrics": data,
        "copy": {
            "headline": "Why This Needs A Benchmark",
            "subtitle": "Public study records explain what was measured; benchmark folds test which signals transfer to a new mission.",
            "reader_rule": "Study evidence becomes transfer evidence only after the split is explicit.",
            "scope": "Benchmarking is an evaluation layer on top of OSDR records.",
        },
        "source_notes": [
            "docs/V1_PAPER_CONTENT.md: 24 OSD studies, 17 missions, approximately 450 binary samples.",
            "evidence_stack_manifest.json: 6 tissue tasks, 22 LOMO folds, 10 mission labels.",
        ],
        "automatic_metrics": image_metrics(OUT_PATH),
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"png": str(OUT_PATH), "grayscale": str(GRAY_PATH), "manifest": str(MANIFEST_PATH)}, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build slide 46 asset: stressor interactions dominate combined effects."""

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
OUT_DIR = ASSET_ROOT / "decompose_interactions"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "stressor_interactions_dominate_combined_effects_premium.png"
GRAY_PATH = OUT_DIR / "stressor_interactions_dominate_combined_effects_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "decompose_interactions_manifest.json"
QA_NOTE = OUT_DIR / "DECOMPOSE_INTERACTIONS_ASSET_VISUAL_QA.md"

V8_PILLARS_MANIFEST = ASSET_ROOT / "v8_pillars" / "v8_pillars_manifest.json"
V8_FIGURE = ROOT / "v8" / "figures" / "Figure2_Stressor_Decomposition.png"
FACTORIAL = ROOT / "v8" / "decompose" / "evaluation" / "factorial_flight_decomposition.json"
ICP = ROOT / "v8" / "causal" / "evaluation" / "icp_dag_summary.json"

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
    "orange": "#E6A15D",
    "violet": "#B39DFF",
    "rose": "#E17882",
    "red": "#F17C88",
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
    "title": load_font(84, True),
    "subtitle": load_font(37),
    "section": load_font(41, True),
    "h2": load_font(36, True),
    "h3": load_font(30, True),
    "body": load_font(27),
    "body_bold": load_font(27, True),
    "small": load_font(24),
    "small_bold": load_font(24, True),
    "tiny": load_font(21),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "metric": load_font(70, True),
    "metric2": load_font(54, True),
    "formula": load_font(31, True),
    "axis": load_font(18, True),
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
    leading: int = 8,
) -> int:
    x, y = xy
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(224, int(draw.textlength(value, font=F["tiny_bold"]) + 76))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def arrow(draw: ImageDraw.ImageDraw, x1: int, y1: int, x2: int, y2: int, color: str, width: int = 4) -> None:
    draw.line((x1, y1, x2, y2), fill=color, width=width)
    pts = [
        (x2, y2),
        (x2 - 20, y2 - 11),
        (x2 - 20, y2 + 11),
    ]
    draw.polygon(pts, fill=color)


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def pct(value: float) -> int:
    return int(round(value * 100))


def load_data() -> dict[str, object]:
    factorial = load_json(FACTORIAL)
    icp = load_json(ICP)
    spleen = factorial["spleen"]["per_flight_tissue"]["thymus"]["variance_attribution_top200"]
    skin = factorial["skin_analog"]["per_flight_tissue"]["skin"]["variance_attribution_top200"]
    hze = factorial["hze_endocrine"]["per_flight_tissue"]["thymus"]["variance_attribution_top200"]
    brain_sig = factorial["brain"]["n_sig_p05"]
    icp_rows = [
        (name, float(stats["mean"]), int(stats["count"]))
        for name, stats in icp["stressor_aggregate_icp"].items()
    ]
    icp_rows.sort(key=lambda row: row[1], reverse=True)
    return {
        "surfaces": [
            (
                "Spleen -> thymus",
                [("HLU", spleen["HLU_frac"], COLORS["blue"]), ("IR", spleen["IR_frac"], COLORS["amber"]), ("HLUxIR", spleen["HLUxIR_frac"], COLORS["rose"])],
                "mouse analog",
            ),
            (
                "Skin -> skin",
                [("HLU", skin["HLU_frac"], COLORS["blue"]), ("IR", skin["IR_frac"], COLORS["amber"]), ("HLUxIR", skin["HLUxIR_frac"], COLORS["rose"])],
                "mouse analog",
            ),
            (
                "HZE -> thymus",
                [("Sex", hze["Sex_frac"], COLORS["teal"]), ("HZE", hze["HZE_frac"], COLORS["amber"]), ("SexxHZE", hze["SexxHZE_frac"], COLORS["rose"])],
                "endocrine analog",
            ),
        ],
        "brain_sig": brain_sig,
        "icp": icp_rows,
    }


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


def draw_factorial_tile(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], label: str, fill: str, sub: str) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 20, fill, "#30425B", 2)
    text(draw, ((x1 + x2) / 2, y1 + 32), label, F["tiny_bold"], COLORS["text"], "mm")
    text(draw, ((x1 + x2) / 2, y1 + 70), sub, F["micro_bold"], COLORS["muted"], "mm")


def draw_method_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 610, 1210, 1462
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 36), "What DECOMPOSE does", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 92),
        "Factorial analog studies split one expression response into main stressor axes and a combined-response surface.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )

    gx, gy = x1 + 76, y1 + 228
    tile_w, tile_h = 175, 106
    draw_factorial_tile(draw, (gx, gy, gx + tile_w, gy + tile_h), "0", COLORS["panel3"], "reference")
    draw_factorial_tile(draw, (gx + 220, gy, gx + 220 + tile_w, gy + tile_h), "HLU", "#19324C", "axis A")
    draw_factorial_tile(draw, (gx, gy + 148, gx + tile_w, gy + 148 + tile_h), "IR", "#40301B", "axis B")
    draw_factorial_tile(draw, (gx + 220, gy + 148, gx + 220 + tile_w, gy + 148 + tile_h), "HLU+IR", "#3D2431", "combined")
    draw.line((gx + 86, gy + 124, gx + 86, gy + 142), fill=COLORS["dim"], width=3)
    draw.line((gx + 196, gy + 53, gx + 216, gy + 53), fill=COLORS["dim"], width=3)
    draw.line((gx + 196, gy + 201, gx + 216, gy + 201), fill=COLORS["dim"], width=3)
    draw.line((gx + 306, gy + 124, gx + 306, gy + 142), fill=COLORS["dim"], width=3)
    arrow(draw, gx + 440, gy + 128, x2 - 428, gy + 128, COLORS["dim"], 4)

    rx, ry = x2 - 390, y1 + 232
    lanes = [
        ("main HLU", COLORS["blue"]),
        ("main radiation", COLORS["amber"]),
        ("time or sex", COLORS["teal"]),
        ("interaction surface", COLORS["rose"]),
    ]
    for i, (label, color) in enumerate(lanes):
        y = ry + i * 82
        rounded(draw, (rx, y, x2 - 74, y + 58), 17, COLORS["panel2"], color, 2)
        text(draw, (rx + 24, y + 17), label, F["tiny_bold"], COLORS["text"])

    rounded(draw, (x1 + 64, y1 + 590, x2 - 64, y1 + 704), 24, COLORS["panel2"], COLORS["violet"], 2)
    text(draw, (x1 + 96, y1 + 624), "expression response", F["formula"], COLORS["text"])
    text(draw, (x1 + 438, y1 + 624), "=", F["formula"], COLORS["dim"])
    text(draw, (x1 + 484, y1 + 624), "main axes + interaction surface", F["formula"], COLORS["muted"])

    rounded(draw, (x1 + 64, y2 - 128, x2 - 64, y2 - 42), 24, "#122234", COLORS["teal"], 2)
    paragraph(
        draw,
        (x1 + 96, y2 - 104),
        "Takeaway: the interaction term measures the extra shape created by the combined condition.",
        F["small"],
        x2 - x1 - 192,
        COLORS["muted"],
        6,
    )


def draw_stacked_bar(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    width: int,
    height: int,
    parts: list[tuple[str, float, str]],
) -> None:
    rounded(draw, (x, y, x + width, y + height), 16, COLORS["panel3"], "#2A394D", 1)
    cursor = x
    for i, (label, value, color) in enumerate(parts):
        seg_w = int(round(width * value))
        if i == len(parts) - 1:
            seg_w = x + width - cursor
        box = (cursor, y, cursor + seg_w, y + height)
        rounded(draw, box, 14, color)
        fill = COLORS["ink"] if color in {COLORS["amber"], COLORS["teal"], COLORS["green"]} else COLORS["text"]
        segment_label = f"{label} {pct(value)}%"
        if seg_w >= 150:
            text(draw, (cursor + seg_w / 2, y + height / 2), segment_label, F["tiny_bold"], fill, "mm")
        elif seg_w >= 104:
            text(draw, (cursor + seg_w / 2, y + height / 2), f"{pct(value)}%", F["tiny_bold"], fill, "mm")
        else:
            text(draw, (cursor + seg_w / 2, y - 14), f"{pct(value)}%", F["micro_bold"], color, "mm")
        cursor += seg_w


def draw_surface_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x1, y1, x2, y2 = 1265, 610, 2575, 1462
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 36), "Interaction share is large", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 92),
        "Top responsive genes often assign a large variance share to the combined-stressor term.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )
    chips = [("spleen", COLORS["blue"]), ("skin", COLORS["teal"]), ("brain/time", COLORS["violet"]), ("HZE endocrine", COLORS["amber"])]
    cx, cy = x1 + 44, y1 + 180
    for label, color in chips:
        w = int(draw.textlength(label, font=F["micro_bold"]) + 42)
        rounded(draw, (cx, cy, cx + w, cy + 40), 13, COLORS["panel2"], color, 2)
        text(draw, (cx + 21, cy + 12), label, F["micro_bold"], COLORS["text"])
        cx += w + 14

    bar_x = x1 + 380
    bar_w = 665
    base_y = y1 + 290
    rows = data["surfaces"]
    for i, (label, parts, note) in enumerate(rows):
        y = base_y + i * 148
        text(draw, (x1 + 68, y + 2), label, F["small_bold"], COLORS["text"])
        text(draw, (x1 + 68, y + 38), note, F["micro_bold"], COLORS["dim"])
        draw_stacked_bar(draw, bar_x, y, bar_w, 64, parts)
        top_part = max(parts, key=lambda part: part[1])
        rounded(draw, (bar_x + bar_w + 24, y - 3, x2 - 62, y + 67), 18, COLORS["panel2"], top_part[2], 2)
        text(draw, (bar_x + bar_w + 52, y + 15), "largest", F["micro_bold"], COLORS["dim"])
        text(draw, (bar_x + bar_w + 52, y + 42), top_part[0], F["tiny_bold"], COLORS["text"])

    rounded(draw, (x1 + 64, y2 - 128, x2 - 64, y2 - 42), 24, "#122234", COLORS["rose"], 2)
    paragraph(
        draw,
        (x1 + 96, y2 - 104),
        "Readout: combined-stressor terms explain 40-62% of selected response surfaces.",
        F["small"],
        x2 - x1 - 192,
        COLORS["muted"],
        6,
    )


def draw_icp_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x1, y1, x2, y2 = 2630, 610, 3720, 1462
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 36), "Causal-stability readout", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 92),
        "ICP ranks stressor axes by how stable their gene associations stay across environments.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )

    chart_x, chart_y = x1 + 255, y1 + 206
    bar_w, bar_h = 610, 36
    max_v = 0.56
    show_order = ["T", "IRxT", "HLUxIRxT", "HLUxT", "HLU", "IR", "HLUxIR", "combined_flight"]
    lookup = {name: (value, count) for name, value, count in data["icp"]}
    for i, name in enumerate(show_order):
        value, count = lookup[name]
        y = chart_y + i * 53
        label = "flight ref." if name == "combined_flight" else name
        color = COLORS["dim"] if name == "combined_flight" else (COLORS["teal"] if "T" in name else COLORS["blue"])
        if "x" in name and name != "combined_flight":
            color = COLORS["rose"]
        text(draw, (x1 + 70, y + 5), label, F["tiny_bold"], COLORS["text"])
        rounded(draw, (chart_x, y, chart_x + bar_w, y + bar_h), 13, COLORS["panel3"], "#2A394D", 1)
        fill_w = int(bar_w * (value / max_v))
        rounded(draw, (chart_x, y, chart_x + fill_w, y + bar_h), 13, color)
        text(draw, (chart_x + fill_w + 16, y + 5), f"{value:.3f}", F["tiny_bold"], COLORS["muted"])
        text(draw, (x2 - 72, y + 6), f"n={count:,}", F["micro_bold"], COLORS["dim"], "ra")

    rounded(draw, (x1 + 70, y2 - 188, x1 + 420, y2 - 20), 26, COLORS["panel2"], COLORS["teal"], 2)
    text(draw, (x1 + 104, y2 - 158), "top mean ICP", F["micro_bold"], COLORS["dim"])
    text(draw, (x1 + 104, y2 - 122), "T 0.540", F["metric2"], COLORS["text"])
    text(draw, (x1 + 104, y2 - 56), "time-rich axis", F["tiny_bold"], COLORS["muted"])

    brain = data["brain_sig"]
    rounded(draw, (x1 + 455, y2 - 188, x2 - 70, y2 - 20), 26, "#122234", COLORS["violet"], 2)
    text(draw, (x1 + 488, y2 - 158), "brain analog context", F["micro_bold"], COLORS["dim"])
    text(draw, (x1 + 488, y2 - 118), f"T genes {brain['T']:,}", F["small_bold"], COLORS["text"])
    text(draw, (x1 + 488, y2 - 76), f"HLUxIRxT genes {brain['HLUxIRxT']:,}", F["tiny_bold"], COLORS["muted"])


def draw_flow_strip(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1530, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "DECOMPOSE decoder", F["h2"], COLORS["text"])
    steps = [
        ("Factorial analogs", "HLU, radiation, time, sex", COLORS["blue"]),
        ("Gene-level model", "fit main + interaction terms", COLORS["amber"]),
        ("Variance share", "assign top-gene signal", COLORS["rose"]),
        ("ICP stability", "rank stable axes", COLORS["teal"]),
        ("Causal map", "carry axes into v8 follow-up", COLORS["green"]),
    ]
    node_w, gap = 570, 62
    start_x, y = x1 + 460, y1 + 152
    for i, (title, desc, color) in enumerate(steps):
        nx = start_x + i * (node_w + gap)
        rounded(draw, (nx, y, nx + node_w, y + 94), 20, COLORS["panel2"], color, 2)
        text(draw, (nx + 26, y + 17), title, F["small_bold"], COLORS["text"])
        text(draw, (nx + 26, y + 55), desc, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            arrow(draw, nx + node_w + 9, y + 47, nx + node_w + gap - 14, y + 47, COLORS["dim"], 4)
    text(draw, (x2 - 80, y1 + 72), "DECOMPOSE", F["metric2"], COLORS["rose"], "ra")
    text(draw, (x2 - 80, y1 + 124), "main axes plus interaction surfaces", F["small_bold"], COLORS["muted"], "ra")


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["rose"])
    paragraph(
        draw,
        (160, 1978),
        "Combined-stressor terms carry a large share of selected response-surface variance and rank stable axes for follow-up analysis.",
        F["small"],
        2520,
        COLORS["muted"],
        8,
    )
    text(draw, (3580, 1960), "Inputs", F["micro_bold"], COLORS["dim"], "ra")
    text(draw, (3580, 1992), "v8 Figure 2 + factorial decomposition JSON", F["tiny"], COLORS["muted"], "ra")
    text(draw, (3580, 2024), "ICP DAG summary", F["tiny"], COLORS["muted"], "ra")


def build() -> None:
    data = load_data()
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 46 | ACT 5 | DECOMPOSE", F["kicker"], COLORS["rose"])
    bx = 1910
    bx = badge(draw, bx, 56, "ANALOGS", "4 programs", COLORS["blue"])
    bx = badge(draw, bx, 56, "INTERACTION", "40-62%", COLORS["rose"])
    bx = badge(draw, bx, 56, "TOP ICP", "T 0.540", COLORS["teal"])
    badge(draw, bx, 56, "LAYER", "v8 Figure 2", COLORS["amber"])

    text(draw, (120, 330), "Stressor Interactions Dominate Combined Effects", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "DECOMPOSE uses factorial analog studies to separate main stressor axes from combined-response surfaces.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_method_panel(draw)
    draw_surface_panel(draw, data)
    draw_icp_panel(draw, data)
    draw_flow_strip(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    ImageOps.colorize(gray, black="#080D14", white="#F3F7FC").save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(gray)

    surfaces = data["surfaces"]
    metrics = {
        "spleen_thymus_interaction_pct": pct(surfaces[0][1][2][1]),
        "skin_skin_interaction_pct": pct(surfaces[1][1][2][1]),
        "hze_endocrine_thymus_interaction_pct": pct(surfaces[2][1][2][1]),
        "top_icp_mean": 0.5401,
        "brain_T_genes": data["brain_sig"]["T"],
        "brain_HLUxIRxT_genes": data["brain_sig"]["HLUxIRxT"],
    }
    manifest = {
        "title": "Stressor Interactions Dominate Combined Effects",
        "readout": "Combined-stressor terms carry a large variance share in selected response surfaces and identify stable follow-up axes.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "source_manifests": {"v8_pillars": str(V8_PILLARS_MANIFEST.relative_to(ROOT))},
        "source_files": {
            "v8_figure2": str(V8_FIGURE.relative_to(ROOT)),
            "factorial_decomposition": str(FACTORIAL.relative_to(ROOT)),
            "icp_dag_summary": str(ICP.relative_to(ROOT)),
        },
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": metrics,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    QA_NOTE.write_text(
        "\n".join(
            [
                "# DECOMPOSE Interactions Asset Visual QA",
                "",
                "Slide 46 rebuilds v8 Figure 2 into a presentation-scale interaction readout.",
                "",
                "Checks to review:",
                "- Method panel explains main axes and interaction surface in one read.",
                "- Stacked bars preserve clear labels for each variance segment.",
                "- ICP panel keeps row labels, values, and n counts separated.",
                "- Bottom DECOMPOSE reading strip has stable arrow spacing.",
                "- Grayscale version preserves the interaction and ICP hierarchy.",
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

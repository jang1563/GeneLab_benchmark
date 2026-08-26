#!/usr/bin/env python3
"""Build slide 47 asset: perturbation hits prioritize follow-up axes."""

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
OUT_DIR = ASSET_ROOT / "intervene_prioritization"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "perturbation_hits_prioritize_followup_axes_premium.png"
GRAY_PATH = OUT_DIR / "perturbation_hits_prioritize_followup_axes_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "intervene_prioritization_manifest.json"
QA_NOTE = OUT_DIR / "INTERVENE_PRIORITIZATION_ASSET_VISUAL_QA.md"

V8_PILLARS_MANIFEST = ASSET_ROOT / "v8_pillars" / "v8_pillars_manifest.json"
EVAL = ROOT / "v8" / "intervene" / "evaluation"
SIGNATURES = ROOT / "v8" / "intervene" / "signatures" / "signatures_manifest.json"
LINCS = EVAL / "lincs_summary.json"
SCORES = EVAL / "multi_tissue_drug_scores.json"
MATRIX = EVAL / "multi_tissue_drug_matrix.csv"
TRIAGE = EVAL / "safety_triage.csv"
CRISPR = EVAL / "crispr_orthog_summary.json"

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
    "title": load_font(82, True),
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
    pts = [(x2, y2), (x2 - 20, y2 - 11), (x2 - 20, y2 + 11)]
    draw.polygon(pts, fill=color)


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def short_name(value: str) -> str:
    names = {
        "QUINACRINE HYDROCHLORIDE": "Quinacrine HCl",
        "Dorsomorphin dihydrochloride": "Dorsomorphin",
        "ALBENDAZOLE": "Albendazole",
    }
    return names.get(value, value)


def axis_label(candidate: str) -> str:
    mapping = {
        "CGP-60474": "CDK cell-cycle axis",
        "QUINACRINE HYDROCHLORIDE": "lysosome / immune axis",
        "AZD-5438": "CDK9-supported axis",
        "geldanamycin": "HSP90 proteostasis axis",
        "mitoxantrone": "TOP2 stress axis",
        "AT-7519": "CDK9-supported axis",
        "Dorsomorphin dihydrochloride": "AMPK/BMP metabolic axis",
        "ALBENDAZOLE": "microtubule thymus axis",
    }
    return mapping.get(candidate, "follow-up axis")


def load_candidates() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with TRIAGE.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            row["candidate_short"] = short_name(str(row["candidate"]))
            row["axis"] = axis_label(str(row["candidate"]))
            row["n_tissues_int"] = int(row["n_tissues"])
            row["mean_rev_float"] = float(row["mean_rev"])
            row["min_rev_float"] = float(row["min_rev"])
            row["max_rev_float"] = float(row["max_rev"])
            row["pareto_bool"] = str(row["pareto_front"]) == "True"
            rows.append(row)
    rows.sort(key=lambda item: (bool(item["pareto_bool"]), int(item["n_tissues_int"]), float(item["mean_rev_float"])), reverse=True)
    return rows


def load_data() -> dict[str, object]:
    signatures = load_json(SIGNATURES)
    scores = load_json(SCORES)
    crispr = load_json(CRISPR)
    lincs = load_json(LINCS)
    candidates = load_candidates()
    kidney_support = crispr["tissues"]["kidney"]["validated_drugs"]
    return {
        "signatures": signatures,
        "scores": scores,
        "lincs": lincs,
        "crispr": crispr,
        "candidates": candidates,
        "kidney_support": kidney_support,
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


def draw_reversal_pair(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, left: str, right: str, color: str) -> None:
    rounded(draw, (x, y, x + 330, y + 86), 20, COLORS["panel2"], color, 2)
    text(draw, (x + 24, y + 18), label, F["micro_bold"], COLORS["dim"])
    text(draw, (x + 24, y + 50), left, F["tiny_bold"], COLORS["text"])
    arrow(draw, x + 165, y + 52, x + 218, y + 52, color, 4)
    text(draw, (x + 236, y + 50), right, F["tiny_bold"], COLORS["text"])


def draw_method_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x1, y1, x2, y2 = 120, 610, 1210, 1462
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 36), "What INTERVENE does", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 92),
        "Signature reversal asks which perturbation profiles move flight-up and flight-down genes in the opposite direction.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )
    top_n = data["signatures"]["top_n"]
    tissues = len(data["signatures"]["tissues"])

    rounded(draw, (x1 + 70, y1 + 228, x1 + 420, y1 + 396), 26, "#122234", COLORS["blue"], 2)
    text(draw, (x1 + 104, y1 + 260), "flight signatures", F["micro_bold"], COLORS["dim"])
    text(draw, (x1 + 104, y1 + 308), f"{tissues} tissues", F["metric2"], COLORS["text"])
    text(draw, (x1 + 104, y1 + 366), f"top-{top_n} up + down", F["tiny_bold"], COLORS["muted"])

    draw_reversal_pair(draw, x1 + 500, y1 + 230, "reversal logic", "flight-up", "pert-down", COLORS["rose"])
    draw_reversal_pair(draw, x1 + 500, y1 + 340, "paired logic", "flight-down", "pert-up", COLORS["teal"])

    steps = [
        ("1", "Export per-tissue ranked genes", COLORS["blue"]),
        ("2", "Query L1000 chemical profiles", COLORS["amber"]),
        ("3", "Score multi-tissue breadth", COLORS["rose"]),
        ("4", "Check CRISPR KO target echoes", COLORS["teal"]),
    ]
    sy = y1 + 474
    for i, (num, body, color) in enumerate(steps):
        y = sy + i * 60
        rounded(draw, (x1 + 74, y, x2 - 74, y + 48), 16, COLORS["panel2"], color, 2)
        rounded(draw, (x1 + 94, y + 8, x1 + 132, y + 42), 11, color)
        text(draw, (x1 + 113, y + 25), num, F["tiny_bold"], COLORS["ink"], "mm")
        text(draw, (x1 + 154, y + 13), body, F["tiny_bold"], COLORS["text"])

    rounded(draw, (x1 + 64, y2 - 128, x2 - 64, y2 - 42), 24, "#122234", COLORS["green"], 2)
    paragraph(
        draw,
        (x1 + 96, y2 - 104),
        "Perturbation hits become pathway axes for follow-up experiments.",
        F["small"],
        x2 - x1 - 192,
        COLORS["muted"],
        6,
    )


def scatter_xy(x: float, y: float, box: tuple[int, int, int, int]) -> tuple[int, int]:
    x1, y1, x2, y2 = box
    xmin, xmax = 0.055, 0.095
    ymin, ymax = 0.055, 0.155
    px = x1 + int((x - xmin) / (xmax - xmin) * (x2 - x1))
    py = y2 - int((y - ymin) / (ymax - ymin) * (y2 - y1))
    return px, py


def draw_pareto_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x1, y1, x2, y2 = 1265, 610, 2575, 1462
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 36), "Multi-tissue scoring narrows the list", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 92),
        "A Pareto view keeps both average reversal and weakest-tissue reversal visible.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )

    metrics = [
        ("unique hits", str(data["scores"]["n_drugs_total"]), COLORS["blue"]),
        ("multi-tissue", str(data["scores"]["n_drugs_multi_tissue"]), COLORS["teal"]),
        ("Pareto front", str(data["scores"]["n_pareto_front"]), COLORS["rose"]),
    ]
    mx, my = x1 + 54, y1 + 178
    for i, (label, value, color) in enumerate(metrics):
        xx = mx + i * 282
        rounded(draw, (xx, my, xx + 246, my + 86), 21, COLORS["panel2"], color, 2)
        text(draw, (xx + 24, my + 18), label, F["micro_bold"], COLORS["dim"])
        text(draw, (xx + 24, my + 60), value, F["metric2"], COLORS["text"], "lm")

    chart = (x1 + 120, y1 + 345, x2 - 120, y2 - 160)
    cx1, cy1, cx2, cy2 = chart
    draw.rectangle(chart, fill=COLORS["panel3"], outline="#2A394D", width=2)
    for i in range(1, 5):
        xx = cx1 + i * (cx2 - cx1) // 5
        draw.line((xx, cy1, xx, cy2), fill=rgba(COLORS["grid"], 95), width=1)
    for i in range(1, 5):
        yy = cy1 + i * (cy2 - cy1) // 5
        draw.line((cx1, yy, cx2, yy), fill=rgba(COLORS["grid"], 95), width=1)

    text(draw, ((cx1 + cx2) / 2, cy2 + 36), "minimum observed tissue reversal", F["tiny_bold"], COLORS["muted"], "mm")
    text(draw, (cx1 - 54, (cy1 + cy2) / 2), "mean reversal", F["tiny_bold"], COLORS["muted"], "mm")
    text(draw, (cx1, cy2 + 8), "0.055", F["axis"], COLORS["dim"])
    text(draw, (cx2, cy2 + 8), "0.095", F["axis"], COLORS["dim"], "ra")
    text(draw, (cx1 - 12, cy2), "0.055", F["axis"], COLORS["dim"], "ra")
    text(draw, (cx1 - 12, cy1), "0.155", F["axis"], COLORS["dim"], "ra")

    candidates = data["candidates"]
    label_points: list[tuple[str, int, int, str]] = []
    for row in candidates:
        x = float(row["min_rev_float"])
        y = float(row["mean_rev_float"])
        px, py = scatter_xy(x, y, chart)
        n_tissues = int(row["n_tissues_int"])
        radius = 17 + n_tissues * 7
        color = COLORS["rose"] if row["pareto_bool"] else (COLORS["blue"] if n_tissues >= 4 else COLORS["amber"])
        draw.ellipse((px - radius, py - radius, px + radius, py + radius), fill=rgba(color, 215), outline=color, width=3)
        if row["pareto_bool"] or row["candidate_short"] in {"AZD-5438", "geldanamycin"}:
            label_points.append((str(row["candidate_short"]), px, py, color))

    label_positions = {
        "AZD-5438": (84, -88),
        "geldanamycin": (96, 24),
        "CGP-60474": (54, -62),
        "Quinacrine HCl": (-270, 46),
    }
    for label, px, py, color in label_points:
        dx, dy = label_positions[label]
        lx, ly = px + dx, py + dy
        w = int(draw.textlength(label, font=F["tiny_bold"]) + 34)
        draw.line((px, py, lx + 12, ly + 19), fill=color, width=3)
        rounded(draw, (lx, ly, lx + w, ly + 40), 13, COLORS["panel2"], color, 2)
        text(draw, (lx + 17, ly + 10), label, F["tiny_bold"], COLORS["text"])

    rounded(draw, (x1 + 70, y2 - 116, x2 - 70, y2 - 42), 22, "#122234", COLORS["rose"], 2)
    paragraph(
        draw,
        (x1 + 104, y2 - 96),
        "Pareto front highlights CGP-60474 and Quinacrine HCl; AZD-5438 and HSP90 hits add broad tissue-axis context.",
        F["small"],
        x2 - x1 - 208,
        COLORS["muted"],
        6,
    )


def draw_axis_row(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    title: str,
    evidence: str,
    color: str,
    tissues: int,
) -> None:
    rounded(draw, (x, y, x + w, y + 74), 20, COLORS["panel2"], color, 2)
    text(draw, (x + 24, y + 16), title, F["tiny_bold"], COLORS["text"])
    text(draw, (x + 24, y + 47), evidence, F["micro_bold"], COLORS["muted"])
    rounded(draw, (x + w - 126, y + 17, x + w - 28, y + 56), 13, color)
    text(draw, (x + w - 77, y + 36), f"{tissues}T", F["tiny_bold"], COLORS["ink"], "mm")


def draw_crispr_grid(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    tissues = ["thymus", "gastroc.", "skin", "eye", "liver", "kidney"]
    cell_w, cell_h = 78, 54
    text(draw, (x, y - 34), "CRISPR KO echo", F["tiny_bold"], COLORS["text"])
    text(draw, (x + 480, y - 34), "CDK9", F["tiny_bold"], COLORS["rose"], "ra")
    for i, tissue in enumerate(tissues):
        tx = x + i * cell_w
        text(draw, (tx + cell_w / 2, y), tissue, F["axis"], COLORS["dim"], "mm")
        fill = COLORS["rose"] if tissue == "kidney" else COLORS["panel3"]
        outline = COLORS["rose"] if tissue == "kidney" else "#2A394D"
        rounded(draw, (tx, y + 24, tx + cell_w - 8, y + 24 + cell_h), 10, fill, outline, 2)
        value = "2" if tissue == "kidney" else "0"
        value_fill = COLORS["text"] if tissue == "kidney" else COLORS["dim"]
        text(draw, (tx + (cell_w - 8) / 2, y + 24 + cell_h / 2), value, F["small_bold"], value_fill, "mm")


def draw_target_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x1, y1, x2, y2 = 2630, 610, 3720, 1462
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 36), "Hits become target-axis priorities", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 92),
        "Chemical and CRISPR evidence are read together as pathway-axis triage for the next experiment.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )

    rows = [
        ("CDK cell-cycle", "4-tissue breadth + CDK9 echo", COLORS["rose"], 4),
        ("HSP90 proteostasis", "4-tissue chemical reversal", COLORS["blue"], 4),
        ("AMPK/BMP metabolic", "gastrocnemius + liver signal", COLORS["amber"], 3),
        ("mTOR / MEK renal", "kidney stress-gene axis", COLORS["teal"], 2),
    ]
    yy = y1 + 224
    for i, row in enumerate(rows):
        draw_axis_row(draw, x1 + 62, yy + i * 96, x2 - x1 - 124, *row)

    draw_crispr_grid(draw, x1 + 82, y1 + 658)

    support = data["kidney_support"]["AZD-5438"][0]
    rounded(draw, (x1 + 610, y1 + 625, x2 - 62, y2 - 42), 24, "#122234", COLORS["violet"], 2)
    text(draw, (x1 + 646, y1 + 660), "example readout", F["micro_bold"], COLORS["dim"])
    text(draw, (x1 + 646, y1 + 694), "CDK9 KO", F["metric2"], COLORS["text"])
    text(draw, (x1 + 646, y1 + 750), "kidney UP genes", F["tiny_bold"], COLORS["muted"])
    text(draw, (x1 + 646, y1 + 780), f"7 overlaps | adj-p {support['adj_p']:.3f}", F["micro_bold"], COLORS["muted"])


def draw_flow_strip(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1530, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "INTERVENE Prioritization Flow", F["h2"], COLORS["text"])
    steps = [
        ("Flight DGE", "per-tissue ranked genes", COLORS["blue"]),
        ("Signatures", "top-150 up and down", COLORS["teal"]),
        ("Reversal query", "chemical profiles", COLORS["amber"]),
        ("Pareto triage", "mean plus weak tissue", COLORS["rose"]),
        ("CRISPR echo", "target-axis check", COLORS["violet"]),
        ("Follow-up axes", "next experiment axes", COLORS["green"]),
    ]
    node_w, gap = 490, 34
    start_x, y = x1 + 280, y1 + 152
    for i, (title, desc, color) in enumerate(steps):
        nx = start_x + i * (node_w + gap)
        rounded(draw, (nx, y, nx + node_w, y + 94), 20, COLORS["panel2"], color, 2)
        text(draw, (nx + 24, y + 17), title, F["small_bold"], COLORS["text"])
        text(draw, (nx + 24, y + 55), desc, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            arrow(draw, nx + node_w + 6, y + 47, nx + node_w + gap - 10, y + 47, COLORS["dim"], 4)
    text(draw, (x2 - 80, y1 + 72), "INTERVENE", F["metric2"], COLORS["amber"], "ra")
    text(draw, (x2 - 80, y1 + 124), "perturbation evidence to target axes", F["small_bold"], COLORS["muted"], "ra")


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["amber"])
    paragraph(
        draw,
        (160, 1978),
        "INTERVENE concentrates 215 perturbation hits into a short list of pathway axes: CDK/cell-cycle, HSP90/proteostasis, AMPK/BMP metabolism, and renal mTOR/MEK.",
        F["small"],
        3180,
        COLORS["muted"],
        8,
    )


def build() -> None:
    data = load_data()
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 47 | ACT 5 | INTERVENE", F["kicker"], COLORS["amber"])
    bx = 1840
    bx = badge(draw, bx, 56, "SIGNATURES", "6 tissues", COLORS["blue"])
    bx = badge(draw, bx, 56, "HITS", f"{data['scores']['n_drugs_total']} perturbagens", COLORS["amber"])
    bx = badge(draw, bx, 56, "MULTI-TISSUE", str(data["scores"]["n_drugs_multi_tissue"]), COLORS["teal"])
    badge(draw, bx, 56, "PARETO", f"{data['scores']['n_pareto_front']} front axes", COLORS["rose"])

    text(draw, (120, 330), "Perturbation Hits Prioritize Follow-Up Axes", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "INTERVENE turns tissue flight signatures into chemical and CRISPR perturbation evidence for ranked pathway axes.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_method_panel(draw, data)
    draw_pareto_panel(draw, data)
    draw_target_panel(draw, data)
    draw_flow_strip(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    ImageOps.colorize(gray, black="#080D14", white="#F3F7FC").save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(gray)

    manifest = {
        "title": "Perturbation Hits Prioritize Follow-Up Axes",
        "claim": "INTERVENE prioritizes perturbation evidence into pathway axes for follow-up experiments.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "source_manifests": {"v8_pillars": str(V8_PILLARS_MANIFEST.relative_to(ROOT))},
        "source_files": {
            "signatures_manifest": str(SIGNATURES.relative_to(ROOT)),
            "lincs_summary": str(LINCS.relative_to(ROOT)),
            "multi_tissue_scores": str(SCORES.relative_to(ROOT)),
            "multi_tissue_matrix": str(MATRIX.relative_to(ROOT)),
            "triage_table": str(TRIAGE.relative_to(ROOT)),
            "crispr_summary": str(CRISPR.relative_to(ROOT)),
        },
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": {
            "signature_tissues": len(data["signatures"]["tissues"]),
            "signature_top_n": data["signatures"]["top_n"],
            "unique_hits": data["scores"]["n_drugs_total"],
            "multi_tissue_hits": data["scores"]["n_drugs_multi_tissue"],
            "pareto_front": data["scores"]["n_pareto_front"],
            "crispr_supported_compound_target_pairs": 2,
            "cdk9_overlap_genes": 7,
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    QA_NOTE.write_text(
        "\n".join(
            [
                "# INTERVENE Prioritization Asset Visual QA",
                "",
                "Slide 47 rebuilds v8 Figure 4 into a presentation-scale perturbation-prioritization proof object.",
                "",
                "Checks to review:",
                "- Reversal logic explains flight-up and flight-down directions in one read.",
                "- Pareto scatter labels separate front candidates and broad tissue-axis context.",
                "- Target-axis rows and CRISPR grid stay readable at presentation scale.",
                "- Bottom INTERVENE reading strip has stable arrow spacing.",
                "- Grayscale version preserves Pareto and target-axis hierarchy.",
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

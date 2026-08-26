#!/usr/bin/env python3
"""Build slide 44 asset: v8 has four incubator pillars."""

from __future__ import annotations

import csv
import json
import math
import re
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
OUT_DIR = ASSET_ROOT / "v8_pillars"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "v8_has_four_incubator_pillars_premium.png"
GRAY_PATH = OUT_DIR / "v8_has_four_incubator_pillars_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "v8_pillars_manifest.json"
QA_NOTE = OUT_DIR / "V8_PILLARS_ASSET_VISUAL_QA.md"

V8_TRANSITION_MANIFEST = ASSET_ROOT / "v8_transition" / "v8_transition_manifest.json"
V8_SUMMARY = ROOT / "v8" / "RESULTS_SUMMARY.md"
V8_BRIDGE = ROOT / "v8" / "bridge" / "evaluation" / "supervised_conservation.json"
V8_SPECIES = ROOT / "v8" / "bridge" / "evaluation" / "species_transfer_nes.json"
V8_TISSUE_NES = ROOT / "v8" / "bridge" / "evaluation" / "tissue_nes_spearman.json"
V8_FACTORIAL = ROOT / "v8" / "decompose" / "evaluation" / "factorial_flight_decomposition.json"
V8_MARS = ROOT / "v8" / "decompose" / "evaluation" / "mars_saturation_summary.json"
V8_OFFLINE = ROOT / "v8" / "intervene" / "evaluation" / "offline_reversal_summary.json"
V8_MULTI_TISSUE = ROOT / "v8" / "intervene" / "evaluation" / "multi_tissue_drug_scores.json"
V8_PARETO = ROOT / "v8" / "intervene" / "evaluation" / "pareto_front.csv"
V8_CRISPR = ROOT / "v8" / "intervene" / "evaluation" / "crispr_orthog_summary.json"
V8_CAUSAL = ROOT / "v8" / "causal" / "evaluation" / "icp_dag_summary.json"
V8_DAG_EDGES = ROOT / "v8" / "causal" / "evaluation" / "dag_edges.csv"

COLORS = {
    "bg": "#0B1119",
    "bg2": "#10141F",
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
    "title": load_font(86, True),
    "subtitle": load_font(37),
    "section": load_font(42, True),
    "h2": load_font(36, True),
    "h3": load_font(29, True),
    "body": load_font(27),
    "body_bold": load_font(27, True),
    "small": load_font(24),
    "small_bold": load_font(24, True),
    "tiny": load_font(21),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "metric": load_font(58, True),
    "metric2": load_font(46, True),
    "label": load_font(16, True),
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
    pts = [(x2, y2), (x2 - 20, y2 - 11), (x2 - 20, y2 + 11)] if x2 >= x1 else [(x2, y2), (x2 + 20, y2 - 11), (x2 + 20, y2 + 11)]
    draw.polygon(pts, fill=color)


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_interaction_range() -> str:
    summary = V8_SUMMARY.read_text(encoding="utf-8")
    match = re.search(r"Interaction dominance:\s*([0-9]+)[–-]([0-9]+)%", summary)
    return f"{match.group(1)}-{match.group(2)}%" if match else "44-61%"


def count_pareto_rows() -> int:
    with V8_PARETO.open(newline="", encoding="utf-8") as handle:
        return max(0, sum(1 for _ in csv.DictReader(handle)))


def top_dag_edges(limit: int = 6) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with V8_DAG_EDGES.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "from": row["node_from"].replace("stressor:", ""),
                    "to": row["node_to"].replace("tissue:", ""),
                    "score": float(row["icp_score"]),
                    "n": int(row["n_genes"]),
                }
            )
    return rows[:limit]


def load_data() -> dict[str, object]:
    bridge = load_json(V8_BRIDGE)
    species = load_json(V8_SPECIES)
    tissue_nes = load_json(V8_TISSUE_NES)
    factorial = load_json(V8_FACTORIAL)
    mars = load_json(V8_MARS)
    offline = load_json(V8_OFFLINE)
    multi_tissue = load_json(V8_MULTI_TISSUE)
    crispr = load_json(V8_CRISPR)
    causal = load_json(V8_CAUSAL)
    transition = load_json(V8_TRANSITION_MANIFEST)
    return {
        "bridge": bridge,
        "species": species,
        "tissue_nes": tissue_nes,
        "factorial": factorial,
        "mars": mars,
        "offline": offline,
        "multi_tissue": multi_tissue,
        "crispr": crispr,
        "causal": causal,
        "transition": transition,
        "interaction_range": parse_interaction_range(),
        "pareto_count": count_pareto_rows(),
        "dag_edges": top_dag_edges(),
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


def metric_color(value: float) -> str:
    if value >= 0.35:
        return COLORS["green"]
    if value >= 0.0:
        return COLORS["teal"]
    if value >= -0.35:
        return COLORS["amber"]
    return COLORS["rose"]


def draw_bridge_visual(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object], color: str) -> None:
    x1, y1, x2, y2 = box
    pivot = data["tissue_nes"]["i4_pivot"]
    tissues = ["gastrocnemius", "kidney", "skin", "liver", "thymus", "eye"]
    comps = list(pivot.keys())[:9]
    cell = 22
    gx, gy = x1 + 36, y1 + 42
    for r, tissue in enumerate(tissues):
        for c, comp in enumerate(comps):
            value = float(pivot[comp][tissue])
            rounded(draw, (gx + c * (cell + 6), gy + r * (cell + 6), gx + c * (cell + 6) + cell, gy + r * (cell + 6) + cell), 4, metric_color(value), None, 1)
    text(draw, (gx, gy + 196), "mouse NES by I4 compartment", F["micro_bold"], COLORS["muted"])
    ax = x1 + 328
    ay = y1 + 42
    bar_w = 320
    base = data["bridge"]["rf"]["cv_mean_auroc_base"]
    aug = data["bridge"]["rf"]["cv_mean_auroc_aug"]
    for i, (label, value, fill) in enumerate([("base", base, COLORS["dim"]), ("+ mouse NES", aug, color)]):
        y = ay + i * 74
        rounded(draw, (ax, y, ax + bar_w, y + 34), 14, COLORS["panel3"], "#2A394D", 1)
        rounded(draw, (ax, y, ax + int(bar_w * value), y + 34), 14, fill)
        text(draw, (ax, y + 44), label, F["micro_bold"], COLORS["muted"])
        text(draw, (ax + bar_w, y + 44), f"{value:.3f}", F["micro_bold"], COLORS["text"], "ra")
    text(draw, (ax, ay + 166), f"+{data['bridge']['rf']['delta_aug_minus_base']['mean']:.3f} AUROC lift", F["tiny_bold"], color)


def draw_decompose_visual(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object], color: str) -> None:
    x1, y1, x2, y2 = box
    cx, cy = x1 + 242, y1 + 120
    nodes = [("HLU", x1 + 88, y1 + 58, COLORS["blue"]), ("IR", x1 + 396, y1 + 58, COLORS["teal"]), ("T", x1 + 242, y1 + 218, COLORS["amber"])]
    for label, nx, ny, fill in nodes:
        draw.line((nx, ny, cx, cy), fill=rgba(fill, 170), width=5)
        draw.ellipse((nx - 46, ny - 46, nx + 46, ny + 46), fill=rgba(fill, 45), outline=fill, width=3)
        text(draw, (nx, ny - 13), label, F["tiny_bold"], COLORS["text"], "mm")
        text(draw, (nx, ny + 15), "axis", F["micro_bold"], COLORS["muted"], "mm")
    draw.ellipse((cx - 64, cy - 64, cx + 64, cy + 64), fill=rgba(color, 55), outline=color, width=4)
    text(draw, (cx, cy - 10), "interaction", F["tiny_bold"], COLORS["text"], "mm")
    text(draw, (cx, cy + 20), data["interaction_range"], F["tiny_bold"], color, "mm")
    bx, by = x1 + 548, y1 + 62
    text(draw, (bx, by), "Mars stressor vector", F["micro_bold"], COLORS["muted"])
    doses = [("IR", data["mars"]["linear_dose"]["IR"], COLORS["teal"]), ("T", data["mars"]["linear_dose"]["T"], COLORS["amber"]), ("HLU", data["mars"]["linear_dose"]["HLU"], COLORS["blue"])]
    for i, (label, value, fill) in enumerate(doses):
        y = by + 40 + i * 48
        rounded(draw, (bx, y, bx + 230, y + 28), 12, COLORS["panel3"], "#2A394D", 1)
        rounded(draw, (bx, y, bx + int(min(1.0, value / 13.0) * 230), y + 28), 12, fill)
        text(draw, (bx + 246, y + 4), f"{label} {value:.1f}x", F["micro_bold"], COLORS["text"])


def draw_intervene_visual(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object], color: str) -> None:
    x1, y1, x2, y2 = box
    funnel = [
        ("signatures", data["offline"]["drugs_scored"], COLORS["blue"]),
        ("multi-tissue", data["multi_tissue"]["n_drugs_multi_tissue"], COLORS["teal"]),
        ("priority axes", data["pareto_count"], color),
    ]
    y = y1 + 44
    widths = [520, 400, 250]
    for i, (label, value, fill) in enumerate(funnel):
        w = widths[i]
        x = x1 + 70 + (520 - w) // 2
        rounded(draw, (x, y, x + w, y + 54), 18, rgba(fill, 45), fill, 2)
        text(draw, (x + 26, y + 15), label, F["tiny_bold"], COLORS["text"])
        text(draw, (x + w - 26, y + 15), f"{value}", F["tiny_bold"], fill, "ra")
        if i < len(funnel) - 1:
            arrow(draw, x1 + 330, y + 62, x1 + 330, y + 96, COLORS["dim"], 4)
        y += 104
    sx, sy = x1 + 610, y1 + 50
    rounded(draw, (sx, sy, sx + 196, sy + 196), 20, COLORS["panel3"], "#2A394D", 2)
    draw.line((sx + 42, sy + 154, sx + 168, sy + 154), fill=COLORS["dim"], width=2)
    draw.line((sx + 42, sy + 154, sx + 42, sy + 34), fill=COLORS["dim"], width=2)
    points = [(82, 92, COLORS["dim"]), (126, 68, color), (108, 118, COLORS["teal"]), (146, 104, color), (70, 126, COLORS["dim"])]
    for px, py, fill in points:
        draw.ellipse((sx + px - 8, sy + py - 8, sx + px + 8, sy + py + 8), fill=fill)
    text(draw, (sx + 54, sy + 166), "mean", F["label"], COLORS["muted"])
    text(draw, (sx + 8, sy + 38), "min", F["label"], COLORS["muted"])


def draw_causal_visual(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object], color: str) -> None:
    x1, y1, x2, y2 = box
    stressors = ["IR", "HLU", "HLUxIR"]
    tissues = ["skin", "thymus", "kidney", "eye"]
    sx = x1 + 72
    tx = x1 + 458
    for i, label in enumerate(stressors):
        y = y1 + 55 + i * 72
        rounded(draw, (sx, y, sx + 146, y + 44), 16, rgba(COLORS["blue"], 35), COLORS["blue"], 2)
        text(draw, (sx + 73, y + 13), label, F["tiny_bold"], COLORS["text"], "ma")
    for i, label in enumerate(tissues):
        y = y1 + 34 + i * 56
        rounded(draw, (tx, y, tx + 170, y + 40), 15, rgba(color, 35), color, 2)
        text(draw, (tx + 85, y + 12), label, F["micro_bold"], COLORS["text"], "ma")
    tissue_y = {"skin": y1 + 54, "thymus": y1 + 110, "kidney": y1 + 166, "eye": y1 + 222}
    stressor_y = {"IR": y1 + 77, "HLU": y1 + 149, "HLUxIR": y1 + 221}
    for edge in data["dag_edges"][:6]:
        source = str(edge["from"])
        target = str(edge["to"])
        if source in stressor_y and target in tissue_y:
            arrow(draw, sx + 146, stressor_y[source], tx, tissue_y[target] + 20, rgba(color, 190), 3)
    rounded(draw, (x1 + 684, y1 + 60, x1 + 808, y1 + 184), 24, rgba(color, 40), color, 3)
    text(draw, (x1 + 746, y1 + 92), f"{data['causal']['n_environments']}", F["metric2"], color, "mm")
    text(draw, (x1 + 746, y1 + 132), "env", F["tiny_bold"], COLORS["text"], "mm")
    text(draw, (x1 + 746, y1 + 160), f"T ICP {data['causal']['stressor_aggregate_icp']['T']['mean']:.3f}", F["micro_bold"], COLORS["muted"], "mm")


def draw_info_row(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str, width: int) -> None:
    text(draw, (x, y), label.upper(), F["micro_bold"], COLORS["dim"])
    paragraph(draw, (x, y + 27), value, F["tiny"], width, COLORS["muted"], 5)
    draw.line((x, y + 84, x + width, y + 84), fill=rgba(color, 120), width=2)


def draw_pillar_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    index: str,
    tag: str,
    title: str,
    question: str,
    visual_fn,
    rows: list[tuple[str, str]],
    metric: str,
    next_slide: str,
    color: str,
    data: dict[str, object],
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 30, COLORS["panel"], "#2B3A50", 2)
    rounded(draw, (x1, y1, x2, y1 + 14), 7, color)
    text(draw, (x1 + 38, y1 + 44), index, F["metric2"], color)
    rounded(draw, (x1 + 138, y1 + 45, x1 + 322, y1 + 91), 15, rgba(color, 45), color, 2)
    text(draw, (x1 + 230, y1 + 59), tag, F["micro_bold"], COLORS["text"], "ma")
    text(draw, (x1 + 38, y1 + 114), title, F["h2"], COLORS["text"])
    paragraph(draw, (x1 + 38, y1 + 164), question, F["small"], x2 - x1 - 76, COLORS["muted"], 7)
    rounded(draw, (x1 + 38, y1 + 258, x2 - 38, y1 + 560), 24, "#0E1724", "#2A394D", 2)
    visual_fn(draw, (x1 + 38, y1 + 258, x2 - 38, y1 + 560), data, color)
    text(draw, (x1 + 42, y1 + 596), metric, F["metric"], color)
    row_y = y1 + 646
    for label, value in rows:
        draw_info_row(draw, x1 + 42, row_y, label, value, color, x2 - x1 - 84)
        row_y += 102
    rounded(draw, (x1 + 42, y2 - 92, x2 - 42, y2 - 38), 18, rgba(color, 35), color, 2)
    text(draw, ((x1 + x2) / 2, y2 - 66), next_slide, F["tiny_bold"], COLORS["text"], "mm")


def draw_pillars(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x0, y0 = 120, 610
    card_w, card_h, gap = 873, 1048, 36
    bridge = data["bridge"]
    modules = [
        (
            "01",
            "BRIDGE",
            "Human Pathway Link",
            "Which mouse pathway signals help predict human spaceflight conservation?",
            draw_bridge_visual,
            [
                ("Evidence intake", f"{bridge['n_pathways']} pathways; {data['species']['n_comparisons']} species comparisons"),
                ("Method move", "augment I4 pathway features with six mouse NES inputs"),
                ("Readout", f"RF AUROC {bridge['rf']['cv_mean_auroc_aug']:.3f}; lift {bridge['rf']['delta_aug_minus_base']['mean']:.3f}"),
            ],
            "AUROC 0.888",
            "Next: slide 45",
            COLORS["blue"],
        ),
        (
            "02",
            "DECOMPOSE",
            "Stressor Interaction Surface",
            "Which stressor combinations create the response surface for Mars-like questions?",
            draw_decompose_visual,
            [
                ("Evidence intake", "HLU, radiation, time, and interaction terms from analog studies"),
                ("Method move", "factorial coefficients become stressor-specific surfaces"),
                ("Readout", f"interaction variance {data['interaction_range']}; Mars vector includes IR 12.9x"),
            ],
            data["interaction_range"],
            "Next: slide 46",
            COLORS["teal"],
        ),
        (
            "03",
            "INTERVENE",
            "Perturbation Priority Queue",
            "Which perturbation axes should become follow-up experiment questions?",
            draw_intervene_visual,
            [
                ("Evidence intake", "tissue signatures queried against LINCS and CRISPR resources"),
                ("Method move", "rank reversal signatures and keep multi-tissue priority axes"),
                ("Readout", f"{data['multi_tissue']['n_drugs_total']} candidates; {data['pareto_count']} Pareto axes"),
            ],
            f"{data['pareto_count']} axes",
            "Next: slide 47",
            COLORS["amber"],
        ),
        (
            "04",
            "CAUSAL",
            "Assumption Map",
            "Which stressor-to-tissue assumptions organize the translational readout?",
            draw_causal_visual,
            [
                ("Evidence intake", f"{data['causal']['n_environments']} environments and {data['causal']['n_genes']:,} genes"),
                ("Method move", "rank stable stressor-to-tissue edges with ICP"),
                ("Readout", f"T mean ICP {data['causal']['stressor_aggregate_icp']['T']['mean']:.3f}; top edge IR to skin"),
            ],
            f"{data['causal']['n_environments']} environments",
            "Next: slide 48",
            COLORS["violet"],
        ),
    ]
    for i, module in enumerate(modules):
        x = x0 + i * (card_w + gap)
        draw_pillar_card(draw, (x, y0, x + card_w, y0 + card_h), *module, data)


def draw_contract_strip(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1710, 3720, 1858
    rounded(draw, (x1, y1, x2, y2), 28, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 36), "Next four slides decoder", F["h2"], COLORS["text"])
    steps = [
        ("Pillar", "module question", COLORS["violet"]),
        ("Input", "evidence intake", COLORS["blue"]),
        ("Method", "analysis move", COLORS["teal"]),
        ("Readout", "bounded result", COLORS["green"]),
        ("Follow-up", "experiment queue", COLORS["amber"]),
    ]
    start_x, node_w, gap = x1 + 850, 455, 48
    y = y1 + 48
    for i, (title, desc, color) in enumerate(steps):
        nx = start_x + i * (node_w + gap)
        rounded(draw, (nx, y, nx + node_w, y + 72), 18, COLORS["panel2"], color, 2)
        text(draw, (nx + 24, y + 13), title, F["tiny_bold"], COLORS["text"])
        text(draw, (nx + 24, y + 42), desc, F["micro"], COLORS["muted"])
        if i < len(steps) - 1:
            arrow(draw, nx + node_w + 8, y + 36, nx + node_w + gap - 14, y + 36, COLORS["dim"], 3)


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "Act 5 map: each v8 pillar converts an established evidence lane into a specific follow-up analysis question.",
        F["small"],
        2470,
        COLORS["muted"],
        8,
    )
    text(draw, (3580, 1960), "Inputs", F["micro_bold"], COLORS["dim"], "ra")
    text(draw, (3580, 1992), "v8 evaluation JSON/CSV", F["tiny"], COLORS["muted"], "ra")
    text(draw, (3580, 2024), "transition manifest + v8 summary", F["tiny"], COLORS["muted"], "ra")


def build() -> None:
    data = load_data()
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 44 | ACT 5 | PILLAR MAP", F["kicker"], COLORS["violet"])
    bx = 2045
    bx = badge(draw, bx, 56, "V8", "4 pillars", COLORS["violet"])
    bx = badge(draw, bx, 56, "BRIDGE", f"AUROC {data['bridge']['rf']['cv_mean_auroc_aug']:.3f}", COLORS["blue"])
    bx = badge(draw, bx, 56, "DECOMPOSE", data["interaction_range"], COLORS["teal"])
    badge(draw, bx, 56, "CAUSAL", f"{data['causal']['n_environments']} env", COLORS["amber"])

    text(draw, (120, 330), "v8 Has Four Incubator Pillars", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "Each pillar turns the established benchmark lane into one translational analysis question.",
        F["subtitle"],
        3000,
        COLORS["muted"],
        8,
    )

    draw_pillars(draw, data)
    draw_contract_strip(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    ImageOps.colorize(gray, black="#080D14", white="#F3F7FC").save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(gray)

    manifest = {
        "title": "v8 Has Four Incubator Pillars",
        "readout": "BRIDGE, DECOMPOSE, INTERVENE, and CAUSAL convert benchmarked evidence into four follow-up analysis questions.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "source_manifests": {
            "v8_transition": str(V8_TRANSITION_MANIFEST.relative_to(ROOT)),
        },
        "source_files": {
            "v8_summary": str(V8_SUMMARY.relative_to(ROOT)),
            "bridge": str(V8_BRIDGE.relative_to(ROOT)),
            "species": str(V8_SPECIES.relative_to(ROOT)),
            "tissue_nes": str(V8_TISSUE_NES.relative_to(ROOT)),
            "factorial": str(V8_FACTORIAL.relative_to(ROOT)),
            "mars": str(V8_MARS.relative_to(ROOT)),
            "offline": str(V8_OFFLINE.relative_to(ROOT)),
            "multi_tissue": str(V8_MULTI_TISSUE.relative_to(ROOT)),
            "pareto": str(V8_PARETO.relative_to(ROOT)),
            "crispr": str(V8_CRISPR.relative_to(ROOT)),
            "causal": str(V8_CAUSAL.relative_to(ROOT)),
            "dag_edges": str(V8_DAG_EDGES.relative_to(ROOT)),
        },
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": {
            "v8_pillars": 4,
            "bridge_rf_auroc_aug": data["bridge"]["rf"]["cv_mean_auroc_aug"],
            "bridge_rf_delta": data["bridge"]["rf"]["delta_aug_minus_base"]["mean"],
            "bridge_pathways": data["bridge"]["n_pathways"],
            "species_comparisons": data["species"]["n_comparisons"],
            "interaction_range": data["interaction_range"],
            "mars_ir_dose": data["mars"]["linear_dose"]["IR"],
            "intervene_candidates": data["multi_tissue"]["n_drugs_total"],
            "intervene_multi_tissue": data["multi_tissue"]["n_drugs_multi_tissue"],
            "pareto_axes": data["pareto_count"],
            "causal_environments": data["causal"]["n_environments"],
            "causal_genes": data["causal"]["n_genes"],
            "causal_t_icp_mean": data["causal"]["stressor_aggregate_icp"]["T"]["mean"],
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    QA_NOTE.write_text(
        "\n".join(
            [
                "# V8 Pillars Asset Visual QA",
                "",
                "Slide 44 is designed as the v8 pillar map before the result-detail slides.",
                "",
                "Checks to review:",
                "- Four pillar cards have equal hierarchy and enough gutters.",
                "- Mini visuals remain legible in thumbnail and full-size views.",
                "- Pillar row labels align across all four cards.",
                "- Bottom reading strip does not collide with pillar cards.",
                "- Grayscale version preserves the four-module distinction.",
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

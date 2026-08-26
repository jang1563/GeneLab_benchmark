#!/usr/bin/env python3
"""Build slide 48 asset: causal maps connect stressor axes to tissue readouts."""

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
OUT_DIR = ASSET_ROOT / "causal_map"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "causal_maps_connect_stressor_axes_to_tissue_readouts_premium.png"
GRAY_PATH = OUT_DIR / "causal_maps_connect_stressor_axes_to_tissue_readouts_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "causal_map_manifest.json"
QA_NOTE = OUT_DIR / "CAUSAL_MAP_ASSET_VISUAL_QA.md"

ICP = ROOT / "v8" / "causal" / "evaluation" / "icp_dag_summary.json"
EDGES = ROOT / "v8" / "causal" / "evaluation" / "dag_edges.csv"
SCORES = ROOT / "v8" / "causal" / "evaluation" / "icp_stressor_pathway_scores.csv"

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


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_edges() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with EDGES.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "from": row["node_from"].replace("stressor:", ""),
                    "to": row["node_to"].replace("tissue:", ""),
                    "icp_score": float(row["icp_score"]),
                    "n_genes": int(row["n_genes"]),
                    "evidence": row["evidence"],
                }
            )
    rows.sort(key=lambda item: float(item["icp_score"]), reverse=True)
    return rows


def load_data() -> dict[str, object]:
    summary = load_json(ICP)
    edges = load_edges()
    return {"summary": summary, "edges": edges}


def stressor_label(value: str) -> str:
    return {
        "HLU": "microgravity",
        "IR": "radiation",
        "HLUxIR": "combined stress",
        "T": "time / isolation",
        "IRxT": "radiation x time",
        "HLUxIRxT": "three-way stress",
        "HLUxT": "microgravity x time",
        "combined_flight": "pooled flight",
    }.get(value, value)


def stressor_color(value: str) -> str:
    if value.startswith("IR"):
        return COLORS["amber"]
    if value.startswith("HLU"):
        return COLORS["blue"] if "x" not in value else COLORS["rose"]
    if value == "T":
        return COLORS["teal"]
    return COLORS["dim"]


def draw_method_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x1, y1, x2, y2 = 120, 610, 1210, 1462
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 36), "What CAUSAL does", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 92),
        "ICP scores how consistently a stressor-gene signal appears across analog and flight settings, then rolls stable genes into tissue edges.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )
    summary = data["summary"]
    cards = [
        ("settings", str(summary["n_environments"]), COLORS["blue"]),
        ("genes", f"{summary['n_genes']:,}", COLORS["teal"]),
        ("edge rows", str(len(data["edges"])), COLORS["rose"]),
    ]
    cx, cy = x1 + 72, y1 + 222
    for i, (label, value, color) in enumerate(cards):
        xx = cx + i * 300
        rounded(draw, (xx, cy, xx + 246, cy + 112), 24, COLORS["panel2"], color, 2)
        text(draw, (xx + 24, cy + 22), label, F["micro_bold"], COLORS["dim"])
        text(draw, (xx + 24, cy + 78), value, F["metric2"], COLORS["text"], "lm")

    flow = [
        ("stressor axis", COLORS["blue"]),
        ("gene stability", COLORS["teal"]),
        ("tissue edge", COLORS["rose"]),
    ]
    fy = y1 + 438
    for i, (label, color) in enumerate(flow):
        xx = x1 + 90 + i * 306
        rounded(draw, (xx, fy, xx + 236, fy + 88), 22, COLORS["panel2"], color, 2)
        text(draw, (xx + 118, fy + 44), label, F["tiny_bold"], COLORS["text"], "mm")
        if i < len(flow) - 1:
            arrow(draw, xx + 248, fy + 44, xx + 292, fy + 44, COLORS["dim"], 4)

    rounded(draw, (x1 + 64, y2 - 230, x2 - 64, y2 - 42), 26, "#122234", COLORS["green"], 2)
    text(draw, (x1 + 96, y2 - 194), "ICP score", F["small_bold"], COLORS["green"])
    paragraph(
        draw,
        (x1 + 96, y2 - 156),
        "Higher ICP means the stressor signal stays more consistent across environments and becomes a heavier map edge.",
        F["small"],
        x2 - x1 - 192,
        COLORS["muted"],
        7,
    )


def draw_edge_network(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x1, y1, x2, y2 = 1265, 610, 2575, 1462
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 36), "Top stressor-to-tissue edges", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 92),
        "The ranked DAG concentrates stressor stability into a small set of tissue readouts.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )

    stressors = {
        "IR": (COLORS["amber"], "radiation"),
        "HLU": (COLORS["blue"], "microgravity"),
        "HLUxIR": (COLORS["rose"], "combined"),
    }
    top_edges = [e for e in data["edges"] if e["from"] in stressors][:5]

    rank_x = x1 + 62
    stressor_x, stressor_w = x1 + 118, 230
    score_x, score_w = x1 + 460, 360
    tissue_x, tissue_w = x2 - 310, 230
    row_y = y1 + 218
    row_gap = 104
    row_h = 74
    max_score = max(float(edge["icp_score"]) for edge in top_edges)

    text(draw, (stressor_x, y1 + 176), "stressor", F["micro_bold"], COLORS["dim"])
    text(draw, (score_x, y1 + 176), "ranked edge", F["micro_bold"], COLORS["dim"])
    text(draw, (tissue_x, y1 + 176), "tissue", F["micro_bold"], COLORS["dim"])

    for i, edge in enumerate(top_edges):
        color, stressor_desc = stressors[edge["from"]]
        y = row_y + i * row_gap
        mid = y + row_h // 2
        tissue_label = "gastroc." if edge["to"] == "gastrocnemius" else str(edge["to"])
        score = float(edge["icp_score"])

        draw.line((stressor_x + stressor_w + 18, mid, score_x - 18, mid), fill=rgba(color, 150), width=4)
        draw.polygon([(score_x - 18, mid), (score_x - 34, mid - 9), (score_x - 34, mid + 9)], fill=rgba(color, 190))
        draw.line((score_x + score_w + 18, mid, tissue_x - 18, mid), fill=rgba(color, 150), width=4)
        draw.polygon([(tissue_x - 18, mid), (tissue_x - 34, mid - 9), (tissue_x - 34, mid + 9)], fill=rgba(color, 190))

        rounded(draw, (rank_x - 26, mid - 26, rank_x + 26, mid + 26), 18, COLORS["panel3"], color, 2)
        text(draw, (rank_x, mid), str(i + 1), F["tiny_bold"], COLORS["text"], "mm")

        rounded(draw, (stressor_x, y, stressor_x + stressor_w, y + row_h), 20, COLORS["panel2"], color, 2)
        text(draw, (stressor_x + 22, y + 17), str(edge["from"]), F["small_bold"], COLORS["text"])
        text(draw, (stressor_x + 22, y + 49), stressor_desc, F["micro_bold"], COLORS["muted"])

        rounded(draw, (score_x, y, score_x + score_w, y + row_h), 20, COLORS["panel3"], color, 2)
        text(draw, (score_x + 22, y + 15), f"{edge['from']} -> {tissue_label}", F["tiny_bold"], COLORS["text"])
        text(draw, (score_x + score_w - 22, y + 15), f"{score:.3f}", F["tiny_bold"], color, "ra")
        bar_x, bar_y = score_x + 22, y + 51
        bar_w = score_w - 44
        rounded(draw, (bar_x, bar_y, bar_x + bar_w, bar_y + 8), 4, "#223046")
        rounded(draw, (bar_x, bar_y, bar_x + int(bar_w * score / max_score), bar_y + 8), 4, color)

        rounded(draw, (tissue_x, y, tissue_x + tissue_w, y + row_h), 20, COLORS["panel2"], COLORS["teal"], 2)
        text(draw, (tissue_x + tissue_w // 2, mid), tissue_label, F["tiny_bold"], COLORS["text"], "mm")

    rounded(draw, (x1 + 70, y2 - 110, x2 - 70, y2 - 42), 22, "#122234", COLORS["teal"], 2)
    paragraph(
        draw,
        (x1 + 104, y2 - 92),
        "Top edges highlight radiation-to-skin, microgravity-to-thymus, and kidney combined-stress links.",
        F["small"],
        x2 - x1 - 208,
        COLORS["muted"],
        6,
    )


def draw_stability_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x1, y1, x2, y2 = 2630, 610, 3720, 1462
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 36), "Stressor stability tiers", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 92),
        "Aggregate ICP ranks the stressor axes before tissue edges are drawn.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )

    order = ["T", "IRxT", "HLUxIRxT", "HLUxT", "HLU", "IR", "HLUxIR", "combined_flight"]
    lookup = data["summary"]["stressor_aggregate_icp"]
    chart_x, chart_y = x1 + 310, y1 + 218
    bar_w, bar_h = 570, 38
    max_v = 0.56
    for i, key in enumerate(order):
        value = float(lookup[key]["mean"])
        count = int(lookup[key]["count"])
        y = chart_y + i * 58
        color = stressor_color(key)
        label = stressor_label(key)
        text(draw, (x1 + 70, y + 6), label, F["tiny_bold"], COLORS["text"])
        rounded(draw, (chart_x, y, chart_x + bar_w, y + bar_h), 13, COLORS["panel3"], "#2A394D", 1)
        fill_w = int(bar_w * value / max_v)
        rounded(draw, (chart_x, y, chart_x + fill_w, y + bar_h), 13, color)
        text(draw, (chart_x + fill_w + 14, y + 6), f"{value:.3f}", F["tiny_bold"], COLORS["muted"])
        text(draw, (x2 - 70, y + 7), f"n={count:,}", F["micro_bold"], COLORS["dim"], "ra")

    rounded(draw, (x1 + 70, y2 - 188, x1 + 420, y2 - 44), 26, COLORS["panel2"], COLORS["teal"], 2)
    text(draw, (x1 + 104, y2 - 158), "top axis", F["micro_bold"], COLORS["dim"])
    text(draw, (x1 + 104, y2 - 118), "T 0.540", F["section"], COLORS["text"])
    text(draw, (x1 + 104, y2 - 68), "time / isolation", F["tiny_bold"], COLORS["muted"])

    top_gene = next(iter(data["summary"]["top20_invariant_genes"].items()))
    rounded(draw, (x1 + 455, y2 - 188, x2 - 70, y2 - 44), 26, "#122234", COLORS["violet"], 2)
    text(draw, (x1 + 488, y2 - 158), "top stable gene", F["micro_bold"], COLORS["dim"])
    text(draw, (x1 + 488, y2 - 118), top_gene[0], F["small_bold"], COLORS["text"])
    text(draw, (x1 + 488, y2 - 80), f"ICP {top_gene[1]:.3f}", F["tiny_bold"], COLORS["muted"])


def draw_flow_strip(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1530, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "CAUSAL Integration Flow", F["h2"], COLORS["text"])
    steps = [
        ("Stressor axes", "microgravity, radiation, time", COLORS["blue"]),
        ("ICP stability", "gene-level consistency", COLORS["teal"]),
        ("Tissue readouts", "edge-weighted response map", COLORS["rose"]),
        ("Human biology", "immune, muscle, organ axes", COLORS["violet"]),
        ("Perturbation axes", "CDK, HSP90, AMPK/BMP", COLORS["green"]),
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
    text(draw, (x2 - 80, y1 + 72), "CAUSAL", F["metric2"], COLORS["teal"], "ra")
    text(draw, (x2 - 80, y1 + 124), "ranked stressor-to-tissue map", F["small_bold"], COLORS["muted"], "ra")


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "CAUSAL condenses stressor evidence into a ranked map: time-related axes lead aggregate stability, while skin, thymus, and kidney anchor the strongest tissue edges.",
        F["small"],
        3180,
        COLORS["muted"],
        8,
    )


def build() -> None:
    data = load_data()
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 48 | ACT 5 | CAUSAL", F["kicker"], COLORS["teal"])
    bx = 1880
    bx = badge(draw, bx, 56, "SETTINGS", str(data["summary"]["n_environments"]), COLORS["blue"])
    bx = badge(draw, bx, 56, "GENES", f"{data['summary']['n_genes']:,}", COLORS["teal"])
    bx = badge(draw, bx, 56, "DAG EDGES", str(len(data["edges"])), COLORS["rose"])
    badge(draw, bx, 56, "TOP ICP", "T 0.540", COLORS["amber"])

    text(draw, (120, 330), "Causal Maps Connect Stressor Axes To Tissue Readouts", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "CAUSAL converts ICP stability into a ranked map from stressor axes to tissue-level biology.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_method_panel(draw, data)
    draw_edge_network(draw, data)
    draw_stability_panel(draw, data)
    draw_flow_strip(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    ImageOps.colorize(gray, black="#080D14", white="#F3F7FC").save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(gray)

    manifest = {
        "title": "Causal Maps Connect Stressor Axes To Tissue Readouts",
        "claim": "CAUSAL converts ICP stability into a ranked stressor-to-tissue map.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "source_files": {
            "icp_summary": str(ICP.relative_to(ROOT)),
            "dag_edges": str(EDGES.relative_to(ROOT)),
            "stressor_pathway_scores": str(SCORES.relative_to(ROOT)),
        },
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": {
            "n_environments": data["summary"]["n_environments"],
            "n_genes": data["summary"]["n_genes"],
            "n_edges": len(data["edges"]),
            "top_edge": data["edges"][0],
            "top_aggregate_icp": data["summary"]["stressor_aggregate_icp"]["T"]["mean"],
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    QA_NOTE.write_text(
        "\n".join(
            [
                "# CAUSAL Map Asset Visual QA",
                "",
                "Slide 48 rebuilds the v8 causal DAG into a presentation-scale ranked map.",
                "",
                "Checks to review:",
                "- Method panel explains ICP edge weighting without crowding.",
                "- Network panel keeps stressor nodes, tissue nodes, and edge labels separated.",
                "- Stability bars preserve row labels, values, and n counts.",
                "- Bottom CAUSAL integration flow has stable arrow spacing.",
                "- Grayscale version preserves edge and bar hierarchy.",
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

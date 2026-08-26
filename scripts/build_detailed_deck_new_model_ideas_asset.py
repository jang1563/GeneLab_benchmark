#!/usr/bin/env python3
"""Build the detailed-deck v7 new-model-ideas comparison asset."""

from __future__ import annotations

import json
import statistics
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
V4_M1 = ROOT / "v4" / "evaluation" / "M1_summary.json"
V7_EVAL = ROOT / "v7" / "evaluation"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "new_model_ideas"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "newer_model_ideas_preserve_the_lesson_premium.png"
GRAY_PATH = OUT_DIR / "newer_model_ideas_preserve_the_lesson_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "new_model_ideas_manifest.json"

TISSUES = ["liver", "gastrocnemius", "kidney", "thymus", "eye", "skin"]
TISSUE_LABELS = {
    "liver": "Liver",
    "gastrocnemius": "Gastro",
    "kidney": "Kidney",
    "thymus": "Thymus",
    "eye": "Eye",
    "skin": "Skin",
}
GRAPH_LABELS = {"wgcna": "WGCNA", "random": "Random", "no_edges": "No-edge"}

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "panel2": "#151F2D",
    "panel3": "#121B28",
    "grid": "#2A3546",
    "axis": "#98A7BA",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#667487",
    "blue": "#66A6E8",
    "navy": "#2D5D9A",
    "sky": "#73A7FF",
    "teal": "#5FD3C4",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "orange": "#E69F00",
    "rose": "#E17882",
    "violet": "#B39DFF",
    "white": "#FFFFFF",
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


def delta_color(delta: float) -> str:
    if delta >= 0:
        return blend("#203547", COLORS["green"], min(1.0, delta / 0.08))
    return blend("#253044", COLORS["rose"], min(1.0, abs(delta) / 0.30))


def score_color(score: float) -> str:
    if score < 0.50:
        return blend("#253044", COLORS["rose"], (0.50 - score) / 0.20)
    if score < 0.70:
        return blend("#253044", COLORS["teal"], (score - 0.50) / 0.20)
    return blend(COLORS["teal"], COLORS["amber"], min(1.0, (score - 0.70) / 0.25))


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
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
    "title": load_font(76, True),
    "subtitle": load_font(35),
    "h2": load_font(42, True),
    "h3": load_font(32, True),
    "body": load_font(29),
    "body_bold": load_font(29, True),
    "small": load_font(25),
    "small_bold": load_font(25, True),
    "tiny": load_font(21),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "stat": load_font(66, True),
    "huge": load_font(90, True),
}


def rounded(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    radius: int,
    fill: str | tuple[int, int, int, int],
    outline: str | tuple[int, int, int, int] | None = None,
    width: int = 1,
) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int | float, int | float],
    value: str,
    font: ImageFont.ImageFont,
    fill: str | tuple[int, int, int] | tuple[int, int, int, int] = COLORS["text"],
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
        else:
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
    fill: str | tuple[int, int, int, int],
    leading: int = 8,
) -> int:
    x, y = xy
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def draw_background(canvas: Image.Image) -> None:
    overlay = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba("#1E2A3A", 28), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba("#1E2A3A", 22), width=1)
    draw.rectangle((0, 0, W, 310), fill=rgba("#121B28", 170))
    draw.rectangle((0, 1855, W, H), fill=rgba("#090E15", 145))
    canvas.alpha_composite(overlay)


def draw_badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(176, int(draw.textlength(value, font=F["tiny_bold"]) + 56))
    rounded(draw, (x, y, x + width, y + 62), 18, rgba("#1B2838", 235), rgba(color, 165), 2)
    text(draw, (x + 18, y + 12), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 18, y + 34), value, F["tiny_bold"], COLORS["text"])
    return x + width + 16


def draw_panel_header(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], label: str, title: str, subtitle: str) -> None:
    x1, y1, x2, _ = box
    text(draw, (x1 + 44, y1 + 34), label.upper(), F["tiny_bold"], COLORS["teal"])
    text(draw, (x1 + 44, y1 + 66), title, F["h2"], COLORS["text"])
    paragraph(draw, (x1 + 44, y1 + 118), subtitle, F["small"], x2 - x1 - 88, COLORS["muted"], 6)


def load_data() -> dict[str, object]:
    m1 = json.loads(V4_M1.read_text())
    pca = {t: float(m1[t]["gene"]["pca_lr"]["auroc"]) for t in TISSUES}
    scprint = {}
    gnn = {}
    for tissue in TISSUES:
        scprint[tissue] = json.loads((V7_EVAL / f"SCPRINT2_{tissue}.json").read_text())
        gnn[tissue] = {}
        for graph in GRAPH_LABELS:
            gnn[tissue][graph] = json.loads((V7_EVAL / f"GNN_{tissue}_{graph}.json").read_text())

    sc_values = [float(scprint[t]["mean_auroc"]) for t in TISSUES]
    pca_values = [pca[t] for t in TISSUES]
    wgcna_values = [float(gnn[t]["wgcna"]["mean_auroc"]) for t in TISSUES]
    graph_best = {
        t: max(GRAPH_LABELS, key=lambda graph: float(gnn[t][graph]["mean_auroc"]))
        for t in TISSUES
    }
    graph_best_counts = {graph: list(graph_best.values()).count(graph) for graph in GRAPH_LABELS}
    return {
        "pca": pca,
        "scprint": scprint,
        "gnn": gnn,
        "pca_mean": statistics.mean(pca_values),
        "sc_mean": statistics.mean(sc_values),
        "wgcna_mean": statistics.mean(wgcna_values),
        "sc_below": sum(1 for t in TISSUES if float(scprint[t]["delta_vs_pca_lr"]) < 0),
        "sc_sig": sum(1 for t in TISSUES if float(scprint[t]["perm_pvalue"]) < 0.05),
        "wgcna_above": sum(1 for t in TISSUES if float(gnn[t]["wgcna"]["delta_vs_pca_lr"]) > 0),
        "wgcna_sig": sum(1 for t in TISSUES if float(gnn[t]["wgcna"]["perm_pvalue"]) < 0.05),
        "graph_best": graph_best,
        "graph_best_counts": graph_best_counts,
    }


def draw_probe_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    number: str,
    title: str,
    color: str,
    steps: list[str],
    caption: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 26, rgba("#0C1420", 225), rgba(color, 150), 2)
    rounded(draw, (x1 + 30, y1 + 28, x1 + 90, y1 + 88), 18, rgba(color, 58), rgba(color, 180), 2)
    text(draw, (x1 + 60, y1 + 45), number, F["h3"], COLORS["text"], "mm")
    text(draw, (x1 + 112, y1 + 34), title, F["h3"], COLORS["text"])
    y = y1 + 118
    for idx, step in enumerate(steps):
        rounded(draw, (x1 + 32, y, x2 - 32, y + 58), 16, rgba("#152130", 230), rgba("#2A3546", 220), 1)
        text(draw, (x1 + 54, y + 18), step, F["small_bold" if idx == len(steps) - 1 else "small"], COLORS["text" if idx == len(steps) - 1 else "muted"])
        if idx < len(steps) - 1:
            draw.line((x1 + 72, y + 61, x1 + 72, y + 83), fill=rgba(color, 190), width=3)
            draw.polygon([(x1 + 64, y + 80), (x1 + 80, y + 80), (x1 + 72, y + 92)], fill=rgba(color, 220))
        y += 92
    paragraph(draw, (x1 + 34, y2 - 116), caption, F["tiny"], x2 - x1 - 68, COLORS["muted"], 6)


def draw_probe_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    box = (120, 330, 1075, 1840)
    rounded(draw, box, 30, COLORS["panel"], rgba("#2A394D", 235), 2)
    draw_panel_header(
        draw,
        box,
        "A. What v7 adds",
        "Two extra probes",
        "Same mission-held-out contract; only the model surface changes.",
    )
    draw_probe_card(
        draw,
        (165, 520, 1030, 1010),
        "1",
        "scPRINT-2 embeddings",
        COLORS["orange"],
        ["Bulk matrix -> AnnData", "scPRINT vocabulary mapping", "256D embedding + classifier"],
        "Question: does a cell foundation embedding improve the bulk task surface?",
    )
    draw_probe_card(
        draw,
        (165, 1060, 1030, 1550),
        "2",
        "GNN topology test",
        COLORS["green"],
        ["Train-fold WGCNA graph", "Random / no-edge controls", "LOMO AUROC comparison"],
        "Question: does graph structure add stable transfer signal beyond tabular ML?",
    )
    rounded(draw, (165, 1600, 1030, 1785), 24, rgba("#111B28", 240), rgba(COLORS["teal"], 110), 2)
    text(draw, (200, 1631), "Readout rule", F["h3"], COLORS["text"])
    paragraph(
        draw,
        (200, 1681),
        "A new architecture earns weight when it wins under the same split, label, and AUROC contract.",
        F["body"],
        770,
        COLORS["muted"],
        8,
    )


def draw_scprint_comparison(draw: ImageDraw.ImageDraw, data: dict[str, object], box: tuple[int, int, int, int]) -> None:
    x1, y1, x2, y2 = box
    pca = data["pca"]
    scprint = data["scprint"]
    text(draw, (x1, y1), "scPRINT-2 vs paired PCA-LR", F["h3"], COLORS["text"])
    text(draw, (x2, y1 + 6), "AUROC", F["tiny_bold"], COLORS["dim"], "ra")
    plot_x1, plot_y1, plot_x2, plot_y2 = x1 + 170, y1 + 72, x2 - 40, y2 - 25
    draw.line((plot_x1, plot_y2, plot_x2, plot_y2), fill=rgba(COLORS["axis"], 130), width=2)
    for tick in [0.30, 0.50, 0.70, 0.90]:
        tx = plot_x1 + (tick - 0.25) / 0.70 * (plot_x2 - plot_x1)
        draw.line((tx, plot_y1, tx, plot_y2), fill=rgba("#2A3546", 145), width=1)
        text(draw, (tx, plot_y2 + 14), f"{tick:.1f}", F["micro"], COLORS["dim"], "ma")
    row_h = (plot_y2 - plot_y1) / len(TISSUES)
    for i, tissue in enumerate(TISSUES):
        y = plot_y1 + row_h * i + row_h * 0.50
        text(draw, (x1 + 6, y), TISSUE_LABELS[tissue], F["small_bold"], COLORS["text"], "lm")
        sc = float(scprint[tissue]["mean_auroc"])
        ref = float(pca[tissue])
        sx = plot_x1 + (sc - 0.25) / 0.70 * (plot_x2 - plot_x1)
        rx = plot_x1 + (ref - 0.25) / 0.70 * (plot_x2 - plot_x1)
        draw.line((sx, y, rx, y), fill=rgba("#8392A7", 145), width=4)
        draw.ellipse((rx - 13, y - 13, rx + 13, y + 13), fill=COLORS["blue"], outline=COLORS["white"], width=1)
        draw.ellipse((sx - 13, y - 13, sx + 13, y + 13), fill=COLORS["orange"], outline=COLORS["white"], width=1)
        text(draw, (rx + 20, y - 18), f"{ref:.2f}", F["micro_bold"], COLORS["blue"])
        text(draw, (sx - 20, y + 12), f"{sc:.2f}", F["micro_bold"], COLORS["orange"], "ra")
    lx = x1 + 8
    ly = y2 - 6
    draw.ellipse((lx, ly - 13, lx + 24, ly + 11), fill=COLORS["blue"])
    text(draw, (lx + 34, ly - 13), "PCA-LR", F["tiny"], COLORS["muted"])
    draw.ellipse((lx + 134, ly - 13, lx + 158, ly + 11), fill=COLORS["orange"])
    text(draw, (lx + 168, ly - 13), "scPRINT-2", F["tiny"], COLORS["muted"])


def draw_gnn_heatmap(draw: ImageDraw.ImageDraw, data: dict[str, object], box: tuple[int, int, int, int]) -> None:
    x1, y1, x2, y2 = box
    gnn = data["gnn"]
    graph_best = data["graph_best"]
    text(draw, (x1, y1), "GNN graph variants vs PCA-LR", F["h3"], COLORS["text"])
    text(draw, (x2, y1 + 6), "delta AUROC", F["tiny_bold"], COLORS["dim"], "ra")
    grid_x1, grid_y1 = x1 + 170, y1 + 82
    cell_w = (x2 - grid_x1 - 18) / len(GRAPH_LABELS)
    cell_h = 68
    for j, graph in enumerate(GRAPH_LABELS):
        cx = grid_x1 + j * cell_w + cell_w * 0.5
        text(draw, (cx, y1 + 48), GRAPH_LABELS[graph], F["tiny_bold"], COLORS["muted"], "ma")
    for i, tissue in enumerate(TISSUES):
        y = grid_y1 + i * (cell_h + 14)
        text(draw, (x1 + 6, y + cell_h / 2), TISSUE_LABELS[tissue], F["small_bold"], COLORS["text"], "lm")
        for j, graph in enumerate(GRAPH_LABELS):
            x = grid_x1 + j * cell_w
            entry = gnn[tissue][graph]
            delta = float(entry["delta_vs_pca_lr"])
            value = float(entry["mean_auroc"])
            fill = delta_color(delta)
            outline = COLORS["white"] if graph_best[tissue] == graph else "#2A3546"
            rounded(draw, (int(x), y, int(x + cell_w - 14), y + cell_h), 14, fill, outline, 2 if graph_best[tissue] == graph else 1)
            text(draw, (x + cell_w * 0.50 - 7, y + 17), f"{delta:+.02f}", F["tiny_bold"], COLORS["text"], "ma")
            text(draw, (x + cell_w * 0.50 - 7, y + 42), f"{value:.2f}", F["micro"], rgba(COLORS["white"], 210), "ma")
    y = y2 - 46
    rounded(draw, (x1, y, x1 + 520, y + 40), 14, blend("#101823", COLORS["green"], 0.18), COLORS["green"], 1)
    text(draw, (x1 + 18, y + 11), "Positive cells are local; topology wins are mixed", F["tiny_bold"], COLORS["muted"])


def draw_evidence_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    box = (1065, 330, 2535, 1840)
    rounded(draw, box, 30, COLORS["panel"], rgba("#2A394D", 235), 2)
    draw_panel_header(
        draw,
        box,
        "B. Same benchmark surface",
        "Paired comparisons keep the result readable",
        "Each row asks whether the new model surface beats the paired tabular baseline.",
    )
    draw_scprint_comparison(draw, data, (1135, 535, 2465, 1050))
    draw.line((1135, 1100, 2465, 1100), fill=rgba("#2A3546", 190), width=2)
    draw_gnn_heatmap(draw, data, (1135, 1148, 2465, 1775))


def stat_tile(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    value: str,
    label: str,
    color: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 22, rgba("#111B28", 240), rgba(color, 120), 2)
    text(draw, (x1 + 26, y1 + 24), value, F["stat"], COLORS["text"])
    paragraph(draw, (x1 + 28, y1 + 102), label, F["small"], x2 - x1 - 56, COLORS["muted"], 6)


def draw_decision_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    box = (2525, 330, 3730, 1840)
    rounded(draw, box, 30, COLORS["panel"], rgba("#2A394D", 235), 2)
    draw_panel_header(
        draw,
        box,
        "C. Decision readout",
        "Task fit remains the bottleneck",
        "The strongest readout remains task fit under mission-held-out scoring.",
    )
    x1, y1, x2, _ = box
    stat_tile(draw, (2585, 535, 2915, 725), "6/6", "scPRINT-2 rows below paired PCA-LR", COLORS["orange"])
    stat_tile(draw, (2945, 535, 3275, 725), "0", "scPRINT-2 local wins vs reference", COLORS["rose"])
    stat_tile(draw, (3305, 535, 3670, 725), "1/6", "WGCNA-GNN rows above PCA-LR", COLORS["green"])

    rounded(draw, (2585, 770, 3670, 1080), 24, rgba("#101723", 245), rgba("#2A3546", 230), 1)
    text(draw, (2625, 805), "Mean AUROC across six tissues", F["h3"], COLORS["text"])
    means = [
        ("PCA-LR", float(data["pca_mean"]), COLORS["blue"]),
        ("WGCNA-GNN", float(data["wgcna_mean"]), COLORS["green"]),
        ("scPRINT-2", float(data["sc_mean"]), COLORS["orange"]),
    ]
    max_v = 0.82
    y = 870
    for label, value, color in means:
        text(draw, (2625, y + 9), label, F["small_bold"], COLORS["text"])
        bar_x1, bar_x2 = 2830, 3575
        rounded(draw, (bar_x1, y, bar_x2, y + 36), 13, rgba("#1B2838", 220), None, 1)
        rounded(draw, (bar_x1, y, int(bar_x1 + (value / max_v) * (bar_x2 - bar_x1)), y + 36), 13, color, None, 1)
        text(draw, (3615, y + 7), f"{value:.3f}", F["small_bold"], COLORS["text"])
        y += 66

    rounded(draw, (2585, 1120, 3670, 1425), 24, rgba("#101723", 245), rgba("#2A3546", 230), 1)
    text(draw, (2625, 1155), "Graph structure check", F["h3"], COLORS["text"])
    counts = data["graph_best_counts"]
    chips = [
        ("WGCNA best", counts["wgcna"], COLORS["green"]),
        ("Random best", counts["random"], COLORS["sky"]),
        ("No-edge best", counts["no_edges"], COLORS["violet"]),
    ]
    x = 2625
    for label, count, color in chips:
        rounded(draw, (x, 1224, x + 305, 1318), 22, blend("#101823", color, 0.20), color, 2)
        text(draw, (x + 24, 1245), f"{count}/6", F["h2"], COLORS["text"])
        text(draw, (x + 24, 1290), label, F["tiny_bold"], COLORS["muted"])
        x += 340
    paragraph(
        draw,
        (2625, 1350),
        "Observed wins are mixed across WGCNA, random, and no-edge variants.",
        F["small"],
        980,
        COLORS["muted"],
        7,
    )

    rounded(draw, (2585, 1470, 3670, 1785), 24, blend("#101823", COLORS["teal"], 0.18), COLORS["teal"], 2)
    text(draw, (2625, 1508), "Slide takeaway", F["h3"], COLORS["text"])
    paragraph(
        draw,
        (2625, 1560),
        "v7 adds credible newer-model probes, but the ranking logic stays the same: robust tabular baselines remain the benchmark floor.",
        F["body_bold"],
        970,
        COLORS["text"],
        9,
    )
    paragraph(
        draw,
        (2625, 1688),
        "Use these rows as extension evidence alongside the core benchmark readout.",
        F["small"],
        970,
        COLORS["muted"],
        7,
    )


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1888, 3730, 2076), 28, rgba("#0B1018", 225), rgba("#263448", 210), 2)
    text(draw, (164, 1926), "Takeaway", F["h3"], COLORS["text"])
    paragraph(
        draw,
        (164, 1979),
        "Paired rows keep newer-model probes anchored to the same mission-held-out AUROC frame as the preceding model slides.",
        F["small"],
        2250,
        COLORS["muted"],
        7,
    )
    text(draw, (3578, 1968), "v7", F["huge"], rgba(COLORS["teal"], 95), "mm")
    text(draw, (3578, 2033), "extension layer", F["tiny_bold"], COLORS["dim"], "mm")


def render() -> Image.Image:
    data = load_data()
    canvas = Image.new("RGBA", (W, H), COLORS["bg"])
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)

    text(draw, (120, 58), "SPACEBIO-BENCH / DETAILED DECK", F["kicker"], COLORS["teal"])
    text(draw, (120, 108), "Newer model ideas preserve the lesson", F["title"], COLORS["text"])
    paragraph(
        draw,
        (122, 205),
        "scPRINT-2 and WGCNA-GNN probes add breadth, while the benchmark contract keeps the comparison anchored.",
        F["subtitle"],
        2300,
        COLORS["muted"],
        6,
    )
    x = 2760
    x = draw_badge(draw, x, 76, "layer", "v7", COLORS["teal"])
    x = draw_badge(draw, x, 76, "tissues", "6", COLORS["blue"])
    x = draw_badge(draw, x, 76, "scPRINT", "6 rows", COLORS["orange"])
    draw_badge(draw, x, 76, "GNN", "18 rows", COLORS["green"])

    draw_probe_panel(draw, data)
    draw_evidence_panel(draw, data)
    draw_decision_panel(draw, data)
    draw_footer(draw)

    manifest = {
        "asset": str(OUT_PATH.relative_to(ROOT)),
        "title": "Newer model ideas preserve the lesson",
        "sources": [
            "v7/evaluation/SCPRINT2_*.json",
            "v7/evaluation/GNN_*_*.json",
            "v4/evaluation/M1_summary.json",
        ],
        "metrics": {
            "pca_lr_mean_auroc": round(float(data["pca_mean"]), 4),
            "scprint2_mean_auroc": round(float(data["sc_mean"]), 4),
            "wgcna_gnn_mean_auroc": round(float(data["wgcna_mean"]), 4),
            "scprint2_rows_below_pca_lr": int(data["sc_below"]),
            "scprint2_significant_rows": int(data["sc_sig"]),
            "wgcna_rows_above_pca_lr": int(data["wgcna_above"]),
            "wgcna_significant_rows": int(data["wgcna_sig"]),
            "gnn_best_graph_counts": data["graph_best_counts"],
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    return canvas


def main() -> None:
    canvas = render()
    rgb = canvas.convert("RGB")
    rgb.save(OUT_PATH, quality=95)
    rgb.convert("L").convert("RGB").save(GRAY_PATH, quality=95)
    print(json.dumps({"asset": str(OUT_PATH.relative_to(ROOT)), "grayscale": str(GRAY_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

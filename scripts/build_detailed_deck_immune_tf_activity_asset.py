#!/usr/bin/env python3
"""Build slide 38 asset: immune and TF activity prioritize follow-up."""

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
OUT_DIR = ASSET_ROOT / "immune_tf_activity"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "immune_tf_activity_prioritize_followup_premium.png"
GRAY_PATH = OUT_DIR / "immune_tf_activity_prioritize_followup_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "immune_tf_activity_manifest.json"

V5_EVAL = ROOT / "v5" / "evaluation"
IMMUNE_FILES = sorted(V5_EVAL.glob("immune_deconv_*.json"))
TF_FILES = sorted(V5_EVAL.glob("tf_activity_*.json"))
SIGNALING = V5_EVAL / "cross_organ_signaling.json"
SOURCE_FIG = ROOT / "v5" / "figures" / "html" / "Fig7_immune_signaling.html"

COLORS = {
    "bg": "#0C111A",
    "bg2": "#091019",
    "header": "#0F1824",
    "footer": "#080D14",
    "panel": "#101823",
    "panel2": "#151F2D",
    "panel3": "#0F1A26",
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
    "red": "#F17C88",
    "ink": "#080D14",
}

TISSUE_COLORS = {
    "skin": "#E17882",
    "thymus": "#B39DFF",
    "kidney": "#8BD17C",
    "liver": "#F4C26B",
    "colon": "#66A6E8",
    "eye": "#5FD3C4",
    "gastrocnemius": "#E69F00",
    "lung": "#A8B4C4",
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
    "h1": load_font(50, True),
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
    "stat": load_font(72, True),
    "stat2": load_font(62, True),
    "axis": load_font(20),
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
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(230, int(draw.textlength(value, font=F["tiny_bold"]) + 76))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def tissue_label(tissue: str) -> str:
    if tissue == "gastrocnemius":
        return "Gastro"
    return tissue.capitalize()


def load_data() -> dict[str, object]:
    immune_rows = []
    immune_by_tissue = {}
    for path in IMMUNE_FILES:
        item = json.loads(path.read_text())
        sig = []
        for name, result in item["cell_types"].items():
            fdr = float(result["fdr_p"])
            if fdr < 0.05:
                sig.append(
                    {
                        "cell_type": name,
                        "cliffs_delta": float(result["cliffs_delta"]),
                        "direction": result["direction"],
                        "fdr_p": fdr,
                    }
                )
        sig.sort(key=lambda row: abs(row["cliffs_delta"]), reverse=True)
        row = {
            "tissue": item["tissue"],
            "n_samples": int(item["n_samples"]),
            "n_cell_types": int(item["n_cell_types"]),
            "n_significant": int(item["n_significant_fdr05"]),
            "significant": sig,
        }
        immune_rows.append(row)
        immune_by_tissue[row["tissue"]] = row
    immune_rows.sort(key=lambda row: (-row["n_significant"], row["tissue"]))

    tf_rows = []
    tf_by_tissue = {}
    for path in TF_FILES:
        item = json.loads(path.read_text())
        sig = []
        for name, result in item["tf_results"].items():
            fdr = float(result["fdr_p"])
            if fdr < 0.05:
                sig.append(
                    {
                        "tf": name,
                        "t_stat": float(result["t_stat_approx"]),
                        "direction": result["direction"],
                        "fdr_p": fdr,
                    }
                )
        sig.sort(key=lambda row: abs(row["t_stat"]), reverse=True)
        row = {
            "tissue": item["tissue"],
            "n_samples": int(item["n_samples"]),
            "n_tfs_tested": int(item["n_tfs_tested"]),
            "n_significant": int(item["n_significant_fdr05"]),
            "significant": sig,
        }
        tf_rows.append(row)
        tf_by_tissue[row["tissue"]] = row
    tf_rows.sort(key=lambda row: (-row["n_significant"], row["tissue"]))

    signaling = json.loads(SIGNALING.read_text())
    active_pairs = []
    for pair in signaling["tissue_pairs"].values():
        for active in pair["active_pairs"]:
            active_pairs.append({"source": pair["source"], "target": pair["target"], **active})

    priority = []
    for tissue in sorted(set(immune_by_tissue) | set(tf_by_tissue)):
        immune = immune_by_tissue[tissue]
        tf = tf_by_tissue[tissue]
        score = immune["n_significant"] * 40 + tf["n_significant"]
        priority.append(
            {
                "tissue": tissue,
                "immune_calls": immune["n_significant"],
                "immune_total": immune["n_cell_types"],
                "tf_calls": tf["n_significant"],
                "tf_total": tf["n_tfs_tested"],
                "score": score,
            }
        )
    priority.sort(key=lambda row: (-row["score"], row["tissue"]))

    return {
        "immune": immune_rows,
        "immune_by_tissue": immune_by_tissue,
        "tf": tf_rows,
        "tf_by_tissue": tf_by_tissue,
        "priority": priority,
        "active_pairs": active_pairs,
        "signaling": signaling,
    }


def draw_background(draw: ImageDraw.ImageDraw) -> None:
    draw.rectangle((0, 0, W, H), fill=COLORS["bg"])
    for y in range(H):
        draw.line((0, y, W, y), fill=blend(COLORS["bg"], COLORS["bg2"], y / H))
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=COLORS["grid"], width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill="#172234", width=1)
    draw.rectangle((0, 0, W, 315), fill=COLORS["header"])
    draw.rectangle((0, 1840, W, H), fill=COLORS["footer"])
    draw.line((0, 315, W, 315), fill="#1E2B3D", width=2)
    draw.line((0, 1840, W, 1840), fill="#1E2B3D", width=2)


def draw_header(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 76), "BIOLOGY / IMMUNE + REGULATORY ACTIVITY", F["kicker"], COLORS["teal"])
    x = 2180
    x = badge(draw, x, 66, "layer", "v5 Fig7", COLORS["blue"])
    x = badge(draw, x, 66, "immune", "skin 6/14", COLORS["rose"])
    x = badge(draw, x, 66, "TF", "241 / 240", COLORS["violet"])
    badge(draw, x, 66, "tissues", "8", COLORS["amber"])


def draw_title(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 382), "Immune And TF Activity Prioritize Follow-Up", F["title"], COLORS["text"])
    paragraph(
        draw,
        (155, 493),
        "v5 separates cell-composition context from regulatory activity so follow-up tissue choices are visible.",
        F["subtitle"],
        2700,
        COLORS["muted"],
        10,
    )


def draw_panel_header(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, step: str, title: str, color: str) -> None:
    rounded(draw, (x, y, x + w, y + 1040), 34, COLORS["panel"], "#29374A", 2)
    rounded(draw, (x + 34, y + 34, x + 98, y + 98), 20, "#172335", color, 2)
    text(draw, (x + 66, y + 51), step, F["h3"], COLORS["text"], "ma")
    text(draw, (x + 120, y + 44), title, F["h2"], COLORS["text"])


def draw_ladder_row(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    label: str,
    value: int,
    total: int,
    max_value: int,
    color: str,
) -> None:
    text(draw, (x, y + 4), label, F["tiny_bold"], COLORS["text"])
    bar_x = x + 190
    bar_w = w - 320
    draw.line((bar_x, y + 24, bar_x + bar_w, y + 24), fill="#263244", width=16)
    fill_w = int(bar_w * max(0.025, value / max_value if max_value else 0))
    draw.line((bar_x, y + 24, bar_x + fill_w, y + 24), fill=color, width=16)
    draw.ellipse((bar_x + fill_w - 12, y + 12, bar_x + fill_w + 12, y + 36), fill=color)
    text(draw, (x + w - 58, y + 3), f"{value}/{total}", F["tiny_bold"], COLORS["text"], "ra")


def draw_chip(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str, width: int) -> None:
    rounded(draw, (x, y, x + width, y + 48), 18, "#172335", color, 2)
    text(draw, (x + 18, y + 12), label, F["tiny_bold"], COLORS["text"])
    text(draw, (x + width - 18, y + 12), value, F["tiny_bold"], color, "ra")


def draw_immune_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, data: dict[str, object]) -> None:
    color = COLORS["rose"]
    immune = data["immune"]
    skin = data["immune_by_tissue"]["skin"]
    draw_panel_header(draw, x, y, w, "1", "Immune deconvolution", color)
    text(draw, (x + 54, y + 145), "6/14", F["stat"], color)
    text(draw, (x + 270, y + 155), "skin cell types", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "Skin concentrates the clearest cell-composition readout, with kidney and thymus adding focused immune context.",
        F["body"],
        w - 108,
        COLORS["muted"],
        9,
    )

    rounded(draw, (x + 54, y + 420, x + w - 54, y + 748), 24, "#0D1520", "#2A394D", 1)
    text(draw, (x + 82, y + 446), "FDR<0.05 cell-type calls", F["small_bold"], COLORS["text"])
    max_sig = max(row["n_significant"] for row in immune)
    for i, row in enumerate(immune):
        draw_ladder_row(
            draw,
            x + 82,
            y + 498 + i * 32,
            w - 164,
            tissue_label(row["tissue"]),
            int(row["n_significant"]),
            int(row["n_cell_types"]),
            max_sig,
            TISSUE_COLORS[row["tissue"]],
        )

    rounded(draw, (x + 54, y + 800, x + w - 54, y + 964), 24, "#172335", "#2A394D", 1)
    text(draw, (x + 82, y + 825), "Skin call set", F["small_bold"], color)
    chip_w = 278
    for i, row in enumerate(skin["significant"][:6]):
        cx = x + 82 + (i % 3) * (chip_w + 18)
        cy = y + 866 + (i // 3) * 56
        direction_color = COLORS["rose"] if row["cliffs_delta"] > 0 else COLORS["blue"]
        name = str(row["cell_type"]).replace("Endothelial cells", "Endothelial").replace("Fibroblasts", "Fibroblast")
        draw_chip(draw, cx, cy, name[:18], f"{row['cliffs_delta']:+.2f}", direction_color, chip_w)


def draw_tf_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, data: dict[str, object]) -> None:
    color = COLORS["violet"]
    tf_rows = data["tf"]
    draw_panel_header(draw, x, y, w, "2", "TF activity", color)
    text(draw, (x + 54, y + 145), "241", F["stat"], color)
    text(draw, (x + 240, y + 155), "skin TFs", F["h3"], COLORS["muted"])
    text(draw, (x + 54, y + 220), "240 thymus TFs", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 282),
        "Regulatory activity creates a second ranking axis, strongest in skin and thymus with kidney and liver as secondary lanes.",
        F["body"],
        w - 108,
        COLORS["muted"],
        9,
    )

    rounded(draw, (x + 54, y + 420, x + w - 54, y + 748), 24, "#0D1520", "#2A394D", 1)
    text(draw, (x + 82, y + 446), "Significant TFs by tissue", F["small_bold"], COLORS["text"])
    max_sig = max(row["n_significant"] for row in tf_rows)
    for i, row in enumerate(tf_rows):
        draw_ladder_row(
            draw,
            x + 82,
            y + 498 + i * 32,
            w - 164,
            tissue_label(row["tissue"]),
            int(row["n_significant"]),
            int(row["n_tfs_tested"]),
            max_sig,
            TISSUE_COLORS[row["tissue"]],
        )

    rounded(draw, (x + 54, y + 800, x + w - 54, y + 964), 24, "#172335", "#2A394D", 1)
    text(draw, (x + 82, y + 825), "Top activity examples", F["small_bold"], color)
    example_tissues = ["skin", "thymus", "kidney", "liver"]
    for i, tissue in enumerate(example_tissues):
        top = [row["tf"] for row in data["tf_by_tissue"][tissue]["significant"][:4]]
        cy = y + 866 + i * 30
        text(draw, (x + 82, cy), tissue_label(tissue), F["tiny_bold"], TISSUE_COLORS[tissue])
        text(draw, (x + 230, cy), " / ".join(top), F["tiny"], COLORS["muted"])


def draw_scatter_axes(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    draw.line((x0, y1, x1, y1), fill="#566275", width=3)
    draw.line((x0, y0, x0, y1), fill="#566275", width=3)
    for tick in range(0, 7, 2):
        px = x0 + int(tick / 6 * (x1 - x0))
        draw.line((px, y1, px, y1 + 10), fill="#566275", width=2)
        text(draw, (px, y1 + 18), str(tick), F["axis"], COLORS["dim"], "ma")
    for tick in [0, 100, 200]:
        py = y1 - int(tick / 250 * (y1 - y0))
        draw.line((x0 - 10, py, x0, py), fill="#566275", width=2)
        text(draw, (x0 - 18, py - 10), str(tick), F["axis"], COLORS["dim"], "ra")
        draw.line((x0, py, x1, py), fill=rgba("#566275", 45), width=1)
    text(draw, ((x0 + x1) / 2, y1 + 34), "immune cell-type calls", F["tiny_bold"], COLORS["muted"], "ma")
    text(draw, (x0 - 48, y0 - 28), "TF calls", F["tiny_bold"], COLORS["muted"])


def scatter_point(box: tuple[int, int, int, int], immune_calls: int, tf_calls: int) -> tuple[int, int]:
    x0, y0, x1, y1 = box
    px = x0 + int(immune_calls / 6 * (x1 - x0))
    py = y1 - int(tf_calls / 250 * (y1 - y0))
    return px, py


def draw_priority_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, data: dict[str, object]) -> None:
    color = COLORS["green"]
    draw_panel_header(draw, x, y, w, "3", "Follow-up priority map", color)
    text(draw, (x + 54, y + 145), "4", F["stat"], color)
    text(draw, (x + 134, y + 155), "primary tissues", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "The priority surface combines immune-call breadth with TF-call breadth, then routes each tissue to the next evidence layer.",
        F["body"],
        w - 108,
        COLORS["muted"],
        9,
    )

    chart = (x + 120, y + 456, x + w - 86, y + 775)
    rounded(draw, (x + 54, y + 420, x + w - 54, y + 878), 24, "#0D1520", "#2A394D", 1)
    draw_scatter_axes(draw, chart)
    for row in data["priority"]:
        tissue = row["tissue"]
        px, py = scatter_point(chart, int(row["immune_calls"]), int(row["tf_calls"]))
        fill = TISSUE_COLORS[tissue]
        r = 19 if tissue in {"skin", "thymus", "kidney", "liver"} else 12
        draw.ellipse((px - r, py - r, px + r, py + r), fill=fill, outline="#0C111A", width=3)
        if tissue in {"skin", "thymus", "kidney", "liver"}:
            label_x = px + 18
            label_y = py - 28
            if tissue == "skin":
                label_x = px - 94
            if tissue == "liver":
                label_x = px + 20
                label_y = py - 8
            text(draw, (label_x, label_y), tissue_label(tissue), F["tiny_bold"], COLORS["text"])
    text(draw, (x + 88, y + 830), "Tissues near the lower-left remain context rows for this layer.", F["tiny"], COLORS["muted"])

    rounded(draw, (x + 54, y + 872, x + w - 54, y + 916), 18, "#172335", COLORS["green"], 2)
    active_pairs = data["active_pairs"]
    if active_pairs:
        pair = active_pairs[0]
        label = f"L-R context: {tissue_label(pair['source'])} -> {tissue_label(pair['target'])} {pair['ligand']}/{pair['receptor']}"
    else:
        label = "L-R context carried as a separate signaling layer"
    text(draw, (x + 82, y + 885), label, F["tiny_bold"], COLORS["text"])

    queue_y = y + 962
    queue = [
        ("Skin", "immune + TF", COLORS["rose"]),
        ("Thymus", "T-cell + TF", COLORS["violet"]),
        ("Kidney", "B/mono + TF", COLORS["green"]),
        ("Liver", "TF follow-up", COLORS["amber"]),
    ]
    card_w = (w - 108 - 3 * 18) // 4
    for i, (head, body, c) in enumerate(queue):
        cx = x + 54 + i * (card_w + 18)
        rounded(draw, (cx, queue_y, cx + card_w, queue_y + 74), 20, "#172335", c, 2)
        text(draw, (cx + 18, queue_y + 13), head, F["small_bold"], COLORS["text"])
        text(draw, (cx + 18, queue_y + 45), body, F["tiny"], COLORS["muted"])


def draw_reader_rule(draw: ImageDraw.ImageDraw) -> None:
    box = (150, 1695, 3690, 1828)
    rounded(draw, box, 30, "#111A27", "#2A394D", 2)
    text(draw, (196, 1736), "Layer decoder", F["h3"], COLORS["teal"])
    steps = [
        ("Composition", "mMCP-counter cell types"),
        ("Regulation", "decoupler TF activity"),
        ("Priority surface", "immune x TF breadth"),
        ("Next layer", "target + biomarker evidence"),
    ]
    x = 700
    for i, (head, body) in enumerate(steps):
        color = [COLORS["rose"], COLORS["violet"], COLORS["green"], COLORS["teal"]][i]
        rounded(draw, (x, 1718, x + 620, 1808), 24, "#172335", color, 2)
        text(draw, (x + 28, 1733), head, F["small_bold"], COLORS["text"])
        text(draw, (x + 28, 1768), body, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            ax = x + 640
            draw.line((ax, 1763, ax + 82, 1763), fill="#6F7E90", width=4)
            draw.polygon([(ax + 82, 1763), (ax + 62, 1752), (ax + 62, 1774)], fill="#6F7E90")
        x += 730


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    source = "Takeaway: immune and TF activity prioritize tissues and axes for follow-up biology."
    scope = "Next: target and biomarker layers translate these axes into candidate follow-up packets."
    paragraph(draw, (150, 1888), source, F["small"], 3440, COLORS["muted"], 7)
    paragraph(draw, (150, 1993), scope, F["small_bold"], 3440, COLORS["text"], 7)


def main() -> None:
    data = load_data()
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)
    draw_background(draw)
    draw_header(draw)
    draw_title(draw)

    panel_y = 625
    panel_w = 1090
    gap = 55
    x1 = 150
    x2 = x1 + panel_w + gap
    x3 = x2 + panel_w + gap
    draw_immune_panel(draw, x1, panel_y, panel_w, data)
    draw_tf_panel(draw, x2, panel_y, panel_w, data)
    draw_priority_panel(draw, x3, panel_y, 1295, data)
    draw_reader_rule(draw)
    draw_footer(draw)

    canvas.save(OUT_PATH, quality=95)
    gray = ImageOps.grayscale(canvas).convert("RGB")
    gray.save(GRAY_PATH, quality=95)

    stat = ImageStat.Stat(gray)
    manifest = {
        "slide": 38,
        "title": "Immune And TF Activity Prioritize Follow-Up",
        "readout": "Skin, thymus, kidney, and liver define the immune/TF follow-up priority surface.",
        "asset": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "sources": [str(path.relative_to(ROOT)) for path in [*IMMUNE_FILES, *TF_FILES, SIGNALING, SOURCE_FIG]],
        "data": {
            "immune_top": [
                {
                    "tissue": row["tissue"],
                    "n_samples": row["n_samples"],
                    "n_cell_types": row["n_cell_types"],
                    "n_significant": row["n_significant"],
                }
                for row in data["immune"][:4]
            ],
            "skin_significant_cell_types": data["immune_by_tissue"]["skin"]["significant"],
            "tf_top": [
                {
                    "tissue": row["tissue"],
                    "n_samples": row["n_samples"],
                    "n_tfs_tested": row["n_tfs_tested"],
                    "n_significant": row["n_significant"],
                    "top_tfs": row["significant"][:6],
                }
                for row in data["tf"][:4]
            ],
            "priority_top": data["priority"][:4],
            "active_lr_pairs": data["active_pairs"],
        },
        "visual_qa": {
            "dimensions": [W, H],
            "mean_grayscale": round(stat.mean[0], 2),
            "stddev_grayscale": round(stat.stddev[0], 2),
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"asset": str(OUT_PATH.relative_to(ROOT)), "grayscale": str(GRAY_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

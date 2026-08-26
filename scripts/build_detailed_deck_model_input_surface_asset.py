#!/usr/bin/env python3
"""Build slide 6 asset: what the model actually sees."""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
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
OUT_DIR = ASSET_ROOT / "model_input_surface"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "what_the_model_actually_sees_premium.png"
GRAY_PATH = OUT_DIR / "what_the_model_actually_sees_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "model_input_surface_manifest.json"

TASK_MANIFEST = ROOT / "v9" / "task_manifests" / "A4_thymus_bulk_lomo.json"
METADATA = ROOT / "processed" / "A_detection" / "thymus" / "thymus_all_missions_metadata.csv"
MATRIX = ROOT / "processed" / "A_detection" / "thymus" / "thymus_all_missions_log2_norm.csv"
PATHWAY_DIR = ROOT / "processed" / "pathway_scores" / "thymus"
METHOD_HARDENING = ASSET_ROOT / "method_hardening" / "method_hardening_manifest.json"

COLORS = {
    "bg": "#0C111A",
    "bg2": "#091019",
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
    "ink": "#080D14",
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
    "h3": load_font(33, True),
    "body": load_font(29),
    "body_bold": load_font(29, True),
    "small": load_font(26),
    "small_bold": load_font(26, True),
    "tiny": load_font(22),
    "tiny_bold": load_font(22, True),
    "micro": load_font(19),
    "micro_bold": load_font(19, True),
    "mono": load_font(24),
    "mono_bold": load_font(24, True),
    "stat": load_font(62, True),
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
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def label_name(raw: str) -> str:
    if raw == "GC":
        return "Ground"
    return raw


def load_pathway_counts() -> dict[str, dict[str, int]]:
    counts: dict[str, dict[str, int]] = defaultdict(dict)
    for path in sorted(PATHWAY_DIR.glob("*_gsva_*.csv")):
        mission, db = path.stem.split("_gsva_")
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.reader(handle)
            header = next(reader)
            n_rows = sum(1 for _ in reader)
        counts[db][mission] = len(header) - 1
        counts[f"{db}_rows"][mission] = n_rows
    return dict(counts)


def matrix_preview(scored_ids: set[str], rows: int = 8, cols: int = 8) -> dict[str, object]:
    samples: list[str] = []
    values: list[list[float]] = []
    genes: list[str] = []
    with MATRIX.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        genes = header[1 : cols + 1]
        for row in reader:
            if row[0] not in scored_ids:
                continue
            samples.append(row[0])
            values.append([float(v) for v in row[1 : cols + 1]])
            if len(samples) >= rows:
                break
    flat = [value for row in values for value in row]
    return {
        "samples": samples,
        "genes": genes,
        "values": values,
        "min": min(flat) if flat else 0.0,
        "max": max(flat) if flat else 1.0,
    }


def pathway_preview(path: Path, rows: int = 7, cols: int = 8) -> dict[str, object]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        features = [name.replace("HALLMARK_", "") for name in header[1 : cols + 1]]
        sample_ids: list[str] = []
        values: list[list[float]] = []
        for row in reader:
            sample_ids.append(row[0])
            values.append([float(v) for v in row[1 : cols + 1]])
            if len(sample_ids) >= rows:
                break
    flat = [value for row in values for value in row]
    return {
        "samples": sample_ids,
        "features": features,
        "values": values,
        "min": min(flat) if flat else -1.0,
        "max": max(flat) if flat else 1.0,
    }


def load_example_score() -> dict[str, object]:
    hardening = json.loads(METHOD_HARDENING.read_text())
    for tissue in hardening["metrics"]["heatmap"]:
        if tissue["tissue"] != "thymus":
            continue
        for feature in tissue["features"]:
            if feature["feature"] == "pathway_kegg":
                return feature
    return {"feature": "pathway_kegg", "model": "pca_lr", "auroc": 0.0, "perm_p": 0.0}


def load_data() -> dict[str, object]:
    manifest = json.loads(TASK_MANIFEST.read_text())
    metadata = read_csv(METADATA)
    scoring_rows = [row for row in metadata if row["label"] in {"Flight", "GC"}]
    scored_ids = {row[""] for row in scoring_rows}
    label_counts = Counter(label_name(row["label"]) for row in scoring_rows)
    mission_counts = Counter(row["mission"] for row in scoring_rows)
    fold = manifest["split"]["folds"][0]
    path_counts = load_pathway_counts()
    hallway = path_counts["hallmark"]
    kegg = path_counts["kegg"]
    return {
        "task_id": manifest["task_id"],
        "tissue": manifest["tissue"],
        "scored_rows": len(scoring_rows),
        "label_counts": dict(label_counts),
        "mission_counts": dict(mission_counts),
        "n_folds": manifest["split"]["n_folds"],
        "genes": int(fold["n_genes_after_var_filter"]),
        "genes_before": int(fold["n_genes_before_var_filter"]),
        "pathway_counts": path_counts,
        "hallmark_features": min(hallway.values()),
        "kegg_min": min(kegg.values()),
        "kegg_max": max(kegg.values()),
        "matrix_preview": matrix_preview(scored_ids),
        "pathway_preview": pathway_preview(PATHWAY_DIR / "RR-6_gsva_hallmark.csv"),
        "example_score": load_example_score(),
    }


def draw_background(draw: ImageDraw.ImageDraw) -> None:
    draw.rectangle((0, 0, W, H), fill=COLORS["bg"])
    for y in range(H):
        draw.line((0, y, W, y), fill=blend(COLORS["bg"], COLORS["bg2"], y / H))
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=COLORS["grid"], width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill="#172234", width=1)
    draw.rectangle((0, 0, W, 315), fill="#0F1824")
    draw.rectangle((0, 1840, W, H), fill=COLORS["ink"])
    draw.line((0, 315, W, 315), fill="#1E2B3D", width=2)
    draw.line((0, 1840, W, 1840), fill="#1E2B3D", width=2)


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(205, int(draw.textlength(value, font=F["tiny_bold"]) + 76))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def draw_header(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    text(draw, (150, 76), "METHOD / INPUT SURFACE", F["kicker"], COLORS["teal"])
    x = 2060
    x = badge(draw, x, 66, "samples", f"{data['scored_rows']} rows", COLORS["blue"])
    x = badge(draw, x, 66, "gene view", f"{data['genes']:,} cols", COLORS["teal"])
    x = badge(draw, x, 66, "pathway view", f"{data['hallmark_features']}+ cols", COLORS["amber"])
    badge(draw, x, 66, "split", f"{data['n_folds']} folds", COLORS["violet"])


def draw_title(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 382), "What The Model Actually Sees", F["title"], COLORS["text"])
    paragraph(
        draw,
        (155, 493),
        "A model receives numeric views of the same samples: gene columns, pathway scores, or compact representations.",
        F["subtitle"],
        1790,
        COLORS["muted"],
        10,
    )


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: str, width: int = 5) -> None:
    x0, y0 = start
    x1, y1 = end
    draw.line((x0, y0, x1, y1), fill=color, width=width)
    draw.polygon([(x1, y1), (x1 - 24, y1 - 14), (x1 - 24, y1 + 14)], fill=color)


def value_color(value: float, lo: float, hi: float) -> tuple[int, int, int, int]:
    if hi <= lo:
        t = 0.5
    else:
        t = (value - lo) / (hi - lo)
    if t < 0.5:
        return rgba(blend(COLORS["blue"], COLORS["teal"], t * 2), 220)
    return rgba(blend(COLORS["teal"], COLORS["amber"], (t - 0.5) * 2), 225)


def draw_heatmap(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    values: list[list[float]],
    lo: float,
    hi: float,
    cell_w: int,
    cell_h: int,
) -> None:
    for r, row in enumerate(values):
        for c, value in enumerate(row):
            box = (
                x + c * cell_w,
                y + r * cell_h,
                x + c * cell_w + cell_w - 4,
                y + r * cell_h + cell_h - 4,
            )
            draw.rounded_rectangle(box, radius=5, fill=value_color(value, lo, hi))


def label_chip(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, color: str, width: int = 132) -> None:
    rounded(draw, (x, y, x + width, y + 42), 13, rgba(color, 42), color, 2)
    text(draw, (x + width / 2, y + 12), label, F["tiny_bold"], COLORS["text"], "ma")


def draw_sample_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x0, y0, x1, y1 = 150, 640, 930, 1585
    rounded(draw, (x0, y0, x1, y1), 24, COLORS["panel"], "#2A394D", 2)
    text(draw, (x0 + 42, y0 + 42), "Same sample rows", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 44, y0 + 100),
        "Each row is one thymus sample with a fixed label and mission.",
        F["body"],
        615,
        COLORS["muted"],
        8,
    )

    stat_y = y0 + 212
    stats = [
        (f"{data['scored_rows']}", "scoring rows", COLORS["blue"]),
        (str(len(data["mission_counts"])), "missions", COLORS["teal"]),
        ("2", "labels", COLORS["amber"]),
    ]
    stat_w = 210
    for i, (value, label, color) in enumerate(stats):
        sx = x0 + 44 + i * (stat_w + 22)
        rounded(draw, (sx, stat_y, sx + stat_w, stat_y + 126), 20, "#121D2B", color, 2)
        text(draw, (sx + 24, stat_y + 18), value, F["stat"], color)
        text(draw, (sx + 24, stat_y + 82), label, F["tiny_bold"], COLORS["muted"])

    table_x, table_y = x0 + 44, y0 + 395
    rounded(draw, (table_x, table_y, x1 - 44, table_y + 410), 18, "#0D1520", "#26364B", 2)
    headers = ["sample row", "label", "mission"]
    col_x = [table_x + 24, table_x + 355, table_x + 510]
    for hx, header in zip(col_x, headers):
        text(draw, (hx, table_y + 24), header.upper(), F["micro_bold"], COLORS["dim"])
    draw.line((table_x + 18, table_y + 62, x1 - 62, table_y + 62), fill="#26364B", width=2)

    label_counts = data["label_counts"]
    rows = [
        ("row 01", "Flight", "RR-6", COLORS["blue"]),
        ("row 02", "Ground", "RR-6", COLORS["teal"]),
        ("row 03", "Flight", "RR-9", COLORS["blue"]),
        ("row 04", "Ground", "RR-9", COLORS["teal"]),
        ("...", f"{label_counts.get('Flight', 0)} Flight", "4-fold", COLORS["amber"]),
        ("...", f"{label_counts.get('Ground', 0)} Ground", "LOMO", COLORS["violet"]),
    ]
    yy = table_y + 86
    for sample, label, mission, color in rows:
        draw.line((table_x + 18, yy + 45, x1 - 62, yy + 45), fill="#182436", width=1)
        text(draw, (col_x[0], yy + 8), sample, F["small_bold"], COLORS["text"])
        label_chip(draw, col_x[1], yy, label, color, width=132)
        text(draw, (col_x[2], yy + 8), mission, F["small_bold"], COLORS["muted"])
        yy += 54

    rounded(draw, (x0 + 44, y1 - 150, x1 - 44, y1 - 44), 18, "#111B28", COLORS["teal"], 2)
    text(draw, (x0 + 74, y1 - 126), "ML translation", F["tiny_bold"], COLORS["teal"])
    paragraph(
        draw,
        (x0 + 74, y1 - 94),
        "Feature = column. View = exposed columns.",
        F["small"],
        590,
        COLORS["text"],
        6,
    )


def draw_gene_visual(draw: ImageDraw.ImageDraw, x: int, y: int, data: dict[str, object]) -> None:
    preview = data["matrix_preview"]
    values = preview["values"][:6]
    draw_heatmap(draw, x, y, values, float(preview["min"]), float(preview["max"]), 31, 24)
    text(draw, (x, y + 158), "Actb  Gapdh  H2-Ab1  Cd3d  ...", F["micro_bold"], COLORS["muted"])


def draw_pathway_visual(draw: ImageDraw.ImageDraw, x: int, y: int, data: dict[str, object]) -> None:
    preview = data["pathway_preview"]
    draw_heatmap(draw, x, y, preview["values"][:6], float(preview["min"]), float(preview["max"]), 33, 24)
    text(draw, (x, y + 158), "E2F  G2M  IFN  OXPHOS  ...", F["micro_bold"], COLORS["muted"])
    node_y = y + 184
    for i, color in enumerate([COLORS["teal"], COLORS["amber"], COLORS["green"], COLORS["violet"]]):
        cx = x + 18 + i * 54
        draw.ellipse((cx, node_y, cx + 28, node_y + 28), fill=rgba(color, 190), outline=color, width=2)
        if i:
            draw.line((cx - 26, node_y + 14, cx, node_y + 14), fill="#5C6D82", width=3)


def draw_compressed_visual(draw: ImageDraw.ImageDraw, x: int, y: int) -> None:
    rounded(draw, (x, y - 10, x + 312, y + 138), 18, "#0D1520", "#26364B", 2)
    for i in range(12):
        bx = x + 25 + i * 23
        h = 30 + ((i * 31) % 74)
        color = [COLORS["violet"], COLORS["blue"], COLORS["teal"], COLORS["amber"]][i % 4]
        rounded(draw, (bx, y + 106 - h, bx + 14, y + 106), 5, rgba(color, 210))
    for i in range(3):
        draw.ellipse((x + 214 + i * 28, y + 16 + i * 19, x + 238 + i * 28, y + 40 + i * 19), fill=rgba(COLORS["violet"], 170))
    text(draw, (x + 24, y + 112), "PCA / encoder", F["micro_bold"], COLORS["muted"])


def draw_lane(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    number: str,
    title: str,
    stat: str,
    body: str,
    color: str,
    visual_fn,
    data: dict[str, object] | None = None,
) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 24, COLORS["panel2"], "#2A394D", 2)
    rounded(draw, (x0 + 24, y0 + 26, x0 + 92, y0 + 94), 18, rgba(color, 48), color, 2)
    text(draw, (x0 + 58, y0 + 44), number, F["h3"], color, "ma")
    text(draw, (x0 + 118, y0 + 31), title, F["h3"], COLORS["text"])
    rounded(draw, (x1 - 386, y0 + 28, x1 - 34, y0 + 86), 15, "#101A27", color, 2)
    text(draw, (x1 - 210, y0 + 45), stat, F["small_bold"], COLORS["text"], "ma")
    paragraph(draw, (x0 + 118, y0 + 83), body, F["small"], 690, COLORS["muted"], 6)
    if data is None:
        visual_fn(draw, x0 + 1035, y0 + 86)
    else:
        visual_fn(draw, x0 + 1035, y0 + 86, data)


def draw_feature_lanes(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x0, y0, x1, y1 = 1030, 640, 2675, 1585
    rounded(draw, (x0, y0, x1, y1), 24, "#0E1722", "#243247", 2)
    text(draw, (x0 + 38, y0 + 35), "Choose a feature view", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 40, y0 + 91),
        "The row identity stays fixed; the columns change what signal the learner can use.",
        F["body"],
        1080,
        COLORS["muted"],
        7,
    )

    lanes = [
        (
            (x0 + 38, y0 + 175, x1 - 38, y0 + 410),
            "1",
            "Gene expression view",
            f"67 x {data['genes']:,}",
            "Raw normalized gene columns. Classical models can read every selected mouse gene as a numeric feature.",
            COLORS["blue"],
            draw_gene_visual,
            data,
        ),
        (
            (x0 + 38, y0 + 445, x1 - 38, y0 + 680),
            "2",
            "Pathway score view",
            f"50 Hallmark | {data['kegg_min']}-{data['kegg_max']} KEGG",
            "GSVA summaries turn many genes into biological-program columns such as cell cycle or immune signaling.",
            COLORS["amber"],
            draw_pathway_visual,
            data,
        ),
        (
            (x0 + 38, y0 + 715, x1 - 38, y0 + 910),
            "3",
            "Compressed input view",
            "compact vector",
            "PCA, encoder, or adapter inputs keep the sample row while changing the coordinate system the model reads.",
            COLORS["violet"],
            draw_compressed_visual,
            None,
        ),
    ]

    for lane in lanes:
        draw_lane(draw, *lane)


def draw_probability_bar(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], label: str, value: float, color: str) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 14, "#0B121B", "#243247", 1)
    fill_w = int((x1 - x0 - 8) * value)
    rounded(draw, (x0 + 4, y0 + 4, x0 + 4 + fill_w, y1 - 4), 11, rgba(color, 190))
    text(draw, (x0 + 18, y0 + 13), label, F["tiny_bold"], COLORS["text"])
    text(draw, (x1 - 20, y0 + 13), f"{value:.2f}", F["tiny_bold"], COLORS["text"], "ra")


def draw_model_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x0, y0, x1, y1 = 2780, 640, 3690, 1585
    rounded(draw, (x0, y0, x1, y1), 24, COLORS["panel"], "#2A394D", 2)
    text(draw, (x0 + 42, y0 + 42), "One scoring gate", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 44, y0 + 100),
        "The evaluator stores which view entered the classifier before reading the metric.",
        F["body"],
        735,
        COLORS["muted"],
        8,
    )

    model_box = (x0 + 92, y0 + 214, x1 - 92, y0 + 392)
    rounded(draw, model_box, 26, "#111C2A", COLORS["teal"], 3)
    text(draw, ((model_box[0] + model_box[2]) / 2, model_box[1] + 34), "classifier", F["h1"], COLORS["text"], "ma")
    paragraph(
        draw,
        (model_box[0] + 58, model_box[1] + 90),
        "learns a Flight vs Ground separator from training missions",
        F["small"],
        model_box[2] - model_box[0] - 116,
        COLORS["muted"],
        6,
    )

    draw_probability_bar(draw, (x0 + 92, y0 + 456, x1 - 92, y0 + 512), "Flight probability", 0.73, COLORS["blue"])
    draw_probability_bar(draw, (x0 + 92, y0 + 528, x1 - 92, y0 + 584), "Ground probability", 0.27, COLORS["teal"])

    score = data["example_score"]
    score_box = (x0 + 60, y0 + 590, x1 - 60, y1 - 48)
    rounded(draw, score_box, 22, "#0D1520", "#26364B", 2)
    text(draw, (score_box[0] + 34, score_box[1] + 30), "Example score row", F["h3"], COLORS["text"])
    text(draw, (score_box[0] + 34, score_box[1] + 70), "stores task + split + feature view + model + metric", F["tiny_bold"], COLORS["amber"])
    rows = [
        ("task_id", "A4_thymus_bulk_lomo"),
        ("feature_view", str(score["feature"])),
        ("model", str(score["model"])),
        ("metric", f"AUROC {float(score['auroc']):.3f}"),
    ]
    yy = score_box[1] + 112
    for label, value in rows:
        rounded(draw, (score_box[0] + 34, yy - 8, score_box[2] - 34, yy + 36), 13, "#121D2B", "#223149", 1)
        text(draw, (score_box[0] + 56, yy + 3), label, F["tiny_bold"], COLORS["dim"])
        text(draw, (score_box[0] + 252, yy + 2), value, F["small_bold"], COLORS["text"])
        yy += 49


def draw_connectors(draw: ImageDraw.ImageDraw) -> None:
    lane_centers = [932, 1202, 1462]
    targets = [858, 990, 1118]
    colors = [COLORS["blue"], COLORS["amber"], COLORS["violet"]]
    for y, target_y, color in zip(lane_centers, targets, colors):
        arrow(draw, (930, y), (1068, y), color, 5)
        arrow(draw, (2630, y), (2780, target_y), color, 4)


def draw_reader_rule(draw: ImageDraw.ImageDraw) -> None:
    y0, y1 = 1625, 1816
    rounded(draw, (150, y0, 3690, y1), 24, "#0E1722", "#243247", 2)
    text(draw, (190, y0 + 35), "Reading path", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (190, y0 + 91),
        "Before comparing models, read the sample split and the feature view. The metric travels with those fields attached.",
        F["body"],
        880,
        COLORS["muted"],
        6,
    )
    steps = [
        ("1", "lock sample rows", "same task ID and mission-held-out split", COLORS["blue"]),
        ("2", "name feature view", "gene, pathway, combined, or compact input", COLORS["amber"]),
        ("3", "read metric", "AUROC / macro-F1 / balanced accuracy", COLORS["green"]),
    ]
    x = 1250
    for number, title, body, color in steps:
        rounded(draw, (x, y0 + 36, x + 730, y1 - 36), 22, "#111C2A", color, 2)
        rounded(draw, (x + 30, y0 + 66, x + 88, y0 + 124), 16, rgba(color, 50), color, 2)
        text(draw, (x + 59, y0 + 80), number, F["h3"], color, "ma")
        text(draw, (x + 112, y0 + 58), title, F["small_bold"], COLORS["text"])
        paragraph(draw, (x + 112, y0 + 96), body, F["small"], 540, COLORS["muted"], 6)
        x += 770


def draw_footer(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    text(draw, (150, 1902), "TAKEAWAY", F["micro_bold"], COLORS["dim"])
    paragraph(
        draw,
        (150, 1934),
        "The model sees a declared input surface: gene matrix, pathway scores, combined table, compact embedding, or prompt summary.",
        F["small"],
        1560,
        COLORS["muted"],
        6,
    )
    text(draw, (2170, 1902), "NEXT", F["micro_bold"], COLORS["dim"])
    paragraph(
        draw,
        (2170, 1934),
        "Performance slides use the same task, split, feature-view, model, and metric fields as the score row.",
        F["small"],
        1480,
        COLORS["muted"],
        6,
    )


def main() -> None:
    data = load_data()
    img = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(img, "RGBA")
    draw_background(draw)
    draw_header(draw, data)
    draw_title(draw)
    draw_sample_panel(draw, data)
    draw_feature_lanes(draw, data)
    draw_model_panel(draw, data)
    draw_connectors(draw)
    draw_reader_rule(draw)
    draw_footer(draw, data)

    img.save(OUT_PATH, quality=95)
    gray = ImageOps.grayscale(img).convert("RGB")
    gray.save(GRAY_PATH, quality=95)
    stat = ImageStat.Stat(gray)

    manifest = {
        "asset": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "source_files": [
            str(TASK_MANIFEST.relative_to(ROOT)),
            str(METADATA.relative_to(ROOT)),
            str(MATRIX.relative_to(ROOT)),
            str(PATHWAY_DIR.relative_to(ROOT)),
            str(METHOD_HARDENING.relative_to(ROOT)),
        ],
        "metrics": {
            "task_id": data["task_id"],
            "scored_rows": data["scored_rows"],
            "genes_after_var_filter": data["genes"],
            "hallmark_features": data["hallmark_features"],
            "kegg_features_min": data["kegg_min"],
            "kegg_features_max": data["kegg_max"],
            "example_score": data["example_score"],
            "grayscale_mean": stat.mean,
            "grayscale_stddev": stat.stddev,
        },
        "scope": "Method primer defining sample rows, feature views, classifier input, and score-row fields for the detailed deck.",
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"asset": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

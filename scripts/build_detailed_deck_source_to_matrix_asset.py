#!/usr/bin/env python3
"""Build slide 5 asset: from study record to sample-by-feature matrix."""

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
OUT_DIR = ASSET_ROOT / "source_to_matrix"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "from_source_record_to_matrix_premium.png"
GRAY_PATH = OUT_DIR / "from_source_record_to_matrix_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "source_to_matrix_manifest.json"

TASK_MANIFEST = ROOT / "v9" / "task_manifests" / "A4_thymus_bulk_lomo.json"
METADATA = ROOT / "processed" / "A_detection" / "thymus" / "thymus_all_missions_metadata.csv"
MATRIX = ROOT / "processed" / "A_detection" / "thymus" / "thymus_all_missions_log2_norm.csv"
SOURCE_INDEX = ROOT / "v9" / "source_inventory.csv"
REFERENCE_SCENE = ROOT / "output" / "premium_methods_scenes" / "methods_data_to_evaluation_overview.png"

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
    "stat": load_font(64, True),
    "mono": load_font(24),
    "mono_bold": load_font(24, True),
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


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sample_matrix_preview(path: Path, rows: int = 7, cols: int = 9) -> dict[str, object]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        genes = header[1 : cols + 1]
        samples: list[str] = []
        values: list[list[float]] = []
        for row in reader:
            samples.append(row[0])
            values.append([float(v) for v in row[1 : cols + 1]])
            if len(samples) >= rows:
                break
    flat = [v for row in values for v in row]
    return {"genes": genes, "samples": samples, "values": values, "min": min(flat), "max": max(flat)}


def load_data() -> dict[str, object]:
    manifest = json.loads(TASK_MANIFEST.read_text())
    metadata = read_csv(METADATA)
    source_rows = read_csv(SOURCE_INDEX)
    source_map = {row["source_id"]: row for row in source_rows}

    mission_counts = Counter(row["mission"] for row in metadata)
    label_counts = Counter(row["label"] for row in metadata)
    mission_osd: dict[str, Counter[str]] = defaultdict(Counter)
    scored_rows = [row for row in metadata if row["label"] in {"Flight", "GC"}]
    scored_mission_counts = Counter(row["mission"] for row in scored_rows)
    for row in metadata:
        mission_osd[row["mission"]][row["osd_id"]] += 1

    fold = manifest["split"]["folds"][0]
    genes_before = int(fold["n_genes_before_var_filter"])
    genes_after = int(fold["n_genes_after_var_filter"])
    excluded_context_rows = int(fold["excluded_bc_ag_samples"])
    matrix_preview = sample_matrix_preview(MATRIX)

    source_cards = []
    for osd_id in ["OSD-244", "OSD-289", "OSD-421"]:
        source = source_map.get(osd_id, {})
        missions = sorted([m for m, counts in mission_osd.items() if counts.get(osd_id, 0)])
        source_cards.append(
            {
                "source_id": osd_id,
                "mission": " + ".join(missions),
                "rows": sum(mission_osd[m][osd_id] for m in missions),
                "scored_rows": sum(scored_mission_counts[m] for m in missions),
                "tissue": source.get("tissue", "thymus"),
                "assay": source.get("assay_modality", "bulk_rna_seq").replace("_", " "),
            }
        )

    return {
        "manifest": manifest,
        "metadata": metadata,
        "source_cards": source_cards,
        "raw_sample_rows": len(metadata),
        "scored_rows": len(scored_rows),
        "excluded_context_rows": excluded_context_rows,
        "mission_counts": dict(mission_counts),
        "label_counts": dict(label_counts),
        "genes_before": genes_before,
        "genes_after": genes_after,
        "matrix_preview": matrix_preview,
        "folds": manifest["split"]["folds"],
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
    text(draw, (150, 76), "METHOD / DATA SURFACE", F["kicker"], COLORS["teal"])
    x = 2045
    x = badge(draw, x, 66, "studies", f"{len(data['source_cards'])} records", COLORS["blue"])
    x = badge(draw, x, 66, "aligned", f"{data['raw_sample_rows']} rows", COLORS["teal"])
    x = badge(draw, x, 66, "scored", f"{data['scored_rows']} rows", COLORS["green"])
    badge(draw, x, 66, "features", f"{data['genes_after']:,} genes", COLORS["violet"])


def draw_title(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 382), "From Study Record To Matrix", F["title"], COLORS["text"])
    paragraph(
        draw,
        (155, 490),
        "Public OSDR records become benchmark input when sample metadata and expression values are aligned into one sample-by-feature table.",
        F["subtitle"],
        2020,
        COLORS["muted"],
        10,
    )


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: str) -> None:
    x1, y1 = start
    x2, y2 = end
    draw.line((x1, y1, x2 - 20, y2), fill=color, width=4)
    draw.polygon([(x2, y2), (x2 - 24, y2 - 13), (x2 - 24, y2 + 13)], fill=color)


def source_card(draw: ImageDraw.ImageDraw, x: int, y: int, card: dict[str, object], color: str) -> None:
    rounded(draw, (x, y, x + 725, y + 126), 20, blend(COLORS["panel2"], color, 0.08), color, 2)
    text(draw, (x + 22, y + 18), str(card["source_id"]), F["h3"], COLORS["text"])
    text(draw, (x + 22, y + 58), str(card["mission"]), F["small_bold"], color)
    text(draw, (x + 22, y + 91), f"{card['rows']} aligned rows / {card['scored_rows']} scored", F["tiny"], COLORS["muted"])
    rounded(draw, (x + 560, y + 28, x + 688, y + 96), 16, "#0D1520", "#2A394D", 1)
    text(draw, (x + 624, y + 49), "bulk", F["micro_bold"], COLORS["muted"], anchor="ma")
    text(draw, (x + 624, y + 73), "RNA-seq", F["micro_bold"], COLORS["muted"], anchor="ma")


def draw_source_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    panel = (150, 680, 1030, 1595)
    rounded(draw, panel, 30, COLORS["panel"], "#2A394D", 2)
    text(draw, (200, 730), "A. Study records stay attached", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (200, 785),
        "Each study row carries mission, tissue, assay, species, and public accession context.",
        F["small"],
        760,
        COLORS["muted"],
        7,
    )
    colors = [COLORS["teal"], COLORS["blue"], COLORS["amber"]]
    for i, card in enumerate(data["source_cards"]):
        source_card(draw, 215, 910 + i * 148, card, colors[i])

    rounded(draw, (215, 1382, 940, 1518), 20, "#111C29", COLORS["green"], 2)
    text(draw, (242, 1410), "Sample metadata", F["small_bold"], COLORS["green"])
    paragraph(
        draw,
        (242, 1446),
        "sample id, mission, label, strain, duration, tissue, assay",
        F["tiny"],
        640,
        COLORS["muted"],
        5,
    )


def cell_color(value: float, lo: float, hi: float) -> str:
    if hi <= lo:
        t = 0.5
    else:
        t = (value - lo) / (hi - lo)
    if t < 0.5:
        return blend("#12314B", COLORS["blue"], t * 1.4)
    return blend(COLORS["blue"], COLORS["amber"], (t - 0.5) * 1.15)


def draw_matrix_crop(draw: ImageDraw.ImageDraw, x: int, y: int, data: dict[str, object], cell: int = 48) -> None:
    values = data["values"]
    lo = float(data["min"])
    hi = float(data["max"])
    genes = data["genes"]
    samples = data["samples"]
    rounded(draw, (x, y, x + 735, y + 500), 22, "#0E1723", "#2A394D", 2)
    text(draw, (x + 28, y + 28), "Actual log2-normalized matrix crop", F["small_bold"], COLORS["teal"])
    text(draw, (x + 28, y + 64), "rows are samples; columns are mouse genes", F["tiny"], COLORS["muted"])
    grid_x = x + 155
    grid_y = y + 140
    for j, gene in enumerate(genes):
        text(draw, (grid_x + j * cell + cell / 2, y + 116), gene[-4:], F["micro_bold"], COLORS["dim"], anchor="mm")
    for i, sample in enumerate(samples):
        short = sample.split("_")[-1][:7]
        text(draw, (x + 28, grid_y + i * cell + cell / 2 - 8), short, F["micro"], COLORS["muted"])
        for j, value in enumerate(values[i]):
            fill = cell_color(float(value), lo, hi)
            rounded(draw, (grid_x + j * cell, grid_y + i * cell, grid_x + j * cell + 38, grid_y + i * cell + 38), 8, fill, None, 0)
    text(draw, (x + 28, y + 458), f"value range in crop: {lo:.1f} to {hi:.1f}", F["micro_bold"], COLORS["dim"])


def filter_bar(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, kept: int, total: int, color: str) -> None:
    rounded(draw, (x, y, x + 615, y + 104), 18, "#111C29", "#26384E", 1)
    text(draw, (x + 22, y + 16), label.upper(), F["micro_bold"], color)
    text(draw, (x + 22, y + 47), f"{kept:,} / {total:,}", F["h3"], COLORS["text"])
    bar_x = x + 260
    bar_y = y + 45
    draw.line((bar_x, bar_y, bar_x + 315, bar_y), fill="#2B394C", width=18)
    draw.line((bar_x, bar_y, bar_x + int(315 * kept / max(total, 1)), bar_y), fill=color, width=18)
    text(draw, (bar_x, y + 72), "kept for scoring surface", F["micro"], COLORS["muted"])


def draw_matrix_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    panel = (1080, 680, 2675, 1595)
    rounded(draw, panel, 30, COLORS["panel"], "#2A394D", 2)
    text(draw, (1130, 730), "B. Metadata and expression align", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (1130, 785),
        "A task joins sample annotations to normalized expression values, then keeps the rows and genes used for scoring.",
        F["small"],
        1390,
        COLORS["muted"],
        7,
    )

    draw_matrix_crop(draw, 1130, 890, data["matrix_preview"])
    filter_bar(draw, 1925, 930, "sample rows", int(data["scored_rows"]), int(data["raw_sample_rows"]), COLORS["green"])
    filter_bar(draw, 1925, 1065, "gene features", int(data["genes_after"]), int(data["genes_before"]), COLORS["violet"])
    rounded(draw, (1925, 1205, 2540, 1332), 18, blend(COLORS["panel2"], COLORS["amber"], 0.08), COLORS["amber"], 2)
    text(draw, (1950, 1230), "Context rows", F["small_bold"], COLORS["amber"])
    paragraph(
        draw,
        (1950, 1266),
        f"{data['excluded_context_rows']} BC/AG rows remain in metadata; the scoring surface uses Flight/Ground rows.",
        F["tiny"],
        555,
        COLORS["muted"],
        5,
    )

    rounded(draw, (1925, 1378, 2540, 1495), 18, blend(COLORS["panel2"], COLORS["teal"], 0.08), COLORS["teal"], 2)
    text(draw, (1950, 1404), "Matrix orientation", F["small_bold"], COLORS["teal"])
    text(draw, (1950, 1442), "samples x mouse-gene features", F["body_bold"], COLORS["text"])


def draw_table_icon(draw: ImageDraw.ImageDraw, x: int, y: int, rows: int, cols: int, color: str) -> None:
    cell = 28
    rounded(draw, (x, y, x + cols * cell + 22, y + rows * cell + 22), 18, "#0C1520", color, 2)
    for i in range(rows):
        for j in range(cols):
            t = ((i + 1) * (j + 2)) % 9 / 8
            fill = blend("#173149", color, t * 0.8)
            rounded(draw, (x + 12 + j * cell, y + 12 + i * cell, x + 12 + j * cell + 18, y + 12 + i * cell + 18), 5, fill)


def draw_output_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    panel = (2725, 680, 3690, 1595)
    rounded(draw, panel, 30, COLORS["panel"], "#2A394D", 2)
    text(draw, (2775, 730), "C. Benchmark input surface", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (2775, 785),
        "The resulting table is the shared input surface used by classical ML and by later feature-view transformations.",
        F["small"],
        825,
        COLORS["muted"],
        7,
    )

    draw_table_icon(draw, 2815, 910, 8, 10, COLORS["teal"])
    text(draw, (3185, 925), "Scoring matrix", F["h3"], COLORS["text"])
    stat_rows = [
        ("rows", f"{data['scored_rows']} samples", COLORS["green"]),
        ("columns", f"{data['genes_after']:,} genes", COLORS["violet"]),
        ("labels", "Flight / Ground", COLORS["amber"]),
        ("task", "A4 thymus bulk LOMO", COLORS["blue"]),
    ]
    yy = 980
    for label, value, color in stat_rows:
        text(draw, (3185, yy), label.upper(), F["micro_bold"], color)
        text(draw, (3320, yy), value, F["tiny_bold"], COLORS["text"])
        yy += 46

    rounded(draw, (2790, 1230, 3638, 1378), 22, blend(COLORS["panel2"], COLORS["blue"], 0.08), COLORS["blue"], 2)
    text(draw, (2820, 1260), "Next slide", F["small_bold"], COLORS["blue"])
    paragraph(
        draw,
        (2820, 1296),
        "The same sample table can be viewed as genes, pathways, or compressed model inputs.",
        F["tiny"],
        760,
        COLORS["muted"],
        5,
    )

    rounded(draw, (2790, 1422, 3638, 1518), 18, "#111C29", COLORS["green"], 2)
    text(draw, (2820, 1450), "Reader use", F["small_bold"], COLORS["green"])
    text(draw, (3010, 1450), "Matrix construction comes before model scoring.", F["tiny_bold"], COLORS["text"])


def stage_chip(draw: ImageDraw.ImageDraw, x: int, y: int, number: str, label: str, detail: str, color: str) -> int:
    w = 618
    rounded(draw, (x, y, x + w, y + 120), 20, blend(COLORS["panel2"], color, 0.08), color, 2)
    rounded(draw, (x + 22, y + 30, x + 82, y + 90), 16, blend(COLORS["panel2"], color, 0.22), color, 2)
    text(draw, (x + 52, y + 63), number, F["h3"], COLORS["text"], anchor="mm")
    text(draw, (x + 108, y + 26), label, F["small_bold"], COLORS["text"])
    text(draw, (x + 108, y + 63), detail, F["tiny"], COLORS["muted"])
    return x + w


def draw_reader_rule(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (150, 1635, 3690, 1815), 24, COLORS["panel"], "#2A394D", 2)
    stages = [
        ("1", "Study", "OSD record", COLORS["blue"]),
        ("2", "Samples", "metadata rows", COLORS["teal"]),
        ("3", "Expression", "normalized counts", COLORS["amber"]),
        ("4", "Task filter", "rows and genes kept", COLORS["violet"]),
        ("5", "Matrix", "sample x feature table", COLORS["green"]),
    ]
    x = 190
    for i, stage in enumerate(stages):
        end_x = stage_chip(draw, x, 1665, *stage)
        if i < len(stages) - 1:
            arrow(draw, (end_x + 8, 1725), (end_x + 60, 1725), COLORS["dim"])
        x = end_x + 82


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (150, 1886, 3690, 2072), 20, COLORS["panel"], "#2A394D", 2)
    text(draw, (190, 1938), "Takeaway", F["tiny_bold"], COLORS["blue"])
    paragraph(
        draw,
        (360, 1934),
        "Study rows, sample metadata, and expression values become one scoring matrix before any model is compared.",
        F["small"],
        2840,
        COLORS["muted"],
        6,
    )
    text(draw, (190, 2005), "Next", F["tiny_bold"], COLORS["amber"])
    paragraph(
        draw,
        (360, 2000),
        "The next slide separates gene, pathway, combined, compact, and prompt input surfaces.",
        F["small"],
        2840,
        COLORS["muted"],
        6,
    )
    text(draw, (3525, 2008), "5", F["h2"], COLORS["teal"])


def build() -> None:
    data = load_data()
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image)
    draw_background(draw)
    draw_header(draw, data)
    draw_title(draw)
    draw_source_panel(draw, data)
    arrow(draw, (1038, 1135), (1075, 1135), COLORS["teal"])
    draw_matrix_panel(draw, data)
    arrow(draw, (2684, 1135), (2720, 1135), COLORS["teal"])
    draw_output_panel(draw, data)
    draw_reader_rule(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=95)
    gray = ImageOps.grayscale(image).convert("RGB")
    gray.save(GRAY_PATH, quality=95)

    stat = ImageStat.Stat(ImageOps.grayscale(image))
    manifest = {
        "title": "From Study Record To Matrix",
        "outputs": {
            "png": str(OUT_PATH.relative_to(ROOT)),
            "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        },
        "sources": [
            str(TASK_MANIFEST.relative_to(ROOT)),
            str(METADATA.relative_to(ROOT)),
            str(MATRIX.relative_to(ROOT)),
            str(SOURCE_INDEX.relative_to(ROOT)),
            str(REFERENCE_SCENE.relative_to(ROOT)),
        ],
        "summary": {
            "source_records": len(data["source_cards"]),
            "aligned_sample_rows": data["raw_sample_rows"],
            "scored_sample_rows": data["scored_rows"],
            "context_rows": data["excluded_context_rows"],
            "genes_before_filter": data["genes_before"],
            "genes_after_filter": data["genes_after"],
        },
        "automatic_metrics": {
            "size": [W, H],
            "mode": image.mode,
            "grayscale_mean": round(stat.mean[0], 2),
            "grayscale_stddev": round(stat.stddev[0], 2),
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "png": str(OUT_PATH.relative_to(ROOT)),
                "grayscale": str(GRAY_PATH.relative_to(ROOT)),
                "manifest": str(MANIFEST_PATH.relative_to(ROOT)),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    build()

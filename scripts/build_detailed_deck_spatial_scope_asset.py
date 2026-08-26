#!/usr/bin/env python3
"""Build slide 36 asset: spatial weak-signal cases sharpen the readout."""

from __future__ import annotations

import json
import statistics
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
OUT_DIR = ASSET_ROOT / "spatial_scope"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "spatial_weak_signal_cases_define_scope_premium.png"
GRAY_PATH = OUT_DIR / "spatial_weak_signal_cases_define_scope_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "spatial_scope_manifest.json"

F3A = ROOT / "v3" / "evaluation" / "F3a_visium_classification.json"
F3B = ROOT / "v3" / "evaluation" / "F3b_visium_svg.json"
F3D = ROOT / "v3" / "evaluation" / "F3d_visium_cross_resolution.json"
SOURCE_FIG = ROOT / "v3" / "figures" / "v3_Fig2_spatial_overview.html"

COLORS = {
    "bg": "#0C111A",
    "bg2": "#091019",
    "header": "#0F1824",
    "footer": "#080D14",
    "panel": "#101823",
    "panel2": "#151F2D",
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

COND_COLORS = {"FLT": "#F17C88", "GC": "#66A6E8"}


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


def pill(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, color: str, font: ImageFont.ImageFont = F["tiny_bold"]) -> int:
    width = int(draw.textlength(label, font=font) + 42)
    rounded(draw, (x, y, x + width, y + 42), 18, "#172335", color, 2)
    text(draw, (x + 21, y + 11), label, font, COLORS["text"])
    return x + width + 14


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(230, int(draw.textlength(value, font=F["tiny_bold"]) + 76))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def load_data() -> dict[str, object]:
    f3a = json.loads(F3A.read_text())
    f3b = json.loads(F3B.read_text())
    f3d = json.loads(F3D.read_text())

    spot_counts = f3a["spot_counts"]
    mapping = f3a["section_mapping"]
    sections = []
    for section, meta in mapping.items():
        sections.append(
            {
                "section": section,
                "animal": meta["animal"],
                "condition": meta["condition"],
                "slide": section.split("_")[1],
                "n_spots": int(spot_counts[section]["n_spots"]),
                "n_genes": int(spot_counts[section]["n_genes"]),
                "cluster_id": int(meta["cluster_id"]),
            }
        )
    sections.sort(key=lambda row: (row["slide"], row["section"]))

    svg_counts = []
    for section, count in f3b["section_svg_counts"].items():
        svg_counts.append(
            {
                "section": section,
                "condition": mapping[section]["condition"],
                "slide": section.split("_")[1],
                "count": int(count),
                "spots": int(spot_counts[section]["n_spots"]),
            }
        )
    svg_counts.sort(key=lambda row: (row["slide"], row["section"]))

    return {
        "dataset": f3a["dataset"],
        "mission": f3a["mission"],
        "tissue": f3a["tissue"],
        "sections": sections,
        "n_sections": len(sections),
        "n_animals": len({row["animal"] for row in sections}),
        "n_spots_total": sum(row["n_spots"] for row in sections),
        "spot_range": [min(row["n_spots"] for row in sections), max(row["n_spots"] for row in sections)],
        "median_spots": statistics.median(row["n_spots"] for row in sections),
        "n_genes": int(f3a["n_genes"]),
        "section_level": f3a["section_level"],
        "animal_level": f3a["animal_level"],
        "bulk_auroc": f3d["bulk_auroc"],
        "spatial_cross": f3d["spatial_auroc"],
        "svg_counts": svg_counts,
        "svg_summary": f3b["svg_count_comparison"],
        "svg_overlap": f3b["svg_overlap"],
        "specific": f3b["condition_specific_svgs"],
        "pc1_batch_percent": 42.5,
        "bulk_deg_count": f3d["deg_overlap"]["bulk_deg_count"],
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
    text(draw, (150, 76), "BIOLOGY / MODALITY READOUT", F["kicker"], COLORS["teal"])
    x = 2180
    x = badge(draw, x, 66, "dataset", "OSD-352 RR-3", COLORS["blue"])
    x = badge(draw, x, 66, "assay", "Visium brain", COLORS["green"])
    x = badge(draw, x, 66, "surface", "12 sections", COLORS["amber"])
    badge(draw, x, 66, "readout", "weak-signal map", COLORS["violet"])


def draw_title(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 382), "Spatial Weak-Signal Cases Sharpen The Readout", F["title"], COLORS["text"])
    paragraph(
        draw,
        (155, 493),
        "OSD-352 brain Visium has real spatial structure; FLT/GC separation remains modest at the current sampling depth.",
        F["subtitle"],
        2300,
        COLORS["muted"],
        10,
    )


def draw_panel_header(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, step: str, title: str, color: str) -> None:
    rounded(draw, (x, y, x + w, y + 1010), 34, COLORS["panel"], "#29374A", 2)
    rounded(draw, (x + 34, y + 34, x + 98, y + 98), 20, "#172335", color, 2)
    text(draw, (x + 66, y + 51), step, F["h3"], COLORS["text"], "ma")
    text(draw, (x + 120, y + 44), title, F["h2"], COLORS["text"])


def draw_input_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, data: dict[str, object]) -> None:
    draw_panel_header(draw, x, y, w, "1", "Visium input surface", COLORS["teal"])

    text(draw, (x + 54, y + 142), f"{data['n_spots_total']:,}", F["stat"], COLORS["text"])
    text(draw, (x + 56, y + 220), "spots across 12 sections", F["h3"], COLORS["muted"])
    pill(draw, x + 56, y + 278, "3 FLT + 3 GC animals", COLORS["rose"])
    pill(draw, x + 340, y + 278, f"{data['n_genes']:,} genes", COLORS["blue"])

    paragraph(
        draw,
        (x + 54, y + 352),
        "The spatial task has enough observed tissue area to inspect section quality before interpreting the classifier readout.",
        F["body"],
        w - 108,
        COLORS["muted"],
        9,
    )

    section_y = y + 520
    section_w = 235
    section_h = 94
    gap_x = 27
    slides = ["158", "159", "304"]
    sections_by_slide = {slide: [row for row in data["sections"] if row["slide"] == slide] for slide in slides}
    min_spot, max_spot = data["spot_range"]
    for i, slide in enumerate(slides):
        sx = x + 60 + i * (section_w + gap_x)
        text(draw, (sx, section_y - 42), f"slide {slide}", F["small_bold"], COLORS["text"])
        for j, row in enumerate(sections_by_slide[slide]):
            yy = section_y + j * 88
            color = COND_COLORS[row["condition"]]
            intensity = (row["n_spots"] - min_spot) / (max_spot - min_spot)
            fill = blend("#111A27", color, 0.28 + 0.42 * intensity)
            rounded(draw, (sx, yy, sx + section_w, yy + section_h - 18), 18, fill, color, 2)
            short = row["section"].replace("Sample_", "")
            text(draw, (sx + 18, yy + 18), short, F["tiny_bold"], COLORS["text"])
            text(draw, (sx + 18, yy + 47), row["condition"], F["micro_bold"], color)
            text(draw, (sx + 118, yy + 47), f"{row['n_spots']:,} spots", F["micro"], COLORS["muted"])

    rounded(draw, (x + 54, y + 905, x + w - 54, y + 964), 20, "#172335", "#2A394D", 1)
    text(draw, (x + 82, y + 922), "Quality cue", F["tiny_bold"], COLORS["teal"])
    text(draw, (x + 222, y + 922), f"spot range {min_spot:,}-{max_spot:,}; median {data['median_spots']:,.0f}", F["tiny"], COLORS["muted"])


def draw_auroc_bar(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    label: str,
    value: float,
    p_value: float,
    color: str,
    axis_x0: int,
    axis_x1: int,
) -> None:
    text(draw, (x, y + 6), label, F["small_bold"], COLORS["text"])
    draw.line((axis_x0, y + 30, axis_x1, y + 30), fill="#263244", width=16)
    chance_x = axis_x0 + int(0.5 * (axis_x1 - axis_x0))
    draw.line((chance_x, y + 4, chance_x, y + 58), fill="#78879A", width=2)
    fill_x = axis_x0 + int(max(0.01, value) * (axis_x1 - axis_x0))
    draw.line((axis_x0, y + 30, fill_x, y + 30), fill=color, width=16)
    draw.ellipse((fill_x - 17, y + 13, fill_x + 17, y + 47), fill=color)
    text(draw, (axis_x1 + 26, y + 4), f"{value:.3f}", F["small_bold"], COLORS["text"])
    text(draw, (axis_x1 + 26, y + 34), f"p={p_value:.1f}", F["tiny"], COLORS["dim"])


def draw_readout_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, data: dict[str, object]) -> None:
    draw_panel_header(draw, x, y, w, "2", "Classifier readout", COLORS["amber"])

    section = data["section_level"]
    animal = data["animal_level"]
    bulk = data["bulk_auroc"]
    text(draw, (x + 54, y + 145), f"{section['auroc']:.3f}", F["stat2"], COLORS["amber"])
    text(draw, (x + 248, y + 155), "section AUROC", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "The classifier stays modest while exact permutation p-values remain broad at six animals.",
        F["body"],
        w - 108,
        COLORS["muted"],
        9,
    )

    axis_x0 = x + 360
    axis_x1 = x + w - 160
    text(draw, (axis_x0, y + 402), "AUROC guide: 0.5 overlap, 1.0 separation", F["tiny"], COLORS["dim"])
    rows = [
        ("Section-level", float(section["auroc"]), float(section["p_exact"]), COLORS["amber"]),
        ("Animal-level", float(animal["auroc"]), float(animal["p_exact"]), COLORS["blue"]),
        ("Companion bulk", float(bulk["auroc"]), float(bulk["p_exact"]), COLORS["rose"]),
    ]
    for i, row in enumerate(rows):
        draw_auroc_bar(draw, x + 64, y + 452 + i * 120, row[0], row[1], row[2], row[3], axis_x0, axis_x1)

    rounded(draw, (x + 54, y + 845, x + w - 54, y + 964), 24, "#172335", "#2A394D", 2)
    text(draw, (x + 82, y + 868), "Readout status", F["small_bold"], COLORS["amber"])
    paragraph(
        draw,
        (x + 82, y + 906),
        "Spatial structure is visible, but brain Visium FLT/GC separation remains modest at the current sample size.",
        F["tiny"],
        w - 160,
        COLORS["muted"],
        5,
    )


def draw_svg_dotplot(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, h: int, data: dict[str, object]) -> None:
    svg_counts = data["svg_counts"]
    max_count = max(row["count"] for row in svg_counts)
    axis_x0 = x + 78
    axis_x1 = x + w - 70
    baseline_y = y + h - 65
    draw.line((axis_x0, baseline_y, axis_x1, baseline_y), fill="#566275", width=2)
    text(draw, (axis_x0, y + 4), "SVG counts per section", F["tiny_bold"], COLORS["text"])
    for i, row in enumerate(svg_counts):
        px = axis_x0 + int((i + 0.5) / len(svg_counts) * (axis_x1 - axis_x0))
        py = baseline_y - int(row["count"] / max_count * (h - 105))
        color = COND_COLORS[row["condition"]]
        draw.line((px, baseline_y, px, py), fill=rgba(color, 120), width=5)
        draw.ellipse((px - 14, py - 14, px + 14, py + 14), fill=color)
        text(draw, (px, baseline_y + 12), row["section"].replace("Sample_", "").split("_")[1], F["micro"], COLORS["dim"], "ma")
    flt_mean = data["svg_summary"]["flt_mean"]
    gc_mean = data["svg_summary"]["gc_mean"]
    text(draw, (x + 82, y + h - 28), f"FLT mean {flt_mean:,.0f}", F["tiny_bold"], COND_COLORS["FLT"])
    text(draw, (x + 320, y + h - 28), f"GC mean {gc_mean:,.0f}", F["tiny_bold"], COND_COLORS["GC"])
    text(draw, (x + 550, y + h - 28), f"p={data['svg_summary']['mannwhitney_p']:.3f}", F["tiny"], COLORS["muted"])


def draw_overlap_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, data: dict[str, object]) -> None:
    overlap = data["svg_overlap"]
    specific = data["specific"]
    cx1, cx2, cy = x + 245, x + 370, y + 125
    r = 94
    draw.ellipse((cx1 - r, cy - r, cx1 + r, cy + r), fill=rgba(COND_COLORS["FLT"], 95), outline=COND_COLORS["FLT"], width=4)
    draw.ellipse((cx2 - r, cy - r, cx2 + r, cy + r), fill=rgba(COND_COLORS["GC"], 95), outline=COND_COLORS["GC"], width=4)
    text(draw, (x + 80, y + 6), "Spatial overlap", F["tiny_bold"], COLORS["text"])
    text(draw, (cx1 - 62, cy - 11), f"{specific['n_flt_specific']}", F["h3"], COLORS["text"])
    text(draw, (cx2 + 34, cy - 11), f"{specific['n_gc_specific']}", F["h3"], COLORS["text"])
    text(draw, ((cx1 + cx2) / 2, cy - 14), f"{overlap['n_intersection']:,}", F["h3"], COLORS["text"], "mm")
    text(draw, (x + 96, y + 223), "FLT-specific", F["micro_bold"], COND_COLORS["FLT"])
    text(draw, (x + 366, y + 223), "GC-specific", F["micro_bold"], COND_COLORS["GC"])


def draw_structure_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, data: dict[str, object]) -> None:
    draw_panel_header(draw, x, y, w, "3", "Spatial structure remains", COLORS["violet"])

    text(draw, (x + 54, y + 145), f"{data['svg_overlap']['jaccard']:.3f}", F["stat2"], COLORS["violet"])
    text(draw, (x + 258, y + 155), "SVG overlap", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "Spatial features are present; overlap and slide structure help interpret the modest classifier readout.",
        F["body"],
        w - 108,
        COLORS["muted"],
        9,
    )

    rounded(draw, (x + 54, y + 410, x + w - 54, y + 658), 24, "#0D1520", "#2A394D", 1)
    draw_svg_dotplot(draw, x + 54, y + 425, w - 108, 218, data)

    rounded(draw, (x + 54, y + 700, x + w - 54, y + 964), 24, "#0D1520", "#2A394D", 1)
    draw_overlap_panel(draw, x + 54, y + 720, w - 108, data)
    rounded(draw, (x + w - 394, y + 794, x + w - 93, y + 910), 22, "#172335", "#2A394D", 1)
    text(draw, (x + w - 365, y + 819), "PC1 cue", F["tiny_bold"], COLORS["amber"])
    text(draw, (x + w - 365, y + 854), f"{data['pc1_batch_percent']:.1f}% slide batch", F["small_bold"], COLORS["text"])


def draw_reader_rule(draw: ImageDraw.ImageDraw) -> None:
    box = (150, 1695, 3690, 1828)
    rounded(draw, box, 30, "#111A27", "#2A394D", 2)
    text(draw, (196, 1736), "Reading path", F["h3"], COLORS["teal"])
    steps = [
        ("Sample surface", "animals, sections, spots"),
        ("Classifier readout", "AUROC + exact p-values"),
        ("Spatial structure", "SVG overlap + batch cue"),
        ("Synthesis", "modest-signal map"),
    ]
    x = 690
    for i, (head, body) in enumerate(steps):
        color = [COLORS["teal"], COLORS["amber"], COLORS["violet"], COLORS["green"]][i]
        rounded(draw, (x, 1718, x + 630, 1808), 24, "#172335", color, 2)
        text(draw, (x + 28, 1733), head, F["small_bold"], COLORS["text"])
        text(draw, (x + 28, 1768), body, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            ax = x + 650
            draw.line((ax, 1763, ax + 82, 1763), fill="#6F7E90", width=4)
            draw.polygon([(ax + 82, 1763), (ax + 62, 1752), (ax + 62, 1774)], fill="#6F7E90")
        x += 735


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    takeaway = "Takeaway: spatial structure is visible, but FLT/GC separation stays modest at this sample size."
    next_line = "Temporal and single-cell context complete the handoff into systems-biology interpretation."
    paragraph(draw, (150, 1905), takeaway, F["small_bold"], 3440, COLORS["text"], 7)
    paragraph(draw, (150, 1993), next_line, F["small"], 3440, COLORS["muted"], 7)


def main() -> None:
    data = load_data()
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)
    draw_background(draw)
    draw_header(draw)
    draw_title(draw)

    panel_y = 625
    panel_w = 1160
    gap = 55
    xs = [150, 150 + panel_w + gap, 150 + (panel_w + gap) * 2]
    draw_input_panel(draw, xs[0], panel_y, panel_w, data)
    draw_readout_panel(draw, xs[1], panel_y, panel_w, data)
    draw_structure_panel(draw, xs[2], panel_y, panel_w, data)
    draw_reader_rule(draw)
    draw_footer(draw)

    canvas.save(OUT_PATH, quality=95)
    gray = ImageOps.grayscale(canvas).convert("RGB")
    gray.save(GRAY_PATH, quality=95)

    stat = ImageStat.Stat(gray)
    manifest = {
        "slide": 36,
        "title": "Spatial Weak-Signal Cases Sharpen The Readout",
        "claim": "OSD-352 brain Visium shows visible spatial structure with modest FLT/GC classifier separation at current sampling depth.",
        "asset": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "sources": [str(F3A.relative_to(ROOT)), str(F3B.relative_to(ROOT)), str(F3D.relative_to(ROOT)), str(SOURCE_FIG.relative_to(ROOT))],
        "data": {
            "n_animals": data["n_animals"],
            "n_sections": data["n_sections"],
            "n_spots_total": data["n_spots_total"],
            "n_genes": data["n_genes"],
            "section_auroc": data["section_level"]["auroc"],
            "animal_auroc": data["animal_level"]["auroc"],
            "bulk_auroc": data["bulk_auroc"]["auroc"],
            "flt_svg_mean": data["svg_summary"]["flt_mean"],
            "gc_svg_mean": data["svg_summary"]["gc_mean"],
            "svg_jaccard": data["svg_overlap"]["jaccard"],
            "pc1_batch_percent": data["pc1_batch_percent"],
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

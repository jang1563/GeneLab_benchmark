#!/usr/bin/env python3
"""Build slide 37 asset: systems biology adds interpretation."""

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
OUT_DIR = ASSET_ROOT / "systems_biology"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "systems_biology_adds_interpretation_premium.png"
GRAY_PATH = OUT_DIR / "systems_biology_adds_interpretation_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "systems_biology_manifest.json"

V5_EVAL = ROOT / "v5" / "evaluation"
IMMUNE_FILES = sorted(V5_EVAL.glob("immune_deconv_*.json"))
TF_FILES = sorted(V5_EVAL.glob("tf_activity_*.json"))
METABOLIC_FILES = sorted(V5_EVAL.glob("metabolic_flux_*.json"))
DRUG_TARGETS = V5_EVAL / "drug_targets.json"
BIOMARKERS = V5_EVAL / "consensus_biomarker_panel.json"
SOURCE_FIGS = [
    ROOT / "v5" / "figures" / "html" / "Fig7_immune_signaling.html",
    ROOT / "v5" / "figures" / "html" / "Fig8_metabolism_drugs.html",
    ROOT / "v5" / "figures" / "html" / "Fig9_consensus_panel.html",
]

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

LANE_COLORS = {
    "immune": COLORS["rose"],
    "tf": COLORS["violet"],
    "metabolic": COLORS["amber"],
    "target": COLORS["green"],
    "biomarker": COLORS["teal"],
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
    "h2": load_font(39, True),
    "h3": load_font(32, True),
    "body": load_font(28),
    "body_bold": load_font(28, True),
    "small": load_font(25),
    "small_bold": load_font(25, True),
    "tiny": load_font(21),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "stat": load_font(72, True),
    "stat2": load_font(62, True),
    "stat3": load_font(54, True),
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


def load_data() -> dict[str, object]:
    immune_rows = []
    for path in IMMUNE_FILES:
        item = json.loads(path.read_text())
        immune_rows.append(
            {
                "tissue": item["tissue"],
                "n_samples": int(item["n_samples"]),
                "n_cell_types": int(item["n_cell_types"]),
                "n_significant": int(item["n_significant_fdr05"]),
            }
        )
    immune_rows.sort(key=lambda row: (-row["n_significant"], row["tissue"]))

    tf_rows = []
    for path in TF_FILES:
        item = json.loads(path.read_text())
        tf_rows.append(
            {
                "tissue": item["tissue"],
                "n_samples": int(item["n_samples"]),
                "n_tfs_tested": int(item["n_tfs_tested"]),
                "n_significant": int(item["n_significant_fdr05"]),
            }
        )
    tf_rows.sort(key=lambda row: (-row["n_significant"], row["tissue"]))

    metabolic_rows = []
    for path in METABOLIC_FILES:
        item = json.loads(path.read_text())
        if item.get("status") != "success":
            continue
        metabolic_rows.append(
            {
                "tissue": item["tissue"],
                "gene_coverage_pct": float(item["gene_coverage_pct"]),
                "n_genes_mapped": int(item["n_genes_mapped"]),
                "n_model_genes": int(item["n_model_genes"]),
                "n_subsystems": int(item["n_subsystems"]),
            }
        )
    metabolic_rows.sort(key=lambda row: row["gene_coverage_pct"], reverse=True)

    drug = json.loads(DRUG_TARGETS.read_text())
    biomarkers = json.loads(BIOMARKERS.read_text())
    validation = [
        {"tissue": tissue, "auroc": float(result["auroc"])}
        for tissue, result in biomarkers["panel_validation"].items()
    ]
    validation.sort(key=lambda row: row["auroc"], reverse=True)

    top_genes = []
    seen = set()
    for row in biomarkers["panel"]:
        gene = row["human_symbol"] or row["gene"]
        gene = gene.upper()
        if gene in seen:
            continue
        seen.add(gene)
        top_genes.append(gene)
        if len(top_genes) == 4:
            break

    return {
        "immune": immune_rows,
        "tf": tf_rows,
        "metabolic": metabolic_rows,
        "drug": drug,
        "biomarkers": biomarkers,
        "validation": validation,
        "top_genes": top_genes,
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
    text(draw, (150, 76), "BIOLOGY / INTERPRETATION LAYER", F["kicker"], COLORS["teal"])
    x = 2280
    x = badge(draw, x, 66, "layer", "v5", COLORS["blue"])
    x = badge(draw, x, 66, "lanes", "5", COLORS["green"])
    x = badge(draw, x, 66, "tissues", "8", COLORS["amber"])
    badge(draw, x, 66, "readout", "interpretation", COLORS["violet"])


def draw_title(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 382), "Systems Biology Adds Interpretation", F["title"], COLORS["text"])
    paragraph(
        draw,
        (155, 493),
        "v5 turns benchmark readouts into immune context, TF activity, metabolic flux, target mapping, and biomarker panels.",
        F["subtitle"],
        2500,
        COLORS["muted"],
        10,
    )


def draw_lane_header(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, step: str, title: str, color: str) -> None:
    rounded(draw, (x, y, x + w, y + 1040), 34, COLORS["panel"], "#29374A", 2)
    rounded(draw, (x + 34, y + 34, x + 96, y + 96), 20, "#172335", color, 2)
    text(draw, (x + 65, y + 50), step, F["h3"], COLORS["text"], "ma")
    text(draw, (x + 118, y + 42), title, F["h2"], COLORS["text"])


def draw_mini_bar(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    label: str,
    value: float,
    max_value: float,
    color: str,
    width: int,
    value_label: str | None = None,
) -> None:
    text(draw, (x, y + 3), label, F["tiny_bold"], COLORS["text"])
    bar_x = x + 190
    draw.line((bar_x, y + 21, bar_x + width, y + 21), fill="#263244", width=15)
    fill_w = int(width * max(0.02, min(1.0, value / max_value)))
    draw.line((bar_x, y + 21, bar_x + fill_w, y + 21), fill=color, width=15)
    draw.ellipse((bar_x + fill_w - 11, y + 10, bar_x + fill_w + 11, y + 32), fill=color)
    label_text = value_label or f"{value:g}"
    text(draw, (bar_x + width + 20, y + 2), label_text, F["tiny_bold"], COLORS["text"])


def draw_immune_lane(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, rows: list[dict[str, object]]) -> None:
    color = LANE_COLORS["immune"]
    draw_lane_header(draw, x, y, w, "1", "Immune context", color)
    text(draw, (x + 54, y + 148), "6/14", F["stat"], color)
    text(draw, (x + 270, y + 158), "skin cell types", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "Deconvolution localizes the strongest immune remodeling signal and keeps tissue context visible.",
        F["body"],
        w - 108,
        COLORS["muted"],
        8,
    )
    rounded(draw, (x + 54, y + 420, x + w - 54, y + 730), 24, "#0D1520", "#2A394D", 1)
    text(draw, (x + 82, y + 446), "FDR<0.05 cell types", F["small_bold"], COLORS["text"])
    max_sig = max(row["n_significant"] for row in rows)
    for i, row in enumerate(rows[:4]):
        label = str(row["tissue"]).capitalize()
        value = int(row["n_significant"])
        total = int(row["n_cell_types"])
        draw_mini_bar(draw, x + 82, y + 500 + i * 50, label, value, max_sig, color, 250, f"{value}/{total}")
    rounded(draw, (x + 54, y + 820, x + w - 54, y + 950), 24, "#172335", "#2A394D", 1)
    text(draw, (x + 84, y + 846), "Readout", F["tiny_bold"], color)
    paragraph(draw, (x + 84, y + 882), "Skin is the high-context immune lane; kidney and thymus each add two FDR calls.", F["tiny"], w - 170, COLORS["muted"], 5)


def draw_tf_lane(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, rows: list[dict[str, object]]) -> None:
    color = LANE_COLORS["tf"]
    draw_lane_header(draw, x, y, w, "2", "TF activity", color)
    text(draw, (x + 54, y + 148), "241", F["stat"], color)
    text(draw, (x + 232, y + 158), "skin TFs", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "Regulatory activity adds a second interpretable axis beside immune composition.",
        F["body"],
        w - 108,
        COLORS["muted"],
        8,
    )
    rounded(draw, (x + 54, y + 420, x + w - 54, y + 730), 24, "#0D1520", "#2A394D", 1)
    text(draw, (x + 82, y + 446), "Significant TFs by tissue", F["small_bold"], COLORS["text"])
    max_sig = max(row["n_significant"] for row in rows)
    for i, row in enumerate(rows[:4]):
        label = str(row["tissue"]).capitalize()
        value = int(row["n_significant"])
        draw_mini_bar(draw, x + 82, y + 500 + i * 50, label, value, max_sig, color, 250, f"{value}")
    rounded(draw, (x + 54, y + 820, x + w - 54, y + 950), 24, "#172335", "#2A394D", 1)
    text(draw, (x + 84, y + 846), "Readout", F["tiny_bold"], color)
    paragraph(draw, (x + 84, y + 882), "Skin and thymus are nearly tied; kidney and liver remain strong secondary lanes.", F["tiny"], w - 170, COLORS["muted"], 5)


def draw_metabolic_lane(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, rows: list[dict[str, object]]) -> None:
    color = LANE_COLORS["metabolic"]
    draw_lane_header(draw, x, y, w, "3", "Metabolic flux", color)
    text(draw, (x + 54, y + 148), "6", F["stat"], color)
    text(draw, (x + 130, y + 158), "tissues solved", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "E-Flux and pFBA transform expression differences into pathway-scale flux hypotheses.",
        F["body"],
        w - 108,
        COLORS["muted"],
        8,
    )
    rounded(draw, (x + 54, y + 420, x + w - 54, y + 730), 24, "#0D1520", "#2A394D", 1)
    text(draw, (x + 82, y + 446), "Gene coverage", F["small_bold"], COLORS["text"])
    min_cov = min(row["gene_coverage_pct"] for row in rows)
    max_cov = max(row["gene_coverage_pct"] for row in rows)
    dot_y = y + 560
    axis_x0 = x + 100
    axis_x1 = x + w - 100
    draw.line((axis_x0, dot_y, axis_x1, dot_y), fill="#263244", width=4)
    for row in rows:
        cov = float(row["gene_coverage_pct"])
        px = axis_x0 + int((cov - min_cov) / (max_cov - min_cov) * (axis_x1 - axis_x0))
        draw.ellipse((px - 17, dot_y - 17, px + 17, dot_y + 17), fill=color)
    text(draw, (axis_x0, dot_y + 36), f"{min_cov:.1f}%", F["tiny_bold"], COLORS["muted"])
    text(draw, (axis_x1, dot_y + 36), f"{max_cov:.1f}%", F["tiny_bold"], COLORS["muted"], "ra")
    text(draw, (x + 82, y + 655), "103 subsystems per solved tissue", F["tiny_bold"], COLORS["text"])
    rounded(draw, (x + 54, y + 820, x + w - 54, y + 950), 24, "#172335", "#2A394D", 1)
    text(draw, (x + 84, y + 846), "Readout", F["tiny_bold"], color)
    paragraph(draw, (x + 84, y + 882), "Coverage spans 84.6-92.9%, supporting a compact cross-tissue flux layer.", F["tiny"], w - 170, COLORS["muted"], 5)


def draw_target_lane(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, drug: dict[str, object]) -> None:
    color = LANE_COLORS["target"]
    draw_lane_header(draw, x, y, w, "4", "Target mapping", color)
    text(draw, (x + 54, y + 148), "271", F["stat"], color)
    text(draw, (x + 236, y + 158), "DGIdb genes", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "Target mapping separates candidate genes, human-mapped symbols, and interaction evidence.",
        F["body"],
        w - 108,
        COLORS["muted"],
        8,
    )
    rounded(draw, (x + 54, y + 420, x + w - 54, y + 730), 24, "#0D1520", "#2A394D", 1)
    steps = [
        ("mouse", int(drug["n_target_genes_mouse"])),
        ("human", int(drug["n_mapped_to_human"])),
        ("DGIdb", int(drug["n_genes_with_dgidb_interactions"])),
    ]
    max_value = steps[0][1]
    bar_max_w = w - 255
    yy = y + 470
    for i, (label, value) in enumerate(steps):
        bar_w = int(bar_max_w * value / max_value)
        rounded(draw, (x + 82, yy + i * 78, x + 82 + bar_w, yy + i * 78 + 46), 18, rgba(color, 90 + i * 35), color, 2)
        text(draw, (x + 106, yy + 11 + i * 78), f"{value:,}", F["small_bold"], COLORS["text"])
        text(draw, (x + w - 96, yy + 11 + i * 78), label, F["tiny_bold"], COLORS["muted"], "ra")
    text(draw, (x + 82, y + 690), f"{int(drug['n_total_dgidb_interactions']):,} total interactions", F["tiny_bold"], COLORS["text"])
    rounded(draw, (x + 54, y + 820, x + w - 54, y + 950), 24, "#172335", "#2A394D", 1)
    text(draw, (x + 84, y + 846), "Readout", F["tiny_bold"], color)
    paragraph(draw, (x + 84, y + 882), "The lane makes downstream evidence auditable before follow-up prioritization.", F["tiny"], w - 170, COLORS["muted"], 5)


def draw_biomarker_lane(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    biomarkers: dict[str, object],
    validation: list[dict[str, object]],
    top_genes: list[str],
) -> None:
    color = LANE_COLORS["biomarker"]
    draw_lane_header(draw, x, y, w, "5", "Biomarker panel", color)
    text(draw, (x + 54, y + 148), f"{int(biomarkers['panel_size'])}", F["stat"], color)
    text(draw, (x + 154, y + 158), "genes", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "A compact panel summarizes recurrent signals and checks where the panel has the clearest readout.",
        F["body"],
        w - 108,
        COLORS["muted"],
        8,
    )
    rounded(draw, (x + 54, y + 420, x + w - 54, y + 730), 24, "#0D1520", "#2A394D", 1)
    text(draw, (x + 82, y + 446), "Panel validation AUROC", F["small_bold"], COLORS["text"])
    max_auc = 1.0
    for i, row in enumerate(validation[:4]):
        label = str(row["tissue"]).replace("gastrocnemius", "Gastro").capitalize()
        value = float(row["auroc"])
        draw_mini_bar(draw, x + 82, y + 500 + i * 50, label, value, max_auc, color, 250, f"{value:.3f}")
    gene_text = " / ".join(top_genes)
    text(draw, (x + 82, y + 690), gene_text, F["tiny_bold"], COLORS["muted"])
    rounded(draw, (x + 54, y + 820, x + w - 54, y + 950), 24, "#172335", "#2A394D", 1)
    text(draw, (x + 84, y + 846), "Readout", F["tiny_bold"], color)
    paragraph(draw, (x + 84, y + 882), "Gastrocnemius and liver give the strongest compact-panel validation readouts.", F["tiny"], w - 170, COLORS["muted"], 5)


def draw_lane_connectors(draw: ImageDraw.ImageDraw) -> None:
    y = 585
    xs = [150 + i * 714 for i in range(5)]
    for i in range(4):
        x0 = xs[i] + 650
        x1 = xs[i + 1] - 24
        draw.line((x0, y, x1, y), fill="#566275", width=4)
        draw.polygon([(x1, y), (x1 - 20, y - 11), (x1 - 20, y + 11)], fill="#566275")
    text(draw, (150, 566), "Benchmark readout flows into interpretation lanes", F["tiny_bold"], COLORS["muted"])


def draw_reader_rule(draw: ImageDraw.ImageDraw) -> None:
    box = (150, 1695, 3690, 1828)
    rounded(draw, box, 30, "#111A27", "#2A394D", 2)
    text(draw, (196, 1736), "Layer decoder", F["h3"], COLORS["teal"])
    steps = [
        ("Benchmark signal", "model score + task context"),
        ("Biology lanes", "immune, TF, flux, target, panel"),
        ("Cross-lane convergence", "shared tissues and genes"),
        ("Follow-up priorities", "ranked questions"),
    ]
    x = 700
    for i, (head, body) in enumerate(steps):
        color = [COLORS["blue"], COLORS["violet"], COLORS["amber"], COLORS["green"]][i]
        rounded(draw, (x, 1718, x + 620, 1808), 24, "#172335", color, 2)
        text(draw, (x + 28, 1733), head, F["small_bold"], COLORS["text"])
        text(draw, (x + 28, 1768), body, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            ax = x + 640
            draw.line((ax, 1763, ax + 82, 1763), fill="#6F7E90", width=4)
            draw.polygon([(ax + 82, 1763), (ax + 62, 1752), (ax + 62, 1774)], fill="#6F7E90")
        x += 730


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    source = "Takeaway: systems biology turns benchmark readouts into immune, TF, metabolic, target, and biomarker lanes."
    scope = "Next: immune/TF activity and target/biomarker slides unpack the highest-signal lanes."
    paragraph(draw, (150, 1888), source, F["small"], 3440, COLORS["muted"], 7)
    paragraph(draw, (150, 1993), scope, F["small_bold"], 3440, COLORS["text"], 7)


def main() -> None:
    data = load_data()
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)
    draw_background(draw)
    draw_header(draw)
    draw_title(draw)
    draw_lane_connectors(draw)

    panel_y = 625
    panel_w = 660
    gap = 54
    xs = [150 + i * (panel_w + gap) for i in range(5)]
    draw_immune_lane(draw, xs[0], panel_y, panel_w, data["immune"])
    draw_tf_lane(draw, xs[1], panel_y, panel_w, data["tf"])
    draw_metabolic_lane(draw, xs[2], panel_y, panel_w, data["metabolic"])
    draw_target_lane(draw, xs[3], panel_y, panel_w, data["drug"])
    draw_biomarker_lane(draw, xs[4], panel_y, panel_w, data["biomarkers"], data["validation"], data["top_genes"])
    draw_reader_rule(draw)
    draw_footer(draw)

    canvas.save(OUT_PATH, quality=95)
    gray = ImageOps.grayscale(canvas).convert("RGB")
    gray.save(GRAY_PATH, quality=95)

    stat = ImageStat.Stat(gray)
    manifest = {
        "slide": 37,
        "title": "Systems Biology Adds Interpretation",
        "readout": "v5 turns benchmark readouts into five biology interpretation lanes.",
        "asset": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "sources": [str(path.relative_to(ROOT)) for path in [*IMMUNE_FILES, *TF_FILES, *METABOLIC_FILES, DRUG_TARGETS, BIOMARKERS, *SOURCE_FIGS]],
        "data": {
            "immune_top": data["immune"][:4],
            "tf_top": data["tf"][:4],
            "metabolic_tissues": len(data["metabolic"]),
            "metabolic_coverage_range": [
                min(row["gene_coverage_pct"] for row in data["metabolic"]),
                max(row["gene_coverage_pct"] for row in data["metabolic"]),
            ],
            "metabolic_subsystems": data["metabolic"][0]["n_subsystems"],
            "target_genes_mouse": data["drug"]["n_target_genes_mouse"],
            "target_genes_human": data["drug"]["n_mapped_to_human"],
            "dgidb_genes": data["drug"]["n_genes_with_dgidb_interactions"],
            "dgidb_interactions": data["drug"]["n_total_dgidb_interactions"],
            "panel_size": data["biomarkers"]["panel_size"],
            "panel_validation_top4": data["validation"][:4],
            "panel_top_genes": data["top_genes"],
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

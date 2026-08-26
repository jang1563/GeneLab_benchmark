#!/usr/bin/env python3
"""Build slide 41 asset: translation details shape the readout."""

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
OUT_DIR = ASSET_ROOT / "translation_details"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "translation_details_define_claim_scope_premium.png"
GRAY_PATH = OUT_DIR / "translation_details_define_claim_scope_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "translation_details_manifest.json"

V6_EVAL = ROOT / "v6" / "evaluation"
SOURCES = {
    "gene": V6_EVAL / "V6_A_gene_conservation.json",
    "pathway": V6_EVAL / "V6_B_pathway_conservation.json",
    "transfer": V6_EVAL / "V6_C_cross_species_transfer.json",
    "biomarker": V6_EVAL / "V6_D_biomarker_validation.json",
    "tf": V6_EVAL / "V6_E_tf_conservation.json",
    "target": V6_EVAL / "V6_F_drug_target_validation.json",
}
SOURCE_FIG = ROOT / "v6" / "figures" / "html" / "v6_FigS8_translation_detail.html"

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
    "eye": "#5FD3C4",
    "thymus": "#B39DFF",
    "kidney": "#8BD17C",
    "skin": "#E17882",
    "liver": "#F4C26B",
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
    return tissue.capitalize()


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_data() -> dict[str, object]:
    gene = load_json(SOURCES["gene"])
    pathway = load_json(SOURCES["pathway"])
    transfer = load_json(SOURCES["transfer"])
    biomarker = load_json(SOURCES["biomarker"])
    tf = load_json(SOURCES["tf"])
    target = load_json(SOURCES["target"])

    common_genes = [int(row["n_common_genes"]) for row in transfer["analyses"]["pre_vs_post"].values()]
    gene_sets = [
        ("SHAP union", gene["enrichment_results"]["shap_union_top100"]),
        ("WGCNA hubs", gene["enrichment_results"]["wgcna_hubs"]),
        ("Panel genes", gene["enrichment_results"]["biomarker_panel"]),
    ]
    gene_rows = [
        {
            "label": label,
            "overlap": int(row["n_overlap_drr"]),
            "expected": float(row["expected_overlap"]),
            "fold": float(row["fold_enrichment"]),
        }
        for label, row in gene_sets
    ]
    pathway_rows = [
        {
            "tissue": tissue,
            "rho": float(pathway["per_tissue_correlations"][tissue]["spearman_rho"]),
        }
        for tissue in ["eye", "thymus", "kidney", "skin", "liver"]
    ]
    detected_genes = []
    seen = set()
    for row in biomarker["gene_results"]:
        symbol = (row.get("human_gene_matched") or row.get("human_symbol_stored") or "").upper()
        if row.get("detected_in_cfrna") and symbol and symbol not in seen:
            seen.add(symbol)
            detected_genes.append(symbol)
        if len(detected_genes) >= 6:
            break
    return {
        "gene": gene,
        "pathway": pathway,
        "transfer": transfer,
        "biomarker": biomarker,
        "tf": tf,
        "target": target,
        "common_min": min(common_genes),
        "common_max": max(common_genes),
        "gene_rows": gene_rows,
        "pathway_rows": pathway_rows,
        "detected_genes": detected_genes,
    }


def background() -> Image.Image:
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image)
    for y in range(H):
        t = y / H
        draw.line((0, y, W, y), fill=blend(COLORS["bg"], COLORS["bg2"], t))
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 25), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 22), width=1)
    draw.rectangle((0, 0, W, 260), fill=COLORS["header"])
    draw.rectangle((0, 1900, W, H), fill=COLORS["footer"])
    return image


def draw_decoder_strip(draw: ImageDraw.ImageDraw) -> None:
    items = [
        ("Ortholog map", "align mouse genes to human gene symbols", COLORS["blue"]),
        ("cfRNA", "human blood transcript readout", COLORS["green"]),
        ("NES", "pathway-level ranked enrichment score", COLORS["teal"]),
        ("Priority tier", "target-class follow-up queue", COLORS["amber"]),
    ]
    x, y = 120, 520
    gap = 32
    w = (3600 - gap * 3) // 4
    for title, desc, color in items:
        rounded(draw, (x, y, x + w, y + 128), 22, COLORS["panel"], "#26364B", 2)
        draw.ellipse((x + 28, y + 32, x + 84, y + 88), fill=color)
        text(draw, (x + 56, y + 60), title[0], F["small_bold"], COLORS["ink"], "mm")
        text(draw, (x + 108, y + 30), title, F["h3"], COLORS["text"])
        paragraph(draw, (x + 108, y + 71), desc, F["tiny"], w - 142, COLORS["muted"], 5)
        x += w + gap


def draw_metric_pill(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str, width: int) -> None:
    rounded(draw, (x, y, x + width, y + 58), 16, COLORS["panel3"], "#2A394D", 1)
    text(draw, (x + 18, y + 11), label, F["micro_bold"], COLORS["dim"])
    text(draw, (x + width - 18, y + 14), value, F["tiny_bold"], color, "ra")


def draw_dots(draw: ImageDraw.ImageDraw, x: int, y: int, n_total: int, n_fill: int, color: str) -> None:
    for i in range(n_total):
        row, col = divmod(i, 10)
        cx = x + col * 31
        cy = y + row * 31
        fill = color if i < n_fill else "#253244"
        outline = color if i < n_fill else "#37465A"
        draw.ellipse((cx, cy, cx + 18, cy + 18), fill=fill, outline=outline, width=2)


def draw_gene_mini_bars(draw: ImageDraw.ImageDraw, x: int, y: int, rows: list[dict[str, object]], width: int) -> None:
    axis_x = x + 170
    bar_w = width - 250
    for i, row in enumerate(rows):
        yy = y + i * 36
        color = [COLORS["teal"], COLORS["violet"], COLORS["amber"]][i]
        text(draw, (x, yy), str(row["label"]), F["micro_bold"], COLORS["text"])
        draw.line((axis_x, yy + 10, axis_x + bar_w, yy + 10), fill="#243246", width=10)
        fill_w = max(2, int(bar_w * min(float(row["fold"]), 1.15) / 1.15))
        draw.line((axis_x, yy + 10, axis_x + fill_w, yy + 10), fill=color, width=10)
        text(draw, (axis_x + bar_w + 18, yy - 3), f"{row['overlap']} / {row['expected']:.1f}", F["micro_bold"], color)


def draw_pathway_mini_axis(draw: ImageDraw.ImageDraw, x: int, y: int, rows: list[dict[str, object]], width: int) -> None:
    label_w = 112
    axis_mid = x + label_w + (width - label_w - 80) // 2
    scale_w = (width - label_w - 122) // 2
    draw.line((axis_mid, y - 8, axis_mid, y + 142), fill=COLORS["dim"], width=2)
    for i, row in enumerate(rows):
        yy = y + i * 31
        tissue = str(row["tissue"])
        rho = float(row["rho"])
        color = TISSUE_COLORS[tissue]
        text(draw, (x, yy - 4), tissue_label(tissue), F["micro_bold"], COLORS["text"])
        end_x = axis_mid + int(rho / 0.5 * scale_w)
        draw.line((axis_mid, yy + 7, end_x, yy + 7), fill=color, width=12)
        text(draw, (end_x + 14 if rho >= 0 else end_x - 14, yy - 6), f"{rho:+.3f}", F["micro_bold"], color, None if rho >= 0 else "ra")


def draw_tier_stack(draw: ImageDraw.ImageDraw, x: int, y: int, tiers: dict[str, int], width: int) -> None:
    colors = {"A": COLORS["teal"], "B": COLORS["green"], "C": COLORS["blue"], "D": COLORS["dim"]}
    fixed = {"A": 48, "B": 58, "D": 48}
    c_width = width - fixed["A"] - fixed["B"] - fixed["D"] - 18
    cursor = x
    for tier, segment_w in [("A", fixed["A"]), ("B", fixed["B"]), ("C", c_width), ("D", fixed["D"])]:
        color = colors[tier]
        rounded(draw, (cursor, y, cursor + segment_w, y + 42), 14, color, None, 0)
        if segment_w > 82:
            text(draw, (cursor + 16, y + 11), f"{tier}: {tiers[tier]}", F["micro_bold"], COLORS["ink"])
        cursor += segment_w + 6


def draw_gate_row(
    draw: ImageDraw.ImageDraw,
    y: int,
    number: str,
    title: str,
    question: str,
    color: str,
    draw_metrics,
) -> None:
    x1, x2 = 120, 2680
    rounded(draw, (x1, y, x2, y + 166), 22, COLORS["panel"], "#26364B", 2)
    draw.ellipse((x1 + 30, y + 45, x1 + 92, y + 107), fill=color)
    text(draw, (x1 + 61, y + 65), number, F["tiny_bold"], COLORS["ink"], "mm")
    text(draw, (x1 + 118, y + 34), title, F["h3"], COLORS["text"])
    paragraph(draw, (x1 + 118, y + 80), question, F["tiny"], 560, COLORS["muted"], 5)
    draw.line((x1 + 760, y + 22, x1 + 760, y + 144), fill="#2B3A4E", width=2)
    draw_metrics(x1 + 798, y + 30)


def draw_gate_ledger(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    rows_y = 720
    gap = 18
    row_h = 166

    def metrics_1(x: int, y: int) -> None:
        draw_metric_pill(draw, x, y, "ortholog map", f"{data['transfer']['ortholog_map_size']:,}", COLORS["blue"], 330)
        draw_metric_pill(draw, x + 360, y, "human cfRNA", f"{data['transfer']['cfRNA_genes']} genes", COLORS["green"], 330)
        draw_metric_pill(draw, x + 720, y, "classifier common", f"{data['common_min']}-{data['common_max']} genes", COLORS["teal"], 410)
        paragraph(draw, (x, y + 78), "This gate defines the measurable cross-species gene space.", F["tiny"], 1040, COLORS["muted"], 5)

    def metrics_2(x: int, y: int) -> None:
        draw_metric_pill(draw, x, y, "human DRR universe", f"{data['gene']['drr_in_universe']} / {data['gene']['universe_size']:,}", COLORS["blue"], 420)
        draw_gene_mini_bars(draw, x + 465, y + 6, data["gene_rows"], 640)

    def metrics_3(x: int, y: int) -> None:
        draw_metric_pill(draw, x, y, "matched pathways", f"{data['pathway']['n_human_pathways']} Hallmark", COLORS["teal"], 360)
        draw_metric_pill(draw, x + 390, y, "mean rank rho", f"{float(data['pathway']['mean_rho']):.3f}", COLORS["teal"], 310)
        draw_pathway_mini_axis(draw, x + 735, y + 8, data["pathway_rows"], 430)

    def metrics_4(x: int, y: int) -> None:
        biomarker = data["biomarker"]
        draw_metric_pill(draw, x, y, "panel detected", f"{biomarker['n_detected_in_cfrna']}/{biomarker['panel_size']}", COLORS["green"], 300)
        draw_metric_pill(draw, x + 326, y, "DE threshold", f"{biomarker['n_de_fdr05']}/{biomarker['panel_size']}", COLORS["amber"], 300)
        draw_dots(draw, x + 655, y + 7, int(biomarker["panel_size"]), int(biomarker["n_detected_in_cfrna"]), COLORS["green"])
        genes = "  ".join(data["detected_genes"][:6])
        text(draw, (x, y + 87), genes, F["tiny_bold"], COLORS["text"])

    def metrics_5(x: int, y: int) -> None:
        tf = data["tf"]
        target = data["target"]
        draw_metric_pill(draw, x, y, "TF lane", f"{tf['n_human_tfs']} TFs | {tf['n_sig_human_tfs']} signals", COLORS["violet"], 470)
        draw_metric_pill(draw, x + 500, y, "target genes", f"{target['n_in_cfrna']}/{target['n_drug_target_genes']}", COLORS["green"], 300)
        draw_tier_stack(draw, x + 830, y + 8, {k: int(v) for k, v in target["tier_counts"].items()}, 315)
        text(draw, (x, y + 86), f"TF mean rho {float(tf['mean_rho']):.3f} | target tiers A {target['tier_counts']['A']}, B {target['tier_counts']['B']}, C {target['tier_counts']['C']}, D {target['tier_counts']['D']}", F["tiny_bold"], COLORS["muted"])

    gate_specs = [
        ("01", "Ortholog + cfRNA Availability", "Which mouse genes can be named and measured on the human side?", COLORS["blue"], metrics_1),
        ("02", "Gene-Set Overlap", "Does a mouse-selected gene set land in the human DRR universe?", COLORS["amber"], metrics_2),
        ("03", "Pathway NES Conservation", "Do pathway ranks move together across matched tissues?", COLORS["teal"], metrics_3),
        ("04", "Biomarker Panel Readback", "Are the proposed panel genes visible in astronaut cfRNA?", COLORS["green"], metrics_4),
        ("05", "TF + Target-Class Queue", "Which regulatory and target layers create follow-up packets?", COLORS["violet"], metrics_5),
    ]
    for i, spec in enumerate(gate_specs):
        draw_gate_row(draw, rows_y + i * (row_h + gap), *spec)


def draw_claim_decoder(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 2740, 720, 3720, 1622
    rounded(draw, (x1, y1, x2, y2), 26, COLORS["panel"], "#26364B", 2)
    text(draw, (x1 + 42, y1 + 36), "Layer decoder", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 42, y1 + 92),
        "Translation strength is assigned by layer, so each measurement answers a narrower question.",
        F["body"],
        x2 - x1 - 84,
        COLORS["muted"],
        8,
    )
    cards = [
        ("Map", "measurable gene space", COLORS["blue"]),
        ("Overlap", "gene-level context", COLORS["amber"]),
        ("Conserve", "pathway bridge", COLORS["teal"]),
        ("Read back", "human observability", COLORS["green"]),
        ("Queue", "follow-up prioritization", COLORS["violet"]),
    ]
    y = y1 + 218
    for i, (title, desc, color) in enumerate(cards):
        rounded(draw, (x1 + 42, y, x2 - 42, y + 92), 20, COLORS["panel2"], color, 2)
        draw.ellipse((x1 + 68, y + 25, x1 + 112, y + 69), fill=color)
        text(draw, (x1 + 90, y + 37), f"{i + 1}", F["micro_bold"], COLORS["ink"], "mm")
        text(draw, (x1 + 138, y + 18), title, F["small_bold"], COLORS["text"])
        text(draw, (x1 + 138, y + 53), desc, F["tiny"], COLORS["muted"])
        if i < len(cards) - 1:
            ax = x1 + 90
            draw.line((ax, y + 94, ax, y + 108), fill=COLORS["dim"], width=3)
            draw.polygon([(ax - 8, y + 108), (ax + 8, y + 108), (ax, y + 120)], fill=COLORS["dim"])
        y += 112

    callout_y = y2 - 100
    rounded(draw, (x1 + 42, callout_y, x2 - 42, callout_y + 70), 20, "#102133", "#33506B", 2)
    text(draw, (x1 + 70, callout_y + 14), "Interpretation move", F["small_bold"], COLORS["teal"])
    text(draw, (x1 + 70, callout_y + 45), "keep layer, statistic, and interpretation aligned", F["tiny_bold"], COLORS["text"])


def draw_reader_flow(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1668, 3720, 1846
    rounded(draw, (x1, y1, x2, y2), 28, "#101927", "#26364B", 2)
    text(draw, (x1 + 38, y1 + 32), "Slide role", F["small_bold"], COLORS["muted"])
    text(draw, (x1 + 260, y1 + 35), "Slide 40 gives the overview", F["small_bold"], COLORS["blue"])
    draw.line((x1 + 805, y1 + 55, x1 + 920, y1 + 55), fill=COLORS["dim"], width=4)
    draw.polygon([(x1 + 920, y1 + 44), (x1 + 920, y1 + 66), (x1 + 945, y1 + 55)], fill=COLORS["dim"])
    text(draw, (x1 + 985, y1 + 35), "Slide 41 shows the layers", F["small_bold"], COLORS["teal"])
    draw.line((x1 + 1505, y1 + 55, x1 + 1620, y1 + 55), fill=COLORS["dim"], width=4)
    draw.polygon([(x1 + 1620, y1 + 44), (x1 + 1620, y1 + 66), (x1 + 1645, y1 + 55)], fill=COLORS["dim"])
    text(draw, (x1 + 1685, y1 + 35), "Slide 42 turns layers into synthesis", F["small_bold"], COLORS["green"])
    paragraph(
        draw,
        (x1 + 260, y1 + 96),
        "The detailed deck should let the audience separate measurable transfer, pathway conservation, human readback, and prioritization before the final interpretation slide.",
        F["small"],
        3180,
        COLORS["muted"],
        8,
    )


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, "#0F1824", "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "Use layer-specific measurements to separate measurable transfer, pathway conservation, human readback, and follow-up prioritization.",
        F["small"],
        2260,
        COLORS["muted"],
        8,
    )
    text(draw, (3580, 1960), "41", F["h2"], COLORS["teal"], "ra")
    text(draw, (3580, 2010), "translation detail", F["tiny_bold"], COLORS["muted"], "ra")


def build() -> None:
    data = load_data()
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 41 | ACT 4 | TRANSLATION", F["kicker"], COLORS["teal"])
    bx = 2170
    bx = badge(draw, bx, 56, "DETAIL", "v6 FigS8", COLORS["teal"])
    bx = badge(draw, bx, 56, "ORTHOLOGS", f"{data['transfer']['ortholog_map_size']:,}", COLORS["blue"])
    bx = badge(draw, bx, 56, "cfRNA", f"{data['transfer']['cfRNA_genes']} genes", COLORS["green"])
    badge(draw, bx, 56, "LAYERS", "5", COLORS["amber"])

    text(draw, (120, 330), "Translation Details Shape The Readout", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "Each mouse-to-human readout passes through a different measurement layer before it becomes interpretable.",
        F["subtitle"],
        2700,
        COLORS["muted"],
        8,
    )

    draw_decoder_strip(draw)
    draw_gate_ledger(draw, data)
    draw_claim_decoder(draw)
    draw_reader_flow(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    gray_rgb = ImageOps.colorize(gray, black="#080D14", white="#F3F7FC")
    gray_rgb.save(GRAY_PATH, quality=96)

    stat = ImageStat.Stat(gray)
    manifest = {
        "title": "Translation Details Shape The Readout",
        "claim": "Layer-specific measurements separate measurable transfer, pathway conservation, human readback, and follow-up prioritization.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "source_figure": str(SOURCE_FIG.relative_to(ROOT)),
        "source_json": {key: str(path.relative_to(ROOT)) for key, path in SOURCES.items()},
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": {
            "ortholog_map_size": data["transfer"]["ortholog_map_size"],
            "cfrna_genes": data["transfer"]["cfRNA_genes"],
            "drr_genes": data["gene"]["drr_in_universe"],
            "pathway_mean_rho": data["pathway"]["mean_rho"],
            "biomarker_detected": data["biomarker"]["n_detected_in_cfrna"],
            "biomarker_panel": data["biomarker"]["panel_size"],
            "target_detected": data["target"]["n_in_cfrna"],
            "target_genes": data["target"]["n_drug_target_genes"],
            "tf_mean_rho": data["tf"]["mean_rho"],
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    build()

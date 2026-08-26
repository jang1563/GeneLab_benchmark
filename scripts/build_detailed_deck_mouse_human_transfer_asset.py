#!/usr/bin/env python3
"""Build slide 40 asset: pathways carry mouse-to-human transfer."""

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
OUT_DIR = ASSET_ROOT / "mouse_human_transfer"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "pathways_carry_mouse_to_human_transfer_premium.png"
GRAY_PATH = OUT_DIR / "pathways_carry_mouse_to_human_transfer_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "mouse_human_transfer_manifest.json"

V6_EVAL = ROOT / "v6" / "evaluation"
SOURCES = {
    "gene": V6_EVAL / "V6_A_gene_conservation.json",
    "pathway": V6_EVAL / "V6_B_pathway_conservation.json",
    "transfer": V6_EVAL / "V6_C_cross_species_transfer.json",
    "biomarker": V6_EVAL / "V6_D_biomarker_validation.json",
    "tf": V6_EVAL / "V6_E_tf_conservation.json",
    "target": V6_EVAL / "V6_F_drug_target_validation.json",
}
SOURCE_FIG = ROOT / "v6" / "figures" / "html" / "v6_Fig10_human_translation.html"

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
    "gastrocnemius": "#E69F00",
    "eye": "#5FD3C4",
    "kidney": "#8BD17C",
    "thymus": "#B39DFF",
    "skin": "#E17882",
    "liver": "#F4C26B",
    "lung": "#A8B4C4",
    "colon": "#66A6E8",
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


def fmt_p(value: float) -> str:
    if value < 0.001:
        return f"p={value:.1e}"
    return f"p={value:.3f}"


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_data() -> dict[str, object]:
    gene = load_json(SOURCES["gene"])
    pathway = load_json(SOURCES["pathway"])
    transfer = load_json(SOURCES["transfer"])
    biomarker = load_json(SOURCES["biomarker"])
    tf = load_json(SOURCES["tf"])
    target = load_json(SOURCES["target"])

    enrich = gene["enrichment_results"]
    gene_rows = [
        {
            "label": "SHAP union",
            "fold": float(enrich["shap_union_top100"]["fold_enrichment"]),
            "overlap": int(enrich["shap_union_top100"]["n_overlap_drr"]),
            "expected": float(enrich["shap_union_top100"]["expected_overlap"]),
        },
        {
            "label": "WGCNA hubs",
            "fold": float(enrich["wgcna_hubs"]["fold_enrichment"]),
            "overlap": int(enrich["wgcna_hubs"]["n_overlap_drr"]),
            "expected": float(enrich["wgcna_hubs"]["expected_overlap"]),
        },
        {
            "label": "Panel genes",
            "fold": float(enrich["biomarker_panel"]["fold_enrichment"]),
            "overlap": int(enrich["biomarker_panel"]["n_overlap_drr"]),
            "expected": float(enrich["biomarker_panel"]["expected_overlap"]),
        },
    ]

    transfer_order = ["gastrocnemius", "eye", "liver", "thymus", "skin", "kidney"]
    transfer_rows = []
    for tissue in transfer_order:
        row = transfer["analyses"]["pre_vs_post"][tissue]
        transfer_rows.append(
            {
                "tissue": tissue,
                "auroc": float(row["transfer_auroc"]),
                "mouse": float(row["mouse_internal_auroc"]),
                "n_common": int(row["n_common_genes"]),
                "p": float(row["perm_p"]),
            }
        )

    pathway_rows = []
    pathway_order = ["eye", "thymus", "kidney", "skin", "liver"]
    for tissue in pathway_order:
        row = pathway["per_tissue_correlations"][tissue]
        pathway_rows.append(
            {
                "tissue": tissue,
                "rho": float(row["spearman_rho"]),
                "p": float(row["p_value"]),
                "n": int(row["n_matched_pathways"]),
                "significant": bool(row["significant"]),
            }
        )

    top_human_nes = sorted(pathway["human_nes"].items(), key=lambda item: abs(float(item[1])), reverse=True)[:4]
    top_human_labels = [label.replace("HALLMARK_", "").replace("_", " ").title() for label, _ in top_human_nes]

    target_tiers = {key: int(value) for key, value in target["tier_counts"].items()}
    tier_a = [row["gene"] for row in target["tier_A_validated"]]
    tier_b = [row["gene"] for row in target["tier_B_promising"][:4]]

    return {
        "gene": gene,
        "gene_rows": gene_rows,
        "transfer": transfer,
        "transfer_rows": transfer_rows,
        "pathway": pathway,
        "pathway_rows": pathway_rows,
        "top_human_labels": top_human_labels,
        "biomarker": biomarker,
        "tf": tf,
        "target": target,
        "target_tiers": target_tiers,
        "target_genes": tier_a + tier_b,
    }


def background() -> Image.Image:
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image)
    for y in range(H):
        t = y / H
        color = blend(COLORS["bg"], COLORS["bg2"], t)
        draw.line((0, y, W, y), fill=color)
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 25), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 22), width=1)
    draw.rectangle((0, 0, W, 260), fill=COLORS["header"])
    draw.rectangle((0, 1900, W, H), fill=COLORS["footer"])
    return image


def panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], title: str, stat: str, stat_label: str, color: str) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 24, COLORS["panel"], "#253247", 2)
    rounded(draw, (x1 + 22, y1 + 22, x2 - 22, y1 + 128), 20, COLORS["panel2"], None, 0)
    text(draw, (x1 + 46, y1 + 43), title, F["h2"], COLORS["text"])
    text(draw, (x2 - 46, y1 + 32), stat, F["stat2"], color, "ra")
    text(draw, (x2 - 46, y1 + 92), stat_label, F["tiny_bold"], COLORS["muted"], "ra")


def draw_metric_chip(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str, width: int) -> None:
    rounded(draw, (x, y, x + width, y + 82), 18, COLORS["panel3"], "#29384D", 1)
    text(draw, (x + 20, y + 16), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 20, y + 44), value, F["small_bold"], color)


def draw_gene_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x1, y1, x2, _ = box
    panel(draw, box, "Gene-Level Transfer Lane", str(data["gene"]["drr_in_universe"]), "human DRR genes", COLORS["blue"])
    y = y1 + 156
    paragraph(
        draw,
        (x1 + 46, y),
        "Gene overlap and ortholog classifiers define the direct-transfer readout before higher-level biology layers.",
        F["body"],
        x2 - x1 - 92,
        COLORS["muted"],
        8,
    )

    chart_x = x1 + 48
    chart_y = y1 + 300
    chart_w = x2 - x1 - 96
    text(draw, (chart_x, chart_y - 44), "DRR enrichment fold", F["h3"], COLORS["text"])
    axis_x = chart_x + 245
    bar_w = chart_w - 360
    draw.line((axis_x, chart_y, axis_x + bar_w, chart_y), fill=COLORS["grid"], width=2)
    one_x = axis_x + int(bar_w * 1.0 / 1.15)
    draw.line((one_x, chart_y - 18, one_x, chart_y + 198), fill=COLORS["dim"], width=2)
    text(draw, (one_x, chart_y - 44), "1.0x guide", F["micro_bold"], COLORS["dim"], "mm")
    for i, row in enumerate(data["gene_rows"]):
        yy = chart_y + 36 + i * 62
        text(draw, (chart_x, yy + 7), row["label"], F["small_bold"], COLORS["text"])
        draw.line((axis_x, yy + 22, axis_x + bar_w, yy + 22), fill="#1D2A3B", width=18)
        width = max(2, int(bar_w * min(float(row["fold"]), 1.15) / 1.15))
        color = [COLORS["teal"], COLORS["violet"], COLORS["amber"]][i]
        draw.line((axis_x, yy + 22, axis_x + width, yy + 22), fill=color, width=18)
        text(draw, (axis_x + width + 18, yy + 9), f"{row['fold']:.2f}x", F["small_bold"], color)
        text(draw, (x2 - 54, yy + 11), f"{row['overlap']} / exp {row['expected']:.1f}", F["tiny"], COLORS["muted"], "ra")

    lower_y = y1 + 570
    text(draw, (chart_x, lower_y - 40), "Ortholog classifier AUROC", F["h3"], COLORS["text"])
    text(draw, (x2 - 54, lower_y - 36), "389-398 common genes", F["tiny"], COLORS["muted"], "ra")
    auroc_axis_x = chart_x + 210
    auroc_bar_w = chart_w - 372
    top = lower_y + 6
    draw.line((auroc_axis_x + int(auroc_bar_w * 0.5), top - 22, auroc_axis_x + int(auroc_bar_w * 0.5), top + 310), fill=COLORS["dim"], width=2)
    text(draw, (auroc_axis_x + int(auroc_bar_w * 0.5), top - 47), "0.50 guide", F["micro_bold"], COLORS["dim"], "mm")
    for i, row in enumerate(data["transfer_rows"]):
        yy = top + i * 50
        tissue = row["tissue"]
        color = TISSUE_COLORS[tissue]
        text(draw, (chart_x, yy + 5), tissue_label(tissue), F["tiny_bold"], COLORS["text"])
        draw.line((auroc_axis_x, yy + 18, auroc_axis_x + auroc_bar_w, yy + 18), fill="#1D2A3B", width=15)
        width = max(3, int(auroc_bar_w * float(row["auroc"])))
        draw.line((auroc_axis_x, yy + 18, auroc_axis_x + width, yy + 18), fill=color, width=15)
        text(draw, (x2 - 54, yy + 3), f"{row['auroc']:.3f}", F["tiny_bold"], color, "ra")


def draw_pathway_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x1, y1, x2, _ = box
    pathway = data["pathway"]
    panel(draw, box, "Pathway Bridge", f"{float(pathway['mean_rho']):.3f}", "mean Spearman rho", COLORS["teal"])
    y = y1 + 156
    paragraph(
        draw,
        (x1 + 46, y),
        "Fifty Hallmark pathway NES profiles provide a direction-aware transfer surface across matched tissues.",
        F["body"],
        x2 - x1 - 92,
        COLORS["muted"],
        8,
    )

    chart_x = x1 + 64
    chart_y = y1 + 330
    chart_w = x2 - x1 - 128
    axis_mid = chart_x + chart_w // 2
    scale_w = chart_w // 2 - 130
    text(draw, (chart_x, chart_y - 52), "Mouse vs human pathway rank correlation", F["h3"], COLORS["text"])
    draw.line((axis_mid - scale_w, chart_y - 6, axis_mid + scale_w, chart_y - 6), fill=COLORS["grid"], width=2)
    for tick in [-0.5, 0.0, 0.5]:
        tx = axis_mid + int(tick / 0.5 * scale_w)
        draw.line((tx, chart_y - 22, tx, chart_y + 350), fill=COLORS["grid"] if tick else COLORS["dim"], width=2)
        text(draw, (tx, chart_y - 48), f"{tick:.1f}", F["micro_bold"], COLORS["dim"], "mm")
    for i, row in enumerate(data["pathway_rows"]):
        yy = chart_y + 34 + i * 70
        tissue = row["tissue"]
        rho = float(row["rho"])
        color = TISSUE_COLORS[tissue]
        text(draw, (chart_x, yy + 2), tissue_label(tissue), F["small_bold"], COLORS["text"])
        end_x = axis_mid + int(rho / 0.5 * scale_w)
        draw.line((axis_mid, yy + 18, end_x, yy + 18), fill=color, width=22)
        draw.ellipse((end_x - 10, yy + 8, end_x + 10, yy + 28), fill=color)
        label_x = end_x + 18 if rho >= 0 else end_x - 18
        anchor = None if rho >= 0 else "ra"
        text(draw, (label_x, yy + 3), f"rho {rho:+.3f}", F["small_bold"], color, anchor)

    callout_y = y1 + 746
    rounded(draw, (x1 + 46, callout_y, x2 - 46, callout_y + 118), 20, COLORS["panel3"], "#29384D", 1)
    text(draw, (x1 + 72, callout_y + 23), "Top human NES signals", F["small_bold"], COLORS["teal"])
    labels = data["top_human_labels"]
    pill_x = x1 + 72
    for label in labels:
        short = label.replace("Interferon Alpha Response", "IFN-alpha").replace("Androgen Response", "Androgen")
        width = int(draw.textlength(short, font=F["tiny_bold"]) + 46)
        rounded(draw, (pill_x, callout_y + 66, pill_x + width, callout_y + 101), 16, "#1B2B3B", "#2D4056", 1)
        text(draw, (pill_x + 23, callout_y + 74), short, F["tiny_bold"], COLORS["text"])
        pill_x += width + 12


def draw_dots(draw: ImageDraw.ImageDraw, x: int, y: int, n_total: int, n_fill: int, color: str) -> None:
    for i in range(n_total):
        row, col = divmod(i, 10)
        cx = x + col * 34
        cy = y + row * 34
        fill = color if i < n_fill else "#263244"
        outline = color if i < n_fill else "#38475A"
        draw.ellipse((cx, cy, cx + 20, cy + 20), fill=fill, outline=outline, width=2)


def draw_target_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x1, y1, x2, _ = box
    target = data["target"]
    biomarker = data["biomarker"]
    panel(draw, box, "Human cfRNA Readback", f"{target['n_in_cfrna']}/{target['n_drug_target_genes']}", "target genes detected", COLORS["green"])
    y = y1 + 156
    paragraph(
        draw,
        (x1 + 46, y),
        "Panel genes and target-class tiers convert translation evidence into auditable follow-up packets.",
        F["body"],
        x2 - x1 - 92,
        COLORS["muted"],
        8,
    )

    chip_y = y1 + 306
    draw_metric_chip(draw, x1 + 46, chip_y, "Biomarker panel", f"{biomarker['n_detected_in_cfrna']}/{biomarker['panel_size']} detected", COLORS["green"], 360)
    draw_metric_chip(draw, x1 + 430, chip_y, "Panel DE calls", f"{biomarker['n_de_fdr05']}/{biomarker['panel_size']} at FDR<0.05", COLORS["amber"], 360)
    draw_dots(draw, x1 + 818, chip_y + 10, int(biomarker["panel_size"]), int(biomarker["n_detected_in_cfrna"]), COLORS["green"])

    tier_y = y1 + 480
    text(draw, (x1 + 46, tier_y - 42), "Target-class evidence tiers", F["h3"], COLORS["text"])
    tiers = data["target_tiers"]
    total = sum(tiers.values())
    tier_colors = {"A": COLORS["teal"], "B": COLORS["green"], "C": COLORS["blue"], "D": COLORS["dim"]}
    stack_x = x1 + 46
    stack_y = tier_y + 14
    stack_w = x2 - x1 - 92
    cursor = stack_x
    min_w = {"A": 54, "B": 64, "C": 0, "D": 50}
    remaining_w = stack_w
    segments = []
    for tier in ["A", "B", "D"]:
        seg_w = min_w[tier]
        segments.append((tier, seg_w))
        remaining_w -= seg_w
    segments.insert(2, ("C", max(180, remaining_w)))
    for tier, seg_w in segments:
        color = tier_colors[tier]
        rounded(draw, (cursor, stack_y, cursor + seg_w - 6, stack_y + 58), 16, color, None, 0)
        if seg_w > 70:
            text(draw, (cursor + 16, stack_y + 16), f"{tier}: {tiers[tier]}", F["tiny_bold"], COLORS["ink"])
        cursor += seg_w
    text(draw, (stack_x, stack_y + 76), f"{total} target genes across A/B/C/D evidence tiers", F["tiny"], COLORS["muted"])

    genes_y = y1 + 644
    text(draw, (x1 + 46, genes_y - 34), "Example target genes", F["small_bold"], COLORS["text"])
    gx = x1 + 46
    for gene in data["target_genes"][:7]:
        width = int(draw.textlength(gene, font=F["tiny_bold"]) + 42)
        rounded(draw, (gx, genes_y, gx + width, genes_y + 42), 16, COLORS["panel3"], "#33445C", 1)
        text(draw, (gx + 21, genes_y + 11), gene, F["tiny_bold"], COLORS["text"])
        gx += width + 12
        if gx > x2 - 220:
            break

    tf_y = y1 + 758
    rounded(draw, (x1 + 46, tf_y, x2 - 46, tf_y + 108), 20, COLORS["panel3"], "#29384D", 1)
    text(draw, (x1 + 72, tf_y + 22), "TF activity lane", F["small_bold"], COLORS["violet"])
    tf = data["tf"]
    tf_label = f"{tf['n_human_tfs']} human TFs | {tf['n_sig_human_tfs']} human signals | mean rho {float(tf['mean_rho']):.3f}"
    text(draw, (x1 + 72, tf_y + 64), tf_label, F["tiny_bold"], COLORS["text"])


def draw_reader_flow(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1694, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 28, "#101927", "#26364B", 2)
    text(draw, (x1 + 38, y1 + 28), "Takeaway", F["small_bold"], COLORS["muted"])
    steps = [
        ("01", "Gene lane", "overlap + AUROC"),
        ("02", "Pathway bridge", "NES rank transfer"),
        ("03", "cfRNA readback", "panel detection"),
        ("04", "Target queue", "tiered follow-up"),
    ]
    node_w = 690
    gap = 70
    start_x = x1 + 318
    node_y = y1 + 42
    for i, (num, title, desc) in enumerate(steps):
        nx = start_x + i * (node_w + gap)
        color = [COLORS["blue"], COLORS["teal"], COLORS["green"], COLORS["amber"]][i]
        rounded(draw, (nx, node_y, nx + node_w, node_y + 82), 20, COLORS["panel2"], color, 2)
        draw.ellipse((nx + 22, node_y + 20, nx + 64, node_y + 62), fill=color)
        text(draw, (nx + 43, node_y + 31), num, F["micro_bold"], COLORS["ink"], "mm")
        text(draw, (nx + 86, node_y + 14), title, F["small_bold"], COLORS["text"])
        text(draw, (nx + 86, node_y + 48), desc, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            ax = nx + node_w + 14
            ay = node_y + 41
            draw.line((ax, ay, ax + gap - 28, ay), fill=COLORS["dim"], width=3)
            draw.polygon([(ax + gap - 28, ay - 9), (ax + gap - 28, ay + 9), (ax + gap - 12, ay)], fill=COLORS["dim"])


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, "#0F1824", "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "Mouse signals prioritize human follow-up through pathway and target-class evidence; direct-gene transfer remains one evidence lane.",
        F["small"],
        2080,
        COLORS["muted"],
        8,
    )
    text(draw, (3580, 1960), "Inputs", F["micro_bold"], COLORS["dim"], "ra")
    text(draw, (3580, 1992), "v6/evaluation V6_A-F JSON", F["tiny"], COLORS["muted"], "ra")
    text(draw, (3580, 2024), "v6 Fig10 human translation", F["tiny"], COLORS["muted"], "ra")


def build() -> None:
    data = load_data()
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 40 | ACT 4 | TRANSLATION", F["kicker"], COLORS["teal"])
    bx = 2290
    bx = badge(draw, bx, 56, "LAYER", "v6 Fig10", COLORS["teal"])
    bx = badge(draw, bx, 56, "PATHWAYS", f"{data['pathway']['n_human_pathways']} Hallmark", COLORS["blue"])
    bx = badge(draw, bx, 56, "cfRNA", f"{data['transfer']['cfRNA_genes']} genes", COLORS["green"])
    badge(draw, bx, 56, "TARGETS", f"{data['target']['n_in_cfrna']}/{data['target']['n_drug_target_genes']}", COLORS["amber"])

    text(draw, (120, 330), "Pathways Carry Mouse-To-Human Transfer", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "Direct gene signals, pathway NES, and human cfRNA readback separate the translation layers.",
        F["subtitle"],
        2550,
        COLORS["muted"],
        8,
    )

    panel_y = 610
    panel_h = 1048
    panel_w = 1160
    gap = 60
    draw_gene_panel(draw, (120, panel_y, 120 + panel_w, panel_y + panel_h), data)
    draw_pathway_panel(draw, (120 + panel_w + gap, panel_y, 120 + panel_w * 2 + gap, panel_y + panel_h), data)
    draw_target_panel(draw, (120 + panel_w * 2 + gap * 2, panel_y, 120 + panel_w * 3 + gap * 2, panel_y + panel_h), data)

    draw_reader_flow(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    gray_rgb = ImageOps.colorize(gray, black="#080D14", white="#F3F7FC")
    gray_rgb.save(GRAY_PATH, quality=96)

    stat = ImageStat.Stat(gray)
    manifest = {
        "title": "Pathways Carry Mouse-To-Human Transfer",
        "readout": "Translation reads as layered evidence: direct gene transfer, pathway NES, cfRNA readback, and target-class queues.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "source_figure": str(SOURCE_FIG.relative_to(ROOT)),
        "source_json": {key: str(path.relative_to(ROOT)) for key, path in SOURCES.items()},
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": {
            "drr_genes": data["gene"]["drr_in_universe"],
            "pathway_mean_rho": data["pathway"]["mean_rho"],
            "n_human_pathways": data["pathway"]["n_human_pathways"],
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

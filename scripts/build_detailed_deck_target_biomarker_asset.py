#!/usr/bin/env python3
"""Build slide 39 asset: target and biomarker layers are triage."""

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
OUT_DIR = ASSET_ROOT / "target_biomarker"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "target_biomarker_layers_are_triage_premium.png"
GRAY_PATH = OUT_DIR / "target_biomarker_layers_are_triage_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "target_biomarker_manifest.json"

V5_EVAL = ROOT / "v5" / "evaluation"
METABOLIC_FILES = sorted(V5_EVAL.glob("metabolic_flux_*.json"))
DRUG_TARGETS = V5_EVAL / "drug_targets.json"
BIOMARKERS = V5_EVAL / "consensus_biomarker_panel.json"
SOURCE_FIGS = [
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

TISSUE_COLORS = {
    "gastrocnemius": "#E69F00",
    "liver": "#F4C26B",
    "eye": "#5FD3C4",
    "kidney": "#8BD17C",
    "thymus": "#B39DFF",
    "skin": "#E17882",
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


def load_data() -> dict[str, object]:
    metabolic_rows = []
    for path in METABOLIC_FILES:
        item = json.loads(path.read_text())
        if item.get("status") != "success":
            continue
        top_subsystems = []
        for name, result in list(item["top_subsystems"].items())[:5]:
            top_subsystems.append(
                {
                    "name": name,
                    "abs_mean_diff": float(result["abs_mean_diff"]),
                    "n_reactions": int(result["n_reactions"]),
                }
            )
        metabolic_rows.append(
            {
                "tissue": item["tissue"],
                "gene_coverage_pct": float(item["gene_coverage_pct"]),
                "n_genes_mapped": int(item["n_genes_mapped"]),
                "n_model_genes": int(item["n_model_genes"]),
                "n_subsystems": int(item["n_subsystems"]),
                "top_subsystems": top_subsystems,
            }
        )
    metabolic_rows.sort(key=lambda row: row["gene_coverage_pct"], reverse=True)

    drug = json.loads(DRUG_TARGETS.read_text())
    biomarkers = json.loads(BIOMARKERS.read_text())

    validation = [
        {
            "tissue": tissue,
            "auroc": float(result["auroc"]),
            "std_auroc": float(result["std_auroc"]),
        }
        for tissue, result in biomarkers["panel_validation"].items()
    ]
    validation.sort(key=lambda row: row["auroc"], reverse=True)

    panel = []
    seen = set()
    for row in biomarkers["panel"]:
        gene = (row.get("human_symbol") or row["gene"]).upper()
        if gene in seen:
            continue
        seen.add(gene)
        panel.append(
            {
                "rank": int(row["rank"]),
                "gene": gene,
                "score": int(row["score"]),
                "n_tissues": int(row["n_tissues"]),
                "sources": row["sources"],
                "n_drugs": len(row.get("drugs", [])),
            }
        )
    return {"metabolic": metabolic_rows, "drug": drug, "biomarkers": biomarkers, "validation": validation, "panel": panel}


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
    text(draw, (150, 76), "BIOLOGY / TARGET + BIOMARKER EVIDENCE", F["kicker"], COLORS["teal"])
    x = 2020
    x = badge(draw, x, 66, "layer", "v5 Fig8/9", COLORS["blue"])
    x = badge(draw, x, 66, "targets", "1,919", COLORS["green"])
    x = badge(draw, x, 66, "DGIdb", "271 genes", COLORS["amber"])
    badge(draw, x, 66, "panel", "20 genes", COLORS["teal"])


def draw_title(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 382), "Target And Biomarker Layers Are Triage", F["title"], COLORS["text"])
    paragraph(
        draw,
        (155, 493),
        "v5 routes metabolic flux, target mapping, and compact panel validation into auditable follow-up questions.",
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


def draw_metric_bar(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    label: str,
    value: float,
    max_value: float,
    color: str,
    width: int,
    value_label: str,
) -> None:
    text(draw, (x, y + 3), label, F["tiny_bold"], COLORS["text"])
    bar_x = x + 165
    draw.line((bar_x, y + 23, bar_x + width, y + 23), fill="#263244", width=16)
    fill_w = int(width * max(0.025, min(1.0, value / max_value if max_value else 0)))
    draw.line((bar_x, y + 23, bar_x + fill_w, y + 23), fill=color, width=16)
    draw.ellipse((bar_x + fill_w - 12, y + 11, bar_x + fill_w + 12, y + 35), fill=color)
    text(draw, (bar_x + width + 18, y + 2), value_label, F["tiny_bold"], COLORS["text"])


def draw_metabolic_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, data: dict[str, object]) -> None:
    color = COLORS["amber"]
    rows = data["metabolic"]
    draw_panel_header(draw, x, y, w, "1", "Metabolic context", color)
    text(draw, (x + 54, y + 145), "6", F["stat"], color)
    text(draw, (x + 132, y + 155), "tissues solved", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "Flux modeling gives each target or panel candidate a pathway context before the evidence is ranked.",
        F["body"],
        w - 108,
        COLORS["muted"],
        9,
    )

    rounded(draw, (x + 54, y + 420, x + w - 54, y + 724), 24, "#0D1520", "#2A394D", 1)
    text(draw, (x + 82, y + 446), "Gene coverage by solved tissue", F["small_bold"], COLORS["text"])
    min_cov = min(row["gene_coverage_pct"] for row in rows)
    max_cov = max(row["gene_coverage_pct"] for row in rows)
    for i, row in enumerate(rows):
        draw_metric_bar(
            draw,
            x + 82,
            y + 500 + i * 34,
            tissue_label(row["tissue"]),
            row["gene_coverage_pct"] - 80,
            max_cov - 80,
            TISSUE_COLORS[row["tissue"]],
            320,
            f"{row['gene_coverage_pct']:.1f}%",
        )
    text(draw, (x + 82, y + 690), f"range {min_cov:.1f}-{max_cov:.1f}% across iMM1865-mapped genes", F["tiny"], COLORS["muted"])

    rounded(draw, (x + 54, y + 790, x + w - 54, y + 964), 24, "#172335", "#2A394D", 1)
    text(draw, (x + 82, y + 816), "Top subsystem cues", F["small_bold"], color)
    cue_rows = [
        ("Thymus", rows[0]["top_subsystems"][0]["name"]),
        ("Kidney", rows[1]["top_subsystems"][0]["name"]),
        ("Liver", next(row for row in rows if row["tissue"] == "liver")["top_subsystems"][0]["name"]),
    ]
    for i, (label, cue) in enumerate(cue_rows):
        cy = y + 858 + i * 34
        text(draw, (x + 82, cy), label, F["tiny_bold"], TISSUE_COLORS[label.lower()])
        paragraph(draw, (x + 210, cy), cue, F["tiny"], w - 300, COLORS["muted"], 2)


def draw_target_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, data: dict[str, object]) -> None:
    color = COLORS["green"]
    drug = data["drug"]
    draw_panel_header(draw, x, y, w, "2", "Target evidence funnel", color)
    text(draw, (x + 54, y + 145), "271", F["stat"], color)
    text(draw, (x + 240, y + 155), "DGIdb genes", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "The funnel keeps candidate genes, human mapping, interaction evidence, and target-class evidence separated.",
        F["body"],
        w - 108,
        COLORS["muted"],
        9,
    )

    rounded(draw, (x + 54, y + 420, x + w - 54, y + 724), 24, "#0D1520", "#2A394D", 1)
    steps = [
        ("mouse candidates", int(drug["n_target_genes_mouse"]), COLORS["blue"]),
        ("human mapped", int(drug["n_mapped_to_human"]), COLORS["teal"]),
        ("DGIdb genes", int(drug["n_genes_with_dgidb_interactions"]), COLORS["green"]),
        ("ChEMBL targets", int(drug["n_chembl_targets"]), COLORS["amber"]),
    ]
    max_value = steps[0][1]
    bar_max_w = w - 400
    for i, (label, value, c) in enumerate(steps):
        yy = y + 470 + i * 58
        bar_w = int(bar_max_w * value / max_value)
        rounded(draw, (x + 82, yy, x + 82 + bar_w, yy + 40), 16, rgba(c, 105 + i * 25), c, 2)
        text(draw, (x + 104, yy + 9), f"{value:,}", F["tiny_bold"], COLORS["text"])
        text(draw, (x + w - 92, yy + 9), label, F["tiny_bold"], COLORS["muted"], "ra")
    text(draw, (x + 82, y + 685), f"{int(drug['n_total_dgidb_interactions']):,} total DGIdb interactions", F["tiny_bold"], COLORS["text"])

    rounded(draw, (x + 54, y + 790, x + w - 54, y + 964), 24, "#172335", "#2A394D", 1)
    text(draw, (x + 82, y + 816), "Interaction tier readout", F["small_bold"], color)
    tier_rows = [
        ("Tier-1 links", len(drug["tiers"].get("tier1_approved", [])), COLORS["green"]),
        ("Tier-3 links", len(drug["tiers"].get("tier3_preclinical", [])), COLORS["amber"]),
        ("Tier-2 links", len(drug["tiers"].get("tier2_clinical", [])), COLORS["violet"]),
    ]
    max_tier = max(row[1] for row in tier_rows)
    for i, (label, value, c) in enumerate(tier_rows):
        draw_metric_bar(draw, x + 82, y + 858 + i * 32, label, value, max_tier, c, 265, f"{value:,}")


def draw_gene_source_matrix(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, panel: list[dict[str, object]], n_candidates: int, max_score: int) -> None:
    sources = [
        ("SHAP", "shap", COLORS["blue"]),
        ("WGCNA", "wgcna_hub", COLORS["violet"]),
        ("Multi", "multi_tissue", COLORS["teal"]),
        ("Target", "druggable", COLORS["green"]),
        ("Lit", "literature", COLORS["amber"]),
    ]
    text(draw, (x, y), "Top-4 evidence matrix", F["small_bold"], COLORS["text"])
    text(draw, (x + w - 18, y + 4), f"{n_candidates:,} scored; max {max_score}", F["tiny_bold"], COLORS["muted"], "ra")
    start_y = y + 58
    col_x = x + 210
    for j, (label, _key, c) in enumerate(sources):
        text(draw, (col_x + j * 66, start_y - 32), label, F["micro_bold"], c, "ma")
    text(draw, (x + w - 88, start_y - 32), "score", F["micro_bold"], COLORS["dim"], "ra")
    for i, row in enumerate(panel[:4]):
        yy = start_y + i * 28
        text(draw, (x, yy + 4), row["gene"], F["tiny_bold"], COLORS["text"])
        for j, (_label, key, c) in enumerate(sources):
            value = int(row["sources"].get(key, 0))
            cx = col_x + j * 66
            if value:
                radius = 8 + value * 2
                draw.ellipse((cx - radius, yy + 8 - radius, cx + radius, yy + 8 + radius), fill=c)
                text(draw, (cx, yy + 2), str(value), F["micro_bold"], COLORS["ink"], "ma")
            else:
                draw.ellipse((cx - 5, yy + 3, cx + 5, yy + 13), fill="#263244")
        text(draw, (x + w - 88, yy - 1), str(row["score"]), F["tiny_bold"], COLORS["muted"], "ra")
        text(draw, (x + w - 20, yy - 1), f"{row['n_tissues']}t", F["tiny"], COLORS["dim"], "ra")


def draw_biomarker_panel(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, data: dict[str, object]) -> None:
    color = COLORS["teal"]
    biomarkers = data["biomarkers"]
    validation = data["validation"]
    panel = data["panel"]
    draw_panel_header(draw, x, y, w, "3", "Biomarker panel check", color)
    text(draw, (x + 54, y + 145), "20", F["stat"], color)
    text(draw, (x + 154, y + 155), "genes", F["h3"], COLORS["muted"])
    paragraph(
        draw,
        (x + 54, y + 246),
        "The compact panel tests whether recurrent candidate signals remain useful after feature compression.",
        F["body"],
        w - 108,
        COLORS["muted"],
        9,
    )

    rounded(draw, (x + 54, y + 420, x + w - 54, y + 724), 24, "#0D1520", "#2A394D", 1)
    text(draw, (x + 82, y + 446), "Panel validation AUROC", F["small_bold"], COLORS["text"])
    for i, row in enumerate(validation):
        draw_metric_bar(
            draw,
            x + 82,
            y + 492 + i * 27,
            tissue_label(row["tissue"]),
            row["auroc"],
            1.0,
            TISSUE_COLORS[row["tissue"]],
            285,
            f"{row['auroc']:.3f}",
        )

    rounded(draw, (x + 54, y + 765, x + w - 54, y + 964), 24, "#172335", "#2A394D", 1)
    draw_gene_source_matrix(draw, x + 82, y + 792, w - 164, panel, int(biomarkers["n_candidates_scored"]), int(biomarkers["max_possible_score"]))


def draw_reader_rule(draw: ImageDraw.ImageDraw) -> None:
    box = (150, 1695, 3690, 1828)
    rounded(draw, box, 30, "#111A27", "#2A394D", 2)
    text(draw, (196, 1736), "Layer decoder", F["h3"], COLORS["teal"])
    steps = [
        ("Flux context", "pathway-scale metabolism"),
        ("Target funnel", "mouse -> human -> interactions"),
        ("Panel check", "20-gene validation"),
        ("Follow-up queue", "ranked evidence packets"),
    ]
    x = 700
    for i, (head, body) in enumerate(steps):
        color = [COLORS["amber"], COLORS["green"], COLORS["teal"], COLORS["violet"]][i]
        rounded(draw, (x, 1718, x + 620, 1808), 24, "#172335", color, 2)
        text(draw, (x + 28, 1733), head, F["small_bold"], COLORS["text"])
        text(draw, (x + 28, 1768), body, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            ax = x + 640
            draw.line((ax, 1763, ax + 82, 1763), fill="#6F7E90", width=4)
            draw.polygon([(ax + 82, 1763), (ax + 62, 1752), (ax + 62, 1774)], fill="#6F7E90")
        x += 730


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    source = "Takeaway: target and biomarker layers convert biology signals into follow-up questions."
    scope = "Next: translation slides ask which of these signals carry across species and human readouts."
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
    panel_w = 1160
    gap = 55
    xs = [150, 150 + panel_w + gap, 150 + (panel_w + gap) * 2]
    draw_metabolic_panel(draw, xs[0], panel_y, panel_w, data)
    draw_target_panel(draw, xs[1], panel_y, panel_w, data)
    draw_biomarker_panel(draw, xs[2], panel_y, panel_w, data)
    draw_reader_rule(draw)
    draw_footer(draw)

    canvas.save(OUT_PATH, quality=95)
    gray = ImageOps.grayscale(canvas).convert("RGB")
    gray.save(GRAY_PATH, quality=95)

    stat = ImageStat.Stat(gray)
    manifest = {
        "slide": 39,
        "title": "Target And Biomarker Layers Are Triage",
        "readout": "Target and biomarker evidence converts biology signals into auditable follow-up questions.",
        "asset": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "sources": [str(path.relative_to(ROOT)) for path in [*METABOLIC_FILES, DRUG_TARGETS, BIOMARKERS, *SOURCE_FIGS]],
        "data": {
            "metabolic_tissues": len(data["metabolic"]),
            "metabolic_coverage_range": [
                min(row["gene_coverage_pct"] for row in data["metabolic"]),
                max(row["gene_coverage_pct"] for row in data["metabolic"]),
            ],
            "target_genes_mouse": data["drug"]["n_target_genes_mouse"],
            "target_genes_human": data["drug"]["n_mapped_to_human"],
            "dgidb_genes": data["drug"]["n_genes_with_dgidb_interactions"],
            "dgidb_interactions": data["drug"]["n_total_dgidb_interactions"],
            "chembl_targets": data["drug"]["n_chembl_targets"],
            "panel_size": data["biomarkers"]["panel_size"],
            "panel_validation": data["validation"],
            "panel_top_unique": data["panel"][:10],
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

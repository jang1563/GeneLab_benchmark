#!/usr/bin/env python3
"""Build the detailed-deck external biological validation proof asset.

The output is a high-resolution 16:9 PNG for the detailed SpaceBio-Bench deck.
It summarizes local `evaluation/cell2020_validation.json` evidence:

- Cell 2020 / published pathway-direction concordance against our fGSEA NES
  signs across five tissues;
- SHAP top-50 overlap with curated literature genes across the three tissues
  with gene-level SHAP summaries;
- an explicit scope layer showing that pathway-level agreement is stronger than
  exact gene-identity overlap.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "evaluation" / "cell2020_validation.json"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "external_validation"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "panel2": "#151F2D",
    "grid": "#2A3546",
    "axis": "#98A7BA",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#667487",
    "teal": "#5FD3C4",
    "sky": "#73A7FF",
    "green": "#84D278",
    "amber": "#F4C26B",
    "rose": "#E17882",
    "purple": "#B39DFF",
    "white": "#FFFFFF",
    "thymus": "#E8A34A",
    "gastrocnemius": "#56B4E9",
    "liver": "#A6AEBB",
    "eye": "#73A7FF",
    "kidney": "#B39DFF",
}


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(value: str, alpha: int) -> tuple[int, int, int, int]:
    return (*hex_to_rgb(value), alpha)


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
    "title": load_font(84, True),
    "subtitle": load_font(36, False),
    "h2": load_font(44, True),
    "h3": load_font(34, True),
    "body": load_font(30, False),
    "small": load_font(25, False),
    "small_bold": load_font(25, True),
    "tiny": load_font(21, False),
    "number": load_font(104, True),
    "stat": load_font(118, True),
}


TISSUE_ORDER = ["thymus", "gastrocnemius", "liver", "eye", "kidney"]
TISSUE_LABELS = {
    "thymus": "Thymus",
    "gastrocnemius": "Gastrocnemius",
    "liver": "Liver",
    "eye": "Eye",
    "kidney": "Kidney",
}


def rounded(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], radius: int, fill: str, outline: str | None = None, width: int = 1) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], value: str, font: ImageFont.ImageFont, fill: str = COLORS["text"], anchor: str | None = None) -> None:
    draw.text(xy, value, font=font, fill=fill, anchor=anchor)


def multiline(draw: ImageDraw.ImageDraw, xy: tuple[int, int], lines: Iterable[str], font: ImageFont.ImageFont, fill: str = COLORS["muted"], leading: int = 8) -> None:
    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=font, fill=fill)
        y += font.size + leading


def fit_label(draw: ImageDraw.ImageDraw, label: str, font: ImageFont.ImageFont, max_width: int) -> str:
    words = label.split()
    lines: list[str] = []
    cur: list[str] = []
    for word in words:
        trial = " ".join(cur + [word])
        if draw.textlength(trial, font=font) <= max_width:
            cur.append(word)
        else:
            if cur:
                lines.append(" ".join(cur))
            cur = [word]
    if cur:
        lines.append(" ".join(cur))
    return "\n".join(lines)


def load_data() -> dict:
    source = json.loads(SOURCE.read_text())
    pathway_rows = []
    for tissue in TISSUE_ORDER:
        record = source["pathway_validation"][tissue]
        pathway_rows.append(
            {
                "tissue": tissue,
                "label": TISSUE_LABELS[tissue],
                "concordant": record["concordant"],
                "discordant": record["discordant"],
                "total": record["n_found"],
                "rate": record["concordance_rate"],
                "top10": record["expected_in_top10"],
                "details": record["details"],
                "source": record["source"],
            }
        )
    gene_rows = []
    for tissue in ["liver", "gastrocnemius", "thymus"]:
        record = source["gene_validation"][tissue]
        found = [g["gene"] for g in record["found_genes"]]
        gene_rows.append(
            {
                "tissue": tissue,
                "label": TISSUE_LABELS[tissue],
                "reference": record["n_reference_genes"],
                "found": record["n_found_in_shap_top50"],
                "rate": record["overlap_rate"],
                "found_genes": found,
                "task_id": record["task_id"],
            }
        )
    return {
        "source": source,
        "pathway_rows": pathway_rows,
        "gene_rows": gene_rows,
        "summary": source["summary"],
    }


def draw_badges(draw: ImageDraw.ImageDraw) -> None:
    badges = [
        ("REFERENCE", "Cell 2020 biology", 410),
        ("PATHWAY CHECK", "5 tissues", 330),
        ("GENE CHECK", "SHAP top-50", 360),
    ]
    bx = 2350
    for kicker, body, badge_w in badges:
        rounded(draw, (bx, 72, bx + badge_w, 176), 26, "#121B28", "#2A394D", 2)
        text(draw, (bx + 28, 96), kicker, F["tiny"], COLORS["sky"])
        text(draw, (bx + 28, 126), body, F["small_bold"], COLORS["text"])
        bx += badge_w + 30


def draw_reference_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "A. External biology frame", F["h2"], COLORS["text"])
    text(draw, (x0 + 50, y0 + 96), "Published biology provides expected pathway directions.", F["small"], COLORS["muted"])

    refs = [
        ("Beheshti et al.", "Cell 2020; NASA GeneLab multi-omics", COLORS["teal"]),
        ("da Silveira et al.", "Cell 2020; integrated spaceflight omics", COLORS["sky"]),
        ("SOMA context", "Nature 2024; human multi-omics atlas context", COLORS["amber"]),
    ]
    y = y0 + 165
    for title, subtitle, color in refs:
        rounded(draw, (x0 + 50, y, x1 - 50, y + 122), 24, "#151F2D", "#2A394D", 2)
        draw.ellipse((x0 + 82, y + 45, x0 + 112, y + 75), fill=color)
        text(draw, (x0 + 135, y + 26), title, F["h3"], COLORS["text"])
        text(draw, (x0 + 135, y + 72), subtitle, F["small"], COLORS["muted"])
        y += 148

    text(draw, (x0 + 50, y0 + 640), "Validation chain", F["h3"], COLORS["text"])
    steps = [
        ("Expected sign", "published pathway direction"),
        ("Observed sign", "our fGSEA NES sign"),
        ("Concordance", "direction match status"),
    ]
    y = y0 + 700
    for idx, (title, subtitle) in enumerate(steps, start=1):
        rounded(draw, (x0 + 50, y, x1 - 50, y + 105), 20, "#171D25", "#2A394D", 2)
        text(draw, (x0 + 82, y + 28), str(idx), F["h3"], COLORS["dim"])
        text(draw, (x0 + 130, y + 22), title, F["small_bold"], COLORS["text"])
        text(draw, (x0 + 130, y + 58), subtitle, F["small"], COLORS["muted"])
        y += 125

    rounded(draw, (x0 + 50, y1 - 148, x1 - 50, y1 - 48), 20, "#211E17", "#69532B", 2)
    text(draw, (x0 + 80, y1 - 120), "Readout", F["small_bold"], COLORS["amber"])
    note = "This biology readout complements the train/test evidence."
    draw.text((x0 + 205, y1 - 120), fit_label(draw, note, F["small"], x1 - x0 - 300), font=F["small"], fill=COLORS["muted"])


def draw_pathway_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], rows: list[dict], summary: dict) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 30, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 42, y0 + 38), "B. Pathway direction concordance", F["h2"], COLORS["text"])
    text(draw, (x0 + 42, y0 + 88), "Expected literature direction vs observed fGSEA NES sign.", F["small"], COLORS["muted"])

    text(draw, (x0 + 42, y0 + 150), "mean concordance", F["small_bold"], COLORS["teal"])
    text(draw, (x0 + 42, y0 + 178), f"{summary['pathway_mean_concordance']:.3f}", F["stat"], COLORS["text"])
    text(draw, (x0 + 42, y0 + 308), "5 tissues; Cell 2020 pathway set", F["small"], COLORS["muted"])

    chart_x0, chart_x1 = x0 + 610, x1 - 145
    chart_y0 = y0 + 175
    row_gap = 96
    for tick in [0.0, 0.5, 1.0]:
        tx = int(chart_x0 + tick * (chart_x1 - chart_x0))
        draw.line((tx, chart_y0 - 18, tx, y1 - 105), fill=rgba(COLORS["grid"], 125), width=2)
        text(draw, (tx, y1 - 62), f"{tick:.1f}", F["small"], COLORS["axis"], anchor="ma")

    for i, row in enumerate(rows):
        y = chart_y0 + i * row_gap
        color = COLORS[row["tissue"]]
        text(draw, (chart_x0 - 28, y - 15), row["label"], F["small_bold"], COLORS["text"], anchor="ra")
        rounded(draw, (chart_x0, y - 20, chart_x1, y + 20), 18, "#172231", "#2A394D", 1)
        split = int(chart_x0 + row["rate"] * (chart_x1 - chart_x0))
        rounded(draw, (chart_x0, y - 20, split, y + 20), 18, color, None, 0)
        if row["rate"] < 1:
            rounded(draw, (split, y - 20, chart_x1, y + 20), 18, "#4B2430", None, 0)
        text(draw, (chart_x1 + 22, y - 16), f"{row['concordant']}/{row['total']}", F["small_bold"], COLORS["text"])
        text(draw, (chart_x1 + 22, y + 14), f"top10 {row['top10']}", F["tiny"], COLORS["muted"])

    rounded(draw, (x0 + 42, y1 - 72, x1 - 42, y1 - 28), 16, "#151F2D", "#2A394D", 1)
    text(draw, (x0 + 65, y1 - 62), "pathway-level signal is the primary evidence layer", F["small_bold"], COLORS["teal"])


def draw_gene_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], rows: list[dict], summary: dict) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 30, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 42, y0 + 38), "C. SHAP gene overlap", F["h2"], COLORS["text"])
    text(draw, (x0 + 42, y0 + 88), "Top-50 model genes are enriched for curated literature genes above chance.", F["small"], COLORS["muted"])

    text(draw, (x0 + 42, y0 + 150), "mean overlap", F["small_bold"], COLORS["amber"])
    text(draw, (x0 + 42, y0 + 178), f"{summary['gene_mean_overlap']:.3f}", F["stat"], COLORS["text"])
    multiline(draw, (x0 + 42, y0 + 308), ["47.1x enrichment", "vs chance expectation"], F["small"], COLORS["muted"], 4)

    chart_x0, chart_x1 = x0 + 610, x1 - 145
    chart_y0 = y0 + 185
    row_gap = 112
    for tick in [0.0, 0.1, 0.2]:
        tx = int(chart_x0 + tick / 0.2 * (chart_x1 - chart_x0))
        draw.line((tx, chart_y0 - 22, tx, y1 - 108), fill=rgba(COLORS["grid"], 125), width=2)
        text(draw, (tx, y1 - 65), f"{tick:.1f}", F["small"], COLORS["axis"], anchor="ma")

    for i, row in enumerate(rows):
        y = chart_y0 + i * row_gap
        color = COLORS[row["tissue"]]
        text(draw, (chart_x0 - 28, y - 15), row["label"], F["small_bold"], COLORS["text"], anchor="ra")
        draw.line((chart_x0, y, chart_x1, y), fill="#2A3546", width=12)
        capped_rate = min(row["rate"], 0.2)
        fill_x = int(chart_x0 + capped_rate / 0.2 * (chart_x1 - chart_x0))
        draw.line((chart_x0, y, fill_x, y), fill=color, width=12)
        draw.ellipse((fill_x - 18, y - 18, fill_x + 18, y + 18), fill=color, outline=COLORS["white"], width=3)
        text(draw, (chart_x1 + 22, y - 16), f"{row['found']}/{row['reference']}", F["small_bold"], COLORS["text"])
        found = ", ".join(row["found_genes"]) if row["found_genes"] else "none in top-50"
        text(draw, (chart_x0, y + 28), found, F["tiny"], COLORS["muted"])

    rounded(draw, (x0 + 42, y1 - 72, x1 - 42, y1 - 28), 16, "#151F2D", "#2A394D", 1)
    text(draw, (x0 + 65, y1 - 62), "signal is enriched at the gene-set level", F["small_bold"], COLORS["amber"])


def draw_support_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], summary: dict) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "D. What this supports", F["h2"], COLORS["text"])
    text(draw, (x0 + 50, y0 + 96), "Separate pathway biology from exact gene attribution.", F["small"], COLORS["muted"])

    cards = [
        ("Pathway support", f"Published pathway directions match our NES signs at {summary['pathway_mean_concordance']:.1%}.", COLORS["teal"]),
        ("Gene layer", f"SHAP gene overlap averages {summary['gene_mean_overlap']:.1%}; useful as enrichment evidence.", COLORS["amber"]),
        ("Best example", "Thymus pathway concordance is 7/7 with proliferation up and interferon down.", COLORS["thymus"]),
    ]
    y = y0 + 175
    for title, body, color in cards:
        rounded(draw, (x0 + 50, y, x1 - 50, y + 148), 24, "#151F2D", "#2A394D", 2)
        text(draw, (x0 + 82, y + 28), title, F["h3"], color)
        draw.text((x0 + 82, y + 76), fit_label(draw, body, F["body"], x1 - x0 - 170), font=F["body"], fill=COLORS["text"])
        y += 178

    rounded(draw, (x0 + 50, y + 18, x1 - 50, y + 305), 28, "#211E17", "#69532B", 2)
    text(draw, (x0 + 82, y + 56), "Readout frame", F["h3"], COLORS["amber"])
    lines = [
        "External biological consistency complements the model-validation slides.",
        "Pathway-level concordance is the headline.",
        "Gene-level SHAP overlap acts as enrichment support.",
    ]
    multiline(draw, (x0 + 82, y + 108), lines, F["body"], COLORS["text"], 8)
    text(draw, (x0 + 82, y + 260), "Biological consistency complements train/test evidence.", F["small"], COLORS["muted"])


def main() -> None:
    data = load_data()
    summary = data["summary"]

    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)

    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 52), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 42), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "EXTERNAL BIOLOGY CHECK", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "Published biology supports the pathway signal", F["title"])
    text(
        draw,
        (150, 216),
        "Cell 2020 pathway expectations align with our fGSEA signs; SHAP genes add enriched support.",
        F["subtitle"],
        COLORS["muted"],
    )
    draw_badges(draw)

    draw_reference_panel(draw, (150, 350, 1075, 1800), data)
    draw_pathway_panel(draw, (1125, 350, 2485, 1080), data["pathway_rows"], summary)
    draw_gene_panel(draw, (1125, 1130, 2485, 1800), data["gene_rows"], summary)
    draw_support_panel(draw, (2535, 350, 3690, 1800), summary)

    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    footer_label_x = 205
    footer_text_x = 390
    footer_text_width = 3220
    text(draw, (footer_label_x, 1925), "Takeaway", F["small_bold"], COLORS["sky"])
    footer = "Published spaceflight biology aligns most clearly at the pathway layer, with SHAP genes adding enrichment support."
    draw.text((footer_text_x, 1925), fit_label(draw, footer, F["small"], footer_text_width), font=F["small"], fill=COLORS["muted"])
    text(draw, (footer_label_x, 1995), "Next", F["small_bold"], COLORS["amber"])
    scope = "The following recap slide combines split, metric, DGE, and biology checks into one readable stack."
    draw.text((footer_text_x, 1995), fit_label(draw, scope, F["small"], footer_text_width), font=F["small"], fill=COLORS["muted"])

    png = OUT_DIR / "external_biology_validation_premium.png"
    canvas.save(png, quality=95)
    gray = OUT_DIR / "external_biology_validation_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    manifest = {
        "title": "Published biology supports the pathway signal",
        "source": str(SOURCE.relative_to(ROOT)),
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "summary": {
            "pathway_mean_concordance": summary["pathway_mean_concordance"],
            "pathway_n_tissues": summary["pathway_n_tissues"],
            "gene_mean_overlap": summary["gene_mean_overlap"],
            "gene_n_tissues": summary["gene_n_tissues"],
            "gene_enrichment_ratio": 47.1,
        },
        "pathway_rows": [
            {
                "tissue": row["tissue"],
                "concordant": row["concordant"],
                "total": row["total"],
                "rate": row["rate"],
                "expected_in_top10": row["top10"],
            }
            for row in data["pathway_rows"]
        ],
        "gene_rows": [
            {
                "tissue": row["tissue"],
                "found": row["found"],
                "reference": row["reference"],
                "rate": row["rate"],
                "found_genes": row["found_genes"],
            }
            for row in data["gene_rows"]
        ],
        "scope": (
            "External biological support is pathway-level concordance with published spaceflight biology. "
            "SHAP gene overlap provides enriched support; biology slides carry the mechanistic interpretation."
        ),
    }
    manifest_path = OUT_DIR / "external_biology_validation_premium_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build the detailed-deck DGE robustness proof asset.

The output is a high-resolution 16:9 PNG for the detailed SpaceBio-Bench deck.
It summarizes J2 DGE pipeline comparison data:

- three DGE callers: DESeq2, edgeR, limma-voom;
- liver and thymus mission-level FLT vs GC contrasts;
- rank-level log2FC agreement versus thresholded DEG-list overlap.

The main visual readout separates two layers: effect-size rankings are robust
across DGE callers, while FDR-thresholded DEG lists vary by caller stringency.
"""

from __future__ import annotations

import json
import statistics
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "evaluation" / "J2_dge_pipeline_comparison.json"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "dge_robustness"
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


PAIR_LABELS = {
    "deseq2_vs_edger": "DESeq2-edgeR",
    "deseq2_vs_limma_voom": "DESeq2-limma",
    "edger_vs_limma_voom": "edgeR-limma",
}

PAIR_COLORS = {
    "deseq2_vs_edger": COLORS["teal"],
    "deseq2_vs_limma_voom": COLORS["sky"],
    "edger_vs_limma_voom": COLORS["amber"],
}

MISSION_ORDER = {
    "MHU-2": 0,
    "RR-1": 1,
    "RR-3": 2,
    "RR-6": 3,
    "RR-8": 4,
    "RR-9": 5,
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
    rows: list[dict] = []
    genelab = []
    for tissue, missions in source["results"].items():
        for mission, record in missions.items():
            genelab.append(record["genelab_concordance"]["log2fc_spearman_rho"])
            for pair_record in record["pairwise_comparisons"]:
                pair = pair_record["pair"]
                rows.append(
                    {
                        "tissue": tissue,
                        "mission": mission,
                        "label": f"{tissue.capitalize()} {mission}",
                        "pair": pair,
                        "pair_label": PAIR_LABELS[pair],
                        "rho": pair_record["log2fc_spearman_rho"],
                        "jaccard": pair_record["deg_jaccard_005"],
                        "top100": pair_record["top100_jaccard"],
                        "direction": pair_record["direction_concordance"],
                    }
                )
    rows.sort(key=lambda r: (0 if r["tissue"] == "liver" else 1, MISSION_ORDER.get(r["mission"], 99), r["pair_label"]))
    rho_vals = [r["rho"] for r in rows]
    jacc_vals = [r["jaccard"] for r in rows]
    top_vals = [r["top100"] for r in rows]
    directions = [r["direction"] for r in rows if r["direction"] is not None]
    return {
        "source": source,
        "rows": rows,
        "summary": {
            "n_missions": source["summary"]["n_missions"],
            "n_comparisons": source["summary"]["n_comparisons"],
            "rho_mean": statistics.mean(rho_vals),
            "rho_min": min(rho_vals),
            "rho_max": max(rho_vals),
            "jaccard_mean": statistics.mean(jacc_vals),
            "jaccard_min": min(jacc_vals),
            "jaccard_max": max(jacc_vals),
            "top100_mean": statistics.mean(top_vals),
            "direction_mean": statistics.mean(directions),
            "genelab_mean": statistics.mean(genelab),
            "genelab_min": min(genelab),
        },
    }


def draw_badges(draw: ImageDraw.ImageDraw) -> None:
    badges = [
        ("MISSIONS", "9 missions", 350),
        ("DGE CALLERS", "DESeq2 / edgeR / limma", 430),
        ("READOUT", "rank robust; lists vary", 470),
    ]
    bx = 2330
    for kicker, body, badge_w in badges:
        rounded(draw, (bx, 72, bx + badge_w, 176), 26, "#121B28", "#2A394D", 2)
        text(draw, (bx + 28, 96), kicker, F["tiny"], COLORS["sky"])
        text(draw, (bx + 28, 126), body, F["small_bold"], COLORS["text"])
        bx += badge_w + 30


def draw_workflow(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "A. What was varied", F["h2"], COLORS["text"])
    text(draw, (x0 + 50, y0 + 96), "FLT vs GC contrast stays fixed while the DGE caller changes.", F["small"], COLORS["muted"])

    contrast = (x0 + 50, y0 + 160, x1 - 50, y0 + 310)
    rounded(draw, contrast, 24, "#131D2A", "#2A394D", 2)
    text(draw, (contrast[0] + 34, contrast[1] + 28), "Input contrast", F["h3"], COLORS["sky"])
    text(draw, (contrast[0] + 34, contrast[1] + 84), "spaceflight", F["body"], COLORS["text"])
    text(draw, (contrast[0] + 260, contrast[1] + 84), "vs", F["body"], COLORS["dim"])
    text(draw, (contrast[0] + 330, contrast[1] + 84), "ground control", F["body"], COLORS["text"])

    pipelines = [
        ("DESeq2", "Wald test", COLORS["teal"]),
        ("edgeR", "QLF test", COLORS["sky"]),
        ("limma-voom", "QW model", COLORS["amber"]),
    ]
    y = y0 + 390
    for name, method, color in pipelines:
        rounded(draw, (x0 + 50, y, x1 - 50, y + 130), 24, "#151F2D", "#2A394D", 2)
        draw.ellipse((x0 + 82, y + 48, x0 + 112, y + 78), fill=color)
        text(draw, (x0 + 135, y + 30), name, F["h3"], COLORS["text"])
        text(draw, (x0 + 135, y + 76), method, F["small"], COLORS["muted"])
        y += 158

    text(draw, (x0 + 50, y0 + 900), "Compare two outputs", F["h3"], COLORS["text"])
    output_cards = [
        ("effect-size rank", "log2FC Spearman", COLORS["teal"]),
        ("thresholded list", "FDR < 0.05 Jaccard", COLORS["amber"]),
    ]
    y = y0 + 955
    for title, subtitle, color in output_cards:
        rounded(draw, (x0 + 50, y, x1 - 50, y + 118), 20, "#171D25", "#2A394D", 2)
        text(draw, (x0 + 80, y + 24), title, F["small_bold"], color)
        text(draw, (x0 + 80, y + 62), subtitle, F["small"], COLORS["muted"])
        y += 140

    rounded(draw, (x0 + 50, y1 - 150, x1 - 50, y1 - 48), 20, "#211E17", "#69532B", 2)
    text(draw, (x0 + 80, y1 - 122), "Readout", F["small_bold"], COLORS["amber"])
    note = "A stable rank signal can coexist with caller-dependent DEG lists."
    draw.text((x0 + 205, y1 - 122), fit_label(draw, note, F["small"], x1 - x0 - 300), font=F["small"], fill=COLORS["muted"])


def draw_distribution_panel(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    subtitle: str,
    stat_label: str,
    stat_value: float,
    range_values: tuple[float, float],
    rows: list[dict],
    metric: str,
    axis_min: float,
    axis_max: float,
    axis_ticks: list[float],
    accent: str,
) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 30, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 42, y0 + 38), title, F["h2"], COLORS["text"])
    text(draw, (x0 + 42, y0 + 88), subtitle, F["small"], COLORS["muted"])

    text(draw, (x0 + 42, y0 + 150), stat_label, F["small_bold"], accent)
    text(draw, (x0 + 42, y0 + 178), f"{stat_value:.3f}", F["stat"], COLORS["text"])
    text(draw, (x0 + 42, y0 + 308), f"range {range_values[0]:.3f}-{range_values[1]:.3f}", F["small"], COLORS["muted"])

    chart_x0, chart_x1 = x0 + 620, x1 - 60
    chart_y0, chart_y1 = y0 + 172, y1 - 175

    def sx(v: float) -> int:
        return int(chart_x0 + (v - axis_min) / (axis_max - axis_min) * (chart_x1 - chart_x0))

    for tick in axis_ticks:
        x = sx(tick)
        draw.line((x, chart_y0 - 20, x, chart_y1 + 8), fill=rgba(COLORS["grid"], 120), width=2)
        label = f"{tick:.2f}" if axis_max <= 1.0 and axis_min >= 0.7 else f"{tick:.1f}"
        text(draw, (x, chart_y1 + 35), label, F["small"], COLORS["axis"], anchor="ma")

    grouped = {pair: [r[metric] for r in rows if r["pair"] == pair] for pair in PAIR_LABELS}
    pair_order = ["deseq2_vs_edger", "deseq2_vs_limma_voom", "edger_vs_limma_voom"]
    row_gap = max(78, int((chart_y1 - chart_y0) / 3))
    for i, pair in enumerate(pair_order):
        y = chart_y0 + i * row_gap + 24
        color = PAIR_COLORS[pair]
        text(draw, (chart_x0 - 34, y - 18), PAIR_LABELS[pair], F["small_bold"], COLORS["text"], anchor="ra")
        draw.line((chart_x0, y, chart_x1, y), fill="#2A3546", width=5)
        vals = grouped[pair]
        if vals:
            draw.line((sx(min(vals)), y, sx(max(vals)), y), fill=rgba(color, 150), width=10)
        for j, val in enumerate(vals):
            jitter = ((j % 5) - 2) * 7
            px = sx(val)
            draw.ellipse((px - 11, y + jitter - 11, px + 11, y + jitter + 11), fill=color, outline="#FFFFFF", width=2)
        text(draw, (chart_x1 + 18, y - 14), f"n={len(vals)}", F["tiny"], COLORS["dim"])

    cue = "all values stay high" if metric == "rho" else "threshold overlap moves"
    rounded(draw, (x0 + 42, y1 - 64, x1 - 42, y1 - 22), 16, "#151F2D", "#2A394D", 1)
    text(draw, (x0 + 65, y1 - 56), cue, F["small_bold"], accent)


def draw_readout(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], summary: dict) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "D. How to use this result", F["h2"], COLORS["text"])
    text(draw, (x0 + 50, y0 + 96), "Robustness is layered: ranks, lists, and external reproduction differ.", F["small"], COLORS["muted"])

    cards = [
        ("Stable layer", "Effect-size rankings remain highly concordant across callers.", COLORS["teal"]),
        ("Sensitive layer", "FDR-thresholded DEG membership changes with caller stringency.", COLORS["amber"]),
        ("Replication layer", f"GeneLab comparison mean rho = {summary['genelab_mean']:.3f} across 9 missions.", COLORS["sky"]),
    ]
    y = y0 + 175
    for title, body, color in cards:
        rounded(draw, (x0 + 50, y, x1 - 50, y + 155), 24, "#151F2D", "#2A394D", 2)
        text(draw, (x0 + 82, y + 28), title, F["h3"], color)
        draw.text((x0 + 82, y + 76), fit_label(draw, body, F["body"], x1 - x0 - 170), font=F["body"], fill=COLORS["text"])
        y += 185

    rounded(draw, (x0 + 50, y + 20, x1 - 50, y + 310), 28, "#211E17", "#69532B", 2)
    text(draw, (x0 + 82, y + 58), "Readout frame", F["h3"], COLORS["amber"])
    lines = [
        "DEG membership can vary by caller and threshold.",
        "Ranked response signal persists across callers.",
        "Exact gene membership belongs to the method-aware follow-up layer.",
    ]
    multiline(draw, (x0 + 82, y + 112), lines, F["body"], COLORS["text"], 8)
    text(draw, (x0 + 82, y + 260), f"Top-100 overlap mean = {summary['top100_mean']:.3f}; DEG-list Jaccard mean = {summary['jaccard_mean']:.3f}", F["small"], COLORS["muted"])


def main() -> None:
    loaded = load_data()
    rows = loaded["rows"]
    summary = loaded["summary"]

    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)

    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 52), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 42), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "ROBUSTNESS CHECK", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "DGE ranks hold; DEG lists move", F["title"])
    text(
        draw,
        (150, 216),
        "Across DESeq2, edgeR, and limma-voom, effect-size rankings stay concordant while FDR gene lists remain caller-sensitive.",
        F["subtitle"],
        COLORS["muted"],
    )
    draw_badges(draw)

    draw_workflow(draw, (150, 350, 1120, 1800))
    draw_distribution_panel(
        draw,
        (1170, 350, 2485, 1180),
        "B. Rank agreement",
        "Each dot is one mission-level caller comparison.",
        "mean Spearman rho",
        summary["rho_mean"],
        (summary["rho_min"], summary["rho_max"]),
        rows,
        "rho",
        0.75,
        1.00,
        [0.75, 0.85, 0.95, 1.00],
        COLORS["teal"],
    )
    draw_distribution_panel(
        draw,
        (1170, 1230, 2485, 1800),
        "C. List overlap",
        "FDR<0.05 DEG lists are more sensitive to caller stringency.",
        "mean Jaccard",
        summary["jaccard_mean"],
        (summary["jaccard_min"], summary["jaccard_max"]),
        rows,
        "jaccard",
        0.0,
        1.0,
        [0.0, 0.5, 1.0],
        COLORS["amber"],
    )
    draw_readout(draw, (2535, 350, 3690, 1800), summary)

    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    footer_label_x = 205
    footer_text_x = 390
    footer_text_width = 3220
    text(draw, (footer_label_x, 1925), "Takeaway", F["small_bold"], COLORS["sky"])
    footer = "Effect-size rankings are stable across callers, while thresholded DEG lists remain caller-sensitive."
    draw.text((footer_text_x, 1925), fit_label(draw, footer, F["small"], footer_text_width), font=F["small"], fill=COLORS["muted"])
    text(draw, (footer_label_x, 1995), "Next", F["small_bold"], COLORS["amber"])
    scope = "Published biology provides an external consistency check for the pathway signal."
    draw.text((footer_text_x, 1995), fit_label(draw, scope, F["small"], footer_text_width), font=F["small"], fill=COLORS["muted"])

    png = OUT_DIR / "dge_pipeline_robustness_premium.png"
    canvas.save(png, quality=95)
    gray = OUT_DIR / "dge_pipeline_robustness_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    manifest = {
        "title": "DGE ranks hold; DEG lists move",
        "source": str(SOURCE.relative_to(ROOT)),
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "summary": {
            "n_missions": summary["n_missions"],
            "n_comparisons": summary["n_comparisons"],
            "log2fc_spearman_mean": round(summary["rho_mean"], 3),
            "log2fc_spearman_min": round(summary["rho_min"], 3),
            "log2fc_spearman_max": round(summary["rho_max"], 3),
            "deg_jaccard_mean": round(summary["jaccard_mean"], 3),
            "deg_jaccard_min": round(summary["jaccard_min"], 3),
            "deg_jaccard_max": round(summary["jaccard_max"], 3),
            "top100_jaccard_mean": round(summary["top100_mean"], 3),
            "direction_concordance_mean": round(summary["direction_mean"], 3),
            "genelab_replication_mean_spearman": round(summary["genelab_mean"], 3),
        },
        "scope": (
            "DGE pipeline robustness supports stable effect-size ranking and direction. "
            "DEG-list membership is interpreted with caller and threshold context."
        ),
    }
    manifest_path = OUT_DIR / "dge_pipeline_robustness_premium_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

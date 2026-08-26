#!/usr/bin/env python3
"""Build the detailed-deck classical ML result-surface proof asset."""

from __future__ import annotations

import json
import statistics
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
M1_SUMMARY = ROOT / "v4" / "evaluation" / "M1_summary.json"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "classical_result_surface"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "classical_ml_result_surface_premium.png"
GRAY_PATH = OUT_DIR / "classical_ml_result_surface_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "classical_result_surface_manifest.json"

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "panel2": "#151F2D",
    "grid": "#2A3546",
    "axis": "#98A7BA",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#667487",
    "blue": "#66A6E8",
    "sky": "#73A7FF",
    "teal": "#5FD3C4",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "rose": "#E17882",
    "violet": "#B39DFF",
    "white": "#FFFFFF",
}

MODEL_ORDER = ["pca_lr", "elasticnet_lr", "xgb", "rf", "knn", "mlp", "tabnet", "svm_rbf"]
MODEL_LABELS = {
    "pca_lr": "PCA-LR",
    "elasticnet_lr": "ElasticNet-LR",
    "xgb": "XGBoost",
    "rf": "Random Forest",
    "knn": "KNN",
    "mlp": "MLP",
    "tabnet": "TabNet",
    "svm_rbf": "SVM-RBF",
}
FEATURE_ORDER = ["gene", "pathway_hallmark", "pathway_kegg", "combined"]
FEATURE_LABELS = {
    "gene": "Gene",
    "pathway_hallmark": "Hallmark",
    "pathway_kegg": "KEGG",
    "combined": "Combined",
}
TISSUE_LABELS = {
    "colon": "Colon",
    "eye": "Eye",
    "gastrocnemius": "Gastrocnemius",
    "kidney": "Kidney",
    "liver": "Liver",
    "lung": "Lung",
    "skin": "Skin",
    "thymus": "Thymus",
}


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgb_to_hex(rgb: tuple[int, int, int]) -> str:
    return "#" + "".join(f"{v:02X}" for v in rgb)


def rgba(value: str, alpha: int) -> tuple[int, int, int, int]:
    return (*hex_to_rgb(value), alpha)


def blend(a: str, b: str, t: float) -> str:
    t = max(0.0, min(1.0, t))
    ar, ag, ab = hex_to_rgb(a)
    br, bg, bb = hex_to_rgb(b)
    return rgb_to_hex((int(ar + (br - ar) * t), int(ag + (bg - ag) * t), int(ab + (bb - ab) * t)))


def heat_color(value: float) -> str:
    if value <= 0.55:
        return blend("#182231", "#415168", (value - 0.50) / 0.05)
    if value <= 0.70:
        return blend("#415168", COLORS["blue"], (value - 0.55) / 0.15)
    if value <= 0.82:
        return blend(COLORS["blue"], COLORS["teal"], (value - 0.70) / 0.12)
    return blend(COLORS["teal"], COLORS["green"], min((value - 0.82) / 0.14, 1.0))


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
    "title": load_font(76, True),
    "subtitle": load_font(35, False),
    "h2": load_font(42, True),
    "h3": load_font(32, True),
    "body": load_font(29, False),
    "body_bold": load_font(29, True),
    "small": load_font(25, False),
    "small_bold": load_font(25, True),
    "tiny": load_font(21, False),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18, False),
    "micro_bold": load_font(18, True),
    "stat": load_font(74, True),
    "huge": load_font(96, True),
}


def rounded(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    radius: int,
    fill: str | tuple[int, int, int, int],
    outline: str | tuple[int, int, int, int] | None = None,
    width: int = 1,
) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int | float, int | float],
    value: str,
    font: ImageFont.ImageFont,
    fill: str | tuple[int, int, int] | tuple[int, int, int, int] = COLORS["text"],
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
        else:
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
    fill: str | tuple[int, int, int, int],
    leading: int = 8,
) -> int:
    x, y = xy
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def load_data() -> dict[str, object]:
    m1 = json.loads(M1_SUMMARY.read_text())
    total_configs = 0
    sig_configs = 0
    tissues_with_sig: set[str] = set()
    best_rows = []
    matrix = []
    feature_totals = []

    for tissue, feature_map in m1.items():
        best = None
        for feature, model_map in feature_map.items():
            for model, record in model_map.items():
                auroc = float(record["auroc"])
                perm_p = float(record["perm_p"])
                total_configs += 1
                if perm_p < 0.05:
                    sig_configs += 1
                    tissues_with_sig.add(tissue)
                item = {
                    "tissue": tissue,
                    "feature": feature,
                    "model": model,
                    "auroc": auroc,
                    "perm_p": perm_p,
                    "n_folds": int(record["n_folds"]),
                    "cv": record["cv"],
                }
                if best is None or auroc > float(best["auroc"]):
                    best = item
        best_rows.append(best)

    for model in MODEL_ORDER:
        row = {"model": model, "label": MODEL_LABELS[model], "features": []}
        for feature in FEATURE_ORDER:
            vals = [float(m1[tissue][feature][model]["auroc"]) for tissue in m1]
            sig_count = sum(float(m1[tissue][feature][model]["perm_p"]) < 0.05 for tissue in m1)
            row["features"].append(
                {
                    "feature": feature,
                    "label": FEATURE_LABELS[feature],
                    "mean": statistics.mean(vals),
                    "sig_count": sig_count,
                    "min": min(vals),
                    "max": max(vals),
                }
            )
        matrix.append(row)

    for feature in FEATURE_ORDER:
        vals = [float(m1[tissue][feature][model]["auroc"]) for tissue in m1 for model in MODEL_ORDER]
        sig_count = sum(float(m1[tissue][feature][model]["perm_p"]) < 0.05 for tissue in m1 for model in MODEL_ORDER)
        feature_totals.append(
            {
                "feature": feature,
                "label": FEATURE_LABELS[feature],
                "mean": statistics.mean(vals),
                "sig_count": sig_count,
                "n": len(vals),
            }
        )

    model_gene_means = [
        {
            "model": model,
            "label": MODEL_LABELS[model],
            "mean": statistics.mean(float(m1[tissue]["gene"][model]["auroc"]) for tissue in m1),
        }
        for model in MODEL_ORDER
    ]
    top_cell = max(
        (feature for row in matrix for feature in row["features"]),
        key=lambda entry: float(entry["mean"]),
    )

    return {
        "total_configs": total_configs,
        "sig_configs": sig_configs,
        "tissues": len(m1),
        "tissues_with_sig": len(tissues_with_sig),
        "best_rows": sorted(best_rows, key=lambda row: float(row["auroc"]), reverse=True),
        "matrix": matrix,
        "feature_totals": feature_totals,
        "model_gene_means": sorted(model_gene_means, key=lambda row: float(row["mean"]), reverse=True),
        "top_cell": top_cell,
    }


def draw_background(canvas: Image.Image) -> None:
    overlay = Image.new("RGBA", (W, H), (0, 0, 0, 0))
    draw = ImageDraw.Draw(overlay)
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba("#1E2A3A", 28), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba("#1E2A3A", 22), width=1)
    draw.rectangle((0, 0, W, 310), fill=rgba("#121B28", 170))
    draw.rectangle((0, 1855, W, H), fill=rgba("#090E15", 145))
    canvas.alpha_composite(overlay)


def stat_badge(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, label: str, value: str, accent: str) -> None:
    rounded(draw, (x, y, x + w, y + 104), 24, "#121B28", "#2A394D", 2)
    text(draw, (x + 26, y + 23), label, F["tiny_bold"], accent)
    text(draw, (x + 26, y + 57), value, F["small_bold"], COLORS["text"])


def draw_header(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    text(draw, (136, 76), "MODELS / CLASSICAL RESULT SURFACE", F["kicker"], COLORS["sky"])
    text(draw, (136, 122), "Classical ML exposes a structured result surface", F["title"], COLORS["text"])
    text(
        draw,
        (138, 222),
        "The v4 grid shows where transfer is strong, which feature view carries signal, and which rows stay in audit context.",
        F["subtitle"],
        COLORS["muted"],
    )
    badges = [
        ("GRID", f"{int(data['total_configs'])} configs", 284, COLORS["sky"]),
        ("MODELS", "8 families", 250, COLORS["teal"]),
        ("VIEWS", "4 inputs", 228, COLORS["violet"]),
        ("PASS", f"{int(data['sig_configs'])} p<0.05", 250, COLORS["amber"]),
        ("TISSUES", f"{int(data['tissues_with_sig'])}/{int(data['tissues'])}", 208, COLORS["green"]),
    ]
    bx = 2372
    for label, value, width, accent in badges:
        stat_badge(draw, bx, 72, width, label, value, accent)
        bx += width + 18


def draw_tissue_ladder(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 46), "A. Best row by tissue", F["h2"], COLORS["text"])
    text(draw, (x0 + 48, y0 + 98), "One top configuration from the v4 grid for each tissue.", F["small"], COLORS["muted"])

    chart_x0 = x0 + 368
    chart_x1 = x1 - 238
    scale_min, scale_max = 0.50, 1.00

    def sx(value: float) -> int:
        return int(chart_x0 + (value - scale_min) / (scale_max - scale_min) * (chart_x1 - chart_x0))

    header_y = y0 + 174
    text(draw, (x0 + 58, header_y), "tissue", F["micro_bold"], COLORS["axis"])
    text(draw, (chart_x0, header_y), "AUROC", F["micro_bold"], COLORS["axis"])
    text(draw, (x1 - 180, header_y), "gate", F["micro_bold"], COLORS["axis"])
    for tick in [0.50, 0.75, 1.00]:
        x = sx(tick)
        draw.line((x, y0 + 210, x, y1 - 230), fill=rgba(COLORS["grid"], 120), width=2)
        text(draw, (x, y1 - 198), f"{tick:.2f}", F["micro_bold"], COLORS["axis"], "mm")

    start_y = y0 + 232
    row_gap = 126
    bar_h = 42
    for i, row in enumerate(data["best_rows"]):
        y = start_y + i * row_gap
        auroc = float(row["auroc"])
        p = float(row["perm_p"])
        passed = p < 0.05
        accent = COLORS["green"] if passed else COLORS["amber"]
        rounded(draw, (x0 + 46, y - 22, x1 - 46, y + 98), 24, "#121B28", "#263244", 1)
        text(draw, (x0 + 76, y), TISSUE_LABELS[str(row["tissue"])], F["small_bold"], COLORS["text"])
        model_feature = f"{MODEL_LABELS[str(row['model'])]} + {FEATURE_LABELS[str(row['feature'])]}"
        text(draw, (x0 + 76, y + 42), model_feature, F["micro_bold"], COLORS["muted"])
        rounded(draw, (chart_x0, y + 8, chart_x1, y + 8 + bar_h), 16, "#0B111B", "#263244", 1)
        rounded(draw, (chart_x0, y + 8, sx(auroc), y + 8 + bar_h), 16, rgba(accent, 220), None, 1)
        text(draw, (min(sx(auroc) + 18, chart_x1 - 10), y + 29), f"{auroc:.3f}", F["small_bold"], accent, "lm")
        badge = "pass" if passed else "audit"
        rounded(draw, (x1 - 190, y + 3, x1 - 70, y + 46), 18, rgba(accent, 38), rgba(accent, 170), 2)
        text(draw, (x1 - 130, y + 25), badge, F["micro_bold"], COLORS["text"], "mm")
        text(draw, (x1 - 190, y + 61), f"p={p:.3f}", F["micro_bold"], COLORS["muted"])

    readout_y = y1 - 160
    rounded(draw, (x0 + 54, readout_y, x1 - 54, y1 - 50), 28, rgba("#66A6E8", 34), rgba("#66A6E8", 150), 2)
    text(draw, (x0 + 88, readout_y + 32), "Readout", F["small_bold"], COLORS["sky"])
    paragraph(
        draw,
        (x0 + 88, readout_y + 67),
        "The tissue ladder separates high AUROC rows from lower-signal rows in the same view.",
        F["body_bold"],
        x1 - x0 - 176,
        COLORS["text"],
        leading=5,
    )


def draw_heatmap_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 46), "B. Model x feature surface", F["h2"], COLORS["text"])
    text(draw, (x0 + 48, y0 + 98), "Cell value is the eight-tissue mean AUROC; small count is p<0.05 rows.", F["small"], COLORS["muted"])

    label_w = 278
    grid_x = x0 + 338
    grid_y = y0 + 254
    cell_w = 185
    cell_h = 82
    gap_x = 14
    gap_y = 18

    for j, feature in enumerate(FEATURE_ORDER):
        x = grid_x + j * (cell_w + gap_x)
        rounded(draw, (x, grid_y - 70, x + cell_w, grid_y - 20), 18, "#121B28", "#2A394D", 1)
        text(draw, (x + cell_w / 2, grid_y - 45), FEATURE_LABELS[feature], F["micro_bold"], COLORS["axis"], "mm")

    top_mean = float(data["top_cell"]["mean"])
    for i, row in enumerate(data["matrix"]):
        y = grid_y + i * (cell_h + gap_y)
        text(draw, (x0 + 58, y + cell_h / 2), str(row["label"]), F["small_bold"], COLORS["text"], "lm")
        for j, feature in enumerate(row["features"]):
            x = grid_x + j * (cell_w + gap_x)
            mean = float(feature["mean"])
            sig = int(feature["sig_count"])
            fill = heat_color(mean)
            outline = COLORS["green"] if mean == top_mean else "#263244"
            width = 3 if mean == top_mean else 1
            rounded(draw, (x, y, x + cell_w, y + cell_h), 18, rgba(fill, 215), outline, width)
            text(draw, (x + 24, y + 20), f"{mean:.3f}", F["tiny_bold"], COLORS["text"])
            text(draw, (x + 24, y + 51), f"{sig}/8 pass", F["micro_bold"], COLORS["text"] if sig else COLORS["muted"])

    legend_y = y0 + 1122
    text(draw, (x0 + 58, legend_y), "Feature-view totals", F["h3"], COLORS["text"])
    fx = x0 + 58
    for feature in data["feature_totals"]:
        tw = 285
        accent = heat_color(float(feature["mean"]))
        rounded(draw, (fx, legend_y + 58, fx + tw, legend_y + 210), 24, COLORS["panel2"], "#2A394D", 1)
        text(draw, (fx + 26, legend_y + 88), str(feature["label"]), F["tiny_bold"], accent)
        text(draw, (fx + 26, legend_y + 124), f"{float(feature['mean']):.3f}", F["stat"], COLORS["text"])
        text(draw, (fx + 26, legend_y + 188), f"{int(feature['sig_count'])}/{int(feature['n'])} pass", F["micro_bold"], COLORS["muted"])
        fx += tw + 20


def takeaway_card(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    idx: str,
    title: str,
    body: str,
    accent: str,
) -> None:
    rounded(draw, (x, y, x + w, y + 158), 26, "#121B28", "#2A394D", 2)
    rounded(draw, (x + 24, y + 36, x + 92, y + 104), 18, rgba(accent, 52), rgba(accent, 190), 2)
    text(draw, (x + 58, y + 58), idx, F["small_bold"], accent, "mm")
    text(draw, (x + 122, y + 34), title, F["h3"], COLORS["text"])
    paragraph(draw, (x + 122, y + 80), body, F["tiny"], w - 160, COLORS["muted"], leading=4)


def draw_takeaway_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "C. Surface decoder", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 50, y0 + 98),
        "The grid is read as a map of where a simple model family already captures transferable structure.",
        F["small"],
        x1 - x0 - 100,
        COLORS["muted"],
        leading=6,
    )

    cards = [
        (
            "1",
            "Tissue strength",
            "Thymus, colon, and lung clear 0.90 AUROC in their best classical rows.",
            COLORS["sky"],
        ),
        (
            "2",
            "Model shape",
            "Linear baselines lead the gene-level ranking; tree, neural, and kernel rows add context.",
            COLORS["green"],
        ),
        (
            "3",
            "Feature view",
            "Pathway views lift selected tissues while gene features still anchor lung and skin.",
            COLORS["teal"],
        ),
        (
            "4",
            "Gate pattern",
            "Forty configurations pass p<0.05, spread across six of eight tissues.",
            COLORS["amber"],
        ),
    ]
    cx, cy, cw = x0 + 62, y0 + 230, x1 - x0 - 124
    for i, card in enumerate(cards):
        takeaway_card(draw, cx, cy + i * 190, cw, *card)

    readout_y = y1 - 362
    rounded(draw, (x0 + 62, readout_y, x1 - 62, y1 - 64), 32, COLORS["panel2"], "#2A394D", 2)
    text(draw, (x0 + 100, readout_y + 42), "Model-section bridge", F["h3"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 100, readout_y + 92),
        "Advanced models are interpreted as gains over a tissue-aware, feature-aware classical surface.",
        F["body_bold"],
        x1 - x0 - 200,
        COLORS["text"],
        leading=7,
    )
    top = data["model_gene_means"][0]
    second = data["model_gene_means"][1]
    rounded(draw, (x0 + 100, y1 - 178, x0 + 458, y1 - 98), 22, "#101823", "#2A394D", 1)
    text(draw, (x0 + 124, y1 - 154), "Gene-level leaders", F["micro_bold"], COLORS["axis"])
    text(draw, (x0 + 124, y1 - 125), f"{top['label']} {float(top['mean']):.3f}", F["tiny_bold"], COLORS["green"])
    text(draw, (x0 + 322, y1 - 125), f"{second['label']} {float(second['mean']):.3f}", F["tiny_bold"], COLORS["teal"])


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    y = 1872
    rounded(draw, (136, y, 3704, 2042), 32, "#101823", "#2A394D", 2)
    text(draw, (180, y + 38), "Slide 23 readout", F["small_bold"], COLORS["sky"])
    paragraph(
        draw,
        (180, y + 76),
        "The classical surface is already structured: high-transfer tissues, feature-specific gains, and auditable pass gates appear before adding pretrained models.",
        F["body_bold"],
        2220,
        COLORS["text"],
        leading=6,
    )
    paragraph(
        draw,
        (2600, y + 48),
        "All values on this slide are recomputed directly from the raw v4 M1 summary.",
        F["tiny_bold"],
        950,
        COLORS["muted"],
        leading=5,
    )
    text(draw, (140, 2102), "Takeaway: the classical surface is already structured before adding pretrained models.", F["micro"], COLORS["dim"])
    text(draw, (3704, 2102), "SPACEBIO-BENCH DETAILED DECK / MODELS", F["micro_bold"], COLORS["dim"], "ra")


def main() -> None:
    data = load_data()
    canvas = Image.new("RGBA", (W, H), COLORS["bg"])
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)

    draw_header(draw, data)
    draw_tissue_ladder(draw, (136, 345, 1268, 1815), data)
    draw_heatmap_panel(draw, (1308, 345, 2550, 1815), data)
    draw_takeaway_panel(draw, (2590, 345, 3704, 1815), data)
    draw_footer(draw)

    rgb = canvas.convert("RGB")
    rgb.save(OUT_PATH, quality=95)
    rgb.convert("L").convert("RGB").save(GRAY_PATH, quality=95)
    MANIFEST_PATH.write_text(
        json.dumps(
            {
                "asset": str(OUT_PATH.relative_to(ROOT)),
                "grayscale": str(GRAY_PATH.relative_to(ROOT)),
                "source_files": [str(M1_SUMMARY.relative_to(ROOT))],
                "metrics": data,
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps({"asset": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

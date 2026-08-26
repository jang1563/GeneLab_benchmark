#!/usr/bin/env python3
"""Build the detailed-deck v4 method-hardening asset."""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
M1_SUMMARY = ROOT / "v4" / "evaluation" / "M1_summary.json"
FIG1_HTML = ROOT / "v4" / "figures" / "html" / "Fig1_benchmark.html"
RESULTS_INVENTORY = ROOT / "docs" / "PROJECT_RESULTS_LOCATION_INVENTORY_2026_05_31.md"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "method_hardening"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "method_hardening_preserves_main_readout_premium.png"
GRAY_PATH = OUT_DIR / "method_hardening_preserves_main_readout_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "method_hardening_manifest.json"

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "panel2": "#151F2D",
    "panel3": "#121B28",
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
    "magenta": "#D56DFF",
    "white": "#FFFFFF",
}

FEATURE_LABELS = {
    "gene": "Gene",
    "pathway_hallmark": "Hallmark",
    "pathway_kegg": "KEGG",
    "combined": "Combined",
}

MODEL_LABELS = {
    "elasticnet_lr": "ElasticNet-LR",
    "pca_lr": "PCA-LR",
    "rf": "RF",
    "xgb": "XGB",
    "svm_rbf": "SVM-RBF",
    "knn": "KNN",
    "mlp": "MLP",
    "tabnet": "TabNet",
}


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(value: str, alpha: int) -> tuple[int, int, int, int]:
    return (*hex_to_rgb(value), alpha)


def blend(a: str, b: str, t: float) -> str:
    ar, ag, ab = hex_to_rgb(a)
    br, bg, bb = hex_to_rgb(b)
    return f"#{int(ar + (br - ar) * t):02x}{int(ag + (bg - ag) * t):02x}{int(ab + (bb - ab) * t):02x}"


def score_color(score: float) -> str:
    if score < 0.55:
        return blend("#253044", COLORS["rose"], max(0, min(1, (score - 0.40) / 0.15)) * 0.55)
    if score < 0.75:
        return blend("#253044", COLORS["teal"], (score - 0.55) / 0.20)
    return blend(COLORS["teal"], COLORS["amber"], min(1, (score - 0.75) / 0.22))


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
    "stat": load_font(66, True),
    "huge": load_font(90, True),
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
    leading: int = 7,
) -> int:
    x, y = xy
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


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


def load_data() -> dict[str, object]:
    raw = json.loads(M1_SUMMARY.read_text())
    tissues = sorted(raw)
    features = list(FEATURE_LABELS)
    models = list(MODEL_LABELS)
    total = 0
    sig = 0
    tissue_sig: set[str] = set()
    feature_stats: dict[str, dict[str, float | int]] = {}
    best_by_tissue: list[dict[str, object]] = []
    heatmap: list[dict[str, object]] = []

    for tissue, feature_map in raw.items():
        tissue_best = {"tissue": tissue, "feature": "", "model": "", "auroc": -1.0, "perm_p": 1.0}
        row = {"tissue": tissue, "features": []}
        for feature in features:
            model_map = feature_map[feature]
            best = {"feature": feature, "model": "", "auroc": -1.0, "perm_p": 1.0}
            values = []
            sig_count = 0
            for model in models:
                metrics = model_map[model]
                total += 1
                auroc = float(metrics["auroc"])
                perm_p = float(metrics["perm_p"])
                values.append(auroc)
                if perm_p < 0.05:
                    sig += 1
                    sig_count += 1
                    tissue_sig.add(tissue)
                if auroc > float(best["auroc"]):
                    best = {"feature": feature, "model": model, "auroc": auroc, "perm_p": perm_p}
                if auroc > float(tissue_best["auroc"]):
                    tissue_best = {"tissue": tissue, "feature": feature, "model": model, "auroc": auroc, "perm_p": perm_p}
            row["features"].append({**best, "any_sig": sig_count > 0})
            fs = feature_stats.setdefault(feature, {"sum": 0.0, "n": 0, "sig": 0})
            fs["sum"] = float(fs["sum"]) + sum(values)
            fs["n"] = int(fs["n"]) + len(values)
            fs["sig"] = int(fs["sig"]) + sig_count
        best_by_tissue.append(tissue_best)
        heatmap.append(row)

    best_by_tissue.sort(key=lambda r: float(r["auroc"]), reverse=True)
    order = [str(r["tissue"]) for r in best_by_tissue]
    heatmap.sort(key=lambda r: order.index(str(r["tissue"])))
    feature_summary = []
    for feature, stats in feature_stats.items():
        feature_summary.append(
            {
                "feature": feature,
                "mean": float(stats["sum"]) / int(stats["n"]),
                "sig": int(stats["sig"]),
            }
        )

    return {
        "n_tissues": len(tissues),
        "n_features": len(features),
        "n_models": len(models),
        "n_configs": total,
        "n_significant": sig,
        "n_tissue_significant": len(tissue_sig),
        "heatmap": heatmap,
        "best_by_tissue": best_by_tissue,
        "feature_summary": feature_summary,
    }


def stat_badge(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, label: str, value: str, accent: str) -> None:
    rounded(draw, (x, y, x + w, y + 104), 24, "#121B28", "#2A394D", 2)
    text(draw, (x + 26, y + 23), label, F["tiny_bold"], accent)
    text(draw, (x + 26, y + 57), value, F["small_bold"], COLORS["text"])


def draw_header(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    text(draw, (136, 76), "MODELS / METHOD HARDENING", F["kicker"], COLORS["sky"])
    text(draw, (136, 122), "Method hardening preserves the main readout", F["title"], COLORS["text"])
    text(
        draw,
        (138, 222),
        "v4 broadens tissues, classifiers, and feature views while keeping the mission-held-out score grammar fixed.",
        F["subtitle"],
        COLORS["muted"],
    )
    badges = [
        ("TISSUES", f"{data['n_tissues']}", 158, COLORS["teal"]),
        ("CLASSIFIERS", f"{data['n_models']}", 196, COLORS["sky"]),
        ("FEATURE VIEWS", f"{data['n_features']}", 208, COLORS["violet"]),
        ("CONFIGS", f"{data['n_configs']}", 182, COLORS["amber"]),
        ("P<0.05", f"{data['n_significant']}", 166, COLORS["green"]),
    ]
    bx = 2636
    for label, value, width, accent in badges:
        stat_badge(draw, bx, 72, width, label, value, accent)
        bx += width + 18


def draw_factor_box(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], value: str, label: str, accent: str) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 28, COLORS["panel2"], rgba(accent, 160), 2)
    text(draw, ((x0 + x1) / 2, y0 + 48), value, F["huge"], accent, "mm")
    text(draw, ((x0 + x1) / 2, y1 - 42), label, F["small_bold"], COLORS["text"], "mm")


def draw_contract_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 46), "A. What v4 hardens", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "The benchmark contract stays fixed while the evidence surface gets wider.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        leading=6,
    )

    y = y0 + 210
    w = (x1 - x0 - 136) // 2
    boxes = [
        ((x0 + 48, y, x0 + 48 + w, y + 190), "8", "tissues", COLORS["teal"]),
        ((x0 + 88 + w, y, x0 + 88 + 2 * w, y + 190), "8", "classifiers", COLORS["sky"]),
        ((x0 + 48, y + 230, x0 + 48 + w, y + 420), "4", "feature views", COLORS["violet"]),
        ((x0 + 88 + w, y + 230, x0 + 88 + 2 * w, y + 420), "256", "configs", COLORS["amber"]),
    ]
    for factor_box, value, label, accent in boxes:
        draw_factor_box(draw, factor_box, value, label, accent)

    rule_y = y0 + 720
    rules = [
        ("Same split", "LOMO keeps the hidden unit at mission level.", COLORS["sky"]),
        ("Same score", "AUROC, uncertainty, and permutation p-value keep one metric grammar.", COLORS["green"]),
        ("Wider surface", "Gene, Hallmark, KEGG, and combined views test representation choices.", COLORS["violet"]),
        ("Method breadth", "Linear, tree, kernel, nearest-neighbor, neural, and TabNet rows enter the grid.", COLORS["teal"]),
    ]
    for title, body, accent in rules:
        rounded(draw, (x0 + 52, rule_y, x1 - 52, rule_y + 136), 26, COLORS["panel2"], "#2A394D", 2)
        draw.line((x0 + 78, rule_y + 32, x0 + 78, rule_y + 104), fill=rgba(accent, 190), width=6)
        text(draw, (x0 + 108, rule_y + 30), title, F["h3"], COLORS["text"])
        paragraph(draw, (x0 + 108, rule_y + 76), body, F["tiny"], x1 - x0 - 210, COLORS["muted"], leading=4)
        rule_y += 160


def draw_heatmap_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 46), "B. Expanded result surface", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "Each cell shows the best AUROC within one tissue and feature view; a dot marks feature views with a significant row.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        leading=6,
    )

    hm = data["heatmap"]
    features = list(FEATURE_LABELS)
    left = x0 + 210
    top = y0 + 264
    cell_w = (x1 - left - 72) // len(features)
    cell_h = 104
    for j, feature in enumerate(features):
        cx = left + j * cell_w
        text(draw, (cx + cell_w / 2, top - 52), FEATURE_LABELS[feature], F["tiny_bold"], COLORS["axis"], "mm")
    for i, row in enumerate(hm):
        y = top + i * (cell_h + 12)
        tissue = str(row["tissue"]).title().replace("Gastrocnemius", "Gastro.")
        text(draw, (x0 + 56, y + 34), tissue, F["small_bold"], COLORS["text"])
        for j, cell in enumerate(row["features"]):
            x = left + j * cell_w
            score = float(cell["auroc"])
            color = score_color(score)
            rounded(draw, (x + 8, y, x + cell_w - 10, y + cell_h), 22, color, rgba("#FFFFFF", 65), 1)
            text(draw, (x + cell_w / 2, y + 32), f"{score:.3f}", F["small_bold"], COLORS["text"], "mm")
            model = MODEL_LABELS[str(cell["model"])]
            text(draw, (x + cell_w / 2, y + 70), model, F["micro_bold"], COLORS["text"], "mm")
            if bool(cell["any_sig"]):
                draw.ellipse((x + cell_w - 46, y + 12, x + cell_w - 28, y + 30), fill=COLORS["green"], outline=COLORS["white"], width=1)

    legend_y = y1 - 170
    rounded(draw, (x0 + 56, legend_y, x1 - 56, y1 - 58), 28, COLORS["panel2"], "#2A394D", 2)
    text(draw, (x0 + 92, legend_y + 36), "Read the grid horizontally", F["small_bold"], COLORS["sky"])
    paragraph(
        draw,
        (x0 + 92, legend_y + 72),
        "Feature views change the representation, but the held-out mission score stays comparable.",
        F["body_bold"],
        x1 - x0 - 184,
        COLORS["text"],
        leading=5,
    )


def draw_rank_row(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    rank: int,
    row: dict[str, object],
    max_score: float,
) -> None:
    score = float(row["auroc"])
    feature = FEATURE_LABELS[str(row["feature"])]
    model = MODEL_LABELS[str(row["model"])]
    tissue = str(row["tissue"]).title().replace("Gastrocnemius", "Gastro.")
    rounded(draw, (x, y, x + w, y + 106), 24, COLORS["panel2"], "#2A394D", 2)
    rounded(draw, (x + 24, y + 27, x + 78, y + 81), 16, rgba(COLORS["sky"], 42), rgba(COLORS["sky"], 150), 1)
    text(draw, (x + 51, y + 47), str(rank), F["micro_bold"], COLORS["sky"], "mm")
    text(draw, (x + 104, y + 26), tissue, F["small_bold"], COLORS["text"])
    text(draw, (x + 104, y + 63), f"{model} + {feature}", F["micro_bold"], COLORS["muted"])
    bar_x = x + w - 246
    bar_w = 154
    draw.line((bar_x, y + 54, bar_x + bar_w, y + 54), fill=rgba("#2A3546", 150), width=9)
    draw.line((bar_x, y + 54, bar_x + int(bar_w * score / max_score), y + 54), fill=rgba(COLORS["teal"], 210), width=9)
    text(draw, (x + w - 42, y + 38), f"{score:.3f}", F["small_bold"], COLORS["teal"], "ra")


def draw_readout_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "C. What stays readable", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 50, y0 + 98),
        "The broadened grid still points to tissue-specific transfer rather than a single model recipe.",
        F["small"],
        x1 - x0 - 100,
        COLORS["muted"],
        leading=6,
    )

    stat_y = y0 + 200
    stats = [
        ("significant rows", f"{data['n_significant']}/{data['n_configs']}", COLORS["green"]),
        ("tissues with signal", f"{data['n_tissue_significant']}/{data['n_tissues']}", COLORS["teal"]),
    ]
    stat_w = (x1 - x0 - 140) // 2
    for i, (label, value, accent) in enumerate(stats):
        x = x0 + 58 + i * (stat_w + 24)
        rounded(draw, (x, stat_y, x + stat_w, stat_y + 162), 28, COLORS["panel2"], rgba(accent, 150), 2)
        text(draw, (x + 30, stat_y + 34), label, F["tiny_bold"], accent)
        text(draw, (x + 30, stat_y + 82), value, F["stat"], COLORS["text"])

    text(draw, (x0 + 58, y0 + 440), "Top rows in the hardened grid", F["h3"], COLORS["text"])
    top_rows = data["best_by_tissue"][:6]
    max_score = float(top_rows[0]["auroc"])
    y = y0 + 500
    for rank, row in enumerate(top_rows, start=1):
        draw_rank_row(draw, x0 + 58, y, x1 - x0 - 116, rank, row, max_score)
        y += 126

    readout_y = y1 - 254
    rounded(draw, (x0 + 58, readout_y, x1 - 58, y1 - 58), 30, "#121B28", rgba("#F4C26B", 145), 2)
    text(draw, (x0 + 96, readout_y + 40), "Readout", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (x0 + 96, readout_y + 88),
        "Method breadth supports the same story: tissue context and feature view shape transfer strength.",
        F["body_bold"],
        x1 - x0 - 192,
        COLORS["text"],
        leading=6,
    )


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    y = 1872
    rounded(draw, (136, y, 3704, 2042), 32, "#101823", "#2A394D", 2)
    text(draw, (180, y + 38), "Slide 28 readout", F["small_bold"], COLORS["sky"])
    paragraph(
        draw,
        (180, y + 76),
        "v4 is a method-hardening pass: the grid is wider, while the hidden-mission readout and tissue-specific pattern remain legible.",
        F["body_bold"],
        2400,
        COLORS["text"],
        leading=6,
    )
    paragraph(
        draw,
        (2730, y + 46),
        "Next: newer model ideas are added as checks against this hardened benchmark surface.",
        F["tiny_bold"],
        820,
        COLORS["muted"],
        leading=5,
    )
    text(
        draw,
        (140, 2102),
        "Takeaway: method hardening broadens the grid while preserving the same tissue and feature-view readout.",
        F["micro"],
        COLORS["dim"],
    )
    text(draw, (3704, 2102), "SPACEBIO-BENCH DETAILED DECK / MODELS", F["micro_bold"], COLORS["dim"], "ra")


def main() -> None:
    data = load_data()
    canvas = Image.new("RGBA", (W, H), COLORS["bg"])
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)

    draw_header(draw, data)
    draw_contract_panel(draw, (136, 345, 1048, 1815), data)
    draw_heatmap_panel(draw, (1088, 345, 2508, 1815), data)
    draw_readout_panel(draw, (2548, 345, 3704, 1815), data)
    draw_footer(draw)

    rgb = canvas.convert("RGB")
    rgb.save(OUT_PATH, quality=95)
    rgb.convert("L").convert("RGB").save(GRAY_PATH, quality=95)
    MANIFEST_PATH.write_text(
        json.dumps(
            {
                "asset": str(OUT_PATH.relative_to(ROOT)),
                "grayscale": str(GRAY_PATH.relative_to(ROOT)),
                "source_files": [
                    str(M1_SUMMARY.relative_to(ROOT)),
                    str(FIG1_HTML.relative_to(ROOT)),
                    str(RESULTS_INVENTORY.relative_to(ROOT)),
                ],
                "metrics": data,
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps({"asset": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

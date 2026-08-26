#!/usr/bin/env python3
"""Build the detailed-deck strong-baselines proof asset."""

from __future__ import annotations

import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
M1_SUMMARY = ROOT / "v4" / "evaluation" / "M1_summary.json"
GENEFORMER_SUMMARY = ROOT / "evaluation" / "geneformer_mouse_gf_all_tissues_summary.json"
SCGPT_SUMMARY = ROOT / "evaluation" / "scgpt_whole_human_all_tissues_summary.json"
V9_BULK_BASELINES = ROOT / "v9" / "reports" / "bulk_lomo_baseline_summary.csv"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "strong_baselines"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "strong_baselines_make_model_claims_meaningful_premium.png"
GRAY_PATH = OUT_DIR / "strong_baselines_make_model_claims_meaningful_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "strong_baselines_manifest.json"

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "panel2": "#151F2D",
    "panel3": "#0F1722",
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
    model_order = ["pca_lr", "elasticnet_lr", "xgb", "rf", "knn", "mlp", "tabnet", "svm_rbf"]
    model_labels = {
        "pca_lr": "PCA-LR",
        "elasticnet_lr": "ElasticNet-LR",
        "xgb": "XGBoost",
        "rf": "Random Forest",
        "knn": "KNN",
        "mlp": "MLP",
        "tabnet": "TabNet",
        "svm_rbf": "SVM-RBF",
    }
    gene_means = []
    total_configs = 0
    sig_configs = 0
    tissues_with_sig: set[str] = set()
    for tissue, feature_map in m1.items():
        for feature_name, model_map in feature_map.items():
            for model, record in model_map.items():
                total_configs += 1
                perm_p = float(record["perm_p"])
                if perm_p < 0.05:
                    sig_configs += 1
                    tissues_with_sig.add(tissue)
    for model in model_order:
        vals = [float(m1[t]["gene"][model]["auroc"]) for t in m1]
        gene_means.append(
            {
                "model": model,
                "label": model_labels[model],
                "mean": statistics.mean(vals),
                "min": min(vals),
                "max": max(vals),
                "n": len(vals),
            }
        )
    gene_means = sorted(gene_means, key=lambda row: float(row["mean"]), reverse=True)

    geneformer = json.loads(GENEFORMER_SUMMARY.read_text())
    scgpt = json.loads(SCGPT_SUMMARY.read_text())
    gf_tasks = geneformer["tissues"]
    scgpt_tasks = scgpt["tasks"]
    gf_wins = sum(float(row["geneformer_mean_auroc"]) > float(row["baseline_mean_auroc"]) for row in gf_tasks.values())
    scgpt_wins = sum(float(row["mean_auroc"]) > float(row["baseline_auroc"]) for row in scgpt_tasks.values())

    v9_by_baseline: dict[str, list[float]] = defaultdict(list)
    with V9_BULK_BASELINES.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row["status"] == "evaluated" and row["auroc"]:
                v9_by_baseline[row["baseline_id"]].append(float(row["auroc"]))
    v9_rows = sum(len(vals) for vals in v9_by_baseline.values())
    v9_summary = {
        baseline_id: {
            "n": len(vals),
            "mean": statistics.mean(vals),
            "min": min(vals),
            "max": max(vals),
        }
        for baseline_id, vals in sorted(v9_by_baseline.items())
    }

    return {
        "gene_means": gene_means,
        "total_configs": total_configs,
        "sig_configs": sig_configs,
        "tissues": len(m1),
        "tissues_with_sig": len(tissues_with_sig),
        "baseline_mean": float(geneformer["overall"]["baseline_mean"]),
        "geneformer_mean": float(geneformer["overall"]["geneformer_mean"]),
        "geneformer_wins": gf_wins,
        "geneformer_n": int(geneformer["overall"]["n_compared"]),
        "scgpt_mean": float(scgpt["overall_mean_auroc"]),
        "scgpt_baseline_mean": float(scgpt["baseline_mean_auroc"]),
        "scgpt_delta": float(scgpt["delta_vs_baseline"]),
        "scgpt_wins": scgpt_wins,
        "scgpt_n": len(scgpt_tasks),
        "v9_rows": v9_rows,
        "v9_summary": v9_summary,
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
    text(draw, (136, 76), "MODELS / BASELINE FLOOR", F["kicker"], COLORS["sky"])
    text(draw, (136, 122), "Strong baselines make model comparisons meaningful", F["title"], COLORS["text"])
    text(
        draw,
        (138, 222),
        "Advanced models inherit the same task contract, then earn interpretation by clearing a strong classical floor.",
        F["subtitle"],
        COLORS["muted"],
    )
    badges = [
        ("GRID", f"{int(data['total_configs'])} configs", 284, COLORS["sky"]),
        ("TISSUES", f"{int(data['tissues'])} tasks", 250, COLORS["teal"]),
        ("PASS", f"{int(data['sig_configs'])} p<0.05", 250, COLORS["amber"]),
        ("TOP FLOOR", f"{float(data['gene_means'][0]['mean']):.3f}", 250, COLORS["green"]),
    ]
    bx = 2650
    for label, value, width, accent in badges:
        stat_badge(draw, bx, 72, width, label, value, accent)
        bx += width + 20


def step_card(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    idx: str,
    title: str,
    body: str,
    accent: str,
) -> None:
    rounded(draw, (x, y, x + w, y + 150), 26, "#121B28", "#2A394D", 2)
    rounded(draw, (x + 24, y + 34, x + 92, y + 102), 18, rgba(accent, 52), rgba(accent, 190), 2)
    text(draw, (x + 58, y + 55), idx, F["small_bold"], accent, "mm")
    text(draw, (x + 122, y + 34), title, F["h3"], COLORS["text"])
    paragraph(draw, (x + 122, y + 79), body, F["tiny"], w - 160, COLORS["muted"], leading=4)


def draw_contract_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 46), "A. Baseline contract", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "A model comparison starts from the same hidden-mission task, preprocessing frame, and score grammar.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        leading=6,
    )

    steps = [
        ("1", "Task is fixed", "Mission-held-out split defines train and test before modeling.", COLORS["sky"]),
        ("2", "Input surface is fixed", "Gene, pathway, or embedding views are compared inside the same task.", COLORS["teal"]),
        ("3", "Training choices stop early", "Transforms and feature selection use training missions only.", COLORS["amber"]),
        ("4", "Classical grid sets floor", "PCA-LR, ElasticNet, tree, kernel, and neural baselines run first.", COLORS["green"]),
        ("5", "Advanced models inherit it", "Pretrained models are read against the matched baseline floor.", COLORS["violet"]),
    ]
    sx, sy, sw = x0 + 56, y0 + 220, x1 - x0 - 112
    for i, (idx, title, body, accent) in enumerate(steps):
        cy = sy + i * 183
        step_card(draw, sx, cy, sw, idx, title, body, accent)
        if i < len(steps) - 1:
            cx = sx + sw // 2
            draw.line((cx, cy + 154, cx, cy + 178), fill=rgba("#98A7BA", 120), width=3)
            draw.polygon([(cx - 10, cy + 174), (cx + 10, cy + 174), (cx, cy + 188)], fill=rgba("#98A7BA", 160))

    readout_y = y1 - 178
    rounded(draw, (x0 + 56, readout_y, x1 - 56, y1 - 50), 30, rgba("#66A6E8", 34), rgba("#66A6E8", 150), 2)
    text(draw, (x0 + 90, readout_y + 33), "Readout", F["small_bold"], COLORS["sky"])
    paragraph(
        draw,
        (x0 + 90, readout_y + 70),
        "The baseline floor carries the model comparison.",
        F["body_bold"],
        x1 - x0 - 180,
        COLORS["text"],
        leading=5,
    )


def draw_classical_chart(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "B. Classical floor is nontrivial", F["h2"], COLORS["text"])
    text(draw, (x0 + 50, y0 + 98), "Gene-level mean AUROC across the eight tissue tasks.", F["small"], COLORS["muted"])

    chart_x0, chart_y0 = x0 + 82, y0 + 214
    chart_x1, chart_y1 = x1 - 66, y0 + 1080
    label_w = 300
    axis_x0 = chart_x0 + label_w
    axis_x1 = chart_x1 - 102
    min_v, max_v = 0.50, 0.82

    def sx(value: float) -> int:
        return int(axis_x0 + (value - min_v) / (max_v - min_v) * (axis_x1 - axis_x0))

    for tick in [0.50, 0.60, 0.70, 0.80]:
        x = sx(tick)
        draw.line((x, chart_y0 - 22, x, chart_y1 + 16), fill=rgba(COLORS["grid"], 150), width=2)
        text(draw, (x, chart_y1 + 42), f"{tick:.2f}", F["micro_bold"], COLORS["axis"], "mm")
    draw.line((axis_x0, chart_y1 + 8, axis_x1, chart_y1 + 8), fill=rgba(COLORS["axis"], 120), width=2)
    text(draw, (sx(0.50) + 16, chart_y0 - 28), "chance", F["micro_bold"], COLORS["dim"])

    rows = data["gene_means"]
    bar_h = 56
    gap = 39
    colors = [COLORS["sky"], COLORS["teal"], COLORS["amber"], COLORS["green"], COLORS["blue"], COLORS["violet"], COLORS["rose"], COLORS["dim"]]
    for i, row in enumerate(rows):
        y = chart_y0 + i * (bar_h + gap)
        label = str(row["label"])
        value = float(row["mean"])
        color = colors[i]
        text(draw, (chart_x0, y + bar_h // 2), label, F["small_bold"], COLORS["text"], "lm")
        rounded(draw, (axis_x0, y, axis_x1, y + bar_h), 18, "#0D141F", "#263244", 1)
        rounded(draw, (axis_x0, y, sx(value), y + bar_h), 18, rgba(color, 210), None, 1)
        text(draw, (sx(value) + 22, y + bar_h // 2), f"{value:.3f}", F["small_bold"], color, "lm")

    callout_y = y0 + 1150
    rounded(draw, (x0 + 54, callout_y, x1 - 54, y1 - 64), 32, COLORS["panel2"], "#2A394D", 2)
    text(draw, (x0 + 92, callout_y + 42), "Result surface", F["h3"], COLORS["text"])
    tiles = [
        ("configured comparisons", f"{int(data['total_configs'])}", COLORS["sky"]),
        ("pass permutation gate", f"{int(data['sig_configs'])}", COLORS["amber"]),
        ("tissues with signal", f"{int(data['tissues_with_sig'])}/{int(data['tissues'])}", COLORS["green"]),
    ]
    tx = x0 + 92
    for label, value, accent in tiles:
        tw = 345
        rounded(draw, (tx, callout_y + 110, tx + tw, callout_y + 250), 24, "#101823", "#2A394D", 1)
        text(draw, (tx + 28, callout_y + 139), label, F["tiny_bold"], accent)
        text(draw, (tx + 28, callout_y + 177), value, F["stat"], COLORS["text"])
        tx += tw + 28


def metric_bar(
    draw: ImageDraw.ImageDraw,
    x0: int,
    y: int,
    w: int,
    label: str,
    value: float,
    accent: str,
    note: str,
    scale_min: float = 0.45,
    scale_max: float = 0.80,
) -> None:
    def sx(v: float) -> int:
        return int(x0 + (v - scale_min) / (scale_max - scale_min) * w)

    row_h = 128
    rounded(draw, (x0 - 12, y, x0 + w + 12, y + row_h), 26, "#101823", "#2A394D", 1)
    text(draw, (x0 + 28, y + 30), label, F["small_bold"], COLORS["text"])
    text(draw, (x0 + w - 26, y + 30), f"{value:.3f}", F["small_bold"], accent, "ra")
    track_y = y + 82
    rounded(draw, (x0 + 28, track_y, x0 + w - 28, track_y + 18), 9, "#0B111B", None, 1)
    inner_x0 = x0 + 28
    inner_w = w - 56

    def inner_sx(v: float) -> int:
        return int(inner_x0 + (v - scale_min) / (scale_max - scale_min) * inner_w)

    value_x = max(inner_x0 + 4, min(inner_sx(value), inner_x0 + inner_w))
    rounded(draw, (inner_x0, track_y, value_x, track_y + 18), 9, rgba(accent, 215), None, 1)
    text(draw, (x0 + 28, y + 105), note, F["micro_bold"], COLORS["muted"])


def draw_model_comparison(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "C. Baseline changes interpretation", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 50, y0 + 98),
        "The next model slides ask where pretraining clears a matched floor.",
        F["small"],
        x1 - x0 - 100,
        COLORS["muted"],
        leading=6,
    )

    chart_x, chart_y, chart_w = x0 + 66, y0 + 250, x1 - x0 - 132
    rows = [
        ("Matched baseline", float(data["baseline_mean"]), COLORS["green"], "same six tissue tasks"),
        ("scGPT whole-human", float(data["scgpt_mean"]), COLORS["amber"], f"wins {int(data['scgpt_wins'])}/{int(data['scgpt_n'])}; delta {float(data['scgpt_delta']):+.3f}"),
        ("Geneformer", float(data["geneformer_mean"]), COLORS["rose"], f"wins {int(data['geneformer_wins'])}/{int(data['geneformer_n'])}; delta {float(data['geneformer_mean']) - float(data['baseline_mean']):+.3f}"),
    ]
    for i, (label, value, accent, note) in enumerate(rows):
        metric_bar(draw, chart_x, chart_y + i * 168, chart_w, label, value, accent, note)

    scale_min, scale_max = 0.45, 0.80
    inner_x0 = chart_x + 28
    inner_w = chart_w - 56

    def inner_sx(v: float) -> int:
        return int(inner_x0 + (v - scale_min) / (scale_max - scale_min) * inner_w)

    chance_x = inner_sx(0.50)
    floor_x = inner_sx(float(data["baseline_mean"]))
    draw.line((chance_x, chart_y + 72, chance_x, chart_y + 468), fill=rgba(COLORS["axis"], 120), width=2)
    text(draw, (chance_x + 8, chart_y + 472), "0.50 chance", F["micro_bold"], COLORS["axis"])
    draw.line((floor_x, chart_y + 72, floor_x, chart_y + 468), fill=rgba(COLORS["green"], 170), width=3)
    text(draw, (floor_x - 10, chart_y + 472), "floor", F["micro_bold"], COLORS["green"], "ra")

    rule_y = y0 + 865
    rounded(draw, (x0 + 66, rule_y, x1 - 66, y1 - 66), 32, COLORS["panel2"], "#2A394D", 2)
    text(draw, (x0 + 104, rule_y + 44), "Interpretation rule", F["h3"], COLORS["text"])
    rules = [
        ("1", "Beat chance", "Score shows separability exists.", COLORS["sky"]),
        ("2", "Clear the floor", "Matched classical baseline earns model-level interpretation.", COLORS["green"]),
        ("3", "Localize wins", "Tissue-specific gains become follow-up hypotheses.", COLORS["violet"]),
    ]
    ry = rule_y + 118
    for idx, title, body, accent in rules:
        rounded(draw, (x0 + 104, ry, x0 + 172, ry + 68), 18, rgba(accent, 46), rgba(accent, 180), 2)
        text(draw, (x0 + 138, ry + 35), idx, F["tiny_bold"], accent, "mm")
        text(draw, (x0 + 200, ry + 1), title, F["small_bold"], COLORS["text"])
        paragraph(draw, (x0 + 200, ry + 36), body, F["tiny"], x1 - x0 - 340, COLORS["muted"], leading=3)
        ry += 122


def draw_footer(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    y = 1872
    rounded(draw, (136, y, 3704, 2042), 32, "#101823", "#2A394D", 2)
    text(draw, (180, y + 38), "Slide 22 readout", F["small_bold"], COLORS["sky"])
    paragraph(
        draw,
        (180, y + 76),
        "A model comparison is meaningful when the same hidden-mission contract already has a strong, auditable baseline floor.",
        F["body_bold"],
        2050,
        COLORS["text"],
        leading=6,
    )
    summary = f"v9 public scaffold: {int(data['v9_rows'])} evaluated baseline rows across nearest-centroid, PCA-LR, and L2 logistic surfaces."
    paragraph(draw, (2460, y + 40), summary, F["tiny_bold"], 1110, COLORS["muted"], leading=5)
    sources = "Takeaway: classical baselines define the floor that advanced model rows need to clear."
    text(draw, (140, 2102), sources, F["micro"], COLORS["dim"])
    text(draw, (3704, 2102), "SPACEBIO-BENCH DETAILED DECK / MODELS", F["micro_bold"], COLORS["dim"], "ra")


def main() -> None:
    data = load_data()
    canvas = Image.new("RGBA", (W, H), COLORS["bg"])
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)

    draw_header(draw, data)
    draw_contract_panel(draw, (136, 345, 1138, 1815))
    draw_classical_chart(draw, (1184, 345, 2518, 1815), data)
    draw_model_comparison(draw, (2564, 345, 3704, 1815), data)
    draw_footer(draw, data)

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
                    str(GENEFORMER_SUMMARY.relative_to(ROOT)),
                    str(SCGPT_SUMMARY.relative_to(ROOT)),
                    str(V9_BULK_BASELINES.relative_to(ROOT)),
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

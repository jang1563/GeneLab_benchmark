#!/usr/bin/env python3
"""Build the detailed-deck model-tier result asset."""

from __future__ import annotations

import glob
import json
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
GENEFORMER_SUMMARY = ROOT / "evaluation" / "geneformer_mouse_gf_all_tissues_summary.json"
SCGPT_SUMMARY = ROOT / "evaluation" / "scgpt_whole_human_all_tissues_summary.json"
RESULTS_SUMMARY = ROOT / "evaluation" / "RESULTS_SUMMARY.md"
CANONICAL_RESULTS = ROOT / "docs" / "CANONICAL_RESULTS_V7_1.md"
EVAL_GLOB = str(ROOT / "evaluation" / "submission_*_zeroshot_A*_eval.json")
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "model_tier_result"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "model_tier_result_task_fit_beats_scale_premium.png"
GRAY_PATH = OUT_DIR / "model_tier_result_task_fit_beats_scale_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "model_tier_result_manifest.json"

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
    "stat": load_font(66, True),
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


def provider_from_eval_path(path: Path) -> str:
    return path.stem[len("submission_") :].split("_zeroshot_")[0]


def load_data() -> dict[str, object]:
    gf = json.loads(GENEFORMER_SUMMARY.read_text())
    scgpt = json.loads(SCGPT_SUMMARY.read_text())

    llm_scores: dict[str, list[float]] = {}
    for path_s in glob.glob(EVAL_GLOB):
        path = Path(path_s)
        data = json.loads(path.read_text())
        provider = provider_from_eval_path(path)
        llm_scores.setdefault(provider, []).append(float(data["summary"]["mean_auroc"]))
    llm_means = {provider: sum(values) / len(values) for provider, values in llm_scores.items()}

    shared_rows = [
        {"label": "PCA-LR / classical ref", "surface": "tabular matrix", "score": float(gf["overall"]["baseline_mean"]), "color": COLORS["blue"]},
        {"label": "scGPT whole-human", "surface": "pretrained expression", "score": float(scgpt["overall_mean_auroc"]), "color": COLORS["teal"]},
        {"label": "Gemini-2.5-Flash", "surface": "prompt diagnostic", "score": llm_means["Gemini-2.5-Flash"], "color": COLORS["amber"]},
        {"label": "Llama-3.3-70B", "surface": "prompt diagnostic", "score": llm_means["Llama-3.3-70B"], "color": COLORS["violet"]},
        {"label": "Mouse-Geneformer", "surface": "pretrained expression", "score": float(gf["overall"]["geneformer_mean"]), "color": COLORS["rose"]},
        {"label": "DeepSeek-V3", "surface": "prompt diagnostic", "score": llm_means["DeepSeek-V3"], "color": COLORS["magenta"]},
    ]

    tissue_order = ["liver", "kidney", "skin", "eye", "thymus", "gastrocnemius"]
    task_by_tissue = {row["tissue"]: task for task, row in gf["tissues"].items()}
    deltas = []
    for tissue in tissue_order:
        task = task_by_tissue[tissue]
        sc_task = scgpt["tasks"][task]
        gf_task = gf["tissues"][task]
        deltas.append(
            {
                "tissue": tissue,
                "scgpt": float(sc_task["delta_vs_baseline"]),
                "geneformer": float(gf_task["delta"]),
                "baseline": float(gf_task["baseline_mean_auroc"]),
                "scgpt_auroc": float(sc_task["mean_auroc"]),
                "geneformer_auroc": float(gf_task["geneformer_mean_auroc"]),
            }
        )

    return {
        "shared_rows": shared_rows,
        "deltas": deltas,
        "n_tasks": int(scgpt["n_tasks"]),
        "scgpt_folds": int(scgpt["n_folds_total"]),
        "geneformer_folds": sum(int(row["geneformer_n_folds"]) for row in gf["tissues"].values()),
        "llm_evals": len(glob.glob(EVAL_GLOB)),
        "scgpt_local_gains": sum(1 for row in deltas if row["scgpt"] > 0),
        "geneformer_local_gains": sum(1 for row in deltas if row["geneformer"] > 0),
    }


def stat_badge(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, label: str, value: str, accent: str) -> None:
    rounded(draw, (x, y, x + w, y + 104), 24, "#121B28", "#2A394D", 2)
    text(draw, (x + 26, y + 23), label, F["tiny_bold"], accent)
    text(draw, (x + 26, y + 57), value, F["small_bold"], COLORS["text"])


def draw_header(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    text(draw, (136, 76), "MODELS / TIER RESULT", F["kicker"], COLORS["sky"])
    text(draw, (136, 122), "Task fit beats scale on the shared benchmark", F["title"], COLORS["text"])
    text(
        draw,
        (138, 222),
        "The model tier label is interpreted through the input surface, matched classical floor, and tissue-local deltas.",
        F["subtitle"],
        COLORS["muted"],
    )
    badges = [
        ("SHARED SURFACE", f"{data['n_tasks']} tasks", 246, COLORS["teal"]),
        ("FM FOLDS", f"{data['scgpt_folds']} + {data['geneformer_folds']}", 214, COLORS["violet"]),
        ("LLM EVALS", f"{data['llm_evals']} files", 198, COLORS["amber"]),
        ("LOCAL scGPT", f"{data['scgpt_local_gains']}/6 gains", 218, COLORS["green"]),
        ("REFERENCE", "0.758", 190, COLORS["blue"]),
    ]
    bx = 2440
    for label, value, width, accent in badges:
        stat_badge(draw, bx, 72, width, label, value, accent)
        bx += width + 18


def fmt3(value: float) -> str:
    return str(Decimal(str(value)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))


def draw_rule_card(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    w: int,
    title: str,
    body: str,
    accent: str,
    idx: str,
) -> None:
    rounded(draw, (x, y, x + w, y + 178), 28, COLORS["panel2"], "#2A394D", 2)
    rounded(draw, (x + 28, y + 42, x + 100, y + 114), 19, rgba(accent, 48), rgba(accent, 180), 2)
    text(draw, (x + 64, y + 68), idx, F["small_bold"], accent, "mm")
    text(draw, (x + 130, y + 36), title, F["h3"], COLORS["text"])
    paragraph(draw, (x + 130, y + 84), body, F["tiny"], w - 174, COLORS["muted"], leading=4)


def draw_comparison_rules(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 46), "A. Compare the right surfaces", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "The slide has two reading levels: shared rows for tier comparison, local deltas for tissue behavior.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        leading=6,
    )

    cards = [
        ("Shared row", "Six-task means compare classical ML, expression-pretrained models, and prompt diagnostics on the same task set.", COLORS["sky"], "1"),
        ("Local delta", "Each tissue asks whether the pretrained expression row clears its matched classical floor.", COLORS["teal"], "2"),
        ("Input surface", "Bulk matrix, pretrained expression tokens, and gene-text prompts carry different information into the gate.", COLORS["violet"], "3"),
        ("Reading order", "Read the score surface first; then interpret what the model tier contributes.", COLORS["green"], "4"),
    ]
    y = y0 + 220
    for title, body, accent, idx in cards:
        draw_rule_card(draw, x0 + 56, y, x1 - x0 - 112, title, body, accent, idx)
        y += 208

    readout_y = y1 - 238
    rounded(draw, (x0 + 56, readout_y, x1 - 56, y1 - 58), 30, "#121B28", rgba("#F4C26B", 145), 2)
    text(draw, (x0 + 92, readout_y + 38), "Result lens", F["h3"], COLORS["amber"])
    paragraph(
        draw,
        (x0 + 92, readout_y + 88),
        "Scale helps when representation, assay surface, species mapping, and task size align.",
        F["body_bold"],
        x1 - x0 - 184,
        COLORS["text"],
        leading=6,
    )


def map_x(value: float, x0: int, w: int, lo: float, hi: float) -> int:
    value = max(lo, min(hi, value))
    return int(x0 + (value - lo) / (hi - lo) * w)


def draw_shared_chart(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], rows: list[dict[str, object]]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 48, y0 + 46), "B. Shared six-task readout", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 48, y0 + 98),
        "Scores are mean AUROC / task score on the shared six-tissue surface.",
        F["small"],
        x1 - x0 - 96,
        COLORS["muted"],
        leading=6,
    )
    chart_x0 = x0 + 420
    chart_y0 = y0 + 245
    chart_w = x1 - chart_x0 - 112
    row_gap = 154
    lo, hi = 0.40, 0.82

    for val, label in [(0.40, "0.40"), (0.50, "chance"), (0.60, "0.60"), (0.70, "0.70"), (0.80, "0.80")]:
        px = map_x(val, chart_x0, chart_w, lo, hi)
        draw.line((px, chart_y0 - 54, px, chart_y0 + row_gap * 5 + 54), fill=rgba("#98A7BA", 115 if val == 0.5 else 58), width=2 if val == 0.5 else 1)
        text(draw, (px, chart_y0 - 86), label, F["micro_bold"], COLORS["axis"], "mm")

    surface_colors = {
        "tabular matrix": COLORS["blue"],
        "pretrained expression": COLORS["teal"],
        "prompt diagnostic": COLORS["amber"],
    }
    for i, row in enumerate(rows):
        y = chart_y0 + i * row_gap
        label = str(row["label"])
        surface = str(row["surface"])
        score = float(row["score"])
        color = str(row["color"])
        text(draw, (x0 + 56, y - 18), label, F["small_bold"], COLORS["text"])
        text(draw, (x0 + 56, y + 18), surface, F["micro_bold"], surface_colors.get(surface, COLORS["muted"]))
        px = map_x(score, chart_x0, chart_w, lo, hi)
        draw.line((chart_x0, y + 2, chart_x0 + chart_w, y + 2), fill=rgba("#2A3546", 150), width=10)
        draw.line((chart_x0, y + 2, px, y + 2), fill=rgba(color, 155), width=10)
        draw.ellipse((px - 22, y - 20, px + 22, y + 24), fill=rgba(color, 225), outline=COLORS["white"], width=2)
        text(draw, (chart_x0 + chart_w + 36, y - 18), fmt3(score), F["small_bold"], color)

    band_y = y1 - 178
    rounded(draw, (x0 + 56, band_y, x1 - 56, y1 - 58), 28, COLORS["panel2"], "#2A394D", 2)
    text(draw, (x0 + 92, band_y + 36), "Pattern", F["small_bold"], COLORS["sky"])
    paragraph(
        draw,
        (x0 + 92, band_y + 72),
        "The classical reference remains highest on this shared surface; prompt diagnostics cluster around the chance line.",
        F["body_bold"],
        x1 - x0 - 184,
        COLORS["text"],
        leading=6,
    )


def draw_delta_chart(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], deltas: list[dict[str, object]]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "C. Tissue-local deltas", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 50, y0 + 98),
        "Dots show AUROC change versus the matched classical floor for each tissue.",
        F["small"],
        x1 - x0 - 100,
        COLORS["muted"],
        leading=6,
    )

    chart_x0 = x0 + 256
    chart_y0 = y0 + 250
    chart_w = x1 - chart_x0 - 110
    row_gap = 142
    lo, hi = -0.56, 0.08
    zero_x = map_x(0.0, chart_x0, chart_w, lo, hi)

    for val, label in [(-0.50, "-0.50"), (-0.25, "-0.25"), (0.0, "matched floor")]:
        px = map_x(val, chart_x0, chart_w, lo, hi)
        draw.line((px, chart_y0 - 60, px, chart_y0 + row_gap * 5 + 50), fill=rgba("#98A7BA", 135 if val == 0 else 58), width=2 if val == 0 else 1)
        text(draw, (px, chart_y0 - 90), label, F["micro_bold"], COLORS["axis"], "mm")

    text(draw, (chart_x0 + 12, chart_y0 + row_gap * 6 + 34), "AUROC change vs matched classical baseline", F["tiny_bold"], COLORS["axis"])

    for i, row in enumerate(deltas):
        y = chart_y0 + i * row_gap
        tissue = str(row["tissue"]).title().replace("Gastrocnemius", "Gastro.")
        sc = float(row["scgpt"])
        gf = float(row["geneformer"])
        sc_x = map_x(sc, chart_x0, chart_w, lo, hi)
        gf_x = map_x(gf, chart_x0, chart_w, lo, hi)
        text(draw, (x0 + 62, y - 18), tissue, F["small_bold"], COLORS["text"])
        draw.line((chart_x0, y + 4, chart_x0 + chart_w, y + 4), fill=rgba("#2A3546", 128), width=8)
        draw.line((gf_x, y + 4, sc_x, y + 4), fill=rgba("#98A7BA", 82), width=5)
        draw.ellipse((gf_x - 16, y - 12, gf_x + 16, y + 20), fill=rgba(COLORS["rose"], 220), outline=COLORS["white"], width=2)
        draw.ellipse((sc_x - 18, y - 14, sc_x + 18, y + 22), fill=rgba(COLORS["teal"], 230), outline=COLORS["white"], width=2)
        if sc > 0:
            label_x = min(sc_x + 28, x1 - 160)
            text(draw, (label_x, y - 18), f"+{fmt3(sc)}", F["micro_bold"], COLORS["green"])

    legend_y = y1 - 250
    rounded(draw, (x0 + 58, legend_y, x1 - 58, y1 - 58), 30, COLORS["panel2"], "#2A394D", 2)
    draw.ellipse((x0 + 92, legend_y + 42, x0 + 126, legend_y + 76), fill=COLORS["teal"])
    text(draw, (x0 + 144, legend_y + 43), "scGPT", F["tiny_bold"], COLORS["text"])
    draw.ellipse((x0 + 292, legend_y + 42, x0 + 326, legend_y + 76), fill=COLORS["rose"])
    text(draw, (x0 + 344, legend_y + 43), "Mouse-Geneformer", F["tiny_bold"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 92, legend_y + 96),
        "Local gains appear in liver and kidney for scGPT; the reference floor remains essential for tier comparisons.",
        F["body_bold"],
        x1 - x0 - 184,
        COLORS["text"],
        leading=6,
    )


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    y = 1872
    rounded(draw, (136, y, 3704, 2042), 32, "#101823", "#2A394D", 2)
    text(draw, (180, y + 38), "Slide 26 readout", F["small_bold"], COLORS["sky"])
    paragraph(
        draw,
        (180, y + 76),
        "Pretraining and scale improve selected tissues; matched classical floors and input-surface accounting keep the comparison calibrated.",
        F["body_bold"],
        2360,
        COLORS["text"],
        leading=6,
    )
    paragraph(
        draw,
        (2700, y + 46),
        "Next: explain why bulk RNA-seq is a hard adaptation surface for single-cell-pretrained models.",
        F["tiny_bold"],
        850,
        COLORS["muted"],
        leading=5,
    )
    text(
        draw,
        (140, 2102),
        "Takeaway: model tier matters less than task fit and the matched baseline floor on this benchmark surface.",
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
    draw_comparison_rules(draw, (136, 345, 1056, 1815), data)
    draw_shared_chart(draw, (1096, 345, 2466, 1815), data["shared_rows"])
    draw_delta_chart(draw, (2506, 345, 3704, 1815), data["deltas"])
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
                    str(RESULTS_SUMMARY.relative_to(ROOT)),
                    str(CANONICAL_RESULTS.relative_to(ROOT)),
                    str(GENEFORMER_SUMMARY.relative_to(ROOT)),
                    str(SCGPT_SUMMARY.relative_to(ROOT)),
                    "evaluation/submission_*_zeroshot_A*_eval.json",
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

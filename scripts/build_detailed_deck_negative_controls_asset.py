#!/usr/bin/env python3
"""Build the detailed-deck negative-controls proof asset."""

from __future__ import annotations

import json
import statistics
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
NC1_SUMMARY = ROOT / "evaluation" / "NC1_permutation_summary.json"
NC2_HK = ROOT / "evaluation" / "NC2_housekeeping_summary.json"
V4_EVAL = ROOT / "v4" / "evaluation"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "negative_controls"
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
    "title": load_font(78, True),
    "subtitle": load_font(36, False),
    "h2": load_font(42, True),
    "h3": load_font(33, True),
    "body": load_font(29, False),
    "body_bold": load_font(29, True),
    "small": load_font(25, False),
    "small_bold": load_font(25, True),
    "tiny": load_font(21, False),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18, False),
    "micro_bold": load_font(18, True),
    "stat": load_font(76, True),
    "huge": load_font(98, True),
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


def quantile(vals: list[float], q: float) -> float:
    vals = sorted(vals)
    if not vals:
        return 0.0
    pos = (len(vals) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(vals) - 1)
    frac = pos - lo
    return vals[lo] * (1 - frac) + vals[hi] * frac


def load_v4_controls() -> dict[str, dict[str, object]]:
    groups: dict[str, list[dict[str, object]]] = {"NC1_shuffled": [], "NC2_random": []}
    for path in V4_EVAL.glob("NC*.json"):
        record = json.loads(path.read_text())
        experiment = str(record["experiment"])
        if experiment in groups:
            groups[experiment].append(record)
    summary: dict[str, dict[str, object]] = {}
    for experiment, records in groups.items():
        vals = [float(r["mean_auroc"]) for r in records]
        ps = [float(r["mean_perm_pvalue"]) for r in records]
        summary[experiment] = {
            "records": records,
            "values": vals,
            "n": len(records),
            "mean": statistics.mean(vals),
            "median": statistics.median(vals),
            "min": min(vals),
            "q1": quantile(vals, 0.25),
            "q3": quantile(vals, 0.75),
            "max": max(vals),
            "mean_perm_p": statistics.mean(ps),
            "sig_count": sum(p < 0.05 for p in ps),
        }
    return summary


def load_data() -> dict[str, object]:
    control_summary = load_v4_controls()
    nc1 = json.loads(NC1_SUMMARY.read_text())
    hk = json.loads(NC2_HK.read_text())

    a_entries = {entry["task_id"]: entry for entry in nc1["entries"] if entry["category"] == "A"}
    examples = [
        {"label": "Thymus", "score": a_entries["A4"]["observed"], "p": a_entries["A4"]["perm_p"], "state": "passes"},
        {"label": "Gastrocnemius", "score": a_entries["A2"]["observed"], "p": a_entries["A2"]["perm_p"], "state": "passes"},
        {"label": "Liver", "score": a_entries["A1"]["observed"], "p": a_entries["A1"]["perm_p"], "state": "boundary"},
        {"label": "Kidney", "score": a_entries["A3"]["observed"], "p": a_entries["A3"]["perm_p"], "state": "boundary"},
    ]
    hk_results = hk["results"]
    hk_pass = sum(1 for row in hk_results if row["verdict"] == "PASS")
    hk_warn = sum(1 for row in hk_results if row["verdict"] == "WARN")
    return {
        "controls": control_summary,
        "permutation_summary": nc1["summary"],
        "examples": examples,
        "housekeeping": {
            "n_hk_mapped": hk["n_hk_mapped"],
            "n_hk_symbols": hk["n_hk_symbols"],
            "pass": hk_pass,
            "audit": hk_warn,
            "rows": hk_results,
        },
    }


def stat_badge(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, label: str, value: str, accent: str) -> None:
    rounded(draw, (x, y, x + w, y + 104), 26, "#121B28", "#2A394D", 2)
    text(draw, (x + 28, y + 24), label, F["tiny_bold"], accent)
    text(draw, (x + 28, y + 56), value, F["small_bold"], COLORS["text"])


def draw_header(draw: ImageDraw.ImageDraw) -> None:
    badges = [
        ("LABEL SHUFFLE", "64 runs", 335, COLORS["sky"]),
        ("RANDOM FEATURES", "64 runs", 355, COLORS["teal"]),
        ("PERM GATE", "28 rows", 270, COLORS["amber"]),
        ("HK AUDIT", "55 genes", 260, COLORS["green"]),
    ]
    bx = 2305
    for label, value, width, accent in badges:
        stat_badge(draw, bx, 72, width, label, value, accent)
        bx += width + 26


def draw_control_family_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "A. Control families", F["h2"], COLORS["text"])
    text(draw, (x0 + 50, y0 + 96), "Each control tests a different route to score inflation.", F["small"], COLORS["muted"])

    rows = [
        ("Label shuffle", "Keep matrix and split; scramble labels.", "Expected: chance", COLORS["sky"]),
        ("Random features", "Keep model and split; remove biology.", "Expected: chance", COLORS["teal"]),
        ("Housekeeping genes", "Use stable genes as a residual-structure audit.", "Expected: lower signal", COLORS["green"]),
        ("Permutation / CI gate", "Promote only rows that survive the null gate.", "Expected: selective pass", COLORS["amber"]),
    ]
    y = y0 + 185
    for idx, (title, body, expected, color) in enumerate(rows, start=1):
        rounded(draw, (x0 + 56, y, x1 - 56, y + 210), 26, "#131D2A", "#2A394D", 2)
        draw.ellipse((x0 + 88, y + 36, x0 + 152, y + 100), fill=rgba(color, 55), outline=color, width=4)
        text(draw, (x0 + 120, y + 54), str(idx), F["small_bold"], COLORS["text"], anchor="mm")
        text(draw, (x0 + 182, y + 34), title, F["h3"], COLORS["text"])
        paragraph(draw, (x0 + 182, y + 82), body, F["small"], x1 - x0 - 280, COLORS["muted"], 7)
        rounded(draw, (x0 + 182, y + 146, x0 + 520, y + 188), 18, rgba(color, 32), color, 1)
        text(draw, (x0 + 204, y + 157), expected, F["tiny_bold"], COLORS["text"])
        y += 246

    rule = (x0 + 56, y1 - 168, x1 - 56, y1 - 52)
    rounded(draw, rule, 24, "#14241F", "#315E51", 2)
    text(draw, (rule[0] + 32, rule[1] + 28), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (rule[0] + 32, rule[1] + 68),
        "Trust rises when null controls cluster near chance and residual controls route to audit.",
        F["small"],
        rule[2] - rule[0] - 70,
        COLORS["text"],
        7,
    )


def x_scale(v: float, x0: int, x1: int) -> int:
    return int(x0 + (v - 0.25) / (0.90 - 0.25) * (x1 - x0))


def draw_distribution_row(
    draw: ImageDraw.ImageDraw,
    x0: int,
    y: int,
    x1: int,
    label: str,
    summary: dict[str, object],
    color: str,
) -> None:
    vals = sorted(float(v) for v in summary["values"])
    axis_y = y + 92
    text(draw, (x0, y), label, F["h3"], COLORS["text"])
    text(draw, (x0, y + 42), f"{int(summary['n'])} runs | mean AUROC {float(summary['mean']):.3f}", F["small"], COLORS["muted"])
    text(draw, (x1 - 10, y + 14), f"{int(summary['sig_count'])}/{int(summary['n'])}", F["stat"], color, anchor="ra")
    text(draw, (x1 - 10, y + 88), "p<0.05", F["tiny"], COLORS["muted"], anchor="ra")

    draw.line((x0, axis_y, x1, axis_y), fill="#26354A", width=3)
    chance_x = x_scale(0.50, x0, x1)
    draw.line((chance_x, axis_y - 60, chance_x, axis_y + 68), fill=rgba(COLORS["amber"], 190), width=4)
    text(draw, (chance_x, axis_y + 78), "0.50 chance", F["tiny"], COLORS["amber"], anchor="ma")

    for i, value in enumerate(vals):
        px = x_scale(value, x0, x1)
        py = axis_y - 23 + (i % 7) * 8
        draw.ellipse((px - 4, py - 4, px + 4, py + 4), fill=rgba(color, 135))

    low = x_scale(float(summary["min"]), x0, x1)
    high = x_scale(float(summary["max"]), x0, x1)
    q1 = x_scale(float(summary["q1"]), x0, x1)
    q3 = x_scale(float(summary["q3"]), x0, x1)
    median = x_scale(float(summary["median"]), x0, x1)
    mean = x_scale(float(summary["mean"]), x0, x1)
    draw.line((low, axis_y - 18, high, axis_y - 18), fill=color, width=5)
    rounded(draw, (q1, axis_y - 44, q3, axis_y + 8), 16, rgba(color, 55), color, 3)
    draw.line((median, axis_y - 48, median, axis_y + 12), fill=COLORS["white"], width=5)
    draw.ellipse((mean - 13, axis_y + 24, mean + 13, axis_y + 50), fill=color, outline=COLORS["white"], width=2)
    text(draw, (mean, axis_y + 60), "mean", F["micro"], COLORS["muted"], anchor="ma")


def draw_distribution_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "B. Shuffled and random inputs return to chance", F["h2"], COLORS["text"])
    text(draw, (x0 + 50, y0 + 96), "Across v4 control sweeps, 0/64 pass the permutation gate.", F["small"], COLORS["muted"])

    axis_x0, axis_x1 = x0 + 90, x1 - 105
    for tick in [0.25, 0.50, 0.75, 0.90]:
        px = x_scale(tick, axis_x0, axis_x1)
        draw.line((px, y0 + 220, px, y0 + 880), fill="#243247", width=2)
        text(draw, (px, y0 + 900), f"{tick:.2f}", F["tiny"], COLORS["axis"], anchor="ma")
    text(draw, ((axis_x0 + axis_x1) // 2, y0 + 940), "Mean AUROC across control runs", F["small"], COLORS["muted"], anchor="ma")

    controls = data["controls"]
    draw_distribution_row(draw, axis_x0, y0 + 190, axis_x1, "Label-shuffled controls", controls["NC1_shuffled"], COLORS["sky"])
    draw_distribution_row(draw, axis_x0, y0 + 555, axis_x1, "Random-feature controls", controls["NC2_random"], COLORS["teal"])

    perm = data["permutation_summary"]
    gate_box = (x0 + 70, y1 - 395, x1 - 70, y1 - 70)
    rounded(draw, gate_box, 28, "#131D2A", "#2A394D", 2)
    text(draw, (gate_box[0] + 34, gate_box[1] + 28), "Permutation / CI gate", F["h3"], COLORS["text"])
    text(draw, (gate_box[0] + 34, gate_box[1] + 72), "The gate keeps result calls selective across benchmark rows.", F["small"], COLORS["muted"])
    total = int(perm["n_total_entries"])
    passed = int(perm["n_significant"])
    filtered = total - passed
    bx0, by0, bw, bh = gate_box[0] + 40, gate_box[1] + 142, gate_box[2] - gate_box[0] - 80, 66
    rounded(draw, (bx0, by0, bx0 + bw, by0 + bh), 30, "#0D141F", "#26354A", 2)
    passed_w = int(bw * passed / total)
    rounded(draw, (bx0, by0, bx0 + passed_w, by0 + bh), 30, COLORS["green"], None, 0)
    draw.rectangle((bx0 + passed_w - 8, by0, bx0 + passed_w + 8, by0 + bh), fill="#0D141F")
    text(draw, (bx0 + 24, by0 + 18), f"{passed} pass", F["small_bold"], "#081018")
    text(draw, (bx0 + passed_w + 24, by0 + 18), f"{filtered} held for audit", F["small_bold"], COLORS["muted"])
    text(draw, (gate_box[0] + 40, gate_box[1] + 238), f"{perm['n_permutation_tested']} permutation-tested rows + {perm['n_ci_based_only']} CI-only transfer rows", F["small"], COLORS["muted"])


def draw_readout_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict[str, object]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 46), "C. Readout rule", F["h2"], COLORS["text"])
    text(draw, (x0 + 50, y0 + 96), "Controls turn score interpretation into a visible gate.", F["small"], COLORS["muted"])

    examples = data["examples"]
    card = (x0 + 60, y0 + 170, x1 - 60, y0 + 540)
    rounded(draw, card, 28, "#14241F", "#315E51", 2)
    text(draw, (card[0] + 34, card[1] + 30), "Signal rows survive the gate", F["h3"], COLORS["green"])
    y = card[1] + 96
    for ex in [r for r in examples if r["state"] == "passes"]:
        rounded(draw, (card[0] + 34, y, card[2] - 34, y + 86), 20, "#101823", "#315E51", 1)
        text(draw, (card[0] + 60, y + 21), ex["label"], F["small_bold"], COLORS["text"])
        text(draw, (card[0] + 260, y + 21), f"AUROC {float(ex['score']):.3f}", F["small_bold"], COLORS["green"])
        text(draw, (card[0] + 505, y + 21), f"perm p={float(ex['p']):.3f}", F["small"], COLORS["muted"])
        y += 105

    card2 = (x0 + 60, y0 + 585, x1 - 60, y0 + 935)
    rounded(draw, card2, 28, "#211E17", "#69532B", 2)
    text(draw, (card2[0] + 34, card2[1] + 30), "Audit rows stay visible", F["h3"], COLORS["amber"])
    y = card2[1] + 96
    for ex in [r for r in examples if r["state"] == "boundary"]:
        rounded(draw, (card2[0] + 34, y, card2[2] - 34, y + 86), 20, "#101823", "#69532B", 1)
        text(draw, (card2[0] + 60, y + 21), ex["label"], F["small_bold"], COLORS["text"])
        text(draw, (card2[0] + 260, y + 21), f"AUROC {float(ex['score']):.3f}", F["small_bold"], COLORS["amber"])
        text(draw, (card2[0] + 505, y + 21), f"perm p={float(ex['p']):.3f}", F["small"], COLORS["muted"])
        y += 105

    hk = data["housekeeping"]
    card3 = (x0 + 60, y0 + 980, x1 - 60, y1 - 55)
    rounded(draw, card3, 28, "#131D2A", "#2A394D", 2)
    text(draw, (card3[0] + 34, card3[1] + 30), "Housekeeping panel routes residuals to audit", F["h3"], COLORS["sky"])
    text(draw, (card3[0] + 34, card3[1] + 86), f"{hk['n_hk_mapped']} mapped genes | {hk['pass']} pass | {hk['audit']} audit flags", F["small_bold"], COLORS["text"])
    paragraph(
        draw,
        (card3[0] + 34, card3[1] + 136),
        "Stable genes can expose batch-like structure. The readout stays strongest when full gene signal separates from this audit panel.",
        F["small"],
        card3[2] - card3[0] - 70,
        COLORS["muted"],
        7,
    )


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    text(draw, (205, 1925), "Takeaway", F["small_bold"], COLORS["sky"])
    source = (
        "Null controls should fall near chance, while housekeeping panels route residual structure into audit instead of interpretation."
    )
    paragraph(draw, (390, 1925), source, F["small"], 3180, COLORS["muted"], 6)
    text(draw, (205, 1995), "Readout", F["small_bold"], COLORS["teal"])
    readout = "Label-shuffled and random-feature controls cluster near chance; housekeeping features mark residual structure for audit."
    paragraph(draw, (390, 1995), readout, F["small"], 3180, COLORS["muted"], 6)


def main() -> None:
    data = load_data()

    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)

    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 48), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 38), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "CORE BENCHMARK CONTROL", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "Negative controls anchor the readout", F["title"])
    text(
        draw,
        (150, 214),
        "Shuffled labels and random features return toward chance, while housekeeping genes route residual structure to audit.",
        F["subtitle"],
        COLORS["muted"],
    )
    draw_header(draw)

    draw_control_family_panel(draw, (150, 350, 1110, 1800))
    draw_distribution_panel(draw, (1170, 350, 2600, 1800), data)
    draw_readout_panel(draw, (2660, 350, 3690, 1800), data)
    draw_footer(draw)

    png = OUT_DIR / "negative_controls_anchor_readout_premium.png"
    canvas.save(png, quality=95)

    gray = OUT_DIR / "negative_controls_anchor_readout_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    controls = data["controls"]
    manifest = {
        "title": "Negative controls anchor the readout",
        "sources": [
            str(NC1_SUMMARY.relative_to(ROOT)),
            str(NC2_HK.relative_to(ROOT)),
            "v4/evaluation/NC1_shuffled_*.json",
            "v4/evaluation/NC2_random_*.json",
        ],
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "control_summary": {
            "label_shuffled": {
                "n": controls["NC1_shuffled"]["n"],
                "mean_auroc": round(float(controls["NC1_shuffled"]["mean"]), 3),
                "median_auroc": round(float(controls["NC1_shuffled"]["median"]), 3),
                "range": [round(float(controls["NC1_shuffled"]["min"]), 3), round(float(controls["NC1_shuffled"]["max"]), 3)],
                "mean_perm_p": round(float(controls["NC1_shuffled"]["mean_perm_p"]), 3),
                "sig_count_p_lt_0_05": controls["NC1_shuffled"]["sig_count"],
            },
            "random_features": {
                "n": controls["NC2_random"]["n"],
                "mean_auroc": round(float(controls["NC2_random"]["mean"]), 3),
                "median_auroc": round(float(controls["NC2_random"]["median"]), 3),
                "range": [round(float(controls["NC2_random"]["min"]), 3), round(float(controls["NC2_random"]["max"]), 3)],
                "mean_perm_p": round(float(controls["NC2_random"]["mean_perm_p"]), 3),
                "sig_count_p_lt_0_05": controls["NC2_random"]["sig_count"],
            },
        },
        "permutation_gate_summary": data["permutation_summary"],
        "housekeeping_summary": data["housekeeping"],
    }
    manifest_path = OUT_DIR / "negative_controls_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

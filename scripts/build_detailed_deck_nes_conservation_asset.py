#!/usr/bin/env python3
"""Build the detailed-deck NES conservation proof asset.

The output is a high-resolution 16:9 PNG intended for the detailed
SpaceBio-Bench deck. It reads the canonical NES conservation table from
evaluation/NES_conservation_vs_transfer.json and draws a consulting-style proof
slide with:

- all tissue points shown;
- the primary 5-tissue rank set with gastrocnemius shown as a coverage note;
- a compact 4-tissue core check focused on thymus, eye, liver, and kidney;
- explicit scope language.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "evaluation" / "NES_conservation_vs_transfer.json"
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "nes_conservation"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160

COLORS = {
    "bg": "#0C111A",
    "panel": "#121A26",
    "panel2": "#172231",
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
    "gray": "#8C98A8",
    "dark_gray": "#536173",
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
    "rank": load_font(27, True),
}


def rank(vals: list[float]) -> list[float]:
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + j + 2) / 2
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def pearson(xs: list[float], ys: list[float]) -> float:
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy)


def spearman(xs: list[float], ys: list[float]) -> float:
    return pearson(rank(xs), rank(ys))


def rounded(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], radius: int, fill: str, outline: str | None = None, width: int = 1) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], value: str, font: ImageFont.ImageFont, fill: str = COLORS["text"], anchor: str | None = None) -> None:
    draw.text(xy, value, font=font, fill=fill, anchor=anchor)


def multiline(draw: ImageDraw.ImageDraw, xy: tuple[int, int], lines: Iterable[str], font: ImageFont.ImageFont, fill: str = COLORS["muted"], leading: int = 10) -> None:
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


def linear_fit(xs: list[float], ys: list[float]) -> tuple[float, float]:
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    denom = sum((x - mx) ** 2 for x in xs)
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / denom
    intercept = my - slope * mx
    return slope, intercept


def main() -> None:
    source = json.loads(SOURCE.read_text())
    raw = source["data"]
    order = ["thymus", "eye", "skin", "liver", "kidney", "gastrocnemius"]
    data = []
    for tissue in order:
        entry = raw[tissue]
        if tissue == "skin":
            note = "NES coverage uses 3 fGSEA sources: RR-6, MHU-2 dorsal, and MHU-2 femoral. Transfer AUROC uses merged MHU-2 plus RR-7."
        elif tissue == "gastrocnemius":
            note = "Gastrocnemius has 2 fGSEA missions in the NES coverage table."
        else:
            note = ""
        data.append(
            {
                "tissue": tissue,
                "label": "Gastrocnemius" if tissue == "gastrocnemius" else tissue.capitalize(),
                "x": float(entry["nes_mean_r"]),
                "y": float(entry["transfer_auroc"]),
                "n_nes": int(entry["n_missions_nes"]),
                "note": note,
            }
        )

    primary = [d for d in data if d["tissue"] != "gastrocnemius"]
    strict = [d for d in data if d["tissue"] not in {"gastrocnemius", "skin"}]
    primary_r = spearman([d["x"] for d in primary], [d["y"] for d in primary])
    strict_r = spearman([d["x"] for d in strict], [d["y"] for d in strict])
    all_r = spearman([d["x"] for d in data], [d["y"] for d in data])
    slope, intercept = linear_fit([d["x"] for d in primary], [d["y"] for d in primary])

    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)

    # Background grid and subtle horizon lines.
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 58), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 46), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "CORE BENCHMARK MECHANISM", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "Conserved pathways predict transfer", F["title"])
    text(
        draw,
        (150, 216),
        "Tissues with more consistent pathway activity across missions tend to support better mission-held-out transfer.",
        F["subtitle"],
        COLORS["muted"],
    )

    # Top badges.
    badges = [
        ("ALL POINTS SHOWN", "6 tissues", 360),
        ("PRIMARY RANK SET", "5 tissues", 360),
        ("READOUT", "transfer association", 450),
    ]
    bx = 2490
    for kicker, body, badge_w in badges:
        rounded(draw, (bx, 72, bx + badge_w, 176), 26, "#121B28", "#2A394D", 2)
        text(draw, (bx + 28, 96), kicker, F["tiny"], COLORS["sky"])
        text(draw, (bx + 28, 126), body, F["small_bold"], COLORS["text"])
        bx += badge_w + 30

    # Main panels.
    plot_box = (150, 350, 2290, 1800)
    side_box = (2380, 350, 3690, 1800)
    rounded(draw, plot_box, 34, "#101823", "#29374A", 2)
    rounded(draw, side_box, 34, "#101823", "#29374A", 2)

    # Scatter plot area.
    px0, py0, px1, py1 = 310, 505, 2150, 1640
    draw.rectangle((px0, py0, px1, py1), fill="#0E1520")
    xmin, xmax = 0.0, 0.65
    ymin, ymax = 0.50, 0.90

    def sx(v: float) -> int:
        return int(px0 + (v - xmin) / (xmax - xmin) * (px1 - px0))

    def sy(v: float) -> int:
        return int(py1 - (v - ymin) / (ymax - ymin) * (py1 - py0))

    # Axes and grid.
    for val in [0.0, 0.15, 0.30, 0.45, 0.60]:
        x = sx(val)
        draw.line((x, py0, x, py1), fill="#253246", width=2)
        text(draw, (x, py1 + 28), f"{val:.2f}", F["small"], COLORS["axis"], anchor="ma")
    for val in [0.50, 0.60, 0.70, 0.80, 0.90]:
        y = sy(val)
        draw.line((px0, y, px1, y), fill="#253246", width=2)
        text(draw, (px0 - 32, y), f"{val:.2f}", F["small"], COLORS["axis"], anchor="rm")

    draw.line((px0, py1, px1, py1), fill=COLORS["axis"], width=3)
    draw.line((px0, py0, px0, py1), fill=COLORS["axis"], width=3)
    text(draw, ((px0 + px1) // 2, 1740), "Pathway activity conservation across missions (mean NES Spearman r)", F["body"], COLORS["muted"], anchor="ma")
    # Rotated y-axis label.
    y_label = Image.new("RGBA", (700, 70), (0, 0, 0, 0))
    yd = ImageDraw.Draw(y_label)
    yd.text((350, 10), "Cross-mission transfer AUROC", font=F["body"], fill=COLORS["muted"], anchor="ma")
    y_label = y_label.rotate(90, expand=True)
    canvas.paste(y_label, (165, 820), y_label)

    text(draw, (310, 415), "A. Tissue-level relationship", F["h2"], COLORS["text"])
    text(draw, (310, 465), "Each dot is a tissue; transfer is Category B PCA-LR AUROC.", F["small"], COLORS["muted"])

    # Trend line for primary 5 tissues.
    x_start, x_end = 0.03, 0.63
    y_start, y_end = intercept + slope * x_start, intercept + slope * x_end
    draw.line((sx(x_start), sy(y_start), sx(x_end), sy(y_end)), fill=COLORS["teal"], width=7)
    draw.line((sx(x_start), sy(y_start) + 8, sx(x_end), sy(y_end) + 8), fill=rgba(COLORS["teal"], 80), width=12)

    tissue_color = {
        "thymus": COLORS["teal"],
        "eye": COLORS["sky"],
        "skin": COLORS["amber"],
        "liver": COLORS["gray"],
        "kidney": COLORS["gray"],
        "gastrocnemius": COLORS["rose"],
    }
    offsets = {
        "thymus": (-155, -80),
        "eye": (-95, -78),
        "skin": (70, -60),
        "liver": (70, -88),
        "kidney": (92, 44),
        "gastrocnemius": (70, -82),
    }
    for item in data:
        x, y = sx(item["x"]), sy(item["y"])
        color = tissue_color[item["tissue"]]
        if item["tissue"] == "gastrocnemius":
            draw.ellipse((x - 24, y - 24, x + 24, y + 24), outline=color, width=7)
            draw.ellipse((x - 10, y - 10, x + 10, y + 10), fill=color)
        else:
            draw.ellipse((x - 26, y - 26, x + 26, y + 26), fill=color, outline="#FFFFFF", width=4)
        dx, dy = offsets[item["tissue"]]
        lx, ly = x + dx, y + dy
        label = f"{item['label']}\nNES {item['x']:.3f} | AUROC {item['y']:.3f}"
        label_w = 340
        line_end_x = lx - 16 if dx > 0 else lx + label_w
        draw.line((x, y, line_end_x, ly + 30), fill=rgba(color, 130), width=2)
        rounded(draw, (lx - 16, ly - 12, lx + label_w, ly + 84), 14, "#172231", "#2C3C52", 1)
        draw.text((lx, ly), label, font=F["tiny"], fill=COLORS["text"])

    # Side panel.
    text(draw, (2460, 425), "B. Read the rank signal", F["h2"])
    text(draw, (2460, 475), "The main relationship is rank-level and tissue-level.", F["small"], COLORS["muted"])

    rounded(draw, (2460, 550, 3075, 820), 28, "#162232", "#36516A", 2)
    text(draw, (2495, 590), "PRIMARY RANK SET", F["small_bold"], COLORS["teal"])
    text(draw, (2495, 620), f"rho = {primary_r:.2f}", F["number"], COLORS["text"])
    multiline(
        draw,
        (2498, 726),
        ["5 tissues", "gastrocnemius shown", "2-mission NES coverage"],
        F["small"],
        COLORS["muted"],
        6,
    )

    rounded(draw, (3105, 550, 3665, 820), 28, "#162232", "#36516A", 2)
    text(draw, (3140, 590), "CORE CHECK", F["small_bold"], COLORS["sky"])
    text(draw, (3140, 620), f"rho = {strict_r:.2f}", F["number"], COLORS["text"])
    multiline(
        draw,
        (3143, 726),
        ["4 tissues", "thymus / eye", "liver / kidney"],
        F["small"],
        COLORS["muted"],
        6,
    )

    # Rank bridge.
    bridge_box = (2460, 890, 3610, 1325)
    rounded(draw, bridge_box, 28, "#131D2A", "#2A394D", 2)
    text(draw, (2495, 925), "Rank bridge", F["h3"], COLORS["text"])
    text(draw, (2495, 968), "Conservation rank mostly tracks transfer rank.", F["small"], COLORS["muted"])
    left_x, right_x = 2690, 3385
    top_y, gap = 1045, 58
    primary_sorted_x = sorted(primary, key=lambda d: d["x"], reverse=True)
    primary_sorted_y = sorted(primary, key=lambda d: d["y"], reverse=True)
    y_by_left = {d["tissue"]: top_y + i * gap for i, d in enumerate(primary_sorted_x)}
    y_by_right = {d["tissue"]: top_y + i * gap for i, d in enumerate(primary_sorted_y)}
    text(draw, (left_x, 1000), "Pathway conservation", F["tiny"], COLORS["teal"], anchor="ma")
    text(draw, (right_x, 1000), "Transfer AUROC", F["tiny"], COLORS["sky"], anchor="ma")
    for i, item in enumerate(primary_sorted_x):
        y = y_by_left[item["tissue"]]
        text(draw, (left_x - 185, y), f"{i+1}", F["rank"], COLORS["dim"], anchor="mm")
        text(draw, (left_x - 150, y - 16), item["label"], F["small_bold"], COLORS["text"])
        text(draw, (left_x - 150, y + 12), f"{item['x']:.3f}", F["tiny"], COLORS["muted"])
    for i, item in enumerate(primary_sorted_y):
        y = y_by_right[item["tissue"]]
        text(draw, (right_x + 185, y), f"{i+1}", F["rank"], COLORS["dim"], anchor="mm")
        text(draw, (right_x + 18, y - 16), item["label"], F["small_bold"], COLORS["text"])
        text(draw, (right_x + 18, y + 12), f"{item['y']:.3f}", F["tiny"], COLORS["muted"])
    for item in primary:
        color = tissue_color[item["tissue"]]
        draw.line((left_x + 120, y_by_left[item["tissue"]], right_x - 120, y_by_right[item["tissue"]]), fill=rgba(color, 210), width=5)
        draw.ellipse((left_x + 108, y_by_left[item["tissue"]] - 10, left_x + 128, y_by_left[item["tissue"]] + 10), fill=color)
        draw.ellipse((right_x - 128, y_by_right[item["tissue"]] - 10, right_x - 108, y_by_right[item["tissue"]] + 10), fill=color)

    # Readout card.
    rounded(draw, (2460, 1370, 3610, 1640), 28, "#211E17", "#69532B", 2)
    text(draw, (2495, 1410), "Readout frame", F["h3"], COLORS["amber"])
    scope_lines = [
        "This supports a biology hypothesis:",
        "pathway-level consistency helps explain gene-level transfer.",
        "Mechanistic interpretation belongs to follow-up biology slides.",
    ]
    multiline(draw, (2495, 1462), scope_lines, F["body"], COLORS["text"], 8)
    text(draw, (2495, 1592), f"All-point Spearman across six tissues: rho = {all_r:.2f}", F["small"], COLORS["muted"])

    # Footer.
    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    footer_label_x = 205
    footer_text_x = 390
    footer_text_width = 3220
    text(draw, (footer_label_x, 1925), "Takeaway", F["small_bold"], COLORS["sky"])
    footer = (
        "Pathway conservation and transfer move together in the strongest tissues. "
        "NES is fGSEA normalized enrichment score; conservation is mean pairwise Spearman r across missions."
    )
    draw.text((footer_text_x, 1925), fit_label(draw, footer, F["small"], footer_text_width), font=F["small"], fill=COLORS["muted"])
    text(draw, (footer_label_x, 1995), "Readout", F["small_bold"], COLORS["teal"])
    readout = "High conservation / high transfer: thymus. Near-zero conservation / low transfer: liver and kidney. Gastrocnemius contributes a 2-mission NES point."
    draw.text((footer_text_x, 1995), fit_label(draw, readout, F["small"], footer_text_width), font=F["small"], fill=COLORS["muted"])

    png = OUT_DIR / "nes_conservation_predicts_transfer_premium.png"
    canvas.save(png, quality=95)

    # QA grayscale and a simple manifest.
    gray = OUT_DIR / "nes_conservation_predicts_transfer_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    manifest = {
        "title": "Conserved pathways predict transfer",
        "source": str(SOURCE.relative_to(ROOT)),
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "data": data,
        "statistics": {
            "spearman_all_6": round(all_r, 3),
            "spearman_primary_5_rank_set": round(primary_r, 3),
            "spearman_core_4_rank_set": round(strict_r, 3),
        },
        "scope": (
            "Association supports the benchmark-mechanism hypothesis that pathway-level "
            "response consistency helps explain gene-level transfer. Mechanistic interpretation belongs to follow-up biology slides."
        ),
        "gastrocnemius_note": "Gastrocnemius has 2 fGSEA missions in the NES coverage table.",
    }
    manifest_path = OUT_DIR / "nes_conservation_predicts_transfer_premium_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

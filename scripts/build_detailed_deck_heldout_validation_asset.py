#!/usr/bin/env python3
"""Build the detailed-deck held-out validation proof asset.

The output is a high-resolution 16:9 PNG for the detailed SpaceBio-Bench deck.
It shows two independent held-out checks:

- RR-23 thymus: train on the original thymus mission set, test on a reserved
  RR-23 mission.
- RR-7 skin: train on MHU-2 + RR-6, test on RR-7.

The slide is designed to make the split visible before the score, and to keep
sample sizes, confidence intervals, and p-values visible.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
    / "heldout_validation"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "panel2": "#151F2D",
    "grid": "#2A3546",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#667487",
    "teal": "#5FD3C4",
    "sky": "#73A7FF",
    "green": "#84D278",
    "amber": "#F4C26B",
    "rose": "#E17882",
    "thymus": "#E8A34A",
    "skin": "#5EC69C",
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
    "number": load_font(118, True),
    "mission": load_font(42, True),
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
    a4 = json.loads((ROOT / "evaluation" / "A4_holdout_results.json").read_text())
    a5 = json.loads((ROOT / "evaluation" / "A5_holdout_results.json").read_text())
    a4_lomo = json.loads((ROOT / "evaluation" / "A4_baseline_results.json").read_text())
    a5_lomo = json.loads((ROOT / "evaluation" / "A5_baseline_results.json").read_text())
    return {
        "thymus": {
            "tissue": "Thymus",
            "mission": "RR-23",
            "duration": "~30 days",
            "source": "OSD-515",
            "color": COLORS["thymus"],
            "headline": a4["lr"],
            "models": [
                ("LR ElasticNet", a4["lr"]),
                ("Random Forest", a4["rf"]),
                ("PCA-50 + LR", a4["pca_lr"]),
            ],
            "lomo_label": "LOMO context",
            "lomo_value": a4_lomo["pca_lr"]["mean_auroc"],
            "lomo_model": "A4 PCA-LR mean",
            "train": ["MHU-1", "MHU-2", "RR-6", "RR-9"],
            "test": ["RR-23"],
            "interpretation": "Reserved held-out mission confirms thymus generalization beyond the LOMO loop.",
        },
        "skin": {
            "tissue": "Skin",
            "mission": "RR-7",
            "duration": "~75 days",
            "source": "OSD-254",
            "color": COLORS["skin"],
            "headline": a5["lr"],
            "models": [
                ("LR ElasticNet", a5["lr"]),
                ("PCA-50 + LR", a5["pca_lr"]),
                ("Random Forest", a5["rf"]),
            ],
            "lomo_label": "LOMO context",
            "lomo_value": a5_lomo["lr"]["mean_auroc"],
            "lomo_model": "A5 LR mean",
            "train": ["MHU-2", "RR-6"],
            "test": ["RR-7"],
            "interpretation": "Second tissue and longest skin mission preserve high-signal generalization.",
        },
    }


def fmt_p(p: float) -> str:
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def draw_mission_tokens(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    tokens: list[str],
    color: str,
    heldout: bool = False,
) -> int:
    cur_x = x
    for token in tokens:
        w = int(draw.textlength(token, font=F["small_bold"])) + 42
        fill = "#1A2635" if not heldout else "#2B1820"
        outline = "#33465D" if not heldout else COLORS["rose"]
        rounded(draw, (cur_x, y, cur_x + w, y + 52), 18, fill, outline, 2)
        dot = color if not heldout else COLORS["rose"]
        draw.ellipse((cur_x + 18, y + 19, cur_x + 30, y + 31), fill=dot)
        text(draw, (cur_x + 38, y + 15), token, F["small_bold"], COLORS["text"])
        cur_x += w + 16
    return cur_x


def draw_score_axis(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], record: dict, color: str, lomo_value: float, lomo_label: str) -> None:
    x0, y0, x1, y1 = box
    axis_y = y0 + 95
    xmin, xmax = 0.50, 1.00

    def sx(v: float) -> int:
        return int(x0 + (v - xmin) / (xmax - xmin) * (x1 - x0))

    draw.line((x0, axis_y, x1, axis_y), fill="#536173", width=4)
    for v in [0.5, 0.7, 0.9, 1.0]:
        x = sx(v)
        draw.line((x, axis_y - 14, x, axis_y + 14), fill="#8490A1", width=3)
        text(draw, (x, axis_y + 30), f"{v:.1f}", F["small"], COLORS["muted"], anchor="ma")
    text(draw, (sx(0.5), axis_y - 44), "chance", F["tiny"], COLORS["dim"], anchor="ma")

    ci_low = record["ci_lower"]
    ci_high = record["ci_upper"]
    auroc = record["auroc"]
    draw.line((sx(ci_low), axis_y, sx(ci_high), axis_y), fill=color, width=12)
    draw.line((sx(ci_low), axis_y - 28, sx(ci_low), axis_y + 28), fill=color, width=5)
    draw.line((sx(ci_high), axis_y - 28, sx(ci_high), axis_y + 28), fill=color, width=5)
    # Held-out diamond
    x = sx(auroc)
    diamond = [(x, axis_y - 42), (x + 42, axis_y), (x, axis_y + 42), (x - 42, axis_y)]
    draw.polygon(diamond, fill=color, outline=COLORS["white"])
    # LOMO context marker
    lx = sx(lomo_value)
    draw.ellipse((lx - 16, axis_y + 74, lx + 16, axis_y + 106), fill="#728094", outline="#D8DEE8", width=3)
    text(draw, (lx, axis_y + 124), lomo_label, F["tiny"], COLORS["muted"], anchor="ma")


def draw_model_ladder(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, models: list[tuple[str, dict]], color: str) -> None:
    rounded(draw, (x, y, x + w, y + 230), 24, "#121C29", "#2A394D", 2)
    text(draw, (x + 28, y + 24), "Model check", F["h3"])
    bar_x = x + 250
    bar_w = w - 310
    for i, (label, rec) in enumerate(models):
        row_y = y + 72 + i * 48
        text(draw, (x + 28, row_y + 6), label, F["small"], COLORS["muted"])
        draw.line((bar_x, row_y + 18, bar_x + bar_w, row_y + 18), fill="#2A3546", width=10)
        fill_w = int(bar_w * max(0, min(1, (rec["auroc"] - 0.5) / 0.5)))
        draw.line((bar_x, row_y + 18, bar_x + fill_w, row_y + 18), fill=color if i == 0 else "#728094", width=10)
        text(draw, (bar_x + bar_w + 18, row_y + 5), f"{rec['auroc']:.3f}", F["small_bold"], COLORS["text"])


def draw_card(draw: ImageDraw.ImageDraw, x: int, y: int, w: int, h: int, item: dict) -> None:
    color = item["color"]
    rounded(draw, (x, y, x + w, y + h), 34, COLORS["panel"], "#2A394D", 2)
    text(draw, (x + 54, y + 54), item["tissue"].upper(), F["kicker"], color)
    text(draw, (x + 54, y + 106), f"{item['mission']} held-out", F["mission"], COLORS["text"])
    text(draw, (x + 54, y + 156), f"{item['source']} | {item['duration']}", F["small"], COLORS["muted"])

    # Score block
    score_x = x + w - 460
    text(draw, (score_x, y + 56), "AUROC", F["small_bold"], COLORS["muted"])
    text(draw, (score_x, y + 82), f"{item['headline']['auroc']:.3f}", F["number"], COLORS["text"])
    ci = f"95% CI [{item['headline']['ci_lower']:.3f}, {item['headline']['ci_upper']:.3f}]"
    text(draw, (score_x, y + 220), ci, F["small"], COLORS["muted"])
    text(draw, (score_x, y + 252), f"perm p={fmt_p(item['headline']['perm_pvalue'])}", F["small_bold"], color)

    # Split map
    split_y = y + 250
    text(draw, (x + 54, split_y), "Train missions", F["small_bold"], COLORS["sky"])
    tx = draw_mission_tokens(draw, x + 54, split_y + 36, item["train"], color, False)
    text(draw, (tx + 30, split_y + 50), "->", F["h3"], COLORS["dim"])
    text(draw, (tx + 90, split_y), "Hidden test", F["small_bold"], COLORS["rose"])
    draw_mission_tokens(draw, tx + 90, split_y + 36, item["test"], color, True)

    # Sample and feature chips.
    chips_y = y + 390
    chips = [
        (f"train n={item['headline']['n_train']}", COLORS["sky"]),
        (f"test n={item['headline']['n_test']}", color),
        (f"{item['headline']['n_flight_test']} FLT / {item['headline']['n_ground_test']} GC", COLORS["amber"]),
        (f"{item['headline']['n_genes']:,} genes", COLORS["muted"]),
    ]
    cur_x = x + 54
    for label, c in chips:
        ww = int(draw.textlength(label, font=F["small_bold"])) + 36
        rounded(draw, (cur_x, chips_y, cur_x + ww, chips_y + 48), 16, "#151F2D", "#2A394D", 2)
        text(draw, (cur_x + 18, chips_y + 13), label, F["small_bold"], c)
        cur_x += ww + 14

    # Axis and model ladder.
    draw_score_axis(
        draw,
        (x + 100, y + 520, x + w - 100, y + 760),
        item["headline"],
        color,
        item["lomo_value"],
        item["lomo_label"],
    )
    text(draw, (x + 100, y + 772), f"Held-out diamond; gray dot = {item['lomo_model']} ({item['lomo_value']:.3f})", F["small"], COLORS["muted"])
    draw_model_ladder(draw, x + 54, y + 870, w - 108, item["models"], color)

    # Claim strip
    rounded(draw, (x + 54, y + h - 138, x + w - 54, y + h - 48), 20, "#171D25", "#2A394D", 2)
    text(draw, (x + 82, y + h - 108), "Interpretation", F["small_bold"], color)
    draw.text((x + 285, y + h - 110), fit_label(draw, item["interpretation"], F["small_bold"], w - 470), font=F["small_bold"], fill=COLORS["muted"])


def main() -> None:
    data = load_data()
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)

    # Background grid.
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 52), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 42), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "GENERALIZATION CHECK", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "Held-out missions confirm the signal", F["title"])
    text(
        draw,
        (150, 216),
        "Two independent test missions support that high-signal tissues generalize beyond the primary LOMO loop.",
        F["subtitle"],
        COLORS["muted"],
    )

    # Header badges.
    badges = [
        ("TEST UNIT", "entire mission"),
        ("TISSUES", "thymus + skin"),
        ("READOUT", "AUROC > 0.88"),
    ]
    bx = 2500
    for kicker, body in badges:
        rounded(draw, (bx, 72, bx + 360, 176), 26, "#121B28", "#2A394D", 2)
        text(draw, (bx + 28, 96), kicker, F["tiny"], COLORS["sky"])
        text(draw, (bx + 28, 126), body, F["small_bold"], COLORS["text"])
        bx += 390

    draw_card(draw, 150, 350, 1710, 1460, data["thymus"])
    draw_card(draw, 1980, 350, 1710, 1460, data["skin"])

    # Bottom scope note.
    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    footer_label_x = 205
    footer_text_x = 390
    footer_text_width = 3220
    text(draw, (footer_label_x, 1925), "Takeaway", F["small_bold"], COLORS["sky"])
    footer = "Reserved thymus and skin missions preserve the high-transfer pattern outside the main mission-held-out runs."
    draw.text((footer_text_x, 1925), fit_label(draw, footer, F["small"], footer_text_width), font=F["small"], fill=COLORS["muted"])
    text(draw, (footer_label_x, 1995), "Next", F["small_bold"], COLORS["amber"])
    scope = "Negative controls then test whether the readout drops toward chance where it should."
    draw.text((footer_text_x, 1995), fit_label(draw, scope, F["small"], footer_text_width), font=F["small"], fill=COLORS["muted"])

    png = OUT_DIR / "heldout_missions_confirm_signal_premium.png"
    canvas.save(png, quality=95)
    gray = OUT_DIR / "heldout_missions_confirm_signal_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    manifest = {
        "title": "Held-out missions confirm the signal",
        "sources": [
            "evaluation/A4_holdout_results.json",
            "evaluation/A5_holdout_results.json",
            "evaluation/A4_baseline_results.json",
            "evaluation/A5_baseline_results.json",
        ],
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "data": {
            key: {
                "tissue": item["tissue"],
                "mission": item["mission"],
                "duration": item["duration"],
                "source": item["source"],
                "headline_lr": {
                    "auroc": round(item["headline"]["auroc"], 3),
                    "ci_lower": round(item["headline"]["ci_lower"], 3),
                    "ci_upper": round(item["headline"]["ci_upper"], 3),
                    "perm_pvalue": item["headline"]["perm_pvalue"],
                    "n_train": item["headline"]["n_train"],
                    "n_test": item["headline"]["n_test"],
                    "n_flight_test": item["headline"]["n_flight_test"],
                    "n_ground_test": item["headline"]["n_ground_test"],
                    "n_genes": item["headline"]["n_genes"],
                },
                "lomo_context_value": round(item["lomo_value"], 3),
                "lomo_context_model": item["lomo_model"],
            }
            for key, item in data.items()
        },
        "scope": "Independent held-out mission checks support reactive-tissue transfer; companion analyses cover other tissues and modalities.",
    }
    manifest_path = OUT_DIR / "heldout_missions_confirm_signal_premium_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

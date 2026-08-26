#!/usr/bin/env python3
"""Build slide 49 asset: v8 turns benchmark evidence into follow-up workstreams."""

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
OUT_DIR = ASSET_ROOT / "v8_summary"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "v8_turns_benchmark_evidence_into_followup_workstreams_premium.png"
GRAY_PATH = OUT_DIR / "v8_turns_benchmark_evidence_into_followup_workstreams_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "v8_summary_manifest.json"
QA_NOTE = OUT_DIR / "V8_SUMMARY_ASSET_VISUAL_QA.md"

WORKSTREAMS = [
    {
        "key": "BRIDGE",
        "title": "Mouse pathways improve human prediction",
        "desc": "Pathway NES vectors add cross-species signal to supervised prediction.",
        "metrics": ["RF AUROC +0.175", "6 mouse NES vectors"],
        "color": "#66A6E8",
        "asset": ASSET_ROOT / "bridge_human_prediction" / "mouse_pathways_improve_human_prediction_premium.png",
        "manifest": ASSET_ROOT / "bridge_human_prediction" / "bridge_human_prediction_manifest.json",
    },
    {
        "key": "DECOMPOSE",
        "title": "Stressor interactions shape combined effects",
        "desc": "Factorial analog surfaces separate main axes from interaction-heavy responses.",
        "metrics": ["40-62% interaction", "ICP stability readout"],
        "color": "#E17882",
        "asset": ASSET_ROOT / "decompose_interactions" / "stressor_interactions_dominate_combined_effects_premium.png",
        "manifest": ASSET_ROOT / "decompose_interactions" / "decompose_interactions_manifest.json",
    },
    {
        "key": "INTERVENE",
        "title": "Perturbation hits prioritize axes",
        "desc": "Chemical reversal and CRISPR support turn tissue signatures into follow-up axes.",
        "metrics": ["215 unique hits", "2 CRISPR-supported pairs"],
        "color": "#8BD17C",
        "asset": ASSET_ROOT / "intervene_prioritization" / "perturbation_hits_prioritize_followup_axes_premium.png",
        "manifest": ASSET_ROOT / "intervene_prioritization" / "intervene_prioritization_manifest.json",
    },
    {
        "key": "CAUSAL",
        "title": "Stressor axes rank tissue readouts",
        "desc": "ICP stability becomes a ranked map from stressor axes to tissue biology.",
        "metrics": ["24 DAG edges", "9 environments"],
        "color": "#5FD3C4",
        "asset": ASSET_ROOT / "causal_map" / "causal_maps_connect_stressor_axes_to_tissue_readouts_premium.png",
        "manifest": ASSET_ROOT / "causal_map" / "causal_map_manifest.json",
    },
]

COLORS = {
    "bg": "#0B1119",
    "bg2": "#111721",
    "header": "#101826",
    "footer": "#090E15",
    "panel": "#111B28",
    "panel2": "#172233",
    "panel3": "#0F1825",
    "grid": "#263245",
    "text": "#F4F7FB",
    "muted": "#AAB6C6",
    "dim": "#687789",
    "blue": "#66A6E8",
    "teal": "#5FD3C4",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "violet": "#B39DFF",
    "rose": "#E17882",
    "ink": "#081018",
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
    "title": load_font(82, True),
    "subtitle": load_font(37),
    "section": load_font(40, True),
    "h2": load_font(34, True),
    "h3": load_font(30, True),
    "body": load_font(26),
    "body_bold": load_font(26, True),
    "small": load_font(23),
    "small_bold": load_font(23, True),
    "tiny": load_font(20),
    "tiny_bold": load_font(20, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "metric": load_font(48, True),
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


def wrap_lines(draw: ImageDraw.ImageDraw, value: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
    words = value.split()
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
    leading: int = 7,
) -> int:
    x, y = xy
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(224, int(draw.textlength(value, font=F["tiny_bold"]) + 78))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def arrow(draw: ImageDraw.ImageDraw, x1: int, y1: int, x2: int, y2: int, color: str, width: int = 4) -> None:
    draw.line((x1, y1, x2, y2), fill=color, width=width)
    draw.polygon([(x2, y2), (x2 - 20, y2 - 11), (x2 - 20, y2 + 11)], fill=color)


def background() -> Image.Image:
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image)
    for y in range(H):
        t = y / H
        draw.line((0, y, W, y), fill=blend(COLORS["bg"], COLORS["bg2"], t))
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 24), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 20), width=1)
    draw.rectangle((0, 0, W, 260), fill=COLORS["header"])
    draw.rectangle((0, 1900, W, H), fill=COLORS["footer"])
    return image


def fit_cover(path: Path, size: tuple[int, int]) -> Image.Image:
    source = Image.open(path).convert("RGB")
    target_w, target_h = size
    scale = max(target_w / source.width, target_h / source.height)
    resized = source.resize((int(source.width * scale), int(source.height * scale)), Image.Resampling.LANCZOS)
    left = (resized.width - target_w) // 2
    top = (resized.height - target_h) // 2
    return resized.crop((left, top, left + target_w, top + target_h))


def paste_thumbnail(canvas: Image.Image, draw: ImageDraw.ImageDraw, path: Path, box: tuple[int, int, int, int]) -> None:
    x1, y1, x2, y2 = box
    thumb = fit_cover(path, (x2 - x1, y2 - y1))
    overlay = Image.new("RGBA", thumb.size, (5, 9, 14, 70))
    thumb = Image.alpha_composite(thumb.convert("RGBA"), overlay).convert("RGB")
    mask = Image.new("L", thumb.size, 0)
    mdraw = ImageDraw.Draw(mask)
    mdraw.rounded_rectangle((0, 0, thumb.width, thumb.height), radius=20, fill=255)
    canvas.paste(thumb, (x1, y1), mask)
    rounded(draw, box, 20, (0, 0, 0, 0), "#2A394D", 2)


def draw_mini_visual(draw: ImageDraw.ImageDraw, key: str, box: tuple[int, int, int, int], color: str) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 20, COLORS["panel3"], "#2A394D", 2)

    if key == "BRIDGE":
        text(draw, (x1 + 24, y1 + 22), "mouse NES", F["micro_bold"], COLORS["muted"])
        text(draw, (x2 - 24, y1 + 22), "human prediction", F["micro_bold"], COLORS["muted"], "ra")
        tissues = ["thymus", "skin", "muscle", "liver", "kidney", "eye"]
        for i, tissue in enumerate(tissues):
            col, row = divmod(i, 3)
            cx = x1 + 24 + col * 126
            cy = y1 + 62 + row * 46
            rounded(draw, (cx, cy, cx + 104, cy + 32), 10, COLORS["panel2"], color, 1)
            text(draw, (cx + 52, cy + 16), tissue, F["micro_bold"], COLORS["text"], "mm")
        arrow(draw, x1 + 300, y1 + 136, x2 - 182, y1 + 136, COLORS["dim"], 4)
        rounded(draw, (x2 - 160, y1 + 94, x2 - 26, y1 + 178), 18, "#122234", color, 2)
        text(draw, (x2 - 93, y1 + 126), "AUROC", F["micro_bold"], COLORS["dim"], "mm")
        text(draw, (x2 - 93, y1 + 154), "+0.175", F["small_bold"], COLORS["text"], "mm")
        bar_x, bar_y = x1 + 26, y2 - 58
        rounded(draw, (bar_x, bar_y, x2 - 26, bar_y + 14), 7, "#253247")
        rounded(draw, (bar_x, bar_y, bar_x + 250, bar_y + 14), 7, COLORS["dim"])
        rounded(draw, (bar_x, bar_y + 24, bar_x + 390, bar_y + 38), 7, color)

    elif key == "DECOMPOSE":
        labels = [("spleen -> thymus", 0.62), ("HZE endocrine", 0.44), ("skin response", 0.40)]
        for i, (label, value) in enumerate(labels):
            y = y1 + 62 + i * 58
            text(draw, (x1 + 26, y + 4), label, F["micro_bold"], COLORS["text"])
            bx = x1 + 192
            rounded(draw, (bx, y, x2 - 28, y + 28), 9, "#253247")
            fill_w = int((x2 - 28 - bx) * value)
            rounded(draw, (bx, y, bx + fill_w, y + 28), 9, color)
            text(draw, (x2 - 38, y + 4), f"{int(value * 100)}%", F["micro_bold"], COLORS["text"], "ra")
        rounded(draw, (x1 + 28, y2 - 62, x2 - 28, y2 - 26), 12, "#122234", color, 1)
        text(draw, (x1 + 48, y2 - 52), "interaction-heavy response surfaces", F["micro_bold"], COLORS["muted"])

    elif key == "INTERVENE":
        text(draw, (x1 + 26, y1 + 24), "signature", F["micro_bold"], COLORS["muted"])
        for i in range(9):
            x = x1 + 30 + (i % 3) * 34
            y = y1 + 62 + (i // 3) * 34
            rounded(draw, (x, y, x + 20, y + 20), 5, color if i in {1, 4, 7} else COLORS["panel2"])
        funnel = [(x1 + 176, y1 + 54), (x1 + 310, y1 + 54), (x1 + 270, y1 + 128), (x1 + 216, y1 + 128)]
        draw.polygon(funnel, fill=rgba(color, 80), outline=color)
        text(draw, (x1 + 243, y1 + 92), "rank", F["micro_bold"], COLORS["text"], "mm")
        arrow(draw, x1 + 322, y1 + 92, x2 - 184, y1 + 92, COLORS["dim"], 4)
        rounded(draw, (x2 - 160, y1 + 54, x2 - 28, y1 + 130), 18, "#122234", color, 2)
        text(draw, (x2 - 94, y1 + 82), "215", F["small_bold"], COLORS["text"], "mm")
        text(draw, (x2 - 94, y1 + 110), "hits", F["micro_bold"], COLORS["muted"], "mm")
        chips = [("chemical", COLORS["green"]), ("CRISPR", COLORS["violet"])]
        for i, (label, chip_color) in enumerate(chips):
            cx = x1 + 74 + i * 190
            cy = y2 - 62
            rounded(draw, (cx, cy, cx + 146, cy + 34), 12, COLORS["panel2"], chip_color, 1)
            text(draw, (cx + 73, cy + 17), label, F["micro_bold"], COLORS["text"], "mm")

    elif key == "CAUSAL":
        rows = [("IR", "0.440", "skin", COLORS["amber"]), ("HLU", "0.439", "thymus", COLORS["blue"]), ("HLUxIR", "0.410", "kidney", COLORS["rose"])]
        for i, (stressor, score, tissue, row_color) in enumerate(rows):
            y = y1 + 48 + i * 58
            rounded(draw, (x1 + 28, y, x1 + 128, y + 34), 11, COLORS["panel2"], row_color, 1)
            text(draw, (x1 + 78, y + 17), stressor, F["micro_bold"], COLORS["text"], "mm")
            arrow(draw, x1 + 140, y + 17, x1 + 214, y + 17, COLORS["dim"], 3)
            rounded(draw, (x1 + 226, y, x1 + 336, y + 34), 11, COLORS["panel3"], row_color, 1)
            text(draw, (x1 + 281, y + 17), score, F["micro_bold"], row_color, "mm")
            arrow(draw, x1 + 348, y + 17, x2 - 138, y + 17, COLORS["dim"], 3)
            rounded(draw, (x2 - 126, y, x2 - 28, y + 34), 11, COLORS["panel2"], COLORS["teal"], 1)
            text(draw, (x2 - 77, y + 17), tissue, F["micro_bold"], COLORS["text"], "mm")


def draw_workstream_card(
    canvas: Image.Image,
    draw: ImageDraw.ImageDraw,
    item: dict[str, object],
    box: tuple[int, int, int, int],
) -> None:
    x1, y1, x2, y2 = box
    color = str(item["color"])
    rounded(draw, box, 28, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 28, y1 + 24, x2 - 28, y1 + 24), fill=color, width=5)

    text(draw, (x1 + 34, y1 + 46), str(item["key"]), F["micro_bold"], color)
    draw_mini_visual(draw, str(item["key"]), (x1 + 52, y1 + 98, x1 + 556, y2 - 52), color)

    tx = x1 + 610
    text(draw, (tx, y1 + 52), str(item["title"]), F["h3"], COLORS["text"])
    paragraph(draw, (tx, y1 + 96), str(item["desc"]), F["small"], x2 - tx - 38, COLORS["muted"], 6)

    chip_y = y1 + 210
    chip_w = (x2 - tx - 56) // 2
    for i, metric in enumerate(item["metrics"]):
        cx = tx + i * (chip_w + 28)
        rounded(draw, (cx, chip_y, cx + chip_w, chip_y + 68), 18, COLORS["panel2"], color, 2)
        text(draw, (cx + 20, chip_y + 22), str(metric), F["tiny_bold"], COLORS["text"])

    rounded(draw, (tx, y2 - 82, x2 - 38, y2 - 34), 16, "#122234", color, 1)
    text(draw, (tx + 20, y2 - 67), "output", F["micro_bold"], COLORS["dim"])
    output_text = {
        "BRIDGE": "cross-species prediction handle",
        "DECOMPOSE": "interaction-aware response handle",
        "INTERVENE": "perturbation-priority handle",
        "CAUSAL": "stressor-to-tissue map handle",
    }[str(item["key"])]
    text(draw, (tx + 98, y2 - 67), output_text, F["micro_bold"], COLORS["muted"])


def draw_handoff(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1516, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Act 6 Handoff", F["h2"], COLORS["text"])
    text(draw, (x2 - 80, y1 + 48), "reproducible benchmark system", F["small_bold"], COLORS["teal"], "ra")

    steps = [
        ("v8 workstreams", "BRIDGE, DECOMPOSE, INTERVENE, CAUSAL", COLORS["violet"]),
        ("task registry", "typed task records", COLORS["blue"]),
        ("metric profiles", "shared evaluation rules", COLORS["teal"]),
        ("evaluator runs", "repeatable score surface", COLORS["green"]),
        ("release package", "reviewable benchmark artifact", COLORS["amber"]),
    ]
    node_w, gap = 585, 58
    start_x, node_y = x1 + 180, y1 + 150
    for i, (title, desc, color) in enumerate(steps):
        nx = start_x + i * (node_w + gap)
        rounded(draw, (nx, node_y, nx + node_w, node_y + 96), 20, COLORS["panel2"], color, 2)
        text(draw, (nx + 26, node_y + 17), title, F["small_bold"], COLORS["text"])
        text(draw, (nx + 26, node_y + 55), desc, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            arrow(draw, nx + node_w + 9, node_y + 48, nx + node_w + gap - 15, node_y + 48, COLORS["dim"], 4)


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "v8 bridges benchmark scores to biology workstreams; Act 6 packages those workstreams into a reproducible benchmark system.",
        F["small"],
        3180,
        COLORS["muted"],
        8,
    )


def build() -> None:
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 49 | ACT 5 | v8 SUMMARY", F["kicker"], COLORS["teal"])
    bx = 1960
    bx = badge(draw, bx, 56, "WORKSTREAMS", "4", COLORS["violet"])
    bx = badge(draw, bx, 56, "RESULT CARDS", "4", COLORS["blue"])
    bx = badge(draw, bx, 56, "NEXT ACT", "PLATFORM", COLORS["teal"])
    badge(draw, bx, 56, "OUTPUT", "WORKSTREAM MAP", COLORS["amber"])

    text(draw, (120, 330), "v8 Turns Benchmark Evidence Into Follow-Up Workstreams", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "The four v8 pillars become reusable analysis handles and a platform handoff.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    card_w, card_h = 1770, 390
    positions = [
        (120, 610, 120 + card_w, 610 + card_h),
        (1950, 610, 1950 + card_w, 610 + card_h),
        (120, 1060, 120 + card_w, 1060 + card_h),
        (1950, 1060, 1950 + card_w, 1060 + card_h),
    ]
    for item, box in zip(WORKSTREAMS, positions, strict=True):
        draw_workstream_card(image, draw, item, box)

    draw_handoff(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    ImageOps.colorize(gray, black="#080D14", white="#F3F7FC").save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(gray)
    manifest = {
        "title": "v8 Turns Benchmark Evidence Into Follow-Up Workstreams",
        "claim": "v8 condenses four proof objects into reusable follow-up workstreams and an Act 6 platform handoff.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "input_assets": [str(Path(item["asset"]).relative_to(ROOT)) for item in WORKSTREAMS],
        "input_manifests": [str(Path(item["manifest"]).relative_to(ROOT)) for item in WORKSTREAMS],
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "visible_copy_policy": "Visible copy uses audience-facing labels and takeaway copy.",
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    QA_NOTE.write_text(
        "\n".join(
            [
                "# v8 Summary Asset Visual QA",
                "",
                "Slide 49 summarizes the four v8 proof-object slides and hands the audience to Act 6.",
                "",
                "Checks to review:",
                "- Four proof cards read independently at presentation scale.",
                "- Thumbnail crops do not hide card titles or metric chips.",
                "- Act 6 handoff arrows stay in gutters and do not touch text.",
                "- Footer uses audience-facing takeaway copy only.",
                "- Grayscale version preserves the four-card hierarchy.",
                "",
                f"Final asset: `{OUT_PATH.relative_to(ROOT)}`",
                f"Grayscale asset: `{GRAY_PATH.relative_to(ROOT)}`",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"output": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    build()

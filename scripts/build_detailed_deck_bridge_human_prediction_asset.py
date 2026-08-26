#!/usr/bin/env python3
"""Build slide 45 asset: mouse pathways improve human prediction."""

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
OUT_DIR = ASSET_ROOT / "bridge_human_prediction"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "mouse_pathways_improve_human_prediction_premium.png"
GRAY_PATH = OUT_DIR / "mouse_pathways_improve_human_prediction_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "bridge_human_prediction_manifest.json"
QA_NOTE = OUT_DIR / "BRIDGE_HUMAN_PREDICTION_ASSET_VISUAL_QA.md"

V8_PILLARS_MANIFEST = ASSET_ROOT / "v8_pillars" / "v8_pillars_manifest.json"
V8_FIGURE = ROOT / "v8" / "figures" / "Figure1_Species_Transfer.png"
V8_BRIDGE = ROOT / "v8" / "bridge" / "evaluation" / "supervised_conservation.json"
V8_SPECIES = ROOT / "v8" / "bridge" / "evaluation" / "species_transfer_nes.json"
V8_TISSUE_NES = ROOT / "v8" / "bridge" / "evaluation" / "tissue_nes_spearman.json"

COLORS = {
    "bg": "#0B1119",
    "bg2": "#10141F",
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
    "orange": "#E6A15D",
    "violet": "#B39DFF",
    "rose": "#E17882",
    "red": "#F17C88",
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
    "title": load_font(86, True),
    "subtitle": load_font(37),
    "section": load_font(42, True),
    "h2": load_font(36, True),
    "h3": load_font(30, True),
    "body": load_font(27),
    "body_bold": load_font(27, True),
    "small": load_font(24),
    "small_bold": load_font(24, True),
    "tiny": load_font(21),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "metric": load_font(68, True),
    "metric2": load_font(54, True),
    "heat": load_font(24, True),
    "axis": load_font(18, True),
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
    leading: int = 8,
) -> int:
    x, y = xy
    for line in wrap_lines(draw, body, font, max_width):
        text(draw, (x, y), line, font, fill)
        y += font.size + leading
    return y


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(224, int(draw.textlength(value, font=F["tiny_bold"]) + 76))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def arrow(draw: ImageDraw.ImageDraw, x1: int, y1: int, x2: int, y2: int, color: str, width: int = 4) -> None:
    draw.line((x1, y1, x2, y2), fill=color, width=width)
    pts = [(x2, y2), (x2 - 20, y2 - 11), (x2 - 20, y2 + 11)] if x2 >= x1 else [(x2, y2), (x2 + 20, y2 - 11), (x2 + 20, y2 + 11)]
    draw.polygon(pts, fill=color)


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_data() -> dict[str, object]:
    bridge = load_json(V8_BRIDGE)
    species = load_json(V8_SPECIES)
    tissue = load_json(V8_TISSUE_NES)
    pillars = load_json(V8_PILLARS_MANIFEST)
    return {"bridge": bridge, "species": species, "tissue": tissue, "pillars": pillars}


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


def heat_color(value: float) -> str:
    value = max(-0.8, min(0.8, value))
    if value >= 0:
        return blend("#F8E7C6", "#D45559", value / 0.8)
    return blend("#D9ECF5", "#2469B2", abs(value) / 0.8)


def draw_species_heatmap(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x1, y1, x2, y2 = 120, 610, 1360, 1460
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 36), "Mouse NES carries a tissue signature", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 92),
        "The input signal is a pathway-level correlation map between mouse tissues and I4 cell compartments.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )
    pivot = data["tissue"]["i4_pivot"]
    tissues = ["gastrocnemius", "kidney", "skin", "liver", "thymus", "eye"]
    comps = ["B", "CD14_Mono", "CD16_Mono", "CD4_T", "CD8_T", "DC", "NK", "other_T", "other"]
    labels = ["B", "CD14", "CD16", "CD4", "CD8", "DC", "NK", "T+", "other"]
    cell_w, cell_h = 86, 70
    gx, gy = x1 + 300, y1 + 210
    for c, label in enumerate(labels):
        text(draw, (gx + c * cell_w + cell_w / 2, gy - 34), label, F["axis"], COLORS["muted"], "mm")
    for r, tissue in enumerate(tissues):
        label = "gastroc." if tissue == "gastrocnemius" else tissue
        text(draw, (gx - 26, gy + r * cell_h + cell_h / 2), label, F["tiny_bold"], COLORS["text"], "rm")
        for c, comp in enumerate(comps):
            value = float(pivot[f"{comp}|Immediately Post-flight"][tissue])
            box = (gx + c * cell_w, gy + r * cell_h, gx + c * cell_w + cell_w - 6, gy + r * cell_h + cell_h - 6)
            rounded(draw, box, 8, heat_color(value))
            fill = COLORS["text"] if abs(value) > 0.32 else COLORS["ink"]
            text(draw, ((box[0] + box[2]) / 2, (box[1] + box[3]) / 2), f"{value:+.2f}", F["heat"], fill, "mm")
    text(draw, (gx, gy + len(tissues) * cell_h + 18), "I4 post-flight compartments", F["small_bold"], COLORS["muted"])
    scale_x, scale_y = x1 + 104, y2 - 110
    for i, (label, fill) in enumerate([("direction A", "#2469B2"), ("near zero", "#D9ECF5"), ("direction B", "#D45559")]):
        rounded(draw, (scale_x + i * 160, scale_y, scale_x + i * 160 + 120, scale_y + 28), 10, fill)
        text(draw, (scale_x + i * 160 + 60, scale_y + 42), label, F["micro_bold"], COLORS["muted"], "mm")
    top_pair = data["tissue"]["top10_i4_pairs"][0]
    rounded(draw, (x1 + 44, y2 - 220, x2 - 44, y2 - 145), 18, COLORS["panel2"], COLORS["blue"], 2)
    text(draw, (x1 + 70, y2 - 199), "Top I4 pair", F["micro_bold"], COLORS["dim"])
    text(draw, (x1 + 70, y2 - 169), f"gastrocnemius -> other: r {top_pair['spearman_r']:.3f}", F["tiny_bold"], COLORS["text"])


def x_for_auc(x0: int, width: int, value: float) -> int:
    return x0 + int((value - 0.45) / 0.50 * width)


def draw_auc_bar(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    label: str,
    base: dict[str, float],
    aug: dict[str, float],
    color: str,
    display_base: float,
    display_aug: float,
    height: int = 118,
) -> None:
    axis_x, axis_w = x + 238, 640
    rounded(draw, (axis_x, y + 35, axis_x + axis_w, y + 63), 12, COLORS["panel3"], "#2A394D", 1)
    for tick in [0.5, 0.7, 0.9]:
        tx = x_for_auc(axis_x, axis_w, tick)
        draw.line((tx, y + 20, tx, y + 82), fill=rgba(COLORS["dim"], 120), width=2)
        text(draw, (tx, y + 88), f"{tick:.1f}", F["micro_bold"], COLORS["dim"], "ma")
    bx = x_for_auc(axis_x, axis_w, base["mean"])
    ax = x_for_auc(axis_x, axis_w, aug["mean"])
    b1 = x_for_auc(axis_x, axis_w, base["ci_low"])
    b2 = x_for_auc(axis_x, axis_w, base["ci_high"])
    a1 = x_for_auc(axis_x, axis_w, aug["ci_low"])
    a2 = x_for_auc(axis_x, axis_w, aug["ci_high"])
    draw.line((b1, y + 35, b2, y + 35), fill=COLORS["muted"], width=5)
    draw.line((a1, y + 63, a2, y + 63), fill=color, width=5)
    draw.ellipse((bx - 12, y + 23, bx + 12, y + 47), fill=COLORS["muted"])
    draw.ellipse((ax - 14, y + 49, ax + 14, y + 77), fill=color)
    text(draw, (x, y + 21), label, F["h3"], COLORS["text"])
    text(draw, (x, y + 60), f"{display_base:.3f} -> {display_aug:.3f}", F["small_bold"], color)
    arrow(draw, bx + 24, y + 50, ax - 24, y + 50, color, 4)


def draw_prediction_lift(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x1, y1, x2, y2 = 1415, 610, 2625, 1460
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 36), "Adding mouse NES raises human prediction", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 92),
        "BRIDGE augments 8 I4 pathway features with 6 mouse tissue NES vectors, then predicts conserved human pathway response.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )
    bridge = data["bridge"]
    rf_base = bridge["rf"]["bootstrap_auroc_base"]
    rf_aug = bridge["rf"]["bootstrap_auroc_aug"]
    lr_base = bridge["logreg"]["bootstrap_auroc_base"]
    lr_aug = bridge["logreg"]["bootstrap_auroc_aug"]
    draw_auc_bar(
        draw,
        x1 + 62,
        y1 + 210,
        "Random forest",
        rf_base,
        rf_aug,
        COLORS["blue"],
        bridge["rf"]["cv_mean_auroc_base"],
        bridge["rf"]["cv_mean_auroc_aug"],
    )
    draw_auc_bar(
        draw,
        x1 + 62,
        y1 + 380,
        "Logistic reg.",
        lr_base,
        lr_aug,
        COLORS["teal"],
        bridge["logreg"]["cv_mean_auroc_base"],
        bridge["logreg"]["cv_mean_auroc_aug"],
    )
    rounded(draw, (x1 + 64, y1 + 590, x2 - 64, y1 + 748), 24, COLORS["panel2"], COLORS["blue"], 2)
    text(draw, (x1 + 94, y1 + 622), "Primary readout", F["micro_bold"], COLORS["dim"])
    text(draw, (x1 + 94, y1 + 660), f"+{bridge['rf']['delta_aug_minus_base']['mean']:.3f} AUROC", F["metric"], COLORS["blue"])
    text(
        draw,
        (x2 - 94, y1 + 660),
        f"95% CI {bridge['rf']['delta_aug_minus_base']['ci_low']:.3f} to {bridge['rf']['delta_aug_minus_base']['ci_high']:.3f}",
        F["small_bold"],
        COLORS["muted"],
        "ra",
    )
    text(draw, (x1 + 94, y1 + 724), "Displayed AUROC uses cross-validation mean; interval uses bootstrap resampling.", F["tiny"], COLORS["muted"])


def draw_tissue_contribution(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x1, y1, x2, y2 = 2680, 610, 3720, 1460
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 36), "Tissue features add complementary cues", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 92),
        "Ablation asks how much the logistic model changes when each mouse tissue NES feature is removed.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )
    ablation = data["bridge"]["logreg"]["ablation_leave_one_mouse_out"]
    rows = []
    for feature, stats in ablation.items():
        tissue = feature.replace("mouse_", "").replace("_NES", "")
        drop = max(0.0, -float(stats["delta_vs_full"]))
        rows.append((tissue, drop))
    rows.sort(key=lambda item: item[1], reverse=True)
    chart_x, chart_y = x1 + 84, y1 + 250
    max_drop = 0.021
    for i, (tissue, drop) in enumerate(rows):
        y = chart_y + i * 72
        label = "gastroc." if tissue == "gastrocnemius" else tissue
        text(draw, (chart_x, y + 10), label, F["tiny_bold"], COLORS["text"])
        rounded(draw, (chart_x + 170, y + 10, x2 - 120, y + 40), 12, COLORS["panel3"], "#2A394D", 1)
        width = int((x2 - 120 - (chart_x + 170)) * (drop / max_drop)) if drop > 0 else 0
        fill = COLORS["amber"] if drop >= 0.015 else COLORS["teal"]
        if width:
            rounded(draw, (chart_x + 170, y + 10, chart_x + 170 + width, y + 40), 12, fill)
        value = f"{drop:.3f}" if drop else "stable"
        text(draw, (x2 - 84, y + 9), value, F["tiny_bold"], fill if drop else COLORS["dim"], "ra")
    rounded(draw, (x1 + 64, y2 - 242, x2 - 64, y2 - 78), 24, COLORS["panel2"], COLORS["green"], 2)
    text(draw, (x1 + 96, y2 - 210), "Takeaway", F["micro_bold"], COLORS["dim"])
    paragraph(
        draw,
        (x1 + 96, y2 - 178),
        "The classifier sees a six-tissue pathway profile, so performance can improve even when individual tissue directions differ.",
        F["small"],
        x2 - x1 - 192,
        COLORS["muted"],
        8,
    )


def draw_flow_strip(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1522, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "BRIDGE decoder", F["h2"], COLORS["text"])
    steps = [
        ("Human features", "8 I4 pathway inputs", COLORS["blue"]),
        ("Mouse prior", "6 tissue NES vectors", COLORS["teal"]),
        ("Classifier", "RF + logistic readout", COLORS["violet"]),
        ("Validation", "5-fold CV + bootstrap", COLORS["green"]),
        ("Readout lane", "pathway-conservation prediction", COLORS["amber"]),
    ]
    node_w, gap = 570, 62
    start_x, y = x1 + 460, y1 + 152
    for i, (title, desc, color) in enumerate(steps):
        nx = start_x + i * (node_w + gap)
        rounded(draw, (nx, y, nx + node_w, y + 94), 20, COLORS["panel2"], color, 2)
        text(draw, (nx + 26, y + 17), title, F["small_bold"], COLORS["text"])
        text(draw, (nx + 26, y + 55), desc, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            arrow(draw, nx + node_w + 9, y + 47, nx + node_w + gap - 14, y + 47, COLORS["dim"], 4)
    text(draw, (x2 - 80, y1 + 72), "BRIDGE", F["metric2"], COLORS["blue"], "ra")
    text(draw, (x2 - 80, y1 + 124), "mouse -> human pathway prior", F["small_bold"], COLORS["muted"], "ra")


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "BRIDGE readout: mouse tissue NES improves supervised human pathway-conservation prediction.",
        F["small"],
        2450,
        COLORS["muted"],
        8,
    )
    text(draw, (3580, 1960), "Inputs", F["micro_bold"], COLORS["dim"], "ra")
    text(draw, (3580, 1992), "v8 Figure 1 + bridge evaluation JSON", F["tiny"], COLORS["muted"], "ra")
    text(draw, (3580, 2024), "species-transfer NES + tissue NES summary", F["tiny"], COLORS["muted"], "ra")


def build() -> None:
    data = load_data()
    image = background()
    draw = ImageDraw.Draw(image)
    bridge = data["bridge"]
    species = data["species"]

    text(draw, (120, 72), "SLIDE 45 | ACT 5 | BRIDGE", F["kicker"], COLORS["blue"])
    bx = 2040
    bx = badge(draw, bx, 56, "PATHWAYS", f"{bridge['n_pathways']} supervised", COLORS["blue"])
    bx = badge(draw, bx, 56, "MOUSE NES", "6 tissues", COLORS["teal"])
    bx = badge(draw, bx, 56, "RF AUROC", f"{bridge['rf']['cv_mean_auroc_aug']:.3f}", COLORS["green"])
    badge(draw, bx, 56, "I4", f"mean r {species['by_mission']['Inspiration4']['mean_spearman_r']:.3f}", COLORS["amber"])

    text(draw, (120, 330), "Mouse Pathways Improve Human Prediction", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "BRIDGE turns mouse tissue pathway NES into a supervised prior for human pathway-conservation prediction.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_species_heatmap(draw, data)
    draw_prediction_lift(draw, data)
    draw_tissue_contribution(draw, data)
    draw_flow_strip(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    ImageOps.colorize(gray, black="#080D14", white="#F3F7FC").save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(gray)

    manifest = {
        "title": "Mouse Pathways Improve Human Prediction",
        "readout": "Adding six mouse-tissue NES vectors improves supervised human pathway-conservation prediction.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "source_manifests": {"v8_pillars": str(V8_PILLARS_MANIFEST.relative_to(ROOT))},
        "source_files": {
            "v8_figure1": str(V8_FIGURE.relative_to(ROOT)),
            "supervised_conservation": str(V8_BRIDGE.relative_to(ROOT)),
            "species_transfer_nes": str(V8_SPECIES.relative_to(ROOT)),
            "tissue_nes_spearman": str(V8_TISSUE_NES.relative_to(ROOT)),
        },
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": {
            "supervised_pathways": bridge["n_pathways"],
            "positive_pathways": bridge["n_positive"],
            "mouse_tissue_features": len(bridge["mouse_features"]),
            "rf_base_auroc": bridge["rf"]["cv_mean_auroc_base"],
            "rf_aug_auroc": bridge["rf"]["cv_mean_auroc_aug"],
            "rf_delta": bridge["rf"]["delta_aug_minus_base"]["mean"],
            "rf_delta_ci_low": bridge["rf"]["delta_aug_minus_base"]["ci_low"],
            "rf_delta_ci_high": bridge["rf"]["delta_aug_minus_base"]["ci_high"],
            "logreg_base_auroc": bridge["logreg"]["cv_mean_auroc_base"],
            "logreg_aug_auroc": bridge["logreg"]["cv_mean_auroc_aug"],
            "logreg_delta": bridge["logreg"]["delta_aug_minus_base"]["mean"],
            "i4_mean_spearman": species["by_mission"]["Inspiration4"]["mean_spearman_r"],
            "i4_max_spearman": species["by_mission"]["Inspiration4"]["max_spearman_r"],
            "tissue_nes_pairs": data["tissue"]["n_pairs"],
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    QA_NOTE.write_text(
        "\n".join(
            [
                "# BRIDGE Human Prediction Asset Visual QA",
                "",
                "Slide 45 rebuilds v8 Figure 1 into a presentation-scale proof object.",
                "",
                "Checks to review:",
                "- Heatmap labels and values remain legible at presentation scale.",
                "- AUROC lift panel clearly shows base-to-augmented movement and CI text.",
                "- Tissue contribution bars do not imply treatment or operational guidance.",
                "- Bottom BRIDGE reading strip has stable arrow spacing.",
                "- Grayscale version preserves heatmap and AUROC contrast.",
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

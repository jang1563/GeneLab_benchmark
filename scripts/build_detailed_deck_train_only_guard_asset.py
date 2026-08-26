#!/usr/bin/env python3
"""Build the detailed-deck train-only leakage guard methodology asset."""

from __future__ import annotations

import csv
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
    / "train_only_guard"
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
    "teal": "#5FD3C4",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "red": "#E16D73",
    "violet": "#9C7CFF",
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
    "title": load_font(80, True),
    "subtitle": load_font(36, False),
    "h2": load_font(44, True),
    "h3": load_font(34, True),
    "body": load_font(30, False),
    "body_bold": load_font(30, True),
    "small": load_font(27, False),
    "small_bold": load_font(27, True),
    "tiny": load_font(23, False),
    "tiny_bold": load_font(23, True),
    "stat": load_font(86, True),
    "mission": load_font(38, True),
}


def rounded(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], radius: int, fill: str, outline: str | None = None, width: int = 1) -> None:
    draw.rounded_rectangle(box, radius=radius, fill=fill, outline=outline, width=width)


def text(draw: ImageDraw.ImageDraw, xy: tuple[int, int], value: str, font: ImageFont.ImageFont, fill: str = COLORS["text"], anchor: str | None = None) -> None:
    draw.text(xy, value, font=font, fill=fill, anchor=anchor)


def wrap_lines(draw: ImageDraw.ImageDraw, label: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
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
    return lines


def multiline(draw: ImageDraw.ImageDraw, xy: tuple[int, int], lines: Iterable[str], font: ImageFont.ImageFont, fill: str, leading: int = 8) -> int:
    x, y = xy
    for line in lines:
        draw.text((x, y), line, font=font, fill=fill)
        y += font.size + leading
    return y


def paragraph(draw: ImageDraw.ImageDraw, xy: tuple[int, int], body: str, font: ImageFont.ImageFont, max_width: int, fill: str, leading: int = 8) -> int:
    return multiline(draw, xy, wrap_lines(draw, body, font, max_width), font, fill, leading)


def arrow(draw: ImageDraw.ImageDraw, start: tuple[int, int], end: tuple[int, int], color: str, width: int = 5) -> None:
    x0, y0 = start
    x1, y1 = end
    draw.line((x0, y0, x1, y1), fill=color, width=width)
    if x1 >= x0:
        draw.polygon([(x1, y1), (x1 - 22, y1 - 13), (x1 - 22, y1 + 13)], fill=color)
    else:
        draw.polygon([(x1, y1), (x1 + 22, y1 - 13), (x1 + 22, y1 + 13)], fill=color)


def load_data() -> dict:
    fold_dir = ROOT / "tasks" / "A4_thymus_lomo" / "fold_RR-9_test"
    fold_info = json.loads((fold_dir / "fold_info.json").read_text())
    selected_genes = [line.strip() for line in (fold_dir / "selected_genes.txt").read_text().splitlines() if line.strip()]
    with (fold_dir / "train_X.csv").open() as handle:
        train_cols = len(next(csv.reader(handle)))
        train_rows = sum(1 for _ in handle)
    with (fold_dir / "test_X.csv").open() as handle:
        test_cols = len(next(csv.reader(handle)))
        test_rows = sum(1 for _ in handle)
    return {
        "fold_dir": fold_dir,
        "fold_info": fold_info,
        "selected_gene_count": len(selected_genes),
        "train_rows": train_rows,
        "test_rows": test_rows,
        "train_feature_cols": train_cols - 1,
        "test_feature_cols": test_cols - 1,
        "first_genes": selected_genes[:8],
    }


def draw_badges(draw: ImageDraw.ImageDraw, data: dict) -> None:
    fold_info = data["fold_info"]
    badges = [
        ("EXAMPLE FOLD", "A4 thymus RR-9", 360),
        ("FROZEN FEATURES", f"{data['selected_gene_count']:,} genes", 370),
        ("HIDDEN USED", "0 for choices", 310),
    ]
    bx = 2310
    for kicker, body, badge_w in badges:
        rounded(draw, (bx, 72, bx + badge_w, 176), 26, "#121B28", "#2A394D", 2)
        text(draw, (bx + 28, 96), kicker, F["tiny"], COLORS["teal"])
        text(draw, (bx + 28, 126), body, F["small_bold"], COLORS["text"])
        bx += badge_w + 30
    text(draw, (3498, 184), f"hidden: {fold_info['test_mission']}", F["tiny"], COLORS["muted"], anchor="ra")


def draw_mission_node(draw: ImageDraw.ImageDraw, cx: int, cy: int, label: str, color: str, *, hidden: bool = False) -> None:
    fill = "#151F2D" if not hidden else "#241719"
    outline = color
    rounded(draw, (cx - 92, cy - 55, cx + 92, cy + 55), 32, fill, outline, 3)
    text(draw, (cx, cy - 20), label, F["mission"], COLORS["text"], anchor="ma")
    if hidden:
        text(draw, (cx, cy + 18), "hidden", F["tiny_bold"], COLORS["red"], anchor="ma")


def draw_bad_choice_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 45, y0 + 45), "A. Where leakage enters", F["h2"], COLORS["text"])
    paragraph(draw, (x0 + 45, y0 + 98), "The hidden mission stays outside preprocessing and feature choice.", F["small"], x1 - x0 - 90, COLORS["muted"], 7)

    rail_y = y0 + 330
    draw_mission_node(draw, x0 + 150, rail_y, "MHU-1", COLORS["teal"])
    draw_mission_node(draw, x0 + 355, rail_y, "MHU-2", COLORS["teal"])
    draw_mission_node(draw, x0 + 560, rail_y, "RR-6", COLORS["teal"])
    draw_mission_node(draw, x0 + 150, rail_y + 240, "RR-9", COLORS["red"], hidden=True)

    choose_box = (x0 + 380, rail_y + 155, x1 - 55, rail_y + 350)
    rounded(draw, choose_box, 26, "#241719", COLORS["red"], 3)
    text(draw, (choose_box[0] + 34, choose_box[1] + 32), "choose features", F["h3"], COLORS["text"])
    paragraph(draw, (choose_box[0] + 34, choose_box[1] + 82), "If RR-9 helps choose features, the score is contaminated.", F["small"], choose_box[2] - choose_box[0] - 68, COLORS["muted"], 7)

    arrow(draw, (x0 + 653, rail_y), (choose_box[0] - 25, choose_box[1] + 64), COLORS["teal"], 5)
    arrow(draw, (x0 + 242, rail_y + 240), (choose_box[0] - 25, choose_box[1] + 132), COLORS["red"], 5)
    text(draw, (x0 + 310, rail_y + 300), "bad path", F["small_bold"], COLORS["red"], anchor="ma")

    rounded(draw, (x0 + 45, y1 - 430, x1 - 45, y1 - 250), 24, "#151F2D", "#2A394D", 2)
    text(draw, (x0 + 78, y1 - 394), "Contamination rule", F["h3"], COLORS["amber"])
    paragraph(draw, (x0 + 78, y1 - 342), "Feature selection, scaling, and model fitting are training-side decisions.", F["small"], x1 - x0 - 156, COLORS["text"], 7)

    rounded(draw, (x0 + 45, y1 - 182, x1 - 45, y1 - 62), 20, "#101823", "#2A394D", 2)
    text(draw, (x0 + 78, y1 - 150), "Correct rule", F["small_bold"], COLORS["teal"])
    text(draw, (x0 + 78, y1 - 110), "RR-9 enters only after choices are frozen.", F["small"], COLORS["muted"])


def draw_step_card(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], title: str, body: str, color: str) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 24, "#151F2D", "#2A394D", 2)
    text(draw, (x0 + 24, y0 + 22), title, F["small_bold"], color)
    paragraph(draw, (x0 + 24, y0 + 66), body, F["small"], x1 - x0 - 48, COLORS["text"], 5)


def draw_pipeline_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict) -> None:
    x0, y0, x1, y1 = box
    fold_info = data["fold_info"]
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 45), "B. Correct fold pipeline", F["h2"], COLORS["text"])
    paragraph(draw, (x0 + 50, y0 + 98), "Training missions define the feature list and model; RR-9 is scored only after those choices are fixed.", F["small"], x1 - x0 - 100, COLORS["muted"], 7)

    train_y = y0 + 340
    hidden_y = y0 + 730
    train_nodes = [x0 + 160, x0 + 350, x0 + 540]
    draw.line((train_nodes[0] + 100, train_y, train_nodes[-1] + 102, train_y), fill=rgba(COLORS["teal"], 140), width=7)
    for cx, mission in zip(train_nodes, fold_info["train_missions"]):
        draw_mission_node(draw, cx, train_y, mission, COLORS["teal"])
    draw_mission_node(draw, x0 + 175, hidden_y, fold_info["test_mission"], COLORS["red"], hidden=True)

    steps = [
        ((x0 + 700, y0 + 250, x0 + 900, y0 + 445), "choose", "train variance only", COLORS["blue"]),
        ((x0 + 940, y0 + 250, x0 + 1140, y0 + 445), "freeze", f"{fold_info['n_genes_after_var_filter']:,} genes", COLORS["teal"]),
        ((x0 + 1180, y0 + 250, x0 + 1360, y0 + 445), "fit", "train labels only", COLORS["violet"]),
    ]
    for step_box, title, body, color in steps:
        draw_step_card(draw, step_box, title, body, color)
    arrow(draw, (train_nodes[-1] + 105, train_y), (steps[0][0][0] - 30, train_y), COLORS["teal"], 5)
    for prev, nxt in zip(steps, steps[1:]):
        arrow(draw, (prev[0][2] + 8, train_y), (nxt[0][0] - 10, train_y), COLORS["axis"], 4)

    guard_x = x0 + 1415
    draw.line((guard_x, y0 + 190, guard_x, y1 - 240), fill=COLORS["red"], width=6)
    text(draw, (guard_x, y0 + 160), "guard", F["small_bold"], COLORS["red"], anchor="ma")

    rounded(draw, (guard_x - 255, train_y + 92, guard_x - 35, train_y + 168), 18, "#101823", "#2A394D", 2)
    text(draw, (guard_x - 145, train_y + 108), "frozen choices", F["tiny_bold"], COLORS["text"], anchor="ma")
    text(draw, (guard_x - 145, train_y + 135), "features + fit params", F["tiny"], COLORS["muted"], anchor="ma")
    arrow(draw, (steps[-1][0][2] + 6, train_y), (guard_x - 25, train_y), COLORS["axis"], 5)

    draw.line((x0 + 260, hidden_y, guard_x - 48, hidden_y), fill=rgba(COLORS["red"], 130), width=6)
    text(draw, (x0 + 430, hidden_y + 34), "waits outside feature choice", F["small_bold"], COLORS["red"])
    score_box = (guard_x + 50, hidden_y - 68, x1 - 60, hidden_y + 68)
    arrow(draw, (guard_x + 30, hidden_y), (score_box[0] - 18, hidden_y), COLORS["red"], 5)
    rounded(draw, score_box, 24, "#241719", COLORS["red"], 3)
    text(draw, ((score_box[0] + score_box[2]) / 2, hidden_y - 33), "score once", F["small_bold"], COLORS["text"], anchor="ma")
    text(draw, ((score_box[0] + score_box[2]) / 2, hidden_y + 4), "RR-9 n=20", F["tiny_bold"], COLORS["muted"], anchor="ma")

    evidence_y = y0 + 1035
    rounded(draw, (x0 + 70, evidence_y, x1 - 70, evidence_y + 260), 28, "#151F2D", "#2A394D", 2)
    text(draw, (x0 + 105, evidence_y + 35), "What freezes before scoring", F["h3"], COLORS["teal"])
    x = x0 + 110
    bars = [
        ("all measured genes", fold_info["n_genes_before_var_filter"], COLORS["dim"], 0.92),
        ("train-variance selected", fold_info["n_genes_after_var_filter"], COLORS["teal"], 0.69),
    ]
    bar_y = evidence_y + 105
    for label, value, color, frac in bars:
        text(draw, (x, bar_y - 4), label, F["small_bold"], COLORS["text"])
        rounded(draw, (x + 360, bar_y, x + 950, bar_y + 34), 16, "#101823", "#2A394D", 1)
        rounded(draw, (x + 360, bar_y, x + 360 + int(590 * frac), bar_y + 34), 16, color, None, 0)
        text(draw, (x + 980, bar_y - 6), f"{value:,}", F["small_bold"], COLORS["text"])
        bar_y += 76
    text(draw, (x0 + 105, evidence_y + 220), "Selected genes are then applied to both train_X and test_X without recomputing on RR-9.", F["small"], COLORS["muted"])

    rounded(draw, (x0 + 70, y1 - 125, x1 - 70, y1 - 45), 20, "#101823", "#2A394D", 2)
    text(draw, (x0 + 100, y1 - 98), "Takeaway", F["small_bold"], COLORS["amber"])
    text(draw, (x0 + 285, y1 - 98), "Held-out status depends on freezing features, scaling, and fit before RR-9 enters.", F["small"], COLORS["muted"])


def draw_audit_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict) -> None:
    x0, y0, x1, y1 = box
    fold_info = data["fold_info"]
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 45, y0 + 45), "C. Fold-file audit hooks", F["h2"], COLORS["text"])
    paragraph(draw, (x0 + 45, y0 + 98), "The example fold exposes the split, frozen feature list, and matching train/test matrices.", F["small"], x1 - x0 - 90, COLORS["muted"], 7)

    stat_y = y0 + 180
    stats = [
        (str(fold_info["n_train"]), "train samples"),
        (str(fold_info["n_test"]), "hidden samples"),
        (f"{data['selected_gene_count']:,}", "frozen genes"),
    ]
    sx = x0 + 45
    for val, label in stats:
        stat_w = 300
        rounded(draw, (sx, stat_y, sx + stat_w, stat_y + 150), 24, "#151F2D", "#2A394D", 2)
        val_x = int(sx + (stat_w - draw.textlength(val, font=F["stat"])) / 2)
        label_x = int(sx + (stat_w - draw.textlength(label, font=F["tiny_bold"])) / 2)
        text(draw, (val_x, stat_y + 12), val, F["stat"], COLORS["text"])
        text(draw, (label_x, stat_y + 112), label, F["tiny_bold"], COLORS["muted"])
        sx += stat_w + 24

    rows = [
        ("fold_info.json", "missions + before/after gene counts"),
        ("selected_genes.txt", f"{data['selected_gene_count']:,} train-selected feature IDs"),
        ("train_X.csv", f"{data['train_rows']} samples x {data['train_feature_cols']:,} genes"),
        ("test_X.csv", f"{data['test_rows']} samples x same frozen genes"),
        ("Category A rule", "logFC features excluded; label-leakage guard"),
    ]
    row_y = y0 + 405
    for name, detail in rows:
        rounded(draw, (x0 + 45, row_y, x1 - 45, row_y + 112), 18, "#151F2D", "#2A394D", 1)
        rounded(draw, (x0 + 75, row_y + 30, x0 + 125, row_y + 80), 15, "#101823", COLORS["teal"], 2)
        text(draw, (x0 + 100, row_y + 40), "ok", F["tiny_bold"], COLORS["teal"], anchor="ma")
        text(draw, (x0 + 150, row_y + 22), name, F["small_bold"], COLORS["text"])
        text(draw, (x0 + 150, row_y + 62), detail, F["tiny"], COLORS["muted"])
        row_y += 126

    y = y1 - 300
    rounded(draw, (x0 + 45, y, x1 - 45, y + 247), 28, "#211E17", "#69532B", 2)
    text(draw, (x0 + 78, y + 36), "Readout frame", F["h3"], COLORS["amber"])
    lines = [
        "Supports leakage control for fold construction.",
        "Performance and mechanism are evaluated in companion slides.",
        "The same guard applies to every fold and feature view.",
    ]
    multiline(draw, (x0 + 78, y + 90), lines, F["small"], COLORS["text"], 9)


def main() -> None:
    data = load_data()
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)

    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 48), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 42), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "LEAKAGE GUARD", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "Training choices stop before the hidden mission", F["title"])
    text(draw, (150, 216), "Within each fold, feature selection, scaling, and fitting are learned on training missions only.", F["subtitle"], COLORS["muted"])
    draw_badges(draw, data)

    draw_bad_choice_panel(draw, (150, 350, 905, 1800), data)
    draw_pipeline_panel(draw, (945, 350, 2600, 1800), data)
    draw_audit_panel(draw, (2635, 350, 3690, 1800), data)

    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    text(draw, (205, 1925), "Takeaway", F["small_bold"], COLORS["blue"])
    source = "Feature selection, scaling, and fitting happen on training missions before the hidden mission is scored."
    paragraph(draw, (390, 1925), source, F["small"], 3140, COLORS["muted"], 7)
    text(draw, (205, 1995), "Next", F["small_bold"], COLORS["amber"])
    scope = "With leakage control established, the metric primer shows AUROC, uncertainty, and permutation context."
    paragraph(draw, (390, 1995), scope, F["small"], 3140, COLORS["muted"], 7)

    png = OUT_DIR / "train_only_feature_guard_premium.png"
    canvas.save(png, quality=95)
    gray = OUT_DIR / "train_only_feature_guard_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    manifest = {
        "title": "Training choices stop before the hidden mission",
        "source_files": [
            "README.md",
            "tasks/README.md",
            "tasks/A4_thymus_lomo/fold_RR-9_test/fold_info.json",
            "tasks/A4_thymus_lomo/fold_RR-9_test/selected_genes.txt",
            "tasks/A4_thymus_lomo/fold_RR-9_test/train_X.csv",
            "tasks/A4_thymus_lomo/fold_RR-9_test/test_X.csv",
        ],
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "summary": {
            "example_fold": "A4_thymus_lomo/fold_RR-9_test",
            "train_missions": data["fold_info"]["train_missions"],
            "hidden_mission": data["fold_info"]["test_mission"],
            "n_train": data["fold_info"]["n_train"],
            "n_test": data["fold_info"]["n_test"],
            "genes_before_var_filter": data["fold_info"]["n_genes_before_var_filter"],
            "genes_after_var_filter": data["fold_info"]["n_genes_after_var_filter"],
            "selected_gene_count": data["selected_gene_count"],
            "train_feature_cols": data["train_feature_cols"],
            "test_feature_cols": data["test_feature_cols"],
        },
        "scope": (
            "The asset shows leakage-control protocol evidence for train-only feature selection and fitting. "
            "Companion slides evaluate model performance, mechanism, and biological interpretation."
        ),
    }
    manifest_path = OUT_DIR / "train_only_feature_guard_premium_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

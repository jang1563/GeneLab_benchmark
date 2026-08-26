#!/usr/bin/env python3
"""Build the detailed-deck mission-held-out methodology proof asset.

The slide teaches the evaluation split for mixed ML / biology audiences:
the hidden unit is an entire mission held out from training-side choices.
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
    / "mission_holdout"
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
    "liver": "#A6AEBB",
    "gastrocnemius": "#56B4E9",
    "kidney": "#B39DFF",
    "thymus": "#E8A34A",
    "skin": "#8BD17C",
    "eye": "#73A7FF",
}

TASKS = [
    ("A1_liver_lomo", "Liver", "liver"),
    ("A2_gastrocnemius_lomo", "Gastrocnemius", "gastrocnemius"),
    ("A3_kidney_lomo", "Kidney", "kidney"),
    ("A4_thymus_lomo", "Thymus", "thymus"),
    ("A5_skin_lomo", "Skin", "skin"),
    ("A6_eye_lomo", "Eye", "eye"),
]


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
    "body_bold": load_font(30, True),
    "small": load_font(27, False),
    "small_bold": load_font(27, True),
    "tiny": load_font(23, False),
    "tiny_bold": load_font(23, True),
    "stat": load_font(90, True),
    "mission": load_font(42, True),
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
    draw.polygon([(x1, y1), (x1 - 22, y1 - 13), (x1 - 22, y1 + 13)], fill=color)


def load_data() -> dict:
    rows = []
    unique_missions: set[str] = set()
    total_folds = 0
    for task_id, label, color_key in TASKS:
        data = json.loads((ROOT / "tasks" / task_id / "task_info.json").read_text())
        folds = data["folds"]
        test_missions = [fold["test_mission"] for fold in folds]
        for fold in folds:
            unique_missions.add(fold["test_mission"])
            unique_missions.update(fold["train_missions"])
        total_folds += data["n_folds_generated"]
        rows.append(
            {
                "task_id": task_id.split("_")[0],
                "label": label,
                "color": COLORS[color_key],
                "n_missions": data["n_missions"],
                "n_folds": data["n_folds_generated"],
                "n_samples": sum(fold["n_test"] for fold in folds),
                "test_missions": test_missions,
                "test_min": min(fold["n_test"] for fold in folds),
                "test_max": max(fold["n_test"] for fold in folds),
            }
        )
    thymus = json.loads((ROOT / "tasks" / "A4_thymus_lomo" / "task_info.json").read_text())
    example = next(fold for fold in thymus["folds"] if fold["test_mission"] == "RR-9")
    return {
        "rows": rows,
        "unique_missions": sorted(unique_missions),
        "total_folds": total_folds,
        "example": example,
    }


def draw_badges(draw: ImageDraw.ImageDraw, data: dict) -> None:
    badges = [
        ("CORE TASKS", "6 tissues", 290),
        ("PUBLIC FOLDS", f"{data['total_folds']} LOMO", 330),
        ("SPLIT UNIT", "mission", 300),
    ]
    bx = 2460
    for kicker, body, badge_w in badges:
        rounded(draw, (bx, 72, bx + badge_w, 176), 26, "#121B28", "#2A394D", 2)
        text(draw, (bx + 28, 96), kicker, F["tiny"], COLORS["teal"])
        text(draw, (bx + 28, 126), body, F["small_bold"], COLORS["text"])
        bx += badge_w + 30


def draw_random_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 45, y0 + 45), "A. Why mission units matter", F["h2"], COLORS["text"])
    paragraph(draw, (x0 + 45, y0 + 98), "Random sample splits can put the same mission on both sides of evaluation.", F["small"], x1 - x0 - 90, COLORS["muted"], 7)

    missions = [("M1", COLORS["teal"]), ("M2", COLORS["blue"]), ("M3", COLORS["amber"])]
    y = y0 + 235
    for idx, (mission, color) in enumerate(missions):
        mx = x0 + 90 + idx * 175
        rounded(draw, (mx, y, mx + 130, y + 92), 36, "#151F2D", color, 3)
        text(draw, (mx + 65, y + 26), mission, F["h3"], COLORS["text"], anchor="ma")
        for k in range(4):
            draw.ellipse((mx + 18 + k * 25, y + 61, mx + 32 + k * 25, y + 75), fill=color)
    text(draw, (x0 + 110, y + 140), "train", F["small_bold"], COLORS["teal"])
    text(draw, (x0 + 410, y + 140), "test", F["small_bold"], COLORS["red"])
    draw.line((x0 + 260, y + 135, x0 + 370, y + 135), fill=COLORS["axis"], width=4)
    draw.line((x0 + 315, y + 110, x0 + 315, y + 165), fill=COLORS["axis"], width=4)
    text(draw, (x0 + 315, y + 205), "same mission can leak context", F["small_bold"], COLORS["red"], anchor="ma")

    leak_box = (x0 + 45, y + 300, x1 - 45, y + 650)
    rounded(draw, leak_box, 24, "#151F2D", "#2A394D", 2)
    text(draw, (leak_box[0] + 32, leak_box[1] + 30), "What leaks", F["h3"], COLORS["text"])
    leak_rows = [
        ("mission handling", "same background"),
        ("batch / protocol", "same context"),
        ("expression baseline", "shortcut signal"),
    ]
    row_y = leak_box[1] + 88
    for label, note in leak_rows:
        rounded(draw, (leak_box[0] + 32, row_y, leak_box[0] + 92, row_y + 50), 16, "#241719", COLORS["red"], 2)
        text(draw, (leak_box[0] + 62, row_y + 12), "X", F["small_bold"], COLORS["red"], anchor="ma")
        text(draw, (leak_box[0] + 115, row_y + 2), label, F["small_bold"], COLORS["text"])
        text(draw, (leak_box[0] + 115, row_y + 32), note, F["tiny"], COLORS["muted"])
        row_y += 78

    rounded(draw, (x0 + 45, y1 - 270, x1 - 45, y1 - 122), 24, "#211E17", "#69532B", 2)
    text(draw, (x0 + 78, y1 - 238), "Problem", F["h3"], COLORS["amber"])
    paragraph(draw, (x0 + 78, y1 - 190), "The model may learn mission context instead of transferable spaceflight biology.", F["small"], x1 - x0 - 156, COLORS["text"], 7)

    rounded(draw, (x0 + 45, y1 - 88, x1 - 45, y1 - 38), 16, "#101823", "#2A394D", 1)
    text(draw, (x0 + 70, y1 - 76), "Use mission-held-out folds instead.", F["small_bold"], COLORS["teal"])


def draw_mission_node(draw: ImageDraw.ImageDraw, cx: int, cy: int, label: str, color: str, mode: str, n: int | None = None) -> None:
    fill = "#151F2D" if mode == "train" else "#241719"
    outline = color if mode == "train" else COLORS["red"]
    rounded(draw, (cx - 100, cy - 64, cx + 100, cy + 64), 38, fill, outline, 4)
    text(draw, (cx, cy - 25), label, F["mission"], COLORS["text"], anchor="ma")
    if n is not None:
        text(draw, (cx, cy + 22), f"n={n}", F["small_bold"], COLORS["muted"], anchor="ma")


def draw_protocol_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], example: dict) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 50, y0 + 45), "B. One fold hides one mission", F["h2"], COLORS["text"])
    paragraph(draw, (x0 + 50, y0 + 98), "Example: A4 thymus fold_RR-9_test trains on three missions and scores RR-9 only after fitting.", F["small"], x1 - x0 - 100, COLORS["muted"], 7)

    rail_y = y0 + 365
    train_xs = [x0 + 190, x0 + 445, x0 + 700]
    draw.line((train_xs[0] + 110, rail_y, train_xs[-1] + 125, rail_y), fill=rgba(COLORS["teal"], 180), width=8)
    for cx, mission in zip(train_xs, example["train_missions"]):
        draw_mission_node(draw, cx, rail_y, mission, COLORS["teal"], "train")

    fit_box = (x0 + 910, rail_y - 88, x0 + 1205, rail_y + 88)
    rounded(draw, fit_box, 26, "#151F2D", "#2A394D", 2)
    text(draw, (fit_box[0] + 35, fit_box[1] + 34), "fit choices", F["h3"], COLORS["text"])
    text(draw, (fit_box[0] + 35, fit_box[1] + 82), "features + model", F["small"], COLORS["muted"])
    arrow(draw, (train_xs[-1] + 125, rail_y), (fit_box[0] - 30, rail_y), COLORS["teal"], 5)

    wall_x = x0 + 1300
    draw.line((wall_x, y0 + 235, wall_x, y1 - 235), fill=COLORS["red"], width=6)
    text(draw, (wall_x, y0 + 205), "holdout gate", F["small_bold"], COLORS["red"], anchor="ma")

    test_cx = x0 + 1515
    draw_mission_node(draw, test_cx, rail_y, example["test_mission"], COLORS["red"], "test", example["n_test"])
    arrow(draw, (fit_box[2] + 30, rail_y), (wall_x - 45, rail_y), COLORS["axis"], 5)
    arrow(draw, (wall_x + 45, rail_y), (test_cx - 120, rail_y), COLORS["red"], 5)

    train_card = (x0 + 70, y0 + 650, x0 + 785, y0 + 885)
    hidden_card = (x0 + 835, y0 + 650, x1 - 70, y0 + 885)
    rounded(draw, train_card, 24, "#151F2D", "#2A394D", 2)
    rounded(draw, hidden_card, 24, "#211E17", "#69532B", 2)
    text(draw, (train_card[0] + 35, train_card[1] + 32), "Training side", F["h3"], COLORS["teal"])
    paragraph(draw, (train_card[0] + 35, train_card[1] + 82), f"{example['n_train']} samples; {example['train_label_distribution']['Flight']} Flight / {example['train_label_distribution']['Ground']} Ground. Variance filter and model fit happen here.", F["small"], train_card[2] - train_card[0] - 70, COLORS["text"], 6)
    text(draw, (hidden_card[0] + 35, hidden_card[1] + 32), "Hidden mission", F["h3"], COLORS["amber"])
    paragraph(draw, (hidden_card[0] + 35, hidden_card[1] + 82), f"{example['test_mission']} stays outside all fitting; {example['n_test']} samples are scored once after the training-side choices are fixed.", F["small"], hidden_card[2] - hidden_card[0] - 70, COLORS["text"], 6)

    rounded(draw, (x0 + 70, y1 - 125, x1 - 70, y1 - 45), 20, "#101823", "#2A394D", 2)
    text(draw, (x0 + 100, y1 - 98), "Takeaway", F["small_bold"], COLORS["amber"])
    text(draw, (x0 + 285, y1 - 98), "Frozen choices precede hidden-mission scoring.", F["small"], COLORS["muted"])


def draw_inventory_panel(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], data: dict) -> None:
    x0, y0, x1, y1 = box
    rounded(draw, box, 34, COLORS["panel"], "#29374A", 2)
    text(draw, (x0 + 45, y0 + 45), "C. Fold inventory", F["h2"], COLORS["text"])
    paragraph(draw, (x0 + 45, y0 + 98), "Each tissue task repeats the same rule: one mission becomes the test fold.", F["small"], x1 - x0 - 90, COLORS["muted"], 7)

    stat_y = y0 + 165
    stats = [("6", "core tissues"), (str(data["total_folds"]), "LOMO folds"), (str(len(data["unique_missions"])), "mission labels")]
    sx = x0 + 45
    stat_w = 295
    for val, label in stats:
        rounded(draw, (sx, stat_y, sx + stat_w, stat_y + 148), 24, "#151F2D", "#2A394D", 2)
        val_x = int(sx + (stat_w - draw.textlength(val, font=F["stat"])) / 2)
        label_x = int(sx + (stat_w - draw.textlength(label, font=F["tiny_bold"])) / 2)
        text(draw, (val_x, stat_y + 12), val, F["stat"], COLORS["text"])
        text(draw, (label_x, stat_y + 112), label, F["tiny_bold"], COLORS["muted"])
        sx += stat_w + 30

    row_y = y0 + 360
    row_h = 112
    row_gap = 124
    for row in data["rows"]:
        rounded(draw, (x0 + 45, row_y, x1 - 45, row_y + row_h), 20, "#151F2D", "#2A394D", 1)
        draw.ellipse((x0 + 72, row_y + 40, x0 + 102, row_y + 70), fill=row["color"])
        text(draw, (x0 + 125, row_y + 22), row["label"], F["small_bold"], COLORS["text"])
        text(draw, (x0 + 125, row_y + 62), f"{row['n_folds']} folds; n_test {row['test_min']}-{row['test_max']}", F["tiny"], COLORS["muted"])

        tile_x = x0 + 435
        for i, mission in enumerate(row["test_missions"][:6]):
            tile_w = max(72, int(draw.textlength(mission, font=F["tiny_bold"]) + 28))
            tx = tile_x
            rounded(draw, (tx, row_y + 29, tx + tile_w, row_y + 82), 16, "#101823", row["color"], 2)
            text(draw, (tx + tile_w / 2, row_y + 44), mission, F["tiny_bold"], COLORS["text"], anchor="ma")
            tile_x += tile_w + 10
        row_y += row_gap

    y = y1 - 285
    rounded(draw, (x0 + 45, y, x1 - 45, y + 232), 28, "#211E17", "#69532B", 2)
    text(draw, (x0 + 78, y + 36), "Readout frame", F["h3"], COLORS["amber"])
    lines = [
        "Evaluation unit: one hidden mission per fold.",
        "Readout: transfer across mission context.",
        "Mechanism and RR-23 / RR-7 evidence appear in companion slides.",
    ]
    multiline(draw, (x0 + 78, y + 90), lines, F["small"], COLORS["text"], 8)


def main() -> None:
    data = load_data()
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(canvas)

    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 48), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 42), width=1)
    draw.rectangle((0, 0, W, 250), fill="#0B1018")

    text(draw, (150, 74), "MISSION-HELD-OUT PROTOCOL", F["kicker"], COLORS["teal"])
    text(draw, (150, 122), "The test set is an entire mission", F["title"])
    text(draw, (150, 216), "Each fold hides one mission so the score asks whether the signal transfers across mission context.", F["subtitle"], COLORS["muted"])
    draw_badges(draw, data)

    draw_random_panel(draw, (150, 350, 910, 1800))
    draw_protocol_panel(draw, (955, 350, 2600, 1800), data["example"])
    draw_inventory_panel(draw, (2645, 350, 3690, 1800), data)

    footer_box = (150, 1880, 3690, 2045)
    rounded(draw, footer_box, 24, "#101823", "#29374A", 2)
    text(draw, (205, 1925), "Takeaway", F["small_bold"], COLORS["blue"])
    source = "The test unit is a whole mission, so the score measures transfer across mission context."
    paragraph(draw, (390, 1925), source, F["small"], 3140, COLORS["muted"], 7)
    text(draw, (205, 1995), "Next", F["small_bold"], COLORS["amber"])
    scope = "After the split is fixed, the next slide checks that train-side choices stop before hidden-mission scoring."
    paragraph(draw, (390, 1995), scope, F["small"], 3140, COLORS["muted"], 7)

    png = OUT_DIR / "mission_heldout_protocol_premium.png"
    canvas.save(png, quality=95)
    gray = OUT_DIR / "mission_heldout_protocol_premium_grayscale.png"
    Image.open(png).convert("L").convert("RGB").save(gray, quality=95)

    manifest = {
        "title": "The test set is an entire mission",
        "source_files": [
            "README.md",
            "tasks/README.md",
            "tasks/A1_liver_lomo/task_info.json",
            "tasks/A2_gastrocnemius_lomo/task_info.json",
            "tasks/A3_kidney_lomo/task_info.json",
            "tasks/A4_thymus_lomo/task_info.json",
            "tasks/A5_skin_lomo/task_info.json",
            "tasks/A6_eye_lomo/task_info.json",
        ],
        "outputs": {
            "png": str(png.relative_to(ROOT)),
            "grayscale": str(gray.relative_to(ROOT)),
        },
        "summary": {
            "core_tissue_tasks": 6,
            "public_lomo_folds": data["total_folds"],
            "unique_mission_labels": len(data["unique_missions"]),
            "example_fold": "A4_thymus_lomo/fold_RR-9_test",
            "example_train_missions": data["example"]["train_missions"],
            "example_hidden_mission": data["example"]["test_mission"],
            "example_n_train": data["example"]["n_train"],
            "example_n_test": data["example"]["n_test"],
        },
        "scope": (
            "Mission-held-out evaluation hides an entire mission per fold. "
            "It controls random-split leakage and establishes benchmark transfer evidence; companion slides carry mechanism evidence."
        ),
    }
    manifest_path = OUT_DIR / "mission_heldout_protocol_premium_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png.relative_to(ROOT)), "grayscale": str(gray.relative_to(ROOT)), "manifest": str(manifest_path.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Build slide 51 asset: manifest, evaluator, and run record."""

from __future__ import annotations

import csv
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
OUT_DIR = ASSET_ROOT / "platform_object_chain"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "manifest_evaluator_run_record_make_scores_rerunnable_premium.png"
GRAY_PATH = OUT_DIR / "manifest_evaluator_run_record_make_scores_rerunnable_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "platform_object_chain_manifest.json"
QA_NOTE = OUT_DIR / "PLATFORM_OBJECT_CHAIN_ASSET_VISUAL_QA.md"

TASK_INDEX = ROOT / "v9" / "task_manifest_index.csv"
TASK_EXAMPLE = ROOT / "v9" / "task_manifests" / "A1_liver_bulk_lomo.json"

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
    "section": load_font(41, True),
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
    "metric": load_font(58, True),
    "mono": load_font(22),
    "mono_bold": load_font(22, True),
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


def arrow(draw: ImageDraw.ImageDraw, x1: int, y1: int, x2: int, y2: int, color: str, width: int = 5) -> None:
    draw.line((x1, y1, x2, y2), fill=color, width=width)
    draw.polygon([(x2, y2), (x2 - 24, y2 - 14), (x2 - 24, y2 + 14)], fill=color)


def down_arrow(draw: ImageDraw.ImageDraw, x: int, y1: int, y2: int, color: str, width: int = 4) -> None:
    draw.line((x, y1, x, y2), fill=color, width=width)
    draw.polygon([(x, y2), (x - 11, y2 - 20), (x + 11, y2 - 20)], fill=color)


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


def load_summary() -> dict[str, object]:
    rows = list(csv.DictReader(TASK_INDEX.open(encoding="utf-8")))
    task = json.loads(TASK_EXAMPLE.read_text(encoding="utf-8"))
    metric_ids = sorted(set(";".join(row["metric_ids"] for row in rows).split(";")))
    return {
        "task_rows": len(rows),
        "folds_total": sum(int(row["n_folds"]) for row in rows),
        "metric_ids": metric_ids,
        "task": task,
    }


def draw_json_lines(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    lines: list[tuple[str, str]],
    max_key_w: int,
    color: str,
) -> int:
    yy = y
    for key, value in lines:
        text(draw, (x, yy), key, F["mono_bold"], color)
        text(draw, (x + max_key_w, yy), value, F["mono"], COLORS["muted"])
        yy += 38
    return yy


def draw_object_panel(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    subtitle: str,
    color: str,
    lines: list[tuple[str, str]],
    output_title: str,
    output_body: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 30, COLORS["panel"], "#2B3A50", 2)
    draw.line((x1 + 32, y1 + 26, x2 - 32, y1 + 26), fill=color, width=5)
    text(draw, (x1 + 42, y1 + 52), title, F["section"], COLORS["text"])
    paragraph(draw, (x1 + 44, y1 + 104), subtitle, F["small"], x2 - x1 - 88, COLORS["muted"], 6)

    code_y = y1 + 210
    rounded(draw, (x1 + 50, code_y, x2 - 50, y2 - 174), 24, COLORS["panel3"], "#2A394D", 1)
    draw_json_lines(draw, x1 + 82, code_y + 34, lines, 230, color)

    rounded(draw, (x1 + 50, y2 - 130, x2 - 50, y2 - 44), 22, "#122234", color, 2)
    text(draw, (x1 + 82, y2 - 102), output_title, F["tiny_bold"], color)
    text(draw, (x1 + 82, y2 - 68), output_body, F["tiny"], COLORS["muted"])


def draw_main_flow(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    task = data["task"]
    metric_ids = data["metric_ids"]
    boxes = {
        "manifest": (120, 610, 1210, 1462),
        "evaluator": (1375, 610, 2465, 1462),
        "record": (2630, 610, 3720, 1462),
    }

    manifest_lines = [
        ("task_id", task["task_id"]),
        ("split", "leave-one-mission-out"),
        ("folds", f"{task['split']['n_folds']} mission folds"),
        ("metrics", "macro_f1, AUROC, mission discrim."),
        ("output", "predictions.csv"),
    ]
    evaluator_lines = [
        ("function", "evaluate_submission()"),
        ("input", "task_manifest + predictions.csv"),
        ("validates", "sample_id, labels, scores"),
        ("computes", ", ".join(metric_ids[:3])),
        ("status", "evaluated / invalid"),
    ]
    record_lines = [
        ("metrics", "metrics.json"),
        ("run", "run_manifest.json"),
        ("stores", "task_id, command, status"),
        ("inputs", "fingerprinted files"),
        ("checks", "metrics + validation"),
    ]

    draw_object_panel(
        draw,
        boxes["manifest"],
        "Task Manifest",
        "The task contract fixes the split, metric set, and prediction format before scoring.",
        COLORS["blue"],
        manifest_lines,
        "fixed contract",
        "one task identity, one split, one metric set",
    )
    draw_object_panel(
        draw,
        boxes["evaluator"],
        "Evaluator",
        "The scoring function validates a submission and computes the declared metrics.",
        COLORS["teal"],
        evaluator_lines,
        "scoring surface",
        "same task plus same predictions gives the same report shape",
    )
    draw_object_panel(
        draw,
        boxes["record"],
        "Run Record",
        "The report files preserve what ran and what scores were produced.",
        COLORS["amber"],
        record_lines,
        "repeatable report",
        "metrics and execution details travel with the result",
    )

    mid_y = 1036
    arrow(draw, boxes["manifest"][2] + 35, mid_y, boxes["evaluator"][0] - 35, mid_y, COLORS["dim"], 5)
    arrow(draw, boxes["evaluator"][2] + 35, mid_y, boxes["record"][0] - 35, mid_y, COLORS["dim"], 5)


def draw_bottom_contract(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1530, 3720, 1848
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Rerunnable Score Path", F["h2"], COLORS["text"])
    text(draw, (x2 - 80, y1 + 48), "manifest + evaluator + run record", F["small_bold"], COLORS["teal"], "ra")

    steps = [
        ("task_manifest.json", "split + metrics", COLORS["blue"]),
        ("predictions.csv", "submitted scores", COLORS["violet"]),
        ("evaluate_submission()", "validation + metrics", COLORS["teal"]),
        ("metrics.json", "score report", COLORS["green"]),
        ("run_manifest.json", "run record", COLORS["amber"]),
    ]
    node_w, gap = 575, 70
    start_x, node_y = x1 + 190, y1 + 150
    for i, (title, desc, color) in enumerate(steps):
        nx = start_x + i * (node_w + gap)
        rounded(draw, (nx, node_y, nx + node_w, node_y + 96), 20, COLORS["panel2"], color, 2)
        text(draw, (nx + 26, node_y + 17), title, F["small_bold"], COLORS["text"])
        text(draw, (nx + 26, node_y + 55), desc, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            arrow(draw, nx + node_w + 10, node_y + 48, nx + node_w + gap - 16, node_y + 48, COLORS["dim"], 4)


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "A score is rerunnable when the task contract, evaluator, metrics report, and run record travel together.",
        F["small"],
        3180,
        COLORS["muted"],
        8,
    )


def build() -> None:
    data = load_summary()
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 51 | ACT 6 | OBJECT CHAIN", F["kicker"], COLORS["teal"])
    bx = 2050
    bx = badge(draw, bx, 56, "TASKS", str(data["task_rows"]), COLORS["blue"])
    bx = badge(draw, bx, 56, "FOLDS", str(data["folds_total"]), COLORS["teal"])
    bx = badge(draw, bx, 56, "METRICS", str(len(data["metric_ids"])), COLORS["green"])
    badge(draw, bx, 56, "RUN FILES", "2", COLORS["amber"])

    text(draw, (120, 330), "Manifest, Evaluator, And Run Record Make Scores Rerunnable", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "Each benchmark score is tied to a task contract, a scoring function, and report files that preserve the run.",
        F["subtitle"],
        3180,
        COLORS["muted"],
        8,
    )

    draw_main_flow(draw, data)
    draw_bottom_contract(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    ImageOps.colorize(gray, black="#080D14", white="#F3F7FC").save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(gray)

    manifest = {
        "title": "Manifest, Evaluator, And Run Record Make Scores Rerunnable",
        "claim": "The platform object chain ties every score to a task contract, evaluator, metrics report, and run record.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "data_inputs": {
            "task_index": str(TASK_INDEX.relative_to(ROOT)),
            "task_example": str(TASK_EXAMPLE.relative_to(ROOT)),
            "run_record_code": "spacebio_bench/reports/run_manifest.py",
            "evaluator_code": "spacebio_bench/evaluate.py",
        },
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": {
            "task_rows": data["task_rows"],
            "folds_total": data["folds_total"],
            "metric_count": len(data["metric_ids"]),
            "metric_ids": data["metric_ids"],
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    QA_NOTE.write_text(
        "\n".join(
            [
                "# Platform Object Chain Asset Visual QA",
                "",
                "Slide 51 explains the three core objects that make benchmark scores rerunnable.",
                "",
                "Checks to review:",
                "- Task Manifest, Evaluator, and Run Record panels read independently.",
                "- Arrows stay between panels and do not cross text.",
                "- JSON-like snippets remain legible at presentation scale.",
                "- Bottom path previews the file sequence without crowding.",
                "- Grayscale version preserves panel hierarchy.",
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

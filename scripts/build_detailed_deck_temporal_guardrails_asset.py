#!/usr/bin/env python3
"""Build slide 34 asset: temporal labels need context."""

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
OUT_DIR = ASSET_ROOT / "temporal_guardrails"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "temporal_labels_need_guardrails_premium.png"
GRAY_PATH = OUT_DIR / "temporal_labels_need_guardrails_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "temporal_guardrails_manifest.json"

TEMPORAL_SUMMARY = ROOT / "v2" / "evaluation" / "T_temporal_summary.json"
SOURCE_FIGURE = ROOT / "v2" / "figures" / "Fig1_temporal.html"

COLORS = {
    "bg": "#0C111A",
    "bg2": "#091019",
    "panel": "#101823",
    "panel2": "#151F2D",
    "panel3": "#0F1A26",
    "grid": "#263244",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#667487",
    "blue": "#66A6E8",
    "teal": "#5FD3C4",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "orange": "#E69F00",
    "violet": "#B39DFF",
    "rose": "#E17882",
    "ink": "#080D14",
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
    "title": load_font(84, True),
    "subtitle": load_font(37),
    "h1": load_font(50, True),
    "h2": load_font(40, True),
    "h3": load_font(32, True),
    "body": load_font(29),
    "body_bold": load_font(29, True),
    "small": load_font(25),
    "small_bold": load_font(25, True),
    "tiny": load_font(21),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "stat": load_font(66, True),
    "mono": load_font(24),
    "mono_bold": load_font(24, True),
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


def wrap_lines(draw: ImageDraw.ImageDraw, label: str, font: ImageFont.ImageFont, max_width: int) -> list[str]:
    words = label.split()
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


def load_data() -> dict[str, object]:
    temporal = json.loads(TEMPORAL_SUMMARY.read_text())
    rr8 = temporal["T1"]["T1b_RR-8_liver"]
    rr6 = temporal["T1"]["T1a_RR-6_liver"]
    rr8_conditions = rr8["conditions"]
    rr6_conditions = rr6["conditions"]
    t2 = temporal["T2"]
    t3 = temporal["T3"]["aging_hypothesis"]
    old = temporal["T3"]["subtasks"]["T3d_OLD_gene"]
    yng = temporal["T3"]["subtasks"]["T3d_YNG_gene"]
    return {
        "preservation": {
            "rr8": {
                "FLT": rr8_conditions["FLT_gene"]["auroc"],
                "GC": rr8_conditions["GC_gene"]["auroc"],
                "BSL": rr8_conditions["BSL_gene"]["auroc"],
                "excess": rr8["interpretation"]["excess_signal"],
                "n_flt": rr8_conditions["FLT_gene"]["n_total"],
                "n_gc": rr8_conditions["GC_gene"]["n_total"],
                "n_bsl": rr8_conditions["BSL_gene"]["n_total"],
            },
            "rr6": {
                "FLT": rr6_conditions["FLT_gene"]["auroc"],
                "GC": rr6_conditions["GC_gene"]["auroc"],
                "BSL": rr6_conditions["BSL_gene"]["auroc"],
                "excess": rr6["interpretation"]["excess_signal"],
            },
        },
        "recovery": {
            "RR-6": {
                "pca_ratio": t2["T2_RR-6"]["pca_recovery_ratio"],
                "n_recovering": t2["T2_RR-6"]["pathway_recovery_summary"]["n_recovering"],
                "n_total": t2["T2_RR-6"]["pathway_recovery_summary"]["n_total"],
                "iss_t_prob": t2["T2_RR-6"]["classification_shift"]["mean_flt_isst_flight_prob"],
                "lar_prob": t2["T2_RR-6"]["classification_shift"]["mean_flt_lar_flight_prob"],
            },
            "RR-8": {
                "pca_ratio": t2["T2_RR-8"]["pca_recovery_ratio"],
                "n_recovering": t2["T2_RR-8"]["pathway_recovery_summary"]["n_recovering"],
                "n_total": t2["T2_RR-8"]["pathway_recovery_summary"]["n_total"],
                "iss_t_prob": t2["T2_RR-8"]["classification_shift"]["mean_flt_isst_flight_prob"],
                "lar_prob": t2["T2_RR-8"]["classification_shift"]["mean_flt_lar_flight_prob"],
            },
        },
        "age": {
            "old_auroc": old["auroc"],
            "yng_auroc": yng["auroc"],
            "delta": t3["delta"],
            "old_n": old["n_flight"] + old["n_gc"],
            "yng_n": yng["n_flight"] + yng["n_gc"],
        },
    }


def draw_background(draw: ImageDraw.ImageDraw) -> None:
    draw.rectangle((0, 0, W, H), fill=COLORS["bg"])
    for y in range(H):
        draw.line((0, y, W, y), fill=blend(COLORS["bg"], COLORS["bg2"], y / H))
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=COLORS["grid"], width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill="#172234", width=1)
    draw.rectangle((0, 0, W, 315), fill="#0F1824")
    draw.rectangle((0, 1840, W, H), fill=COLORS["ink"])
    draw.line((0, 315, W, 315), fill="#1E2B3D", width=2)
    draw.line((0, 1840, W, 1840), fill="#1E2B3D", width=2)


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(205, int(draw.textlength(value, font=F["tiny_bold"]) + 76))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def draw_header(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 76), "BIOLOGY / TEMPORAL CONTEXT", F["kicker"], COLORS["teal"])
    x = 2180
    x = badge(draw, x, 66, "dataset", "v2 temporal", COLORS["blue"])
    x = badge(draw, x, 66, "tissue", "liver", COLORS["green"])
    x = badge(draw, x, 66, "missions", "RR-6 / RR-8", COLORS["amber"])
    badge(draw, x, 66, "layers", "T1-T3", COLORS["violet"])


def draw_title(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 382), "Temporal Labels Need Context", F["title"], COLORS["text"])
    paragraph(
        draw,
        (155, 493),
        "ISS-T and LAR labels carry preservation, return timing, and age context. Separating the layers makes timepoint biology easier to read.",
        F["subtitle"],
        2050,
        COLORS["muted"],
        10,
    )


def draw_bar(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    label: str,
    value: float,
    color: str,
    max_w: int = 620,
    height: int = 44,
) -> None:
    text(draw, (x, y + 8), label, F["small_bold"], COLORS["muted"])
    bar_x = x + 115
    rounded(draw, (bar_x, y, bar_x + max_w, y + height), 12, "#0C131D", "#25364B", 1)
    rounded(draw, (bar_x + 4, y + 4, bar_x + 4 + int((max_w - 8) * value), y + height - 4), 9, rgba(color, 215))
    text(draw, (bar_x + max_w + 26, y + 8), f"{value:.3f}", F["small_bold"], COLORS["text"])


def draw_panel_title(draw: ImageDraw.ImageDraw, x: int, y: int, number: str, title: str, color: str) -> None:
    rounded(draw, (x, y, x + 66, y + 66), 18, rgba(color, 46), color, 2)
    text(draw, (x + 33, y + 18), number, F["h3"], color, "ma")
    text(draw, (x + 90, y + 10), title, F["h2"], COLORS["text"])


def draw_preservation_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x0, y0, x1, y1 = 150, 640, 1230, 1585
    rounded(draw, (x0, y0, x1, y1), 24, COLORS["panel"], "#2A394D", 2)
    draw_panel_title(draw, x0 + 42, y0 + 40, "1", "Preservation layer", COLORS["rose"])
    paragraph(
        draw,
        (x0 + 132, y0 + 96),
        "Control rows separate strongly by ISS-T versus LAR timing, so the label readout includes handling signal.",
        F["body"],
        805,
        COLORS["muted"],
        8,
    )

    rr8 = data["preservation"]["rr8"]  # type: ignore[index]
    text(draw, (x0 + 54, y0 + 230), "RR-8 liver AUROC by row type", F["small_bold"], COLORS["text"])
    draw_bar(draw, x0 + 62, y0 + 292, "FLT", float(rr8["FLT"]), COLORS["blue"])
    draw_bar(draw, x0 + 62, y0 + 362, "GC", float(rr8["GC"]), COLORS["rose"])
    draw_bar(draw, x0 + 62, y0 + 432, "BSL", float(rr8["BSL"]), COLORS["amber"])

    rounded(draw, (x0 + 58, y0 + 540, x1 - 58, y0 + 706), 20, "#111B28", COLORS["rose"], 2)
    text(draw, (x0 + 86, y0 + 572), "Readout", F["tiny_bold"], COLORS["rose"])
    text(draw, (x0 + 86, y0 + 615), "GC 0.973 > FLT 0.930", F["h3"], COLORS["text"])
    paragraph(
        draw,
        (x0 + 86, y0 + 660),
        "The preservation/handling axis is visible in controls.",
        F["small"],
        800,
        COLORS["muted"],
        6,
    )

    decoder = "ISS-T = terminal collection on ISS  |  LAR = live animal return  |  BSL = baseline"
    rounded(draw, (x0 + 58, y1 - 120, x1 - 58, y1 - 44), 18, "#0D1520", "#26364B", 2)
    text(draw, (x0 + 86, y1 - 93), decoder, F["tiny_bold"], COLORS["muted"])


def draw_recovery_row(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    mission: str,
    row: dict[str, float | int],
    color: str,
) -> None:
    text(draw, (x, y + 14), mission, F["h3"], color)
    rail_x0, rail_x1 = x + 160, x + 690
    rail_y = y + 34
    draw.line((rail_x0, rail_y, rail_x1, rail_y), fill="#3A4B62", width=5)
    dots = [(rail_x0, "ISS-T", COLORS["rose"]), ((rail_x0 + rail_x1) // 2, "LAR", color), (rail_x1, "baseline", COLORS["green"])]
    for dx, label, dot_color in dots:
        draw.ellipse((dx - 18, rail_y - 18, dx + 18, rail_y + 18), fill=dot_color, outline=dot_color, width=2)
        text(draw, (dx, rail_y + 34), label, F["micro_bold"], COLORS["dim"], "ma")
    text(draw, (x + 160, y + 106), "PCA ratio", F["tiny_bold"], COLORS["dim"])
    text(draw, (x + 285, y + 92), f"{float(row['pca_ratio']):.3f}", F["h2"], color)
    text(draw, (x + 485, y + 106), "pathways", F["tiny_bold"], COLORS["dim"])
    text(draw, (x + 610, y + 94), f"{int(row['n_recovering'])}/{int(row['n_total'])}", F["h3"], COLORS["text"])
    text(draw, (x + 160, y + 156), "probability shift", F["tiny_bold"], COLORS["dim"])
    text(draw, (x + 350, y + 145), f"{float(row['iss_t_prob']):.3f} -> {float(row['lar_prob']):.3f}", F["small_bold"], COLORS["muted"])


def draw_recovery_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x0, y0, x1, y1 = 1380, 640, 2460, 1585
    rounded(draw, (x0, y0, x1, y1), 24, COLORS["panel"], "#2A394D", 2)
    draw_panel_title(draw, x0 + 42, y0 + 40, "2", "Recovery projection", COLORS["amber"])
    paragraph(
        draw,
        (x0 + 132, y0 + 96),
        "LAR profiles shift toward baseline along a projection axis. The direction and pathway fraction differ by mission.",
        F["body"],
        805,
        COLORS["muted"],
        8,
    )

    recovery = data["recovery"]  # type: ignore[index]
    draw_recovery_row(draw, x0 + 64, y0 + 260, "RR-6", recovery["RR-6"], COLORS["teal"])  # type: ignore[index]
    draw_recovery_row(draw, x0 + 64, y0 + 525, "RR-8", recovery["RR-8"], COLORS["amber"])  # type: ignore[index]

    rounded(draw, (x0 + 58, y1 - 130, x1 - 58, y1 - 44), 18, "#111B28", COLORS["amber"], 2)
    text(draw, (x0 + 86, y1 - 102), "Use as descriptive recovery signal attached to mission context.", F["small_bold"], COLORS["text"])
    text(draw, (x0 + 86, y1 - 68), "PCA ratio and pathway count travel together.", F["tiny_bold"], COLORS["muted"])


def draw_age_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    x0, y0, x1, y1 = 2610, 640, 3690, 1585
    rounded(draw, (x0, y0, x1, y1), 24, COLORS["panel"], "#2A394D", 2)
    draw_panel_title(draw, x0 + 42, y0 + 40, "3", "Age / context layer", COLORS["violet"])
    paragraph(
        draw,
        (x0 + 132, y0 + 96),
        "Age group changes response strength in RR-8 liver, adding a second context field to temporal interpretation.",
        F["body"],
        805,
        COLORS["muted"],
        8,
    )

    age = data["age"]  # type: ignore[index]
    chart_x, chart_y = x0 + 78, y0 + 280
    draw_bar(draw, chart_x, chart_y, "OLD", float(age["old_auroc"]), COLORS["violet"], 610, 52)
    draw_bar(draw, chart_x, chart_y + 92, "YNG", float(age["yng_auroc"]), COLORS["blue"], 610, 52)

    rounded(draw, (x0 + 78, y0 + 505, x1 - 78, y0 + 760), 22, "#111B28", COLORS["violet"], 2)
    text(draw, (x0 + 112, y0 + 540), "RR-8 liver spaceflight detection", F["small_bold"], COLORS["muted"])
    text(draw, (x0 + 112, y0 + 588), f"OLD - YNG delta {float(age['delta']):.3f}", F["h2"], COLORS["text"])
    text(draw, (x0 + 112, y0 + 646), f"OLD n={int(age['old_n'])} | YNG n={int(age['yng_n'])}", F["small_bold"], COLORS["muted"])

    for i, label in enumerate(["preservation", "return timing", "age context"]):
        px = x0 + 92 + i * 292
        rounded(draw, (px, y1 - 134, px + 248, y1 - 48), 18, "#0D1520", [COLORS["rose"], COLORS["amber"], COLORS["violet"]][i], 2)
        text(draw, (px + 124, y1 - 101), label, F["tiny_bold"], COLORS["text"], "ma")


def draw_reader_rule(draw: ImageDraw.ImageDraw) -> None:
    y0, y1 = 1625, 1816
    rounded(draw, (150, y0, 3690, y1), 24, "#0E1722", "#243247", 2)
    text(draw, (190, y0 + 36), "Reading path", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (190, y0 + 94),
        "Temporal labels become context fields for biology interpretation. Separate collection handling, return timing, and age/context first.",
        F["body"],
        940,
        COLORS["muted"],
        6,
    )
    steps = [
        ("1", "preservation", "controls establish handling signal", COLORS["rose"]),
        ("2", "recovery", "LAR-to-baseline projection", COLORS["amber"]),
        ("3", "context", "age group changes response strength", COLORS["violet"]),
    ]
    x = 1290
    for number, title, body, color in steps:
        rounded(draw, (x, y0 + 40, x + 710, y1 - 38), 22, "#111C2A", color, 2)
        rounded(draw, (x + 30, y0 + 68, x + 88, y0 + 126), 16, rgba(color, 50), color, 2)
        text(draw, (x + 59, y0 + 82), number, F["h3"], color, "ma")
        text(draw, (x + 112, y0 + 58), title, F["small_bold"], COLORS["text"])
        paragraph(draw, (x + 112, y0 + 96), body, F["small"], 520, COLORS["muted"], 6)
        x += 750


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 1902), "TAKEAWAY", F["micro_bold"], COLORS["dim"])
    paragraph(
        draw,
        (150, 1934),
        "Preservation, return timing, and age context should travel with temporal labels before biology interpretation.",
        F["small"],
        1550,
        COLORS["muted"],
        6,
    )
    text(draw, (2170, 1902), "NEXT", F["micro_bold"], COLORS["dim"])
    paragraph(
        draw,
        (2170, 1934),
        "Single-cell pilots, spatial weak-signal cases, and systems biology add the next interpretation layers.",
        F["small"],
        1480,
        COLORS["muted"],
        6,
    )


def main() -> None:
    data = load_data()
    img = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(img, "RGBA")
    draw_background(draw)
    draw_header(draw)
    draw_title(draw)
    draw_preservation_panel(draw, data)
    draw_recovery_panel(draw, data)
    draw_age_panel(draw, data)
    draw_reader_rule(draw)
    draw_footer(draw)

    img.save(OUT_PATH, quality=95)
    gray = ImageOps.grayscale(img).convert("RGB")
    gray.save(GRAY_PATH, quality=95)
    stat = ImageStat.Stat(gray)
    manifest = {
        "asset": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "source_files": [str(TEMPORAL_SUMMARY.relative_to(ROOT)), str(SOURCE_FIGURE.relative_to(ROOT))],
        "metrics": {
            "rr8_preservation_auroc": data["preservation"]["rr8"],
            "recovery": data["recovery"],
            "age": data["age"],
            "grayscale_mean": stat.mean,
            "grayscale_stddev": stat.stddev,
        },
        "takeaway": "Temporal context slide for preservation, recovery projection, and age/context labels.",
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"asset": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

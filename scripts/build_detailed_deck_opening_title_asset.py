#!/usr/bin/env python3
"""Build slide 1 opening title asset for the detailed SpaceBio-Bench deck."""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFilter, ImageFont, ImageOps, ImageStat


ROOT = Path(__file__).resolve().parent.parent
ASSET_ROOT = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
)
OUT_DIR = ASSET_ROOT / "opening_title"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "spacebiobench_opening_title_premium.png"
GRAY_PATH = OUT_DIR / "spacebiobench_opening_title_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "spacebiobench_opening_title_manifest.json"

EVIDENCE_STACK = ASSET_ROOT / "evidence_stack" / "evidence_stack_manifest.json"
METHOD_HARDENING = ASSET_ROOT / "method_hardening" / "method_hardening_manifest.json"

COLORS = {
    "bg": "#0C111A",
    "bg2": "#091019",
    "panel": "#101823",
    "panel2": "#151F2D",
    "panel3": "#1A2534",
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
    "red": "#F06E7D",
    "white": "#FFFFFF",
}


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


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
    "title": load_font(142, True),
    "title_small": load_font(72, True),
    "subtitle": load_font(42),
    "h1": load_font(54, True),
    "h2": load_font(42, True),
    "h3": load_font(32, True),
    "body": load_font(30),
    "body_bold": load_font(30, True),
    "small": load_font(25),
    "small_bold": load_font(25, True),
    "tiny": load_font(21),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "stat": load_font(78, True),
    "chip": load_font(28, True),
}


def rounded(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    radius: int,
    fill: str,
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
    leading: int = 9,
) -> int:
    x, y = xy
    for block in body.splitlines() or [""]:
        if not block:
            y += font.size + leading
            continue
        for line in wrap_lines(draw, block, font, max_width):
            text(draw, (x, y), line, font, fill)
            y += font.size + leading
    return y


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


def load_metrics() -> dict[str, object]:
    evidence = load_json(EVIDENCE_STACK)["metrics"]
    method = load_json(METHOD_HARDENING)["metrics"]
    return {
        "tissues": evidence["mission_tasks"],
        "lomo_folds": evidence["lomo_folds"],
        "mission_labels": evidence["mission_labels"],
        "model_configs": method["n_configs"],
        "significant_rows": method["n_significant"],
        "control_rows": evidence["control_rows"],
        "dge_rho": evidence["dge_log2fc_spearman_mean"],
    }


def draw_background(draw: ImageDraw.ImageDraw) -> None:
    draw.rectangle((0, 0, W, H), fill=COLORS["bg"])
    for y in range(H):
        t = y / H
        color = blend(COLORS["bg"], COLORS["bg2"], t)
        draw.line((0, y, W, y), fill=color)

    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=COLORS["grid"], width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill="#172234", width=1)

    draw.rectangle((0, 0, W, 330), fill="#0F1824")
    draw.rectangle((0, 1840, W, H), fill="#080D14")
    draw.line((0, 330, W, 330), fill="#1E2B3D", width=2)
    draw.line((0, 1840, W, 1840), fill="#1E2B3D", width=2)


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(210, int(draw.textlength(value, font=F["tiny_bold"]) + 72))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def draw_header(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    text(draw, (150, 82), "DETAILED TECHNICAL WALKTHROUGH", F["kicker"], COLORS["teal"])
    x = 2355
    x = badge(draw, x, 72, "tasks", f"{data['tissues']} tissues", COLORS["teal"])
    x = badge(draw, x, 72, "folds", f"{data['lomo_folds']} LOMO", COLORS["blue"])
    badge(draw, x, 72, "configs", f"{data['model_configs']} configs", COLORS["amber"])


def draw_title_block(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 410), "SpaceBio-Bench", F["title"], COLORS["text"])
    paragraph(
        draw,
        (160, 580),
        "Mission-held-out evaluation for spaceflight biology under domain shift.",
        F["subtitle"],
        1260,
        COLORS["muted"],
        11,
    )
    rounded(draw, (160, 745, 1385, 1002), 26, "#101A27", "#2D3D53", 2)
    text(draw, (205, 792), "Thesis", F["tiny_bold"], COLORS["teal"])
    paragraph(
        draw,
        (205, 835),
        "Benchmarking spaceflight biology requires evidence that survives a hidden mission.",
        F["h1"],
        1085,
        COLORS["text"],
        10,
    )
    rounded(draw, (160, 1040, 1385, 1245), 22, blend("#101823", COLORS["amber"], 0.11), COLORS["amber"], 2)
    text(draw, (205, 1082), "Reader question", F["tiny_bold"], COLORS["amber"])
    paragraph(
        draw,
        (205, 1125),
        "Can a model trained on known missions recognize a biology signal in a new mission?",
        F["body_bold"],
        1085,
        COLORS["text"],
        8,
    )


def draw_lane_chip(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], label: str, sub: str, color: str) -> None:
    rounded(draw, box, 18, "#142033", color, 2)
    x1, y1, _, _ = box
    text(draw, (x1 + 24, y1 + 20), label, F["chip"], COLORS["text"])
    text(draw, (x1 + 24, y1 + 58), sub, F["tiny"], COLORS["muted"])


def arrow(draw: ImageDraw.ImageDraw, x1: int, y1: int, x2: int, y2: int, color: str) -> None:
    draw.line((x1, y1, x2, y2), fill=color, width=5)
    draw.polygon([(x2, y2), (x2 - 22, y2 - 13), (x2 - 22, y2 + 13)], fill=color)


def draw_benchmark_contract(draw: ImageDraw.ImageDraw) -> None:
    panel = (1500, 430, 3685, 1298)
    rounded(draw, panel, 34, COLORS["panel"], "#2A394D", 2)
    text(draw, (1550, 478), "Benchmark contract", F["h2"], COLORS["text"])
    text(draw, (1550, 531), "Same biological question, scored only after the mission is hidden.", F["small"], COLORS["muted"])

    y = 700
    source = (1570, y, 1995, y + 180)
    task = (2150, y, 2575, y + 180)
    score = (2730, y, 3155, y + 180)
    claim = (3310, y, 3635, y + 180)
    draw_lane_chip(draw, source, "Public records", "OSDR study tables", COLORS["blue"])
    draw_lane_chip(draw, task, "Fixed task", "input + label + split", COLORS["teal"])
    draw_lane_chip(draw, score, "Hidden mission", "scored once", COLORS["amber"])
    draw_lane_chip(draw, claim, "Readout status", "context attached", COLORS["green"])
    arrow(draw, 2015, y + 90, 2130, y + 90, COLORS["teal"])
    arrow(draw, 2595, y + 90, 2710, y + 90, COLORS["amber"])
    arrow(draw, 3175, y + 90, 3290, y + 90, COLORS["green"])

    lane_y = 980
    rounded(draw, (1570, lane_y, 3635, lane_y + 250), 26, "#0F1722", "#24364D", 2)
    text(draw, (1615, lane_y + 34), "What has to stay visible throughout the deck", F["h3"], COLORS["text"])
    rows = [
        ("Split", "the test unit is a whole mission", COLORS["teal"]),
        ("Input", "genes, pathways, embeddings, or text prompts", COLORS["violet"]),
        ("Metric", "AUROC plus uncertainty defines the readout", COLORS["blue"]),
        ("Layering", "benchmark readout stays separate from mechanism and translation", COLORS["amber"]),
    ]
    x = 1615
    for label, body, color in rows:
        rounded(draw, (x, lane_y + 95, x + 470, lane_y + 218), 18, blend("#121C29", color, 0.10), color, 2)
        text(draw, (x + 24, lane_y + 130), label, F["tiny_bold"], color)
        paragraph(draw, (x + 24, lane_y + 165), body, F["tiny"], 405, COLORS["muted"], 3)
        x += 500


def draw_proof_tiles(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    tiles = [
        ("1", "Method first", f"{data['lomo_folds']} mission-held-out folds; train choices stop before scoring.", COLORS["teal"]),
        ("2", "Comparable model surfaces", f"{data['model_configs']} tested configurations under one task contract.", COLORS["violet"]),
        ("3", "Scores become readout status", f"{data['control_rows']} controls plus DGE and biology checks guide interpretation.", COLORS["amber"]),
    ]
    x = 160
    for idx, title, body, color in tiles:
        rounded(draw, (x, 1390, x + 1225, 1748), 28, COLORS["panel"], "#2A394D", 2)
        rounded(draw, (x + 35, 1428, x + 105, 1498), 20, blend("#101823", color, 0.22), color, 2)
        text(draw, (x + 70, 1464), idx, F["h3"], COLORS["text"], anchor="mm")
        text(draw, (x + 135, 1433), title, F["h2"], COLORS["text"])
        paragraph(draw, (x + 135, 1492), body, F["body"], 995, COLORS["muted"], 8)
        x += 1270


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (150, 1886, 3690, 2068), 22, "#111B28", "#26384E", 2)
    text(draw, (190, 1928), "Takeaway", F["tiny_bold"], COLORS["blue"])
    paragraph(
        draw,
        (360, 1928),
        "SpaceBio-Bench asks whether biology signals survive an unseen mission, not just a familiar study split.",
        F["small"],
        2570,
        COLORS["muted"],
        7,
    )
    text(draw, (190, 1992), "Next", F["tiny_bold"], COLORS["amber"])
    paragraph(
        draw,
        (360, 1992),
        "The opening section turns study records into fixed tasks, hidden-mission folds, and interpretable score rows.",
        F["small"],
        2750,
        COLORS["muted"],
        7,
    )
    text(draw, (3510, 1998), "1", F["h2"], COLORS["teal"])


def render() -> Image.Image:
    data = load_metrics()
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image)
    draw_background(draw)
    draw_header(draw, data)
    draw_title_block(draw)
    draw_benchmark_contract(draw)
    draw_proof_tiles(draw, data)
    draw_footer(draw)
    return image


def image_metrics(path: Path) -> dict[str, object]:
    image = Image.open(path).convert("RGB")
    gray = ImageOps.grayscale(image)
    stat = ImageStat.Stat(gray)
    edges = gray.filter(ImageFilter.FIND_EDGES)
    edge_stat = ImageStat.Stat(edges)
    return {
        "width_px": image.width,
        "height_px": image.height,
        "luma_mean": round(stat.mean[0], 2),
        "luma_stddev": round(stat.stddev[0], 2),
        "edge_mean": round(edge_stat.mean[0], 2),
        "nonblank_pass": stat.stddev[0] >= 16.0 and edge_stat.mean[0] >= 2.0,
    }


def main() -> None:
    image = render()
    image.save(OUT_PATH, quality=95)
    gray = ImageOps.grayscale(image).convert("RGB")
    gray.save(GRAY_PATH, quality=95)
    metrics = load_metrics()
    manifest = {
        "title": "SpaceBio-Bench",
        "slide": 1,
        "status": "ready",
        "role": "opening thesis title scene",
        "created": "2026-06-11",
        "outputs": {
            "png": str(OUT_PATH.relative_to(ROOT)),
            "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        },
        "metrics": metrics,
        "copy": {
            "headline": "SpaceBio-Bench",
            "subtitle": "Mission-held-out evaluation for spaceflight biology under domain shift.",
            "thesis": "Benchmarking spaceflight biology requires evidence that survives a hidden mission.",
            "reader_question": "Can a model trained on known missions recognize a biology signal in a new mission?",
            "readout_frame": "Benchmark transfer readout under tested task contracts; mechanism and translation appear as companion layers.",
        },
        "source_manifests": [
            str(EVIDENCE_STACK.relative_to(ROOT)),
            str(METHOD_HARDENING.relative_to(ROOT)),
        ],
        "automatic_metrics": image_metrics(OUT_PATH),
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"png": str(OUT_PATH), "grayscale": str(GRAY_PATH), "manifest": str(MANIFEST_PATH)}, indent=2))


if __name__ == "__main__":
    main()

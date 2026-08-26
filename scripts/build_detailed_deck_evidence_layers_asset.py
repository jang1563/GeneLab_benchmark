#!/usr/bin/env python3
"""Build slide 3 asset: evidence layers organize the story."""

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
OUT_DIR = ASSET_ROOT / "evidence_layers"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "evidence_layers_set_claim_strength_premium.png"
GRAY_PATH = OUT_DIR / "evidence_layers_set_claim_strength_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "evidence_layers_manifest.json"

EVIDENCE_STACK = ASSET_ROOT / "evidence_stack" / "evidence_stack_manifest.json"
METHOD_HARDENING = ASSET_ROOT / "method_hardening" / "method_hardening_manifest.json"
EXTERNAL_VALIDATION = ASSET_ROOT / "external_validation" / "external_biology_validation_premium_manifest.json"

COLORS = {
    "bg": "#0C111A",
    "bg2": "#091019",
    "panel": "#101823",
    "panel2": "#151F2D",
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
    "stat": load_font(64, True),
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
    leading: int = 8,
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


def load_metrics() -> dict[str, object]:
    evidence = json.loads(EVIDENCE_STACK.read_text())["metrics"]
    method = json.loads(METHOD_HARDENING.read_text())["metrics"]
    external = json.loads(EXTERNAL_VALIDATION.read_text())["summary"]
    return {
        "tissues": evidence["mission_tasks"],
        "lomo_folds": evidence["lomo_folds"],
        "mission_labels": evidence["mission_labels"],
        "model_configs": method["n_configs"],
        "control_rows": evidence["control_rows"],
        "dge_rho": evidence["dge_log2fc_spearman_mean"],
        "pathway_concordance": external["pathway_mean_concordance"],
        "gene_enrichment": external["gene_enrichment_ratio"],
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
    draw.rectangle((0, 1840, W, H), fill="#080D14")
    draw.line((0, 315, W, 315), fill="#1E2B3D", width=2)
    draw.line((0, 1840, W, 1840), fill="#1E2B3D", width=2)


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(205, int(draw.textlength(value, font=F["tiny_bold"]) + 72))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def draw_header(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    text(draw, (150, 76), "OPENING / STORY LAYERS", F["kicker"], COLORS["teal"])
    x = 2265
    x = badge(draw, x, 66, "benchmark", f"{data['lomo_folds']} folds", COLORS["teal"])
    x = badge(draw, x, 66, "controls", f"{data['control_rows']} rows", COLORS["amber"])
    badge(draw, x, 66, "methods", f"{data['model_configs']} configs", COLORS["violet"])


def draw_title(draw: ImageDraw.ImageDraw) -> None:
    text(draw, (150, 382), "Evidence Layers Organize The Story", F["title"], COLORS["text"])
    paragraph(
        draw,
        (155, 490),
        "The same project carries benchmark results, biological support, translation hypotheses, and platform readiness as separate evidence layers.",
        F["subtitle"],
        1740,
        COLORS["muted"],
        10,
    )


def draw_stack_row(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    number: str,
    label: str,
    question: str,
    body: str,
    color: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 22, blend("#101823", color, 0.10), color, 2)
    rounded(draw, (x1 + 22, y1 + 26, x1 + 88, y1 + 92), 18, blend("#101823", color, 0.22), color, 2)
    text(draw, (x1 + 55, y1 + 60), number, F["h3"], COLORS["text"], anchor="mm")
    text(draw, (x1 + 116, y1 + 20), label, F["tiny_bold"], color)
    text(draw, (x1 + 116, y1 + 54), question, F["h3"], COLORS["text"])
    text(draw, (x1 + 116, y1 + 96), body, F["micro"], COLORS["muted"])


def draw_claim_stack(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    panel = (150, 680, 2110, 1620)
    rounded(draw, panel, 30, COLORS["panel"], "#2A394D", 2)
    text(draw, (200, 730), "A. Evidence ladder", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (200, 785),
        "Each layer adds a different kind of confidence to the readout.",
        F["small"],
        1640,
        COLORS["muted"],
        7,
    )
    rows = [
        (
            "1",
            "Data layer",
            "What is available?",
            "OSDR records, metadata, expression matrices, DGE outputs.",
            COLORS["blue"],
        ),
        (
            "2",
            "Task layer",
            "What is scored?",
            "Input, label, split, metric, and readout are fixed.",
            COLORS["violet"],
        ),
        (
            "3",
            "Benchmark readout",
            "What transfers?",
            f"{data['tissues']} tissues, {data['lomo_folds']} hidden-mission folds, AUROC plus uncertainty.",
            COLORS["teal"],
        ),
        (
            "4",
            "Robustness layer",
            "How stable is the readout?",
            f"{data['control_rows']} controls, DGE rho {data['dge_rho']:.3f}, method grid checks.",
            COLORS["amber"],
        ),
        (
            "5",
            "Biology support",
            "How does the readout align?",
            f"Pathway concordance {data['pathway_concordance']:.3f}; SHAP enrichment {data['gene_enrichment']:.1f}x.",
            COLORS["green"],
        ),
    ]
    y = 900
    for row in rows:
        draw_stack_row(draw, (205, y, 2050, y + 118), *row)
        y += 132


def small_metric(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], value: str, label: str, color: str) -> None:
    x1, y1, x2, _ = box
    rounded(draw, box, 18, "#142033", color, 2)
    text(draw, (x1 + 24, y1 + 20), value, F["stat"], COLORS["text"])
    paragraph(draw, (x1 + 24, y1 + 92), label, F["tiny"], x2 - x1 - 48, COLORS["muted"], 3)


def layer_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    label: str,
    title: str,
    body: str,
    color: str,
    chips: list[str],
) -> None:
    x1, y1, x2, _ = box
    rounded(draw, box, 24, blend("#101823", color, 0.10), color, 2)
    text(draw, (x1 + 30, y1 + 26), label, F["tiny_bold"], color)
    text(draw, (x1 + 30, y1 + 64), title, F["h3"], COLORS["text"])
    paragraph(draw, (x1 + 30, y1 + 108), body, F["tiny"], x2 - x1 - 60, COLORS["muted"], 4)
    chip_x = x1 + 30
    chip_y = y1 + 184
    for chip in chips:
        chip_w = int(draw.textlength(chip, font=F["micro_bold"]) + 34)
        rounded(draw, (chip_x, chip_y, chip_x + chip_w, chip_y + 34), 12, "#111B28", color, 1)
        text(draw, (chip_x + 17, chip_y + 9), chip, F["micro_bold"], COLORS["muted"])
        chip_x += chip_w + 10


def draw_layer_map(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    panel = (2190, 680, 3690, 1620)
    rounded(draw, panel, 30, COLORS["panel"], "#2A394D", 2)
    text(draw, (2240, 730), "B. Where v1-v9 belong", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (2240, 785),
        "Version labels are project history; the slide deck reads them as evidence layers.",
        F["small"],
        1240,
        COLORS["muted"],
        7,
    )
    layer_card(
        draw,
        (2240, 900, 3630, 1150),
        "v1-v7",
        "Benchmark and robustness evidence",
        "Mission-held-out scores, controls, method hardening, model probes, and biological support form the core benchmark spine.",
        COLORS["teal"],
        ["transfer", "controls", "DGE", "models"],
    )
    layer_card(
        draw,
        (2240, 1185, 2928, 1458),
        "v8",
        "Translation hypotheses",
        "BRIDGE, DECOMPOSE, INTERVENE, and CAUSAL organize follow-up ideas after benchmark evidence is visible.",
        COLORS["amber"],
        ["BRIDGE", "DECOMPOSE", "INTERVENE", "CAUSAL"],
    )
    layer_card(
        draw,
        (2945, 1185, 3630, 1458),
        "v9",
        "Platform readiness",
        "Manifests, evaluator, run records, and data packaging make the benchmark auditable and reusable.",
        COLORS["blue"],
        ["manifest", "evaluator", "run record", "package"],
    )
    rounded(draw, (2240, 1500, 3630, 1568), 16, blend("#101823", COLORS["green"], 0.10), COLORS["green"], 2)
    text(draw, (2268, 1519), "Takeaway", F["tiny_bold"], COLORS["green"])
    text(draw, (2438, 1519), "Use the strongest layer that actually exists for each result.", F["tiny"], COLORS["muted"])


def draw_metric_strip(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    y = 1662
    metrics = [
        (f"{data['tissues']}", "tissue tasks", COLORS["teal"]),
        (f"{data['lomo_folds']}", "hidden-mission folds", COLORS["blue"]),
        (f"{data['model_configs']}", "tested model configs", COLORS["violet"]),
        (f"{data['control_rows']}", "control rows", COLORS["amber"]),
        (f"{data['pathway_concordance']:.3f}", "pathway concordance", COLORS["green"]),
    ]
    x = 175
    for value, label, color in metrics:
        small_metric(draw, (x, y, x + 680, y + 132), value, label, color)
        x += 705


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (150, 1886, 3690, 2068), 22, "#111B28", "#26384E", 2)
    text(draw, (190, 1928), "Takeaway", F["tiny_bold"], COLORS["green"])
    paragraph(
        draw,
        (360, 1928),
        "The deck moves from fixed benchmark tasks to conserved biology, follow-up analyses, and a portable benchmark package.",
        F["small"],
        2760,
        COLORS["muted"],
        7,
    )
    text(draw, (3510, 1998), "3", F["h2"], COLORS["teal"])


def render() -> Image.Image:
    data = load_metrics()
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image)
    draw_background(draw)
    draw_header(draw, data)
    draw_title(draw)
    draw_claim_stack(draw, data)
    draw_layer_map(draw, data)
    draw_metric_strip(draw, data)
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
    ImageOps.grayscale(image).convert("RGB").save(GRAY_PATH, quality=95)
    data = load_metrics()
    manifest = {
        "title": "Evidence Layers Organize The Story",
        "slide": 3,
        "status": "ready",
        "role": "v1-v9 evidence-layer story map",
        "created": "2026-06-11",
        "outputs": {
            "png": str(OUT_PATH.relative_to(ROOT)),
            "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        },
        "metrics": data,
        "copy": {
            "headline": "Evidence Layers Organize The Story",
            "subtitle": "The same project carries benchmark results, biological support, translation hypotheses, and platform readiness as separate evidence layers.",
            "audience_read": "Use the strongest layer that actually exists for each result.",
            "takeaway": "The deck moves from fixed benchmark tasks to conserved biology, follow-up analyses, and a portable benchmark package.",
        },
        "source_manifests": [
            str(EVIDENCE_STACK.relative_to(ROOT)),
            str(METHOD_HARDENING.relative_to(ROOT)),
            str(EXTERNAL_VALIDATION.relative_to(ROOT)),
        ],
        "automatic_metrics": image_metrics(OUT_PATH),
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"png": str(OUT_PATH), "grayscale": str(GRAY_PATH), "manifest": str(MANIFEST_PATH)}, indent=2))


if __name__ == "__main__":
    main()

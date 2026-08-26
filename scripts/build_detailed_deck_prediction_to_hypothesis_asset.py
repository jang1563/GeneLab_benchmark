#!/usr/bin/env python3
"""Build slide 42 asset: prediction becomes hypothesis through layered readouts."""

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
OUT_DIR = ASSET_ROOT / "prediction_hypothesis"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "prediction_becomes_hypothesis_through_evidence_gates_premium.png"
GRAY_PATH = OUT_DIR / "prediction_becomes_hypothesis_through_evidence_gates_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "prediction_hypothesis_manifest.json"

CORE_MANIFEST = ASSET_ROOT / "core_takeaway" / "core_takeaway_manifest.json"
EVIDENCE_MANIFEST = ASSET_ROOT / "evidence_stack" / "evidence_stack_manifest.json"
TRANSLATION_MANIFEST = ASSET_ROOT / "translation_details" / "translation_details_manifest.json"
SYSTEMS_QA = ASSET_ROOT / "systems_biology" / "SYSTEMS_BIOLOGY_ASSET_VISUAL_QA.md"
V6_EVAL = ROOT / "v6" / "evaluation"
V6_TARGET = V6_EVAL / "V6_F_drug_target_validation.json"

COLORS = {
    "bg": "#0C111A",
    "bg2": "#091019",
    "header": "#0F1824",
    "footer": "#080D14",
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
    "red": "#F17C88",
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
    "stat": load_font(72, True),
    "stat2": load_font(62, True),
    "axis": load_font(20),
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


def badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(230, int(draw.textlength(value, font=F["tiny_bold"]) + 76))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_data() -> dict[str, object]:
    core = load_json(CORE_MANIFEST)["metrics"]
    evidence = load_json(EVIDENCE_MANIFEST)["metrics"]
    translation = load_json(TRANSLATION_MANIFEST)["metrics"]
    target = load_json(V6_TARGET)
    return {
        "core": core,
        "evidence": evidence,
        "translation": translation,
        "target": target,
    }


def background() -> Image.Image:
    image = Image.new("RGB", (W, H), COLORS["bg"])
    draw = ImageDraw.Draw(image)
    for y in range(H):
        t = y / H
        draw.line((0, y, W, y), fill=blend(COLORS["bg"], COLORS["bg2"], t))
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill=rgba(COLORS["grid"], 25), width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill=rgba(COLORS["grid"], 22), width=1)
    draw.rectangle((0, 0, W, 260), fill=COLORS["header"])
    draw.rectangle((0, 1900, W, H), fill=COLORS["footer"])
    return image


def draw_status_badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, color: str) -> None:
    rounded(draw, (x, y, x + 250, y + 54), 18, "#102133", color, 2)
    text(draw, (x + 22, y + 15), label, F["tiny_bold"], COLORS["text"])


def draw_stage_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    number: str,
    title: str,
    question: str,
    supported_readout: str,
    color: str,
    proof_rows: list[tuple[str, str]],
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 26, COLORS["panel"], "#26364B", 2)
    draw.ellipse((x1 + 34, y1 + 30, x1 + 96, y1 + 92), fill=color)
    text(draw, (x1 + 65, y1 + 51), number, F["tiny_bold"], COLORS["ink"], "mm")
    text(draw, (x1 + 122, y1 + 32), title, F["h2"], COLORS["text"])
    paragraph(draw, (x1 + 122, y1 + 86), question, F["small"], x2 - x1 - 160, COLORS["muted"], 7)
    rounded(draw, (x1 + 38, y1 + 176, x2 - 38, y1 + 274), 20, COLORS["panel2"], color, 2)
    text(draw, (x1 + 66, y1 + 199), "Supported readout", F["micro_bold"], COLORS["dim"])
    paragraph(draw, (x1 + 66, y1 + 228), supported_readout, F["tiny_bold"], x2 - x1 - 120, COLORS["text"], 5)
    y = y1 + 318
    for label, value in proof_rows:
        rounded(draw, (x1 + 38, y, x2 - 38, y + 62), 16, COLORS["panel3"], "#2A394D", 1)
        text(draw, (x1 + 60, y + 14), label, F["micro_bold"], COLORS["dim"])
        text(draw, (x2 - 60, y + 14), value, F["tiny_bold"], color, "ra")
        y += 76


def draw_claim_ladder(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    core = data["core"]
    evidence = data["evidence"]
    translation = data["translation"]
    target = data["target"]
    x0, y0 = 120, 610
    card_w, card_h, gap = 830, 770, 55
    specs = [
        (
            "01",
            "Prediction Score",
            "Does a model generalize across held-out missions?",
            "Benchmark signal is supported for tested mouse bulk tasks.",
            COLORS["blue"],
            [
                ("mission tasks", f"{core['mission_tasks']}"),
                ("LOMO folds", f"{core['lomo_folds']}"),
                ("top tissue AUROC", f"{core['top_tissue']} {core['top_tissue_auroc']:.2f}"),
                ("held-out checks", f"thymus {core['heldout_thymus_auroc']:.3f} | skin {core['heldout_skin_auroc']:.3f}"),
            ],
        ),
        (
            "02",
            "Robust Pattern",
            "How does the signal behave across controls and method variation?",
            "Model result becomes a stable comparison surface across matched configurations.",
            COLORS["teal"],
            [
                ("model configs", f"{evidence['model_configs']}"),
                ("significant rows", f"{evidence['significant_model_rows']}"),
                ("control rows", f"{evidence['control_rows']} | {evidence['significant_control_rows']} sig"),
                ("DGE rank stability", f"rho {evidence['dge_log2fc_spearman_mean']:.3f}"),
            ],
        ),
        (
            "03",
            "Biology Interpretation",
            "Which biological layer explains or contextualizes the score?",
            "Immune, TF, metabolic, target, and biomarker lanes add interpretation.",
            COLORS["green"],
            [
                ("interpretation lanes", "5"),
                ("pathway support", f"{evidence['pathway_mean_concordance']:.3f} concordance"),
                ("gene support", f"{evidence['gene_enrichment_ratio']:.1f}x enrichment"),
                ("target triage", f"{target['tier_counts']['A']} A | {target['tier_counts']['B']} B tiers"),
            ],
        ),
        (
            "04",
            "Translation + Hypothesis",
            "Which part of the biology carries into human readback?",
            "Human follow-up stays at pathway, readback, and target-queue level.",
            COLORS["amber"],
            [
                ("ortholog map", f"{translation['ortholog_map_size']:,}"),
                ("pathway bridge", f"mean rho {translation['pathway_mean_rho']:.3f}"),
                ("panel readback", f"{translation['biomarker_detected']}/{translation['biomarker_panel']} detected"),
                ("target readback", f"{translation['target_detected']}/{translation['target_genes']} detected"),
            ],
        ),
    ]
    for i, spec in enumerate(specs):
        x = x0 + i * (card_w + gap)
        draw_stage_card(draw, (x, y0, x + card_w, y0 + card_h), *spec)
        if i < len(specs) - 1:
            ax = x + card_w + 12
            ay = y0 + 382
            draw.line((ax, ay, ax + gap - 24, ay), fill=COLORS["dim"], width=4)
            draw.polygon([(ax + gap - 24, ay - 12), (ax + gap - 24, ay + 12), (ax + gap - 4, ay)], fill=COLORS["dim"])


def draw_reader_rule(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1435, 3720, 1690
    rounded(draw, (x1, y1, x2, y2), 30, "#101927", "#26364B", 2)
    text(draw, (x1 + 42, y1 + 36), "Reading path", F["small_bold"], COLORS["muted"])
    rules = [
        ("Score", "supports benchmark result", COLORS["blue"]),
        ("Robustness", "stabilizes the readout", COLORS["teal"]),
        ("Biology", "adds interpretation lanes", COLORS["green"]),
        ("Translation", "narrows human readback", COLORS["amber"]),
        ("Hypothesis", "queues follow-up tests", COLORS["violet"]),
    ]
    node_w, gap = 620, 34
    start_x = x1 + 310
    y = y1 + 84
    for i, (title, desc, color) in enumerate(rules):
        nx = start_x + i * (node_w + gap)
        rounded(draw, (nx, y, nx + node_w, y + 102), 20, COLORS["panel2"], color, 2)
        text(draw, (nx + 26, y + 18), title, F["small_bold"], COLORS["text"])
        text(draw, (nx + 26, y + 58), desc, F["tiny"], COLORS["muted"])
        if i < len(rules) - 1:
            ax = nx + node_w + 7
            ay = y + 51
            draw.line((ax, ay, ax + gap - 16, ay), fill=COLORS["dim"], width=3)
            draw.polygon([(ax + gap - 16, ay - 8), (ax + gap - 16, ay + 8), (ax + gap - 3, ay)], fill=COLORS["dim"])


def draw_transition_panel(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1732, 3720, 1858
    rounded(draw, (x1, y1, x2, y2), 26, COLORS["panel"], "#26364B", 2)
    text(draw, (x1 + 42, y1 + 33), "Act 4 close", F["small_bold"], COLORS["green"])
    text(draw, (x1 + 280, y1 + 33), "The deck now has a disciplined route from prediction to interpretation.", F["small_bold"], COLORS["text"])
    text(draw, (x1 + 280, y1 + 75), "Act 5 can introduce v8 as incubation because the benchmark and biology layers are already visible.", F["small"], COLORS["muted"])
    draw_status_badge(draw, x2 - 915, y1 + 36, "benchmark readout", COLORS["blue"])
    draw_status_badge(draw, x2 - 635, y1 + 36, "biology layer", COLORS["green"])
    draw_status_badge(draw, x2 - 355, y1 + 36, "hypothesis queue", COLORS["amber"])


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, "#0F1824", "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "Act 4 synthesis: score, robustness, biology, translation, and hypothesis stay in separate layers.",
        F["small"],
        2200,
        COLORS["muted"],
        8,
    )
    text(draw, (3580, 1960), "42", F["h2"], COLORS["green"], "ra")
    text(draw, (3580, 2010), "Act 4 synthesis", F["tiny_bold"], COLORS["muted"], "ra")


def build() -> None:
    data = load_data()
    image = background()
    draw = ImageDraw.Draw(image)

    core = data["core"]
    translation = data["translation"]
    text(draw, (120, 72), "SLIDE 42 | ACT 4 | SYNTHESIS", F["kicker"], COLORS["green"])
    bx = 2180
    bx = badge(draw, bx, 56, "TASKS", f"{core['mission_tasks']} tissues", COLORS["blue"])
    bx = badge(draw, bx, 56, "ROBUSTNESS", f"{data['evidence']['model_configs']} configs", COLORS["teal"])
    bx = badge(draw, bx, 56, "BIOLOGY", "5 lanes", COLORS["green"])
    badge(draw, bx, 56, "TRANSLATION", f"{translation['biomarker_detected']}/{translation['biomarker_panel']} panel", COLORS["amber"])

    text(draw, (120, 330), "Prediction Becomes Hypothesis Through Layered Readouts", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 432),
        "The biology story is interpretable when each statement stays in the layer supported by its measurements.",
        F["subtitle"],
        2820,
        COLORS["muted"],
        8,
    )

    draw_claim_ladder(draw, data)
    draw_reader_rule(draw)
    draw_transition_panel(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    gray_rgb = ImageOps.colorize(gray, black="#080D14", white="#F3F7FC")
    gray_rgb.save(GRAY_PATH, quality=96)

    stat = ImageStat.Stat(gray)
    manifest = {
        "title": "Prediction Becomes Hypothesis Through Layered Readouts",
        "claim": "Act 4 closes by separating benchmark score, robustness, biology interpretation, translation readback, and hypothesis generation.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "source_manifests": {
            "core_takeaway": str(CORE_MANIFEST.relative_to(ROOT)),
            "evidence_stack": str(EVIDENCE_MANIFEST.relative_to(ROOT)),
            "translation_details": str(TRANSLATION_MANIFEST.relative_to(ROOT)),
            "systems_biology_qa": str(SYSTEMS_QA.relative_to(ROOT)),
        },
        "source_json": {"v6_target": str(V6_TARGET.relative_to(ROOT))},
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": {
            "mission_tasks": core["mission_tasks"],
            "lomo_folds": core["lomo_folds"],
            "model_configs": data["evidence"]["model_configs"],
            "significant_control_rows": data["evidence"]["significant_control_rows"],
            "biology_lanes": 5,
            "pathway_mean_rho": translation["pathway_mean_rho"],
            "biomarker_detected": translation["biomarker_detected"],
            "biomarker_panel": translation["biomarker_panel"],
            "target_detected": translation["target_detected"],
            "target_genes": translation["target_genes"],
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUT_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    build()

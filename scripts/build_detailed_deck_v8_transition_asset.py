#!/usr/bin/env python3
"""Build slide 43 asset: v8 starts after the core readout."""

from __future__ import annotations

import csv
import json
import re
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
OUT_DIR = ASSET_ROOT / "v8_transition"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "v8_starts_after_the_evidence_gates_premium.png"
GRAY_PATH = OUT_DIR / "v8_starts_after_the_evidence_gates_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "v8_transition_manifest.json"
QA_NOTE = OUT_DIR / "V8_TRANSITION_ASSET_VISUAL_QA.md"

CORE_MANIFEST = ASSET_ROOT / "core_takeaway" / "core_takeaway_manifest.json"
EVIDENCE_MANIFEST = ASSET_ROOT / "evidence_stack" / "evidence_stack_manifest.json"
TRANSLATION_MANIFEST = ASSET_ROOT / "translation_details" / "translation_details_manifest.json"
V8_SUMMARY = ROOT / "v8" / "RESULTS_SUMMARY.md"
V8_BRIDGE = ROOT / "v8" / "bridge" / "evaluation" / "supervised_conservation.json"
V8_SPECIES = ROOT / "v8" / "bridge" / "evaluation" / "species_transfer_nes.json"
V8_OFFLINE = ROOT / "v8" / "intervene" / "evaluation" / "offline_reversal_summary.json"
V8_PARETO = ROOT / "v8" / "intervene" / "evaluation" / "pareto_front.csv"
V8_CRISPR = ROOT / "v8" / "intervene" / "evaluation" / "crispr_orthog_summary.json"
V8_CAUSAL = ROOT / "v8" / "causal" / "evaluation" / "icp_dag_summary.json"

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
    "title": load_font(88, True),
    "subtitle": load_font(38),
    "section": load_font(46, True),
    "h2": load_font(38, True),
    "h3": load_font(31, True),
    "body": load_font(28),
    "body_bold": load_font(28, True),
    "small": load_font(24),
    "small_bold": load_font(24, True),
    "tiny": load_font(21),
    "tiny_bold": load_font(21, True),
    "micro": load_font(18),
    "micro_bold": load_font(18, True),
    "metric": load_font(62, True),
    "metric2": load_font(48, True),
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
    width = max(230, int(draw.textlength(value, font=F["tiny_bold"]) + 76))
    rounded(draw, (x, y, x + width, y + 72), 18, "#172335", color, 2)
    text(draw, (x + 22, y + 13), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 22, y + 39), value, F["tiny_bold"], COLORS["text"])
    return x + width + 18


def arrow(draw: ImageDraw.ImageDraw, x1: int, y1: int, x2: int, y2: int, color: str, width: int = 5) -> None:
    draw.line((x1, y1, x2, y2), fill=color, width=width)
    if x2 >= x1:
        pts = [(x2, y2), (x2 - 24, y2 - 13), (x2 - 24, y2 + 13)]
    else:
        pts = [(x2, y2), (x2 + 24, y2 - 13), (x2 + 24, y2 + 13)]
    draw.polygon(pts, fill=color)


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_interaction_range() -> str:
    summary = V8_SUMMARY.read_text(encoding="utf-8")
    match = re.search(r"Interaction dominance:\s*([0-9]+)[–-]([0-9]+)%", summary)
    if match:
        return f"{match.group(1)}-{match.group(2)}%"
    return "44-61%"


def count_pareto_rows() -> int:
    if not V8_PARETO.exists():
        return 2
    with V8_PARETO.open(newline="", encoding="utf-8") as handle:
        return max(0, sum(1 for _ in csv.DictReader(handle)))


def crispr_support_count(data: dict[str, object]) -> int:
    count = 0
    for tissue in data.get("tissues", {}).values():
        if isinstance(tissue, dict):
            count += len(tissue.get("validated_drugs", {}))
    return count


def load_data() -> dict[str, object]:
    core = load_json(CORE_MANIFEST)["metrics"]
    evidence = load_json(EVIDENCE_MANIFEST)["metrics"]
    translation = load_json(TRANSLATION_MANIFEST)["metrics"]
    bridge = load_json(V8_BRIDGE)
    species = load_json(V8_SPECIES)
    offline = load_json(V8_OFFLINE)
    crispr = load_json(V8_CRISPR)
    causal = load_json(V8_CAUSAL)
    return {
        "core": core,
        "evidence": evidence,
        "translation": translation,
        "bridge": bridge,
        "species": species,
        "offline": offline,
        "crispr": crispr,
        "causal": causal,
        "interaction_range": parse_interaction_range(),
        "pareto_count": count_pareto_rows(),
        "crispr_support_count": crispr_support_count(crispr),
    }


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


def draw_gate_row(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    step: str,
    title: str,
    metric: str,
    body: str,
    color: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 20, COLORS["panel2"], "#2A394D", 2)
    draw.ellipse((x1 + 26, y1 + 25, x1 + 82, y1 + 81), fill=color)
    text(draw, (x1 + 54, y1 + 53), step, F["micro_bold"], COLORS["ink"], "mm")
    text(draw, (x1 + 108, y1 + 22), title, F["h3"], COLORS["text"])
    text(draw, (x2 - 34, y1 + 25), metric, F["small_bold"], color, "ra")
    paragraph(draw, (x1 + 108, y1 + 66), body, F["tiny"], x2 - x1 - 155, COLORS["muted"], 6)


def draw_evidence_gates(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    core = data["core"]
    evidence = data["evidence"]
    translation = data["translation"]
    x1, y1, x2, y2 = 120, 610, 1445, 1438
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 35), "Core readout is established", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 96),
        "The core deck has already separated score, robustness, biology, and translation readback.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )
    rows = [
        (
            "01",
            "Benchmark score",
            f"{core['mission_tasks']} tasks | {core['lomo_folds']} folds",
            f"{evidence['model_configs']} model configurations read through one held-out mission contract.",
            COLORS["blue"],
        ),
        (
            "02",
            "Robustness surface",
            f"{evidence['significant_control_rows']} control hits",
            "Control and method-hardening views keep the score lane auditable.",
            COLORS["teal"],
        ),
        (
            "03",
            "Biology interpretation",
            "5 lanes",
            "Immune, TF, metabolic, target, and biomarker layers explain where the signal lives.",
            COLORS["green"],
        ),
        (
            "04",
            "Translation readback",
            f"{translation['biomarker_detected']}/{translation['biomarker_panel']} panel",
            f"{translation['ortholog_map_size']:,} orthologs and {translation['target_detected']}/{translation['target_genes']} target genes define human-facing readback.",
            COLORS["amber"],
        ),
    ]
    top = y1 + 205
    for i, row in enumerate(rows):
        draw_gate_row(draw, (x1 + 38, top + i * 144, x2 - 38, top + 118 + i * 144), *row)


def draw_threshold(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 1518, 666, 1904, 1388
    rounded(draw, (x1, y1, x2, y2), 34, "#101A28", "#34465E", 2)
    draw.rectangle((x1 + 64, y1 + 104, x1 + 102, y2 - 104), fill=COLORS["blue"])
    draw.rectangle((x2 - 102, y1 + 104, x2 - 64, y2 - 104), fill=COLORS["amber"])
    rounded(draw, (x1 + 78, y1 + 76, x2 - 78, y1 + 136), 18, "#18263A", "#3D516D", 2)
    rounded(draw, (x1 + 78, y2 - 136, x2 - 78, y2 - 76), 18, "#18263A", "#3D516D", 2)
    text(draw, ((x1 + x2) / 2, y1 + 202), "CORE", F["metric2"], COLORS["text"], "mm")
    text(draw, ((x1 + x2) / 2, y1 + 272), "READOUT", F["metric2"], COLORS["text"], "mm")
    text(draw, ((x1 + x2) / 2, y1 + 342), "ESTABLISHED", F["small_bold"], COLORS["teal"], "mm")
    rounded(draw, (x1 + 58, y1 + 410, x2 - 58, y1 + 486), 18, "#18263A", COLORS["blue"], 2)
    text(draw, ((x1 + x2) / 2, y1 + 436), "benchmark readout", F["tiny_bold"], COLORS["text"], "mm")
    text(draw, ((x1 + x2) / 2, y1 + 462), "visible", F["micro_bold"], COLORS["muted"], "mm")
    rounded(draw, (x1 + 58, y1 + 510, x2 - 58, y1 + 586), 18, "#18263A", COLORS["amber"], 2)
    text(draw, ((x1 + x2) / 2, y1 + 536), "translation layers", F["tiny_bold"], COLORS["text"], "mm")
    text(draw, ((x1 + x2) / 2, y1 + 562), "visible", F["micro_bold"], COLORS["muted"], "mm")
    rounded(draw, (x1 + 58, y2 - 108, x2 - 58, y2 - 54), 16, "#18263A", COLORS["green"], 2)
    text(draw, ((x1 + x2) / 2, y2 - 82), "ready for v8 incubation", F["micro_bold"], COLORS["text"], "mm")
    arrow(draw, 1444, 1018, 1518, 1018, COLORS["dim"], 5)
    arrow(draw, 1904, 1018, 2025, 1018, COLORS["dim"], 5)


def draw_module_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    tag: str,
    title: str,
    metric: str,
    input_line: str,
    readout_line: str,
    color: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 26, COLORS["panel"], "#2A394D", 2)
    rounded(draw, (x1 + 30, y1 + 28, x1 + 178, y1 + 78), 16, rgba(color, 40), color, 2)
    text(draw, (x1 + 104, y1 + 53), tag, F["tiny_bold"], COLORS["text"], "mm")
    text(draw, (x1 + 204, y1 + 28), title, F["h2"], COLORS["text"])
    text(draw, (x1 + 32, y1 + 105), metric, F["metric2"], color)
    rounded(draw, (x1 + 34, y1 + 172, x2 - 34, y1 + 222), 15, COLORS["panel2"], "#2A394D", 1)
    text(draw, (x1 + 58, y1 + 187), "Input", F["micro_bold"], COLORS["dim"])
    text(draw, (x1 + 158, y1 + 185), input_line, F["tiny_bold"], COLORS["muted"])
    rounded(draw, (x1 + 34, y1 + 236, x2 - 34, y1 + 282), 15, COLORS["panel3"], "#2A394D", 1)
    text(draw, (x1 + 58, y1 + 250), "Readout", F["micro_bold"], COLORS["dim"])
    text(draw, (x1 + 158, y1 + 248), readout_line, F["tiny"], COLORS["muted"])


def draw_v8_modules(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    bridge = data["bridge"]
    species = data["species"]
    offline = data["offline"]
    causal = data["causal"]
    x1, y1, x2, y2 = 2025, 610, 3720, 1438
    rounded(draw, (x1, y1, x2, y2), 30, "#101A28", "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 35), "v8 incubation modules", F["section"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 44, y1 + 96),
        "Each module starts from a bounded evidence lane and converts it into the next analysis question.",
        F["small"],
        x2 - x1 - 88,
        COLORS["muted"],
        7,
    )
    card_w, card_h = 790, 296
    gap_x, gap_y = 38, 28
    top = y1 + 205
    left = x1 + 38
    modules = [
        (
            "BRIDGE",
            "Human Link",
            f"AUROC {bridge['rf']['cv_mean_auroc_aug']:.3f}",
            f"{bridge['n_pathways']} pathways | {species['n_comparisons']} comparisons",
            "pathway-conservation prior",
            COLORS["blue"],
        ),
        (
            "DECOMPOSE",
            "Stressor Surface",
            data["interaction_range"],
            "HLU, radiation, time, and interaction terms",
            "regime-change questions",
            COLORS["teal"],
        ),
        (
            "INTERVENE",
            "Perturbation Priority",
            f"{data['pareto_count']} Pareto axes",
            f"{offline['drugs_scored']} scored candidates | CRISPR support",
            "follow-up experiment queue",
            COLORS["amber"],
        ),
        (
            "CAUSAL",
            "Assumption Map",
            f"{causal['n_environments']} environments",
            f"T ICP {causal['stressor_aggregate_icp']['T']['mean']:.3f}",
            "DAG-organized analysis map",
            COLORS["violet"],
        ),
    ]
    for i, module in enumerate(modules):
        col, row = i % 2, i // 2
        x = left + col * (card_w + gap_x)
        y = top + row * (card_h + gap_y)
        draw_module_card(draw, (x, y, x + card_w, y + card_h), *module)


def draw_reader_move(draw: ImageDraw.ImageDraw) -> None:
    x1, y1, x2, y2 = 120, 1500, 3720, 1818
    rounded(draw, (x1, y1, x2, y2), 30, COLORS["panel"], "#2B3A50", 2)
    text(draw, (x1 + 42, y1 + 38), "Act 5 move", F["h2"], COLORS["text"])
    paragraph(
        draw,
        (x1 + 42, y1 + 91),
        "Read v8 as translational incubation: the deck moves from what the benchmark supports to what follow-up analysis should ask next.",
        F["small"],
        1600,
        COLORS["muted"],
        7,
    )
    steps = [
        ("Score", "mission-held-out", COLORS["blue"]),
        ("Robustness", "controls + methods", COLORS["teal"]),
        ("Biology", "interpretation lanes", COLORS["green"]),
        ("Translation", "human readback", COLORS["amber"]),
        ("Incubation", "v8 modules", COLORS["violet"]),
    ]
    node_w, gap = 560, 60
    start_x, y = x1 + 450, y1 + 185
    for i, (title, desc, color) in enumerate(steps):
        nx = start_x + i * (node_w + gap)
        rounded(draw, (nx, y, nx + node_w, y + 92), 20, COLORS["panel2"], color, 2)
        text(draw, (nx + 26, y + 18), title, F["small_bold"], COLORS["text"])
        text(draw, (nx + 26, y + 54), desc, F["tiny"], COLORS["muted"])
        if i < len(steps) - 1:
            arrow(draw, nx + node_w + 10, y + 46, nx + node_w + gap - 12, y + 46, COLORS["dim"], 4)
    text(draw, (x2 - 84, y1 + 74), "Why here?", F["metric2"], COLORS["violet"], "ra")
    text(draw, (x2 - 84, y1 + 130), "because the core readout is already visible", F["small_bold"], COLORS["muted"], "ra")


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1910, 3720, 2076), 24, COLORS["header"], "#253247", 2)
    text(draw, (160, 1940), "Takeaway", F["small_bold"], COLORS["teal"])
    paragraph(
        draw,
        (160, 1978),
        "Act 5 transition: v8 converts benchmarked evidence into translational follow-up analysis modules.",
        F["small"],
        2450,
        COLORS["muted"],
        8,
    )
    text(draw, (3580, 1960), "43", F["h2"], COLORS["violet"], "ra")
    text(draw, (3580, 2010), "v8 transition", F["tiny_bold"], COLORS["muted"], "ra")


def build() -> None:
    data = load_data()
    image = background()
    draw = ImageDraw.Draw(image)

    text(draw, (120, 72), "SLIDE 43 | ACT 5 | V8 TRANSITION", F["kicker"], COLORS["violet"])
    bx = 2180
    bx = badge(draw, bx, 56, "READOUT", "ready through 42", COLORS["green"])
    bx = badge(draw, bx, 56, "V8", "4 modules", COLORS["violet"])
    bx = badge(draw, bx, 56, "BRIDGE", f"AUROC {data['bridge']['rf']['cv_mean_auroc_aug']:.3f}", COLORS["blue"])
    badge(draw, bx, 56, "V8 INPUTS", "evaluations", COLORS["amber"])

    text(draw, (120, 330), "v8 Starts After The Core Readout", F["title"], COLORS["text"])
    paragraph(
        draw,
        (124, 434),
        "Benchmarking establishes the core readout; v8 uses that readout to incubate translational questions.",
        F["subtitle"],
        3100,
        COLORS["muted"],
        8,
    )

    draw_evidence_gates(draw, data)
    draw_threshold(draw)
    draw_v8_modules(draw, data)
    draw_reader_move(draw)
    draw_footer(draw)

    image.save(OUT_PATH, quality=96)
    gray = ImageOps.grayscale(image)
    ImageOps.colorize(gray, black="#080D14", white="#F3F7FC").save(GRAY_PATH, quality=96)
    stat = ImageStat.Stat(gray)

    manifest = {
        "title": "v8 Starts After The Core Readout",
        "claim": "Act 5 begins after score, robustness, biology, and translation layers have established the core readout.",
        "output": str(OUT_PATH.relative_to(ROOT)),
        "grayscale": str(GRAY_PATH.relative_to(ROOT)),
        "source_manifests": {
            "core_takeaway": str(CORE_MANIFEST.relative_to(ROOT)),
            "evidence_stack": str(EVIDENCE_MANIFEST.relative_to(ROOT)),
            "translation_details": str(TRANSLATION_MANIFEST.relative_to(ROOT)),
        },
        "source_files": {
            "v8_summary": str(V8_SUMMARY.relative_to(ROOT)),
            "bridge": str(V8_BRIDGE.relative_to(ROOT)),
            "species": str(V8_SPECIES.relative_to(ROOT)),
            "offline_reversal": str(V8_OFFLINE.relative_to(ROOT)),
            "pareto_front": str(V8_PARETO.relative_to(ROOT)),
            "crispr": str(V8_CRISPR.relative_to(ROOT)),
            "causal": str(V8_CAUSAL.relative_to(ROOT)),
        },
        "size": [W, H],
        "gray_mean": stat.mean[0],
        "gray_stddev": stat.stddev[0],
        "metrics": {
            "mission_tasks": data["core"]["mission_tasks"],
            "lomo_folds": data["core"]["lomo_folds"],
            "model_configs": data["evidence"]["model_configs"],
            "significant_control_rows": data["evidence"]["significant_control_rows"],
            "ortholog_map_size": data["translation"]["ortholog_map_size"],
            "biomarker_detected": data["translation"]["biomarker_detected"],
            "biomarker_panel": data["translation"]["biomarker_panel"],
            "target_detected": data["translation"]["target_detected"],
            "target_genes": data["translation"]["target_genes"],
            "v8_modules": 4,
            "bridge_rf_auroc_aug": data["bridge"]["rf"]["cv_mean_auroc_aug"],
            "bridge_rf_delta": data["bridge"]["rf"]["delta_aug_minus_base"]["mean"],
            "bridge_pathways": data["bridge"]["n_pathways"],
            "species_comparisons": data["species"]["n_comparisons"],
            "interaction_range": data["interaction_range"],
            "offline_drugs_scored": data["offline"]["drugs_scored"],
            "pareto_axes": data["pareto_count"],
            "crispr_supported_drugs": data["crispr_support_count"],
            "causal_environments": data["causal"]["n_environments"],
            "causal_t_icp_mean": data["causal"]["stressor_aggregate_icp"]["T"]["mean"],
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    QA_NOTE.write_text(
        "\n".join(
            [
                "# V8 Transition Asset Visual QA",
                "",
                "Slide 43 is designed as an Act 5 transition rather than a result-detail slide.",
                "",
                "Checks to review:",
                "- Header badges remain inside the safe margin.",
                "- Left core readout, center threshold, and right v8 modules have clear separation.",
                "- Arrows do not touch text labels or module cards.",
                "- Reader-move strip reads at presentation scale.",
                "- Grayscale version preserves contrast for blue, teal, green, amber, and violet lanes.",
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

#!/usr/bin/env python3
"""Build the detailed-deck evidence-stack summary asset."""

from __future__ import annotations

import json
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent.parent
ASSET_ROOT = (
    ROOT
    / "outputs"
    / "019e8bd1-3b42-7821-9e1c-beebfa8e2ece"
    / "presentations"
    / "spacebiobench-detailed-deck"
    / "assets"
)
OUT_DIR = ASSET_ROOT / "evidence_stack"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "evidence_stack_turns_scores_into_claim_status_premium.png"
GRAY_PATH = OUT_DIR / "evidence_stack_turns_scores_into_claim_status_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "evidence_stack_manifest.json"

MANIFESTS = {
    "mission": ASSET_ROOT / "mission_holdout" / "mission_heldout_protocol_premium_manifest.json",
    "guard": ASSET_ROOT / "train_only_guard" / "train_only_feature_guard_premium_manifest.json",
    "metric": ASSET_ROOT / "metric_primer" / "metric_primer_auroc_uncertainty_premium_manifest.json",
    "method": ASSET_ROOT / "method_hardening" / "method_hardening_manifest.json",
    "new_models": ASSET_ROOT / "new_model_ideas" / "new_model_ideas_manifest.json",
    "controls": ASSET_ROOT / "negative_controls" / "negative_controls_manifest.json",
    "dge": ASSET_ROOT / "dge_robustness" / "dge_pipeline_robustness_premium_manifest.json",
    "biology": ASSET_ROOT / "external_validation" / "external_biology_validation_premium_manifest.json",
}

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "panel2": "#151F2D",
    "panel3": "#121B28",
    "grid": "#2A3546",
    "axis": "#98A7BA",
    "text": "#F3F7FC",
    "muted": "#A8B4C4",
    "dim": "#667487",
    "blue": "#66A6E8",
    "sky": "#73A7FF",
    "teal": "#5FD3C4",
    "green": "#8BD17C",
    "amber": "#F4C26B",
    "orange": "#E69F00",
    "rose": "#E17882",
    "violet": "#B39DFF",
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
    "title": load_font(76, True),
    "subtitle": load_font(35),
    "h2": load_font(42, True),
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
    "huge": load_font(90, True),
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
    fill: str | tuple[int, int, int] = COLORS["text"],
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
        else:
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


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


def load_data() -> dict[str, object]:
    manifests = {key: load_json(path) for key, path in MANIFESTS.items()}
    mission = manifests["mission"]["summary"]
    guard = manifests["guard"]["summary"]
    metric = manifests["metric"]["decision_rule"]
    method = manifests["method"]["metrics"]
    new_models = manifests["new_models"]["metrics"]
    controls = manifests["controls"]["control_summary"]
    dge = manifests["dge"]["summary"]
    biology = manifests["biology"]["summary"]
    return {
        "mission_tasks": mission["core_tissue_tasks"],
        "lomo_folds": mission["public_lomo_folds"],
        "mission_labels": mission["unique_mission_labels"],
        "train_cols": guard["train_feature_cols"],
        "test_cols": guard["test_feature_cols"],
        "mean_gate": metric["mean_auroc"],
        "ci_gate": metric["ci_lower"],
        "p_gate": metric["perm_p"],
        "configs": method["n_configs"],
        "sig_rows": method["n_significant"],
        "sig_tissues": method["n_tissue_significant"],
        "pca_mean": new_models["pca_lr_mean_auroc"],
        "sc_mean": new_models["scprint2_mean_auroc"],
        "wgcna_mean": new_models["wgcna_gnn_mean_auroc"],
        "control_rows": controls["label_shuffled"]["n"] + controls["random_features"]["n"],
        "control_sig": controls["label_shuffled"]["sig_count_p_lt_0_05"] + controls["random_features"]["sig_count_p_lt_0_05"],
        "dge_rho": dge["log2fc_spearman_mean"],
        "dge_deg_jaccard": dge["deg_jaccard_mean"],
        "pathway_concordance": biology["pathway_mean_concordance"],
        "pathway_tissues": biology["pathway_n_tissues"],
        "gene_enrichment": biology["gene_enrichment_ratio"],
    }


def draw_background(canvas: Image.Image) -> None:
    draw = ImageDraw.Draw(canvas)
    for x in range(0, W, 160):
        draw.line((x, 0, x, H), fill="#142033", width=1)
    for y in range(0, H, 160):
        draw.line((0, y, W, y), fill="#111B2A", width=1)
    draw.rectangle((0, 0, W, 310), fill="#101824")
    draw.rectangle((0, 1855, W, H), fill="#090E15")


def draw_badge(draw: ImageDraw.ImageDraw, x: int, y: int, label: str, value: str, color: str) -> int:
    width = max(180, int(draw.textlength(value, font=F["tiny_bold"]) + 56))
    rounded(draw, (x, y, x + width, y + 62), 18, "#1B2838", color, 2)
    text(draw, (x + 18, y + 12), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 18, y + 34), value, F["tiny_bold"], COLORS["text"])
    return x + width + 16


def draw_panel_header(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], label: str, title: str, subtitle: str) -> None:
    x1, y1, x2, _ = box
    text(draw, (x1 + 44, y1 + 34), label.upper(), F["tiny_bold"], COLORS["teal"])
    text(draw, (x1 + 44, y1 + 66), title, F["h2"], COLORS["text"])
    paragraph(draw, (x1 + 44, y1 + 118), subtitle, F["small"], x2 - x1 - 88, COLORS["muted"], 6)


def draw_reader_panel(draw: ImageDraw.ImageDraw) -> None:
    box = (120, 330, 1010, 1840)
    rounded(draw, box, 30, COLORS["panel"], "#2A394D", 2)
    draw_panel_header(
        draw,
        box,
        "A. Reading path",
        "Read the layers before interpreting",
        "Each layer answers a different audience question before interpretation moves forward.",
    )
    steps = [
        ("1", "Define the hidden unit", "The test set is an entire hidden mission."),
        ("2", "Freeze training choices", "Feature selection and fitting happen before the hidden mission is opened."),
        ("3", "Apply one score gate", "AUROC, CI, and permutation p-value travel together."),
        ("4", "Add companion evidence", "Controls, DGE stability, and biology support set the readout context."),
    ]
    y = 575
    for number, title, body in steps:
        rounded(draw, (170, y, 960, y + 210), 24, "#111B28", "#2A3546", 1)
        rounded(draw, (200, y + 28, 258, y + 86), 18, blend("#111B28", COLORS["teal"], 0.35), COLORS["teal"], 2)
        text(draw, (229, y + 57), number, F["h3"], COLORS["text"], "mm")
        text(draw, (280, y + 32), title, F["h3"], COLORS["text"])
        paragraph(draw, (280, y + 82), body, F["small"], 615, COLORS["muted"], 7)
        y += 245
    rounded(draw, (170, 1615, 960, 1785), 24, blend("#101823", COLORS["amber"], 0.15), COLORS["amber"], 2)
    text(draw, (210, 1650), "Working logic", F["h3"], COLORS["text"])
    paragraph(
        draw,
        (210, 1698),
        "The deck separates score evidence, robustness evidence, biological support, and follow-up hypotheses.",
        F["small_bold"],
        690,
        COLORS["text"],
        8,
    )


def draw_layer(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    idx: int,
    label: str,
    title: str,
    question: str,
    metric: str,
    color: str,
) -> None:
    x1, y1, x2, y2 = box
    fill = blend("#101823", color, 0.14 + idx * 0.01)
    rounded(draw, box, 26, fill, color, 2)
    rounded(draw, (x1 + 28, y1 + 26, x1 + 92, y1 + 90), 18, blend("#101823", color, 0.35), color, 2)
    text(draw, (x1 + 60, y1 + 58), str(idx), F["h3"], COLORS["text"], "mm")
    text(draw, (x1 + 120, y1 + 24), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x1 + 120, y1 + 52), title, F["h3"], COLORS["text"])
    paragraph(draw, (x1 + 120, y1 + 96), question, F["tiny"], 720, COLORS["muted"], 6)
    rounded(draw, (x2 - 385, y1 + 30, x2 - 34, y2 - 30), 20, "#0D1520", "#2A3546", 1)
    paragraph(draw, (x2 - 352, y1 + 54), metric, F["small_bold"], 300, COLORS["text"], 7)


def draw_stack_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    box = (1000, 330, 2665, 1840)
    rounded(draw, box, 30, COLORS["panel"], "#2A394D", 2)
    draw_panel_header(
        draw,
        box,
        "B. Evidence stack",
        "Layers accumulate; they do different work",
        "The stack keeps benchmark evidence separate from robustness checks and biological support.",
    )
    layers = [
        (
            "Task contract",
            "Hidden mission split",
            "Can we read the score as transfer across missions?",
            f"{data['mission_tasks']} tissue tasks\n{data['lomo_folds']} LOMO folds\n{data['mission_labels']} mission labels",
            COLORS["blue"],
        ),
        (
            "Leakage guard",
            "Train-only operations",
            "Do modeling choices stop before the hidden mission?",
            f"{data['train_cols']:,} train columns\n{data['test_cols']:,} test columns\nsame feature surface",
            COLORS["teal"],
        ),
        (
            "Metric gate",
            "AUROC + uncertainty",
            "How does a score become a reproducible readout?",
            f"AUROC {data['mean_gate']}\nCI-low {data['ci_gate']}\nperm p {data['p_gate']}",
            COLORS["amber"],
        ),
        (
            "Model breadth",
            "Baseline floor + extensions",
            "Does the lesson hold across model families?",
            f"{data['configs']} configs\n{data['sig_rows']} significant rows\n{data['sig_tissues']} tissues with signal",
            COLORS["violet"],
        ),
        (
            "Robustness",
            "Controls + DGE stability",
            "Do controls and preprocessing checks support the readout?",
            f"{data['control_rows']} control rows\n{data['control_sig']} significant controls\nDGE rho {data['dge_rho']:.3f}",
            COLORS["green"],
        ),
        (
            "Biology support",
            "Pathway-level concordance",
            "Does external biology align with the pathway signal?",
            f"{data['pathway_concordance']:.3f} concordance\n{data['pathway_tissues']} tissues\n{data['gene_enrichment']:.1f}x gene enrichment",
            COLORS["orange"],
        ),
    ]
    y = 550
    layer_h = 182
    gap = 20
    for i, (label, title, question, metric, color) in enumerate(layers, start=1):
        offset = (i - 1) * 18
        draw_layer(draw, (1055 + offset, y, 2610 - offset, y + layer_h), i, label, title, question, metric, color)
        if i < len(layers):
            cx = 1832
            draw.line((cx, y + layer_h + 2, cx, y + layer_h + gap - 3), fill="#3B4A60", width=4)
            draw.polygon([(cx - 10, y + layer_h + gap - 7), (cx + 10, y + layer_h + gap - 7), (cx, y + layer_h + gap + 8)], fill="#3B4A60")
        y += layer_h + gap


def status_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    title: str,
    body: str,
    tag: str,
    color: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 24, "#111B28", color, 2)
    rounded(draw, (x1 + 26, y1 + 26, x1 + 176, y1 + 64), 15, blend("#101823", color, 0.22), color, 1)
    text(draw, (x1 + 44, y1 + 36), tag.upper(), F["micro_bold"], COLORS["text"])
    text(draw, (x1 + 26, y1 + 86), title, F["h3"], COLORS["text"])
    paragraph(draw, (x1 + 26, y1 + 136), body, F["small"], x2 - x1 - 52, COLORS["muted"], 7)


def draw_claim_panel(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    box = (2655, 330, 3730, 1840)
    rounded(draw, box, 30, COLORS["panel"], "#2A394D", 2)
    draw_panel_header(
        draw,
        box,
        "C. What the stack carries",
        "What the first three acts establish",
        "The stack turns multiple evidence layers into a concise benchmark readout.",
    )
    status_card(
        draw,
        (2705, 545, 3680, 825),
        "Benchmark readout",
        "Mission-level transfer scores are comparable because split, feature selection, scoring, and uncertainty are fixed.",
        "supported",
        COLORS["teal"],
    )
    status_card(
        draw,
        (2705, 865, 3680, 1145),
        "Model-family lesson",
        f"Classical baselines remain the benchmark floor across the {data['configs']}-configuration surface and v7 extension probes.",
        "supported",
        COLORS["violet"],
    )
    status_card(
        draw,
        (2705, 1185, 3680, 1465),
        "Biology interpretation",
        "Pathway concordance and DGE rank stability support a biology layer for follow-up prioritization.",
        "support layer",
        COLORS["orange"],
    )
    rounded(draw, (2705, 1510, 3680, 1785), 28, blend("#101823", COLORS["green"], 0.15), COLORS["green"], 2)
    text(draw, (2745, 1550), "Next slide", F["h3"], COLORS["text"])
    paragraph(
        draw,
        (2745, 1602),
        "Convert the stack into the concise core benchmark takeaway before moving into biology and translation layers.",
        F["body_bold"],
        840,
        COLORS["text"],
        9,
    )


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1888, 3730, 2076), 28, "#0B1018", "#263448", 2)
    text(draw, (164, 1926), "Takeaway", F["h3"], COLORS["text"])
    paragraph(
        draw,
        (164, 1979),
        "A score is easiest to trust when the split, training path, metric gate, controls, model breadth, and biology readout are visible in one stack.",
        F["small"],
        2920,
        COLORS["muted"],
        7,
    )
    text(draw, (3570, 1960), "7-31", F["huge"], blend(COLORS["bg"], COLORS["teal"], 0.85), "mm")
    text(draw, (3570, 2030), "layered readout", F["tiny_bold"], COLORS["dim"], "mm")


def render() -> Image.Image:
    data = load_data()
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)

    text(draw, (120, 58), "SPACEBIO-BENCH / DETAILED DECK", F["kicker"], COLORS["teal"])
    text(draw, (120, 108), "Evidence Stack Makes Scores Readable", F["title"], COLORS["text"])
    paragraph(
        draw,
        (122, 205),
        "Split, leakage guard, score gate, model breadth, robustness, and biology support each answer a separate reader question.",
        F["subtitle"],
        2450,
        COLORS["muted"],
        6,
    )
    x = 2525
    x = draw_badge(draw, x, 76, "tasks", f"{data['mission_tasks']}", COLORS["blue"])
    x = draw_badge(draw, x, 76, "folds", f"{data['lomo_folds']}", COLORS["teal"])
    x = draw_badge(draw, x, 76, "configs", f"{data['configs']}", COLORS["violet"])
    x = draw_badge(draw, x, 76, "DGE rho", f"{data['dge_rho']:.3f}", COLORS["green"])
    draw_badge(draw, x, 76, "bio", f"{data['pathway_concordance']:.3f}", COLORS["orange"])

    draw_reader_panel(draw)
    draw_stack_panel(draw, data)
    draw_claim_panel(draw, data)
    draw_footer(draw)

    manifest = {
        "asset": str(OUT_PATH.relative_to(ROOT)),
        "title": "Evidence Stack Makes Scores Readable",
        "source_manifests": [str(path.relative_to(ROOT)) for path in MANIFESTS.values()],
        "metrics": {
            "mission_tasks": int(data["mission_tasks"]),
            "lomo_folds": int(data["lomo_folds"]),
            "mission_labels": int(data["mission_labels"]),
            "train_feature_columns_example": int(data["train_cols"]),
            "model_configs": int(data["configs"]),
            "significant_model_rows": int(data["sig_rows"]),
            "significant_tissues": int(data["sig_tissues"]),
            "control_rows": int(data["control_rows"]),
            "significant_control_rows": int(data["control_sig"]),
            "dge_log2fc_spearman_mean": float(data["dge_rho"]),
            "pathway_mean_concordance": float(data["pathway_concordance"]),
            "gene_enrichment_ratio": float(data["gene_enrichment"]),
        },
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    return canvas


def main() -> None:
    canvas = render()
    canvas.save(OUT_PATH, quality=95)
    canvas.convert("L").convert("RGB").save(GRAY_PATH, quality=95)
    print(json.dumps({"asset": str(OUT_PATH.relative_to(ROOT)), "grayscale": str(GRAY_PATH.relative_to(ROOT)), "manifest": str(MANIFEST_PATH.relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()

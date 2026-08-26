#!/usr/bin/env python3
"""Build the detailed-deck core benchmark takeaway asset."""

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
OUT_DIR = ASSET_ROOT / "core_takeaway"
OUT_DIR.mkdir(parents=True, exist_ok=True)

W, H = 3840, 2160
OUT_PATH = OUT_DIR / "core_benchmark_takeaway_premium.png"
GRAY_PATH = OUT_DIR / "core_benchmark_takeaway_premium_grayscale.png"
MANIFEST_PATH = OUT_DIR / "core_takeaway_manifest.json"

MANIFESTS = {
    "tissue": ASSET_ROOT / "tissue_hierarchy" / "tissue_transfer_hierarchy_premium_manifest.json",
    "heldout": ASSET_ROOT / "heldout_validation" / "heldout_missions_confirm_signal_premium_manifest.json",
    "evidence_stack": ASSET_ROOT / "evidence_stack" / "evidence_stack_manifest.json",
    "method": ASSET_ROOT / "method_hardening" / "method_hardening_manifest.json",
    "new_models": ASSET_ROOT / "new_model_ideas" / "new_model_ideas_manifest.json",
    "external": ASSET_ROOT / "external_validation" / "external_biology_validation_premium_manifest.json",
}

COLORS = {
    "bg": "#0C111A",
    "panel": "#101823",
    "panel2": "#151F2D",
    "grid": "#2A3546",
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
    "mega": load_font(104, True),
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
    tissue_rows = manifests["tissue"]["rows"]
    top_tissue = tissue_rows[0]
    heldout = manifests["heldout"]["data"]
    stack = manifests["evidence_stack"]["metrics"]
    method = manifests["method"]["metrics"]
    external = manifests["external"]["summary"]
    new_models = manifests["new_models"]["metrics"]
    return {
        "top_tissue": top_tissue["tissue"],
        "top_tissue_auroc": top_tissue["auroc"],
        "thymus_minus_liver": manifests["tissue"]["comparisons"]["thymus_minus_liver"],
        "heldout_thymus": heldout["thymus"]["headline_lr"]["auroc"],
        "heldout_skin": heldout["skin"]["headline_lr"]["auroc"],
        "mission_tasks": stack["mission_tasks"],
        "lomo_folds": stack["lomo_folds"],
        "mission_labels": stack["mission_labels"],
        "model_configs": method["n_configs"],
        "significant_model_rows": method["n_significant"],
        "significant_tissues": method["n_tissue_significant"],
        "control_rows": stack["control_rows"],
        "significant_control_rows": stack["significant_control_rows"],
        "dge_rho": stack["dge_log2fc_spearman_mean"],
        "pathway_concordance": external["pathway_mean_concordance"],
        "pathway_tissues": external["pathway_n_tissues"],
        "gene_enrichment": external["gene_enrichment_ratio"],
        "pca_mean": new_models["pca_lr_mean_auroc"],
        "wgcna_mean": new_models["wgcna_gnn_mean_auroc"],
        "sc_mean": new_models["scprint2_mean_auroc"],
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
    width = max(178, int(draw.textlength(value, font=F["tiny_bold"]) + 56))
    rounded(draw, (x, y, x + width, y + 62), 18, "#1B2838", color, 2)
    text(draw, (x + 18, y + 12), label.upper(), F["micro_bold"], COLORS["dim"])
    text(draw, (x + 18, y + 34), value, F["tiny_bold"], COLORS["text"])
    return x + width + 16


def panel_header(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], label: str, title: str, subtitle: str) -> None:
    x1, y1, x2, _ = box
    text(draw, (x1 + 44, y1 + 34), label.upper(), F["tiny_bold"], COLORS["teal"])
    text(draw, (x1 + 44, y1 + 66), title, F["h2"], COLORS["text"])
    paragraph(draw, (x1 + 44, y1 + 118), subtitle, F["small"], x2 - x1 - 88, COLORS["muted"], 6)


def metric_tile(draw: ImageDraw.ImageDraw, box: tuple[int, int, int, int], value: str, label: str, color: str) -> None:
    x1, y1, x2, _ = box
    rounded(draw, box, 24, "#111B28", color, 2)
    text(draw, (x1 + 28, y1 + 28), value, F["stat"], COLORS["text"])
    paragraph(draw, (x1 + 28, y1 + 104), label, F["small"], x2 - x1 - 56, COLORS["muted"], 7)


def draw_core_claim(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    box = (120, 330, 1235, 1840)
    rounded(draw, box, 30, COLORS["panel"], "#2A394D", 2)
    panel_header(
        draw,
        box,
        "A. Core benchmark takeaway",
        "Transfer is measurable and tissue-structured",
        "The first three acts establish the benchmark readout before interpretation begins.",
    )
    rounded(draw, (175, 545, 1180, 890), 30, blend("#101823", COLORS["teal"], 0.16), COLORS["teal"], 2)
    text(draw, (220, 585), "Core readout", F["h3"], COLORS["text"])
    paragraph(
        draw,
        (220, 640),
        "Spaceflight bulk RNA-seq transfer can be evaluated under mission-held-out scoring, with the strongest signals tracking tissue context and task fit.",
        F["body_bold"],
        900,
        COLORS["text"],
        10,
    )
    paragraph(
        draw,
        (220, 780),
        "This gives the project a benchmark floor for comparing models and a biological map for follow-up.",
        F["small"],
        900,
        COLORS["muted"],
        7,
    )
    metric_tile(draw, (175, 960, 505, 1175), str(data["mission_tasks"]), "core tissue tasks", COLORS["blue"])
    metric_tile(draw, (530, 960, 835, 1175), str(data["lomo_folds"]), "mission-held-out folds", COLORS["teal"])
    metric_tile(draw, (860, 960, 1180, 1175), str(data["model_configs"]), "method configs", COLORS["violet"])
    rounded(draw, (175, 1245, 1180, 1645), 28, "#111B28", "#2A3546", 2)
    text(draw, (220, 1285), "Result spine", F["h3"], COLORS["text"])
    rows = [
        (COLORS["teal"], f"{data['top_tissue']} leads", f"mean AUROC {data['top_tissue_auroc']:.3f}; thymus-liver delta {data['thymus_minus_liver']:.3f}"),
        (COLORS["green"], "Held-out missions confirm signal", f"RR-23 thymus {data['heldout_thymus']:.3f}; RR-7 skin {data['heldout_skin']:.3f}"),
        (COLORS["amber"], "Biology layer aligns", f"pathway concordance {data['pathway_concordance']:.3f}; gene enrichment {data['gene_enrichment']:.1f}x"),
    ]
    y = 1350
    for color, title, body in rows:
        draw.ellipse((220, y + 8, 244, y + 32), fill=color)
        text(draw, (265, y), title, F["small_bold"], COLORS["text"])
        text(draw, (265, y + 34), body, F["tiny"], COLORS["muted"])
        y += 78
    rounded(draw, (175, 1705, 1180, 1785), 22, blend("#101823", COLORS["amber"], 0.14), COLORS["amber"], 2)
    text(draw, (220, 1730), "The benchmark floor anchors model comparison and biology follow-up.", F["small_bold"], COLORS["text"])


def claim_card(
    draw: ImageDraw.ImageDraw,
    box: tuple[int, int, int, int],
    number: str,
    title: str,
    statement: str,
    evidence: str,
    status: str,
    color: str,
) -> None:
    x1, y1, x2, y2 = box
    rounded(draw, box, 28, blend("#101823", color, 0.12), color, 2)
    rounded(draw, (x1 + 30, y1 + 30, x1 + 92, y1 + 92), 18, blend("#101823", color, 0.32), color, 2)
    text(draw, (x1 + 61, y1 + 61), number, F["h3"], COLORS["text"], "mm")
    rounded(draw, (x2 - 245, y1 + 34, x2 - 34, y1 + 72), 16, "#0E1722", color, 1)
    text(draw, (x2 - 223, y1 + 44), status.upper(), F["micro_bold"], COLORS["text"])
    text(draw, (x1 + 120, y1 + 30), title, F["h3"], COLORS["text"])
    paragraph(draw, (x1 + 120, y1 + 84), statement, F["small"], x2 - x1 - 170, COLORS["muted"], 7)
    rounded(draw, (x1 + 120, y2 - 84, x2 - 34, y2 - 24), 18, "#0D1520", "#2A3546", 1)
    paragraph(draw, (x1 + 145, y2 - 68), evidence, F["tiny_bold"], x2 - x1 - 205, COLORS["text"], 5)


def draw_claim_status(draw: ImageDraw.ImageDraw, data: dict[str, object]) -> None:
    box = (1225, 330, 2615, 1840)
    rounded(draw, box, 30, COLORS["panel"], "#2A394D", 2)
    panel_header(
        draw,
        box,
        "B. What the first three acts establish",
        "Three supported statements",
        "Each statement is tied to a visible evidence layer, so the audience can inspect the basis for the readout.",
    )
    claim_card(
        draw,
        (1285, 555, 2555, 890),
        "1",
        "Benchmark readout is disciplined",
        "The evaluation unit, leakage guard, metric gate, and controls make model scores comparable across the benchmark surface.",
        f"{data['mission_tasks']} tasks | {data['lomo_folds']} folds | {data['control_rows']} controls / {data['significant_control_rows']} sig",
        "supported",
        COLORS["teal"],
    )
    claim_card(
        draw,
        (1285, 945, 2555, 1280),
        "2",
        "Model lesson is stable",
        "Tuned classical baselines remain the floor while newer model probes add breadth to the comparison.",
        f"{data['model_configs']} configs | {data['significant_model_rows']} significant rows | PCA mean {data['pca_mean']:.3f}",
        "supported",
        COLORS["violet"],
    )
    claim_card(
        draw,
        (1285, 1335, 2555, 1670),
        "3",
        "Biology support is actionable",
        "DGE rank stability and published pathway concordance provide a follow-up map for mechanism-oriented slides.",
        f"DGE rho {data['dge_rho']:.3f} | pathway concordance {data['pathway_concordance']:.3f}",
        "support layer",
        COLORS["orange"],
    )
    rounded(draw, (1285, 1715, 2555, 1788), 22, blend("#101823", COLORS["green"], 0.15), COLORS["green"], 2)
    text(draw, (1325, 1736), "Core status: benchmark evidence is ready for the biology turn.", F["small_bold"], COLORS["text"])


def draw_bridge(draw: ImageDraw.ImageDraw) -> None:
    box = (2605, 330, 3730, 1840)
    rounded(draw, box, 30, COLORS["panel"], "#2A394D", 2)
    panel_header(
        draw,
        box,
        "C. Bridge to biology",
        "Now ask what the signal means",
        "The next act shifts from scoring models to interpreting biological context.",
    )
    lanes = [
        ("Temporal labels", "preservation and recovery context", COLORS["blue"]),
        ("Single-cell pilots", "cell-type context for bulk signals", COLORS["teal"]),
        ("Spatial / weak-signal cases", "where modality limits are informative", COLORS["amber"]),
        ("Systems biology", "immune, TF, metabolic, target, and biomarker layers", COLORS["green"]),
        ("Translation", "mouse-to-human signal transfer with constraints", COLORS["violet"]),
    ]
    y = 565
    for idx, (title, body, color) in enumerate(lanes, start=1):
        rounded(draw, (2665, y, 3670, y + 165), 24, "#111B28", color, 2)
        rounded(draw, (2700, y + 34, 2750, y + 84), 15, blend("#101823", color, 0.28), color, 2)
        text(draw, (2725, y + 59), str(idx), F["small_bold"], COLORS["text"], "mm")
        text(draw, (2780, y + 35), title, F["h3"], COLORS["text"])
        paragraph(draw, (2780, y + 86), body, F["small"], 800, COLORS["muted"], 7)
        if idx < len(lanes):
            draw.line((3168, y + 168, 3168, y + 198), fill="#3B4A60", width=4)
            draw.polygon([(3158, y + 194), (3178, y + 194), (3168, y + 208)], fill="#3B4A60")
        y += 205
    rounded(draw, (2665, 1625, 3670, 1785), 28, blend("#101823", COLORS["teal"], 0.16), COLORS["teal"], 2)
    text(draw, (2710, 1662), "Biology handoff", F["h3"], COLORS["text"])
    paragraph(
        draw,
        (2710, 1710),
        "We move from whether transfer is measurable to which biological layers explain and prioritize it.",
        F["body_bold"],
        870,
        COLORS["text"],
        9,
    )


def draw_footer(draw: ImageDraw.ImageDraw) -> None:
    rounded(draw, (120, 1888, 3730, 2076), 28, "#0B1018", "#263448", 2)
    text(draw, (164, 1926), "Takeaway", F["h3"], COLORS["text"])
    paragraph(
        draw,
        (164, 1979),
        "The core benchmark section closes with a measurable transfer readout, a stable model lesson, and a biology-ready interpretation path.",
        F["small"],
        2850,
        COLORS["muted"],
        7,
    )
    text(draw, (3568, 1960), "33", F["mega"], blend(COLORS["bg"], COLORS["teal"], 0.85), "mm")
    text(draw, (3568, 2032), "core close", F["tiny_bold"], COLORS["dim"], "mm")


def render() -> Image.Image:
    data = load_data()
    canvas = Image.new("RGB", (W, H), COLORS["bg"])
    draw_background(canvas)
    draw = ImageDraw.Draw(canvas)

    text(draw, (120, 58), "SPACEBIO-BENCH / DETAILED DECK", F["kicker"], COLORS["teal"])
    text(draw, (120, 108), "Core benchmark takeaway", F["title"], COLORS["text"])
    paragraph(
        draw,
        (122, 205),
        "Acts 1-3 establish a mission-held-out benchmark floor and a biology-ready interpretation path.",
        F["subtitle"],
        2300,
        COLORS["muted"],
        6,
    )
    x = 2635
    x = draw_badge(draw, x, 76, "tasks", f"{data['mission_tasks']}", COLORS["blue"])
    x = draw_badge(draw, x, 76, "folds", f"{data['lomo_folds']}", COLORS["teal"])
    x = draw_badge(draw, x, 76, "configs", f"{data['model_configs']}", COLORS["violet"])
    x = draw_badge(draw, x, 76, "DGE rho", f"{data['dge_rho']:.3f}", COLORS["green"])
    draw_badge(draw, x, 76, "bio", f"{data['pathway_concordance']:.3f}", COLORS["orange"])

    draw_core_claim(draw, data)
    draw_claim_status(draw, data)
    draw_bridge(draw)
    draw_footer(draw)

    manifest = {
        "asset": str(OUT_PATH.relative_to(ROOT)),
        "title": "Core benchmark takeaway",
        "source_manifests": [str(path.relative_to(ROOT)) for path in MANIFESTS.values()],
        "metrics": {
            "mission_tasks": int(data["mission_tasks"]),
            "lomo_folds": int(data["lomo_folds"]),
            "model_configs": int(data["model_configs"]),
            "significant_model_rows": int(data["significant_model_rows"]),
            "top_tissue": data["top_tissue"],
            "top_tissue_auroc": float(data["top_tissue_auroc"]),
            "heldout_thymus_auroc": float(data["heldout_thymus"]),
            "heldout_skin_auroc": float(data["heldout_skin"]),
            "dge_log2fc_spearman_mean": float(data["dge_rho"]),
            "pathway_mean_concordance": float(data["pathway_concordance"]),
            "gene_enrichment_ratio": float(data["gene_enrichment"]),
        },
        "readout_status": [
            "benchmark readout supported",
            "model-family lesson supported",
            "biology support layer ready for interpretation",
        ],
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

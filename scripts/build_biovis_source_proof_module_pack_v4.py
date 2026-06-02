#!/usr/bin/env python3
"""Build v0.4 source-proof placeholder modules for SpaceBio-Bench visuals.

v0.4 shifts from biological decoration to evidence-ready frames. Each module
reserves space for a source object, status badges, provenance line, and caveat
before final empirical figures are available.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "assets" / "biovis_source_proof_module_pack_v0_4"
LIGHT_SVG = OUT / "modules" / "svg"
LIGHT_PNG = OUT / "modules" / "png"
DARK_SVG = OUT / "modules_dark" / "svg"
DARK_PNG = OUT / "modules_dark" / "png"
PREVIEW = OUT / "preview"
QA = OUT / "qa"
CREATED = "2026-06-02"

COLORS = {
    "paper": "#F7F4EC",
    "paper2": "#FBFAF7",
    "ink": "#17212B",
    "muted": "#5D6978",
    "rule": "#AEB8C5",
    "blue": "#2D6F9F",
    "sky": "#6BAED6",
    "green": "#178B63",
    "teal": "#1AA090",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "purple": "#6750A4",
    "red": "#B23A3A",
    "deep": "#0D1720",
    "deep2": "#111D27",
    "deep3": "#172636",
    "white": "#F2F6F8",
    "soft": "#B8C4CF",
}


@dataclass(frozen=True)
class ProofModule:
    asset_id: str
    label: str
    module_type: str
    proof_object: str
    intended_use: str
    required_replacement: str
    caution: str


MODULES: list[ProofModule] = [
    ProofModule(
        "source_dataset_record_plate",
        "Source dataset record plate",
        "source-record",
        "dataset accession / record screenshot",
        "OSDR, GeneLab, GEO, or repository record proof",
        "replace placeholder with real accession/source screenshot",
        "do not use as proof unless source id and date are present",
    ),
    ProofModule(
        "expression_matrix_proof_plate",
        "Expression matrix proof plate",
        "matrix-proof",
        "real expression matrix, QC heatmap, or feature table",
        "feature construction and data-processing proof",
        "replace schematic heatmap with real matrix/QC output",
        "do not imply heatmap values are measured while placeholder remains",
    ),
    ProofModule(
        "single_cell_embedding_proof_plate",
        "Single-cell embedding proof plate",
        "single-cell-proof",
        "real embedding or cell-state distribution",
        "single-cell payload and modality proof",
        "replace schematic embedding with real UMAP/t-SNE/spatial map",
        "32px symbols need labels; embedding needs cell count/source metadata",
    ),
    ProofModule(
        "held_out_task_proof_plate",
        "Held-out task proof plate",
        "benchmark-proof",
        "train/test or mission-held-out split diagram",
        "benchmark hygiene, leakage boundary, and evaluation proof",
        "replace placeholder split with actual split manifest or task table",
        "do not use held-out badge without actual withheld evaluation boundary",
    ),
    ProofModule(
        "organoid_evidence_plate",
        "Organoid evidence plate",
        "extension-proof",
        "organoid microscopy, assay panel, or source-derived culture proof",
        "human organoid extension or future validation track",
        "replace placeholder with cited organoid source image or assay output",
        "must keep preliminary/extension caveat unless full benchmark is ready",
    ),
    ProofModule(
        "result_claim_plate",
        "Result claim plate",
        "result-proof",
        "real result plot plus claim-boundary footer",
        "processed result or validated claim surface",
        "replace placeholder plot with final result figure",
        "separate processed result from validated claim",
    ),
]


def c(token: str) -> str:
    return COLORS[token]


def esc(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def palette(dark: bool) -> dict[str, str]:
    if dark:
        return {
            "bg": c("deep"),
            "panel": c("deep2"),
            "panel2": c("deep3"),
            "ink": c("white"),
            "muted": c("soft"),
            "rule": "#617282",
            "slot": "#0F1A24",
            "slot2": "#152435",
        }
    return {
        "bg": "none",
        "panel": c("paper2"),
        "panel2": c("paper"),
        "ink": c("ink"),
        "muted": c("muted"),
        "rule": c("rule"),
        "slot": "#FFFFFF",
        "slot2": "#F1EEE6",
    }


def txt(x: int, y: int, value: str, size: int, color: str, *, weight: int = 400, anchor: str = "start") -> str:
    return (
        f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" fill="{color}" text-anchor="{anchor}" '
        f'letter-spacing="0">{esc(value)}</text>'
    )


def line(x1: int, y1: int, x2: int, y2: int, color: str, *, width: float = 1.4, opacity: float = 1.0, dash: str | None = None) -> str:
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<path d="M{x1} {y1} L{x2} {y2}" stroke="{color}" stroke-width="{width}" '
        f'opacity="{opacity}" stroke-linecap="round"{dash_attr}/>'
    )


def badge(x: int, y: int, label: str, tone: str, pal: dict[str, str]) -> str:
    width = max(110, 42 + len(label) * 10)
    return "\n".join(
        [
            f'<rect x="{x}" y="{y}" width="{width}" height="34" rx="14" fill="{pal["panel"]}" stroke="{c(tone)}" stroke-width="2"/>',
            f'<circle cx="{x + 19}" cy="{y + 17}" r="6" fill="{c(tone)}"/>',
            txt(x + 34, y + 22, label, 13, pal["ink"], weight=700),
        ]
    )


def shell(width: int, height: int, module: ProofModule, pal: dict[str, str], dark: bool, body: str) -> str:
    bg = f'<rect width="{width}" height="{height}" fill="{pal["bg"]}"/>' if dark else '<rect width="100%" height="100%" fill="none"/>'
    meta = json.dumps(
        {
            "kind": "biovis_source_proof_module_v0_4_dark" if dark else "biovis_source_proof_module_v0_4",
            "asset_id": module.asset_id,
            "proof_object": module.proof_object,
            "required_replacement": module.required_replacement,
        },
        sort_keys=True,
    )
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="{esc(module.label)}">
  <metadata>{esc(meta)}</metadata>
  {bg}
  <rect x="18" y="18" width="{width - 36}" height="{height - 36}" rx="24" fill="{pal["panel"]}" stroke="{pal["rule"]}" stroke-width="1.4" opacity="0.98"/>
  {txt(42, 58, module.label, 28, pal["ink"], weight=700)}
  {txt(42, 88, module.intended_use, 14, pal["muted"])}
  {line(42, 112, width - 42, 112, pal["rule"], width=1.2, opacity=0.55)}
{body}
</svg>
'''


def source_dataset_record_plate(module: ProofModule, dark: bool) -> str:
    pal = palette(dark)
    body = [
        f'<rect x="42" y="142" width="670" height="360" rx="18" fill="{pal["slot"]}" stroke="{pal["rule"]}" stroke-width="1.6"/>',
        txt(74, 184, "SOURCE RECORD", 18, pal["muted"], weight=700),
        txt(74, 235, "OSDR / GeneLab / GEO accession", 28, pal["ink"], weight=700),
        txt(74, 280, "Replace with real dataset screenshot", 17, pal["muted"]),
        f'<rect x="74" y="328" width="280" height="34" rx="8" fill="{pal["panel2"]}" stroke="{pal["rule"]}" stroke-width="1"/>',
        txt(92, 351, "source_id: pending", 14, pal["muted"], weight=700),
        f'<rect x="74" y="382" width="590" height="54" rx="10" fill="{pal["panel2"]}" stroke="{pal["rule"]}" stroke-width="1" opacity="0.92"/>',
        txt(94, 416, "sample count / organism / assay / mission date", 16, pal["muted"]),
        f'<rect x="750" y="142" width="290" height="360" rx="18" fill="{pal["panel2"]}" stroke="{pal["rule"]}" stroke-width="1.4"/>',
        txt(778, 184, "Required evidence", 18, pal["ink"], weight=700),
        badge(778, 214, "source", "blue", pal),
        badge(778, 258, "schematic until replaced", "muted", pal),
        badge(778, 302, "caveat", "red", pal),
        txt(778, 378, "Footer must include source URL/date", 14, pal["muted"]),
        txt(778, 414, "No biological claim without record", 14, pal["muted"]),
        txt(42, 548, "Proof object slot: real source record screenshot or accession panel.", 14, pal["muted"]),
    ]
    return shell(1080, 600, module, pal, dark, "\n  ".join(body))


def expression_matrix_proof_plate(module: ProofModule, dark: bool) -> str:
    pal = palette(dark)
    heat = []
    colors = [c("blue"), c("sky"), c("amber"), c("rose"), c("green")]
    for r in range(8):
        for col_i in range(12):
            tone = colors[(r * 3 + col_i * 2) % len(colors)]
            opacity = 0.22 + ((r + col_i) % 5) * 0.10
            heat.append(f'<rect x="{72 + col_i * 38}" y="{162 + r * 32}" width="35" height="29" fill="{tone}" opacity="{opacity:.2f}"/>')
    body = [
        f'<rect x="42" y="142" width="560" height="340" rx="18" fill="{pal["slot"]}" stroke="{pal["rule"]}" stroke-width="1.6"/>',
        *heat,
        txt(72, 438, "Placeholder heatmap. Replace with real matrix/QC.", 13, pal["muted"]),
        f'<rect x="630" y="142" width="408" height="340" rx="18" fill="{pal["panel2"]}" stroke="{pal["rule"]}" stroke-width="1.4"/>',
        txt(662, 184, "Matrix proof checklist", 18, pal["ink"], weight=700),
        badge(662, 214, "processed", "teal", pal),
        badge(662, 258, "checksum", "green", pal),
        badge(662, 302, "source", "blue", pal),
        txt(662, 366, "Needs feature universe, sample axis,", 14, pal["muted"]),
        txt(662, 392, "normalization, and frozen input hash.", 14, pal["muted"]),
        txt(42, 548, "Use only after replacing schematic values with real output.", 14, pal["muted"]),
    ]
    return shell(1080, 600, module, pal, dark, "\n  ".join(body))


def single_cell_embedding_proof_plate(module: ProofModule, dark: bool) -> str:
    pal = palette(dark)
    dots = []
    clusters = [
        (170, 235, c("blue")),
        (260, 310, c("green")),
        (355, 230, c("rose")),
        (430, 330, c("amber")),
        (250, 205, c("purple")),
    ]
    for idx, (cx, cy, tone) in enumerate(clusters):
        dots.append(f'<ellipse cx="{cx}" cy="{cy}" rx="90" ry="46" fill="{tone}" opacity="0.13"/>')
        for j in range(18):
            x = cx + ((j * 37 + idx * 17) % 120) - 60
            y = cy + ((j * 23 + idx * 11) % 66) - 33
            dots.append(f'<circle cx="{x}" cy="{y}" r="4.4" fill="{tone}" opacity="0.76"/>')
    body = [
        f'<rect x="42" y="142" width="560" height="340" rx="18" fill="{pal["slot"]}" stroke="{pal["rule"]}" stroke-width="1.6"/>',
        *dots,
        txt(72, 438, "Replace with real embedding or spatial/cell-state map.", 13, pal["muted"]),
        f'<rect x="630" y="142" width="408" height="340" rx="18" fill="{pal["panel2"]}" stroke="{pal["rule"]}" stroke-width="1.4"/>',
        txt(662, 184, "Required labels", 18, pal["ink"], weight=700),
        txt(662, 232, "species / tissue / assay", 17, pal["ink"], weight=700),
        txt(662, 270, "cell count / filtering", 17, pal["ink"], weight=700),
        txt(662, 308, "batch or mission split", 17, pal["ink"], weight=700),
        badge(662, 356, "prelim", "amber", pal),
        badge(662, 400, "source", "blue", pal),
        txt(42, 548, "Use as modality proof, not as result claim until real embedding is inserted.", 14, pal["muted"]),
    ]
    return shell(1080, 600, module, pal, dark, "\n  ".join(body))


def held_out_task_proof_plate(module: ProofModule, dark: bool) -> str:
    pal = palette(dark)
    body = [
        f'<rect x="42" y="142" width="670" height="340" rx="18" fill="{pal["slot"]}" stroke="{pal["rule"]}" stroke-width="1.6"/>',
        f'<rect x="92" y="218" width="220" height="126" rx="18" fill="{pal["panel2"]}" stroke="{c("blue")}" stroke-width="2.5"/>',
        txt(202, 276, "TRAIN", 26, pal["ink"], weight=700, anchor="middle"),
        txt(202, 310, "fit / select / tune", 14, pal["muted"], anchor="middle"),
        f'<rect x="442" y="218" width="220" height="126" rx="18" fill="{pal["panel2"]}" stroke="{c("green")}" stroke-width="2.5"/>',
        txt(552, 276, "HELD-OUT", 26, pal["ink"], weight=700, anchor="middle"),
        txt(552, 310, "evaluate only", 14, pal["muted"], anchor="middle"),
        line(328, 281, 426, 281, c("red"), width=4.0, opacity=0.95, dash="8 10"),
        txt(377, 258, "no leakage", 14, c("red"), weight=700, anchor="middle"),
        f'<rect x="750" y="142" width="290" height="340" rx="18" fill="{pal["panel2"]}" stroke="{pal["rule"]}" stroke-width="1.4"/>',
        txt(778, 184, "Split proof", 18, pal["ink"], weight=700),
        badge(778, 214, "train-only", "purple", pal),
        badge(778, 258, "held-out", "blue", pal),
        badge(778, 302, "caveat", "red", pal),
        txt(778, 376, "Needs split manifest and", 14, pal["muted"]),
        txt(778, 402, "task construction timestamp.", 14, pal["muted"]),
        txt(42, 548, "Use only with actual train/test or mission-held-out manifest.", 14, pal["muted"]),
    ]
    return shell(1080, 600, module, pal, dark, "\n  ".join(body))


def organoid_evidence_plate(module: ProofModule, dark: bool) -> str:
    pal = palette(dark)
    cells = []
    for idx, (cx, cy, tone) in enumerate([(190, 250, c("purple")), (280, 235, c("green")), (250, 330, c("sky")), (360, 310, c("rose"))]):
        cells.append(f'<circle cx="{cx}" cy="{cy}" r="{58 + idx * 6}" fill="{tone}" opacity="0.20"/>')
        cells.append(f'<circle cx="{cx}" cy="{cy}" r="17" fill="{tone}" opacity="0.46"/>')
    body = [
        f'<rect x="42" y="142" width="560" height="340" rx="18" fill="{pal["slot"]}" stroke="{pal["rule"]}" stroke-width="1.6"/>',
        *cells,
        txt(72, 438, "Replace with cited organoid image, assay, or source panel.", 13, pal["muted"]),
        f'<rect x="630" y="142" width="408" height="340" rx="18" fill="{pal["panel2"]}" stroke="{pal["rule"]}" stroke-width="1.4"/>',
        txt(662, 184, "Organoid extension status", 18, pal["ink"], weight=700),
        badge(662, 214, "prelim", "amber", pal),
        badge(662, 258, "source", "blue", pal),
        badge(662, 302, "caveat", "red", pal),
        txt(662, 366, "Required: culture model, assay,", 14, pal["muted"]),
        txt(662, 392, "source/citation, and benchmark role.", 14, pal["muted"]),
        txt(42, 548, "Keep extension caveat unless organoid benchmark payload is complete.", 14, pal["muted"]),
    ]
    return shell(1080, 600, module, pal, dark, "\n  ".join(body))


def result_claim_plate(module: ProofModule, dark: bool) -> str:
    pal = palette(dark)
    plot = [
        line(88, 424, 548, 424, pal["rule"], width=1.4, opacity=0.8),
        line(88, 188, 88, 424, pal["rule"], width=1.4, opacity=0.8),
        '<path d="M110 386 C180 330, 218 355, 268 286 C320 214, 392 247, 526 184" fill="none" stroke="' + c("green") + '" stroke-width="5" stroke-linecap="round"/>',
        '<path d="M110 408 C186 390, 234 362, 290 336 C368 298, 430 286, 526 260" fill="none" stroke="' + c("blue") + '" stroke-width="4" stroke-linecap="round" opacity="0.76"/>',
    ]
    body = [
        f'<rect x="42" y="142" width="560" height="340" rx="18" fill="{pal["slot"]}" stroke="{pal["rule"]}" stroke-width="1.6"/>',
        *plot,
        txt(100, 174, "RESULT PLOT SLOT", 15, pal["muted"], weight=700),
        txt(72, 460, "Replace with final result; keep processed vs validated labels.", 13, pal["muted"]),
        f'<rect x="630" y="142" width="408" height="340" rx="18" fill="{pal["panel2"]}" stroke="{pal["rule"]}" stroke-width="1.4"/>',
        txt(662, 184, "Claim boundary", 18, pal["ink"], weight=700),
        badge(662, 214, "processed", "teal", pal),
        badge(662, 258, "validated", "green", pal),
        badge(662, 302, "caveat", "red", pal),
        txt(662, 366, "Validated claim requires", 14, pal["muted"]),
        txt(662, 392, "reviewed result and source trace.", 14, pal["muted"]),
        txt(42, 548, "Do not upgrade processed result to validated claim without review.", 14, pal["muted"]),
    ]
    return shell(1080, 600, module, pal, dark, "\n  ".join(body))


BUILDERS = {
    "source_dataset_record_plate": source_dataset_record_plate,
    "expression_matrix_proof_plate": expression_matrix_proof_plate,
    "single_cell_embedding_proof_plate": single_cell_embedding_proof_plate,
    "held_out_task_proof_plate": held_out_task_proof_plate,
    "organoid_evidence_plate": organoid_evidence_plate,
    "result_claim_plate": result_claim_plate,
}


def ensure_dirs() -> None:
    for directory in [OUT, LIGHT_SVG, LIGHT_PNG, DARK_SVG, DARK_PNG, PREVIEW, QA]:
        directory.mkdir(parents=True, exist_ok=True)


def run_rsvg(svg_path: Path, png_path: Path, *, width: int = 1800) -> bool:
    exe = shutil.which("rsvg-convert")
    if not exe:
        return False
    subprocess.run([exe, "-w", str(width), "-o", str(png_path), str(svg_path)], check=True)
    return True


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def font(size: int, *, bold: bool = False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def paste_fit(canvas: Image.Image, image_path: Path, box: tuple[int, int, int, int]) -> None:
    image = Image.open(image_path).convert("RGBA")
    x, y, w, h = box
    scale = min(w / image.width, h / image.height)
    size = (max(1, int(image.width * scale)), max(1, int(image.height * scale)))
    image = image.resize(size, Image.Resampling.LANCZOS)
    canvas.alpha_composite(image, (x + (w - image.width) // 2, y + (h - image.height) // 2))


def draw_wrapped(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    value: str,
    font_obj,
    fill: tuple[int, int, int],
    *,
    max_width: int,
    line_height: int,
    max_lines: int = 3,
) -> int:
    words = value.split()
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if draw.textlength(candidate, font=font_obj) <= max_width:
            current = candidate
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    lines = lines[:max_lines]
    x, y = xy
    for idx, line_value in enumerate(lines):
        draw.text((x, y + idx * line_height), line_value, font=font_obj, fill=fill)
    return y + len(lines) * line_height


def grayscale_copy(src: Path, dst: Path) -> bool:
    with Image.open(src) as image:
        image.convert("LA").convert("RGBA").save(dst)
    return True


def contact_sheet(path: Path, *, dark: bool = False) -> bool:
    width, height = 2400, 1760
    bg = hex_to_rgb(c("deep")) if dark else hex_to_rgb(c("paper"))
    fg = hex_to_rgb(c("white")) if dark else hex_to_rgb(c("ink"))
    muted = hex_to_rgb(c("soft")) if dark else hex_to_rgb(c("muted"))
    card = hex_to_rgb(c("deep2")) if dark else hex_to_rgb(c("paper2"))
    rule = hex_to_rgb("#617282" if dark else c("rule"))
    canvas = Image.new("RGBA", (width, height), bg + (255,))
    draw = ImageDraw.Draw(canvas)
    title = font(40, bold=True)
    body = font(19)
    label = font(24, bold=True)
    small = font(15)
    draw.text((72, 44), f"Source-proof placeholder modules v0.4 {'dark' if dark else 'light'}", font=title, fill=fg)
    draw.text((72, 96), "Evidence-ready frames: proof slot, status badges, source line, and caveat position.", font=body, fill=muted)
    png_root = DARK_PNG if dark else LIGHT_PNG
    for idx, module in enumerate(MODULES):
        col_i = idx % 2
        row_i = idx // 2
        x = 72 + col_i * 1140
        y = 168 + row_i * 470
        draw.rounded_rectangle((x, y, x + 1040, y + 405), radius=18, fill=card, outline=rule, width=1)
        paste_fit(canvas, png_root / f"{module.asset_id}.png", (x + 28, y + 24, 710, 395))
        draw_wrapped(draw, (x + 770, y + 50), module.label, label, fg, max_width=235, line_height=28, max_lines=2)
        draw.text((x + 770, y + 110), module.module_type, font=small, fill=muted)
        draw.text((x + 770, y + 154), "Proof object:", font=small, fill=fg)
        next_y = draw_wrapped(draw, (x + 770, y + 180), module.proof_object, small, muted, max_width=240, line_height=22, max_lines=3)
        draw.text((x + 770, max(next_y + 24, y + 248)), "Replace:", font=small, fill=fg)
        draw_wrapped(draw, (x + 770, max(next_y + 50, y + 274)), module.required_replacement, small, muted, max_width=240, line_height=22, max_lines=3)
    canvas.convert("RGB").save(path)
    return True


def write_manifest(conversion_ok: bool) -> None:
    records = []
    for module in MODULES:
        for variant, svg_root, png_root in [
            ("light", LIGHT_SVG, LIGHT_PNG),
            ("dark", DARK_SVG, DARK_PNG),
        ]:
            records.append(
                {
                    "asset_id": f"{module.asset_id}_{variant}",
                    "base_asset_id": module.asset_id,
                    "asset_type": "source_proof_placeholder_module",
                    "variant": variant,
                    "label": module.label,
                    "module_type": module.module_type,
                    "proof_object": module.proof_object,
                    "intended_use": module.intended_use,
                    "required_replacement": module.required_replacement,
                    "caution": module.caution,
                    "svg": rel(svg_root / f"{module.asset_id}.svg"),
                    "png_preview": rel(png_root / f"{module.asset_id}.png"),
                    "license": "Original project artwork; repository MIT license.",
                }
            )
    manifest = {
        "created": CREATED,
        "version": "0.4",
        "scope": "Source-proof placeholder modules for evidence-ready SpaceBio-Bench slides.",
        "design_rule": "Reserve proof object, status badge, source line, and caveat before final figures are available.",
        "asset_counts": {
            "modules_base": len(MODULES),
            "variants": 2,
            "total": len(records),
        },
        "conversion_ok": conversion_ok,
        "preview": {
            "light_contact_sheet": rel(PREVIEW / "biovis_source_proof_modules_v0_4_light.png"),
            "dark_contact_sheet": rel(PREVIEW / "biovis_source_proof_modules_v0_4_dark.png"),
            "light_grayscale_qa": rel(QA / "biovis_source_proof_modules_v0_4_light_grayscale.png"),
            "dark_grayscale_qa": rel(QA / "biovis_source_proof_modules_v0_4_dark_grayscale.png"),
        },
        "records": records,
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with (OUT / "manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)
    qa = {
        "created": CREATED,
        "automatic_checks": {
            "expected_svg_count": len(records),
            "svg_count": len(list(LIGHT_SVG.glob("*.svg"))) + len(list(DARK_SVG.glob("*.svg"))),
            "png_count": len(list(LIGHT_PNG.glob("*.png"))) + len(list(DARK_PNG.glob("*.png"))),
            "rsvg_conversion_available": conversion_ok,
            "manifest_json_exists": (OUT / "manifest.json").exists(),
            "manifest_csv_exists": (OUT / "manifest.csv").exists(),
            "contact_sheets_exist": (PREVIEW / "biovis_source_proof_modules_v0_4_light.png").exists()
            and (PREVIEW / "biovis_source_proof_modules_v0_4_dark.png").exists(),
            "grayscale_qa_exists": (QA / "biovis_source_proof_modules_v0_4_light_grayscale.png").exists()
            and (QA / "biovis_source_proof_modules_v0_4_dark_grayscale.png").exists(),
        },
        "manual_review": {
            "initial_verdict": "pass with conditions",
            "proof_slot": "pass: each module reserves a visible source/result/proof-object area before final evidence is inserted",
            "claim_boundary": "pass: status badges, source requirements, and caveat language remain attached to the proof surface",
            "density": "conditional pass: one large proof module is appropriate for main slides; multiple proof modules should be review-board or appendix material",
            "visual_notes": "dark variants feel more premium and presentation-ready; light variants are clear but read more like review or appendix boards",
            "next_step": "replace placeholders with actual OSDR/GeneLab/GEO record screenshots, matrix/QC outputs, embeddings, split manifests, and result plots",
        },
    }
    (QA / "biovis_source_proof_modules_v0_4_qa.json").write_text(json.dumps(qa, indent=2) + "\n")


def main() -> None:
    ensure_dirs()
    conversion_ok = shutil.which("rsvg-convert") is not None
    for module in MODULES:
        light_svg = BUILDERS[module.asset_id](module, False)
        dark_svg = BUILDERS[module.asset_id](module, True)
        (LIGHT_SVG / f"{module.asset_id}.svg").write_text(light_svg)
        (DARK_SVG / f"{module.asset_id}.svg").write_text(dark_svg)
        if conversion_ok:
            run_rsvg(LIGHT_SVG / f"{module.asset_id}.svg", LIGHT_PNG / f"{module.asset_id}.png")
            run_rsvg(DARK_SVG / f"{module.asset_id}.svg", DARK_PNG / f"{module.asset_id}.png")
    contact_sheet(PREVIEW / "biovis_source_proof_modules_v0_4_light.png", dark=False)
    contact_sheet(PREVIEW / "biovis_source_proof_modules_v0_4_dark.png", dark=True)
    grayscale_copy(PREVIEW / "biovis_source_proof_modules_v0_4_light.png", QA / "biovis_source_proof_modules_v0_4_light_grayscale.png")
    grayscale_copy(PREVIEW / "biovis_source_proof_modules_v0_4_dark.png", QA / "biovis_source_proof_modules_v0_4_dark_grayscale.png")
    write_manifest(conversion_ok)
    print(json.dumps({"output": rel(OUT), "modules": len(MODULES), "variants": 2, "conversion_ok": conversion_ok}, indent=2))


if __name__ == "__main__":
    main()

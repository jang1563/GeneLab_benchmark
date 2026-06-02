#!/usr/bin/env python3
"""Build premium biological motif assets for SpaceBio-Bench slides.

These are not simple icons. They are transparent SVG motif plates designed to
sit on top of scientific slide scenes as cell, molecular, tissue, organ, and
assay signposts without replacing source-derived evidence.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "assets" / "biovis_motif_pack_v0_2"
SVG_DIR = OUT / "svg"
PNG_DIR = OUT / "png"

W = 960
H = 640

COLORS = {
    "ink": "#17212B",
    "paper": "#FBFAF6",
    "muted": "#5D6978",
    "rule": "#AEB8C5",
    "source_blue": "#286EA6",
    "bio_green": "#15815E",
    "assay_teal": "#159A8A",
    "label_amber": "#A36F13",
    "test_red": "#B91C1C",
    "model_purple": "#6D3EDB",
    "shadow": "#9C9487",
}


@dataclass(frozen=True)
class Motif:
    motif_id: str
    category: str
    label: str
    semantic_use: str
    primary_color: str
    secondary_color: str
    caution: str


MOTIFS = [
    Motif("cell_state_microfield", "cellular", "Cell state microfield", "single-cell state, cell composition, sample biology", "bio_green", "source_blue", "Use as cellular context, not a specific cell type."),
    Motif("mitochondrial_stress_field", "cellular", "Mitochondrial stress field", "organelle stress, metabolism, radiation or oxidative context", "assay_teal", "label_amber", "Use as biological interpretation only unless mitochondria are measured."),
    Motif("organoid_spheroid_field", "cellular", "Organoid spheroid field", "human organoid or 3D culture extension", "model_purple", "bio_green", "Pair with draft/diagnostic caveats for organoid tracks."),
    Motif("single_cell_droplet_field", "cellular", "Single-cell droplet field", "single-cell assay, cell barcode, droplet capture", "source_blue", "assay_teal", "Do not imply runnable scRNA-seq payload readiness."),
    Motif("tissue_contour_section", "tissue", "Tissue contour section", "bulk tissue, histology-like context, tissue-specific result", "bio_green", "rule", "Symbolic tissue surface, not microscopy proof."),
    Motif("epithelial_barrier_field", "tissue", "Epithelial barrier field", "barrier tissue, tissue architecture, lining", "assay_teal", "bio_green", "Generic barrier architecture only."),
    Motif("vascular_microsection", "tissue", "Vascular microsection", "vascular/perfusion/cfRNA origin context", "test_red", "source_blue", "Use red carefully so it is not confused with held-out/test status."),
    Motif("spatial_spot_tissue_field", "tissue", "Spatial spot tissue field", "spatial transcriptomics, tissue spots, location-aware biology", "bio_green", "model_purple", "Use only for spatial/modality explanation or clearly bounded future work."),
    Motif("rna_expression_trace", "molecular", "RNA expression trace", "RNA-seq expression, molecular readout, transcript signal", "assay_teal", "source_blue", "Default molecular signal motif for expression assays."),
    Motif("dna_context_trace", "molecular", "DNA context trace", "gene set, sequence context, feature identity", "source_blue", "model_purple", "Avoid using as RNA-seq expression shorthand unless caption clarifies."),
    Motif("protein_complex_field", "molecular", "Protein complex field", "protein, proteomics, pathway protein-level intuition", "label_amber", "bio_green", "Do not imply direct proteomics without data support."),
    Motif("pathway_program_field", "molecular", "Pathway program field", "pathway summary, biological program, network feature", "bio_green", "label_amber", "Use for summarized programs, not causal pathways by default."),
    Motif("brain_neural_context", "organ", "Brain neural context", "brain, neural organoid, CNS interpretation", "model_purple", "source_blue", "Label clearly if the asset represents organoid rather than organ tissue."),
    Motif("kidney_tubule_context", "organ", "Kidney tubule context", "kidney, renal tissue, tubule-like biology", "bio_green", "label_amber", "Use with tissue labels for specificity."),
    Motif("liver_lobule_context", "organ", "Liver lobule context", "liver, metabolism, lobule-like tissue context", "label_amber", "test_red", "Keep amber/red separated from test-boundary semantics."),
    Motif("sample_to_matrix_assay", "assay", "Sample to matrix assay", "sample collection to expression matrix/task input", "source_blue", "assay_teal", "Do not imply wet-lab protocol steps not in the analysis."),
]


def c(name: str) -> str:
    return COLORS[name]


def ensure_dirs() -> None:
    SVG_DIR.mkdir(parents=True, exist_ok=True)
    PNG_DIR.mkdir(parents=True, exist_ok=True)


def defs(motif: Motif) -> str:
    p = c(motif.primary_color)
    s = c(motif.secondary_color)
    return f"""
  <defs>
    <filter id="soft_shadow" x="-20%" y="-20%" width="140%" height="140%">
      <feDropShadow dx="12" dy="16" stdDeviation="14" flood-color="{c('shadow')}" flood-opacity="0.16"/>
    </filter>
    <radialGradient id="bio_glow" cx="46%" cy="42%" r="64%">
      <stop offset="0%" stop-color="{p}" stop-opacity="0.18"/>
      <stop offset="58%" stop-color="{s}" stop-opacity="0.08"/>
      <stop offset="100%" stop-color="{c('paper')}" stop-opacity="0"/>
    </radialGradient>
    <linearGradient id="surface_fade" x1="0%" y1="0%" x2="100%" y2="100%">
      <stop offset="0%" stop-color="{c('paper')}" stop-opacity="0.78"/>
      <stop offset="100%" stop-color="{c('paper')}" stop-opacity="0.16"/>
    </linearGradient>
  </defs>"""


def field_layer(motif: Motif) -> str:
    p = c(motif.primary_color)
    s = c(motif.secondary_color)
    rule = c("rule")
    lines = [
        '<rect width="960" height="640" fill="none"/>',
        f'<ellipse cx="482" cy="322" rx="356" ry="218" fill="url(#bio_glow)" opacity="0.92"/>',
        f'<path d="M104 318 H856" stroke="{rule}" stroke-width="1.8" opacity="0.18"/>',
        f'<path d="M188 168 V510 M480 114 V546 M774 178 V498" stroke="{rule}" stroke-width="1.4" opacity="0.12"/>',
        f'<path d="M154 426 C305 358, 607 356, 812 432" fill="none" stroke="{p}" stroke-width="2.0" opacity="0.12" stroke-linecap="round"/>',
        f'<path d="M158 230 C330 292, 620 290, 800 218" fill="none" stroke="{s}" stroke-width="1.8" opacity="0.10" stroke-linecap="round"/>',
    ]
    return "\n".join(f"  {line}" for line in lines)


def circle(x: int, y: int, r: int, fill: str, opacity: float, stroke: str = "", sw: float = 0) -> str:
    stroke_attr = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
    return f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" opacity="{opacity:.2f}"{stroke_attr}/>'


def ellipse(x: int, y: int, rx: int, ry: int, fill: str, opacity: float, stroke: str = "", sw: float = 0) -> str:
    stroke_attr = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
    return f'<ellipse cx="{x}" cy="{y}" rx="{rx}" ry="{ry}" fill="{fill}" opacity="{opacity:.2f}"{stroke_attr}/>'


def motif_body(motif: Motif) -> str:
    p = c(motif.primary_color)
    s = c(motif.secondary_color)
    ink = c("ink")
    muted = c("muted")
    rule = c("rule")
    paper = c("paper")
    aid = motif.motif_id
    body: list[str] = []
    common = 'stroke-linecap="round" stroke-linejoin="round"'

    if aid == "cell_state_microfield":
        body += [
            f'<g filter="url(#soft_shadow)">',
            ellipse(405, 310, 176, 114, paper, 0.88, p, 8),
            ellipse(570, 345, 152, 96, paper, 0.78, s, 6),
            circle(386, 315, 42, s, 0.78),
            circle(565, 346, 34, p, 0.56),
            f'<path d="M270 314 C360 255, 540 262, 690 354" fill="none" stroke="{rule}" stroke-width="5" opacity="0.36" {common}/>',
            f'</g>',
        ]
        for x, y, r in [(285, 255, 9), (472, 236, 7), (680, 292, 8), (636, 420, 6), (336, 398, 6)]:
            body.append(circle(x, y, r, p, 0.30))
    elif aid == "mitochondrial_stress_field":
        body += [
            f'<g filter="url(#soft_shadow)">',
            f'<path d="M280 374 C210 272, 310 168, 472 194 C636 220, 721 328, 636 427 C542 535, 360 492, 280 374 Z" fill="{paper}" stroke="{p}" stroke-width="9" {common}/>',
            f'<path d="M310 354 C388 248, 420 460, 498 342 C558 252, 584 418, 660 314" fill="none" stroke="{s}" stroke-width="12" opacity="0.82" {common}/>',
            f'<path d="M316 268 C420 202, 570 220, 650 300" fill="none" stroke="{rule}" stroke-width="5" opacity="0.30" {common}/>',
            f'</g>',
            f'<path d="M672 204 l22 -28 M706 258 l34 -4 M660 250 l20 26" stroke="{c("label_amber")}" stroke-width="6" opacity="0.70" {common}/>',
        ]
    elif aid == "organoid_spheroid_field":
        body.append(f'<g filter="url(#soft_shadow)">')
        body.append(circle(486, 320, 156, paper, 0.92, p, 9))
        cells = [(414, 276, 54, s, 0.38), (486, 252, 50, p, 0.30), (556, 284, 56, s, 0.42), (416, 360, 48, p, 0.34), (512, 360, 68, s, 0.36), (580, 370, 40, p, 0.26)]
        for x, y, r, col, op in cells:
            body.append(circle(x, y, r, col, op))
        body.append(f'<path d="M332 320 C420 214, 573 220, 658 322 C574 450, 414 448, 332 320 Z" fill="none" stroke="{rule}" stroke-width="5" opacity="0.30" {common}/>')
        body.append('</g>')
    elif aid == "single_cell_droplet_field":
        body.append(f'<g filter="url(#soft_shadow)">')
        for idx, (x, y) in enumerate([(322, 260), (410, 220), (505, 250), (598, 218), (638, 314), (546, 370), (438, 356), (348, 330)]):
            col = p if idx % 2 == 0 else s
            body.append(circle(x, y, 39, paper, 0.86, col, 6))
            body.append(circle(x - 7, y + 3, 12, col, 0.55))
        body.append(f'<path d="M280 422 C405 456, 584 456, 704 412" fill="none" stroke="{rule}" stroke-width="5" opacity="0.32" {common}/>')
        body.append('</g>')
    elif aid == "tissue_contour_section":
        body += [
            f'<g filter="url(#soft_shadow)">',
            f'<path d="M252 340 C280 182, 451 138, 618 184 C760 224, 756 394, 624 474 C452 578, 210 500, 252 340 Z" fill="{paper}" stroke="{p}" stroke-width="9" opacity="0.92" {common}/>',
            f'<path d="M304 314 C418 238, 585 254, 704 338" fill="none" stroke="{rule}" stroke-width="6" opacity="0.34" {common}/>',
            f'<path d="M318 402 C430 340, 574 352, 662 422" fill="none" stroke="{s}" stroke-width="6" opacity="0.42" {common}/>',
            f'</g>',
        ]
        for x, y, r, col in [(410, 290, 13, p), (520, 334, 12, s), (608, 282, 10, p), (464, 420, 10, s)]:
            body.append(circle(x, y, r, col, 0.42))
    elif aid == "epithelial_barrier_field":
        body.append(f'<g filter="url(#soft_shadow)">')
        body.append(f'<path d="M246 420 H714" stroke="{ink}" stroke-width="8" opacity="0.44" {common}/>')
        for x, col, h in [(330, p, 150), (430, s, 190), (530, p, 166), (630, s, 142)]:
            body.append(f'<path d="M{x - 48} 410 C{x - 42} {410 - h}, {x + 42} {410 - h}, {x + 48} 410 Z" fill="{paper}" stroke="{col}" stroke-width="8" opacity="0.92" {common}/>')
            body.append(circle(x, 356, 16, col, 0.48))
        body.append(f'<path d="M286 454 C420 478, 568 478, 704 452" fill="none" stroke="{rule}" stroke-width="4" opacity="0.30" {common}/>')
        body.append('</g>')
    elif aid == "vascular_microsection":
        body += [
            f'<g filter="url(#soft_shadow)">',
            circle(490, 324, 156, paper, 0.90, p, 9),
            circle(490, 324, 92, "none", 1.0, s, 16),
            circle(490, 324, 42, paper, 0.86, muted, 6),
            f'<path d="M382 248 C460 186, 584 202, 642 284" fill="none" stroke="{rule}" stroke-width="5" opacity="0.30" {common}/>',
            f'</g>',
            f'<path d="M350 324 H426 M554 324 H630" stroke="{s}" stroke-width="5" opacity="0.52" {common}/>',
        ]
    elif aid == "spatial_spot_tissue_field":
        body.append(f'<g filter="url(#soft_shadow)">')
        body.append(f'<path d="M298 332 C326 198, 488 176, 636 222 C724 250, 714 390, 596 446 C456 512, 262 468, 298 332 Z" fill="{paper}" stroke="{p}" stroke-width="8" opacity="0.86" {common}/>')
        for row, y in enumerate([250, 306, 362, 418]):
            for col, x in enumerate([360, 430, 500, 570, 640]):
                op = 0.32 if (row + col) % 2 else 0.52
                body.append(circle(x, y, 12, s, op))
        body.append(f'<path d="M330 386 C430 316, 560 318, 675 384" fill="none" stroke="{rule}" stroke-width="4" opacity="0.36" {common}/>')
        body.append('</g>')
    elif aid == "rna_expression_trace":
        body.append(f'<g filter="url(#soft_shadow)">')
        body.append(f'<path d="M188 352 C290 154, 430 494, 554 270 C620 150, 704 182, 770 270" fill="none" stroke="{p}" stroke-width="14" {common}/>')
        for x, y, dx, dy in [(268, 294, 58, -38), (372, 316, 54, 36), (488, 332, 58, -44), (612, 244, 58, 34), (704, 246, 54, -28)]:
            body.append(f'<path d="M{x} {y} l{dx} {dy}" stroke="{s}" stroke-width="8" opacity="0.76" {common}/>')
            body.append(circle(x + dx, y + dy, 9, s, 0.78))
        body.append('</g>')
    elif aid == "dna_context_trace":
        body.append(f'<g filter="url(#soft_shadow)">')
        body.append(f'<path d="M328 130 C632 246, 632 394, 328 510" fill="none" stroke="{p}" stroke-width="12" {common}/>')
        body.append(f'<path d="M632 130 C328 246, 328 394, 632 510" fill="none" stroke="{s}" stroke-width="12" {common}/>')
        for y, x1, x2 in [(184, 382, 578), (238, 336, 624), (292, 360, 600), (346, 360, 600), (400, 336, 624), (454, 382, 578)]:
            body.append(f'<path d="M{x1} {y} H{x2}" stroke="{rule}" stroke-width="7" opacity="0.56" {common}/>')
        body.append('</g>')
    elif aid == "protein_complex_field":
        body.append(f'<g filter="url(#soft_shadow)">')
        for x, y, r, col, op in [(390, 286, 70, p, 0.68), (496, 260, 60, s, 0.56), (486, 372, 78, p, 0.40), (598, 352, 52, s, 0.48), (610, 250, 38, p, 0.36)]:
            body.append(circle(x, y, r, col, op))
        body.append(f'<path d="M330 210 C430 120, 642 150, 702 296" fill="none" stroke="{ink}" stroke-width="6" opacity="0.20" {common}/>')
        body.append('</g>')
    elif aid == "pathway_program_field":
        body.append(f'<g filter="url(#soft_shadow)">')
        nodes = [(300, 318, 42, p), (450, 220, 42, s), (450, 416, 42, s), (620, 318, 46, p), (725, 244, 32, s), (725, 392, 32, s)]
        for x, y, r, col in nodes:
            body.append(circle(x, y, r, paper, 0.92, col, 8))
        for x1, y1, x2, y2 in [(342, 302, 410, 246), (342, 334, 410, 390), (492, 246, 580, 302), (492, 390, 580, 334), (666, 300, 696, 258), (666, 336, 696, 378), (342, 318, 574, 318)]:
            body.append(f'<path d="M{x1} {y1} L{x2} {y2}" stroke="{rule}" stroke-width="5" opacity="0.62" {common}/>')
        body.append('</g>')
    elif aid == "brain_neural_context":
        body.append(f'<g filter="url(#soft_shadow)">')
        body.append(f'<path d="M306 360 C186 276, 240 146, 384 166 C444 70, 616 98, 638 220 C770 238, 780 398, 658 458 C618 574, 414 560, 350 450 C322 452, 298 414, 306 360 Z" fill="{paper}" stroke="{p}" stroke-width="8" opacity="0.92" {common}/>')
        for d in ["M398 194 C348 274, 392 334, 488 322", "M520 150 C468 244, 536 298, 636 288", "M372 426 C496 380, 594 408, 656 482"]:
            body.append(f'<path d="{d}" fill="none" stroke="{s}" stroke-width="5" opacity="0.60" {common}/>')
        body.append('</g>')
    elif aid == "kidney_tubule_context":
        body.append(f'<g filter="url(#soft_shadow)">')
        body.append(f'<path d="M554 146 C390 104, 282 230, 284 370 C286 526, 456 586, 610 500 C720 438, 720 306, 632 262 C576 234, 580 190, 628 166 C612 158, 588 152, 554 146 Z" fill="{paper}" stroke="{p}" stroke-width="9" opacity="0.92" {common}/>')
        body.append(f'<path d="M584 244 C494 312, 506 430, 612 482" fill="none" stroke="{s}" stroke-width="12" opacity="0.70" {common}/>')
        body.append(f'<path d="M604 358 H764" stroke="{s}" stroke-width="7" opacity="0.62" {common}/>')
        body.append('</g>')
    elif aid == "liver_lobule_context":
        body.append(f'<g filter="url(#soft_shadow)">')
        body.append(f'<path d="M480 148 L690 268 V438 L480 548 L270 438 V268 Z" fill="{paper}" stroke="{p}" stroke-width="9" opacity="0.92" {common}/>')
        body.append(circle(480, 350, 42, s, 0.72))
        for x2, y2 in [(480, 198), (480, 500), (318, 282), (642, 282), (318, 418), (642, 418)]:
            body.append(f'<path d="M480 350 L{x2} {y2}" stroke="{rule}" stroke-width="5" opacity="0.58" {common}/>')
        body.append('</g>')
    elif aid == "sample_to_matrix_assay":
        body.append(f'<g filter="url(#soft_shadow)">')
        body.append(f'<path d="M246 158 H356 M276 158 V410 C276 486, 326 502, 382 466 C412 446, 416 420, 416 380 V158" fill="{paper}" stroke="{p}" stroke-width="8" {common}/>')
        body.append(f'<path d="M282 346 C324 320, 362 384, 410 346 V424 C410 478, 282 478, 282 424 Z" fill="{s}" opacity="0.55"/>')
        body.append(f'<path d="M458 320 H540" stroke="{rule}" stroke-width="6" opacity="0.62" {common}/>')
        body.append(f'<path d="M520 260 L566 320 L520 380" fill="none" stroke="{rule}" stroke-width="6" opacity="0.62" {common}/>')
        body.append(f'<rect x="590" y="210" width="180" height="220" fill="{paper}" stroke="{muted}" stroke-width="7" opacity="0.92"/>')
        for x in [626, 662, 698, 734]:
            body.append(f'<path d="M{x} 210 V430" stroke="{rule}" stroke-width="3" opacity="0.45"/>')
        for y in [254, 298, 342, 386]:
            body.append(f'<path d="M590 {y} H770" stroke="{rule}" stroke-width="3" opacity="0.45"/>')
        for x, y in [(626, 254), (698, 298), (734, 342), (662, 386)]:
            body.append(f'<rect x="{x - 11}" y="{y - 11}" width="22" height="22" rx="4" fill="{c("model_purple")}" opacity="0.62"/>')
        body.append('</g>')
    else:
        raise ValueError(aid)
    return "\n".join(f"  {line}" for line in body)


def svg_for(motif: Motif) -> str:
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-label="{motif.label}">
  <metadata>SpaceBio-Bench original premium biological motif asset v0.2; motif_id={motif.motif_id}; category={motif.category}; license=repository MIT license.</metadata>
{defs(motif)}
{field_layer(motif)}
{motif_body(motif)}
</svg>
'''


def contact_sheet_svg() -> str:
    cols = 4
    thumb_w = 300
    thumb_h = 200
    cell_w = 365
    cell_h = 310
    pad_x = 54
    pad_y = 150
    rows = (len(MOTIFS) + cols - 1) // cols
    width = cols * cell_w + pad_x * 2
    height = rows * cell_h + pad_y + 70
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="{c("paper")}"/>',
        f'<text x="{pad_x}" y="62" font-family="Arial, Helvetica, sans-serif" font-size="32" font-weight="700" fill="{c("ink")}">SpaceBio-Bench biological motif asset pack v0.2</text>',
        f'<text x="{pad_x}" y="98" font-family="Arial, Helvetica, sans-serif" font-size="17" fill="{c("muted")}">Transparent SVG motif plates with scientific texture, depth, and biological specificity.</text>',
    ]
    for idx, motif in enumerate(MOTIFS):
        col = idx % cols
        row = idx // cols
        x = pad_x + col * cell_w
        y = pad_y + row * cell_h
        parts += [
            f'<rect x="{x - 10}" y="{y - 18}" width="{thumb_w + 28}" height="{thumb_h + 110}" fill="#FFFFFF" fill-opacity="0.46" stroke="{c("rule")}" stroke-width="1" stroke-opacity="0.30"/>',
            f'<image href="svg/{motif.motif_id}.svg" x="{x + 4}" y="{y}" width="{thumb_w}" height="{thumb_h}" preserveAspectRatio="xMidYMid meet"/>',
            f'<text x="{x}" y="{y + 236}" font-family="Arial, Helvetica, sans-serif" font-size="18" font-weight="700" fill="{c("ink")}">{motif.label}</text>',
            f'<text x="{x}" y="{y + 262}" font-family="Arial, Helvetica, sans-serif" font-size="14" fill="{c("muted")}">{motif.category} | {motif.motif_id}</text>',
        ]
    parts.append("</svg>\n")
    return "\n".join(parts)


def usage_sheet_svg() -> str:
    width = 3840
    height = 2160
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="{c("paper")}"/>',
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="{c("rule")}" opacity="0.035"/>',
        f'<text x="180" y="190" font-family="Arial, Helvetica, sans-serif" font-size="82" font-weight="700" fill="{c("ink")}">Biological motif usage test</text>',
        f'<text x="180" y="265" font-family="Arial, Helvetica, sans-serif" font-size="34" fill="{c("muted")}">Full-canvas Z-stack check: biological specificity without card-box decoration.</text>',
        f'<path d="M180 396 H3660 M180 1090 H3660 M1920 396 V1900" stroke="{c("rule")}" stroke-width="2" opacity="0.28"/>',
        f'<path d="M180 1900 H3660" stroke="{c("rule")}" stroke-width="2" opacity="0.22"/>',
    ]
    panels = [
        (
            "source_to_task",
            "Source context becomes a task",
            "sample, RNA signal, expression matrix",
            230,
            465,
            ["sample_to_matrix_assay", "rna_expression_trace", "single_cell_droplet_field"],
        ),
        (
            "feature_abstraction",
            "Molecular features become biology",
            "RNA trace, pathway program, protein context",
            2060,
            465,
            ["rna_expression_trace", "pathway_program_field", "protein_complex_field"],
        ),
        (
            "organoid_neural",
            "Organoid and neural extension",
            "3D culture marker plus neural context",
            230,
            1160,
            ["organoid_spheroid_field", "brain_neural_context", "single_cell_droplet_field"],
        ),
        (
            "tissue_context",
            "Tissue and organ context",
            "tissue section, kidney, liver, spatial spots",
            2060,
            1160,
            ["tissue_contour_section", "kidney_tubule_context", "liver_lobule_context", "spatial_spot_tissue_field"],
        ),
    ]
    for _, title, subtitle, x, y, ids in panels:
        parts.append(f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="42" font-weight="700" fill="{c("ink")}">{title}</text>')
        parts.append(f'<text x="{x}" y="{y + 50}" font-family="Arial, Helvetica, sans-serif" font-size="24" fill="{c("muted")}">{subtitle}</text>')
        parts.append(f'<path d="M{x} {y + 96} H{x + 1420}" stroke="{c("rule")}" stroke-width="2" opacity="0.22"/>')
        if len(ids) == 3:
            placements = [(x + 30, y + 135, 420, 280), (x + 485, y + 145, 420, 280), (x + 940, y + 128, 420, 280)]
        else:
            placements = [(x + 0, y + 140, 350, 234), (x + 350, y + 140, 350, 234), (x + 700, y + 140, 350, 234), (x + 1050, y + 140, 350, 234)]
        for motif_id, (ix, iy, iw, ih) in zip(ids, placements):
            href = f"svg/{motif_id}.svg"
            parts.append(f'<image href="{href}" x="{ix}" y="{iy}" width="{iw}" height="{ih}" preserveAspectRatio="xMidYMid meet"/>')
        parts.append(f'<path d="M{x + 120} {y + 500} C{x + 520} {y + 560}, {x + 960} {y + 560}, {x + 1320} {y + 500}" stroke="{c("rule")}" stroke-width="3" opacity="0.25" fill="none"/>')
    parts.append(f'<text x="180" y="2030" font-family="Arial, Helvetica, sans-serif" font-size="24" fill="{c("muted")}">QA verdict target: motifs should behave as restrained biological context, not as decorative clip-art or figure replacements.</text>')
    parts.append("</svg>\n")
    return "\n".join(parts)


def run_rsvg(svg_path: Path, png_path: Path, width: int | None = None, height: int | None = None) -> bool:
    exe = shutil.which("rsvg-convert")
    if not exe:
        return False
    cmd = [exe]
    if width is not None:
        cmd += ["-w", str(width)]
    if height is not None:
        cmd += ["-h", str(height)]
    cmd += ["-o", str(png_path), str(svg_path)]
    subprocess.run(cmd, check=True)
    return True


def write_manifest(previews: dict[str, str], contact_svg: Path, contact_png: Path | None, usage_svg: Path, usage_png: Path | None) -> None:
    records = []
    for motif in MOTIFS:
        svg = SVG_DIR / f"{motif.motif_id}.svg"
        records.append(
            {
                "motif_id": motif.motif_id,
                "category": motif.category,
                "label": motif.label,
                "semantic_use": motif.semantic_use,
                "primary_color_token": motif.primary_color,
                "secondary_color_token": motif.secondary_color,
                "svg": str(svg.relative_to(ROOT)),
                "png_preview": previews.get(motif.motif_id, ""),
                "license": "Original project artwork; repository MIT license.",
                "caution": motif.caution,
            }
        )
    manifest = {
        "created": "2026-06-02",
        "version": "0.2",
        "scope": "Premium biological motif plates for SpaceBio-Bench slide/deck visuals.",
        "design_rule": "Use as Z2/Z3 biological texture or signpost; do not substitute for measured source evidence.",
        "contact_sheet_svg": str(contact_svg.relative_to(ROOT)),
        "contact_sheet_png": str(contact_png.relative_to(ROOT)) if contact_png else "",
        "usage_sheet_svg": str(usage_svg.relative_to(ROOT)),
        "usage_sheet_png": str(usage_png.relative_to(ROOT)) if usage_png else "",
        "motif_count": len(records),
        "records": records,
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with (OUT / "manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    ensure_dirs()
    previews: dict[str, str] = {}
    for motif in MOTIFS:
        svg = SVG_DIR / f"{motif.motif_id}.svg"
        svg.write_text(svg_for(motif))
        png = PNG_DIR / f"{motif.motif_id}.png"
        if run_rsvg(svg, png, width=1440, height=960):
            previews[motif.motif_id] = str(png.relative_to(ROOT))

    contact_svg = OUT / "biovis_motif_pack_contact_sheet.svg"
    contact_svg.write_text(contact_sheet_svg())
    contact_png = OUT / "biovis_motif_pack_contact_sheet.png"
    contact_png_path = contact_png if run_rsvg(contact_svg, contact_png) else None
    usage_svg = OUT / "biovis_motif_pack_usage_sheet.svg"
    usage_svg.write_text(usage_sheet_svg())
    usage_png = OUT / "biovis_motif_pack_usage_sheet.png"
    usage_png_path = usage_png if run_rsvg(usage_svg, usage_png) else None
    write_manifest(previews, contact_svg, contact_png_path, usage_svg, usage_png_path)

    print(
        json.dumps(
            {
                "output": str(OUT.relative_to(ROOT)),
                "motif_count": len(MOTIFS),
                "contact_sheet_png": str(contact_png.relative_to(ROOT)) if contact_png.exists() else "",
                "usage_sheet_png": str(usage_png.relative_to(ROOT)) if usage_png.exists() else "",
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

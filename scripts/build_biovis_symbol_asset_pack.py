#!/usr/bin/env python3
"""Build an original SpaceBio-Bench biological visual symbol pack.

The pack is SVG-first. PNG previews are generated with rsvg-convert when it is
available, avoiding heavyweight plotting dependencies for simple icon assets.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "assets" / "biovis_symbol_pack_v0_1"
SVG_DIR = OUT / "svg"
PNG_DIR = OUT / "png"

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
}


@dataclass(frozen=True)
class Asset:
    asset_id: str
    category: str
    label: str
    semantic_use: str
    primary_color: str
    secondary_color: str
    caution: str


ASSETS = [
    Asset("cell_nucleus", "cellular", "Cell + nucleus", "single-cell, sample, cellular state", "bio_green", "source_blue", "Do not imply a specific cell type without a label."),
    Asset("mitochondrion", "cellular", "Mitochondrion", "metabolism, stress response, organelle biology", "assay_teal", "label_amber", "Use as organelle shorthand, not measured microscopy evidence."),
    Asset("organoid_sphere", "cellular", "Organoid sphere", "human organoid or 3D culture extension", "model_purple", "bio_green", "Pair with source/caveat text for draft organoid tracks."),
    Asset("tissue_section", "tissue", "Tissue section", "bulk tissue or section-level evidence", "bio_green", "rule", "Avoid using as histology proof unless tied to a real image."),
    Asset("epithelium_layer", "tissue", "Epithelial layer", "barrier tissue, lining, tissue architecture", "assay_teal", "bio_green", "Generic tissue architecture symbol only."),
    Asset("vascular_slice", "tissue", "Vascular slice", "vascular, perfusion, cfRNA origin, transport", "test_red", "source_blue", "Use red only when vascular/test-boundary semantics are not confused."),
    Asset("dna_helix", "molecular", "DNA helix", "gene features, sequence-level context", "source_blue", "model_purple", "Do not use for RNA-seq expression unless caption clarifies."),
    Asset("rna_strand", "molecular", "RNA strand", "RNA-seq, transcript expression, molecular readout", "assay_teal", "source_blue", "Good default for expression assays."),
    Asset("protein_complex", "molecular", "Protein complex", "protein/pathway or proteomics shorthand", "label_amber", "bio_green", "Do not imply direct proteomics unless data support it."),
    Asset("pathway_nodes", "molecular", "Pathway nodes", "pathway summaries, network biology", "bio_green", "label_amber", "Use for summarized programs, not causal pathways by default."),
    Asset("brain_outline", "organ", "Brain", "neural, brain/organoid interpretation", "model_purple", "source_blue", "Keep distinct from brain organoid unless labeled."),
    Asset("lung_outline", "organ", "Lung", "respiration, thoracic tissue, spaceflight physiology", "assay_teal", "source_blue", "Generic organ symbol only."),
    Asset("kidney_outline", "organ", "Kidney", "kidney tissue or renal transfer result", "bio_green", "label_amber", "Use with tissue labels for specificity."),
    Asset("liver_lobule", "organ", "Liver lobule", "liver/metabolism tissue context", "label_amber", "test_red", "Avoid overusing amber/orange as dominant palette."),
    Asset("sample_tube", "assay", "Sample tube", "sample collection or source material", "source_blue", "assay_teal", "Use as sample symbol, not as a wet-lab protocol claim."),
    Asset("expression_matrix", "assay", "Expression matrix", "matrix, task input, measured features", "muted", "model_purple", "Pair with task/source context so it is not a loose table."),
]


def c(name: str) -> str:
    return COLORS[name]


def icon_body(asset: Asset) -> list[str]:
    p = c(asset.primary_color)
    s = c(asset.secondary_color)
    ink = c("ink")
    muted = c("muted")
    rule = c("rule")
    paper = c("paper")
    common = 'stroke-linecap="round" stroke-linejoin="round"'
    aid = asset.asset_id

    if aid == "cell_nucleus":
        return [
            f'<ellipse cx="256" cy="256" rx="165" ry="120" fill="{paper}" stroke="{p}" stroke-width="16" opacity="0.96" {common}/>',
            f'<circle cx="246" cy="256" r="48" fill="{s}" opacity="0.82"/>',
            f'<circle cx="300" cy="220" r="13" fill="{p}" opacity="0.42"/><circle cx="180" cy="286" r="10" fill="{p}" opacity="0.35"/>',
            f'<path d="M143 246 C190 218, 322 214, 373 252" fill="none" stroke="{rule}" stroke-width="7" opacity="0.45" {common}/>',
        ]
    if aid == "mitochondrion":
        return [
            f'<path d="M151 296 C116 225, 170 144, 265 150 C363 157, 414 222, 383 296 C351 370, 231 389, 151 296 Z" fill="{paper}" stroke="{p}" stroke-width="15" {common}/>',
            f'<path d="M177 282 C216 229, 242 332, 277 277 C310 226, 329 319, 363 263" fill="none" stroke="{s}" stroke-width="15" opacity="0.82" {common}/>',
            f'<path d="M175 221 C237 183, 326 190, 365 239" fill="none" stroke="{rule}" stroke-width="7" opacity="0.42" {common}/>',
        ]
    if aid == "organoid_sphere":
        return [
            f'<circle cx="256" cy="256" r="142" fill="{paper}" stroke="{p}" stroke-width="13" opacity="0.97"/>',
            f'<circle cx="215" cy="217" r="45" fill="{s}" opacity="0.45"/><circle cx="290" cy="214" r="40" fill="{p}" opacity="0.33"/>',
            f'<circle cx="204" cy="294" r="39" fill="{p}" opacity="0.38"/><circle cx="290" cy="298" r="50" fill="{s}" opacity="0.35"/>',
            f'<circle cx="254" cy="258" r="28" fill="{ink}" opacity="0.18"/>',
        ]
    if aid == "tissue_section":
        return [
            f'<path d="M112 270 C128 154, 246 105, 372 142 C443 164, 426 286, 355 352 C263 438, 95 386, 112 270 Z" fill="{paper}" stroke="{p}" stroke-width="13" {common}/>',
            f'<path d="M143 281 C210 230, 306 231, 388 284" fill="none" stroke="{rule}" stroke-width="8" opacity="0.42" {common}/>',
            f'<path d="M163 335 C231 300, 308 306, 360 346" fill="none" stroke="{s}" stroke-width="8" opacity="0.55" {common}/>',
            f'<circle cx="206" cy="235" r="15" fill="{p}" opacity="0.45"/><circle cx="278" cy="260" r="14" fill="{s}" opacity="0.58"/><circle cx="333" cy="219" r="12" fill="{p}" opacity="0.35"/>',
        ]
    if aid == "epithelium_layer":
        return [
            f'<path d="M112 334 H400" stroke="{ink}" stroke-width="10" opacity="0.52" {common}/>',
            f'<path d="M126 316 C152 240, 192 240, 218 316 Z" fill="{paper}" stroke="{p}" stroke-width="11" {common}/>',
            f'<path d="M207 316 C233 228, 281 228, 307 316 Z" fill="{paper}" stroke="{s}" stroke-width="11" {common}/>',
            f'<path d="M296 316 C322 240, 362 240, 388 316 Z" fill="{paper}" stroke="{p}" stroke-width="11" {common}/>',
            f'<circle cx="172" cy="286" r="13" fill="{p}" opacity="0.45"/><circle cx="257" cy="280" r="13" fill="{s}" opacity="0.55"/><circle cx="342" cy="286" r="13" fill="{p}" opacity="0.45"/>',
        ]
    if aid == "vascular_slice":
        return [
            f'<circle cx="256" cy="256" r="144" fill="{paper}" stroke="{p}" stroke-width="14"/>',
            f'<circle cx="256" cy="256" r="78" fill="none" stroke="{s}" stroke-width="18" opacity="0.75"/>',
            f'<path d="M185 211 C230 168, 297 169, 337 212" fill="none" stroke="{rule}" stroke-width="7" opacity="0.45" {common}/>',
            f'<circle cx="256" cy="256" r="31" fill="{paper}" stroke="{muted}" stroke-width="7" opacity="0.80"/>',
        ]
    if aid == "dna_helix":
        body = [
            f'<path d="M175 102 C351 176, 351 336, 175 410" fill="none" stroke="{p}" stroke-width="15" {common}/>',
            f'<path d="M337 102 C161 176, 161 336, 337 410" fill="none" stroke="{s}" stroke-width="15" {common}/>',
        ]
        for y, x1, x2 in [(148, 201, 311), (196, 175, 337), (244, 186, 326), (292, 186, 326), (340, 175, 337), (388, 201, 311)]:
            body.append(f'<path d="M{x1} {y} H{x2}" stroke="{rule}" stroke-width="9" opacity="0.58" {common}/>')
        return body
    if aid == "rna_strand":
        body = [f'<path d="M123 298 C173 188, 249 356, 310 227 C341 163, 375 153, 400 203" fill="none" stroke="{p}" stroke-width="17" {common}/>']
        for x, y, dx, dy in [(158, 258, 42, -32), (210, 282, 42, 25), (268, 278, 42, -32), (326, 211, 42, 20), (374, 193, 42, -18)]:
            body.append(f'<path d="M{x} {y} l{dx} {dy}" stroke="{s}" stroke-width="10" opacity="0.76" {common}/><circle cx="{x + dx}" cy="{y + dy}" r="8" fill="{s}" opacity="0.70"/>')
        return body
    if aid == "protein_complex":
        return [
            f'<circle cx="207" cy="235" r="57" fill="{p}" opacity="0.72"/><circle cx="286" cy="223" r="48" fill="{s}" opacity="0.64"/>',
            f'<circle cx="255" cy="303" r="61" fill="{p}" opacity="0.45"/><circle cx="333" cy="297" r="38" fill="{s}" opacity="0.52"/>',
            f'<path d="M178 183 C230 137, 330 151, 376 220" fill="none" stroke="{ink}" stroke-width="8" opacity="0.24" {common}/>',
        ]
    if aid == "pathway_nodes":
        return [
            f'<circle cx="150" cy="255" r="34" fill="{paper}" stroke="{p}" stroke-width="12"/><circle cx="256" cy="174" r="34" fill="{paper}" stroke="{s}" stroke-width="12"/>',
            f'<circle cx="256" cy="336" r="34" fill="{paper}" stroke="{s}" stroke-width="12"/><circle cx="364" cy="255" r="34" fill="{paper}" stroke="{p}" stroke-width="12"/>',
            f'<path d="M184 242 L222 197 M184 268 L222 314 M290 197 L330 242 M290 314 L330 268" stroke="{rule}" stroke-width="9" opacity="0.70" {common}/>',
            f'<path d="M194 255 H320" stroke="{ink}" stroke-width="7" opacity="0.32" {common}/>',
        ]
    if aid == "brain_outline":
        return [
            f'<path d="M154 285 C100 242, 128 153, 204 154 C233 98, 329 105, 348 168 C418 181, 426 273, 369 312 C353 379, 251 395, 216 339 C187 340, 165 319, 154 285 Z" fill="{paper}" stroke="{p}" stroke-width="13" {common}/>',
            f'<path d="M211 174 C188 220, 213 253, 260 249 M285 143 C257 193, 290 229, 342 224 M194 307 C249 290, 302 300, 348 333" fill="none" stroke="{s}" stroke-width="8" opacity="0.62" {common}/>',
        ]
    if aid == "lung_outline":
        return [
            f'<path d="M256 122 V236" stroke="{s}" stroke-width="15" {common}/>',
            f'<path d="M256 218 C210 188, 170 184, 143 235 C112 294, 133 369, 205 377 C238 381, 250 336, 250 263 C250 240, 242 226, 222 214" fill="{paper}" stroke="{p}" stroke-width="13" {common}/>',
            f'<path d="M256 218 C302 188, 342 184, 369 235 C400 294, 379 369, 307 377 C274 381, 262 336, 262 263 C262 240, 270 226, 290 214" fill="{paper}" stroke="{p}" stroke-width="13" {common}/>',
            f'<path d="M214 263 C190 285, 185 318, 196 346 M298 263 C322 285, 327 318, 316 346" fill="none" stroke="{rule}" stroke-width="7" opacity="0.55" {common}/>',
        ]
    if aid == "kidney_outline":
        return [
            f'<path d="M289 117 C212 95, 148 162, 145 257 C142 363, 226 421, 314 379 C369 353, 385 290, 347 254 C316 225, 311 204, 343 177 C333 146, 315 125, 289 117 Z" fill="{paper}" stroke="{p}" stroke-width="14" {common}/>',
            f'<path d="M301 184 C259 223, 259 290, 312 329" fill="none" stroke="{s}" stroke-width="13" opacity="0.72" {common}/><path d="M315 255 H390" stroke="{s}" stroke-width="10" opacity="0.72" {common}/>',
        ]
    if aid == "liver_lobule":
        return [
            f'<path d="M256 119 L376 188 V326 L256 393 L136 326 V188 Z" fill="{paper}" stroke="{p}" stroke-width="13" {common}/>',
            f'<circle cx="256" cy="256" r="34" fill="{s}" opacity="0.70"/>',
            f'<path d="M256 154 V220 M256 292 V358 M167 205 L226 238 M286 274 L345 307 M345 205 L286 238 M226 274 L167 307" stroke="{rule}" stroke-width="8" opacity="0.62" {common}/>',
        ]
    if aid == "sample_tube":
        return [
            f'<path d="M191 124 H321 M212 124 V324 C212 376, 300 376, 300 324 V124" fill="{paper}" stroke="{p}" stroke-width="13" {common}/>',
            f'<path d="M214 284 C243 268, 270 306, 298 286 V331 C298 365, 214 365, 214 331 Z" fill="{s}" opacity="0.62"/>',
            f'<path d="M229 177 H283 M229 220 H283" stroke="{rule}" stroke-width="7" opacity="0.60" {common}/>',
        ]
    if aid == "expression_matrix":
        body = [f'<rect x="126" y="140" width="260" height="232" fill="{paper}" stroke="{p}" stroke-width="12"/>']
        for i in range(1, 6):
            x = 126 + i * 43
            body.append(f'<path d="M{x} 140 V372" stroke="{rule}" stroke-width="5" opacity="0.46"/>')
        for i in range(1, 5):
            y = 140 + i * 46
            body.append(f'<path d="M126 {y} H386" stroke="{rule}" stroke-width="5" opacity="0.46"/>')
        for x, y in [(169, 186), (255, 232), (341, 278), (212, 324), (298, 186)]:
            body.append(f'<rect x="{x - 14}" y="{y - 14}" width="28" height="28" rx="5" fill="{s}" opacity="0.70"/>')
        return body
    raise ValueError(aid)


def svg_for(asset: Asset) -> str:
    body = "\n".join(f"  {line}" for line in icon_body(asset))
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="512" height="512" viewBox="0 0 512 512" role="img" aria-label="{asset.label}">
  <metadata>SpaceBio-Bench original biological symbol asset v0.1; asset_id={asset.asset_id}; category={asset.category}; license=repository MIT license.</metadata>
  <rect width="512" height="512" fill="none"/>
{body}
</svg>
'''


def contact_sheet_svg() -> str:
    cols = 4
    cell_w = 350
    cell_h = 330
    pad_x = 52
    pad_y = 150
    width = cols * cell_w + pad_x * 2
    rows = (len(ASSETS) + cols - 1) // cols
    height = rows * cell_h + pad_y + 70
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="{c("paper")}"/>',
        f'<text x="{pad_x}" y="62" font-family="Arial, Helvetica, sans-serif" font-size="32" font-weight="700" fill="{c("ink")}">SpaceBio-Bench biological symbol asset pack v0.1</text>',
        f'<text x="{pad_x}" y="98" font-family="Arial, Helvetica, sans-serif" font-size="17" fill="{c("muted")}">Original editable SVG primitives for cell, molecular, tissue, organ, and assay visuals.</text>',
    ]
    for idx, asset in enumerate(ASSETS):
        col = idx % cols
        row = idx // cols
        x = pad_x + col * cell_w
        y = pad_y + row * cell_h
        href = f"svg/{asset.asset_id}.svg"
        parts.append(f'<rect x="{x - 10}" y="{y - 16}" width="300" height="280" rx="0" fill="none" stroke="{c("rule")}" stroke-width="1" opacity="0.22"/>')
        parts.append(f'<image href="{href}" x="{x + 38}" y="{y}" width="190" height="190" preserveAspectRatio="xMidYMid meet"/>')
        parts.append(f'<text x="{x}" y="{y + 225}" font-family="Arial, Helvetica, sans-serif" font-size="18" font-weight="700" fill="{c("ink")}">{asset.label}</text>')
        parts.append(f'<text x="{x}" y="{y + 250}" font-family="Arial, Helvetica, sans-serif" font-size="14" fill="{c("muted")}">{asset.category} | {asset.asset_id}</text>')
    parts.append("</svg>\n")
    return "\n".join(parts)


def run_rsvg(svg_path: Path, png_path: Path, size: int | None = None) -> bool:
    exe = shutil.which("rsvg-convert")
    if not exe:
        return False
    cmd = [exe]
    if size is not None:
        cmd += ["-w", str(size), "-h", str(size)]
    cmd += ["-o", str(png_path), str(svg_path)]
    subprocess.run(cmd, check=True)
    return True


def ensure_dirs() -> None:
    SVG_DIR.mkdir(parents=True, exist_ok=True)
    PNG_DIR.mkdir(parents=True, exist_ok=True)


def write_manifest(previews: dict[str, str], contact_sheet_svg_path: Path, contact_sheet_png_path: Path | None) -> None:
    records = []
    for asset in ASSETS:
        svg_path = SVG_DIR / f"{asset.asset_id}.svg"
        records.append(
            {
                "asset_id": asset.asset_id,
                "category": asset.category,
                "label": asset.label,
                "semantic_use": asset.semantic_use,
                "primary_color_token": asset.primary_color,
                "secondary_color_token": asset.secondary_color,
                "svg": str(svg_path.relative_to(ROOT)),
                "png_preview": previews.get(asset.asset_id, ""),
                "license": "Original project artwork; repository MIT license.",
                "caution": asset.caution,
            }
        )
    manifest = {
        "created": "2026-06-02",
        "version": "0.1",
        "scope": "Original biological symbols for SpaceBio-Bench premium slide/deck visuals.",
        "design_rule": "Use as subtle scientific object vocabulary; do not substitute these symbols for source-derived proof images.",
        "contact_sheet_svg": str(contact_sheet_svg_path.relative_to(ROOT)),
        "contact_sheet_png": str(contact_sheet_png_path.relative_to(ROOT)) if contact_sheet_png_path else "",
        "asset_count": len(records),
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
    for asset in ASSETS:
        svg_path = SVG_DIR / f"{asset.asset_id}.svg"
        svg_path.write_text(svg_for(asset))
        png_path = PNG_DIR / f"{asset.asset_id}.png"
        if run_rsvg(svg_path, png_path, size=1024):
            previews[asset.asset_id] = str(png_path.relative_to(ROOT))

    contact_svg = OUT / "biovis_symbol_pack_contact_sheet.svg"
    contact_svg.write_text(contact_sheet_svg())
    contact_png = OUT / "biovis_symbol_pack_contact_sheet.png"
    contact_png_path = contact_png if run_rsvg(contact_svg, contact_png) else None

    write_manifest(previews, contact_svg, contact_png_path)
    print(json.dumps({"output": str(OUT.relative_to(ROOT)), "asset_count": len(ASSETS), "contact_sheet_png": str(contact_png.relative_to(ROOT)) if contact_png.exists() else ""}, indent=2))


if __name__ == "__main__":
    main()

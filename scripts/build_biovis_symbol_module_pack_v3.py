#!/usr/bin/env python3
"""Build reusable biological symbols, badges, and modules for v0.3 visuals.

This pack complements the v0.3 biological motif plates. The motif pack provides
texture and atmosphere; this pack provides small, editable visual vocabulary
for methods, evidence status, species/modality scope, and claim boundaries.
"""

from __future__ import annotations

import csv
import html
import json
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "assets" / "biovis_symbol_module_pack_v0_3"
SYMBOL_SVG = OUT / "symbols" / "svg"
SYMBOL_PNG = OUT / "symbols" / "png"
BADGE_SVG = OUT / "badges" / "svg"
BADGE_PNG = OUT / "badges" / "png"
MODULE_SVG = OUT / "modules" / "svg"
MODULE_PNG = OUT / "modules" / "png"
MODULE_DARK_SVG = OUT / "modules_dark" / "svg"
MODULE_DARK_PNG = OUT / "modules_dark" / "png"
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
}


@dataclass(frozen=True)
class SymbolAsset:
    asset_id: str
    category: str
    label: str
    semantic_use: str
    primary: str
    secondary: str
    allowed_use: str
    prohibited_use: str


@dataclass(frozen=True)
class BadgeAsset:
    asset_id: str
    label: str
    tone: str
    semantic_use: str
    allowed_use: str
    prohibited_use: str


@dataclass(frozen=True)
class ModuleAsset:
    asset_id: str
    label: str
    module_type: str
    semantic_use: str
    allowed_use: str
    prohibited_use: str


SYMBOLS: list[SymbolAsset] = [
    SymbolAsset("gene_locus", "molecular", "Gene locus", "gene-level features, sequence context", "blue", "purple", "gene/feature signpost", "do not imply DNA-seq evidence by itself"),
    SymbolAsset("rna_readout", "molecular", "RNA readout", "RNA-seq and expression measurement", "teal", "blue", "default expression assay marker", "do not use for protein abundance without label"),
    SymbolAsset("protein_program", "molecular", "Protein program", "proteins, complexes, pathway activity", "amber", "green", "protein/pathway shorthand", "do not imply direct proteomics unless supported"),
    SymbolAsset("pathway_network", "molecular", "Pathway network", "pathway-level interpretation", "green", "amber", "program or network summary", "do not present as causal wiring"),
    SymbolAsset("mitochondrial_stress", "cellular", "Mitochondrial stress", "organelle stress and metabolism", "teal", "rose", "stress biology marker", "do not overclaim microscopy evidence"),
    SymbolAsset("cell_state", "cellular", "Cell state", "cell-level state or phenotype", "green", "blue", "single-cell or cell-state marker", "do not imply a specific cell type without label"),
    SymbolAsset("cell_population", "cellular", "Cell population", "cell clusters and population shifts", "blue", "green", "population or embedding marker", "do not use as measured UMAP without data"),
    SymbolAsset("organoid_rosette", "cellular", "Organoid rosette", "organoid or 3D culture extension", "purple", "green", "organoid/culture context", "do not imply organoid validation if it is only future scope"),
    SymbolAsset("tissue_section", "tissue", "Tissue section", "bulk tissue or histology-like context", "amber", "rose", "tissue evidence anchor", "do not present as source histology"),
    SymbolAsset("spatial_spots", "tissue", "Spatial spots", "spatial transcriptomics or spatial assay", "purple", "green", "spatial modality marker", "do not imply spatial data exists without source"),
    SymbolAsset("vascular_context", "tissue", "Vascular context", "perfusion, transport, vascular tissue", "red", "blue", "vascular biology cue", "do not let red imply failure unless intended"),
    SymbolAsset("epithelial_barrier", "tissue", "Epithelial barrier", "barrier or lining tissue", "teal", "green", "barrier tissue cue", "do not imply a specific organ by itself"),
    SymbolAsset("brain_organ", "organ", "Brain organ", "neural tissue or brain organoid context", "purple", "blue", "neural result anchor", "keep distinct from organoid unless labeled"),
    SymbolAsset("lung_organ", "organ", "Lung organ", "lung or respiration biology", "teal", "blue", "lung/tissue anchor", "generic organ symbol only"),
    SymbolAsset("kidney_organ", "organ", "Kidney organ", "renal tissue and kidney result slides", "green", "amber", "kidney/tissue anchor", "generic organ symbol only"),
    SymbolAsset("liver_organ", "organ", "Liver organ", "liver, metabolism, detox biology", "amber", "green", "liver/tissue anchor", "avoid amber-dominant slide palettes"),
    SymbolAsset("sample_tube", "assay", "Sample tube", "sample source and collection", "blue", "teal", "sample/source marker", "not a wet-lab protocol proof"),
    SymbolAsset("bulk_rna_assay", "assay", "Bulk RNA assay", "bulk RNA-seq or transcriptomics", "teal", "blue", "bulk expression lane", "do not use for single-cell payloads"),
    SymbolAsset("single_cell_droplet", "assay", "Single-cell droplet", "single-cell capture and cell barcodes", "blue", "green", "scRNA-seq modality marker", "do not imply actual droplet chemistry"),
    SymbolAsset("proteomics_assay", "assay", "Proteomics assay", "mass-spec or proteomics readout", "amber", "purple", "proteomics modality marker", "do not imply proteomics data without source"),
    SymbolAsset("microscopy_frame", "assay", "Microscopy frame", "imaging, microscopy, histology proof slot", "rose", "blue", "source-image placeholder", "replace with source image for empirical claim"),
    SymbolAsset("omics_matrix", "assay", "Omics matrix", "feature matrix or expression table", "blue", "amber", "matrix/task input marker", "do not use as real heatmap values"),
    SymbolAsset("human_model", "species", "Human model", "human samples, cell lines, organoids", "blue", "purple", "human/model-system marker", "do not imply clinical cohort"),
    SymbolAsset("mouse_model", "species", "Mouse model", "mouse and rodent experiments", "green", "amber", "mouse model marker", "do not use for rat"),
    SymbolAsset("rat_model", "species", "Rat model", "rat model or rodent extension", "amber", "green", "rat model marker", "do not use interchangeably with mouse"),
    SymbolAsset("fly_model", "species", "Fly model", "Drosophila experiments", "purple", "blue", "fly/model organism marker", "do not imply mammalian result"),
    SymbolAsset("worm_model", "species", "Worm model", "C. elegans or worm experiments", "teal", "green", "worm/model organism marker", "do not imply tissue anatomy"),
    SymbolAsset("fish_model", "species", "Fish model", "fish or zebrafish experiments", "blue", "teal", "fish/model organism marker", "do not imply mammal"),
    SymbolAsset("plant_model", "species", "Plant model", "plant space biology experiments", "green", "amber", "plant/model marker", "do not imply animal biology"),
    SymbolAsset("microbe_model", "species", "Microbe model", "microbes, yeast, microbial systems", "rose", "purple", "microbe/model marker", "do not imply pathogen relevance"),
    SymbolAsset("orbit_microgravity", "space-context", "Orbit microgravity", "spaceflight or microgravity context", "blue", "sky", "mission/stressor context", "do not use as measured gravity condition"),
    SymbolAsset("radiation_stressor", "space-context", "Radiation stressor", "radiation or stressor context", "red", "amber", "space stressor marker", "avoid alarm styling unless warranted"),
    SymbolAsset("source_record", "trust", "Source record", "source dataset or record provenance", "blue", "teal", "source provenance marker", "do not imply validation"),
    SymbolAsset("checksum_record", "trust", "Checksum record", "hash, freeze, reproducibility", "green", "blue", "artifact freeze marker", "not a scientific quality claim by itself"),
    SymbolAsset("split_guard", "trust", "Split guard", "train/test split or leakage boundary", "red", "blue", "leakage-control marker", "do not use as generic warning"),
    SymbolAsset("caveat_flag", "trust", "Caveat flag", "limitation, preliminary status, caveat", "amber", "red", "caveat marker", "do not hide major caveats in small icon-only form"),
]

BADGES: list[BadgeAsset] = [
    BadgeAsset("generated_context", "schematic", "muted", "generated context asset", "mark schematic/non-source imagery", "do not use on measured source figures"),
    BadgeAsset("source_proof", "source", "blue", "source-derived evidence", "mark real source-derived proof panels", "do not use for generated art"),
    BadgeAsset("processed_result", "processed", "teal", "processed analysis output", "mark computed outputs and derived results", "do not imply independent validation"),
    BadgeAsset("validated_result", "validated", "green", "validated result", "mark reviewed or validated result surfaces", "do not use for preliminary checks"),
    BadgeAsset("preliminary", "prelim", "amber", "draft/preliminary status", "mark draft or exploratory analysis", "do not hide uncertainty elsewhere"),
    BadgeAsset("caveat", "caveat", "red", "explicit limitation", "mark claim boundary or limitation", "do not use decoratively"),
    BadgeAsset("train_only", "train-only", "purple", "train-only processing", "mark leakage-safe train-only step", "do not use if test information leaks in"),
    BadgeAsset("held_out", "held-out", "blue", "held-out evaluation", "mark withheld test or mission split", "do not use for random internal checks"),
]

MODULES: list[ModuleAsset] = [
    ModuleAsset("biological_scale_ladder", "Biological scale ladder", "scale", "molecule-cell-tissue-organ-organism hierarchy", "methods overview, bridge slides, figure annotations", "do not present as a causal sequence"),
    ModuleAsset("sample_to_feature_stack", "Sample to feature stack", "workflow", "sample to assay to matrix to program to model", "methods pipeline explanation", "do not imply all steps were run for every dataset"),
    ModuleAsset("species_coverage_strip", "Species coverage strip", "scope", "multi-species/model-system scope", "scope slide, dataset coverage, extension planning", "do not imply included analysis without status labels"),
    ModuleAsset("modality_lane_set", "Modality lane set", "scope", "bulk, single-cell, spatial, proteomics, imaging, organoid modalities", "data collection and analysis scope", "do not imply availability for every species"),
    ModuleAsset("trust_status_ribbon", "Trust status ribbon", "trust", "source, checksum, split guard, processed result, caveat", "method QA and artifact governance", "do not use as generic decoration"),
    ModuleAsset("claim_boundary_bar", "Claim boundary bar", "trust", "schematic context vs source proof vs processed result vs validated claim", "claim boundary and footer modules", "do not shrink until labels become unreadable"),
    ModuleAsset("space_bio_assay_lane", "Space biology assay lane", "workflow", "mission/stressor to tissue to assay to feature layer", "SpaceBio-Bench method explanation", "do not imply a temporal biological mechanism"),
]


def col(token: str) -> str:
    return COLORS[token]


def esc(value: str) -> str:
    return html.escape(value, quote=True)


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def common_attrs() -> str:
    return 'stroke-linecap="round" stroke-linejoin="round" vector-effect="non-scaling-stroke"'


def icon_body(asset_id: str, primary: str, secondary: str) -> str:
    p = col(primary)
    s = col(secondary)
    ink = col("ink")
    paper = col("paper2")
    rule = col("rule")
    muted = col("muted")
    c = common_attrs()

    if asset_id == "gene_locus":
        rails = [
            f'<path d="M58 79 C128 49, 128 207, 198 177" fill="none" stroke="{p}" stroke-width="8" {c}/>',
            f'<path d="M198 79 C128 49, 128 207, 58 177" fill="none" stroke="{s}" stroke-width="8" {c}/>',
        ]
        for y, x1, x2 in [(89, 74, 182), (111, 64, 192), (133, 68, 188), (155, 64, 192), (177, 74, 182)]:
            rails.append(f'<path d="M{x1} {y} H{x2}" stroke="{rule}" stroke-width="5" opacity="0.65" {c}/>')
        return "\n".join(rails)
    if asset_id == "rna_readout":
        return "\n".join(
            [
                f'<path d="M42 144 C72 67, 115 188, 154 105 C177 55, 202 62, 216 96" fill="none" stroke="{p}" stroke-width="10" {c}/>',
                f'<path d="M70 129 l25 -18 M107 142 l25 17 M143 126 l25 -21 M183 88 l23 13" stroke="{s}" stroke-width="6" opacity="0.75" {c}/>',
                f'<circle cx="95" cy="111" r="5" fill="{s}"/><circle cx="132" cy="159" r="5" fill="{s}"/><circle cx="168" cy="105" r="5" fill="{s}"/><circle cx="206" cy="101" r="5" fill="{s}"/>',
            ]
        )
    if asset_id == "protein_program":
        return "\n".join(
            [
                f'<circle cx="92" cy="104" r="34" fill="{p}" opacity="0.75"/>',
                f'<circle cx="138" cy="99" r="29" fill="{s}" opacity="0.64"/>',
                f'<circle cx="120" cy="147" r="38" fill="{p}" opacity="0.44"/>',
                f'<circle cx="169" cy="144" r="24" fill="{s}" opacity="0.55"/>',
                f'<path d="M72 72 C111 40, 178 49, 205 93" fill="none" stroke="{ink}" stroke-width="5" opacity="0.24" {c}/>',
            ]
        )
    if asset_id == "pathway_network":
        return "\n".join(
            [
                f'<path d="M83 126 H172 M128 78 V178 M93 93 L164 163 M164 93 L93 163" stroke="{rule}" stroke-width="5" opacity="0.72" {c}/>',
                f'<circle cx="83" cy="126" r="18" fill="{paper}" stroke="{p}" stroke-width="7"/>',
                f'<circle cx="128" cy="78" r="18" fill="{paper}" stroke="{s}" stroke-width="7"/>',
                f'<circle cx="172" cy="126" r="18" fill="{paper}" stroke="{p}" stroke-width="7"/>',
                f'<circle cx="128" cy="178" r="18" fill="{paper}" stroke="{s}" stroke-width="7"/>',
            ]
        )
    if asset_id == "mitochondrial_stress":
        return "\n".join(
            [
                f'<path d="M68 151 C47 105, 79 57, 134 59 C191 62, 219 101, 199 151 C180 201, 111 207, 68 151 Z" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M83 141 C105 104, 119 167, 141 130 C160 99, 171 158, 191 122" fill="none" stroke="{s}" stroke-width="8" opacity="0.86" {c}/>',
                f'<path d="M196 52 l-18 32 h25 l-22 44" fill="none" stroke="{col("red")}" stroke-width="7" opacity="0.82" {c}/>',
            ]
        )
    if asset_id == "cell_state":
        return "\n".join(
            [
                f'<ellipse cx="128" cy="128" rx="85" ry="63" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<circle cx="123" cy="128" r="26" fill="{s}" opacity="0.78"/>',
                f'<circle cx="157" cy="108" r="7" fill="{p}" opacity="0.45"/><circle cx="93" cy="148" r="6" fill="{p}" opacity="0.38"/>',
                f'<path d="M69 125 C98 109, 160 109, 188 128" fill="none" stroke="{rule}" stroke-width="4" opacity="0.55" {c}/>',
            ]
        )
    if asset_id == "cell_population":
        body = []
        for x, y, r, token, opacity in [
            (79, 103, 14, p, 0.75),
            (108, 82, 11, s, 0.70),
            (132, 121, 16, p, 0.48),
            (166, 95, 13, s, 0.64),
            (181, 147, 15, p, 0.55),
            (116, 166, 14, s, 0.58),
            (73, 153, 11, p, 0.44),
        ]:
            body.append(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{token}" opacity="{opacity}"/>')
        body.append(f'<path d="M58 182 C91 203, 164 207, 202 171" fill="none" stroke="{rule}" stroke-width="5" opacity="0.48" {c}/>')
        return "\n".join(body)
    if asset_id == "organoid_rosette":
        return "\n".join(
            [
                f'<circle cx="128" cy="128" r="76" fill="{paper}" stroke="{p}" stroke-width="8"/>',
                f'<circle cx="101" cy="101" r="24" fill="{s}" opacity="0.44"/><circle cx="147" cy="99" r="22" fill="{p}" opacity="0.36"/>',
                f'<circle cx="98" cy="149" r="23" fill="{p}" opacity="0.39"/><circle cx="151" cy="150" r="28" fill="{s}" opacity="0.35"/>',
                f'<circle cx="127" cy="127" r="18" fill="{ink}" opacity="0.16"/>',
            ]
        )
    if asset_id == "tissue_section":
        return "\n".join(
            [
                f'<path d="M51 136 C59 72, 129 47, 196 70 C231 82, 222 153, 184 188 C132 236, 42 201, 51 136 Z" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M70 139 C104 114, 159 113, 204 143" fill="none" stroke="{rule}" stroke-width="5" opacity="0.50" {c}/>',
                f'<path d="M80 171 C117 151, 161 154, 190 177" fill="none" stroke="{s}" stroke-width="5" opacity="0.64" {c}/>',
                f'<circle cx="104" cy="117" r="8" fill="{p}" opacity="0.44"/><circle cx="145" cy="133" r="8" fill="{s}" opacity="0.60"/><circle cx="178" cy="111" r="7" fill="{p}" opacity="0.35"/>',
            ]
        )
    if asset_id == "spatial_spots":
        spots = []
        for x, y, token, op in [
            (78, 96, p, 0.55),
            (111, 82, s, 0.75),
            (146, 93, p, 0.55),
            (181, 115, s, 0.64),
            (92, 135, s, 0.66),
            (132, 134, p, 0.55),
            (171, 155, p, 0.45),
            (111, 172, s, 0.55),
            (149, 184, p, 0.42),
        ]:
            spots.append(f'<circle cx="{x}" cy="{y}" r="9" fill="{token}" opacity="{op}"/>')
        return "\n".join(
            [
                f'<path d="M51 136 C59 72, 129 47, 196 70 C231 82, 222 153, 184 188 C132 236, 42 201, 51 136 Z" fill="{paper}" stroke="{rule}" stroke-width="6" opacity="0.75" {c}/>',
                *spots,
            ]
        )
    if asset_id == "vascular_context":
        return "\n".join(
            [
                f'<circle cx="128" cy="128" r="75" fill="{paper}" stroke="{p}" stroke-width="8"/>',
                f'<circle cx="128" cy="128" r="40" fill="none" stroke="{s}" stroke-width="10" opacity="0.78"/>',
                f'<path d="M90 105 C113 83, 149 83, 169 105" fill="none" stroke="{rule}" stroke-width="4" opacity="0.52" {c}/>',
                f'<circle cx="128" cy="128" r="17" fill="{paper}" stroke="{muted}" stroke-width="4" opacity="0.82"/>',
            ]
        )
    if asset_id == "epithelial_barrier":
        return "\n".join(
            [
                f'<path d="M50 171 H207" stroke="{ink}" stroke-width="6" opacity="0.50" {c}/>',
                f'<path d="M59 160 C73 111, 98 111, 112 160 Z" fill="{paper}" stroke="{p}" stroke-width="7" {c}/>',
                f'<path d="M106 160 C121 103, 150 103, 165 160 Z" fill="{paper}" stroke="{s}" stroke-width="7" {c}/>',
                f'<path d="M157 160 C171 111, 196 111, 210 160 Z" fill="{paper}" stroke="{p}" stroke-width="7" {c}/>',
                f'<circle cx="86" cy="143" r="7" fill="{p}" opacity="0.45"/><circle cx="136" cy="139" r="7" fill="{s}" opacity="0.55"/><circle cx="184" cy="143" r="7" fill="{p}" opacity="0.45"/>',
            ]
        )
    if asset_id == "brain_organ":
        return "\n".join(
            [
                f'<path d="M74 145 C46 124, 61 78, 101 79 C117 50, 167 54, 177 87 C214 94, 218 143, 188 164 C179 198, 127 205, 109 176 C93 177, 80 164, 74 145 Z" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M104 90 C92 115, 106 132, 130 129 M145 72 C130 99, 148 117, 176 114 M96 158 C126 149, 154 154, 178 171" fill="none" stroke="{s}" stroke-width="5" opacity="0.62" {c}/>',
            ]
        )
    if asset_id == "lung_organ":
        return "\n".join(
            [
                f'<path d="M128 61 V121" stroke="{s}" stroke-width="8" {c}/>',
                f'<path d="M128 112 C103 95, 82 96, 67 123 C48 158, 61 199, 101 203 C121 205, 125 180, 125 143 C125 129, 120 119, 109 113" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M128 112 C153 95, 174 96, 189 123 C208 158, 195 199, 155 203 C135 205, 131 180, 131 143 C131 129, 136 119, 147 113" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M106 145 C94 157, 91 176, 97 191 M150 145 C162 157, 165 176, 159 191" fill="none" stroke="{rule}" stroke-width="4" opacity="0.55" {c}/>',
            ]
        )
    if asset_id == "kidney_organ":
        return "\n".join(
            [
                f'<path d="M145 59 C105 48, 72 84, 71 132 C70 188, 116 219, 162 196 C191 182, 199 149, 179 130 C163 114, 160 103, 177 89 C172 72, 160 63, 145 59 Z" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M151 94 C129 116, 130 151, 158 171" fill="none" stroke="{s}" stroke-width="7" opacity="0.76" {c}/>',
                f'<path d="M159 130 H205" stroke="{s}" stroke-width="6" opacity="0.76" {c}/>',
            ]
        )
    if asset_id == "liver_organ":
        return "\n".join(
            [
                f'<path d="M52 132 C60 85, 96 65, 147 75 C194 84, 223 111, 215 148 C207 187, 162 201, 111 191 C72 184, 46 165, 52 132 Z" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M106 107 C138 117, 158 140, 167 177" fill="none" stroke="{s}" stroke-width="7" opacity="0.75" {c}/>',
                f'<circle cx="109" cy="135" r="14" fill="{s}" opacity="0.42"/>',
            ]
        )
    if asset_id == "sample_tube":
        return "\n".join(
            [
                f'<path d="M93 57 H163 M105 57 V165 C105 194, 151 194, 151 165 V57" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M107 144 C122 136, 137 155, 149 145 V169 C149 188, 107 188, 107 169 Z" fill="{s}" opacity="0.66"/>',
                f'<path d="M113 86 H143 M113 109 H143" stroke="{rule}" stroke-width="4" opacity="0.60" {c}/>',
            ]
        )
    if asset_id == "bulk_rna_assay":
        return "\n".join(
            [
                f'<rect x="61" y="73" width="134" height="111" rx="12" fill="{paper}" stroke="{p}" stroke-width="8"/>',
                f'<path d="M84 135 C104 90, 130 160, 154 112 C169 81, 184 86, 194 106" fill="none" stroke="{s}" stroke-width="7" {c}/>',
                f'<path d="M82 161 H176" stroke="{rule}" stroke-width="5" opacity="0.48" {c}/>',
            ]
        )
    if asset_id == "single_cell_droplet":
        return "\n".join(
            [
                f'<path d="M128 52 C168 101, 187 130, 184 158 C181 193, 156 214, 128 214 C100 214, 75 193, 72 158 C69 130, 88 101, 128 52 Z" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<circle cx="128" cy="154" r="28" fill="{s}" opacity="0.58"/>',
                f'<circle cx="128" cy="154" r="11" fill="{p}" opacity="0.72"/>',
                f'<circle cx="103" cy="129" r="6" fill="{s}" opacity="0.58"/><circle cx="155" cy="129" r="6" fill="{s}" opacity="0.58"/>',
            ]
        )
    if asset_id == "proteomics_assay":
        return "\n".join(
            [
                f'<path d="M70 180 H188" stroke="{rule}" stroke-width="6" opacity="0.58" {c}/>',
                f'<path d="M91 180 V78 L128 62 L165 78 V180" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M99 141 H157 M106 118 H150 M114 95 H142" stroke="{s}" stroke-width="6" opacity="0.78" {c}/>',
                f'<circle cx="128" cy="62" r="14" fill="{s}" opacity="0.62"/>',
            ]
        )
    if asset_id == "microscopy_frame":
        return "\n".join(
            [
                f'<rect x="58" y="64" width="140" height="128" rx="12" fill="{paper}" stroke="{p}" stroke-width="8"/>',
                f'<circle cx="107" cy="124" r="22" fill="{s}" opacity="0.40"/><circle cx="150" cy="139" r="28" fill="{p}" opacity="0.31"/>',
                f'<path d="M78 165 C111 143, 157 145, 181 166" fill="none" stroke="{rule}" stroke-width="5" opacity="0.58" {c}/>',
                f'<path d="M77 85 H111 M77 102 H96" stroke="{s}" stroke-width="5" opacity="0.72" {c}/>',
            ]
        )
    if asset_id == "omics_matrix":
        body = [f'<rect x="59" y="65" width="138" height="126" fill="{paper}" stroke="{p}" stroke-width="8"/>']
        for i in range(1, 5):
            x = 59 + i * 28
            body.append(f'<path d="M{x} 65 V191" stroke="{rule}" stroke-width="3" opacity="0.46"/>')
        for i in range(1, 4):
            y = 65 + i * 31
            body.append(f'<path d="M59 {y} H197" stroke="{rule}" stroke-width="3" opacity="0.46"/>')
        for x, y, token in [(87, 96, s), (143, 127, s), (171, 158, s), (115, 158, s), (171, 96, s)]:
            body.append(f'<rect x="{x-8}" y="{y-8}" width="16" height="16" rx="3" fill="{token}" opacity="0.72"/>')
        return "\n".join(body)
    if asset_id == "human_model":
        return "\n".join(
            [
                f'<circle cx="128" cy="88" r="29" fill="{paper}" stroke="{p}" stroke-width="8"/>',
                f'<path d="M83 203 C87 154, 103 130, 128 130 C153 130, 169 154, 173 203 Z" fill="{paper}" stroke="{s}" stroke-width="8" {c}/>',
                f'<path d="M105 162 H151" stroke="{rule}" stroke-width="5" opacity="0.52" {c}/>',
            ]
        )
    if asset_id == "mouse_model":
        return "\n".join(
            [
                f'<path d="M55 146 C79 98, 150 91, 189 123 C203 135, 200 154, 184 166 C142 197, 83 184, 55 146 Z" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M74 150 C42 160, 38 184, 72 190" fill="none" stroke="{s}" stroke-width="6" opacity="0.82" {c}/>',
                f'<path d="M173 120 C185 103, 202 105, 207 123" fill="none" stroke="{s}" stroke-width="5" opacity="0.70" {c}/>',
                f'<path d="M104 171 C126 162, 150 160, 174 166" fill="none" stroke="{rule}" stroke-width="4" opacity="0.42" {c}/>',
            ]
        )
    if asset_id == "rat_model":
        return "\n".join(
            [
                f'<path d="M45 149 C74 103, 151 101, 199 128 C220 140, 217 160, 194 174 C145 203, 75 190, 45 149 Z" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M65 154 C18 170, 20 203, 76 198" fill="none" stroke="{s}" stroke-width="6" opacity="0.82" {c}/>',
                f'<path d="M183 127 C194 113, 209 116, 214 131" fill="none" stroke="{s}" stroke-width="5" opacity="0.64" {c}/>',
                f'<path d="M91 174 C122 163, 160 162, 190 171" fill="none" stroke="{rule}" stroke-width="4" opacity="0.42" {c}/>',
            ]
        )
    if asset_id == "fly_model":
        return "\n".join(
            [
                f'<ellipse cx="101" cy="113" rx="38" ry="22" fill="none" stroke="{s}" stroke-width="6" opacity="0.74" transform="rotate(-32 101 113)"/>',
                f'<ellipse cx="155" cy="113" rx="38" ry="22" fill="none" stroke="{s}" stroke-width="6" opacity="0.74" transform="rotate(32 155 113)"/>',
                f'<path d="M128 70 C146 101, 148 166, 128 203 C108 166, 110 101, 128 70 Z" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M128 90 V184" stroke="{rule}" stroke-width="4" opacity="0.48" {c}/>',
                f'<path d="M96 164 L75 188 M160 164 L181 188" stroke="{rule}" stroke-width="4" opacity="0.62" {c}/>',
            ]
        )
    if asset_id == "worm_model":
        return "\n".join(
            [
                f'<path d="M53 151 C82 92, 132 206, 173 116 C186 88, 204 77, 216 87" fill="none" stroke="{p}" stroke-width="16" {c}/>',
                f'<path d="M74 134 L91 142 M103 144 L119 153 M135 156 L151 164 M168 123 L184 131" stroke="{s}" stroke-width="5" opacity="0.75" {c}/>',
            ]
        )
    if asset_id == "fish_model":
        return "\n".join(
            [
                f'<path d="M62 128 C94 81, 162 78, 199 128 C162 178, 94 175, 62 128 Z" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M199 128 L229 99 V157 Z" fill="{paper}" stroke="{s}" stroke-width="8" {c}/>',
                f'<circle cx="96" cy="119" r="5" fill="{ink}"/>',
                f'<path d="M126 91 C118 113, 119 143, 132 166" fill="none" stroke="{rule}" stroke-width="5" opacity="0.55" {c}/>',
            ]
        )
    if asset_id == "plant_model":
        return "\n".join(
            [
                f'<path d="M128 199 V77" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M128 126 C86 83, 58 87, 45 129 C80 151, 109 147, 128 126 Z" fill="{paper}" stroke="{s}" stroke-width="8" {c}/>',
                f'<path d="M128 104 C168 62, 200 70, 211 110 C177 133, 148 128, 128 104 Z" fill="{paper}" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M128 158 C164 131, 194 142, 203 177 C170 190, 144 184, 128 158 Z" fill="{paper}" stroke="{s}" stroke-width="8" {c}/>',
            ]
        )
    if asset_id == "microbe_model":
        return "\n".join(
            [
                f'<ellipse cx="102" cy="119" rx="43" ry="24" fill="{paper}" stroke="{p}" stroke-width="8" transform="rotate(-24 102 119)"/>',
                f'<ellipse cx="158" cy="146" rx="48" ry="25" fill="{paper}" stroke="{s}" stroke-width="8" transform="rotate(25 158 146)"/>',
                f'<path d="M61 92 l-17 -18 M71 153 l-20 14 M137 86 l18 -19 M200 125 l23 -8 M117 180 l-4 24" stroke="{rule}" stroke-width="5" opacity="0.65" {c}/>',
                f'<circle cx="94" cy="116" r="5" fill="{p}" opacity="0.55"/><circle cx="162" cy="145" r="5" fill="{s}" opacity="0.55"/>',
            ]
        )
    if asset_id == "orbit_microgravity":
        return "\n".join(
            [
                f'<circle cx="128" cy="128" r="42" fill="{paper}" stroke="{p}" stroke-width="8"/>',
                f'<path d="M45 141 C79 71, 178 54, 214 100 C237 130, 206 170, 149 187 C98 202, 48 188, 42 157" fill="none" stroke="{s}" stroke-width="7" opacity="0.78" {c}/>',
                f'<circle cx="211" cy="103" r="9" fill="{s}"/>',
                f'<path d="M101 121 C116 110, 139 111, 155 128" fill="none" stroke="{rule}" stroke-width="4" opacity="0.55" {c}/>',
            ]
        )
    if asset_id == "radiation_stressor":
        return "\n".join(
            [
                f'<circle cx="128" cy="128" r="18" fill="{p}" opacity="0.72"/>',
                f'<path d="M128 55 V92 M128 164 V201 M55 128 H92 M164 128 H201 M76 76 L103 103 M153 153 L180 180 M180 76 L153 103 M103 153 L76 180" stroke="{s}" stroke-width="8" opacity="0.82" {c}/>',
                f'<circle cx="128" cy="128" r="72" fill="none" stroke="{rule}" stroke-width="5" opacity="0.42"/>',
            ]
        )
    if asset_id == "source_record":
        return "\n".join(
            [
                f'<rect x="70" y="50" width="116" height="154" rx="12" fill="{paper}" stroke="{p}" stroke-width="8"/>',
                f'<path d="M94 88 H162 M94 116 H162 M94 144 H139" stroke="{rule}" stroke-width="6" opacity="0.70" {c}/>',
                f'<circle cx="160" cy="177" r="18" fill="{s}" opacity="0.72"/><path d="M151 177 l7 7 l14 -17" fill="none" stroke="{paper}" stroke-width="5" {c}/>',
            ]
        )
    if asset_id == "checksum_record":
        return "\n".join(
            [
                f'<rect x="62" y="70" width="132" height="116" rx="14" fill="{paper}" stroke="{p}" stroke-width="8"/>',
                f'<path d="M90 99 h18 M121 99 h18 M152 99 h18 M90 128 h18 M121 128 h18 M152 128 h18 M90 157 h18 M121 157 h18 M152 157 h18" stroke="{rule}" stroke-width="6" opacity="0.72" {c}/>',
                f'<path d="M68 53 L91 75 M188 53 L165 75" stroke="{s}" stroke-width="7" opacity="0.76" {c}/>',
            ]
        )
    if asset_id == "split_guard":
        return "\n".join(
            [
                f'<rect x="56" y="72" width="144" height="112" rx="14" fill="{paper}" stroke="{p}" stroke-width="8"/>',
                f'<path d="M128 78 V178" stroke="{p}" stroke-width="8" opacity="0.72" {c}/>',
                f'<path d="M83 128 H112 M144 128 H173" stroke="{s}" stroke-width="7" {c}/>',
                f'<path d="M93 106 L112 128 L93 150 M163 106 L144 128 L163 150" fill="none" stroke="{s}" stroke-width="6" {c}/>',
            ]
        )
    if asset_id == "caveat_flag":
        return "\n".join(
            [
                f'<path d="M83 206 V57" stroke="{p}" stroke-width="8" {c}/>',
                f'<path d="M83 65 C117 45, 151 81, 188 59 V135 C151 157, 117 121, 83 141 Z" fill="{paper}" stroke="{s}" stroke-width="8" {c}/>',
                f'<path d="M136 82 V111" stroke="{p}" stroke-width="7" {c}/><circle cx="136" cy="128" r="5" fill="{p}"/>',
            ]
        )
    raise ValueError(f"Unknown asset_id: {asset_id}")


def symbol_group(asset_id: str, x: float, y: float, size: float, primary: str | None = None, secondary: str | None = None) -> str:
    symbol = next(item for item in SYMBOLS if item.asset_id == asset_id)
    p = primary or symbol.primary
    s = secondary or symbol.secondary
    scale = size / 256.0
    return f'<g transform="translate({x:.2f} {y:.2f}) scale({scale:.4f})">\n{icon_body(asset_id, p, s)}\n</g>'


def metadata(kind: str, payload: dict[str, str]) -> str:
    compact = json.dumps({"kind": kind, **payload}, sort_keys=True)
    return f"<metadata>{esc(compact)}</metadata>"


def svg_doc(width: int, height: int, body: str, meta: str, *, background: str = "none") -> str:
    bg = "" if background == "none" else f'<rect width="{width}" height="{height}" fill="{background}"/>'
    return f'''<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img">
  {meta}
  {bg}
{body}
</svg>
'''


def symbol_svg(asset: SymbolAsset) -> str:
    body = "\n".join(
        [
            '  <rect width="256" height="256" fill="none"/>',
            f'  <g>{icon_body(asset.asset_id, asset.primary, asset.secondary)}</g>',
        ]
    )
    return svg_doc(
        256,
        256,
        body,
        metadata(
            "biovis_symbol_v0_3",
            {
                "asset_id": asset.asset_id,
                "label": asset.label,
                "category": asset.category,
                "semantic_use": asset.semantic_use,
            },
        ),
    )


def badge_svg(asset: BadgeAsset) -> str:
    tone = col(asset.tone) if asset.tone in COLORS else col("muted")
    text_color = col("ink")
    fill = col("paper2")
    label = esc(asset.label)
    mark = {
        "generated_context": f'<circle cx="52" cy="56" r="12" fill="none" stroke="{tone}" stroke-width="5"/><path d="M43 56 H61" stroke="{tone}" stroke-width="4" {common_attrs()}/>',
        "source_proof": f'<circle cx="52" cy="56" r="13" fill="{tone}" opacity="0.80"/><path d="M45 56 l6 6 l12 -15" fill="none" stroke="{fill}" stroke-width="4.5" {common_attrs()}/>',
        "processed_result": f'<path d="M38 45 H66 V68 H38 Z" fill="none" stroke="{tone}" stroke-width="5"/><path d="M44 57 H60" stroke="{tone}" stroke-width="4" {common_attrs()}/>',
        "validated_result": f'<path d="M39 58 l10 10 l22 -26" fill="none" stroke="{tone}" stroke-width="6" {common_attrs()}/>',
        "preliminary": f'<path d="M52 39 V68" stroke="{tone}" stroke-width="6" {common_attrs()}/><circle cx="52" cy="76" r="4.5" fill="{tone}"/>',
        "caveat": f'<path d="M52 35 L73 75 H31 Z" fill="none" stroke="{tone}" stroke-width="5" {common_attrs()}/><path d="M52 49 V62" stroke="{tone}" stroke-width="4" {common_attrs()}/><circle cx="52" cy="69" r="3.5" fill="{tone}"/>',
        "train_only": f'<path d="M35 56 H69 M50 42 L69 56 L50 70" fill="none" stroke="{tone}" stroke-width="5" {common_attrs()}/>',
        "held_out": f'<rect x="36" y="39" width="32" height="34" rx="7" fill="none" stroke="{tone}" stroke-width="5"/><path d="M44 39 V31 C44 20 60 20 60 31 V39" fill="none" stroke="{tone}" stroke-width="5" {common_attrs()}/>',
    }[asset.asset_id]
    body = "\n".join(
        [
            f'  <rect x="6" y="8" width="308" height="96" rx="28" fill="{fill}" stroke="{tone}" stroke-width="3.5" opacity="0.98"/>',
            f'  {mark}',
            f'  <text x="92" y="64" font-family="Arial, Helvetica, sans-serif" font-size="34" font-weight="700" fill="{text_color}" letter-spacing="0">{label}</text>',
        ]
    )
    return svg_doc(
        320,
        112,
        body,
        metadata("biovis_badge_v0_3", {"asset_id": asset.asset_id, "label": asset.label, "semantic_use": asset.semantic_use}),
    )


def module_frame(width: int, height: int, title: str, subtitle: str) -> list[str]:
    return [
        f'<rect width="{width}" height="{height}" fill="none"/>',
        f'<text x="28" y="42" font-family="Arial, Helvetica, sans-serif" font-size="24" font-weight="700" fill="{col("ink")}" letter-spacing="0">{esc(title)}</text>',
        f'<text x="28" y="70" font-family="Arial, Helvetica, sans-serif" font-size="13" fill="{col("muted")}" letter-spacing="0">{esc(subtitle)}</text>',
        f'<path d="M28 86 H{width-28}" stroke="{col("rule")}" stroke-width="1.3" opacity="0.45"/>',
    ]


def label_text(x: int, y: int, text: str, *, size: int = 15, weight: int = 700, color: str = "ink", anchor: str = "middle") -> str:
    return f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="{size}" font-weight="{weight}" fill="{col(color)}" text-anchor="{anchor}" letter-spacing="0">{esc(text)}</text>'


def small_note(x: int, y: int, text: str, *, anchor: str = "middle") -> str:
    return f'<text x="{x}" y="{y}" font-family="Arial, Helvetica, sans-serif" font-size="11.5" fill="{col("muted")}" text-anchor="{anchor}" letter-spacing="0">{esc(text)}</text>'


def arrow(x1: int, y1: int, x2: int, y2: int, *, color: str = "rule") -> str:
    return f'<path d="M{x1} {y1} H{x2}" stroke="{col(color)}" stroke-width="2.5" opacity="0.72" {common_attrs()}/><path d="M{x2-9} {y2-7} L{x2} {y2} L{x2-9} {y2+7}" fill="none" stroke="{col(color)}" stroke-width="2.5" opacity="0.72" {common_attrs()}/>'


def module_biological_scale_ladder() -> str:
    width, height = 1280, 360
    parts = module_frame(width, height, "Biological scale ladder", "Use as compact orientation: molecular signal to organism/model scope.")
    steps = [
        ("gene_locus", "Molecule", "gene / RNA / protein"),
        ("cell_state", "Cell", "state / phenotype"),
        ("tissue_section", "Tissue", "bulk or section"),
        ("kidney_organ", "Organ", "organ-specific result"),
        ("human_model", "Model", "species / system"),
    ]
    x0, gap = 112, 242
    y_icon = 126
    for idx, (asset_id, label, note) in enumerate(steps):
        x = x0 + idx * gap
        parts.append(symbol_group(asset_id, x - 54, y_icon, 108))
        parts.append(label_text(x, 268, label))
        parts.append(small_note(x, 291, note))
        if idx < len(steps) - 1:
            parts.append(arrow(x + 79, 176, x + gap - 78, 176))
    return svg_doc(width, height, "\n".join(f"  {p}" for p in parts), metadata("biovis_module_v0_3", {"asset_id": "biological_scale_ladder"}))


def module_sample_to_feature_stack() -> str:
    width, height = 1420, 410
    parts = module_frame(width, height, "Sample to feature stack", "Reusable methods module for explaining how source material becomes model-ready features.")
    steps = [
        ("sample_tube", "Sample", "source material"),
        ("bulk_rna_assay", "Assay", "RNA readout"),
        ("omics_matrix", "Matrix", "feature table"),
        ("pathway_network", "Program", "interpretable layer"),
        ("split_guard", "Model task", "leakage guard"),
    ]
    x0, gap = 118, 276
    for idx, (asset_id, label, note) in enumerate(steps):
        x = x0 + idx * gap
        parts.append(f'<rect x="{x-70}" y="116" width="140" height="138" rx="20" fill="{col("paper2")}" stroke="{col("rule")}" stroke-width="1.2" opacity="0.82"/>')
        parts.append(symbol_group(asset_id, x - 50, 130, 100))
        parts.append(label_text(x, 292, label))
        parts.append(small_note(x, 315, note))
        if idx < len(steps) - 1:
            parts.append(arrow(x + 88, 184, x + gap - 90, 184))
    parts.append(f'<path d="M63 354 H1358" stroke="{col("rule")}" stroke-width="1.4" opacity="0.32"/>')
    parts.append(small_note(63, 380, "Replace schematic icons with source-derived proof plates when the slide makes an empirical claim.", anchor="start"))
    return svg_doc(width, height, "\n".join(f"  {p}" for p in parts), metadata("biovis_module_v0_3", {"asset_id": "sample_to_feature_stack"}))


def module_species_coverage_strip() -> str:
    width, height = 1500, 340
    parts = module_frame(width, height, "Species coverage strip", "Small scope strip for mouse/human core plus extension candidates.")
    species = [
        ("human_model", "human", "core"),
        ("organoid_rosette", "organoid", "extension"),
        ("mouse_model", "mouse", "core"),
        ("rat_model", "rat", "extension"),
        ("fly_model", "fly", "possible"),
        ("worm_model", "worm", "possible"),
        ("fish_model", "fish", "possible"),
        ("plant_model", "plant", "reserve"),
        ("microbe_model", "microbe", "reserve"),
    ]
    x0, gap = 95, 158
    for idx, (asset_id, label, status) in enumerate(species):
        x = x0 + idx * gap
        tone = {"core": "green", "extension": "blue", "possible": "amber", "reserve": "muted"}[status]
        parts.append(f'<circle cx="{x}" cy="156" r="61" fill="{col("paper2")}" stroke="{col(tone)}" stroke-width="3.2" opacity="0.98"/>')
        parts.append(symbol_group(asset_id, x - 41, 115, 82))
        parts.append(label_text(x, 250, label, size=14))
        parts.append(small_note(x, 272, status))
    return svg_doc(width, height, "\n".join(f"  {p}" for p in parts), metadata("biovis_module_v0_3", {"asset_id": "species_coverage_strip"}))


def module_modality_lane_set() -> str:
    width, height = 1380, 430
    parts = module_frame(width, height, "Modality lane set", "Compact data-modality legend for collection, processing, or benchmark-scope slides.")
    lanes = [
        ("bulk_rna_assay", "Bulk RNA", "sample-level expression", "blue"),
        ("single_cell_droplet", "Single-cell", "cell-resolved counts", "green"),
        ("spatial_spots", "Spatial", "section + location", "purple"),
        ("proteomics_assay", "Proteomics", "protein abundance", "amber"),
        ("microscopy_frame", "Imaging", "source proof slot", "rose"),
        ("organoid_rosette", "Organoid", "human model extension", "teal"),
    ]
    y0, lane_h = 108, 48
    for idx, (asset_id, label, note, tone) in enumerate(lanes):
        y = y0 + idx * lane_h
        parts.append(f'<path d="M74 {y+22} H1304" stroke="{col(tone)}" stroke-width="7" opacity="0.16" {common_attrs()}/>')
        parts.append(symbol_group(asset_id, 70, y - 8, 58))
        parts.append(label_text(154, y + 29, label, size=16, anchor="start"))
        parts.append(small_note(342, y + 29, note, anchor="start"))
        parts.append(f'<path d="M580 {y+22} H1220" stroke="{col("rule")}" stroke-width="1.2" stroke-dasharray="3 8" opacity="0.36"/>')
        parts.append(f'<circle cx="1262" cy="{y+22}" r="7" fill="{col(tone)}" opacity="0.64"/>')
    return svg_doc(width, height, "\n".join(f"  {p}" for p in parts), metadata("biovis_module_v0_3", {"asset_id": "modality_lane_set"}))


def module_trust_status_ribbon() -> str:
    width, height = 1400, 360
    parts = module_frame(width, height, "Trust status ribbon", "Reusable governance strip for source provenance, freezes, split hygiene, and caveats.")
    steps = [
        ("source_record", "source", "recorded"),
        ("checksum_record", "freeze", "hash / manifest"),
        ("split_guard", "split", "train-only guard"),
        ("omics_matrix", "result", "processed output"),
        ("caveat_flag", "caveat", "claim boundary"),
    ]
    x0, gap = 124, 276
    for idx, (asset_id, label, note) in enumerate(steps):
        x = x0 + idx * gap
        parts.append(f'<circle cx="{x}" cy="170" r="58" fill="{col("paper2")}" stroke="{col("rule")}" stroke-width="1.2" opacity="0.94"/>')
        parts.append(symbol_group(asset_id, x - 42, 128, 84))
        parts.append(label_text(x, 267, label))
        parts.append(small_note(x, 291, note))
        if idx < len(steps) - 1:
            parts.append(f'<path d="M{x+68} 170 H{x+gap-68}" stroke="{col("rule")}" stroke-width="2.2" opacity="0.60" {common_attrs()}/>')
    parts.append(f'<rect x="36" y="320" width="1328" height="1.2" fill="{col("rule")}" opacity="0.32"/>')
    return svg_doc(width, height, "\n".join(f"  {p}" for p in parts), metadata("biovis_module_v0_3", {"asset_id": "trust_status_ribbon"}))


def module_claim_boundary_bar() -> str:
    width, height = 1360, 300
    parts = module_frame(width, height, "Claim boundary bar", "Footer-scale module: keep schematic context separate from proof and validation.")
    items = [
        ("generated_context", "Generated context", "visual orientation only"),
        ("source_proof", "Source proof", "real source-derived object"),
        ("processed_result", "Processed result", "computed analysis output"),
        ("validated_result", "Validated claim", "reviewed claim surface"),
    ]
    x0, gap = 112, 312
    for idx, (badge_id, title, note) in enumerate(items):
        badge = next(item for item in BADGES if item.asset_id == badge_id)
        tone = col(badge.tone)
        x = x0 + idx * gap
        parts.append(f'<rect x="{x-70}" y="112" width="140" height="54" rx="20" fill="{col("paper2")}" stroke="{tone}" stroke-width="2.4"/>')
        parts.append(f'<text x="{x}" y="147" font-family="Arial, Helvetica, sans-serif" font-size="17" font-weight="700" fill="{col("ink")}" text-anchor="middle" letter-spacing="0">{esc(badge.label)}</text>')
        parts.append(label_text(x, 202, title, size=14))
        parts.append(small_note(x, 225, note))
        if idx < len(items) - 1:
            parts.append(arrow(x + 88, 139, x + gap - 90, 139))
    return svg_doc(width, height, "\n".join(f"  {p}" for p in parts), metadata("biovis_module_v0_3", {"asset_id": "claim_boundary_bar"}))


def module_space_bio_assay_lane() -> str:
    width, height = 1480, 370
    parts = module_frame(width, height, "Space biology assay lane", "Compact narrative lane: mission context to tissue, assay, feature layer, and split guard.")
    steps = [
        ("orbit_microgravity", "Mission", "spaceflight context"),
        ("radiation_stressor", "Stressor", "condition / exposure"),
        ("tissue_section", "Tissue", "sample source"),
        ("bulk_rna_assay", "Assay", "omics readout"),
        ("pathway_network", "Feature layer", "program summary"),
        ("split_guard", "Benchmark", "held-out guard"),
    ]
    x0, gap = 95, 247
    parts.append(f'<path d="M83 176 H1372" stroke="{col("blue")}" stroke-width="10" opacity="0.10" {common_attrs()}/>')
    for idx, (asset_id, label, note) in enumerate(steps):
        x = x0 + idx * gap
        parts.append(f'<circle cx="{x}" cy="176" r="55" fill="{col("paper2")}" stroke="{col("rule")}" stroke-width="1.3" opacity="0.96"/>')
        parts.append(symbol_group(asset_id, x - 39, 137, 78))
        parts.append(label_text(x, 272, label, size=14))
        parts.append(small_note(x, 294, note))
        if idx < len(steps) - 1:
            parts.append(f'<path d="M{x+63} 176 H{x+gap-63}" stroke="{col("rule")}" stroke-width="2.0" opacity="0.55" {common_attrs()}/>')
    return svg_doc(width, height, "\n".join(f"  {p}" for p in parts), metadata("biovis_module_v0_3", {"asset_id": "space_bio_assay_lane"}))


def module_svg(asset: ModuleAsset) -> str:
    builders = {
        "biological_scale_ladder": module_biological_scale_ladder,
        "sample_to_feature_stack": module_sample_to_feature_stack,
        "species_coverage_strip": module_species_coverage_strip,
        "modality_lane_set": module_modality_lane_set,
        "trust_status_ribbon": module_trust_status_ribbon,
        "claim_boundary_bar": module_claim_boundary_bar,
        "space_bio_assay_lane": module_space_bio_assay_lane,
    }
    return builders[asset.asset_id]()


def dark_module_svg(svg_text: str) -> str:
    replacements = {
        COLORS["ink"]: "#F2F6F8",
        COLORS["muted"]: "#B8C4CF",
        COLORS["rule"]: "#617282",
        COLORS["paper"]: "#0D1720",
        COLORS["paper2"]: "#111D27",
        COLORS["blue"]: "#75B7E0",
        COLORS["sky"]: "#9AC9E8",
        COLORS["green"]: "#35B989",
        COLORS["teal"]: "#42C7B8",
        COLORS["amber"]: "#D7A64A",
        COLORS["rose"]: "#D982A3",
        COLORS["purple"]: "#A28AE8",
        COLORS["red"]: "#DF6F6F",
    }
    dark = svg_text
    dark = dark.replace('"biovis_module_v0_3"', '"biovis_module_v0_3_dark"')
    for source, target in replacements.items():
        dark = dark.replace(source, target)
    dark = dark.replace('fill="none"/>', f'fill="{COLORS["deep"]}"/>', 1)
    return dark


def run_rsvg(svg_path: Path, png_path: Path, *, width: int | None = None, height: int | None = None) -> bool:
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


def grayscale_copy(src: Path, dst: Path) -> bool:
    try:
        from PIL import Image
    except ImportError:
        return False
    with Image.open(src) as image:
        gray = image.convert("LA").convert("RGBA")
        gray.save(dst)
    return True


def load_font(size: int, *, bold: bool = False):
    try:
        from PIL import ImageFont
    except ImportError:  # pragma: no cover
        return None
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf" if bold else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/Library/Fonts/Arial.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except OSError:
            continue
    return ImageFont.load_default()


def hex_to_rgb(value: str) -> tuple[int, int, int]:
    value = value.lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def draw_centered(draw, xy: tuple[int, int], text: str, font, fill: str) -> None:
    bbox = draw.textbbox((0, 0), text, font=font)
    width = bbox[2] - bbox[0]
    draw.text((xy[0] - width / 2, xy[1]), text, font=font, fill=hex_to_rgb(col(fill)))


def paste_fit(canvas, image_path: Path, box: tuple[int, int, int, int]) -> None:
    from PIL import Image

    image = Image.open(image_path).convert("RGBA")
    x, y, w, h = box
    scale = min(w / image.width, h / image.height)
    new_size = (max(1, int(image.width * scale)), max(1, int(image.height * scale)))
    image = image.resize(new_size, Image.Resampling.LANCZOS)
    px = x + (w - image.width) // 2
    py = y + (h - image.height) // 2
    canvas.alpha_composite(image, (px, py))


def write_symbol_contact_png(path: Path) -> bool:
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return False
    cols = 6
    cell_w = 244
    cell_h = 242
    pad = 52
    top = 130
    rows = (len(SYMBOLS) + cols - 1) // cols
    width = cols * cell_w + pad * 2
    height = top + rows * cell_h + 70
    canvas = Image.new("RGBA", (width, height), hex_to_rgb(col("paper")) + (255,))
    draw = ImageDraw.Draw(canvas)
    title_font = load_font(34, bold=True)
    body_font = load_font(16)
    label_font = load_font(14, bold=True)
    small_font = load_font(11)
    draw.text((pad, 34), "Biological symbol module pack v0.3: symbols", font=title_font, fill=hex_to_rgb(col("ink")))
    draw.text((pad, 78), "Small editable vectors for biological scale, assay modality, species/model systems, and trust cues.", font=body_font, fill=hex_to_rgb(col("muted")))
    for idx, asset in enumerate(SYMBOLS):
        col_i = idx % cols
        row_i = idx // cols
        x = pad + col_i * cell_w
        y = top + row_i * cell_h
        draw.rounded_rectangle((x, y, x + 196, y + 196), radius=20, fill=hex_to_rgb(col("paper2")), outline=hex_to_rgb(col("rule")), width=1)
        paste_fit(canvas, SYMBOL_PNG / f"{asset.asset_id}.png", (x + 44, y + 18, 108, 108))
        draw_centered(draw, (x + 98, y + 140), asset.label, label_font, "ink")
        draw_centered(draw, (x + 98, y + 164), asset.category, small_font, "muted")
    canvas.convert("RGB").save(path)
    return True


def write_badge_module_contact_png(path: Path) -> bool:
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return False
    width, height = 1800, 2200
    canvas = Image.new("RGBA", (width, height), hex_to_rgb(col("paper")) + (255,))
    draw = ImageDraw.Draw(canvas)
    title_font = load_font(36, bold=True)
    body_font = load_font(17)
    h2_font = load_font(24, bold=True)
    label_font = load_font(21, bold=True)
    note_font = load_font(14)
    small_font = load_font(13)
    draw.text((60, 36), "Biological symbol module pack v0.3: badges and modules", font=title_font, fill=hex_to_rgb(col("ink")))
    draw.text((60, 84), "Reusable components for consulting-grade methods, scope, and evidence-status communication.", font=body_font, fill=hex_to_rgb(col("muted")))
    draw.text((60, 146), "Status badges", font=h2_font, fill=hex_to_rgb(col("ink")))
    for idx, badge in enumerate(BADGES):
        x = 60 + (idx % 4) * 420
        y = 180 + (idx // 4) * 128
        paste_fit(canvas, BADGE_PNG / f"{badge.asset_id}.png", (x, y, 340, 96))

    draw.text((60, 472), "Modules", font=h2_font, fill=hex_to_rgb(col("ink")))
    y = 518
    for module in MODULES:
        draw.rounded_rectangle((60, y, width - 60, y + 208), radius=18, fill=hex_to_rgb(col("paper2")), outline=hex_to_rgb(col("rule")), width=1)
        paste_fit(canvas, MODULE_PNG / f"{module.asset_id}.png", (92, y + 20, 860, 168))
        draw.text((1000, y + 46), module.label, font=label_font, fill=hex_to_rgb(col("ink")))
        draw.text((1000, y + 82), module.semantic_use, font=note_font, fill=hex_to_rgb(col("muted")))
        draw.text((1000, y + 116), f"type: {module.module_type}", font=small_font, fill=hex_to_rgb(col("muted")))
        draw.text((1000, y + 146), f"use: {module.allowed_use}", font=small_font, fill=hex_to_rgb(col("muted")))
        y += 230
    canvas.convert("RGB").save(path)
    return True


def write_dark_module_contact_png(path: Path) -> bool:
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return False
    width, height = 1800, 1820
    canvas = Image.new("RGBA", (width, height), hex_to_rgb(col("deep")) + (255,))
    draw = ImageDraw.Draw(canvas)
    title_font = load_font(36, bold=True)
    body_font = load_font(17)
    label_font = load_font(21, bold=True)
    note_font = load_font(14)
    draw.text((60, 38), "Biological symbol module pack v0.3: dark modules", font=title_font, fill=(242, 246, 248))
    draw.text((60, 86), "Dark-field variants for premium slides that use deep canvas, layered texture, or source-image backgrounds.", font=body_font, fill=(184, 196, 207))
    y = 142
    for module in MODULES:
        draw.rounded_rectangle((60, y, width - 60, y + 210), radius=18, fill=(17, 29, 39), outline=(97, 114, 130), width=1)
        paste_fit(canvas, MODULE_DARK_PNG / f"{module.asset_id}.png", (92, y + 20, 900, 170))
        draw.text((1040, y + 52), module.label, font=label_font, fill=(242, 246, 248))
        draw.text((1040, y + 88), module.semantic_use, font=note_font, fill=(184, 196, 207))
        draw.text((1040, y + 122), f"type: {module.module_type}", font=note_font, fill=(184, 196, 207))
        y += 230
    canvas.convert("RGB").save(path)
    return True


def write_micro_icon_qa_png(path: Path) -> bool:
    try:
        from PIL import Image, ImageDraw
    except ImportError:
        return False
    cols = 4
    row_h = 86
    pad = 44
    cell_w = 410
    card_w = 384
    width = pad * 2 + cols * cell_w
    rows = (len(SYMBOLS) + cols - 1) // cols
    height = 130 + rows * row_h + 64
    canvas = Image.new("RGBA", (width, height), hex_to_rgb(col("paper")) + (255,))
    draw = ImageDraw.Draw(canvas)
    title_font = load_font(30, bold=True)
    body_font = load_font(15)
    label_font = load_font(13, bold=True)
    small_font = load_font(10)
    draw.text((pad, 34), "Micro-size symbol QA", font=title_font, fill=hex_to_rgb(col("ink")))
    draw.text((pad, 72), "Each symbol is shown at 32, 48, and 64 px. Use labels when the 32 px version is ambiguous.", font=body_font, fill=hex_to_rgb(col("muted")))
    for idx, asset in enumerate(SYMBOLS):
        col_i = idx % cols
        row_i = idx // cols
        x = pad + col_i * cell_w
        y = 126 + row_i * row_h
        draw.rounded_rectangle((x, y, x + card_w, y + 64), radius=12, fill=hex_to_rgb(col("paper2")), outline=hex_to_rgb(col("rule")), width=1)
        draw.text((x + 12, y + 13), asset.label, font=label_font, fill=hex_to_rgb(col("ink")))
        draw.text((x + 12, y + 35), asset.category, font=small_font, fill=hex_to_rgb(col("muted")))
        paste_fit(canvas, SYMBOL_PNG / f"{asset.asset_id}.png", (x + 208, y + 16, 32, 32))
        paste_fit(canvas, SYMBOL_PNG / f"{asset.asset_id}.png", (x + 258, y + 8, 48, 48))
        paste_fit(canvas, SYMBOL_PNG / f"{asset.asset_id}.png", (x + 320, y, 64, 64))
    canvas.convert("RGB").save(path)
    return True


def contact_sheet_symbols() -> str:
    cols = 6
    cell_w = 212
    cell_h = 218
    pad = 44
    top = 120
    rows = (len(SYMBOLS) + cols - 1) // cols
    width = cols * cell_w + pad * 2
    height = top + rows * cell_h + 58
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="{col("paper")}"/>',
        f'<text x="{pad}" y="52" font-family="Arial, Helvetica, sans-serif" font-size="30" font-weight="700" fill="{col("ink")}" letter-spacing="0">Biological symbol module pack v0.3: symbols</text>',
        f'<text x="{pad}" y="85" font-family="Arial, Helvetica, sans-serif" font-size="15" fill="{col("muted")}" letter-spacing="0">Small editable vectors for biological scale, assay modality, species/model systems, and trust cues.</text>',
    ]
    for idx, asset in enumerate(SYMBOLS):
        col_i = idx % cols
        row_i = idx // cols
        x = pad + col_i * cell_w
        y = top + row_i * cell_h
        parts.append(f'<rect x="{x}" y="{y}" width="174" height="178" rx="18" fill="{col("paper2")}" stroke="{col("rule")}" stroke-width="1" opacity="0.65"/>')
        parts.append(f'<image href="../symbols/png/{asset.asset_id}.png" x="{x+38}" y="{y+12}" width="98" height="98" preserveAspectRatio="xMidYMid meet"/>')
        parts.append(f'<text x="{x+87}" y="{y+138}" font-family="Arial, Helvetica, sans-serif" font-size="13" font-weight="700" fill="{col("ink")}" text-anchor="middle" letter-spacing="0">{esc(asset.label)}</text>')
        parts.append(f'<text x="{x+87}" y="{y+160}" font-family="Arial, Helvetica, sans-serif" font-size="10.5" fill="{col("muted")}" text-anchor="middle" letter-spacing="0">{esc(asset.category)}</text>')
    parts.append("</svg>\n")
    return "\n".join(parts)


def contact_sheet_badges_modules() -> str:
    width, height = 1600, 1700
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="{col("paper")}"/>',
        f'<text x="54" y="58" font-family="Arial, Helvetica, sans-serif" font-size="32" font-weight="700" fill="{col("ink")}" letter-spacing="0">Biological symbol module pack v0.3: badges and modules</text>',
        f'<text x="54" y="94" font-family="Arial, Helvetica, sans-serif" font-size="16" fill="{col("muted")}" letter-spacing="0">Reusable components for consulting-grade methods, scope, and evidence-status communication.</text>',
        f'<text x="54" y="150" font-family="Arial, Helvetica, sans-serif" font-size="22" font-weight="700" fill="{col("ink")}" letter-spacing="0">Status badges</text>',
    ]
    for idx, badge in enumerate(BADGES):
        x = 54 + (idx % 4) * 374
        y = 174 + (idx // 4) * 140
        parts.append(f'<image href="../badges/png/{badge.asset_id}.png" x="{x}" y="{y}" width="300" height="105" preserveAspectRatio="xMidYMid meet"/>')
    parts.append(f'<text x="54" y="486" font-family="Arial, Helvetica, sans-serif" font-size="22" font-weight="700" fill="{col("ink")}" letter-spacing="0">Modules</text>')
    y = 522
    for module in MODULES:
        parts.append(f'<rect x="54" y="{y}" width="1492" height="142" rx="18" fill="{col("paper2")}" stroke="{col("rule")}" stroke-width="1" opacity="0.72"/>')
        parts.append(f'<image href="../modules/png/{module.asset_id}.png" x="80" y="{y+14}" width="720" height="114" preserveAspectRatio="xMinYMid meet"/>')
        parts.append(f'<text x="838" y="{y+48}" font-family="Arial, Helvetica, sans-serif" font-size="20" font-weight="700" fill="{col("ink")}" letter-spacing="0">{esc(module.label)}</text>')
        parts.append(f'<text x="838" y="{y+76}" font-family="Arial, Helvetica, sans-serif" font-size="13" fill="{col("muted")}" letter-spacing="0">{esc(module.semantic_use)}</text>')
        parts.append(f'<text x="838" y="{y+103}" font-family="Arial, Helvetica, sans-serif" font-size="12" fill="{col("muted")}" letter-spacing="0">type: {esc(module.module_type)}</text>')
        y += 158
    parts.append("</svg>\n")
    return "\n".join(parts)


def ensure_dirs() -> None:
    for directory in [OUT, SYMBOL_SVG, SYMBOL_PNG, BADGE_SVG, BADGE_PNG, MODULE_SVG, MODULE_PNG, MODULE_DARK_SVG, MODULE_DARK_PNG, PREVIEW, QA]:
        directory.mkdir(parents=True, exist_ok=True)


def write_manifest() -> None:
    symbol_records = []
    for asset in SYMBOLS:
        symbol_records.append(
            {
                "asset_id": asset.asset_id,
                "asset_type": "symbol",
                "category": asset.category,
                "label": asset.label,
                "semantic_use": asset.semantic_use,
                "primary_color_token": asset.primary,
                "secondary_color_token": asset.secondary,
                "svg": rel(SYMBOL_SVG / f"{asset.asset_id}.svg"),
                "png_preview": rel(SYMBOL_PNG / f"{asset.asset_id}.png"),
                "allowed_use": asset.allowed_use,
                "prohibited_use": asset.prohibited_use,
                "license": "Original project artwork; repository MIT license.",
            }
        )
    badge_records = []
    for asset in BADGES:
        badge_records.append(
            {
                "asset_id": asset.asset_id,
                "asset_type": "badge",
                "category": "trust-status",
                "label": asset.label,
                "semantic_use": asset.semantic_use,
                "primary_color_token": asset.tone,
                "secondary_color_token": "",
                "svg": rel(BADGE_SVG / f"{asset.asset_id}.svg"),
                "png_preview": rel(BADGE_PNG / f"{asset.asset_id}.png"),
                "allowed_use": asset.allowed_use,
                "prohibited_use": asset.prohibited_use,
                "license": "Original project artwork; repository MIT license.",
            }
        )
    module_records = []
    for asset in MODULES:
        module_records.append(
            {
                "asset_id": asset.asset_id,
                "asset_type": "module",
                "category": asset.module_type,
                "label": asset.label,
                "semantic_use": asset.semantic_use,
                "primary_color_token": "",
                "secondary_color_token": "",
                "svg": rel(MODULE_SVG / f"{asset.asset_id}.svg"),
                "png_preview": rel(MODULE_PNG / f"{asset.asset_id}.png"),
                "allowed_use": asset.allowed_use,
                "prohibited_use": asset.prohibited_use,
                "license": "Original project artwork; repository MIT license.",
            }
        )
    module_dark_records = []
    for asset in MODULES:
        module_dark_records.append(
            {
                "asset_id": f"{asset.asset_id}_dark",
                "asset_type": "module_dark",
                "category": asset.module_type,
                "label": f"{asset.label} dark",
                "semantic_use": f"{asset.semantic_use}; dark-field variant",
                "primary_color_token": "deep",
                "secondary_color_token": "",
                "svg": rel(MODULE_DARK_SVG / f"{asset.asset_id}.svg"),
                "png_preview": rel(MODULE_DARK_PNG / f"{asset.asset_id}.png"),
                "allowed_use": f"{asset.allowed_use}; dark-canvas slides",
                "prohibited_use": asset.prohibited_use,
                "license": "Original project artwork; repository MIT license.",
            }
        )
    records = symbol_records + badge_records + module_records + module_dark_records
    manifest = {
        "created": CREATED,
        "version": "0.3",
        "scope": "Reusable biological symbols, status badges, light modules, and dark-field modules for SpaceBio-Bench v0.3 visual system.",
        "design_rule": "Symbols explain biological or evidence semantics; they must not substitute for source-derived proof imagery.",
        "asset_counts": {
            "symbols": len(SYMBOLS),
            "badges": len(BADGES),
            "modules_light": len(MODULES),
            "modules_dark": len(MODULES),
            "total": len(records),
        },
        "preview": {
            "symbol_contact_sheet_svg": rel(PREVIEW / "biovis_symbol_module_pack_v0_3_symbols.svg"),
            "symbol_contact_sheet_png": rel(PREVIEW / "biovis_symbol_module_pack_v0_3_symbols.png"),
            "badge_module_contact_sheet_svg": rel(PREVIEW / "biovis_symbol_module_pack_v0_3_badges_modules.svg"),
            "badge_module_contact_sheet_png": rel(PREVIEW / "biovis_symbol_module_pack_v0_3_badges_modules.png"),
            "dark_module_contact_sheet_png": rel(PREVIEW / "biovis_symbol_module_pack_v0_3_dark_modules.png"),
            "micro_symbol_qa_png": rel(QA / "biovis_symbol_module_pack_v0_3_micro_symbol_qa.png"),
            "grayscale_symbol_contact_sheet": rel(QA / "biovis_symbol_module_pack_v0_3_symbols_grayscale.png"),
            "grayscale_badge_module_contact_sheet": rel(QA / "biovis_symbol_module_pack_v0_3_badges_modules_grayscale.png"),
        },
        "records": records,
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with (OUT / "manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)


def write_qa(conversion_ok: bool, grayscale_ok: bool) -> None:
    expected = len(SYMBOLS) + len(BADGES) + len(MODULES) * 2
    svg_paths = list(SYMBOL_SVG.glob("*.svg")) + list(BADGE_SVG.glob("*.svg")) + list(MODULE_SVG.glob("*.svg")) + list(MODULE_DARK_SVG.glob("*.svg"))
    png_paths = list(SYMBOL_PNG.glob("*.png")) + list(BADGE_PNG.glob("*.png")) + list(MODULE_PNG.glob("*.png")) + list(MODULE_DARK_PNG.glob("*.png"))
    qa = {
        "created": CREATED,
        "automatic_checks": {
            "expected_asset_count": expected,
            "svg_asset_count": len(svg_paths),
            "png_preview_count": len(png_paths),
            "rsvg_conversion_available": conversion_ok,
            "all_svg_assets_created": len(svg_paths) == expected,
            "all_png_previews_created": len(png_paths) == expected,
            "manifest_json_exists": (OUT / "manifest.json").exists(),
            "manifest_csv_exists": (OUT / "manifest.csv").exists(),
            "grayscale_previews_created": grayscale_ok,
            "dark_module_contact_sheet_exists": (PREVIEW / "biovis_symbol_module_pack_v0_3_dark_modules.png").exists(),
            "micro_symbol_qa_exists": (QA / "biovis_symbol_module_pack_v0_3_micro_symbol_qa.png").exists(),
        },
        "manual_review_template": {
            "symbol_readability": "inspect at 64px, 128px, and contact-sheet scale",
            "micro_symbol_readability": "inspect 32px, 48px, and 64px icon samples before footer use",
            "module_readability": "inspect labels and arrows at slide-footer and half-width module scale",
            "dark_module_readability": "inspect on deep canvas; use dark modules only on dark-field slides",
            "claim_boundary": "symbols and badges are semantic signposts, not proof images",
            "design_risk": "avoid badge clutter and avoid turning modules into dense mini-slides",
        },
    }
    (QA / "biovis_symbol_module_pack_v0_3_qa.json").write_text(json.dumps(qa, indent=2) + "\n")


def main() -> None:
    ensure_dirs()

    conversion_ok = shutil.which("rsvg-convert") is not None

    for asset in SYMBOLS:
        path = SYMBOL_SVG / f"{asset.asset_id}.svg"
        path.write_text(symbol_svg(asset))
        if conversion_ok:
            run_rsvg(path, SYMBOL_PNG / f"{asset.asset_id}.png", width=512, height=512)

    for asset in BADGES:
        path = BADGE_SVG / f"{asset.asset_id}.svg"
        path.write_text(badge_svg(asset))
        if conversion_ok:
            run_rsvg(path, BADGE_PNG / f"{asset.asset_id}.png", width=640)

    for asset in MODULES:
        path = MODULE_SVG / f"{asset.asset_id}.svg"
        light_svg = module_svg(asset)
        path.write_text(light_svg)
        if conversion_ok:
            run_rsvg(path, MODULE_PNG / f"{asset.asset_id}.png", width=1800)
        dark_path = MODULE_DARK_SVG / f"{asset.asset_id}.svg"
        dark_path.write_text(dark_module_svg(light_svg))
        if conversion_ok:
            run_rsvg(dark_path, MODULE_DARK_PNG / f"{asset.asset_id}.png", width=1800)

    symbol_contact = PREVIEW / "biovis_symbol_module_pack_v0_3_symbols.svg"
    symbol_contact.write_text(contact_sheet_symbols())
    badge_module_contact = PREVIEW / "biovis_symbol_module_pack_v0_3_badges_modules.svg"
    badge_module_contact.write_text(contact_sheet_badges_modules())
    symbol_contact_png = PREVIEW / "biovis_symbol_module_pack_v0_3_symbols.png"
    badge_module_contact_png = PREVIEW / "biovis_symbol_module_pack_v0_3_badges_modules.png"
    wrote_symbol_contact = write_symbol_contact_png(symbol_contact_png)
    wrote_badge_module_contact = write_badge_module_contact_png(badge_module_contact_png)
    write_dark_module_contact_png(PREVIEW / "biovis_symbol_module_pack_v0_3_dark_modules.png")
    write_micro_icon_qa_png(QA / "biovis_symbol_module_pack_v0_3_micro_symbol_qa.png")
    if conversion_ok and not wrote_symbol_contact:
        run_rsvg(symbol_contact, PREVIEW / "biovis_symbol_module_pack_v0_3_symbols.png")
    if conversion_ok and not wrote_badge_module_contact:
        run_rsvg(badge_module_contact, PREVIEW / "biovis_symbol_module_pack_v0_3_badges_modules.png")

    grayscale_ok = False
    if (PREVIEW / "biovis_symbol_module_pack_v0_3_symbols.png").exists():
        grayscale_ok = grayscale_copy(PREVIEW / "biovis_symbol_module_pack_v0_3_symbols.png", QA / "biovis_symbol_module_pack_v0_3_symbols_grayscale.png")
    if (PREVIEW / "biovis_symbol_module_pack_v0_3_badges_modules.png").exists():
        grayscale_ok = grayscale_copy(PREVIEW / "biovis_symbol_module_pack_v0_3_badges_modules.png", QA / "biovis_symbol_module_pack_v0_3_badges_modules_grayscale.png") and grayscale_ok

    write_manifest()
    write_qa(conversion_ok, grayscale_ok)

    print(
        json.dumps(
            {
                "output": rel(OUT),
                "symbols": len(SYMBOLS),
                "badges": len(BADGES),
                "modules_light": len(MODULES),
                "modules_dark": len(MODULES),
                "conversion_ok": conversion_ok,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()

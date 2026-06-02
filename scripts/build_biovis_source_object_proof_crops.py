#!/usr/bin/env python3
"""Build real source-object proof crops for v0.4 BioVis modules.

These crops replace schematic placeholders with local, auditable evidence:
actual normalized-count matrix values and manifest-derived split counts.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output" / "biovis_source_object_proof_crops_v0_1"
CROPS = OUT / "proof_crops"
MOCKS = OUT / "production_mocks"
CREATED = "2026-06-02"

COLORS = {
    "deep": "#0D1720",
    "deep2": "#111D27",
    "deep3": "#172636",
    "white": "#F2F6F8",
    "soft": "#B8C4CF",
    "muted": "#6D7D8D",
    "blue": "#2D6F9F",
    "sky": "#6BAED6",
    "green": "#178B63",
    "teal": "#1AA090",
    "amber": "#B4832F",
    "rose": "#B45A7E",
    "red": "#B23A3A",
    "purple": "#6750A4",
}


def rgb(token: str) -> tuple[int, int, int]:
    value = COLORS[token].lstrip("#")
    return tuple(int(value[i : i + 2], 16) for i in (0, 2, 4))


def rgba(token: str, alpha: int) -> tuple[int, int, int, int]:
    return rgb(token) + (alpha,)


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load_json(path: str) -> object:
    with (ROOT / path).open() as handle:
        return json.load(handle)


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


F = {
    "eyebrow": font(20, bold=True),
    "title": font(44, bold=True),
    "subtitle": font(24),
    "h1": font(34, bold=True),
    "h2": font(26, bold=True),
    "body": font(20),
    "small": font(15),
    "tiny": font(12),
}


def ensure() -> None:
    for directory in [OUT, CROPS, MOCKS]:
        directory.mkdir(parents=True, exist_ok=True)


def write_wrapped(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    text: str,
    font_obj,
    fill,
    *,
    width: int,
    line_height: int,
    max_lines: int,
) -> int:
    words = text.split()
    lines: list[str] = []
    current = ""
    for word in words:
        candidate = f"{current} {word}".strip()
        if draw.textlength(candidate, font=font_obj) <= width:
            current = candidate
        else:
            if current:
                lines.append(current)
            current = word
    if current:
        lines.append(current)
    if len(lines) > max_lines:
        lines = lines[:max_lines]
        lines[-1] = lines[-1].rstrip(".") + "..."
    x, y = xy
    for idx, line in enumerate(lines):
        draw.text((x, y + idx * line_height), line, font=font_obj, fill=fill)
    return y + len(lines) * line_height


def compact_sha(value: str) -> str:
    return f"{value[:10]}...{value[-8:]}" if value else "NA"


def badge(draw: ImageDraw.ImageDraw, xy: tuple[int, int], text: str, tone: str, *, fill_token: str = "deep2") -> None:
    x, y = xy
    width = max(104, int(draw.textlength(text, font=F["small"])) + 44)
    draw.rounded_rectangle((x, y, x + width, y + 34), radius=14, fill=rgba(fill_token, 240), outline=rgb(tone), width=2)
    draw.ellipse((x + 14, y + 12, x + 24, y + 22), fill=rgb(tone))
    draw.text((x + 34, y + 9), text, font=F["small"], fill=rgb("white"))


def card(draw: ImageDraw.ImageDraw, xywh: tuple[int, int, int, int], *, alpha: int = 236, outline: str = "deep3", radius: int = 24) -> None:
    x, y, w, h = xywh
    draw.rounded_rectangle((x, y, x + w, y + h), radius=radius, fill=rgba("deep2", alpha), outline=rgba(outline, 220), width=2)


def read_matrix_crop(path: Path, *, rows: int = 42) -> tuple[list[str], list[str], list[list[float]]]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        samples = header[1:]
        genes: list[str] = []
        values: list[list[float]] = []
        for row in reader:
            if len(row) < 2:
                continue
            try:
                numeric = [math.log2(float(value) + 1.0) for value in row[1:]]
            except ValueError:
                continue
            genes.append(row[0].strip('"'))
            values.append(numeric)
            if len(values) >= rows:
                break
    return genes, samples, values


def value_color(value: float, low: float, high: float) -> tuple[int, int, int]:
    if high <= low:
        t = 0.5
    else:
        t = max(0.0, min(1.0, (value - low) / (high - low)))
    stops = [
        (0.0, (20, 35, 51)),
        (0.32, rgb("blue")),
        (0.62, rgb("teal")),
        (0.82, rgb("amber")),
        (1.0, rgb("rose")),
    ]
    for idx in range(len(stops) - 1):
        t0, c0 = stops[idx]
        t1, c1 = stops[idx + 1]
        if t0 <= t <= t1:
            local = (t - t0) / (t1 - t0)
            return tuple(int(c0[channel] + (c1[channel] - c0[channel]) * local) for channel in range(3))
    return stops[-1][1]


def draw_heatmap(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    values: list[list[float]],
    *,
    cell_w: int,
    cell_h: int,
    gap: int = 2,
) -> None:
    flat = [value for row in values for value in row]
    flat_sorted = sorted(flat)
    lo = flat_sorted[int(len(flat_sorted) * 0.04)]
    hi = flat_sorted[int(len(flat_sorted) * 0.96)]
    x0, y0 = xy
    for r, row in enumerate(values):
        row_mean = sum(row) / len(row)
        centered = [value - row_mean for value in row]
        for c, value in enumerate(centered):
            color = value_color(value, lo - row_mean, hi - row_mean)
            x = x0 + c * (cell_w + gap)
            y = y0 + r * (cell_h + gap)
            draw.rounded_rectangle((x, y, x + cell_w, y + cell_h), radius=2, fill=color)


def render_organoid_matrix_crop() -> Path:
    manifest = load_json("v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json")
    matrices = manifest["split"]["expression_matrix_sources"]
    datasets = [
        ("OSD-863", "cortical organoid", ROOT / matrices["OSD-863"]["local_matrix_path"], "sky"),
        ("OSD-871", "dopaminergic organoid", ROOT / matrices["OSD-871"]["local_matrix_path"], "teal"),
    ]
    canvas = Image.new("RGBA", (2200, 1250), rgb("deep") + (255,))
    draw = ImageDraw.Draw(canvas)
    draw.text((70, 54), "Organoid matrix proof", font=F["title"], fill=rgb("white"))
    draw.text((70, 108), "Actual normalized-count matrix crops with manifest-derived counts and hashes.", font=F["subtitle"], fill=rgb("soft"))
    draw.line((70, 158, 2130, 158), fill=rgba("deep3", 255), width=2)

    x_positions = [70, 1130]
    for idx, (source_id, label, matrix_path, tone) in enumerate(datasets):
        x = x_positions[idx]
        entry = matrices[source_id]
        genes, samples, values = read_matrix_crop(matrix_path)
        card(draw, (x, 210, 1000, 830), alpha=220)
        draw.text((x + 36, 244), source_id, font=F["h1"], fill=rgb("white"))
        draw.text((x + 36, 284), label, font=F["body"], fill=rgb("soft"))
        badge(draw, (x + 36, 322), "matrix downloaded", tone)
        badge(draw, (x + 236, 322), "sample-aligned", "green")
        heat_x, heat_y = x + 46, 386
        draw_heatmap(draw, (heat_x, heat_y), values, cell_w=32, cell_h=10, gap=2)
        draw.rectangle((heat_x, heat_y, heat_x + len(samples) * 34 - 2, heat_y + len(values) * 12 - 2), outline=rgba("soft", 120), width=1)
        draw.text((x + 48, 920), f"{entry['n_feature_rows']:,} genes x {entry['n_sample_columns']} samples", font=F["h2"], fill=rgb("white"))
        draw.text((x + 48, 962), f"sha256 {compact_sha(entry['sha256'])}", font=F["body"], fill=rgb("soft"))
        draw.text((x + 48, 1002), Path(entry["local_matrix_path"]).name, font=F["small"], fill=rgb("muted"))
        draw.text((x + 48, 1068), f"Rows shown: {len(genes)} genes; all values are log2(count+1), row-centered for display.", font=F["small"], fill=rgb("soft"))

    card(draw, (70, 1084, 2060, 96), alpha=188, radius=18)
    draw.text((104, 1118), "Claim boundary", font=F["h2"], fill=rgb("white"))
    draw.text((310, 1122), "This proves local matrix presence and sample alignment; it does not freeze the full source payload or validate benchmark performance.", font=F["body"], fill=rgb("soft"))
    out = CROPS / "organoid_matrix_audit_proof.png"
    canvas.convert("RGB").save(out)
    return out


def draw_fold_lane(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int],
    fold: dict,
    *,
    width: int,
    tone: str,
    compact: bool = False,
) -> None:
    x, y = xy
    train_w = int(width * 0.56)
    test_w = int(width * 0.28)
    draw.text((x, y), str(fold["heldout_value"]), font=F["h2"], fill=rgb("white"))
    draw.text((x, y + 34), str(fold["heldout_factor"]).replace("_", " "), font=F["small"], fill=rgb("soft"))
    y0 = y + 72
    draw.rounded_rectangle((x, y0, x + train_w, y0 + 76), radius=16, fill=rgba("deep3", 235), outline=rgb("blue"), width=2)
    draw.rounded_rectangle((x + train_w + 92, y0, x + train_w + 92 + test_w, y0 + 76), radius=16, fill=rgba("deep3", 235), outline=rgb(tone), width=2)
    draw.text((x + 28, y0 + 22), f"TRAIN {fold['n_train']}", font=F["h2"], fill=rgb("white"))
    draw.text((x + train_w + 120, y0 + 22), f"TEST {fold['n_test']}", font=F["h2"], fill=rgb("white"))
    dash_y = y0 + 38
    for xx in range(x + train_w + 20, x + train_w + 78, 14):
        draw.line((xx, dash_y, xx + 7, dash_y), fill=rgb("red"), width=3)
    draw.text((x + train_w + 26, y0 + 10), "hidden", font=F["tiny"], fill=rgb("red"))
    if not compact:
        train_dist = ", ".join([f"{k} {v}" for k, v in fold["train_label_distribution"].items()])
        test_dist = ", ".join([f"{k} {v}" for k, v in fold["test_label_distribution"].items()])
        draw.text((x + 28, y0 + 104), f"train labels: {train_dist}", font=F["small"], fill=rgb("soft"))
        draw.text((x + train_w + 120, y0 + 104), f"test labels: {test_dist}", font=F["small"], fill=rgb("soft"))


def render_osd120_split_crop() -> Path:
    manifest = load_json("v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json")
    primary = manifest["split"]["primary_candidate_folds"]
    secondary = manifest["split"]["secondary_light_treatment_folds"]
    canvas = Image.new("RGBA", (2200, 1250), rgb("deep") + (255,))
    draw = ImageDraw.Draw(canvas)
    draw.text((70, 54), "OSD-120 interaction split proof", font=F["title"], fill=rgb("white"))
    draw.text((70, 108), "Manifest-derived train/test geometry for genotype-by-light Arabidopsis root task.", font=F["subtitle"], fill=rgb("soft"))
    draw.line((70, 158, 2130, 158), fill=rgba("deep3", 255), width=2)
    badge(draw, (70, 186), "OSD-120", "green")
    badge(draw, (190, 186), "within-source diagnostic", "amber")
    badge(draw, (430, 186), "not leave-one-mission-out", "red")

    card(draw, (70, 270, 1280, 780), alpha=216)
    draw.text((110, 310), "Primary candidate: genotype/ecotype holdout", font=F["h1"], fill=rgb("white"))
    y = 378
    for fold in primary:
        draw_fold_lane(draw, (110, y), fold, width=1040, tone="green")
        y += 204

    card(draw, (1410, 270, 720, 780), alpha=216)
    draw.text((1450, 310), "Secondary light-treatment holdout", font=F["h2"], fill=rgb("white"))
    y = 384
    for fold in secondary:
        draw_fold_lane(draw, (1450, y), fold, width=580, tone="amber", compact=True)
        y += 232

    card(draw, (70, 1084, 2060, 96), alpha=188, radius=18)
    draw.text((104, 1118), "Claim boundary", font=F["h2"], fill=rgb("white"))
    draw.text((310, 1122), "This is a within-source interaction diagnostic; it should not be presented as mission-held-out or cross-species generalization.", font=F["body"], fill=rgb("soft"))
    out = CROPS / "osd120_interaction_split_proof.png"
    canvas.convert("RGB").save(out)
    return out


def render_organoid_split_crop() -> Path:
    manifest = load_json("v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json")
    folds = manifest["split"]["candidate_folds"][:4]
    canvas = Image.new("RGBA", (2200, 1250), rgb("deep") + (255,))
    draw = ImageDraw.Draw(canvas)
    draw.text((70, 54), "Organoid blocked split proof", font=F["title"], fill=rgb("white"))
    draw.text((70, 108), "Manifest-derived folds for organoid type and microglia condition.", font=F["subtitle"], fill=rgb("soft"))
    draw.line((70, 158, 2130, 158), fill=rgba("deep3", 255), width=2)
    badge(draw, (70, 186), "OSD-863/871", "sky")
    badge(draw, (220, 186), "single ISS mission pilot", "amber")
    badge(draw, (460, 186), "not mission-held-out", "red")

    for idx, fold in enumerate(folds):
        x = 70 + (idx % 2) * 1060
        y = 290 + (idx // 2) * 360
        card(draw, (x, y, 990, 300), alpha=216)
        draw_fold_lane(draw, (x + 34, y + 36), fold, width=790, tone="teal", compact=False)

    card(draw, (70, 1084, 2060, 96), alpha=188, radius=18)
    draw.text((104, 1118), "Claim boundary", font=F["h2"], fill=rgb("white"))
    draw.text((310, 1122), "Folds are sample-count-backed and useful for method explanation, but single-mission organoid data cannot support LOMO claims.", font=F["body"], fill=rgb("soft"))
    out = CROPS / "organoid_blocked_split_proof.png"
    canvas.convert("RGB").save(out)
    return out


def paste_fit(canvas: Image.Image, path: Path, box: tuple[int, int, int, int], *, opacity: float = 1.0) -> None:
    image = Image.open(path).convert("RGBA")
    x, y, w, h = box
    scale = min(w / image.width, h / image.height)
    image = image.resize((max(1, int(image.width * scale)), max(1, int(image.height * scale))), Image.Resampling.LANCZOS)
    if opacity < 1.0:
        alpha = image.getchannel("A")
        alpha = alpha.point(lambda value: int(value * opacity))
        image.putalpha(alpha)
    canvas.alpha_composite(image, (x + (w - image.width) // 2, y + (h - image.height) // 2))


def render_production_mock(title: str, subtitle: str, crop_path: Path, out_path: Path, *, lane: str, caveat: str) -> Path:
    W, H = 3840, 2160
    canvas = Image.new("RGBA", (W, H), rgb("deep") + (255,))
    draw = ImageDraw.Draw(canvas)
    for x in [420, 1920, 3320]:
        draw.line((x, 150, x, H - 150), fill=rgba("deep3", 125), width=2)
    for y in [340, 1140, 1840]:
        draw.line((150, y, W - 150, y), fill=rgba("deep3", 150), width=2)
    draw.text((160, 112), lane.upper(), font=F["eyebrow"], fill=rgb("sky"))
    draw.text((160, 164), title, font=font(68, bold=True), fill=rgb("white"))
    draw.text((160, 252), subtitle, font=F["subtitle"], fill=rgb("soft"))
    draw.line((160, 318, W - 160, 318), fill=rgba("deep3", 255), width=2)

    card(draw, (180, 410, 2340, 1340), alpha=188, radius=34)
    paste_fit(canvas, crop_path, (245, 470, 2210, 1200))

    card(draw, (2660, 410, 1020, 490), alpha=220, radius=28)
    draw.text((2715, 465), "How to use this proof", font=F["h1"], fill=rgb("white"))
    usage = [
        "Use the crop as a Z2 evidence surface.",
        "Keep accession/source and status visible.",
        "Add interpretation as editable text only.",
        "Do not promote draft diagnostics to validated claims.",
    ]
    y = 535
    for item in usage:
        draw.ellipse((2718, y + 8, 2734, y + 24), fill=rgb("green"))
        write_wrapped(draw, (2755, y), item, F["body"], rgb("soft"), width=810, line_height=28, max_lines=2)
        y += 74

    card(draw, (2660, 970, 1020, 480), alpha=220, radius=28)
    draw.text((2715, 1025), "Status layer", font=F["h1"], fill=rgb("white"))
    badge(draw, (2715, 1088), "source-derived local proof", "green")
    badge(draw, (2715, 1144), "draft boundary visible", "amber")
    badge(draw, (2715, 1200), "no overclaim", "red")
    write_wrapped(draw, (2715, 1280), caveat, F["body"], rgb("soft"), width=830, line_height=30, max_lines=4)

    card(draw, (2660, 1520, 1020, 230), alpha=190, radius=28)
    draw.text((2715, 1575), "Next replacement", font=F["h1"], fill=rgb("white"))
    write_wrapped(draw, (2715, 1630), "Add official OSDR source-record screenshot as a companion proof object, then run small-screen QA.", F["body"], rgb("soft"), width=830, line_height=30, max_lines=3)

    draw.text((160, 1960), "Proof-crop sprint output. Not a final manuscript figure.", font=F["body"], fill=rgb("muted"))
    canvas.convert("RGB").save(out_path)
    return out_path


def render_contact_sheet(paths: list[Path]) -> Path:
    canvas = Image.new("RGBA", (3000, 1900), rgb("deep") + (255,))
    draw = ImageDraw.Draw(canvas)
    draw.text((80, 58), "Real proof crops for v0.4 replacement", font=F["title"], fill=rgb("white"))
    draw.text((80, 112), "Local evidence objects generated from matrix files and task manifests.", font=F["subtitle"], fill=rgb("soft"))
    labels = [
        "Organoid matrix audit proof",
        "OSD-120 interaction split proof",
        "Organoid blocked split proof",
    ]
    positions = [(80, 220), (1540, 220), (80, 1080)]
    for path, label, (x, y) in zip(paths, labels, positions):
        card(draw, (x, y, 1360, 760), alpha=218)
        paste_fit(canvas, path, (x + 34, y + 34, 1292, 628))
        draw.text((x + 34, y + 680), label, font=F["h2"], fill=rgb("white"))
        draw.text((x + 34, y + 716), rel(path), font=F["small"], fill=rgb("muted"))
    card(draw, (1540, 1080, 1360, 760), alpha=170)
    draw.text((1585, 1132), "Visual QA verdict", font=F["h1"], fill=rgb("white"))
    notes = [
        "Matrix proof now contains real normalized count values, not placeholder heatmap blocks.",
        "Split proof is manifest-derived and visually separates train/test boundaries.",
        "OSD-918 remains source-record/blocker only until local h5ad exists.",
        "Use one proof crop per main slide; combine several only in appendix boards.",
    ]
    y = 1208
    for note in notes:
        draw.ellipse((1588, y + 8, 1606, y + 26), fill=rgb("green"))
        write_wrapped(draw, (1628, y), note, F["body"], rgb("soft"), width=1160, line_height=34, max_lines=3)
        y += 132
    out = MOCKS / "03_proof_crop_contact_sheet.png"
    canvas.convert("RGB").save(out)
    return out


def write_manifest(crop_paths: list[Path], mock_paths: list[Path]) -> None:
    manifest = {
        "created": CREATED,
        "scope": "Real source-object proof crops for v0.4 placeholder replacement.",
        "source_inputs": {
            "organoid_task_manifest": "v9/human_organoid/task_manifests/draft_human_organoid_spaceflight.json",
            "organoid_matrices": [
                "v9/human_organoid/matrices/GLDS-716_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
                "v9/human_organoid/matrices/GLDS-720_rna_seq_Normalized_Counts_GLbulkRNAseq.csv",
            ],
            "osd120_task_manifest": "v9/multispecies/interaction_task_manifests/draft_osd120_arabidopsis_root_light_interaction_spaceflight.json",
        },
        "outputs": {
            "proof_crops": [rel(path) for path in crop_paths],
            "production_mocks": [rel(path) for path in mock_paths],
        },
        "claim_boundary": "Generated crops prove local evidence presence and manifest-derived split geometry; they are not final validated manuscript result figures.",
    }
    (OUT / "source_object_proof_crops_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    qa = {
        "created": CREATED,
        "automatic_checks": {
            "crop_count": len(crop_paths),
            "mock_count": len(mock_paths),
            "all_crops_exist": all(path.exists() for path in crop_paths),
            "all_mocks_exist": all(path.exists() for path in mock_paths),
            "manifest_exists": (OUT / "source_object_proof_crops_manifest.json").exists(),
        },
        "manual_review": {
            "initial_verdict": "ready for production-slide replacement test",
            "strongest_main_slide_candidate": "osd120_interaction_split_proof",
            "strongest_appendix_candidate": "organoid_matrix_audit_proof",
            "caveat": "OSDR webpage screenshots still need separate capture before final source-record slides.",
        },
    }
    (OUT / "source_object_proof_crops_qa.json").write_text(json.dumps(qa, indent=2) + "\n")


def main() -> None:
    ensure()
    matrix = render_organoid_matrix_crop()
    osd120 = render_osd120_split_crop()
    organoid_split = render_organoid_split_crop()
    mock1 = render_production_mock(
        "Replace the schematic matrix with real local evidence",
        "Human organoid extension: local normalized-count matrices and alignment metadata",
        matrix,
        MOCKS / "01_dark_organoid_matrix_proof_replacement.png",
        lane="human organoid proof",
        caveat="Matrix presence and sample alignment are source-object evidence; payload freeze and benchmark claims remain draft.",
    )
    mock2 = render_production_mock(
        "Make the held-out split visible before discussing scores",
        "OSD-120 Arabidopsis root task: genotype-by-light interaction split geometry",
        osd120,
        MOCKS / "02_dark_osd120_split_proof_replacement.png",
        lane="multispecies method proof",
        caveat="Within-source interaction diagnostic only; do not describe this as mission-held-out or cross-species generalization.",
    )
    sheet = render_contact_sheet([matrix, osd120, organoid_split])
    crop_paths = [matrix, osd120, organoid_split]
    mock_paths = [mock1, mock2, sheet]
    write_manifest(crop_paths, mock_paths)
    print(json.dumps({"output": rel(OUT), "proof_crops": len(crop_paths), "production_mocks": len(mock_paths)}, indent=2))


if __name__ == "__main__":
    main()

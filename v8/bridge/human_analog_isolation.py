"""Tier 1A — Human isolation-analog transcriptome integration.

Ingests public human-isolation GEO series (Mars500 520-day confinement, Antarctic
Concordia winter-over, HERA/HI-SEAS), fits a 2-way factorial per gene along the
axis available in each dataset, and compares the resulting stressor coefficients
to v8 rodent ICP-validated β's via per-pathway Spearman.

Framework-first: dataset definitions live in CANDIDATES; a GEO series that fails
to download (restricted, withdrawn, or lacking per-sample design metadata) is
skipped with a transparent error rather than crashing the whole run.

Outputs:
  v8/bridge/evaluation/human_analog_betas_{accession}.csv
  v8/bridge/evaluation/rodent_to_human_analog_transfer.json
"""
from __future__ import annotations

import json
import re
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from scipy import stats

warnings.filterwarnings("ignore")

try:
    import GEOparse
except ImportError:
    print("ERROR: GEOparse not installed. Run: pip install GEOparse --break-system-packages")
    sys.exit(1)

REPO_ROOT = Path(__file__).resolve().parents[2]
BRIDGE_DIR = Path(__file__).resolve().parent
EVAL_DIR = BRIDGE_DIR / "evaluation"
CACHE_DIR = BRIDGE_DIR / "geo_cache"
EVAL_DIR.mkdir(parents=True, exist_ok=True)
CACHE_DIR.mkdir(parents=True, exist_ok=True)


# -------- dataset registry --------

@dataclass
class AnalogStudy:
    accession: str            # GEO series accession
    label: str                # short descriptor
    axes: list[str]           # stressor axes to model (e.g. ["isolation", "time"])
    design_fn: Callable[[pd.DataFrame], pd.DataFrame]  # maps GSM metadata -> design columns
    platform_hint: str = ""   # "illumina_array", "affy", "rnaseq_counts"
    notes: str = ""


def _astronaut_iss_design(meta: pd.DataFrame) -> pd.DataFrame:
    """GSE74708 ISS astronauts (hair follicle microarray).

    10 astronauts, 6 timepoints including 2 in-flight. Design axes:
      - inflight (0/1): currently on ISS when sampled
      - time_num: pre-flight (−), inflight (+small), post-flight (+large)
    Extracted from title strings like 'Preflight-180', 'Inflight-90', 'Postflight-0'.
    """
    df = meta.copy()
    title = df["title"].astype(str)
    blob = (title + " | " + df.get("characteristics_ch1", pd.Series("", index=df.index)).astype(str)).str.lower()

    inflight = np.where(blob.str.contains("inflight|in-flight|on-orbit|in orbit"), 1, 0)
    postflight = np.where(blob.str.contains("postflight|post-flight|post flight|r\\+"), 1, 0)
    preflight = np.where(blob.str.contains(r"preflight|pre-flight|pre flight|l-"), 1, 0)

    # Extract numeric days (signed)
    day = np.zeros(len(df))
    for i, t in enumerate(blob):
        m = re.search(r"([lr+\-])[\s-]?(\d+)", t)
        if m:
            sign = -1 if m.group(1) in ("l", "-") else 1
            day[i] = sign * int(m.group(2))

    df["inflight"] = inflight
    df["postflight"] = postflight
    df["day_signed"] = day
    df["inflight_x_day"] = inflight * np.abs(day)
    return df[["inflight", "postflight", "day_signed"]]


def _concordia_design(meta: pd.DataFrame) -> pd.DataFrame:
    """GSE149321 Antarctic Concordia (3233 m, 12-month winter-over, n=13).

    Timepoints: PRE, Day30, Day150, POST. Axes:
      - isolation: 1 during mission (Day30, Day150), 0 otherwise (PRE, POST)
      - time_num: 0=PRE, 1=Day30, 2=Day150, 3=POST
      - iso_x_time: interaction
    """
    df = meta.copy()
    title = df["title"].astype(str).str.lower()
    blob = (title + " | " + df.get("characteristics_ch1", pd.Series("", index=df.index)).astype(str).str.lower())

    t = np.full(len(df), -1)
    t = np.where(blob.str.contains(r"\bpre\b|t0|before"), 0, t)
    t = np.where(blob.str.contains(r"day\s*30|\bd30\b|t1"), 1, t)
    t = np.where(blob.str.contains(r"day\s*150|\bd150\b|t2"), 2, t)
    t = np.where(blob.str.contains(r"\bpost\b|t3|after"), 3, t)

    iso = np.where(np.isin(t, [1, 2]), 1, 0)
    df["isolation"] = iso
    df["time_num"] = t
    df["iso_x_time"] = iso * t
    df = df[df["time_num"] >= 0]  # drop unmatched
    return df[["isolation", "time_num", "iso_x_time"]]


CANDIDATES: list[AnalogStudy] = [
    AnalogStudy("GSE74708", "ISS_astronaut_hair", ["inflight", "day_signed"],
                _astronaut_iss_design, "affy_array",
                notes="Terada 2016: 10 astronauts x 6 timepoints, hair follicle (PMID 27029003)"),
    AnalogStudy("GSE149321", "Concordia_12mo", ["isolation", "time_num"],
                _concordia_design, "agilent_array",
                notes="Pisa 2020: n=13 Concordia winter-over, PBMC, 4 tps (PMID 32518224)"),
]


# -------- GEO ingest --------

def _download_series_matrix(accession: str) -> Path | None:
    """Download GEO series_matrix.txt.gz directly via HTTPS (bypasses GEOparse bugs)."""
    dest = CACHE_DIR / f"{accession}_series_matrix.txt.gz"
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    prefix = accession[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/matrix/{accession}_series_matrix.txt.gz"
    import subprocess
    r = subprocess.run(["curl", "-sL", "-o", str(dest), url], capture_output=True)
    if r.returncode != 0 or not dest.exists() or dest.stat().st_size < 1000:
        return None
    return dest


def _download_platform_annotation(gpl_id: str) -> Path | None:
    """Download GPL annotation SOFT table for probe->symbol mapping."""
    dest = CACHE_DIR / f"{gpl_id}.annot.gz"
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    prefix = gpl_id[:-3] + "nnn" if len(gpl_id) > 5 else gpl_id + "nnn"
    # Try annot file first (has symbol); fall back to data table from SOFT
    urls = [
        f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{prefix}/{gpl_id}/annot/{gpl_id}.annot.gz",
        f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{prefix}/{gpl_id}/soft/{gpl_id}_family.soft.gz",
    ]
    import subprocess
    for url in urls:
        r = subprocess.run(["curl", "-sfL", "-o", str(dest), url], capture_output=True)
        if dest.exists() and dest.stat().st_size > 1000:
            return dest
    return None


def _parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """Parse a GEO series_matrix.txt.gz; returns (expr[probe x sample], meta[sample x fields], gpl_id)."""
    import gzip
    meta_rows: dict[str, list[str]] = {}
    sample_ids: list[str] = []
    expr_lines: list[str] = []
    gpl_id = ""
    in_table = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!Sample_geo_accession"):
                sample_ids = line.split("\t")[1:]
                sample_ids = [s.strip('"') for s in sample_ids]
            elif line.startswith("!Sample_title"):
                meta_rows["title"] = [s.strip('"') for s in line.split("\t")[1:]]
            elif line.startswith("!Sample_source_name_ch1"):
                meta_rows["source"] = [s.strip('"') for s in line.split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                fld = [s.strip('"') for s in line.split("\t")[1:]]
                key = f"char_{len([k for k in meta_rows if k.startswith('char_')])}"
                meta_rows[key] = fld
            elif line.startswith("!Series_platform_id"):
                gpl_id = line.split("\t")[1].strip('"')
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            elif line.startswith("!series_matrix_table_end"):
                in_table = False
            elif in_table:
                expr_lines.append(line)

    if not expr_lines:
        raise RuntimeError("no expression table in series matrix")

    from io import StringIO
    expr = pd.read_csv(StringIO("\n".join(expr_lines)), sep="\t", index_col=0, na_values=["", "NA", "null"])
    expr.index = expr.index.astype(str).str.strip('"')
    expr.columns = [c.strip('"') for c in expr.columns]

    meta = pd.DataFrame(meta_rows, index=sample_ids)
    # Concatenate characteristic cols into a single string for regex matching
    char_cols = [c for c in meta.columns if c.startswith("char_")]
    meta["characteristics_ch1"] = meta[char_cols].fillna("").agg(" | ".join, axis=1) if char_cols else ""

    return expr, meta, gpl_id


def _parse_platform_annotation(path: Path) -> dict[str, str] | None:
    """Parse GPL annot/soft to produce probe_id -> gene_symbol dict."""
    import gzip
    probe2sym: dict[str, str] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        # SOFT format: look for !platform_table_begin
        in_table = False
        header: list[str] = []
        sym_idx = None
        id_idx = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!platform_table_begin") or line.startswith("!annotation_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end") or line.startswith("!annotation_table_end"):
                break
            if in_table:
                parts = line.split("\t")
                if not header:
                    header = [p.strip('"') for p in parts]
                    for i, h in enumerate(header):
                        if h.lower() in ("id", "probe_id", "probeset_id"):
                            id_idx = i
                        if h in ("Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol", "GeneSymbol", "ILMN_Gene"):
                            sym_idx = i
                    if sym_idx is None:
                        for i, h in enumerate(header):
                            if "gene" in h.lower() and "symbol" in h.lower():
                                sym_idx = i
                                break
                    if id_idx is None:
                        id_idx = 0
                    continue
                if sym_idx is None or len(parts) <= max(id_idx, sym_idx):
                    continue
                pid = parts[id_idx].strip('"')
                sym = parts[sym_idx].strip('"').split("///")[0].split(" // ")[0].strip()
                if sym and sym.upper() not in ("", "NA", "---"):
                    probe2sym[pid] = sym.upper()
    return probe2sym if probe2sym else None


def load_geo(accession: str) -> tuple[pd.DataFrame, pd.DataFrame] | None:
    """Return (expression_matrix[gene x sample], sample_metadata) or None on failure."""
    mat_path = _download_series_matrix(accession)
    if mat_path is None:
        print(f"  [{accession}] series matrix download failed")
        return None
    try:
        expr, meta, gpl_id = _parse_series_matrix(mat_path)
    except Exception as e:
        print(f"  [{accession}] parse failed: {e}")
        return None
    print(f"  expr shape: {expr.shape}, GPL: {gpl_id}")

    # Map probes -> gene symbols
    if gpl_id:
        ann_path = _download_platform_annotation(gpl_id)
        if ann_path is not None:
            probe2sym = _parse_platform_annotation(ann_path)
            if probe2sym:
                expr_num = expr.apply(pd.to_numeric, errors="coerce")
                expr_num["_sym"] = expr_num.index.map(probe2sym.get)
                expr_num = expr_num[expr_num["_sym"].notna()]
                before = len(expr_num)
                expr_num = expr_num.groupby("_sym").mean(numeric_only=True)
                print(f"  mapped {before} probes -> {len(expr_num)} unique gene symbols (GPL: {gpl_id})")
                return expr_num, meta
            else:
                print(f"  [{accession}] platform annotation parse empty; keeping probe IDs")
        else:
            print(f"  [{accession}] platform annotation download failed for {gpl_id}")

    expr_num = expr.apply(pd.to_numeric, errors="coerce").dropna(how="all")
    return expr_num, meta


# -------- factorial fit --------

def fit_factorial_human(expr: pd.DataFrame, design: pd.DataFrame) -> pd.DataFrame:
    """OLS per gene: expr ~ 1 + axis1 + axis2 + axis1:axis2 (2-way).

    Returns DataFrame with beta, SE, t, p per design column.
    """
    samples = [s for s in expr.columns if s in design.index]
    if len(samples) < 6:
        raise ValueError(f"too few matched samples: {len(samples)}")
    X = design.loc[samples].astype(float).values
    X = np.column_stack([np.ones(len(samples)), X])
    Y = expr[samples].astype(float).values  # genes x samples
    # log-transform if appears raw
    if np.nanmax(Y) > 50:
        Y = np.log2(np.clip(Y, 1e-3, None) + 1)

    n, p = X.shape
    results = []
    # Least-squares in bulk
    XtX_inv = np.linalg.pinv(X.T @ X)
    B = XtX_inv @ X.T @ Y.T  # (p, n_genes)
    resid = Y.T - X @ B
    sigma2 = (resid ** 2).sum(axis=0) / max(n - p, 1)
    se = np.sqrt(np.outer(np.diag(XtX_inv), sigma2))  # (p, n_genes)
    t = np.divide(B, se, out=np.zeros_like(B), where=se > 0)
    pval = 2 * (1 - stats.t.cdf(np.abs(t), df=max(n - p, 1)))

    cols = ["intercept"] + list(design.columns)
    out = {"gene": expr.index.values}
    for i, c in enumerate(cols):
        out[f"beta_{c}"] = B[i]
        out[f"se_{c}"] = se[i]
        out[f"t_{c}"] = t[i]
        out[f"p_{c}"] = pval[i]
    return pd.DataFrame(out)


# -------- cross-species comparison --------

def load_rodent_betas() -> dict[str, pd.DataFrame]:
    d = REPO_ROOT / "v8/decompose/evaluation"
    out = {}
    for name in ["spleen", "skin_analog", "brain"]:
        p = d / f"factorial_{name}_betas.csv"
        if p.exists():
            df = pd.read_csv(p)
            df["gene"] = df["gene"].astype(str).str.upper()
            out[name] = df
    return out


def compare_to_rodent(human: pd.DataFrame, rodents: dict[str, pd.DataFrame], human_axis: str,
                      rodent_axis: str) -> dict:
    """Spearman correlation of per-gene effect sizes between human & each rodent analog."""
    h_col = f"beta_{human_axis}"
    r_col = f"beta_{rodent_axis}"
    if h_col not in human.columns:
        return {"error": f"{h_col} not in human betas"}
    human = human[["gene", h_col]].dropna()
    human["gene"] = human["gene"].astype(str).str.upper()
    result = {}
    for name, rdf in rodents.items():
        if r_col not in rdf.columns:
            result[name] = {"error": f"{r_col} missing"}
            continue
        merged = human.merge(rdf[["gene", r_col]], on="gene").dropna()
        if len(merged) < 50:
            result[name] = {"n_overlap": len(merged), "note": "too few"}
            continue
        rho, p = stats.spearmanr(merged[h_col], merged[r_col])
        result[name] = {"n_overlap": int(len(merged)), "spearman_r": float(rho), "p": float(p)}
    return result


# -------- main --------

def main() -> None:
    rodents = load_rodent_betas()
    print(f"loaded rodent factorial βs: {list(rodents)}")

    summary = {"rodent_analogs": list(rodents), "human_studies": {}}

    for study in CANDIDATES:
        print(f"\n=== {study.accession} ({study.label}) ===")
        loaded = load_geo(study.accession)
        if loaded is None:
            summary["human_studies"][study.accession] = {"status": "fetch_failed"}
            continue
        expr, meta = loaded
        print(f"  expr shape: {expr.shape}, meta rows: {len(meta)}")

        try:
            design = study.design_fn(meta)
        except Exception as e:
            print(f"  design extraction failed: {e}")
            summary["human_studies"][study.accession] = {"status": "design_failed", "err": str(e)}
            continue

        # Drop constant columns (no contrast -> cannot estimate)
        nunique = design.nunique()
        design = design.loc[:, nunique > 1]
        if design.shape[1] == 0:
            summary["human_studies"][study.accession] = {
                "status": "no_variable_axes",
                "nunique_by_axis": nunique.to_dict(),
            }
            print(f"  no variable design axes (all constant) — skip")
            continue

        try:
            betas = fit_factorial_human(expr, design)
        except Exception as e:
            summary["human_studies"][study.accession] = {"status": "fit_failed", "err": str(e)}
            print(f"  fit failed: {e}")
            continue

        out_p = EVAL_DIR / f"human_analog_betas_{study.accession}.csv"
        betas.to_csv(out_p, index=False)
        print(f"  wrote {out_p} ({len(betas)} genes, axes: {list(design.columns)})")

        # Cross-species: compare first human axis to rodent IR (radiation) + HLU (unloading)
        per_axis = {}
        for h_ax in design.columns:
            per_axis[h_ax] = {
                "vs_rodent_IR": compare_to_rodent(betas, rodents, h_ax, "IR"),
                "vs_rodent_HLU": compare_to_rodent(betas, rodents, h_ax, "HLU"),
                "vs_rodent_HLUxIR": compare_to_rodent(betas, rodents, h_ax, "HLUxIR"),
            }

        summary["human_studies"][study.accession] = {
            "status": "ok",
            "label": study.label,
            "n_genes": int(len(betas)),
            "n_samples": int(expr.shape[1]),
            "axes": list(design.columns),
            "cross_species": per_axis,
        }

    summary_p = EVAL_DIR / "rodent_to_human_analog_transfer.json"
    summary_p.write_text(json.dumps(summary, indent=2))
    print(f"\n✓ wrote {summary_p}")

    # Human-readable top-line
    print("\n=== TOP-LINE ===")
    for acc, info in summary["human_studies"].items():
        if info.get("status") != "ok":
            print(f"  {acc}: {info.get('status')}")
            continue
        best = None
        for h_ax, comps in info["cross_species"].items():
            for r_ax, d in comps.items():
                if isinstance(d, dict) and "spearman_r" in d:
                    for tissue, res in d.items() if isinstance(d, dict) else []:
                        pass  # nested structure
        # flatten search for best rho
        flat = []
        for h_ax, comps in info["cross_species"].items():
            for contrast, tissue_dict in comps.items():
                for tissue, res in tissue_dict.items():
                    if isinstance(res, dict) and "spearman_r" in res:
                        flat.append((abs(res["spearman_r"]), h_ax, contrast, tissue, res))
        flat.sort(reverse=True)
        if flat:
            _, h_ax, contrast, tissue, res = flat[0]
            print(f"  {acc} [{info['label']}]: best ρ={res['spearman_r']:+.3f} "
                  f"({h_ax} vs {contrast}/{tissue}, n={res['n_overlap']}, p={res['p']:.2e})")


if __name__ == "__main__":
    main()

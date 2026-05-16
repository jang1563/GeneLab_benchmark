"""v8 Pillar 1 (BRIDGE) — per-tissue mouse NES vs human compartments.

Refines link_spaceomicsbench.py by building a per-tissue mouse NES profile
(aggregated across missions) and computing the full
mouse-tissue × human-compartment Spearman matrix. Hypothesis: thymus (strongest
v1 transfer) ↔ I4 lymphoid subsets (CD4_T, CD8_T, B, NK) should dominate.

Outputs
-------
  v8/bridge/evaluation/tissue_nes_spearman.csv
  v8/bridge/evaluation/tissue_nes_spearman.json
"""
from __future__ import annotations

import glob
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
FGSEA_ROOT = REPO_ROOT / "processed" / "fgsea"
SOB_ROOT = Path(os.environ.get("SPACEOMICS_ROOT", REPO_ROOT.parent / "SpaceOmicsBench")).expanduser()
I4_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_i4_pbmc.csv"
TWINS_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_NASA_Twins.csv"
OUT_DIR = Path(__file__).resolve().parent / "evaluation"

TISSUES = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]
DBS = ["hallmark", "kegg", "reactome", "mitocarta", "c2cgp", "c5bp"]
MISSION_ALIASES = {"TBD": "OSD-397"}


def _pathway_col_upper(df: pd.DataFrame) -> pd.DataFrame:
    """fgsea writes pathway IDs already in uppercase but ensure consistency."""
    df = df.copy()
    df["pathway"] = df["pathway"].astype(str).str.strip()
    return df


def load_mouse_tissue_nes() -> pd.DataFrame:
    """Return long df: tissue, pathway, mouse_NES, n_missions."""
    frames = []
    for tissue in TISSUES:
        for db in DBS:
            for f in sorted(glob.glob(str(FGSEA_ROOT / tissue / f"*_fgsea_{db}.csv"))):
                base = os.path.basename(f).replace(".csv", "")
                m = re.match(r"(.+)_fgsea_(\w+)", base)
                if not m:
                    continue
                mission = MISSION_ALIASES.get(m.group(1), m.group(1))
                try:
                    df = pd.read_csv(f, usecols=["pathway", "NES", "padj"])
                except Exception:
                    continue
                df = _pathway_col_upper(df)
                df["tissue"] = tissue
                df["mission"] = mission
                df["db"] = db
                frames.append(df)
    long = pd.concat(frames, ignore_index=True)
    long["NES"] = pd.to_numeric(long["NES"], errors="coerce")
    agg = (
        long.dropna(subset=["NES"])
        .groupby(["tissue", "pathway"])
        .agg(mouse_NES=("NES", "mean"), n_missions=("mission", "nunique"))
        .reset_index()
    )
    return agg


def load_human() -> pd.DataFrame:
    frames = []

    i4 = pd.read_csv(I4_CSV)
    i4 = _pathway_col_upper(i4)
    i4["mission"] = "Inspiration4"
    i4["compartment"] = i4["celltype"].astype(str) + "|" + i4["timepoint"].astype(str)
    frames.append(i4[["mission", "compartment", "pathway", "NES"]])

    tw = pd.read_csv(TWINS_CSV)
    tw = tw.rename(columns={"Pathway": "pathway"})
    tw = _pathway_col_upper(tw)
    tw["mission"] = "NASA_Twins"
    tw["compartment"] = tw["CellType"].astype(str) + "|" + tw["timepoint"].astype(str)
    frames.append(tw[["mission", "compartment", "pathway", "NES"]])

    out = pd.concat(frames, ignore_index=True)
    out["NES"] = pd.to_numeric(out["NES"], errors="coerce")
    return out.dropna(subset=["NES"])


def _concordance(x: np.ndarray, y: np.ndarray, t: float = 1.0) -> float:
    mask = (np.abs(x) >= t) & (np.abs(y) >= t)
    return float(np.mean(np.sign(x[mask]) == np.sign(y[mask]))) if mask.sum() else float("nan")


def compute_matrix(mouse: pd.DataFrame, human: pd.DataFrame, min_pathways: int = 20) -> pd.DataFrame:
    rows = []
    mouse_by_tissue = {t: sub[["pathway", "mouse_NES"]] for t, sub in mouse.groupby("tissue")}
    for (mission, compartment), human_sub in human.groupby(["mission", "compartment"]):
        human_sub = human_sub.drop_duplicates("pathway")
        for tissue, mouse_sub in mouse_by_tissue.items():
            merged = mouse_sub.merge(human_sub[["pathway", "NES"]], on="pathway", how="inner")
            n = len(merged)
            if n < min_pathways:
                rows.append({"mouse_tissue": tissue, "human_mission": mission,
                             "human_compartment": compartment, "n_pathways": n,
                             "spearman_r": np.nan, "spearman_p": np.nan,
                             "concordance_1.0": np.nan})
                continue
            rho, p = spearmanr(merged["mouse_NES"], merged["NES"])
            rows.append({
                "mouse_tissue": tissue,
                "human_mission": mission,
                "human_compartment": compartment,
                "n_pathways": int(n),
                "spearman_r": float(rho),
                "spearman_p": float(p),
                "concordance_1.0": _concordance(merged["mouse_NES"].values, merged["NES"].values, 1.0),
            })
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    mouse = load_mouse_tissue_nes()
    human = load_human()
    print(f"mouse: {mouse.shape[0]} tissue×pathway rows, {mouse.tissue.nunique()} tissues, "
          f"{mouse.pathway.nunique()} unique pathways")
    print(f"human: {human.shape[0]} rows, {human.groupby('mission').compartment.nunique().to_dict()}")

    mat = compute_matrix(mouse, human)
    csv_path = OUT_DIR / "tissue_nes_spearman.csv"
    mat.to_csv(csv_path, index=False)

    i4 = mat[mat.human_mission == "Inspiration4"]
    pivot = i4.pivot_table(index="mouse_tissue", columns="human_compartment", values="spearman_r")
    pivot = pivot.reindex(TISSUES)

    summary = {
        "n_pairs": int(len(mat)),
        "mean_r_by_mouse_tissue": {t: float(v) for t, v in mat.groupby("mouse_tissue").spearman_r.mean().items()},
        "i4_best_human_compartment_per_tissue": {
            t: (i4[(i4.mouse_tissue == t) & i4.spearman_r.notna()]
                .nlargest(1, "spearman_r")[["human_compartment", "spearman_r", "n_pathways"]]
                .to_dict(orient="records")[0] if (i4.mouse_tissue == t).any() else None)
            for t in TISSUES
        },
        "top10_i4_pairs": i4.nlargest(10, "spearman_r")[
            ["mouse_tissue", "human_compartment", "spearman_r", "n_pathways", "concordance_1.0"]
        ].to_dict(orient="records"),
        "i4_pivot": pivot.round(3).fillna("NA").to_dict(),
    }
    with open(OUT_DIR / "tissue_nes_spearman.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== Inspiration4 mouse-tissue × human-compartment Spearman r ===")
    print(pivot.round(3).to_string())
    print("\nMean r by mouse tissue (across all human compartments):")
    for t, v in sorted(summary["mean_r_by_mouse_tissue"].items(), key=lambda kv: -kv[1]):
        print(f"  {t:<14} r={v:.3f}")


if __name__ == "__main__":
    main()

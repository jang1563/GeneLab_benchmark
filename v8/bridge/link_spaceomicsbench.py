"""v8 Pillar 1 (BRIDGE) — mouse → human pathway-NES conservation.

Uses SpaceOmicsBench's pre-computed fGSEA NES tables as the human reference
rather than re-running fGSEA from raw SOMA data. This mirrors v1 Fig 4
(NES-conservation r=0.9 predicts Category B transfer) but applied across species.

Inputs
------
  mouse:  SpaceOmicsBench/v2_public/data/processed/gt_conserved_pathways_GeneLab.csv
  human:  SpaceOmicsBench/v2_public/data/processed/gt_conserved_pathways_i4_pbmc.csv
          SpaceOmicsBench/v2_public/data/processed/gt_conserved_pathways_NASA_Twins.csv

Outputs (v8/bridge/evaluation/)
-------
  species_transfer_nes.json   per-comparison Spearman, direction concordance, n_pathways
  species_transfer_nes.csv    flat table for plotting
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
SOB_ROOT = Path(os.environ.get("SPACEOMICS_ROOT", REPO_ROOT.parent / "SpaceOmicsBench")).expanduser()
MOUSE_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_GeneLab.csv"
I4_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_i4_pbmc.csv"
TWINS_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_NASA_Twins.csv"

OUT_DIR = Path(__file__).resolve().parent / "evaluation"


def load_mouse() -> pd.DataFrame:
    df = pd.read_csv(MOUSE_CSV)
    return df[["pathway", "NES", "padj", "info"]].rename(columns={"NES": "mouse_NES", "padj": "mouse_padj"})


def _concordance(x: pd.Series, y: pd.Series, threshold: float = 1.0) -> float:
    """Sign-match rate on pathways where both |NES| >= threshold. NaN if no pathways pass."""
    mask = (x.abs() >= threshold) & (y.abs() >= threshold)
    if mask.sum() == 0:
        return float("nan")
    return float((np.sign(x[mask]) == np.sign(y[mask])).mean())


def compare_one(mouse: pd.DataFrame, human: pd.DataFrame, human_NES_col: str = "NES") -> dict:
    merged = mouse.merge(
        human[["pathway", human_NES_col]].rename(columns={human_NES_col: "human_NES"}),
        on="pathway",
        how="inner",
    )
    n = len(merged)
    if n < 5:
        return {"n_pathways": n, "spearman_r": float("nan"), "spearman_p": float("nan"),
                "concordance_1.0": float("nan"), "concordance_1.5": float("nan")}
    rho, p = spearmanr(merged["mouse_NES"], merged["human_NES"])
    return {
        "n_pathways": int(n),
        "spearman_r": float(rho),
        "spearman_p": float(p),
        "concordance_1.0": _concordance(merged["mouse_NES"], merged["human_NES"], 1.0),
        "concordance_1.5": _concordance(merged["mouse_NES"], merged["human_NES"], 1.5),
    }


def run_i4_pbmc(mouse: pd.DataFrame) -> pd.DataFrame:
    df = pd.read_csv(I4_CSV)
    rows = []
    for (celltype, timepoint), sub in df.groupby(["celltype", "timepoint"]):
        r = compare_one(mouse, sub)
        r.update({"mission": "Inspiration4", "compartment": f"{celltype}|{timepoint}"})
        rows.append(r)
    return pd.DataFrame(rows)


def run_twins(mouse: pd.DataFrame) -> pd.DataFrame:
    df = pd.read_csv(TWINS_CSV)
    NES_col = "NES"
    rows = []
    group_cols = [c for c in ("CellType", "celltype", "Celltype", "timepoint", "ExperimentType") if c in df.columns]
    if not group_cols:
        return pd.DataFrame()
    for keys, sub in df.groupby(group_cols):
        r = compare_one(mouse.rename(columns={"info": "_info"}), sub.rename(columns={"Pathway": "pathway"}), human_NES_col=NES_col)
        tag = "|".join(str(k) for k in (keys if isinstance(keys, tuple) else (keys,)))
        r.update({"mission": "NASA_Twins", "compartment": tag})
        rows.append(r)
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mouse = load_mouse()
    i4 = run_i4_pbmc(mouse)
    twins = run_twins(mouse)
    all_df = pd.concat([i4, twins], ignore_index=True)
    all_df = all_df.sort_values(["mission", "spearman_r"], ascending=[True, False])

    csv_path = out_dir / "species_transfer_nes.csv"
    json_path = out_dir / "species_transfer_nes.json"
    all_df.to_csv(csv_path, index=False)

    summary = {
        "n_mouse_pathways": int(len(mouse)),
        "n_comparisons": int(len(all_df)),
        "by_mission": {
            m: {
                "n": int(sub.shape[0]),
                "mean_spearman_r": float(sub["spearman_r"].mean()),
                "max_spearman_r": float(sub["spearman_r"].max()) if sub["spearman_r"].notna().any() else float("nan"),
                "mean_concordance_1.0": float(sub["concordance_1.0"].mean()),
            }
            for m, sub in all_df.groupby("mission")
        },
        "top5_compartments": all_df.nlargest(5, "spearman_r")[["mission", "compartment", "spearman_r", "n_pathways", "concordance_1.0"]].to_dict(orient="records"),
    }
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"Wrote {csv_path}")
    print(f"Wrote {json_path}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

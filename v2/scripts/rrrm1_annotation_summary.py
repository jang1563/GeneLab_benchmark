#!/usr/bin/env python3
"""
rrrm1_annotation_summary.py
Create cross-tissue summary tables and a stacked bar plot from RRRM-1 broad annotations.
"""

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

BASE = Path(os.environ.get("SCRATCH_DIR", ".")) / "rrrm1_scrna" / "downstream_initial" / "annotations"
TABLEDIR = BASE / "tables"
FIGDIR = BASE / "figures"
SUMMARYDIR = BASE / "summary"


def infer_tissue(path: Path) -> str:
    name = path.name
    for tissue in ["blood", "eye", "muscle", "skin"]:
        if f"_{tissue}_" in name:
            return tissue
    raise ValueError(f"Could not infer tissue from {path}")


def main() -> None:
    SUMMARYDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)

    rows = []
    for path in sorted(TABLEDIR.glob("RRRM1_*_broad_celltype_counts.csv")):
        tissue = infer_tissue(path)
        df = pd.read_csv(path)
        df["tissue"] = tissue
        rows.append(df)

    counts = pd.concat(rows, ignore_index=True)
    counts = counts[["tissue", "broad_celltype", "n_cells"]].sort_values(
        ["tissue", "n_cells"],
        ascending=[True, False],
    )
    counts.to_csv(SUMMARYDIR / "RRRM1_broad_celltype_counts_all_tissues.csv", index=False)

    totals = counts.groupby("tissue", as_index=False)["n_cells"].sum().rename(
        columns={"n_cells": "total_cells"}
    )
    prop = counts.merge(totals, on="tissue", how="left")
    prop["fraction"] = prop["n_cells"] / prop["total_cells"]
    prop["percent"] = prop["fraction"] * 100
    prop.to_csv(SUMMARYDIR / "RRRM1_broad_celltype_proportions_all_tissues.csv", index=False)

    plot_df = prop.pivot(index="tissue", columns="broad_celltype", values="fraction").fillna(0)
    plot_df = plot_df.loc[["blood", "eye", "muscle", "skin"]]

    ax = plot_df.plot(
        kind="bar",
        stacked=True,
        figsize=(10, 5),
        width=0.8,
        colormap="tab20",
    )
    ax.set_ylabel("Fraction of cells")
    ax.set_xlabel("")
    ax.set_title("RRRM-1 Broad Cell Type Composition by Tissue")
    ax.legend(
        title="Broad cell type",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
    )
    plt.tight_layout()
    plt.savefig(FIGDIR / "RRRM1_broad_celltype_stacked_bar.png", dpi=220, bbox_inches="tight")
    plt.close()

    print(f"Saved: {SUMMARYDIR / 'RRRM1_broad_celltype_counts_all_tissues.csv'}")
    print(f"Saved: {SUMMARYDIR / 'RRRM1_broad_celltype_proportions_all_tissues.csv'}")
    print(f"Saved: {FIGDIR / 'RRRM1_broad_celltype_stacked_bar.png'}")


if __name__ == "__main__":
    main()

"""v8 Pillar 1 (BRIDGE) — figure for the mouse-tissue x human-compartment
proteostasis matrix.

Reads v8/bridge/evaluation/proteostasis_matrix.csv (produced by
proteostasis_matrix.py) and renders a two-panel supplementary figure:

  Panel A : Spearman r heatmap on the 11 proteostasis pathways
            (mouse tissues on rows, human compartments on columns, columns
            grouped by source: Inspiration4 PBMC | NASA Twins | Inspiration4
            Skin).
  Panel B : Same layout, on the 14 non-proteostasis pathways from the same
            25-pathway overlap, as a within-overlap biology control.

Cells with n_pathways < 3 are masked. Diverging Red-Blue colormap centered
on 0, vmin=-1, vmax=+1. Cells annotated with the Spearman r and a single
'*' marker when permutation_p < 0.05.

Outputs
-------
  v8/figures/SupplementaryFigure_Proteostasis_Conservation.png
  v8/figures/SupplementaryFigure_Proteostasis_Conservation.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Rectangle

REPO_ROOT = Path(__file__).resolve().parents[2]
MATRIX_CSV = REPO_ROOT / "v8" / "bridge" / "evaluation" / "proteostasis_matrix.csv"
FIG_DIR = REPO_ROOT / "v8" / "figures"

MOUSE_TISSUES_ORDER = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]
MISSION_ORDER = ["Inspiration4_PBMC", "NASA_Twins", "Inspiration4_Skin"]
MISSION_LABEL = {
    "Inspiration4_PBMC": "Inspiration4 PBMC",
    "NASA_Twins": "NASA Twins",
    "Inspiration4_Skin": "i4 Skin",
}

# Compartment-label cleanup for the column tick labels.
COMPARTMENT_RENAMES = {
    "|Immediately Post-flight": "",
    "|Post-flight vs Pre-flight": " (post)",
    "|In-flight vs Pre-flight": " (in)",
    ";Multivariate (PolyA+ and Ribodepleted together)": ";Mv",
    ";Ribodepleted": ";Rb",
    ";PolyA+": ";Pa",
}

PANEL_TITLES = {
    "proteostasis_11": "A   Proteostasis-11 pathway set (translation/proteostasis)",
    "nonproteostasis_14": "B   Non-proteostasis-14 pathway set (within-overlap control)",
}


def _short_compartment(label: str) -> str:
    for k, v in COMPARTMENT_RENAMES.items():
        label = label.replace(k, v)
    return label


def _ordered_columns(df: pd.DataFrame) -> list[tuple[str, str]]:
    """Return [(mission, compartment), ...] in mission-grouped order. Within a
    mission, sort by mean Spearman r across mouse tissues (most positive
    first) so the figure reads as a gradient."""
    cols: list[tuple[str, str]] = []
    for mission in MISSION_ORDER:
        sub = df[df.mission == mission]
        if sub.empty:
            continue
        rank = (
            sub.groupby("compartment")["spearman_r"]
            .mean()
            .sort_values(ascending=False)
            .index.tolist()
        )
        for c in rank:
            cols.append((mission, c))
    return cols


def _pivot(df: pd.DataFrame, value_col: str, cols: list[tuple[str, str]]) -> pd.DataFrame:
    """Pivot to a [mouse_tissue x (mission|compartment)] matrix in the order
    given by cols."""
    df = df.copy()
    df["colkey"] = df["mission"] + "::" + df["compartment"]
    pivot = df.pivot(index="mouse_tissue", columns="colkey", values=value_col)
    pivot = pivot.reindex(MOUSE_TISSUES_ORDER)
    pivot = pivot.reindex(columns=[f"{m}::{c}" for m, c in cols])
    return pivot


def _draw_panel(ax: plt.Axes, df: pd.DataFrame, title: str, *,
                show_xlabels: bool, cbar_ax: plt.Axes | None) -> None:
    cols = _ordered_columns(df)
    r_pivot = _pivot(df, "spearman_r", cols)
    p_pivot = _pivot(df, "permutation_p", cols)
    n_pivot = _pivot(df, "n_pathways", cols)

    # Mask cells with n_pathways < 3 or missing r.
    mask = (n_pivot < 3) | r_pivot.isna()
    data = np.ma.array(r_pivot.values, mask=mask.values)

    norm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad(color="#f0f0f0")
    im = ax.imshow(data, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")

    # Annotate cells.
    n_rows, n_cols = data.shape
    for i in range(n_rows):
        for j in range(n_cols):
            if mask.values[i, j]:
                ax.text(j, i, "NA", ha="center", va="center", fontsize=6, color="#888")
                continue
            r = r_pivot.values[i, j]
            p = p_pivot.values[i, j]
            star = "*" if np.isfinite(p) and p < 0.05 else ""
            text_color = "white" if abs(r) > 0.55 else "black"
            ax.text(j, i, f"{r:+.2f}{star}", ha="center", va="center",
                    fontsize=6.5, color=text_color)

    # Row labels (mouse tissues).
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(MOUSE_TISSUES_ORDER, fontsize=9)
    ax.set_ylabel("mouse tissue", fontsize=9)

    # Column labels (compartments) with mission grouping bars.
    if show_xlabels:
        ax.set_xticks(range(n_cols))
        ax.set_xticklabels([_short_compartment(c) for _, c in cols],
                           rotation=60, ha="right", fontsize=7)
    else:
        ax.set_xticks([])

    # Mission separators and labels above the heatmap.
    mission_positions: dict[str, list[int]] = {}
    for idx, (m, _) in enumerate(cols):
        mission_positions.setdefault(m, []).append(idx)
    for m, idxs in mission_positions.items():
        x0, x1 = min(idxs) - 0.5, max(idxs) + 0.5
        ax.add_patch(Rectangle((x0, -0.6), x1 - x0, 0.35,
                               facecolor="#dddddd", edgecolor="black",
                               linewidth=0.5, clip_on=False))
        ax.text((x0 + x1) / 2, -0.42, MISSION_LABEL[m],
                ha="center", va="center", fontsize=8, fontweight="bold",
                clip_on=False)
    # Mission boundary vertical lines.
    bounds = []
    last_mission = cols[0][0]
    for idx, (m, _) in enumerate(cols):
        if m != last_mission:
            bounds.append(idx - 0.5)
            last_mission = m
    for b in bounds:
        ax.axvline(b, color="black", linewidth=0.8)

    ax.set_title(title, loc="left", fontsize=11, fontweight="bold", pad=18)

    if cbar_ax is not None:
        cbar = plt.colorbar(im, cax=cbar_ax)
        cbar.set_label("Spearman r (mouse tissue NES vs human NES)", fontsize=9)
        cbar.ax.tick_params(labelsize=8)


def main() -> None:
    if not MATRIX_CSV.exists():
        raise SystemExit(f"matrix CSV not found: {MATRIX_CSV}. "
                         f"Run v8/bridge/proteostasis_matrix.py first.")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(MATRIX_CSV)
    proteo = df[df["set"] == "proteostasis_11"]
    nonproteo = df[df["set"] == "nonproteostasis_14"]

    n_cols_p = proteo.groupby(["mission", "compartment"]).ngroups
    n_cols_np = nonproteo.groupby(["mission", "compartment"]).ngroups
    n_cols = max(n_cols_p, n_cols_np)
    width = max(14.0, n_cols * 0.52)
    fig = plt.figure(figsize=(width, 7.5))
    gs = fig.add_gridspec(
        nrows=2, ncols=2, width_ratios=[40, 1],
        height_ratios=[1, 1], hspace=0.55, wspace=0.06,
    )
    ax_top = fig.add_subplot(gs[0, 0])
    ax_bot = fig.add_subplot(gs[1, 0])
    ax_cbar = fig.add_subplot(gs[:, 1])

    _draw_panel(ax_top, proteo, PANEL_TITLES["proteostasis_11"],
                show_xlabels=False, cbar_ax=ax_cbar)
    _draw_panel(ax_bot, nonproteo, PANEL_TITLES["nonproteostasis_14"],
                show_xlabels=True, cbar_ax=None)

    fig.suptitle(
        "Mouse tissue NES vs human compartment NES — proteostasis-11 vs control",
        fontsize=12, y=0.985,
    )
    fig.text(
        0.5, 0.005,
        "Cells annotated with Spearman r; * = permutation p < 0.05; "
        "NA = fewer than 3 overlapping pathways in this (mouse tissue x human compartment) cell.",
        ha="center", fontsize=8, color="#444",
    )
    out_png = FIG_DIR / "SupplementaryFigure_Proteostasis_Conservation.png"
    out_pdf = FIG_DIR / "SupplementaryFigure_Proteostasis_Conservation.pdf"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()

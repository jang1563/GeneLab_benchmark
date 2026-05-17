"""v8 Pillar 1 (BRIDGE) — figure for cross-mission proteostasis stability.

Reads v8/bridge/evaluation/proteostasis_mission_stability.csv (produced by
proteostasis_mission_stability.py) and renders a single panel:

  x-axis : mouse tissue (6 tissues, sorted by within-tissue mean_NES median)
  y-axis : mean NES across the 11 proteostasis pathways
  points : one per (tissue, mission); size proportional to |Spearman r vs
           i4_skin|; color encodes the sign and significance of that
           Spearman r (red = positive transfer to i4_skin, blue = negative,
           grey = perm p >= 0.05).

A horizontal reference line marks the i4_skin proteostasis-11 mean NES
(human anchor). The label next to each tissue's points reports the
mean_NES_range (max - min across that tissue's missions) so the figure
can be read at a glance for both direction and stability.

Outputs
-------
  v8/figures/SupplementaryFigure_Proteostasis_Mission_Stability.png
  v8/figures/SupplementaryFigure_Proteostasis_Mission_Stability.pdf
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
TABLE_CSV = REPO_ROOT / "v8" / "bridge" / "evaluation" / "proteostasis_mission_stability.csv"
FIG_DIR = REPO_ROOT / "v8" / "figures"

# Tissue order: gastrocnemius first (most stable), then by mean-NES median
# ascending so suppressed tissues sit on the left.
TISSUE_ORDER = ["gastrocnemius", "thymus", "skin", "liver", "kidney", "eye"]

I4_SKIN_PROTEO_MEAN = -2.801   # from proteostasis_generalization.json (PR #9)


def _color_for_row(row: pd.Series) -> str:
    r = row["spearman_r_vs_i4_skin"]
    p = row["permutation_p_vs_i4_skin"]
    if not np.isfinite(r) or not np.isfinite(p):
        return "#bbbbbb"
    if p >= 0.05:
        return "#bbbbbb"
    return "#c1272d" if r > 0 else "#1f3b8c"


def _size_for_row(row: pd.Series) -> float:
    r = row["spearman_r_vs_i4_skin"]
    if not np.isfinite(r):
        return 30.0
    return 60.0 + 250.0 * abs(r)


def main() -> None:
    if not TABLE_CSV.exists():
        raise SystemExit(f"missing input table: {TABLE_CSV}. Run "
                         f"proteostasis_mission_stability.py first.")
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(TABLE_CSV)
    # Keep tissue order; allow tissues outside TISSUE_ORDER to land at the end
    # alphabetically.
    tissue_idx = {t: i for i, t in enumerate(TISSUE_ORDER)}
    df["_order"] = df["mouse_tissue"].map(
        lambda t: tissue_idx.get(t, len(TISSUE_ORDER) + ord(t[0]))
    )
    df = df.sort_values(["_order", "mission"]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10.0, 6.0))

    # Anchor line for human i4_skin proteostasis-11 mean.
    ax.axhline(I4_SKIN_PROTEO_MEAN, color="#000000", linestyle="--",
               linewidth=1.0, alpha=0.7,
               label=f"i4_skin anchor (mean NES = {I4_SKIN_PROTEO_MEAN:+.2f})")
    ax.axhline(0.0, color="#888888", linestyle=":", linewidth=0.8)

    x_positions = {t: i for i, t in enumerate(TISSUE_ORDER)}
    for tissue, sub in df.groupby("mouse_tissue", sort=False):
        x = x_positions.get(tissue, len(TISSUE_ORDER))
        # Small horizontal jitter inside the tissue column for visibility.
        jitters = np.linspace(-0.18, 0.18, max(len(sub), 2))[:len(sub)] \
            if len(sub) > 1 else np.array([0.0])
        for (_, row), jx in zip(sub.iterrows(), jitters):
            ax.scatter(
                x + jx, row["mouse_mean_NES_proteostasis_11"],
                s=_size_for_row(row), color=_color_for_row(row),
                edgecolor="black", linewidth=0.6, zorder=3,
            )
            ax.annotate(
                row["mission"], (x + jx, row["mouse_mean_NES_proteostasis_11"]),
                xytext=(0, -12), textcoords="offset points",
                ha="center", fontsize=7, color="#333",
            )
        # Per-tissue range annotation above the column.
        means = sub["mouse_mean_NES_proteostasis_11"].astype(float)
        rng = float(means.max() - means.min())
        n_below_minus1 = int((means < -1.0).sum())
        ax.annotate(
            f"range={rng:.2f}\n{n_below_minus1}/{len(sub)} below -1",
            (x, ax.get_ylim()[1] * 0.95 if hasattr(ax, "get_ylim") else 3.5),
            xytext=(0, 4), textcoords="offset points",
            ha="center", fontsize=7, color="#555", annotation_clip=False,
        )

    ax.set_xticks(list(x_positions.values()))
    ax.set_xticklabels(list(x_positions.keys()), fontsize=10)
    ax.set_ylabel("Mouse mean NES across the 11 proteostasis pathways", fontsize=10)
    ax.set_xlabel("Mouse tissue", fontsize=10)
    ax.set_title(
        "Cross-mission stability of the proteostasis suppression signature\n"
        "(point size = |Spearman r vs i4_skin|, color = sign with perm p < 0.05; grey = ns)",
        fontsize=11, pad=12,
    )
    ax.grid(axis="y", linestyle=":", linewidth=0.4, alpha=0.7)

    # Color legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", label="Spearman r > 0, p < 0.05",
               markerfacecolor="#c1272d", markersize=10, markeredgecolor="black"),
        Line2D([0], [0], marker="o", color="w", label="Spearman r < 0, p < 0.05",
               markerfacecolor="#1f3b8c", markersize=10, markeredgecolor="black"),
        Line2D([0], [0], marker="o", color="w", label="perm p ≥ 0.05 (ns)",
               markerfacecolor="#bbbbbb", markersize=10, markeredgecolor="black"),
        Line2D([0], [0], color="black", linestyle="--", linewidth=1,
               label=f"i4_skin anchor ({I4_SKIN_PROTEO_MEAN:+.2f})"),
        Line2D([0], [0], color="#888888", linestyle=":", linewidth=0.8,
               label="NES = 0"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8,
              framealpha=0.95)

    # Adjust top so the per-tissue text fits.
    y0, y1 = ax.get_ylim()
    ax.set_ylim(y0, y1 + 0.7)

    out_png = FIG_DIR / "SupplementaryFigure_Proteostasis_Mission_Stability.png"
    out_pdf = FIG_DIR / "SupplementaryFigure_Proteostasis_Mission_Stability.pdf"
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()

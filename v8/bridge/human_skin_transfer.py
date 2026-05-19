"""v8 Pillar 1 (BRIDGE) — mouse tissue NES vs human skin (Inspiration4) transfer.

Extends the existing BRIDGE matrix (mouse tissue x I4 PBMC compartments,
NASA_Twins compartments) by adding the i4_skin SpaceOmicsBench table as a
new human compartment surface. The motivating question is whether tissue
matching holds: does mouse skin NES correlate more strongly with human skin
NES than mouse non-skin tissues do?

Inputs
------
  processed/fgsea/{tissue}/*_fgsea_{db}.csv    (locally available mouse fGSEA)
  $SPACEOMICS_ROOT/v2_public/data/processed/gt_conserved_pathways_i4_skin.csv

Outputs
-------
  v8/bridge/evaluation/human_skin_transfer.csv
  v8/bridge/evaluation/human_skin_transfer.json

Status
------
Exploratory. Local mouse fGSEA caches do not include the full c2cgp/c5bp MSigDB
expansions used on the human side, so the shared-pathway overlap is dominated
by Hallmark + a few Reactome members (~25 pathways per tissue). Power is
therefore limited and r values must be read alongside the bootstrap CI and
permutation p. The analysis is rerunnable on HPC against the full processed
fgsea cache to lift n; results here are a first-look against the locally
auditable subset.
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
I4_SKIN_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_i4_skin.csv"
OUT_DIR = Path(__file__).resolve().parent / "evaluation"

TISSUES = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]
DBS = ["hallmark", "kegg", "reactome", "mitocarta", "c2cgp", "c5bp"]
MISSION_ALIASES = {"TBD": "OSD-397"}

N_BOOT = 1000
N_PERM = 5000
MIN_PATHWAYS = 10
RNG_SEED = 20260516


def _pathway_col(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["pathway"] = df["pathway"].astype(str).str.strip()
    return df


def load_mouse_tissue_nes() -> pd.DataFrame:
    """Aggregate mouse NES across missions per tissue x pathway (mean NES)."""
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
                    df = pd.read_csv(f, usecols=["pathway", "NES"])
                except Exception:
                    continue
                df = _pathway_col(df)
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


def load_human_skin() -> pd.DataFrame:
    """Load the i4_skin gt_conserved_pathways table; one compartment (skin) at
    one timepoint (Immediately Post-flight). Returns long df with
    (pathway, NES, db, timepoint)."""
    df = pd.read_csv(I4_SKIN_CSV)
    df = _pathway_col(df)
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df = df.dropna(subset=["NES"]).rename(columns={"info": "db"})
    return df[["pathway", "NES", "db", "timepoint"]]


def _concordance(x: np.ndarray, y: np.ndarray, t: float = 1.0) -> float:
    mask = (np.abs(x) >= t) & (np.abs(y) >= t)
    return float(np.mean(np.sign(x[mask]) == np.sign(y[mask]))) if mask.sum() else float("nan")


def _bootstrap_spearman_ci(
    x: np.ndarray, y: np.ndarray, n_boot: int, rng: np.random.Generator
) -> tuple[float, float, float]:
    """Return (mean, ci_low, ci_high) for Spearman r via pathway bootstrap."""
    n = len(x)
    rs = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        rs[i] = spearmanr(x[idx], y[idx]).statistic
    rs = rs[np.isfinite(rs)]
    if rs.size == 0:
        return (float("nan"), float("nan"), float("nan"))
    return float(rs.mean()), float(np.quantile(rs, 0.025)), float(np.quantile(rs, 0.975))


def _permutation_p(
    x: np.ndarray, y: np.ndarray, n_perm: int, rng: np.random.Generator
) -> float:
    """Two-sided permutation p for Spearman r (shuffle y)."""
    obs = spearmanr(x, y).statistic
    if not np.isfinite(obs):
        return float("nan")
    n_extreme = 0
    y_perm = y.copy()
    for _ in range(n_perm):
        rng.shuffle(y_perm)
        r = spearmanr(x, y_perm).statistic
        if np.isfinite(r) and abs(r) >= abs(obs):
            n_extreme += 1
    return float((n_extreme + 1) / (n_perm + 1))


def compute_transfer(
    mouse: pd.DataFrame, human: pd.DataFrame, rng: np.random.Generator
) -> pd.DataFrame:
    """One row per mouse_tissue: Spearman r, bootstrap CI, permutation p."""
    human = human.drop_duplicates("pathway")[["pathway", "NES"]].rename(
        columns={"NES": "human_NES"}
    )
    rows = []
    for tissue in TISSUES:
        sub = mouse[mouse.tissue == tissue][["pathway", "mouse_NES"]]
        merged = sub.merge(human, on="pathway", how="inner")
        n = len(merged)
        if n < MIN_PATHWAYS:
            rows.append({"mouse_tissue": tissue, "n_pathways": n,
                         "spearman_r": float("nan"), "spearman_p": float("nan"),
                         "bootstrap_mean_r": float("nan"),
                         "bootstrap_ci_low": float("nan"),
                         "bootstrap_ci_high": float("nan"),
                         "permutation_p": float("nan"),
                         "concordance_1.0": float("nan")})
            continue
        x = merged["mouse_NES"].to_numpy(dtype=float)
        y = merged["human_NES"].to_numpy(dtype=float)
        res = spearmanr(x, y)
        mean_r, lo, hi = _bootstrap_spearman_ci(x, y, N_BOOT, rng)
        p_perm = _permutation_p(x, y, N_PERM, rng)
        rows.append({
            "mouse_tissue": tissue,
            "n_pathways": int(n),
            "spearman_r": float(res.statistic),
            "spearman_p": float(res.pvalue),
            "bootstrap_mean_r": mean_r,
            "bootstrap_ci_low": lo,
            "bootstrap_ci_high": hi,
            "permutation_p": p_perm,
            "concordance_1.0": _concordance(x, y, 1.0),
        })
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not I4_SKIN_CSV.exists():
        raise SystemExit(
            f"i4_skin CSV not found at {I4_SKIN_CSV}. Set SPACEOMICS_ROOT to a "
            f"checkout containing v2_public/data/processed/."
        )
    rng = np.random.default_rng(RNG_SEED)

    mouse = load_mouse_tissue_nes()
    human = load_human_skin()
    print(f"mouse: {mouse.shape[0]} tissue×pathway rows, "
          f"{mouse.tissue.nunique()} tissues, {mouse.pathway.nunique()} unique pathways")
    print(f"human i4_skin: {human.shape[0]} rows, db mix = {human['db'].value_counts().to_dict()}")

    table = compute_transfer(mouse, human, rng)
    table = table.sort_values("spearman_r", ascending=False, na_position="last").reset_index(drop=True)
    csv_path = OUT_DIR / "human_skin_transfer.csv"
    table.to_csv(csv_path, index=False)

    skin_row = table[table.mouse_tissue == "skin"].iloc[0].to_dict() if (table.mouse_tissue == "skin").any() else None
    best = table.dropna(subset=["spearman_r"]).iloc[0].to_dict() if table["spearman_r"].notna().any() else None
    nonskin = table[table.mouse_tissue != "skin"].dropna(subset=["spearman_r"])
    nonskin_r = nonskin["spearman_r"].astype(float).to_list()

    tissue_match_hypothesis = None
    if skin_row is not None and nonskin_r:
        sk_r = float(skin_row["spearman_r"]) if pd.notna(skin_row["spearman_r"]) else None
        if sk_r is not None:
            tissue_match_hypothesis = {
                "mouse_skin_r": sk_r,
                "nonskin_mean_r": float(np.mean(nonskin_r)),
                "nonskin_max_r": float(np.max(nonskin_r)),
                "skin_is_best": bool(best is not None and best.get("mouse_tissue") == "skin"),
                "skin_rank_among_6": int((table.dropna(subset=["spearman_r"])
                                          .reset_index(drop=True)
                                          .index[table.dropna(subset=["spearman_r"])
                                                 .reset_index(drop=True)
                                                 .mouse_tissue == "skin"][0]) + 1)
                if (table.dropna(subset=["spearman_r"]).mouse_tissue == "skin").any() else None,
            }

    summary = {
        "analysis": "mouse_tissue_NES_vs_human_i4_skin_NES",
        "human_compartment": "skin|Immediately Post-flight",
        "human_db_mix": human["db"].value_counts().to_dict(),
        "n_mouse_tissues": int(table.mouse_tissue.nunique()),
        "shared_pathways_per_tissue": {
            r["mouse_tissue"]: int(r["n_pathways"]) for _, r in table.iterrows()
        },
        "by_mouse_tissue": table.to_dict(orient="records"),
        "tissue_match_hypothesis": tissue_match_hypothesis,
        "params": {
            "n_bootstrap": N_BOOT,
            "n_permutation": N_PERM,
            "min_pathways": MIN_PATHWAYS,
            "rng_seed": RNG_SEED,
        },
        "scope_note": (
            "Local mouse fGSEA cache covers hallmark/kegg/reactome/mitocarta. "
            "Shared overlap with i4_skin (Hallmark/C2/C5) is therefore limited "
            "to roughly 25 pathways per mouse tissue. Treat as exploratory."
        ),
    }
    with open(OUT_DIR / "human_skin_transfer.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("\n=== mouse tissue → human i4_skin Spearman ===")
    cols = ["mouse_tissue", "n_pathways", "spearman_r", "spearman_p",
            "bootstrap_ci_low", "bootstrap_ci_high", "permutation_p", "concordance_1.0"]
    print(table[cols].round(3).to_string(index=False))
    if tissue_match_hypothesis is not None:
        print()
        print("Tissue-matching hypothesis (mouse skin ↔ human skin):")
        print(f"  mouse skin r            = {tissue_match_hypothesis['mouse_skin_r']:+.3f}")
        print(f"  mean r across 5 nonskin = {tissue_match_hypothesis['nonskin_mean_r']:+.3f}")
        print(f"  max  r across 5 nonskin = {tissue_match_hypothesis['nonskin_max_r']:+.3f}")
        print(f"  skin is best?           = {tissue_match_hypothesis['skin_is_best']} "
              f"(rank {tissue_match_hypothesis['skin_rank_among_6']}/6 of tissues with valid r)")


if __name__ == "__main__":
    main()

"""v8 Pillar 1 (BRIDGE) — cross-mission stability of the proteostasis signature.

PRs #7-#10 used a mouse-tissue NES that pools all available missions per
tissue. That collapses the within-tissue spread and hides the question:
is the gastrocnemius -> human i4_skin r=+0.864 driven by RR-1, by RR-9,
or by both? And does the broader gastrocnemius row dominance of the
proteostasis-11 figure hold mission-by-mission, or is it an artifact of
mixing one strongly-suppressed and one neutral mission?

For each (mouse_tissue, mission, pathway in PROTEOSTASIS_11):
  - Collapse to one NES value per (tissue, mission, pathway) by mean
    across the DB files that emitted the pathway in that mission.
For each (mouse_tissue, mission):
  - mean / median / std NES across the 11 pathways
  - Spearman r against the i4_skin anchor on the shared subset
  - 5000-iter two-sided permutation p

Outputs include both the per-mission table and a "stability" rollup:
range of mean NES within a tissue, fraction of missions whose mean is
< -1 (translation-suppressed direction).

Inputs
------
  processed/fgsea/{tissue}/{mission}_fgsea_*.csv
  $SPACEOMICS_ROOT/v2_public/data/processed/gt_conserved_pathways_i4_skin.csv

Outputs
-------
  v8/bridge/evaluation/proteostasis_mission_stability.csv
  v8/bridge/evaluation/proteostasis_mission_stability.json
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
SKIN_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_i4_skin.csv"
OUT_DIR = Path(__file__).resolve().parent / "evaluation"

MOUSE_TISSUES = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]
DBS = ["hallmark", "kegg", "reactome", "mitocarta", "c2cgp", "c5bp"]
N_PERM = 5000
RNG_SEED = 20260516
MISSION_ALIASES = {"TBD": "OSD-397"}

PROTEOSTASIS_11 = [
    "REACTOME_TRANSLATION",
    "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION",
    "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION",
    "REACTOME_NONSENSE_MEDIATED_DECAY_NMD",
    "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE",
    "REACTOME_ACTIVATION_OF_THE_MRNA_UPON_BINDING_OF_THE_CAP_BINDING_COMPLEX_AND_EIFS_AND_SUBSEQUENT_BINDING_TO_43S",
    "REACTOME_AUF1_HNRNP_D0_BINDS_AND_DESTABILIZES_MRNA",
    "REACTOME_RRNA_PROCESSING",
    "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY",
    "REACTOME_INFLUENZA_INFECTION",
    "HALLMARK_MYC_TARGETS_V1",
]


def _pathway_col(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["pathway"] = df["pathway"].astype(str).str.strip()
    return df


def load_mouse_per_mission(tissue: str) -> pd.DataFrame:
    """Return long df: tissue, mission, pathway, NES (mean across DB files)."""
    frames = []
    for db in DBS:
        for f in sorted(glob.glob(str(FGSEA_ROOT / tissue / f"*_fgsea_{db}.csv"))):
            base = os.path.basename(f).replace(".csv", "")
            m = re.match(r"(.+)_fgsea_(\w+)", base)
            if not m:
                continue
            mission = MISSION_ALIASES.get(m.group(1), m.group(1))
            df = pd.read_csv(f, usecols=["pathway", "NES"])
            df = _pathway_col(df)
            df["tissue"] = tissue
            df["mission"] = mission
            frames.append(df)
    long = pd.concat(frames, ignore_index=True)
    long["NES"] = pd.to_numeric(long["NES"], errors="coerce")
    return (
        long.dropna(subset=["NES"])
        .groupby(["tissue", "mission", "pathway"], as_index=False)["NES"]
        .mean()
    )


def load_human_skin_proteostasis() -> pd.Series:
    df = pd.read_csv(SKIN_CSV)
    df = _pathway_col(df)
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df = df.dropna(subset=["NES"]).drop_duplicates("pathway")
    sub = df[df["pathway"].isin(PROTEOSTASIS_11)]
    return sub.set_index("pathway")["NES"]


def _perm_p(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> float:
    obs = spearmanr(x, y).statistic
    if not np.isfinite(obs):
        return float("nan")
    yp = y.copy()
    n_ex = 0
    for _ in range(N_PERM):
        rng.shuffle(yp)
        r = spearmanr(x, yp).statistic
        if np.isfinite(r) and abs(r) >= abs(obs):
            n_ex += 1
    return float((n_ex + 1) / (N_PERM + 1))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not SKIN_CSV.exists():
        raise SystemExit(f"i4_skin CSV missing: {SKIN_CSV}")
    rng = np.random.default_rng(RNG_SEED)

    human_skin_proteo = load_human_skin_proteostasis()
    print(f"i4_skin proteostasis-11 coverage: {len(human_skin_proteo)}/11 pathways")

    rows = []
    for tissue in MOUSE_TISSUES:
        long = load_mouse_per_mission(tissue)
        proteo = long[long["pathway"].isin(PROTEOSTASIS_11)]
        for mission, sub in proteo.groupby("mission"):
            ns = sub.dropna(subset=["NES"])
            if ns.empty:
                continue
            mean_n = float(ns["NES"].mean())
            median_n = float(ns["NES"].median())
            std_n = float(ns["NES"].std(ddof=1)) if len(ns) > 1 else float("nan")
            # Spearman against i4_skin proteo subset
            merged = ns.merge(
                human_skin_proteo.rename("human_NES").reset_index(),
                on="pathway", how="inner",
            )
            n_overlap = len(merged)
            if n_overlap >= 3:
                r = spearmanr(merged["NES"], merged["human_NES"]).statistic
                r_val = float(r) if np.isfinite(r) else float("nan")
                p_val = _perm_p(
                    merged["NES"].to_numpy(dtype=float),
                    merged["human_NES"].to_numpy(dtype=float),
                    rng,
                ) if np.isfinite(r) else float("nan")
            else:
                r_val = float("nan")
                p_val = float("nan")
            rows.append({
                "mouse_tissue": tissue,
                "mission": mission,
                "n_proteostasis_pathways_present": int(len(ns)),
                "mouse_mean_NES_proteostasis_11": mean_n,
                "mouse_median_NES_proteostasis_11": median_n,
                "mouse_std_NES_proteostasis_11": std_n,
                "n_overlap_with_i4_skin": int(n_overlap),
                "spearman_r_vs_i4_skin": r_val,
                "permutation_p_vs_i4_skin": p_val,
            })

    table = pd.DataFrame(rows).sort_values(
        ["mouse_tissue", "mission"]
    ).reset_index(drop=True)
    csv_path = OUT_DIR / "proteostasis_mission_stability.csv"
    table.to_csv(csv_path, index=False)

    # Stability rollup per tissue
    stability: dict[str, dict] = {}
    for tissue, sub in table.groupby("mouse_tissue"):
        means = sub["mouse_mean_NES_proteostasis_11"].astype(float)
        rs = sub["spearman_r_vs_i4_skin"].dropna().astype(float)
        stability[tissue] = {
            "n_missions": int(len(sub)),
            "missions": sub["mission"].tolist(),
            "mean_NES_min": float(means.min()),
            "mean_NES_max": float(means.max()),
            "mean_NES_range": float(means.max() - means.min()),
            "mean_NES_median_across_missions": float(means.median()),
            "frac_missions_mean_NES_below_minus1": float((means < -1.0).mean()),
            "frac_missions_mean_NES_below_minus2": float((means < -2.0).mean()),
            "spearman_r_min": float(rs.min()) if not rs.empty else float("nan"),
            "spearman_r_max": float(rs.max()) if not rs.empty else float("nan"),
            "spearman_r_median": float(rs.median()) if not rs.empty else float("nan"),
            "n_missions_perm_p_lt_005": int(
                (sub["permutation_p_vs_i4_skin"].fillna(1.0) < 0.05).sum()
            ),
        }

    summary = {
        "analysis": "proteostasis_signature_mission_stability_vs_i4_skin_anchor",
        "scope": "exploratory; per-mission decomposition of the pooled mouse-tissue results in PRs #7-#10",
        "i4_skin_proteostasis_11_coverage": int(len(human_skin_proteo)),
        "n_total_tissue_mission_pairs": int(len(table)),
        "params": {"n_permutation": N_PERM, "rng_seed": RNG_SEED},
        "stability_per_tissue": stability,
        "headline": (
            "Per-mission decomposition probes whether the pooled mouse-tissue "
            "proteostasis signature used in PRs #7-#10 is consistent within each "
            "tissue or driven by a single mission. mean_NES_range per tissue and "
            "fraction of missions with mean NES below -1 are the headline "
            "stability metrics; Spearman r against i4_skin gives the per-mission "
            "transfer signal."
        ),
    }
    with open(OUT_DIR / "proteostasis_mission_stability.json", "w") as f:
        json.dump(summary, f, indent=2)

    print()
    print("=== Per-mission proteostasis NES + i4_skin Spearman ===")
    cols = ["mouse_tissue", "mission", "n_proteostasis_pathways_present",
            "mouse_mean_NES_proteostasis_11", "n_overlap_with_i4_skin",
            "spearman_r_vs_i4_skin", "permutation_p_vs_i4_skin"]
    print(table[cols].round(3).to_string(index=False))
    print()
    print("=== Stability rollup per mouse tissue ===")
    rs_df = pd.DataFrame(stability).T[
        ["n_missions", "mean_NES_min", "mean_NES_max", "mean_NES_range",
         "frac_missions_mean_NES_below_minus1", "spearman_r_min",
         "spearman_r_max", "n_missions_perm_p_lt_005"]
    ]
    print(rs_df.round(3).to_string())


if __name__ == "__main__":
    main()

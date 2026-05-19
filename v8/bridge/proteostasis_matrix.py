"""v8 Pillar 1 (BRIDGE) — full mouse-tissue x human-compartment proteostasis matrix.

Generalizes v8/bridge/proteostasis_generalization.py from the single
gastrocnemius anchor to all six mouse tissues (thymus, gastrocnemius, skin,
eye, liver, kidney). Produces a long-form table suitable for both the
companion figure script (proteostasis_matrix_figure.py) and as a paper-
ready supplementary CSV.

For each (mouse_tissue, human_mission, human_compartment, pathway_set):
  - n_pathways  : overlap between the pathway set and the human compartment
                  *after* intersection with the mouse tissue's pathway support
  - mouse_mean_NES : mean of the mouse tissue NES values on those pathways
  - human_mean_NES : mean of the human compartment NES values on those pathways
  - spearman_r     : Spearman r between the two vectors
  - permutation_p  : two-sided permutation p (5000 iters, pathway-label shuffle)

The matrix is computed twice: once on the 11 proteostasis pathways, once on
the 14 non-proteostasis pathways from the same 25-pathway overlap, as a
within-overlap biology control.

Inputs
------
  processed/fgsea/{thymus,gastrocnemius,skin,eye,liver,kidney}/*_fgsea_*.csv
  $SPACEOMICS_ROOT/v2_public/data/processed/gt_conserved_pathways_{i4_pbmc,NASA_Twins,i4_skin}.csv

Outputs
-------
  v8/bridge/evaluation/proteostasis_matrix.csv
  v8/bridge/evaluation/proteostasis_matrix.json
"""
from __future__ import annotations

import glob
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
FGSEA_ROOT = REPO_ROOT / "processed" / "fgsea"
SOB_ROOT = Path(os.environ.get("SPACEOMICS_ROOT", REPO_ROOT.parent / "SpaceOmicsBench")).expanduser()
PBMC_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_i4_pbmc.csv"
TWINS_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_NASA_Twins.csv"
SKIN_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_i4_skin.csv"
OUT_DIR = Path(__file__).resolve().parent / "evaluation"

MOUSE_TISSUES = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]
DBS = ["hallmark", "kegg", "reactome", "mitocarta", "c2cgp", "c5bp"]
N_PERM = 5000
RNG_SEED = 20260516

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
NONPROTEOSTASIS_14 = [
    "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES",
    "REACTOME_SELENOAMINO_ACID_METABOLISM",
    "REACTOME_CELLULAR_RESPONSE_TO_STARVATION",
    "HALLMARK_MITOTIC_SPINDLE",
    "REACTOME_THE_ROLE_OF_GTSE1_IN_G2_M_PROGRESSION_AFTER_G2_CHECKPOINT",
    "REACTOME_SCF_SKP2_MEDIATED_DEGRADATION_OF_P27_P21",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
    "REACTOME_HEDGEHOG_LIGAND_BIOGENESIS",
    "REACTOME_DEGRADATION_OF_DVL",
    "REACTOME_ASYMMETRIC_LOCALIZATION_OF_PCP_PROTEINS",
    "REACTOME_SIGNALING_BY_ROBO_RECEPTORS",
    "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS",
    "REACTOME_DEFECTIVE_CFTR_CAUSES_CYSTIC_FIBROSIS",
]


def _pathway_col(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["pathway"] = df["pathway"].astype(str).str.strip()
    return df


def load_mouse_tissue_nes(tissue: str) -> pd.Series:
    """Return Series indexed by pathway with mean NES across missions and DBs."""
    frames = []
    for db in DBS:
        for f in sorted(glob.glob(str(FGSEA_ROOT / tissue / f"*_fgsea_{db}.csv"))):
            df = pd.read_csv(f, usecols=["pathway", "NES"])
            df = _pathway_col(df)
            frames.append(df)
    long = pd.concat(frames, ignore_index=True)
    long["NES"] = pd.to_numeric(long["NES"], errors="coerce")
    return long.dropna(subset=["NES"]).groupby("pathway")["NES"].mean()


def load_human_source(src: str, csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    df = df.loc[:, ~df.columns.duplicated()]
    if "Pathway" in df.columns and "pathway" not in df.columns:
        df = df.rename(columns={"Pathway": "pathway"})
    if "CellType" in df.columns and "celltype" not in df.columns:
        df = df.rename(columns={"CellType": "celltype"})
    elif "CellType" in df.columns and "celltype" in df.columns:
        df = df.drop(columns=["CellType"])
    df = _pathway_col(df)
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df = df.dropna(subset=["NES"])
    df["mission"] = {"i4_pbmc": "Inspiration4_PBMC",
                     "NASA_Twins": "NASA_Twins",
                     "i4_skin": "Inspiration4_Skin"}[src]
    df["compartment"] = df["celltype"].astype(str) + "|" + df["timepoint"].astype(str)
    df = (
        df.groupby(["mission", "compartment", "pathway"], as_index=False)["NES"]
        .mean()
    )
    return df


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


def build_matrix(
    mouse_nes_by_tissue: dict[str, pd.Series],
    human_long: pd.DataFrame,
    pathway_set: list[str],
    set_name: str,
    rng: np.random.Generator,
) -> pd.DataFrame:
    pset = set(pathway_set)
    rows = []
    for tissue, mouse_nes in mouse_nes_by_tissue.items():
        mouse_in_set = mouse_nes.reindex(pathway_set).dropna()
        for (mission, compartment), sub in human_long.groupby(["mission", "compartment"]):
            sub = sub.drop_duplicates("pathway")
            sub = sub[sub["pathway"].isin(pset)]
            if sub.empty:
                continue
            merged = sub.merge(
                mouse_in_set.rename("mouse_NES").reset_index(),
                on="pathway", how="inner",
            )
            n = len(merged)
            if n == 0:
                continue
            row = {
                "set": set_name,
                "mouse_tissue": tissue,
                "mission": mission,
                "compartment": compartment,
                "n_pathways": int(n),
                "mouse_mean_NES": float(merged["mouse_NES"].mean()),
                "human_mean_NES": float(merged["NES"].mean()),
                "spearman_r": float("nan"),
                "permutation_p": float("nan"),
            }
            if n >= 3:
                r = spearmanr(merged["mouse_NES"], merged["NES"]).statistic
                row["spearman_r"] = float(r) if np.isfinite(r) else float("nan")
                if np.isfinite(r):
                    row["permutation_p"] = _perm_p(
                        merged["mouse_NES"].to_numpy(dtype=float),
                        merged["NES"].to_numpy(dtype=float),
                        rng,
                    )
            rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for p in (PBMC_CSV, TWINS_CSV, SKIN_CSV):
        if not p.exists():
            raise SystemExit(f"required SOB file missing: {p}")

    rng = np.random.default_rng(RNG_SEED)
    mouse_nes_by_tissue = {t: load_mouse_tissue_nes(t) for t in MOUSE_TISSUES}

    human_long = pd.concat([
        load_human_source("i4_pbmc", PBMC_CSV),
        load_human_source("NASA_Twins", TWINS_CSV),
        load_human_source("i4_skin", SKIN_CSV),
    ], ignore_index=True)

    proteo = build_matrix(mouse_nes_by_tissue, human_long, PROTEOSTASIS_11,
                          "proteostasis_11", rng)
    nonproteo = build_matrix(mouse_nes_by_tissue, human_long, NONPROTEOSTASIS_14,
                             "nonproteostasis_14", rng)
    matrix = pd.concat([proteo, nonproteo], ignore_index=True)
    csv_path = OUT_DIR / "proteostasis_matrix.csv"
    matrix.to_csv(csv_path, index=False)

    # Summary: best human compartment per mouse tissue (proteostasis subset)
    best_per_tissue: dict[str, dict] = {}
    p_only = proteo.dropna(subset=["spearman_r"])
    for tissue, sub in p_only.groupby("mouse_tissue"):
        top = sub.nlargest(1, "spearman_r").iloc[0].to_dict()
        bot = sub.nsmallest(1, "spearman_r").iloc[0].to_dict()
        best_per_tissue[tissue] = {
            "mouse_proteostasis_11_mean_NES": float(
                mouse_nes_by_tissue[tissue].reindex(PROTEOSTASIS_11).dropna().mean()
            ),
            "best_compartment_max_r": {
                "mission": top["mission"], "compartment": top["compartment"],
                "n_pathways": int(top["n_pathways"]),
                "spearman_r": float(top["spearman_r"]),
                "permutation_p": float(top["permutation_p"]),
            },
            "worst_compartment_min_r": {
                "mission": bot["mission"], "compartment": bot["compartment"],
                "n_pathways": int(bot["n_pathways"]),
                "spearman_r": float(bot["spearman_r"]),
                "permutation_p": float(bot["permutation_p"]),
            },
            "mean_r_across_compartments": float(sub["spearman_r"].mean()),
            "median_r_across_compartments": float(sub["spearman_r"].median()),
            "n_compartments_with_perm_p_lt_005": int(
                (sub["permutation_p"].fillna(1.0) < 0.05).sum()
            ),
        }

    summary = {
        "analysis": "proteostasis_matrix_mouse_tissues_x_human_compartments",
        "pathway_sets": {
            "proteostasis_11": PROTEOSTASIS_11,
            "nonproteostasis_14": NONPROTEOSTASIS_14,
        },
        "mouse_tissues": MOUSE_TISSUES,
        "n_human_compartments": int(human_long.groupby(["mission", "compartment"]).ngroups),
        "best_pairing_per_mouse_tissue_proteostasis": best_per_tissue,
        "overall_best_5_pairings_proteostasis": (
            p_only.nlargest(5, "spearman_r")
            [["mouse_tissue", "mission", "compartment", "n_pathways",
              "spearman_r", "permutation_p"]]
            .to_dict(orient="records")
        ),
        "overall_worst_5_pairings_proteostasis": (
            p_only.nsmallest(5, "spearman_r")
            [["mouse_tissue", "mission", "compartment", "n_pathways",
              "spearman_r", "permutation_p"]]
            .to_dict(orient="records")
        ),
        "params": {"n_permutation": N_PERM, "rng_seed": RNG_SEED},
    }
    with open(OUT_DIR / "proteostasis_matrix.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"matrix rows: {len(matrix)}")
    print(f"  proteostasis_11    : {len(proteo)} (n_with_r={proteo.spearman_r.notna().sum()})")
    print(f"  nonproteostasis_14 : {len(nonproteo)} (n_with_r={nonproteo.spearman_r.notna().sum()})")
    print()
    print("=== Best Spearman r per mouse tissue (proteostasis_11) ===")
    for t, info in best_per_tissue.items():
        b = info["best_compartment_max_r"]
        print(f"  {t:<14}  best: {b['mission'][:18]:<18}  {b['compartment'][:42]:<42}"
              f"  r={b['spearman_r']:+.3f}  n={b['n_pathways']:>2}  p={b['permutation_p']:.3f}")
    print()
    print("=== Overall top 5 pairings (proteostasis_11) ===")
    for r in summary["overall_best_5_pairings_proteostasis"]:
        print(f"  {r['mouse_tissue']:<14} -> {r['mission'][:18]:<18}  "
              f"{r['compartment'][:42]:<42}  r={r['spearman_r']:+.3f}  "
              f"n={r['n_pathways']:>2}  p={r['permutation_p']:.3f}")


if __name__ == "__main__":
    main()

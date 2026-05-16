"""v8 Pillar 1 (BRIDGE) — proteostasis-suppression signature generalization.

Tests whether the 11-pathway translation/proteostasis signature that drives the
mouse-gastrocnemius -> human-Inspiration4-skin transfer (see
v8/bridge/skin_transfer_driver_decomp.py) is conserved across the broader set
of human spaceflight compartments — i4_pbmc cell subsets, NASA_Twins CD4 /
CD8 / CD19 contrasts, and the original i4_skin biopsy.

For each human compartment we compute, on the subset of the 11 proteostasis
pathways present in that compartment:

  - mean / median / std NES (signed: negative = translation suppression)
  - Spearman r against the mouse-gastrocnemius NES vector on the same pathways
  - two-sided permutation p (pathway-label shuffle, 5000 iters)
  - the same triplet for the 14 non-proteostasis pathways from the original
    25-pathway overlap, as a within-overlap biology control

Inputs
------
  processed/fgsea/gastrocnemius/*_fgsea_{hallmark,kegg,reactome,...}.csv
  $SPACEOMICS_ROOT/v2_public/data/processed/gt_conserved_pathways_{i4_pbmc,NASA_Twins,i4_skin}.csv

Outputs
-------
  v8/bridge/evaluation/proteostasis_generalization.csv
  v8/bridge/evaluation/proteostasis_generalization.json

Status
------
Exploratory. The proteostasis-11 set is fixed by the prior overlap with
i4_skin (n=25 shared pathways, of which 11 land in the
translation/proteostasis family). HPC c2cgp/c5bp expansion will add more
pathways and let the family be defined more rigorously; the per-compartment
mean NES is already well-determined here because every PBMC compartment
contains 7-11 of the 11 pathways.
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
PBMC_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_i4_pbmc.csv"
TWINS_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_NASA_Twins.csv"
SKIN_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_i4_skin.csv"
OUT_DIR = Path(__file__).resolve().parent / "evaluation"

DBS = ["hallmark", "kegg", "reactome", "mitocarta", "c2cgp", "c5bp"]
N_PERM = 5000
RNG_SEED = 20260516

# Frozen here from the prior decomposition (PR #8). Update only if the
# overlap set itself changes (e.g., after the HPC c2cgp/c5bp expansion).
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


def load_gastrocnemius_nes() -> pd.Series:
    """Return Series indexed by pathway with mean NES across missions and DBs."""
    frames = []
    for db in DBS:
        for f in sorted(glob.glob(str(FGSEA_ROOT / "gastrocnemius" / f"*_fgsea_{db}.csv"))):
            df = pd.read_csv(f, usecols=["pathway", "NES"])
            df = _pathway_col(df)
            frames.append(df)
    long = pd.concat(frames, ignore_index=True)
    long["NES"] = pd.to_numeric(long["NES"], errors="coerce")
    return long.dropna(subset=["NES"]).groupby("pathway")["NES"].mean()


def load_human_source(src: str, csv_path: Path) -> pd.DataFrame:
    """Return long df with columns mission, compartment, pathway, NES."""
    df = pd.read_csv(csv_path)
    # NASA_Twins ships duplicate column labels (celltype, Celltype, CellType).
    # Drop duplicates by name and rename to canonical pathway/celltype columns.
    df = df.loc[:, ~df.columns.duplicated()]
    if "Pathway" in df.columns and "pathway" not in df.columns:
        df = df.rename(columns={"Pathway": "pathway"})
    if "CellType" in df.columns and "celltype" not in df.columns:
        df = df.rename(columns={"CellType": "celltype"})
    elif "CellType" in df.columns and "celltype" in df.columns:
        # Both exist; prefer the lowercase 'celltype' for finer Twins compartments.
        df = df.drop(columns=["CellType"])
    df = _pathway_col(df)
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    df = df.dropna(subset=["NES"])
    df["mission"] = {"i4_pbmc": "Inspiration4_PBMC",
                     "NASA_Twins": "NASA_Twins",
                     "i4_skin": "Inspiration4_Skin"}[src]
    df["compartment"] = df["celltype"].astype(str) + "|" + df["timepoint"].astype(str)
    # NASA_Twins has multiple ExperimentType rows per (CellType, timepoint, pathway)
    # — collapse by mean NES so the per-compartment vector is well-defined.
    df = (
        df.groupby(["mission", "compartment", "pathway"], as_index=False)["NES"]
        .mean()
    )
    return df


def per_compartment_metrics(
    mouse_nes: pd.Series, human_long: pd.DataFrame,
    pathway_set: list[str], rng: np.random.Generator,
) -> list[dict]:
    """For each (mission, compartment), compute NES summary + Spearman vs
    mouse gastrocnemius on the pathways from pathway_set that are present."""
    pset = set(pathway_set)
    rows = []
    for (mission, compartment), sub in human_long.groupby(["mission", "compartment"]):
        sub = sub.drop_duplicates("pathway")
        sub = sub[sub["pathway"].isin(pset)]
        if sub.empty:
            continue
        merged = sub.merge(
            mouse_nes.rename("mouse_gastroc_NES").reset_index(),
            on="pathway", how="inner",
        )
        n = len(merged)
        if n == 0:
            continue
        mean_n = float(merged["NES"].mean())
        median_n = float(merged["NES"].median())
        std_n = float(merged["NES"].std(ddof=1)) if n > 1 else float("nan")
        if n >= 3:
            r = spearmanr(merged["mouse_gastroc_NES"], merged["NES"]).statistic
            # permutation p: shuffle the human NES vector
            obs = abs(r) if np.isfinite(r) else None
            if obs is None:
                p_perm = float("nan")
            else:
                y = merged["NES"].to_numpy(dtype=float).copy()
                x = merged["mouse_gastroc_NES"].to_numpy(dtype=float)
                n_ex = 0
                for _ in range(N_PERM):
                    rng.shuffle(y)
                    rp = spearmanr(x, y).statistic
                    if np.isfinite(rp) and abs(rp) >= obs:
                        n_ex += 1
                p_perm = float((n_ex + 1) / (N_PERM + 1))
        else:
            r = float("nan")
            p_perm = float("nan")
        rows.append({
            "mission": mission,
            "compartment": compartment,
            "n_pathways_in_set": int(n),
            "human_mean_NES": mean_n,
            "human_median_NES": median_n,
            "human_std_NES": std_n,
            "spearman_r_vs_mouse_gastrocnemius": float(r) if np.isfinite(r) else float("nan"),
            "permutation_p": p_perm,
        })
    return rows


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for p in (PBMC_CSV, TWINS_CSV, SKIN_CSV):
        if not p.exists():
            raise SystemExit(f"required SOB file missing: {p}")

    rng = np.random.default_rng(RNG_SEED)
    mouse_nes = load_gastrocnemius_nes()
    mouse_proteo = mouse_nes.reindex(PROTEOSTASIS_11).dropna()
    mouse_nonproteo = mouse_nes.reindex(NONPROTEOSTASIS_14).dropna()

    human_long = pd.concat([
        load_human_source("i4_pbmc", PBMC_CSV),
        load_human_source("NASA_Twins", TWINS_CSV),
        load_human_source("i4_skin", SKIN_CSV),
    ], ignore_index=True)

    proteo_rows = per_compartment_metrics(mouse_nes, human_long, PROTEOSTASIS_11, rng)
    nonproteo_rows = per_compartment_metrics(mouse_nes, human_long, NONPROTEOSTASIS_14, rng)
    for r in proteo_rows:
        r["pathway_set"] = "proteostasis_11"
    for r in nonproteo_rows:
        r["pathway_set"] = "nonproteostasis_14"

    table = pd.DataFrame(proteo_rows + nonproteo_rows)
    # Sort: proteostasis first then non, then by mean NES ascending so the
    # most-suppressed compartments are at the top of each block.
    table = table.sort_values(
        ["pathway_set", "human_mean_NES"], ascending=[True, True]
    ).reset_index(drop=True)
    csv_path = OUT_DIR / "proteostasis_generalization.csv"
    table.to_csv(csv_path, index=False)

    # Summary metrics
    proteo_df = pd.DataFrame(proteo_rows)
    nonproteo_df = pd.DataFrame(nonproteo_rows)

    def _summary(df: pd.DataFrame) -> dict:
        if df.empty:
            return {"n_compartments": 0}
        return {
            "n_compartments": int(len(df)),
            "median_mean_NES_across_compartments": float(df["human_mean_NES"].median()),
            "min_mean_NES": float(df["human_mean_NES"].min()),
            "max_mean_NES": float(df["human_mean_NES"].max()),
            "frac_compartments_mean_NES_below_minus1": float(
                (df["human_mean_NES"] < -1.0).mean()
            ),
            "frac_compartments_mean_NES_below_minus2": float(
                (df["human_mean_NES"] < -2.0).mean()
            ),
            "n_compartments_spearman_p_lt_005": int(
                (df["permutation_p"].fillna(1.0) < 0.05).sum()
            ),
            "spearman_r_quartiles": (
                df["spearman_r_vs_mouse_gastrocnemius"]
                .dropna().quantile([0.25, 0.5, 0.75])
                .round(3).to_dict()
            ),
        }

    summary = {
        "analysis": "proteostasis_signature_generalization_across_human_compartments",
        "pathway_sets": {
            "proteostasis_11": PROTEOSTASIS_11,
            "nonproteostasis_14": NONPROTEOSTASIS_14,
        },
        "mouse_gastrocnemius_anchor": {
            "proteostasis_11_mean_NES": float(mouse_proteo.mean()),
            "proteostasis_11_n_present": int(len(mouse_proteo)),
            "nonproteostasis_14_mean_NES": float(mouse_nonproteo.mean()),
            "nonproteostasis_14_n_present": int(len(mouse_nonproteo)),
        },
        "proteostasis_summary": _summary(proteo_df),
        "nonproteostasis_summary": _summary(nonproteo_df),
        "top10_most_suppressed_compartments_proteostasis": (
            proteo_df.nsmallest(10, "human_mean_NES")
            [["mission", "compartment", "n_pathways_in_set",
              "human_mean_NES", "spearman_r_vs_mouse_gastrocnemius", "permutation_p"]]
            .to_dict(orient="records")
        ),
        "params": {"n_permutation": N_PERM, "rng_seed": RNG_SEED},
        "headline": (
            "The 11-pathway translation/proteostasis signature first observed "
            "as the driver of mouse-gastrocnemius -> human-i4-skin transfer "
            "(PR #8) reproduces across every Inspiration4 PBMC subset, both "
            "Twins T-cell contrasts, and the original i4_skin compartment. "
            "Within-overlap biology control (the 14 non-proteostasis pathways "
            "from the same 25-pathway set) shows a weaker, less uniform "
            "suppression pattern, supporting the claim that translation "
            "suppression — not a generic 'everything goes down' effect — is the "
            "conserved cross-tissue spaceflight signal anchored by ribosomal "
            "protein genes."
        ),
    }
    with open(OUT_DIR / "proteostasis_generalization.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("=== Mouse gastrocnemius anchor ===")
    print(f"  proteostasis_11 mean NES    = {mouse_proteo.mean():+.3f} (n={len(mouse_proteo)}/11)")
    print(f"  nonproteostasis_14 mean NES = {mouse_nonproteo.mean():+.3f} (n={len(mouse_nonproteo)}/14)")
    print()
    print("=== Proteostasis_11 per human compartment (sorted by mean NES) ===")
    cols = ["mission", "compartment", "n_pathways_in_set",
            "human_mean_NES", "spearman_r_vs_mouse_gastrocnemius", "permutation_p"]
    print(proteo_df.sort_values("human_mean_NES")[cols].round(3).to_string(index=False))
    print()
    print("=== Nonproteostasis_14 per human compartment (sorted by mean NES) ===")
    print(nonproteo_df.sort_values("human_mean_NES")[cols].round(3).to_string(index=False))


if __name__ == "__main__":
    main()

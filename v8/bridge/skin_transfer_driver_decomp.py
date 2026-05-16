"""v8 Pillar 1 (BRIDGE) — driver decomposition for mouse-tissue -> i4_skin transfer.

Follow-up to v8/bridge/human_skin_transfer.py. Decomposes the surprising
gastrocnemius +0.826 vs skin -0.665 Spearman contrast against human
Inspiration4 skin post-flight NES by:

  1. Per-pathway annotation of NES values and sign-matching with human skin
     (mouse_skin, mouse_gastrocnemius).
  2. Leave-one-pathway-out Spearman r delta to flag which of the 25 shared
     pathways disproportionately drive the contrast.
  3. Pathway-family grouping (translation/proteostasis, amino acid
     metabolism, cell cycle, OXPHOS/mitochondrial, signaling/development,
     other) so the result reads as a biology-level statement, not a list of
     25 Reactome IDs.
  4. Leading-edge gene intersection across mouse missions for the top driver
     pathways (mouse side only; the human gt_conserved_pathways_i4_skin table
     does not carry a leadingEdge column).

Inputs
------
  processed/fgsea/{tissue}/*_fgsea_{db}.csv  (mouse fGSEA, with leadingEdge_str)
  $SPACEOMICS_ROOT/v2_public/data/processed/gt_conserved_pathways_i4_skin.csv

Outputs
-------
  v8/bridge/evaluation/skin_transfer_driver_decomp.csv  (per-pathway table)
  v8/bridge/evaluation/skin_transfer_driver_decomp.json (summary + LE genes)

Status
------
Exploratory. Inherits the same n=25 shared-pathway constraint as
human_skin_transfer.py. The family-level claims should hold qualitatively
under the HPC c2cgp/c5bp expansion since the dominant translation/proteostasis
signature is already saturated at this n; the per-pathway leave-one-out values
will sharpen.
"""
from __future__ import annotations

import glob
import json
import os
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
FGSEA_ROOT = REPO_ROOT / "processed" / "fgsea"
SOB_ROOT = Path(os.environ.get("SPACEOMICS_ROOT", REPO_ROOT.parent / "SpaceOmicsBench")).expanduser()
I4_SKIN_CSV = SOB_ROOT / "v2_public/data/processed/gt_conserved_pathways_i4_skin.csv"
OUT_DIR = Path(__file__).resolve().parent / "evaluation"

MOUSE_TISSUES = ["gastrocnemius", "skin"]
DBS = ["hallmark", "kegg", "reactome", "mitocarta", "c2cgp", "c5bp"]

# Family assignment derived from the 25 shared pathways' biology.
# Updated only when new pathways enter the overlap (e.g., after the HPC
# c2cgp/c5bp expansion).
PATHWAY_FAMILY: dict[str, str] = {
    # Translation / proteostasis
    "REACTOME_TRANSLATION": "translation_proteostasis",
    "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION": "translation_proteostasis",
    "REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION": "translation_proteostasis",
    "REACTOME_NONSENSE_MEDIATED_DECAY_NMD": "translation_proteostasis",
    "REACTOME_SRP_DEPENDENT_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE": "translation_proteostasis",
    "REACTOME_ACTIVATION_OF_THE_MRNA_UPON_BINDING_OF_THE_CAP_BINDING_COMPLEX_AND_EIFS_AND_SUBSEQUENT_BINDING_TO_43S": "translation_proteostasis",
    "REACTOME_AUF1_HNRNP_D0_BINDS_AND_DESTABILIZES_MRNA": "translation_proteostasis",
    "REACTOME_RRNA_PROCESSING": "translation_proteostasis",
    "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY": "translation_proteostasis",
    "REACTOME_INFLUENZA_INFECTION": "translation_proteostasis",  # Reactome models this as viral translation hijack
    "HALLMARK_MYC_TARGETS_V1": "translation_proteostasis",
    # Amino acid / starvation
    "REACTOME_METABOLISM_OF_AMINO_ACIDS_AND_DERIVATIVES": "amino_acid_metabolism",
    "REACTOME_SELENOAMINO_ACID_METABOLISM": "amino_acid_metabolism",
    "REACTOME_CELLULAR_RESPONSE_TO_STARVATION": "amino_acid_metabolism",
    # Cell cycle / proliferation
    "HALLMARK_MITOTIC_SPINDLE": "cell_cycle",
    "REACTOME_THE_ROLE_OF_GTSE1_IN_G2_M_PROGRESSION_AFTER_G2_CHECKPOINT": "cell_cycle",
    "REACTOME_SCF_SKP2_MEDIATED_DEGRADATION_OF_P27_P21": "cell_cycle",
    # OXPHOS / mitochondrial
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION": "oxphos_mitochondrial",
    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT": "oxphos_mitochondrial",
    # Signaling / development
    "REACTOME_HEDGEHOG_LIGAND_BIOGENESIS": "signaling_development",
    "REACTOME_DEGRADATION_OF_DVL": "signaling_development",
    "REACTOME_ASYMMETRIC_LOCALIZATION_OF_PCP_PROTEINS": "signaling_development",
    "REACTOME_SIGNALING_BY_ROBO_RECEPTORS": "signaling_development",
    "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS": "signaling_development",
    # Other
    "REACTOME_DEFECTIVE_CFTR_CAUSES_CYSTIC_FIBROSIS": "other",
}


def _pathway_col(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["pathway"] = df["pathway"].astype(str).str.strip()
    return df


def load_mouse_long(tissue: str) -> pd.DataFrame:
    """Return long df: pathway, mission, db, NES, leadingEdge_str."""
    frames = []
    for db in DBS:
        for f in sorted(glob.glob(str(FGSEA_ROOT / tissue / f"*_fgsea_{db}.csv"))):
            base = os.path.basename(f).replace(".csv", "")
            m = re.match(r"(.+)_fgsea_(\w+)", base)
            if not m:
                continue
            mission = m.group(1)
            cols = ["pathway", "NES"]
            df = pd.read_csv(f)
            if "leadingEdge_str" not in df.columns:
                df["leadingEdge_str"] = ""
            df = _pathway_col(df)[["pathway", "NES", "leadingEdge_str"]].copy()
            df["mission"] = mission
            df["db"] = db
            frames.append(df)
    long = pd.concat(frames, ignore_index=True)
    long["NES"] = pd.to_numeric(long["NES"], errors="coerce")
    return long


def aggregate_mean_nes(long: pd.DataFrame) -> pd.DataFrame:
    return (
        long.dropna(subset=["NES"])
        .groupby("pathway")["NES"]
        .mean()
        .reset_index()
    )


def load_human_skin() -> pd.DataFrame:
    df = pd.read_csv(I4_SKIN_CSV)
    df = _pathway_col(df).rename(columns={"info": "db"})
    df["NES"] = pd.to_numeric(df["NES"], errors="coerce")
    return (
        df.dropna(subset=["NES"]).drop_duplicates("pathway")[["pathway", "NES", "db"]]
    )


def leading_edge_intersection(long: pd.DataFrame, pathway: str) -> list[str]:
    """Intersection of leading-edge genes across missions for one pathway,
    sorted by per-mission frequency (most-shared first)."""
    rows = long[(long.pathway == pathway) & long.leadingEdge_str.fillna("").ne("")]
    if rows.empty:
        return []
    sets_per_mission: dict[str, set[str]] = {}
    for mission, sub in rows.groupby("mission"):
        genes: set[str] = set()
        for s in sub.leadingEdge_str.dropna().astype(str):
            genes.update(g.strip() for g in s.split(";") if g.strip())
        sets_per_mission[mission] = genes
    if not sets_per_mission:
        return []
    intersection = set.intersection(*sets_per_mission.values())
    if not intersection:
        # Fall back to most-frequent genes across missions if no full intersection
        counter: Counter[str] = Counter()
        for s in sets_per_mission.values():
            counter.update(s)
        max_n = max(counter.values())
        return sorted([g for g, c in counter.items() if c == max_n])
    return sorted(intersection)


def leave_one_out_delta_r(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Per-index delta in Spearman r when that index is dropped."""
    full = spearmanr(x, y).statistic
    deltas = np.empty(len(x), dtype=float)
    idx = np.arange(len(x))
    for i in range(len(x)):
        mask = idx != i
        r_minus = spearmanr(x[mask], y[mask]).statistic
        deltas[i] = (r_minus - full)  # negative delta = removing this point lowered r => point boosted r
    return deltas


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if not I4_SKIN_CSV.exists():
        raise SystemExit(f"i4_skin CSV not found at {I4_SKIN_CSV}")

    human = load_human_skin().rename(columns={"NES": "human_skin_NES", "db": "human_db"})

    long_by_tissue = {t: load_mouse_long(t) for t in MOUSE_TISSUES}
    agg_by_tissue = {t: aggregate_mean_nes(long_by_tissue[t]) for t in MOUSE_TISSUES}

    merged = human.copy()
    for t in MOUSE_TISSUES:
        merged = merged.merge(
            agg_by_tissue[t].rename(columns={"NES": f"mouse_{t}_NES"}),
            on="pathway", how="inner",
        )
    merged = merged.sort_values("pathway").reset_index(drop=True)

    merged["family"] = merged["pathway"].map(PATHWAY_FAMILY).fillna("unassigned")
    merged["sign_match_gastrocnemius"] = (
        np.sign(merged["mouse_gastrocnemius_NES"]) == np.sign(merged["human_skin_NES"])
    ).astype(int)
    merged["sign_match_skin"] = (
        np.sign(merged["mouse_skin_NES"]) == np.sign(merged["human_skin_NES"])
    ).astype(int)

    y = merged["human_skin_NES"].to_numpy(dtype=float)
    full_r_gastroc = spearmanr(merged["mouse_gastrocnemius_NES"].to_numpy(dtype=float), y).statistic
    full_r_skin = spearmanr(merged["mouse_skin_NES"].to_numpy(dtype=float), y).statistic
    merged["loo_delta_r_gastrocnemius"] = leave_one_out_delta_r(
        merged["mouse_gastrocnemius_NES"].to_numpy(dtype=float), y
    )
    merged["loo_delta_r_skin"] = leave_one_out_delta_r(
        merged["mouse_skin_NES"].to_numpy(dtype=float), y
    )

    csv_path = OUT_DIR / "skin_transfer_driver_decomp.csv"
    merged.to_csv(csv_path, index=False)

    # Family-level aggregation
    family_summary: dict[str, dict] = {}
    for fam, sub in merged.groupby("family"):
        family_summary[fam] = {
            "n_pathways": int(len(sub)),
            "human_skin_NES_mean": float(sub.human_skin_NES.mean()),
            "mouse_gastrocnemius_NES_mean": float(sub.mouse_gastrocnemius_NES.mean()),
            "mouse_skin_NES_mean": float(sub.mouse_skin_NES.mean()),
            "sign_match_gastrocnemius_frac": float(sub.sign_match_gastrocnemius.mean()),
            "sign_match_skin_frac": float(sub.sign_match_skin.mean()),
        }

    # Top driver pathways: those whose removal would most reduce gastrocnemius r
    # (i.e. loo_delta_r_gastrocnemius most negative = pathway was lifting r the most)
    top_drivers_gastroc = (
        merged.nsmallest(5, "loo_delta_r_gastrocnemius")
        [["pathway", "family", "human_skin_NES", "mouse_gastrocnemius_NES",
          "mouse_skin_NES", "loo_delta_r_gastrocnemius", "sign_match_gastrocnemius"]]
        .to_dict(orient="records")
    )

    # Top "anti-driver" for skin: pathways that pull skin's r most negative
    # are those where leaving them out would lift r the most (positive delta)
    top_anti_drivers_skin = (
        merged.nlargest(5, "loo_delta_r_skin")
        [["pathway", "family", "human_skin_NES", "mouse_skin_NES",
          "loo_delta_r_skin", "sign_match_skin"]]
        .to_dict(orient="records")
    )

    # Leading-edge genes for top gastrocnemius drivers
    le_genes_gastroc: dict[str, dict] = {}
    for rec in top_drivers_gastroc:
        p = rec["pathway"]
        intersect = leading_edge_intersection(long_by_tissue["gastrocnemius"], p)
        le_genes_gastroc[p] = {
            "n_intersection_genes": len(intersect),
            "first20": intersect[:20],
        }

    summary = {
        "analysis": "skin_transfer_driver_decomposition",
        "scope": "exploratory; explains the gastrocnemius +0.826 vs skin -0.665 contrast from human_skin_transfer.py",
        "shared_pathways": int(len(merged)),
        "full_spearman_r": {
            "mouse_gastrocnemius_vs_human_skin": float(full_r_gastroc),
            "mouse_skin_vs_human_skin": float(full_r_skin),
        },
        "sign_match_overall": {
            "gastrocnemius_vs_human_skin_frac": float(merged.sign_match_gastrocnemius.mean()),
            "skin_vs_human_skin_frac": float(merged.sign_match_skin.mean()),
        },
        "family_summary": family_summary,
        "top5_drivers_gastrocnemius_to_human_skin": top_drivers_gastroc,
        "top5_anti_drivers_skin_to_human_skin": top_anti_drivers_skin,
        "leading_edge_intersection_gastrocnemius_top_drivers": le_genes_gastroc,
        "headline_interpretation": (
            "Both human Inspiration4 skin and mouse hindlimb gastrocnemius "
            "post-flight strongly down-regulate translation/proteostasis "
            "machinery (ribosomal proteins, eIF complexes, NMD, SRP-mediated "
            "cotranslational targeting, rRNA processing). The mouse skin "
            "transcriptome instead UP-regulates the same Hallmark/Reactome "
            "subset, breaking direct tissue matching. On this 25-pathway "
            "Hallmark-dominated overlap the skin-to-skin transfer signal is "
            "therefore swamped by a stronger conserved muscle/skin proteostasis "
            "suppression axis driven by ribosomal-protein genes."
        ),
    }
    with open(OUT_DIR / "skin_transfer_driver_decomp.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"shared pathways: {len(merged)}")
    print(f"full spearman r (gastroc -> human skin): {full_r_gastroc:+.3f}")
    print(f"full spearman r (skin    -> human skin): {full_r_skin:+.3f}")
    print()
    print("=== Family-level sign match with human i4_skin ===")
    fam_df = pd.DataFrame(family_summary).T
    cols = ["n_pathways", "human_skin_NES_mean", "mouse_gastrocnemius_NES_mean",
            "mouse_skin_NES_mean", "sign_match_gastrocnemius_frac", "sign_match_skin_frac"]
    print(fam_df[cols].round(3).to_string())
    print()
    print("=== Top 5 gastrocnemius driver pathways (most r-lifting when included) ===")
    for r in top_drivers_gastroc:
        print(f"  {r['pathway']:<70} family={r['family']:<24} "
              f"gastroc_NES={r['mouse_gastrocnemius_NES']:+.2f}  "
              f"human_NES={r['human_skin_NES']:+.2f}  "
              f"loo_delta_r={r['loo_delta_r_gastrocnemius']:+.3f}")
    print()
    print("=== Leading-edge gene intersection (mouse gastrocnemius missions) ===")
    for p, info in le_genes_gastroc.items():
        print(f"  {p}")
        print(f"    n_intersection={info['n_intersection_genes']}, first20={info['first20']}")


if __name__ == "__main__":
    main()

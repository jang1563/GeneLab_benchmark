"""v8 Pillar 3 (INTERVENE) — export per-tissue up/down gene signatures.

Builds LINCS/CMap-ready signatures from GeneLab DESeq2 DGE files. For each
tissue, aggregates across missions using sign-consistent ranking (Stouffer's
or simple mean-Wald), then emits top-N up and top-N down human-symbol genes.

Output format (L2S2/SigCom LINCS-compatible):
  v8/intervene/signatures/{tissue}_up.txt
  v8/intervene/signatures/{tissue}_dn.txt
  v8/intervene/signatures/signatures_manifest.json
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = REPO_ROOT / "data" / "mouse"
SYMBOL_MAP = REPO_ROOT / "processed" / "ensembl_symbol_map.csv"
OUT_DIR = Path(__file__).resolve().parent / "signatures"

TISSUES = ["thymus", "gastrocnemius", "skin", "eye", "liver", "kidney"]
TOP_N = 150  # LINCS L1000 uses 978-landmark signatures; 150 up + 150 dn is standard query size


def _find_de_files(tissue: str) -> list[Path]:
    roots = [Path("."), REPO_ROOT]
    seen: set[Path] = set()
    out: list[Path] = []
    for root in roots:
        for pat in (f"data/mouse/{tissue}/*/GLDS-*_differential_expression_GLbulkRNAseq.csv",
                    f"data/mouse/{tissue}/*/GLDS-*_differential_expression.csv"):
            for p in sorted(root.glob(pat)):
                rp = p.resolve()
                if rp in seen:
                    continue
                seen.add(rp)
                out.append(p)
    return out


def _detect_flight_stat(cols: list[str]) -> str | None:
    stat_cols = [c for c in cols if c.startswith("Stat_")]
    candidates = [c for c in stat_cols if re.search(r"Space.?Flight|FLT", c) and re.search(r"Ground.?Control|GC|Vivarium", c)]
    if not candidates:
        return stat_cols[0] if stat_cols else None
    control_first = [c for c in candidates
                     if re.match(r"Stat_\((GC|Ground|VIV|BSL)", c)]
    return control_first[0] if control_first else candidates[0]


def load_tissue_de(tissue: str) -> pd.DataFrame:
    """Return long df: mission, gene (human symbol uppercase), stat."""
    frames = []
    for f in _find_de_files(tissue):
        mission = f.parent.name
        try:
            df = pd.read_csv(f)
        except Exception as e:
            print(f"  SKIP {f}: {e}")
            continue
        stat_col = _detect_flight_stat(list(df.columns))
        if stat_col is None:
            continue
        sym_col = "SYMBOL" if "SYMBOL" in df.columns else None
        if sym_col is None:
            continue
        out = df[[sym_col, stat_col]].rename(columns={sym_col: "gene", stat_col: "stat"})
        out["gene"] = out["gene"].astype(str).str.upper().str.strip()
        out["stat"] = pd.to_numeric(out["stat"], errors="coerce")
        out = out.dropna(subset=["gene", "stat"])
        out = out[out["gene"].ne("") & out["gene"].ne("NAN")]
        out = out.groupby("gene", as_index=False)["stat"].mean()  # dedupe
        # Normalize contrast direction: flip if control-first so positive = up-in-flight
        # detect_flight_stat returned control-first preferentially → positive stat = control > flight
        # We want UP-in-flight; flip sign.
        if re.match(r"Stat_\((GC|Ground|VIV|BSL)", stat_col):
            out["stat"] = -out["stat"]
        out["mission"] = mission
        frames.append(out)
    if not frames:
        return pd.DataFrame(columns=["gene", "stat", "mission"])
    return pd.concat(frames, ignore_index=True)


def aggregate_signature(long: pd.DataFrame, top_n: int = TOP_N) -> tuple[list[str], list[str], pd.DataFrame]:
    """Aggregate per-mission Wald stats into a tissue-level ranked list.

    Uses simple mean z-score across missions (approximate Stouffer's).
    Requires a gene be measured in at least 2 missions if >=2 missions available.
    """
    if long.empty:
        return [], [], pd.DataFrame()
    n_missions = long["mission"].nunique()
    # z-score each mission separately so contrasts are comparable
    long = long.copy()
    long["z"] = long.groupby("mission")["stat"].transform(lambda s: (s - s.mean()) / (s.std() or 1))
    agg = long.groupby("gene").agg(mean_z=("z", "mean"), n=("mission", "nunique")).reset_index()
    if n_missions >= 2:
        agg = agg[agg["n"] >= 2]
    agg = agg.sort_values("mean_z", ascending=False)
    up = agg.head(top_n)["gene"].tolist()
    dn = agg.tail(top_n).sort_values("mean_z")["gene"].tolist()
    return up, dn, agg


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--top-n", type=int, default=TOP_N)
    ap.add_argument("--out-dir", default=str(OUT_DIR))
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    manifest = {"top_n": args.top_n, "tissues": {}}

    for tissue in TISSUES:
        long = load_tissue_de(tissue)
        if long.empty:
            print(f"[{tissue}] NO data")
            continue
        up, dn, agg = aggregate_signature(long, top_n=args.top_n)
        up_path = out / f"{tissue}_up.txt"
        dn_path = out / f"{tissue}_dn.txt"
        ranked_path = out / f"{tissue}_ranked.csv"
        up_path.write_text("\n".join(up) + "\n")
        dn_path.write_text("\n".join(dn) + "\n")
        agg.to_csv(ranked_path, index=False)
        manifest["tissues"][tissue] = {
            "n_missions": int(long["mission"].nunique()),
            "n_genes_ranked": int(len(agg)),
            "n_up_emitted": len(up),
            "n_dn_emitted": len(dn),
            "top_up": up[:5],
            "top_dn": dn[:5],
            "files": {"up": str(up_path.name), "dn": str(dn_path.name), "ranked": str(ranked_path.name)},
        }
        print(f"[{tissue}] missions={long['mission'].nunique()}  "
              f"genes_ranked={len(agg)}  up[:3]={up[:3]}  dn[:3]={dn[:3]}")

    (out / "signatures_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"\nWrote {out}/signatures_manifest.json")


if __name__ == "__main__":
    main()

# DATA_CATALOG.md — GeneLab_benchmark
**자동 생성**: `catalog_datasets.py` — 2026-02-28 07:18 UTC (수동 업데이트: 2026-03-22 — v4 lung/colon 추가)
**PLAN.md 버전**: v1.6 | **Benchmark 버전**: v4.0

> ⚠️ Status legend:
> ✅ Verified (bulk RNA-seq) | ❌ NOT FOUND | ⚠️ Wrong assay type
> 🔁 DUPLICATE (excluded from Category A) | ⛔ EXCLUDED by design

---

## Candidate Datasets (PLAN.md v0.5)

| OSD ID | Tissue | Mission | Track | Status | Assay | Control Types | Bulk RNA-seq Files | Notes |
|---|---|---|---|---|---|---|---|---|
| OSD-48 | liver | RR-1 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC, VC, BC | 23 | primary |
| OSD-137 | liver | RR-3 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 20 | primary |
| OSD-168 | liver | RR-1+3 | 2a | 🔁 DUPLICATE — Category J only | bulk_rnaseq | GC | 54 | DUPLICATE of OSD-48+137 — Category J only (DESIGN_DECISIONS DD-05) |
| OSD-245 | liver | RR-6 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 358 | primary |
| OSD-379 | liver | RR-8 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 290 | primary |
| OSD-242 | liver | RR-9 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 144 | primary |
| OSD-686 | liver | MHU-2 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | unknown | 26 | primary — GLDS-617 RNA-seq files (uG vs GC vs 1G-centrifuge, n=3/group) |
| OSD-617 | liver | MHU-2 | exc | ⛔ EXCLUDED (design) | unknown | unknown | 0 | EXCLUDE — cytology/estrous cycle only, NOT RNA-seq. RNA-seq is in OSD-686 |
| OSD-173 | liver | STS-135 | qc | ⚠️ FAILS QC (n<3/group) | bulk_rnaseq | GC | 14 | C57BL/6CR (C57BL/6 subline), n=2/group — FAILS QC minimum (n<3/group) |
| OSD-25 | liver | STS-135 | exc | ⛔ EXCLUDED (design) | microarray | unknown | 0 | MICROARRAY confirmed — not RNA-seq. Completely excluded |
| OSD-101 | gastrocnemius | RR-1 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 20 | primary |
| OSD-326 | gastrocnemius | RR-9 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | unknown | 24 | primary |
| OSD-401 | gastrocnemius | RR-5 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | unknown | 16 | primary |
| OSD-104 | soleus | RR-1 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 20 | supplementary — check mission count |
| OSD-638 | soleus | MHU-8 | 2a | ❌ NOT FOUND | unknown | ? | 0 | supplementary — unverified |
| OSD-102 | kidney | RR-1 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 20 | primary |
| OSD-163 | kidney | RR-3 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 26 | primary |
| OSD-253 | kidney | RR-7 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 232 | primary |
| OSD-674 | kidney | ? | 2a | ❌ NOT FOUND | unknown | ? | 0 | unverified — mission TBD |
| OSD-4 | thymus | STS-118 | 2a | ⚠️ MICROARRAY (not Category A) | microarray | unknown | 0 | shuttle mission — verify assay type |
| OSD-244 | thymus | RR-6 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 122 | primary |
| OSD-289 | thymus | MHU-2 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 44 | primary |
| OSD-421 | thymus | RR-9 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 29 | primary |
| OSD-238 | skin | MHU-2 (dorsal) | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | AG, GC, VC | 149 | A5 primary — merged with OSD-239 as mission="MHU-2" |
| OSD-239 | skin | MHU-2 (femoral) | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | AG, GC | 142 | A5 primary — merged with OSD-238 as mission="MHU-2" |
| OSD-240 | skin | RR-5 (dorsal) | 2b | ⚠️ TRACK 2b (BAL-TAL strain) | bulk_rnaseq | GC | 20 | BAL-TAL균주 — Track 2a 제외, Track 2b 후보 |
| OSD-241 | skin | RR-5 (femoral) | 2b | ⚠️ TRACK 2b (BAL-TAL strain) | bulk_rnaseq | GC | 19 | BAL-TAL균주 — Track 2a 제외, Track 2b 후보 |
| OSD-243 | skin | RR-6 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | BC, GC | 322 | A5 primary |
| OSD-254 | skin | RR-7 | 2a | ✅ Verified (bulk RNA-seq, C57BL/6J subset) | bulk_rnaseq | GC, VC, BSL, C3H | 80 | A5 primary — C57BL/6J 30개 추출 (BSL+C3H 50개 제외). OSD-254 = mixed strain study |
| OSD-689 | skin | RR-8 | exc | ⛔ EXCLUDED (design) | single_cell | unknown | 50 | EXCLUDE — scRNA-seq (GLDS-689), not bulk RNA-seq |
| OSD-793 | skin | RRRM-1 | exc | ⛔ EXCLUDED (design) | unknown | unknown | 0 | EXCLUDE — targetRNAseq (NanoString panel), not standard bulk mRNA |
| OSD-100 | eye | RR-1 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 14 | primary |
| OSD-194 | eye | RR-3 | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC | 16 | primary |
| OSD-397 | eye/retina | TBD | 2a | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | unknown | 36 | bulk RNA-seq confirmed (n=9 FLT, n=7 GC). Mission TBD — verify ISA metadata |
| OSD-664 | eye | MHU-8 | exc | ⛔ EXCLUDED (design) | unknown | unknown | 0 | EXCLUDE — western blot + immunofluorescence only, NOT RNA-seq |
| OSD-87 | eye | ? | exc | ⛔ EXCLUDED (design) | microarray | unknown | 0 | EXCLUDE — microarray |
| OSD-295 | soleus | HU | HU | 🔬 HU ANALOG (bulk_rnaseq) | bulk_rnaseq | unknown | 44 | Hindlimb Unloading soleus RNA-seq — B6 candidate |
| OSD-335 | liver | HU+HZE | exc | ⛔ EXCLUDED (design) | bulk_rnaseq | unknown | 80 | EXCLUDE: miRNA-Seq (not mRNA bulk RNA-seq) |
| OSD-270 | heart | RR-3 | exc | ⛔ EXCLUDED (design) | spatial | GC | 15 | EXCLUDE from A: Visium spatial transcriptomics |
| OSD-596 | heart | NG-11 | exc | ⛔ EXCLUDED (design) | unknown | unknown | 0 | heart — verify assay type, ≤2 missions |
| **OSD-248** | **lung** | **RR-6** | **2a** | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC, BC | 58 | **v4 추가** — 20 FLT + 19 GC + 19 BC. BC는 v4에서 ground로 포함 |
| **OSD-247** | **colon** | **RR-6** | **2a** | ✅ Verified (bulk RNA-seq) | bulk_rnaseq | GC, BC | 54 | **v4 추가** — 19 FLT + 17 GC + 18 BC. BC는 v4에서 ground로 포함 |

---

## Summary

- ✅ Verified (bulk RNA-seq): **26** (+2 v4: OSD-247 colon, OSD-248 lung; +2 Track 2b: OSD-240/241; +1 subset: OSD-254)
- ❌ Not Found: **2**
- ⚠️ Wrong assay type / Track 2b: **4** (OSD-4 microarray, OSD-240/241 BAL-TAL)
- ⛔/🔁 Excluded by design: **10**

**A5 Skin 최종 구성** (Track 2a, 3-fold LOMO):
- MHU-2: OSD-238 (dorsal) + OSD-239 (femoral) 병합 → 35 binary samples
- RR-6: OSD-243 → 37 binary samples
- RR-7: OSD-254 C57BL/6J 서브셋 → 30 binary samples

---

## Notes on Excluded Studies

### GLDS-168 / OSD-168 — DUPLICATE
OSD-168 is a reprocessed combination of OSD-48 (RR-1) + OSD-137 (RR-3) liver data.
Including in Category A would cause sample duplication across LOMO train/test folds.
**Use only in Category J (J1): pipeline version comparison (v1 vs v2 vs v3).**
Reference: DESIGN_DECISIONS.md DD-05

### OSD-25 (STS-135) — Track 2b only
Uses BALB/c mouse strain, not C57BL/6J. Excluded from Track 2a (primary analysis).
Track 2b (secondary, supplementary): included as cross-strain experiment.
Reference: DESIGN_DECISIONS.md DD-06

### OSD-240/241 (RR-5 Skin) — Track 2b only
BAL-TAL strain (BALB/c × C57BL/6 cross or similar), NOT C57BL/6J.
Excluded from Track 2a (A5). Future Track 2b cross-strain analysis candidate.
Sample IDs include "BAL-TAL" prefix (confirmed from SampleTable condition strings).

### OSD-254 (RR-7 Skin) — C57BL/6J subset extracted
Mixed strain study: C57BL/6J (n=50 incl. BSL) + C3H/HeJ (n=30).
For Track 2a (A5): extracted C57BL/6J non-BSL samples only = 30 samples (10F + 10GC + 10VC).
Time points: 25 days and 75 days post-launch (both included, treated as same RR-7 mission).

### OSD-270 (Heart, Visium) — NOT Category A
Visium spatial transcriptomics, not bulk RNA-seq. Category F (v2.0) only.
Reference: PLAN.md v0.5 Section 3.2

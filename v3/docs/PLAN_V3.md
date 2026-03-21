# GeneLabBench v3 — 종합 구현 계획 (Rev.3, self-review 반영)

## Context

v1.0 (bulk 6-tissue, 25+ tasks, LOMO, ML/FM/LLM) 및 v2.0 (cross-species E1-E3, scRNA-seq F1-F2, temporal T1-T3) 완료.
v3.0은 **다종 비교 확장, 공간전사체, RRRM-2 scRNA-seq, 방사선, 기존 분석 강화**로 확장.

**User 방향**: "이 방향으로 하되 새로운 데이터 추가 (다른 OSDR, 새로운 종, spatial)하고 기존 분석 확장 부분 강화하자"

### 데이터 검증 결과 (Rev.2에서 수정됨)

| 원래 가정 | 실제 확인 | 영향 |
|----------|----------|------|
| GLDS-73/80/109/117 = RNA-seq | ❌ **전부 microarray** (Affymetrix CEL) | Phase 1 전면 재설계 |
| OSD-588 = Drosophila 우주비행 | ❌ **기생봉 감염 지상 실험** | Drosophila → GLDS-207로 교체 |
| GLDS-113 = C. elegans RNA-seq | ❌ **microarray** (n=8) | C. elegans RNA-seq 별도 확인 필요 |
| GLDS-613~628 = GeoMx spatial | ❌ **일반 bulk RNA-seq** | F4 GeoMx 삭제 |
| — | ✅ **RR-10 cerebellum Visium** 발견 (5+5) | 새로운 primary spatial dataset |
| — | ✅ **GLDS-402~405 RRRM-2 scRNA-seq** 발견 | 새로운 Phase 추가 |
| — | ✅ **GLDS-235~249** 13개 신규 bulk RNA-seq | 기존 분석 확장 가능 |

### Task ID 변경 사항 (v1/v2와 충돌 방지)
- v1 H1/H2/H3 = 사전등록 가설 → v3 방사선 태스크는 **R1/R2/R3**로 명명
- v1 A7(Brain)/A8(Heart) = 미구현 보류 → v3에서 **A7(Lung)/A8(Extended Skin)**으로 재할당
- v2 F4 = GeoMx (데이터 검증 실패) → **F4 retired**, F5부터 사용

### HPC RRRM-2 데이터 (User 제보)
- Cayuga에서 RRRM-2 처리 파일 있음 (User 확인)
- 가능 위치: `/athena/masonlab/scratch/projects/GeneLab/ISS_single_cell/.../GLDS-{402,403,404,405}/` (2022 Seurat 객체)
- Phase 0에서 정확한 경로와 파일 형식 확인 필요

---

## Phase 0: 데이터 카탈로그 + 인프라 준비 (즉시 실행)

### 0.1 OSDR API 검증 (최우선)

Phase 1-5 시작 전, 모든 데이터셋 메타데이터 확인:

| 데이터셋 | 확인 사항 | 상태 |
|----------|----------|------|
| GLDS-207 (Drosophila head) | RNA-seq 여부, sample design (FLT vs GC?), n, strain | ⚠️ 웹 확인 완료, API 상세 필요 |
| GLDS-7/37/38/120/208 (Arabidopsis) | RNA-seq 확인, FLT vs GC design, n, ecotype | ⚠️ 일부 확인, 상세 필요 |
| RR-10 cerebellum Visium | OSD ID 확인 (Nature Comm. 2024 data availability에서), 파일 구조 | ⚠️ **최우선 — Phase 2 blocker** |
| GLDS-402~405 (RRRM-2 scRNA) | Cayuga 처리 파일 경로 + 형식 확인 (Seurat→h5ad 변환 필요 여부) | ⚠️ User 제보, 상세 필요 |
| GLDS-235~249 (신규 bulk) | 각 데이터셋별 조직/미션/n 매핑 (13개 전수) | ⚠️ 일부 확인 |
| GLDS-202/211/237 (방사선+HLU) | 2×2 design 확인, 정확한 dose | ⚠️ 웹 확인 완료 |
| OSD-804 (MicroCT) | 파일 형식, feature 추출 가능성 | ⚠️ 미확인 |
| C. elegans RNA-seq | GLDS-113 외 RNA-seq 존재 여부 (OSDR "C. elegans" + "RNA-seq" 검색) | ❓ **최우선** |

### 0.2 Phase 0 Exit Criteria (Go/No-Go)

| Phase | 최소 조건 | No-Go 시 대안 |
|-------|----------|-------------|
| Phase 1 (다종) | ≥3종 (mouse+human+1) RNA-seq FLT vs GC 확보 | Drosophila 실패 시 Arabidopsis만으로 진행 (3종) |
| Phase 2 (공간) | RR-10 Visium OSD ID 확인 + 데이터 접근 가능 | GLDS-270 heart만으로 축소 (exploratory only) |
| Phase 3 (RRRM-2) | GLDS-402~405 중 ≥2 접근 가능 | RRRM-1 추가 분석으로 대체 |
| Phase 4 (기존 확장) | GLDS-235~249 중 ≥3 FLT vs GC RNA-seq 확인 | v1 데이터만으로 cross-tissue 확장 |
| Phase 5 (방사선) | GLDS-202/211/237 중 ≥2 접근 + 2×2 design 확인 | 조건부 Phase → 삭제 |

### 0.3 신규 스크립트
- `v3/scripts/catalog_v3_datasets.py` — OSDR API 자동 검증
- `v3/docs/DATA_CATALOG_V3.md` — 검증 결과 문서화
- `v3/docs/PLAN_V3.md` — 이 계획 저장

### 0.4 Ortholog Framework 구축
- Ensembl BioMart: mouse↔Drosophila, mouse↔Arabidopsis ortholog 다운로드
- C. elegans: RNA-seq 데이터 확인 후 조건부 추가
- `v3/scripts/build_ortholog_db.py`, `v3/scripts/ortholog_utils.py`

---

## Phase 1: Multi-Species Conservation (Category E 확장) — 핵심 Phase

**복잡도**: L | **의존성**: Phase 0 (데이터 검증)
**핵심 질문**: 우주비행 전사체 반응은 종 간에 보존되는가? 진화적 거리와 어떤 관계인가?

### 사용 가능 데이터 (검증됨)

| 종 | 데이터셋 | 조직 | n (추정) | Pathway DB |
|----|---------|------|---------|-----------|
| Mouse (6 tissues) | v1.0 기존 데이터 | liver/thymus/kidney/eye/skin/gastro | 50-200 | Hallmark/KEGG/Reactome |
| Human | JAXA cfRNA (v2 E1) | blood cfRNA | ~20 | Hallmark |
| Drosophila | **GLDS-207** | head | ~30 | **KEGG** (Hallmark 불가) |
| Arabidopsis | **GLDS-7/37/38/120/208** | seedlings/roots | varies | **KEGG** (Hallmark 불가) |
| C. elegans | **미확인** — Phase 0에서 RNA-seq 존재 확인 | whole organism | ? | KEGG |

### Pathway DB 전략 (중요한 설계 결정)

- **Hallmark**: 포유류 전용 → mouse + human만 사용 가능
- **KEGG**: 범용 (mouse, Drosophila, C. elegans, Arabidopsis 모두 지원) → **다종 비교의 primary DB**
- **GO Biological Process**: 범용이나 너무 광범위 → supplementary

→ **E4는 KEGG 기반 NES concordance** (Hallmark은 mouse-human 쌍에만 supplementary로 보고)
→ **(DD-18 후보)**: Cross-kingdom pathway comparison에서 KEGG를 primary DB로 선택 — Hallmark은 포유류 전용, GO BP는 너무 광범위

### 태스크 정의

| Task | 내용 | 방법 | 평가 |
|------|------|------|------|
| **E4a** | Mouse↔Drosophila KEGG NES concordance | fGSEA(KEGG) per species → Spearman r | perm p, bootstrap CI |
| **E4b** | Mouse↔Arabidopsis KEGG NES concordance | 동일 | cross-kingdom 보존 |
| **E4c** | Mouse↔C. elegans KEGG NES concordance | 조건부 (RNA-seq 확인 시) | |
| **E4d** | 전체 종 pairwise concordance matrix | E4a-c + E1(mouse-human) 통합 | heatmap |
| **E5** | Phylogenetic distance vs concordance | 정성적 분석 (4-5 데이터포인트 → 통계적 검정 불충분, 시각화 중심) | scatter plot |

### 신규 스크립트
- `v3/scripts/fetch_multispecies_osdr.py` — Drosophila/Arabidopsis 데이터 다운로드
- `v3/scripts/preprocess_multispecies.py` — 종별 정규화 (genome annotation 상이)
- `v3/scripts/run_multispecies_fgsea.R` — KEGG 기반 fGSEA (종별 KEGG gene set 사용)
- `v3/scripts/e4_multispecies_nes.py` — 다종 NES concordance (v2 E1 패턴 재사용)

### 검증 기준
- [ ] ≥2 종 쌍에서 NES 상관 유의미 (perm p < 0.05)
- [ ] 각 쌍별 ortholog 커버리지 + pathway 커버리지 보고
- [ ] 예상: mouse-Drosophila > mouse-Arabidopsis (동물-식물 경계)

### 산출물
- `v3/processed/orthologs/{mouse_drosophila,mouse_arabidopsis}.tsv`
- `v3/evaluation/E4{a-d}_multispecies_nes.json`, `E5_phylogenetic.json`
- `v3/figures/v3_Fig1_multispecies_concordance.html`

---

## Phase 2: Spatial Transcriptomics (Category F) — 핵심 Phase

**복잡도**: L-XL | **의존성**: Phase 0 (OSD ID 확인)
**핵심 질문**: 공간적 해상도가 우주비행 반응 탐지를 개선하는가?

### 사용 가능 데이터

| 데이터셋 | 조직 | 미션 | 기술 | n | 비고 |
|----------|------|------|------|---|------|
| **RR-10 cerebellum** | Cerebellum | RR-10 | 10x Visium | 5 FLT + 5 GC + basal/vivarium | ★ Primary (Nature Comm. 2024) |
| **GLDS-270** | Heart | RR-3 | 10x Visium | 3 FLT + 2 GC | Exploratory (극소 n) |

- **F4 GeoMx (GLDS-613~628) 삭제**: 실제 bulk RNA-seq이었음, GeoMx 아님

### 태스크 정의

**F3: RR-10 Cerebellum Visium** (Primary)

| Task | 유형 | 내용 |
|------|------|------|
| F3a | Classification | Spot-level FLT vs GC (spatial PCA-LR) |
| F3b | Spatial pattern | Spatially variable genes: Moran's I (FLT vs GC 차이) |
| F3c | Region analysis | 소뇌 해부학적 영역별 (granular layer, Purkinje, molecular layer) 반응 차이 |
| F3d | Cross-resolution | Spot-level AUROC vs v1 bulk brain/cerebellum AUROC 비교 |

**F3e: GLDS-270 Heart Visium** (Exploratory)

| Task | 유형 | 내용 |
|------|------|------|
| F3e | Exploratory | Heart spatial pattern (n=5, descriptive only, 정량적 벤치마크 부적합) |

**F3e 실행 게이트 — ❌ NO-GO (2026-03-21 최종 결정)**:
- OSDR에 raw FASTQ만 존재, SpaceRanger output 없음. 외부 저장소(GEO 등)에도 processed data 없음.
- n=5 (3 FLT + 2 GC) → C(5,2)=10 permutations, min p=0.1. 정량적 분류 통계적으로 불가.
- Brain Visium (n=6)도 AUROC=0.139 (negative result) → heart(n=5)는 더 약할 것으로 예상.
- SpaceRanger 재처리 비용(reference genome + tissue image + GPU 4-6hr) 대비 가치 부족.
- **Future work**: SpaceRanger output이 OSDR에 추가되면 재검토.

### 신규 인프라
- `v3/scripts/spatial_utils.py` — Visium 로딩, Squidpy spatial analysis, Moran's I
- `v3/scripts/process_visium.py` — Space Ranger output → Scanpy/Squidpy 파이프라인
- Conda env: `scanpy`, `squidpy`, `spatialdata` (Cayuga 설치 필요)

### 검증 기준
- [ ] F3a: Spot-level AUROC > 0.6
- [ ] F3b: ≥20 spatially variable genes FLT vs GC 차이 유의 (FDR < 0.05)
- [ ] F3d: Spot-level AUROC vs bulk AUROC 비교 보고

### 산출물
- `v3/evaluation/F3{a-d}_visium_cerebellum.json`, `F3e_visium_heart.json`
- `v3/figures/v3_Fig2_spatial_overview.html`

---

## Phase 3: RRRM-2 scRNA-seq (Category F 확장) — 핵심 Phase

**복잡도**: L | **의존성**: Phase 0 (데이터 확인)
**핵심 질문**: RRRM-2의 scRNA-seq 결과가 RRRM-1 (v2 F2)과 일치하는가? 재현성 검증.

### 사용 가능 데이터

| 데이터셋 | 조직 | 기술 |
|----------|------|------|
| GLDS-402 | Femur bone marrow | 10x Chromium scRNA-seq |
| GLDS-403 | Humerus bone marrow | 10x Chromium scRNA-seq |
| GLDS-404 | PBMCs | 10x Chromium scRNA-seq |
| GLDS-405 | Spleen | 10x Chromium scRNA-seq |

### v2 F2 (RRRM-1)과의 비교

| v2 F2 (RRRM-1) | v3 (RRRM-2) | 비교 가능 |
|-----------------|-------------|----------|
| Blood (PBMCs) | PBMCs (GLDS-404) | ✅ **직접 비교** |
| Muscle | — | ❌ |
| Eye | — | ❌ |
| Skin | — | ❌ |
| — | Femur bone marrow (GLDS-402) | 🆕 새 조직 |
| — | Humerus bone marrow (GLDS-403) | 🆕 새 조직 |
| — | Spleen (GLDS-405) | 🆕 새 조직 |

### 태스크 정의

| Task | 내용 | 방법 |
|------|------|------|
| **F5a** | Cell-type composition per tissue (4 tissues) | Scanpy clustering + annotation |
| **F5b** | Cell-type fGSEA (Hallmark 50) per cell type | Pseudo-bulk → fGSEA (v2 F2-B 패턴) |
| **F5c** | Cell-type LOAO classifier | PCA-LR per cell type (v2 F2-C 패턴) |
| **F5d** | Cross-mission PBMC concordance | RRRM-1 blood vs RRRM-2 PBMC NES Spearman r (v2 F2-D 패턴 확장) |
| **F5e** | Bone marrow cell-type benchmark | Novel: immune progenitor spaceflight detection |

### 신규 스크립트
- `v3/scripts/fetch_rrrm2_scrna.py` — GLDS-402~405 다운로드
- `v3/scripts/process_rrrm2_scrna.py` — Scanpy QC + clustering + annotation (v2 파이프라인 재사용)
- `v3/scripts/f5_rrrm2_benchmark.py` — F5a-e 통합 벤치마크 (v2 F2 코드 패턴)
- `v3/scripts/f5d_cross_mission_scrna.py` — RRRM-1 vs RRRM-2 concordance

### 검증 기준
- [ ] F5c: ≥3 cell types에서 AUROC > 0.7
- [ ] F5d: RRRM-1↔RRRM-2 PBMC NES r > E1 bulk r=0.352 (scRNA-seq 우위 확인)
- [ ] F5e: Bone marrow immune progenitor 분류 가능 여부

### 산출물
- `v3/evaluation/F5{a-e}_rrrm2_scrna.json`
- `v3/figures/v3_Fig3_rrrm2_scrna.html`

---

## Phase 4: 기존 분석 확장 (Strengthen Existing) — 중요

**복잡도**: M-L | **의존성**: 없음 (독립)
**핵심 질문**: v1.0 벤치마크를 새 조직/미션/모델로 확장하면 결론이 일반화되는가?

### 4.1 신규 조직 + 미션 추가 (GLDS-235~249)

확인된 13개 신규 bulk RNA-seq 데이터셋 (일부, Phase 0에서 전수 확인):

| 데이터셋 | 조직 | 미션 | 추가 가치 |
|----------|------|------|----------|
| GLDS-248 | **Lung** | RR-6 | 🆕 v1에 없는 새 조직 |
| GLDS-235 | Liver | RR | 기존 liver fold 확장 |
| GLDS-240/241 | Dorsal/Femoral Skin | RR-5 | skin fold 확장 |
| GLDS-243 | Dorsal Skin | RR-6 | skin fold 확장 |
| GLDS-249 | Muscle | 치료 개입 | muscle fold 확장 + therapeutic |

**주의**: GLDS-237 (skin)은 **방사선+HLU 2×2 design** → Phase 5 전용. Phase 4 A8에 포함하지 않음.

### 태스크 정의 (v1 A7=Brain, A8=Heart 미구현 → v3에서 재할당)

| Task | 내용 | 평가 | 상세 |
|------|------|------|------|
| **A7** | Lung spaceflight detection | AUROC, bootstrap CI, perm p | GLDS-248 기반, v1 6-tissue LOMO와 동일 방법론 |
| **A8** | Extended skin LOMO (RR-5 fold 추가) | AUROC delta vs v1 skin | GLDS-240/241/243 + 기존 skin 통합 |
| **B_ext** | 7-tissue cross-mission transfer matrix | Transfer AUROC 7×7 matrix | v1 C (4-pair) → v3 전체 매트릭스 확장 (lung 포함) |
| **A_held** | 추가 held-out validation | AUROC, perm p | 신규 미션 데이터를 held-out test로 사용 (v1 A4 thymus, A5 skin 패턴) |

### 4.2 추가 Foundation Model 평가

| 모델 | 유형 | v1/v2 결과 | v3 추가 |
|------|------|-----------|---------|
| PCA-LR | Classical | 0.758 (baseline) | 신규 조직에도 적용 |
| Geneformer | Gene FM | 0.476 | — (v1 완료) |
| scGPT | Gene FM | 0.666 | — (v1.3 완료) |
| **scFoundation** | Gene FM | 미평가 | 🆕 최신 모델 (2024) |
| **UCE** | Gene FM | 미평가 | 🆕 Universal Cell Embeddings |

**운영 상태 (2026-03-20, FINAL)**:
- `B_ext`: ✅ **COMPLETE**. 7×7 transfer matrix (42 pairs), Method A + C, bootstrap CI + permutation p. `B_ext_transfer_matrix.json` (43KB).
- `UCE` (7 tissues, seeded): ✅ **COMPLETE** (job 2707938, `np.random.seed(42)` 적용). 결과:
  - liver 0.459 (p=0.83), kidney 0.489 (p=0.57), thymus **0.632** (p=0.031), gastro 0.578 (p=0.22), eye 0.550 (p=0.30)
  - lung 0.555 (5-fold CV), colon 0.449 (5-fold CV)
  - UCE 비결정적 문제 해결: `eval_data.py`의 `np.random.choice`/`shuffle`이 원인. Seed 패치 후 안정.
- `scFoundation` (7 tissues, full): ✅ **COMPLETE** (job 2707916, CI + perm test 포함). 결과:
  - liver **0.635** [0.548, 0.715] (p=0.001), kidney 0.541 [0.419, 0.658] (p=0.25)
  - thymus 0.487 [0.337, 0.627] (p=0.58), gastro **0.691** [0.498, 0.867] (p=0.035)
  - eye 0.563 [0.375, 0.750] (p=0.26), lung 0.389 (5-fold), colon **0.755** (5-fold)
- **FM 결론**: PCA-LR (liver 0.758) 대비 모든 FM이 열등. scFoundation이 UCE보다 전반적으로 우수. 유의미: scFoundation liver (p=0.001), gastro (p=0.035); UCE thymus (p=0.031). FM은 spaceflight detection에 PCA-LR 대비 이점 없음 → 사전학습 cell atlas 지식이 spaceflight perturbation 탐지에 도움 안 됨.
- 상세 실행 순서와 smoke-test 절차는 `v3/docs/HPC_FM_RUNBOOK.md` 참조.

### 4.3 Cross-Tissue Transfer 확장

v1 Category C: 4 pairs × 3 methods → v3: **전체 7×7 transfer matrix** (directional 42 pairs) + 신규 조직(lung, colon)

### 신규 스크립트
- `v3/scripts/fetch_new_bulk_osdr.py` — GLDS-235~249 다운로드 + 전처리
- `v3/scripts/extend_v1_benchmark.py` — A7/A8/B_ext 태스크 (v1 패턴 재사용)
- `v3/scripts/run_new_fm.py` — scFoundation/UCE fine-tuning (Phase 0에서 모델 가용성 확인)

### 검증 기준
- [x] A7 (Lung): AUROC 0.608 (목표 0.65 미달이나 chance 이상, 의미있는 탐지)
- [x] B_ext: 42 pairs 완료, Method A 평균 ~0.57, v1 4-pair 결과와 비교 가능
- [x] ≥1 new FM이 Geneformer (0.476)보다 개선: **scFoundation gastro 0.691 (p=0.035)**, liver 0.635 (p=0.001), colon 0.755. UCE thymus 0.632 (p=0.031)
- [x] scFoundation full 7 tissues (CI+perm) 완료: job 2707916
- [x] UCE deterministic re-run (seeded, job 2707938) 완료: 비결정적 문제 해결

### 산출물
- `v3/evaluation/A7_lung.json`, `A8_extended_skin.json`, `B_ext_transfer_matrix.json`
- `v3/evaluation/FM_scfoundation.json`, `FM_uce.json`

---

## Phase 5: Radiation × Microgravity Interaction (Category H 재설계) — 조건부

**복잡도**: M | **의존성**: Phase 0 (microarray 처리 결정)

### 문제와 재설계

**원래 계획**: GLDS-73/80/109/117로 dose-response regression → **불가** (전부 microarray)

**사용 가능 RNA-seq 방사선 데이터**:

| 데이터셋 | 조직 | Design | 방사선 | Dose |
|----------|------|--------|--------|------|
| GLDS-202 | Retina + Brain | 2×2 (LDR × HLU) | 0.4 Gy | 단일 |
| GLDS-211 | Spleen | 2×2 (LDR × HLU) | chronic low-dose | 단일 |
| GLDS-237 | Skin | 2×2 (LDR × HLU) | chronic low-dose | 단일 |

→ **Dose-response regression 불가** (단일 선량만). **2×2 interaction analysis**로 재설계:

### 태스크 정의 (재설계, v1 H1-H3 가설 ID와 구분하여 R1-R3으로 명명)

| Task | 내용 | 방법 |
|------|------|------|
| **R1** | 방사선 vs HLU 효과 분리 | 2×2 ANOVA: main effects + interaction |
| **R2** | ISS flight NES vs ground radiation NES concordance | fGSEA NES Spearman r (GLDS-202/211/237 각각) |
| **R3** | Radiation + HLU combined vs ISS flight | 2×2 combined effect가 ISS와 더 유사한가? |

### 선택지: Microarray 포함 여부

| Option | 장점 | 단점 |
|--------|------|------|
| **A: RNA-seq만 (GLDS-202/211/237)** | DD-01 일관성, 기존 파이프라인 재사용 | dose-response 불가, 3개 조직만 |
| **B: Microarray도 포함 (GLDS-73/80/109/117)** | dose-response regression 가능, 다양한 선종 | 새 전처리 파이프라인 필요 (RMA), DD-01 확장 |

→ **User 결정 필요**: 이 선택지를 Phase 0에서 user에게 확인

### 검증 기준
- [ ] R1: Radiation main effect가 HLU main effect와 분리 가능한가?
- [ ] R2: ISS flight NES와 ground radiation NES 간 상관 보고
- [ ] (Option B 시) ≥2 조직에서 dose-response Spearman ρ > 0.5

### 동기: v2 T4에서 ISS 선량 범위(6.45-8.05 mGy)가 dose-response 분석에 불충분함을 확인 → 지상 방사선 데이터 필요

---

## Phase 6: Multimodal Integration (Category I) — stretch goal

**복잡도**: M | **의존성**: Phase 5 (방사선 metadata)

### 태스크

| Task | 모달리티 | Baseline 비교 |
|------|---------|--------------|
| I1 | Gene + body/organ weight | vs A1 gene-only (0.758) |
| I2 | Gene + MicroCT (OSD-804) | vs A1 (bone 조직 한정) |

- **I3 (RadLab+EDA) 축소**: Phase 5 결과에 따라 결정

**I1 실행 게이트 — ❌ NO-GO (2026-03-21 최종 결정)**:
- OSDR SampleTable CSV에 body weight 컬럼 없음. OSDR API에도 phenotype/weight metadata 없음.
- Publication supplementary table에서 수동 큐레이션 필요 → 노동 대비 가치 낮음.
- 매칭 가능한 mission 수 불확실 → LOMO 신뢰성 보장 불가.
- **Future work**: OSDR에 phenotype metadata가 추가되거나, 논문별 수동 큐레이션이 완료되면 재검토.

### 검증 기준
- [ ] ≥1 fusion 전략에서 AUROC > gene-only + 2%p

---

## Stretch Goals (핵심 Phase 완료 후)

| 항목 | Phase | 이유 |
|------|-------|------|
| J6: STAR vs HISAT2 | 7 | 300 GB FASTQ 필요, 낮은 novelty |
| J7: mm10 vs mm39 | 7 | J6와 동시 실행 |
| G1-G6: Microbiome | 8 | 완전히 새로운 도메인 (QIIME2), 높은 인프라 비용 |

---

## 실행 순서 (Dependency Graph)

```
Phase 0 (데이터 검증 + ortholog)
    │
    ├──→ Phase 1 (Multi-Species E4/E5) ──────┐
    │                                          │
    ├──→ Phase 2 (Spatial F3)  ───────────────┤
    │                                          ├──→ 통합 Figure + 논문
    ├──→ Phase 3 (RRRM-2 scRNA F5) ──────────┤
    │                                          │
    ├──→ Phase 4 (기존 확장 A7/A8/FM) ────────┤
    │                                          │
    └──→ Phase 5 (Radiation H, 조건부) ────────┘
              │
              └──→ Phase 6 (Multimodal I, stretch)
```

**Phase 1, 2, 3, 4 모두 Phase 0 완료 후 병렬 시작 가능**

---

## v3 통합 Figure (4 main, 기존 컨벤션 유지)

**명명**: `v3/figures/v3_Fig{N}_*.html` (v1/v2와 구분)
**스타일**: D3.js v7, Okabe-Ito palette, Nature Methods typography (12/10/9/7.5px), self-contained HTML+SVG

| Figure | 내용 | Phase |
|--------|------|-------|
| v3_Fig1 | Multi-species KEGG NES concordance matrix + phylogenetic distance | Phase 1 |
| v3_Fig2 | Spatial Visium: cerebellum spot-level classification + region analysis | Phase 2 |
| v3_Fig3 | RRRM-2 scRNA: cell-type benchmarks + RRRM-1 vs RRRM-2 concordance | Phase 3 |
| v3_Fig4 | 기존 확장: 신규 조직/FM 비교 + multimodal (있을 경우) | Phase 4+6 |

---

## Evaluation JSON 스키마 호환성

- **분류 태스크** (E4, F3a, F5c, A7 등): 기존 스키마 유지 (`auroc`, `ci_lower`, `perm_p`)
- **NES concordance** (E4a-d, F5d, R2): v2 E1 스키마 유지 (`spearman_r`, `bootstrap_ci`, `perm_p`)
- **Regression** (dose-response, Option B 선택 시만): 새 스키마 추가 (`spearman_rho`, `r_squared`, `rmse`)
- **Spatial** (F3b): 새 필드 (`n_spatially_variable_genes`, `morans_i_stats`)

---

## 리스크 평가 (업데이트)

| 리스크 | 확률 | 영향 | 완화 |
|--------|------|------|------|
| GLDS-207 Drosophila n < 6 또는 design 부적합 | 중 | Phase 1 약화 | Arabidopsis 5+ datasets로 보완 |
| RR-10 Visium OSD ID 불확실 | 중 | Phase 2 지연 | Nature Comm. 2024 corresponding author에게 확인 |
| RRRM-2 데이터 비공개/접근 불가 | 저 | Phase 3 차단 | GLDS-402~405 공개 확인됨 |
| scFoundation/UCE Cayuga 설치 실패 | 중 | Phase 4 FM 부분 약화 | PCA-LR baseline으로 신규 조직 평가 (핵심은 데이터) |
| Arabidopsis KEGG coverage < 50 pathways | 중 | E4b 약화 | GO BP로 대체 |
| Squidpy/Space Ranger Cayuga 호환성 | 중 | Phase 2 지연 | scanpy 기반 대안 (squidpy 없이 기본 분석) |

---

## HPC 리소스 추정 (핵심 Phase만)

| Phase | CPU 시간 | GPU 시간 | 저장공간 |
|-------|---------|---------|---------|
| Phase 0 (검증+ortholog) | ~5h | 0 | ~2 GB |
| Phase 1 (다종) | ~30h | 0 | ~15 GB |
| Phase 2 (공간) | ~50h | 0 | ~50 GB |
| Phase 3 (RRRM-2 scRNA) | ~40h | 0 | ~30 GB |
| Phase 4 (기존 확장) | ~30h | ~10h (FM) | ~20 GB |
| Phase 5 (방사선) | ~20h | 0 | ~10 GB |
| **합계 (핵심)** | **~175h** | **~10h** | **~127 GB** |

---

## 과학적 서사 (Scientific Narrative)

v3의 통합 질문: **"GeneLabBench: 종, 해상도, 모달리티를 넘어 우주비행 전사체 벤치마크의 일반화 검증"**

1. **v1**: 한 종(mouse), 한 기술(bulk RNA-seq), 한 시점 → 기초 벤치마크
2. **v2**: 두 종(mouse+human), 시간(temporal), 세포 수준(scRNA-seq) → 차원 확장
3. **v3**: 다종(4종), 공간(Visium), 재현성(RRRM-2), 기존 강화(신규 조직+FM) → **일반화 검증**

핵심 논점:
- 우주비행 반응이 **kingdom boundary**를 넘어 보존되는가? (E4/E5)
- **공간적 위치**가 세포 수준 이상의 정보를 제공하는가? (F3)
- v2 scRNA 결과가 **독립 미션에서 재현**되는가? (F5d)
- v1 결론이 **새 조직/FM에서도 유지**되는가? (Phase 4)

### 음성 결과 해석 (Scientific rigor)

| Phase | 음성 결과 | 해석 |
|-------|----------|------|
| Phase 1 | Mouse-Drosophila concordance ≈ 0 | 우주비행 반응은 포유류 특이적 |
| Phase 2 | Spot-level AUROC < bulk AUROC | 공간적 노이즈가 해상도 이점을 상쇄 |
| Phase 3 | RRRM-1↔RRRM-2 PBMC NES r < 0.3 | v2 결과가 미션 특이적 (재현 실패) |
| Phase 4 | 신규 조직 AUROC < 0.6 | 특정 조직은 우주비행에 덜 반응 |
| Phase 5 | Radiation NES ≉ ISS flight NES | ISS 반응은 방사선 이외 요인이 지배 |

---

## 수정/생성 파일 목록

### Phase 0 (즉시)
| 파일 | 작업 |
|------|------|
| `v3/docs/PLAN_V3.md` | CREATE |
| `v3/docs/DATA_CATALOG_V3.md` | CREATE |
| `v3/scripts/catalog_v3_datasets.py` | CREATE |
| `v3/scripts/build_ortholog_db.py` | CREATE |
| `v3/scripts/ortholog_utils.py` | CREATE |

### Phase 1-5 (각 phase 별 스크립트는 위 섹션 참조)

### 재사용 기존 코드
| 파일 | 재사용 부분 |
|------|-----------|
| `scripts/utils.py` | load_metadata, align_features_with_meta |
| `scripts/evaluate_submission.py` | bootstrap AUROC, permutation p |
| `scripts/generate_tasks.py` | LOMO split 패턴 |
| `v2/scripts/cross_species_nes_comparison.py` | NES concordance 패턴 |
| `v2/scripts/rrrm1_f2*.py` | scRNA-seq 벤치마크 패턴 |

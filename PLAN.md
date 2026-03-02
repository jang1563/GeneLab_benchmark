# GeneLab_benchmark: Comprehensive Space Biology Multi-Modal Benchmark
**최종 수정**: 2026-02-28
**버전**: 0.6 (Pathway Analysis Integration)

---

## 1. 비전 및 목표

### 핵심 비전
NASA OSDR/GeneLab의 마우스 다조직 RNA-seq 데이터를 출발점으로,
우주비행 전사체 서명의 **미션 간·조직 간 일반화 능력**을 체계적으로 평가하는
벤치마크를 단계적으로 구축한다.

### 이중 목적
1. **AI/ML 벤치마크**: 알고리즘의 우주비행 오믹스 학습 및 일반화 능력 비교
2. **파이프라인 표준화 벤치마크**: 분석 파이프라인 간 재현성·일관성 정량 평가

### 협력 체계
- NASA GeneLab AI/ML AWG 프로그램 담당자와 공동연구 진행 중
- 내부 데이터 접근 및 기존 AWG 분석과의 정합성 확보

### Expected Findings (사전 등록 가설)

코드 착수 전 확정된 검증 가설. 데이터 분석은 이 가설을 검증 또는 반증하기 위해 수행.

```
H1 (Primary): 마우스 우주비행 전사체 서명의 미션 간 일반화 능력은 조직별로 유의미하게 다르다.
  간(Liver)은 가장 일관된 서명을 보이고, 근육(Muscle)은 미션 특이적 반응이 지배적이다.
  → 검증: B1 LOMO AUROC matrix에서 Liver 평균 > Muscle 평균 (paired Wilcoxon test)

H2 (Secondary): ML 모델의 cross-mission 전이 실패는 배치효과보다 생물학적 다양성에 기인한다.
  ComBat-seq 보정 후 B1 transfer AUROC 개선 < 5%p.
  → 검증: J3 배치보정 비교에서 보정 전후 B1 AUROC 차이 측정

H3 (Tertiary): 우주비행 전사체 반응은 유전자 수준보다 경로(pathway) 수준에서 더 잘 보존된다.
  Category C: 방법 A(gene-level) < 방법 B(DEG overlap) << 방법 C(pathway-level).
  → 검증: C1~C4 각 task에서 3가지 방법 AUROC 비교
```

### 차별화 포인트 (vs. 기존 벤치마크)

| 축 | SpaceOmicsBench | BPS Benchmark | SOMA Atlas | **이 프로젝트 v1** |
|---|---|---|---|---|
| 종 | 인간 | 마우스 | 인간 | **마우스 (주)** |
| 조직 | 혈액/PBMC/피부 | 간 only | 혈액/다조직 | **8개 조직 × 10+ 미션** |
| 데이터 유형 | 오믹스 | 오믹스 | 오믹스 | **bulk RNA-seq + 미션 메타데이터** |
| ML Task 수 | 21개 | 없음 | 없음 | **~25개** |
| 크로스미션 일반화 | 없음 | 없음 | 없음 | **LOMO full cross-table** |
| 크로스티슈 전이 | 없음 | 없음 | 없음 | **체계적 transfer benchmark** |
| 파이프라인 비교 | 없음 | 없음 | 없음 | **DESeq2 vs edgeR vs limma** |

---

## 2. 모듈식 출판 전략 (버전별 범위)

현재 계획에서 A~J 10개 카테고리는 세 가지 독립적인 과학적 질문에 해당한다.
**각 버전은 독립적으로 출판 가능한 단위**로 설계하고, 인프라와 코드는 누적적으로 재활용한다.

```
┌──────────────────────────────────────────────────────────────┐
│  v1.0  질문: "마우스 우주비행 전사체 서명이                    │
│              미션 간·조직 간 일반화되는가?"                    │
│  범위: Category A + B + C + D (bulk RNA-seq, 마우스)          │
│  Task: ~25개 / GLDS: ~25개 / 조직: 8개 / 미션: 10+           │
│  저널: Genome Biology                                         │
│  기간: 6~9개월                                                │
├──────────────────────────────────────────────────────────────┤
│  v2.0  질문: "우주비행 반응이 종·세포·마이크로바이옴            │
│              수준에서 보존되는가?"                             │
│  범위: Category E + F + G (다종, 단일세포, 미생물)            │
│  저널: Cell Systems / npj Systems Biology                     │
│  전제: v1.0 출판 후, C. elegans 데이터 존재 확인 후           │
├──────────────────────────────────────────────────────────────┤
│  v3.0  질문: "물리적 우주 환경이 오믹스 프로파일과             │
│              어떻게 연결되는가?"                               │
│  범위: Category H (재설계) + I + J                            │
│  저널: Nature Communications                                  │
│  전제: 지상 방사선 실험 OSDR 데이터 확보 후                    │
└──────────────────────────────────────────────────────────────┘
```

> **현재 문서의 초점**: v1.0 설계 완성 및 실행.
> v2.0·v3.0은 섹션 11 "Future Extensions"에 개요만 기록.

---

## 3. 핵심 설계 원칙 (리뷰 반영)

### 3.1 대조군 유형 명시 (Critical)

OSDR 연구에는 세 종류의 지상 대조군이 존재하며, 연구마다 구성이 다르다.
모든 데이터셋에 `control_type` 필드를 반드시 기록하고, 분석을 두 트랙으로 분리한다.

| 대조군 유형 | 정의 | 측정 대상 |
|---|---|---|
| **Ground Control (GC)** | 지구에서 동일 AEM 하드웨어 사용 | 미세중력 효과만 (순수) |
| **Vivarium Control (VC)** | 일반 사육장 케이지 | 미세중력 + 하드웨어 효과 합산 |
| **Basal Control (BC)** | 발사 전 기준값 | 발사 스트레스 포함 |

- **Track 1**: Flight vs GC → 진짜 미세중력 효과 (GC 있는 연구만)
- **Track 2**: Flight vs VC → 전체 우주비행 효과 (VC 있는 연구 전체)
- Track 1이 더 엄밀하지만 샘플 수 적음 → 두 트랙 결과를 나란히 보고

### 3.2 Category A = Bulk RNA-seq 전용

Visium (spatial), scRNA-seq, microarray는 Category A에서 제외.
이 assay 유형들은 각각 Category F (단일세포/공간) 또는 Category J (파이프라인)에만 포함.

```
포함: STAR-aligned, DESeq2 normalized counts (Illumina bulk RNA-seq)
제외: 10X Visium, NanoString GeoMx, scRNA-seq, snRNA-seq, Affymetrix/Agilent microarray
```

### 3.3 Task 포함 최소 기준 (샘플 수)

| 기준 | 값 | 근거 |
|---|---|---|
| 최소 미션 수 | ≥ 3 | LOMO 의미 있으려면 |
| 최소 총 샘플 수 | ≥ 50 | 학습 세트 30+ 확보 |
| 최소 그룹 반복 수 | ≥ 3/그룹/미션 | 통계적 최소 요건 |

→ A8 (Heart, 2 missions)는 supplementary로 이동.

### 3.4 크로스스터디 Feature Representation 정책

다른 GLDS 연구 간의 DESeq2 normalized counts는 직접 비교 불가.
용도에 따라 representation을 달리 한다:

| Task 유형 | Feature Representation | 근거 |
|---|---|---|
| A (단일 조직, 이진 분류) | **log2(DESeq2 normalized counts), per-sample** | ⚠️ LFC 사용 금지 (label leakage: LFC 계산에 flight/ground 레이블 필요) |
| B (크로스미션) | **LFC** (apeglm shrinkage, 학습 미션에서만 계산) | scale-invariant, mission-내 정규화. 테스트 미션 LFC 사전 노출 금지 |
| C (크로스티슈) | Pathway score 또는 LFC (방법별 상이) | 3가지 방법 병렬 비교 (섹션 5 Category C 참조) |
| D (조건 예측) | log2(normalized counts) 또는 LFC | task별 선택, 명시 필수 |
| J (파이프라인 비교) | 원본 counts 사용 | 파이프라인 차이 측정 목적 |

**LFC 계산 방법 확정** (Category B, C에만 적용):
```
Primary:    DESeq2 lfcShrink(type="apeglm")  — n=3/group 소규모에서 분산 안정화
Comparison: DESeq2 lfcShrink(type="ashr")   — Category J 파이프라인 비교용
Forbidden:  raw LFC (shrinkage 없음) — n<5 환경에서 분산 극심
```

### 3.5 Feature Selection 정책

전체 유전자(~22,000개)를 그대로 쓰면 차원의 저주 발생.
다음 순서로 필터링:

1. **Low-expression filter**: 전체 샘플의 20% 미만에서 발현(count > 1)인 유전자 제거
2. **Variance filter**: 하위 25% variance 유전자 추가 제거
3. 최종 feature 수: ~5,000~8,000개 (조직별 상이)
4. Baseline에서 이 필터 사용, 원본(22k) vs 필터링(5~8k) 비교를 J 카테고리에 포함

> ⚠️ **LOMO Leakage 주의사항 (P0)**:
> Feature selection은 반드시 **LOMO 루프 내부**에서 수행. 전체 데이터 기준 필터링은
> holdout 미션 정보가 variance 계산에 노출되어 성능 과대 추정을 유발.
> ```python
> # 잘못된 방법 (현재 구현 주의):
> filtered_genes = variance_filter(ALL_missions_including_holdout)
>
> # 올바른 방법 — generate_tasks.py에서 split_aware_feature_filter()로 구현:
> for train_missions, test_mission in lomo_splits:
>     filtered_genes = variance_filter(train_missions_only)
>     train_X = train_data[filtered_genes]
>     test_X  = test_data[filtered_genes]  # 같은 유전자 세트 적용
> ```

### 3.6 근육 유형 분리

RR-1의 GLDS-99 (EDL), 101 (gastrocnemius), 103 (quadriceps), 104 (soleus)는
**근육 유형별 별도 task** 또는 **근육 유형 다중 클래스 task**로 처리.
Soleus(서행근)와 EDL(속행근)의 spaceflight 반응이 반대 방향일 수 있음.

### 3.7 지상 아날로그(Hindlimb Unloading) 포함

HU(Hindlimb Unloading)는 미세중력 지상 아날로그:
- HU = 미세중력 효과만 (방사선 없음)
- 실제 ISS = 미세중력 + 방사선 + 격리 + 스트레스

OSDR에 HU RNA-seq 연구 다수 존재. Phase 1 데이터 카탈로그 시 수집.
**B 카테고리 확장**: ISS flight vs HU 비교 task 추가 → "진짜 우주비행 특이 서명" 식별.

### 3.8 마우스 계통 처리 정책

OSDR 연구에서 C57BL/6J (대부분 ISS 미션) vs BALB/c (STS-135 등 일부 셔틀 미션) 혼용.
계통 차이는 spaceflight 신호보다 클 수 있어 명시적 분리 필요.

```
Track 2a (Primary): C57BL/6J 미션만 포함
  → Category A, B, C의 주 분석 (논문 main figures)
  → B1 6×6 matrix = C57BL/6J 미션 6개 (GLDS-48, 137, 245, 379, 242, 617)

Track 2b (Secondary): 모든 계통 포함 (C57BL/6J + BALB/c)
  → Category A, B, C의 supplementary 분석
  → BALB/c 포함 task는 "cross-strain" 실험으로 명시적 reframe

D4 (strain prediction AUROC)로 계통 효과 크기 정량화:
  → D4 AUROC가 높을수록 계통이 강한 confounder임을 투명하게 보고
```

각 GLDS의 `mouse_strain` 필드는 `missions.json`의 필수 항목.

### 3.9 Task 난이도 계층화

SpaceOmicsBench 방식과 호환:

| 레벨 | 기준 | 예시 |
|---|---|---|
| **Calibration** | Baseline AUROC > 0.85, 검증 목적 | A1 Liver (많은 미션, 강한 신호) |
| **Standard** | Baseline AUROC 0.65~0.85 | A3 Kidney, D2 Mission ID |
| **Advanced** | Baseline AUROC 0.55~0.65, 의미 있는 신호 | B1 Cross-mission, C1 Cross-tissue |
| **Frontier** | Baseline ≈ random, 학습 가능성 경계 | C4 Universal signature, B3 BION |

---

## 4. 데이터 전략 (v1.0 범위)

### 4.1 마우스 Bulk RNA-seq — 조직 × 미션 매트릭스

Category A 포함 가능 조직 (bulk RNA-seq 확인된 것만, GLDS ID 검증 필요):

| 조직 | RR-1 | RR-3 | RR-6 | RR-9 | MHU-2 | 기타 | 미션 수 | 포함 여부 |
|---|---|---|---|---|---|---|---|---|
| **간 (Liver)** | OSD-48 | OSD-137 | OSD-245 | OSD-242 | **OSD-686** (GLDS-617 RNA-seq) | OSD-379(RR-8) | **6** (Track 2a, C57BL/6J) | **A1 ✓** |
| **Soleus** | OSD-104 | - | - | - | - | ~~GLDS-638(MHU-8)~~ (❌ Not Found) | 1 | **제외 → A7 이후** |
| **Gastrocnemius** | OSD-101 | - | - | OSD-326 | - | OSD-401(RR-5) | 3 | **A2 ✓** |
| **Quadriceps** | GLDS-103 | - | - | - | - | - | 1 | **제외** |
| **EDL** | GLDS-99 | - | - | - | - | - | 1 | **제외** |
| **신장 (Kidney)** | OSD-102 | OSD-163 | - | - | - | OSD-253(RR-7), ~~OSD-674~~ (❌) | **3** | **A3 ✓** (ギリギリ 최소) |
| **흉선 (Thymus)** | - | - | OSD-244 | OSD-421 | OSD-289 | ~~OSD-4~~ (⚠️ microarray) | **3** | **A4 ✓** |
| **피부 (Skin)** | - | - | OSD-243 | - | OSD-238, OSD-239 | ~~OSD-689~~ (⚠️ scRNA) | **2 missions** | **⚠️ 미션 수 부족** |
| **눈/망막 (Eye)** | OSD-100 | OSD-194 | - | - | - | OSD-664 (⚠️ unknown) | **2+?** | **⚠️ 확인 필요** |
| **뇌 (Brain)** | - | - | - | - | - | GLDS-33(STS-128 micro) | 확인 필요 | **보류** |
| **심장 (Heart)** | - | GLDS-596 (NG-11) | - | - | - | - | ≤2 | **Supplementary** |
| **부신 (Adrenal)** | - | GLDS-161 | - | - | - | - | 1 | **제외** |
| **폐 (Lung)** | - | - | GLDS-248 | - | - | - | 1 | **제외** |
| **대장 (Colon)** | - | - | GLDS-247 | - | - | - | 1 | **제외** |

> **첫 번째 작업**: `catalog_datasets.py`로 위 GLDS ID 실제 존재·파일 형식 전수 검증.
> Visium(GLDS-270), scRNA(GLDS-748 등)는 Category A에서 완전 제외.
>
> ✅ **catalog_datasets.py 검증 결과 (2026-02-28)** — 주요 발견:
>
> **OSD-617 수정**: OSD-617은 세포학(에스트러스 사이클) 데이터만 포함. MHU-2 간 RNA-seq는
> **OSD-686**에 있음 (GLDS-617 파일명 사용). 조건: uG(n=3), GC(n=3), 인공중력(n=3).
> Category A1: uG vs GC 비교 사용. 인공중력 그룹 → D6 (인공중력 효과) task로 활용.
>
> **STS-135 간 RNA-seq**: OSD-173 확인. C57BL/6CR 계통, **n=2/group** → QC 탈락 (최소 n=3).
> OSD-25 (STS-135) = Microarray → 완전 제외. STS-135 간 데이터 없음.
>
> **A1 Liver 확정 미션 (6개, Track 2a C57BL/6J)**:
> OSD-48(RR-1), OSD-137(RR-3), OSD-245(RR-6), OSD-379(RR-8), OSD-242(RR-9), OSD-686(MHU-2)
> → B1 6×6 matrix 가능 ✓
>
> ⚠️ **Skin (A5) 위기**: OSD-689(RR-8) = scRNA-seq (❌), OSD-238+239 = MHU-2 (2개가 같은 미션).
> 실제로 2개 미션만 확보 (MHU-2 + RR-6) → 추가 미션 탐색 필요.
>
> ⚠️ **Eye (A6) 확인 필요**: OSD-664(MHU-8) = assay type unknown. 2 confirmed missions만.
> 3개 미션 미충족 시 A6 supplementary로 강등.
>
> ⚠️ **GLDS-168 제외 (중복 샘플)**: GLDS-168은 GLDS-48(RR-1) + GLDS-137(RR-3) 재처리.
> **Category A에서 완전 제외. Category J (J1)에만 사용.**

### 4.2 미션 메타데이터 표준 스키마

```json
{
  "mission_id": "RR-1",
  "vehicle": "SpaceX CRS-4",
  "launch_date": "2014-09-21",
  "duration_days": 37,
  "hardware": "AEM",
  "mouse_strain": "C57BL/6J",
  "mouse_sex": "Female",
  "mouse_age_weeks": 10,
  "n_flight": 10,
  "n_ground_control": 10,
  "n_vivarium_control": 10,
  "control_types_available": ["GC", "VC"],
  "radiation_total_gy": 0.06,
  "radiation_source": "RadLab_DOSTEL1",
  "avg_temperature_c": 26.2,
  "avg_co2_ppm": null,
  "altitude_km": 400,
  "inclination_deg": 51.6,
  "has_hu_analog": false,
  "genelab_studies": {
    "liver": "GLDS-48",
    "gastrocnemius": "GLDS-101",
    "kidney": "GLDS-102",
    "eye": "GLDS-100"
  }
}
```

### 4.3 비(非)오믹스 데이터 — v1.0에서의 역할

v1.0에서 physiology/ALSDA 데이터는 **ML task 입력이 아닌 보조 역할**로만 사용:

| 데이터 유형 | v1.0 역할 | 대표 OSD |
|---|---|---|
| 체중, 장기 무게 | 샘플 품질 확인, 이상값 필터 | 거의 모든 RR |
| 방사선량 (RadLab) | 미션 메타데이터 covariate | RR-1~23 |
| 환경 텔레메트리 (EDA) | 미션 메타데이터 covariate | RR-8 이후 |
| MicroCT (OSD-804) | 생물학적 검증 참고 | RR-1 |
| 보행 데이터 (OSD-478) | 생물학적 검증 참고 | 35일 ISS |
| CBC / CMP | 임상 검증 참고 | OSD-569 (인간) |

→ v2.0·v3.0에서 Category I (멀티모달 통합) task로 발전.

### 4.4 Pathway Enrichment Analysis — 로컬 실행 정책

OSDR은 fGSEA 결과를 제공하지 않음 (DESeq2 DGE까지만 처리).
GeneLab 파이프라인 v.D 이후 DGE에 Wald statistic (`Stat_` 컬럼) 포함 →
이를 fGSEA ranking metric으로 사용. 상세 설계: DESIGN_DECISIONS.md DD-15.

| Gene Set DB | MSigDB Category | Set 수 | 역할 | 적용 |
|---|---|---|---|---|
| **Hallmark** | H | 50 | Primary — Category C Method C, 논문 main figures | C, B해석 |
| **KEGG** | C2:CP:KEGG_MEDICUS | ~658 | Secondary — 대사 경로 중심 (간 분석) | C보조, J비교 |
| **Reactome** | C2:CP:REACTOME | ~1,787 | Secondary — 포괄적 경로 커버리지 | C보조, J비교 |

두 수준의 분석:
1. **fGSEA** (그룹): Flight vs Ground DEG의 pathway enrichment → pathway 선택, 해석
2. **GSVA** (샘플): 개별 샘플의 pathway activity score → ML feature (전 카테고리 적용)

```
scripts/run_fgsea.R              — 그룹 수준 fGSEA (미션×조직별)
scripts/compute_pathway_scores.R — 샘플 수준 GSVA pathway score
scripts/preprocess_pathways.py   — Python 통합 + Category C feature 생성
```

---

## 5. v1.0 Benchmark Task 카탈로그 (~25개)

### Category A: Spaceflight Detection (조직별 이진 분류)

**입력**: log2(DESeq2 size-factor normalized counts), per-sample
**출력**: Flight (1) vs Ground (0)
**Split**: LOMO (Leave-One-Mission-Out)
**평가**: AUROC (primary), AUPRC, macro-F1
**주의**: LFC 사용 금지 — LFC는 flight/ground 레이블 없이 계산 불가 (label leakage)

| Task | 조직 | 미션 수 | 샘플 수(예상) | 난이도 | 대조군 Track | 상태 |
|---|---|---|---|---|---|---|
| **A1** | Liver | **6** (RR-1,3,6,8,9 + MHU-2) | ~90 | Calibration | 1+2 | ✅ 확정 |
| **A2** | Gastrocnemius | 3 (RR-1, RR-5, RR-9) | ~60 | Standard | 2 | ✅ 확정 |
| **A3** | Kidney | 3 (RR-1, RR-3, RR-7) | ~50 | Standard | 1+2 | ✅ 최소 기준 충족 |
| **A4** | Thymus | 3 (RR-6, MHU-2, RR-9) | ~50 | Standard | 2 | ✅ 최소 기준 충족 |
| **A5** | Skin | **2** (MHU-2 dorsal/femoral + RR-6) | ~40 | Advanced | 2 | ⚠️ 미션 부족 → supplementary 강등 검토 |
| **A6** | Eye/Retina | **3** (RR-1, RR-3 + OSD-397 TBD) | ~45 | Advanced | 2 | ✅ OSD-397 bulk RNA-seq 확인. 미션 ID 검증 필요 |
| **A7** | Brain (bulk only) | TBD | TBD | TBD | 확인 후 | 보류 |

> 📊 **MHU-2 특이사항** (OSD-686): 3그룹 존재 (uG, 인공중력 1G, 지상 GC).
> A1에서는 uG vs GC 비교. 인공중력(AG) 그룹 → **D6 task** (인공중력 효과 예측)로 활용.
> 이는 유일한 인공중력 조건 — 독특한 가치.

**생물학적 Ground Truth 검증** (sanity check용, task 정의와 별개):
- A1 (간): ANGPTL4, PCK1 상향조절, 미토콘드리아 유전자 세트 (Cell 2020)
- A2 (근육): MuRF1/TRIM63, MAFbx/FBXO32 (근위축 E3 리가제) 탐지 여부
- A3 (신장): Nrf2 경로 활성화 (OSD-253 참조)

---

### Category B: Cross-Mission Generalization (크로스미션 일반화)

**핵심 질문**: 한 미션으로 학습한 모델이 독립적 미션에서 동작하는가?
**설계 원칙**: LFC feature 사용, 배치효과 제거 없이 평가 → 진정한 일반화 측정
**평가**: Transfer AUROC (학습·테스트 미션 쌍 전체 6×6 행렬)

| Task | 조직 | 설명 | 난이도 |
|---|---|---|---|
| **B1** | Liver | LOMO full 6×6 cross-mission matrix | Advanced |
| **B2** | Liver | ISS missions → STS-135 (shuttle) | Frontier |
| **B3** | Liver | ISS → BION-M1 (biosatellite) | Frontier |
| **B4** | Gastrocnemius | LOMO 3-mission | Advanced |
| **B5** | Kidney | LOMO 4-mission | Advanced |
| **B6** | Liver | ISS flight vs HU analog (지상 아날로그) | Advanced |

> B6는 HU 데이터 catalog 확인 후 포함 여부 결정.

**배치보정 비교 실험** (B1~B5 각각):
- 보정 없음 vs ComBat-seq vs limma::removeBatchEffect vs RUVseq 각각의 transfer AUROC 비교
- Harmony는 단일세포 임베딩용 → bulk RNA-seq에는 비표준. 참고용으로만 포함, 기본 비교 대상 아님
- → Category J (파이프라인) 태스크와 연결

---

### Category C: Cross-Tissue Transfer (크로스티슈 전이)

**핵심 질문**: 한 조직의 spaceflight LFC 서명이 다른 조직에서도 유효한가?
**방법**: Gene rank 또는 LFC 기반. 공통 유전자 세트로 모델 전이.
**설계 명확화**: 같은 미션 내 다른 조직 간 비교 (예: RR-1 liver feature → RR-1 kidney 예측)

| Task | 학습 조직 | 테스트 조직 | 공유 생물학 | 기대 Transfer AUROC | 난이도 |
|---|---|---|---|---|---|
| **C1** | Liver | Kidney | 대사, 산화스트레스 | 0.60~0.70 | Advanced |
| **C2** | Liver | Gastrocnemius | 에너지 대사 | 0.55~0.65 | Advanced |
| **C3** | Liver | Thymus | 면역 (Kupffer cell) | 0.55~0.65 | Advanced |
| **C4** | Thymus | Kidney | 면역-신장 | 0.55~0.60 | Frontier |
| **C5** | 모든 조직 풀 | holdout 조직 | Universal signature | ~random | Frontier |

**3가지 방법 병렬 비교** (각 C task에 대해 모두 실행, H3 검증):

| 방법 | 설명 | 예상 성능 | 역할 |
|---|---|---|---|
| 방법 A (Gene-level covariate shift) | 학습 조직 분류기 가중치 → 테스트 조직에 직접 적용 | 낮음 (조직 특이 유전자 문제) | 하한선 |
| 방법 B (DEG feature transfer) | 학습 조직 top-200 DEG → 테스트 조직에서 같은 유전자로 새 모델 학습 | 중간 | 중간 baseline |
| 방법 C (Pathway-level) [PRIMARY] | 학습 조직 fGSEA top-20 pathway (Hallmark, padj<0.05, \|NES\| 기준) → 테스트 조직 GSVA score → 분류 | 높음 (pathway 보존) | H3 핵심 |

*기대 발견: 방법 A < 방법 B << 방법 C → "spaceflight 반응은 유전자보다 경로 수준에서 더 잘 보존"*

**Gene-set 수준 병행 분석**:
- 조직 쌍별 DEG Jaccard index (top-200 DEG 기준)
- 경로(pathway) 수준 보존도: GSEA NES rank correlation (Spearman r)
- STRING 네트워크에서 공통 허브 유전자 식별

---

### Category D: Condition & Confounder Prediction (조건 예측)

**목적 이중성**: (1) 우주비행 조건을 전사체로 역추정, (2) 주요 confound 변수 정량화

| Task | 입력 | 예측 대상 | 유형 | 생물학적 의미 | 난이도 |
|---|---|---|---|---|---|
| **D1** | Gene expression (normalized counts) | 비행 기간 클래스 (3-class ordinal: ≤14d / 26~40d / ≥75d) | Ordinal classification | 기간-반응 구분 가능성 | Advanced |
| **D2** | Gene expression | 비행 기간 (이산: ≤20d / 21~40d / >40d) | 3-class | 임상 위험 단계 | Standard |
| **D3** | Gene expression | 미션 ID (다중 클래스, 간) | 6-class | 미션별 발현 차이 | Standard |
| **D4** | Gene expression | 마우스 계통 (C57BL/6 vs BALB/c) | Binary | **confounder 정량** | Calibration |
| **D5** | Gene expression | 하드웨어 유형 (AEM vs Shuttle) | Binary | **confounder 정량** | Calibration |
| **D6** | Gene expression | 인공중력 여부 (uG vs 1G-centrifuge vs GC) | 3-class | **MHU-2 인공중력 효과 분리** | Advanced |

> **D6 데이터 출처**: OSD-686 (MHU-2) 유일하게 3그룹 보유 (uG n=3, AG-centrifuge n=3, GC n=3).
> "우주비행 반응 = 미세중력 + 방사선 + 스트레스" 중 미세중력 성분만 분리 가능한 독특한 task.
>
> **D1 설계 근거**: 실제 비행 기간이 3~4 클러스터에 집중 (≤14d/26~40d/≥75d).
> 연속 회귀(Pearson r) 사용 시 클러스터 구조를 무시하고 성능을 과대 평가. Ordinal 3-class로 재설계.
> Primary metric: macro-F1 (ordinal 특성 반영: Spearman r between predicted/true class rank)
>
> **D4, D5 해석 주의**: 계통(D4), 하드웨어(D5)는 미션/기간과 공선성.
> 높은 예측 성능 = confounder 독립 효과가 아닌 "해당 변수의 미션별 배치 차이"의 상한선으로 해석.
> 논문 Methods에 명시 필수: "D4 AUROC ≤ D3 AUROC (BALB/c는 RR-3에서만 사용)"

#### ✅ Category D 구현 결과 (2026-03-01)

**실현 가능성 평가**: D1/D2 (비행 기간) 비실현 — 모든 미션 30-40일 범위. D4 (계통) — thymus만 2계통 보유. D5 (하드웨어) — 큐레이션 필요. → **D3 + D6만 구현.**

| Task | Metric | Gene-level | Pathway-level | p-value | 비고 |
|---|---|---|---|---|---|
| **D3** (Liver 6-class) | macro-F1 | **1.000** [1.00, 1.00] | 0.056 [0.04, 0.08] | <0.001 (gene) / NS (pathway) | 유전자=완벽 미션 분리, 경로=배치 불변 |
| **D6 Liver** (MHU-2) | macro-F1 | **0.886** [0.56, 1.00] | 0.413 NS | 0.002 (gene) | uG/AG/GC 3-class LOO |
| **D6 Thymus** (MHU-2) | macro-F1 | **0.657** [0.33, 0.90] | 0.641 [0.28, 1.00] | 0.037 (gene) / 0.052 (pathway) | 유전자≈경로 |

> **핵심 발견**: D3 gene F1=1.000은 강력한 미션 배치 효과 증거. Pathway F1=0.056은 GSVA Hallmark이
> 배치에 불변(batch-invariant)함을 확인 → 이는 Category B/C에서 pathway 사용의 과학적 근거.
>
> **방법**: D3 = RepeatedStratifiedKFold (5×10), D6 = LOO-CV (n=9, n_per_class=3).
> 스크립트: `scripts/condition_prediction.py`. 결과: `evaluation/D_condition_summary.json`.

---

### Category J: Pipeline Standardization Benchmark (v1.0 포함 subset)

v1.0에서는 핵심 파이프라인 비교 task만 포함:

| Task | 비교 대상 | 평가 지표 |
|---|---|---|
| **J1** | GeneLab pipeline v1 vs v2 vs v3 | DEG Jaccard (top-100/500), AUROC 변화 |
| **J2** | DESeq2 vs edgeR vs limma-voom | Top DEG 일치율, Category A AUROC 영향 |
| **J3** | 배치보정 없음 vs ComBat-seq vs limma::removeBatchEffect vs RUVseq | Category B transfer AUROC 변화 (Harmony는 참고용만: bulk RNA-seq 비표준) |
| **J4** | 전체 유전자(22k) vs 필터링(5~8k) | 연산비용 대비 성능 trade-off |
| **J5** | Gene-level vs Pathway-level (GSVA) feature representation | 전 카테고리 (A,B,C,D) 성능 비교 |

> J5: 전 카테고리에서 gene-level vs pathway-level feature의 AUROC/F1 체계적 비교.
> Hallmark(50 pathways) primary, KEGG/Reactome secondary. DD-15 참조.
>
> J6 (STAR vs HISAT2), J7 (mm10 vs mm39): 원본 FASTQ 재처리 필요 → v3.0으로 이동.

#### ✅ J5 구현 결과 (2026-03-01): Gene vs Pathway 체계적 비교

12개 비교 (A×5, C×4, D×3). 전체: Gene wins 8/12, Pathway wins 4/12.

**Category A (Spaceflight Detection, LOMO AUROC)**:

| 조직 | Gene | Pathway | Diff (P-G) | Winner |
|---|---|---|---|---|
| Liver | 0.670 | 0.574 | -0.10 | Gene |
| Gastrocnemius | 0.824 | 0.688 | -0.14 | Gene |
| **Kidney** | 0.432 | **0.743** | **+0.31** | **Pathway** ★ |
| Thymus | 0.923 | 0.879 | -0.04 | Gene |
| **Eye** | 0.789 | **0.915** | **+0.13** | **Pathway** ★ |

> **핵심 발견 — Kidney Rescue**: Gene-level AUROC 0.43 (chance 이하) → Pathway-level 0.74로 "구조".
> Eye도 0.79→0.92 향상. **Pathway가 gene-level 실패를 보완하는 조직 존재** → H3 부분 지지.

**Category C/D 요약**: C에서 pathway 2/4 승리 (liver→gastro에서 C=0.867>>A=0.563), D에서 gene 3/3 승리 (배치/조건 구분에 유전자 우월).

> **Mean diff (Category A)**: +0.03 (본질적 동등). **전체**: -0.11 (D 카테고리의 배치 효과로 gene 쏠림).
> 스크립트: `scripts/gene_vs_pathway_comparison.py`. 결과: `evaluation/J5_gene_vs_pathway.json`.

---

## 6. 평가 체계

### 6.1 Split 전략

| Split 유형 | 독립성 단위 | 적용 Task | 설명 |
|---|---|---|---|
| **LOMO** (Leave-One-Mission-Out) | Mission (cage-level leakage 구조적으로 없음) | A, B | 미션 단위 holdout. N_mission ≥ 3 |
| **LOTO** (Leave-One-Tissue-Out) | Tissue | C | 조직 단위 holdout |
| **Mission-stratified 80/20** | Mission (같은 미션은 같은 fold) | D | ⚠️ 랜덤 샘플 분리 금지 — mission-level 독립성 보장 |
| **GLDS-stratified 80/20** | GLDS study | J | 동일 GLDS 샘플은 train 또는 test 중 하나만 |
| **Cross-species** | Species | E (v2.0) | 종 단위 완전 분리 |

> ⚠️ **D, J 카테고리 split 주의**: 기존 "80/20 stratified" 설계에서 같은 미션/GLDS의 샘플이
> train/test에 분리될 경우 batch leakage 발생. 반드시 **mission/GLDS 단위 그룹 분리** 적용.

### 6.2 평가 지표

| Task 유형 | Primary | CI 방법 | 모델 간 비교 | Secondary |
|---|---|---|---|---|
| 이진 분류 (A, B, C) | AUROC | Bootstrap 95% CI (n=1000) | DeLong's test | AUPRC, macro-F1, Precision@50 |
| 다중 클래스 (D3~D5) | macro-F1 | Bootstrap 95% CI (n=1000) | Wilcoxon signed-rank | per-class F1, accuracy |
| 회귀 (D1) | Spearman r | Bootstrap 95% CI (n=1000) | — | RMSE, R² |
| Transfer (B) | Transfer AUROC | Bootstrap 95% CI | DeLong's test | Source→Target delta |
| 파이프라인 (J) | DEG Jaccard | Bootstrap 95% CI | Wilcoxon | Rank correlation, AUROC delta |

**보고 형식 표준** (모든 task):
```
AUROC = 0.81 [95% CI: 0.74–0.87], permutation p < 0.001
LOMO fold 평균 ± SD: 0.81 ± 0.06
```

**다중 검정 보정 정책**:
- task 내 모델 비교: DeLong's test + Bonferroni correction
- task 간 비교: Benjamini-Hochberg FDR (q < 0.05)
- Composite score: CI 비중첩 = 유의미한 차이로 판단

### 6.3 점수 정규화 및 Composite Score

```python
# Task별 정규화 (SpaceOmicsBench 호환)
normalized_score = (model_score - random_baseline) / (1.0 - random_baseline)

# 카테고리별 평균
category_score = mean(normalized_scores within category)

# Composite score (카테고리 동일 가중치)
composite = mean([score_A, score_B, score_C, score_D, score_J])
```

**Random baseline — task 유형별 올바른 계산** (P0 수정: metric과 동일 스케일 유지):

| Task 유형 | Primary Metric | Random Baseline | 수식 근거 |
|---|---|---|---|
| 이진 분류 (A, B, C) | AUROC | **0.5** | 랜덤 분류기의 기대 AUROC |
| 다중 클래스 (D3~D5) | macro-F1 | **1/n_classes** (macro-F1 기준) | ≠ 1/K accuracy. macro-F1은 클래스별 F1 평균 |
| 회귀 (D1) | Spearman r | **0.0** | 랜덤 예측의 기대 상관 |
| Transfer (B) | Transfer AUROC | **0.5** | 랜덤 전이 |

> ⚠️ 구버전 오류: "다중 클래스 random baseline = 1/K accuracy"는 primary metric인 macro-F1과
> 스케일이 달라 정규화 수식이 수학적으로 잘못됨. 위 표 기준으로 통일.

### 6.4 난이도별 집계

| 난이도 | Task 수 (예상) | 해석 |
|---|---|---|
| Calibration | 3~4개 | 기준선 확인. 여기서 실패 시 데이터 문제 |
| Standard | 6~8개 | 합리적 접근으로 개선 가능 |
| Advanced | 8~10개 | 의미 있는 개선 가능하지만 어려움 |
| Frontier | 3~5개 | 현재 방법론의 한계 탐색 |

---

## 7. 데이터 품질 필터링

### 7.1 RNA-seq 샘플 수준 QC

| 지표 | 제외 기준 | 근거 |
|---|---|---|
| STAR mapping rate | < 70% | 오염 또는 저품질 |
| rRNA 비율 | > 20% | 라이브러리 실패 |
| Library size (total reads) | < 1,000,000 | 통계적 검정력 부족 |
| 검출 유전자 수 | < 5,000 | 표현 다양성 부족 |
| PCA outlier | > 4 SD from centroid | 기술적 이상값 |
| 샘플 간 상관 | Pearson r > 0.99 (같은 그룹 내 아닌 경우) | 중복 샘플 |

### 7.2 연구(GLDS) 수준 QC

| 조건 | 처리 |
|---|---|
| 생물학적 반복 수 < 3/그룹 | Category A 제외, J만 허용 |
| Flight 또는 Ground 그룹 누락 | Task 제외 |
| 메타데이터 누락 (성별, 계통, 기간) | 수동 큐레이션 또는 제외 |
| Bulk RNA-seq이 없는 조직 | Category A 제외 |
| GeneLab 처리 파이프라인 버전 | 기록 필수 (J 카테고리 활용) |

### 7.3 미션 수준 QC

```python
# Task 포함 기준 (하드 컷오프)
def is_eligible_for_task(missions_list):
    return (
        len(missions_list) >= 3 and          # LOMO 의미 있으려면
        sum(n_samples) >= 50 and             # 학습 세트 30+ 확보
        all(n >= 3 for n in n_per_group)     # 통계적 최소 요건
    )
```

---

## 8. 방사선·환경 메타데이터 활용 (v1.0에서의 역할)

v1.0에서 RadLab / EDA 데이터는 **task 입력이 아닌 메타데이터 covariate**로만 활용.
H 카테고리(방사선-오믹스 통합)는 설계 문제로 **v3.0으로 이동** (아래 참조).

### v1.0 활용 방식

```
1. missions.json에 각 미션의 총 방사선량 (Gy) 기록
2. Category B 결과 분석 시 secondary 변수로:
   - 미션 간 transfer AUROC와 방사선량 차이의 상관 분석
   - "방사선량이 비슷한 미션끼리는 더 잘 전이되는가?"
3. Category D 결과 해석 시:
   - D3 (미션 ID 예측) 성능과 방사선량 차이 매트릭스 비교
```

### v3.0에서 H 카테고리 재설계 방향

현재 H 카테고리의 문제: ISS 내 모든 마우스는 같은 방사선량 → sample-level regression 불가.

v3.0 해결책:
- OSDR에서 **지상 방사선 단독 실험** 데이터 수집 (Fe-56, proton, X-ray 선량별)
- 방사선 선종 × 선량 × 조직의 DEG 패턴 → 진정한 dose-response 모델
- ISS 데이터와 지상 방사선 데이터 비교 → "ISS 특이 vs 방사선 단독" 분리

---

## 9. 프로젝트 디렉토리 구조

```
GeneLab_benchmark/
├── PLAN.md                          # 이 파일
├── README.md                        # 공개용 소개 (v1.0 기준)
├── DATA_CATALOG.md                  # 전체 OSD/GLDS 목록 (자동 생성)
│
├── scripts/
│   ├── catalog_datasets.py          # OSDR API → DATA_CATALOG.md 자동 생성
│   │                                # + GLDS ID 검증 + control_type 확인
│   ├── fetch_osdr.py                # OSDR 데이터 다운로드 (fetch_genelab.py 기반)
│   ├── fetch_radlab.py              # RadLab API → 미션별 방사선량 수집
│   ├── fetch_eda.py                 # EDA API → 환경 텔레메트리 수집
│   ├── build_missions_json.py       # missions.json 통합 생성
│   ├── quality_filter.py            # QC 필터링 (GLDS 수준 + 샘플 수준)
│   ├── preprocess/
│   │   ├── preprocess_mouse_rnaseq.py   # LFC + rank feature 생성
│   │   └── batch_correction.py          # ComBat-seq, Harmony
│   ├── run_fgsea.R                  # 그룹 수준 fGSEA (DD-15)
│   ├── compute_pathway_scores.R     # 샘플 수준 GSVA pathway score (DD-15)
│   ├── preprocess_pathways.py       # Python 통합 + Category C feature 생성
│   ├── generate_tasks.py            # Task JSON + split 파일 생성
│   └── validate_tasks.py            # 스키마 검증
│
├── data/
│   ├── mouse/
│   │   ├── liver/                   # GLDS-48, 137, 168, 242, 245, 379, 617
│   │   ├── gastrocnemius/           # GLDS-101, 326, 401
│   │   ├── kidney/                  # GLDS-102, 163, 253, 674
│   │   ├── thymus/                  # GLDS-4, 244, 289, 421
│   │   ├── skin/                    # GLDS-238, 239, 243, 689
│   │   ├── eye/                     # GLDS-100, 194, 664
│   │   └── [future_v2]/             # brain, heart, bone_marrow (v2.0)
│   └── [future_v2]/                 # human/, celegans/, arabidopsis/, microbiome/
│
├── mission_metadata/
│   ├── missions.json                # 전체 미션 메타데이터
│   ├── radlab_summary.csv           # 미션별 방사선량
│   └── eda_summary.csv              # 미션별 환경 데이터
│
├── processed/                       # 벤치마크용 정제 CSV
│   ├── A_detection/
│   ├── B_cross_mission/
│   ├── C_cross_tissue/
│   ├── D_condition/
│   ├── J_pipeline/
│   ├── fgsea/                       # 그룹 수준 fGSEA 결과 (DD-15)
│   │   ├── {tissue}/                # {mission}_fgsea_{db}.csv
│   │   └── summary/                 # all_fgsea_{db}.csv (전 미션 통합)
│   └── pathway_scores/              # 샘플 수준 GSVA scores (DD-15)
│       └── {tissue}/                # {mission}_gsva_{db}.csv
│
├── tasks/                           # Task 정의 JSON (SpaceOmicsBench 스키마 호환)
├── splits/                          # Train/test split JSON
│
├── evaluation/
│   ├── eval_harness.py
│   ├── metrics.py
│   ├── pipeline_compare.py
│   ├── B_cross_mission_summary.json    # Category B 결과 (5 tissues)
│   ├── C_cross_tissue_summary.json     # Category C 결과 (4 pairs × 3 methods)
│   ├── D_condition_summary.json        # Category D 결과 (D3, D6×2)
│   └── J5_gene_vs_pathway.json         # J5 비교 결과 (12 comparisons)
│
├── baselines/
│   ├── run_baselines.py             # LogReg, RF, XGBoost, LightGBM, Random
│   └── baseline_results.json
│
└── docs/
    ├── DOWNLOAD_INSTRUCTIONS.md
    ├── DATA_INVENTORY.md            # GLDS ID 전수 검증 결과 포함
    ├── PROVENANCE.md
    ├── BIOLOGICAL_GROUND_TRUTH.md  # 검증 기준 유전자 세트 (Cell 2020, SOMA 2024)
    └── DESIGN_DECISIONS.md         # control_type, feature representation 등 결정 사항
```

---

## 10. 개발 로드맵 (v1.0)

### Phase 1: 데이터 인프라 + 간(Liver) 검증

**목표**: GLDS 전수 검증 + 간 데이터 완성 + Task A1, B1 동작 확인

**체크포인트** (3조건 AND 충족 시 Phase 2 진행):
1. A1 LOMO AUROC > 0.70, 95% CI 하한 > 0.60
2. Permutation p < 0.05 (랜덤 레이블 대비 유의미한 신호)
3. Known spaceflight genes (ANGPTL4, PCK1) SHAP top-50 이내

*3조건 중 하나라도 불충족 시: QC 재검토 → 미션별 mapping rate·샘플 수 확인*

```
[ ] catalog_datasets.py
    - OSDR API로 대상 GLDS ID 실제 존재 확인
    - 각 GLDS의 bulk RNA-seq 파일명·형식 확인
    - control_type (GC/VC/BC) 기록
    - → DATA_CATALOG.md + GLDS_verified.json 자동 생성

[ ] fetch_osdr.py
    - 간 GLDS 7개 다운로드 (GLDS-48, 137, 168, 242, 245, 379, 617)
    - DESeq2 normalized counts + sample metadata 우선

[ ] fetch_radlab.py + fetch_eda.py
    - 각 미션 방사선량, 환경 데이터 수집
    - build_missions_json.py → missions.json 생성

[ ] quality_filter.py
    - 샘플 수준 QC (mapping rate, library size, PCA outlier)
    - 미션 수준 QC (n_per_group ≥ 3 확인)
    - QC 리포트 자동 생성

[ ] preprocess_mouse_rnaseq.py
    - Raw counts → log2(DESeq2 normalized counts) per sample → Category A용
    - Raw counts → LFC (apeglm shrinkage, vs mission GC/VC, 두 Track) → Category B/C용
    - Gene rank 계산
    - Low-expression filter 적용 (전체 샘플 기준. variance filter는 LOMO loop 내부에서)
    - 표준화 CSV 저장

[ ] Task A1, B1 구현 + Baseline 실행
    - LOMO split 생성
    - LogReg, RF, LightGBM baseline
    - A1 sanity check: ANGPTL4, PCK1 feature importance 상위?

✅ fGSEA + GSVA 파이프라인 구축 (DD-15)
    - ✅ run_fgsea.R, compute_pathway_scores.R, preprocess_pathways.py 작성
    - ✅ Liver 6개 미션 fGSEA 실행 (Hallmark+KEGG+Reactome)
    - ✅ Liver GSVA score 계산 + gene vs pathway 비교
    - ✅ 다조직 fGSEA + GSVA 확장 (5 tissues × 17 missions × 3 DBs = 51 fGSEA + 54 GSVA)
```

### Phase 2: 다조직 확장 + C, D, J 구현

**목표**: 6개 조직 × 10+ 미션 + 전체 v1.0 task 완성

```
✅ Gastrocnemius, Kidney, Thymus, Eye 데이터 추가 (Skin 제외: 미션 수 부족)
✅ Category B 구현 (cross-mission transfer, 5 tissues, bootstrap CI + permutation)
✅ Category C 구현 (cross-tissue transfer)
    - ✅ 방법 A(gene), B(DEG), C(pathway) 3방법 × 4 tissue pairs
    - ✅ J5: gene-level vs pathway-level 전 카테고리 비교
    - H3 결과: 조건부 지지 (tissue-pair dependent)
✅ Category D 구현 (condition prediction)
    - ✅ D3 (mission ID 6-class, liver): gene F1=1.0 / pathway F1=0.06
    - ✅ D6 (MHU-2 gravity 3-class): liver gene F1=0.886, thymus gene F1=0.657
✅ Category J5 구현 (gene vs pathway 12 comparisons)
[ ] Category J1-J4 구현 (pipeline comparison: DESeq2 vs edgeR vs limma)
[ ] batch_correction.py (ComBat-seq, Harmony) → B 카테고리 재평가
[ ] 전체 baseline 실행 + baseline_results.json 생성
[ ] 난이도 계층화 검증 (예측 AUROC vs 실제 AUROC 비교)
```

### Phase 3: 문서화 + 논문 초안

```
[ ] BIOLOGICAL_GROUND_TRUTH.md 완성
[ ] DESIGN_DECISIONS.md 완성
[ ] README.md 공개용 작성
[ ] 논문 초안 (Genome Biology)
    - Methods: 데이터, task 설계, split, metric
    - Results: baseline 성능, cross-mission generalization 발견
    - Discussion: 어떤 조직/미션 조합이 일반화 가능한가
```

---

## 11. Future Extensions (v2.0 · v3.0 개요)

### v2.0: 다종·단일세포·마이크로바이옴 (별도 논문)

**전제 조건**:
- v1.0 출판 완료
- C. elegans OSDR 데이터 존재·품질 직접 확인
- 인간 데이터 SpaceOmicsBench 연동 방식 확정

**Category E (Cross-Species Conservation)**:
- E1: 마우스 간 LFC → 인간 혈액 LFC 경로 수준 보존도 (Spearman r on GSEA NES)
- E2: C. elegans 유전체 → 마우스 ortholog 공통 DEG (검증 후 설계)
- E3: Arabidopsis ROS/열충격 유전자 → 마우스 상동 유전자 발현
- ※ E1: DEG 직접 비교 대신 **경로(pathway) 수준 rank correlation** 사용 (Direct DEG list comparison은 조직 유형이 달라 부적합)

**Category F (Single-Cell & Spatial)**:
- F1: RRRM-2 뇌 snRNA-seq (GLDS-589) — 세포유형별 spaceflight 반응
- F2: RRRM-1 다조직 scRNA-seq (GLDS-746~762) — 세포 조성 변화 (검증 후)
- F3: RR-3 심장 Visium (GLDS-270) — 공간 영역별 분류
- F4: RR-18 뇌 GeoMx (GLDS-613~628) — 뇌 영역별 분류

**Category G (Microbiome)**:
- G1: ISS 환경 마이크로바이옴 출처 분류 (OSD-572,573)
- G2: 마우스 장 마이크로바이옴 flight vs ground
- G3: 인간 마이크로바이옴 시계열 (OSD-630)
- G6: 마이크로바이옴 × 호스트 전사체 공통 변화 (paired 샘플 확인 필요)

### v3.0: 방사선·멀티모달·파이프라인 완전판 (별도 논문)

**Category H (Radiation-Omics) — 재설계 방향**:
- ISS 데이터만으로는 샘플 내 방사선량이 상수 → H 카테고리 불가
- 지상 방사선 실험 OSDR 데이터 추가 (HZE, proton, X-ray 선종별)
- Task H1: 지상 방사선 선량별 전사체 dose-response (진짜 regression task)
- Task H2: ISS flight vs 선종별 지상 방사선 비교 (ISS 특이 서명 식별)

**Category I (Multi-Modal)**:
- I1: RNA-seq + 체중/장기무게 통합 (AUROC 향상도 측정)
- I2: RNA-seq + MicroCT (RR-1, OSD-804)
- I3: RNA-seq + RadLab + EDA 완전 멀티모달

**Category J 완전판**:
- J5: STAR vs HISAT2 (원본 FASTQ 재처리)
- J6: 참조 게놈 버전 (mm10 vs mm39)

---

## 12. 벤치마크 타겟 모델 & Dataset Freeze 정책

### 12.1 타겟 모델 계층 (3-tier)

**Tier 1 — Classical ML Baselines (Phase 1, 현재 구현)**
- Logistic Regression (ElasticNet), Random Forest, XGBoost, LightGBM
- PCA-LR (low-dimensional baseline)
- 평가 방식: LOMO AUROC, Bootstrap 95% CI, permutation p-value

**Tier 2 — Gene Expression Foundation Models (Phase 2, DD-13)**
- Geneformer (Theodoris et al. 2023, Nature): gene rank token → fine-tune classification
- scGPT (Cell Systems 2024): gene expression → pretraining → downstream task
- 입력 형식: OSDR normalized counts → gene rank order (token sequence)
- 태스크: A (spaceflight detection), B (cross-mission generalization)

**Tier 3 — Text-based LLMs (Phase 2, 신규 추가)**
- GPT-4o, Claude Opus, Llama 3 등 텍스트 LLM
- 입력 형식 (zero-shot/few-shot):
  ```
  System: "You are a bioinformatics expert analyzing spaceflight transcriptomics."
  User: "The following 50 genes are most variable in a mouse {tissue} sample:
         Gene1 (ENSMUSGXXXXXXX), Gene2, ..., Gene50.
         Is this sample from spaceflight (Flight) or ground control (Ground)?
         Answer: Flight or Ground, with confidence 0-1."
  ```
- 평가: zero-shot accuracy vs ML baseline; few-shot (n=3, 5, 10 examples) scaling
- 비교 포인트: LLM의 생물학 도메인 지식이 flight 예측에 도움이 되는가?

**타겟 모델 요약표**:

| Tier | 모델 유형 | 입력 | 구현 단계 |
|------|----------|------|---------|
| 1 | Classical ML | Gene expression matrix | Phase 1 (완료) |
| 2 | Gene FM (Geneformer) | Gene rank tokens | Phase 2 |
| 3 | Text LLM (GPT-4o) | Gene list as text | Phase 2 |

### 12.2 v1.0 Dataset Freeze 정책

**원칙**: v1.0 출판 시 포함되는 미션/데이터는 완전히 고정됨. 이후 새 미션은 v2.0으로.

**v1.0 포함 기준**:
- OSDR 공개 기준일: **2026-03-01** 이전 공개 미션
- 처리 완료 기준: normalized counts + metadata QC 통과
- 최소 요건: n_flight ≥ 3, n_ground ≥ 3 (DD-11 준수)

**현재 확정 v1.0 미션 (Track 2a, C57BL/6J)**:

| 조직 | 미션 | OSD | n | Track | 상태 |
|------|------|-----|---|-------|------|
| Thymus | MHU-2 | OSD-289 | 9 | 2a | ✓ |
| Thymus | RR-6 | OSD-244 | 54 | 2a | ✓ |
| Thymus | RR-9 | OSD-421 | 20 | 2a | ✓ |
| Thymus | RR-23 (OSD-515) | OSD-515 | 27 | 2a | 미다운로드 |
| Gastrocnemius | RR-1, RR-5, RR-9 | 복수 | ~50 | 2a | ✓ |

**v1.0 Held-Out 결정 (2026-03-01 확정)**:
- OSD-515 (RR-23): **v1.0 held-out test set** (labels 비공개) — leaderboard 외부 평가용
- A4 LOMO folds (MHU-1/2, RR-6, RR-9) + OSD-515 held-out = v1.0 thymus 태스크 최종 구성

---

## 13. v1.0 MVP 범위 재정의 (2026-03-01 결정)

**배경**: 원래 PLAN.md (v0.5) Category C/D/J 포함 25+ tasks는 단기 출판에 비현실적. 공개 벤치마크 선례 (MAQC, CAGI)에서 v1.0은 명확한 범위로 시작 후 버전 확장이 표준.

### v1.0 확정 범위 (2026-03-01 업데이트)

| 항목 | v1.0 포함 | v2.0+ 이동 | 비고 |
|------|----------|-----------|------|
| Category A (탐지) | 5 tissues LOMO (A1-A6) | Soleus/EDL 분리 | 완료 |
| Category B (전이 행렬) | 5 tissues 전체 | — | 완료, 5-tissue summary |
| Category C (조직 간 전이) | 4 pairs × 3 methods | — | 완료 |
| Category D (조건 예측) | D3+D4+D5+D6 | — | 완료 (mission/strain/hardware/gravity) |
| Category J | J3 (batch) + J5 (gene vs pathway) | J1/J2/J4 | J3+J5 완료 |
| Negative controls | NC1 (permutation) + NC2 (housekeeping) | — | 완료 |
| External validation | Cell 2020 concordance | SOMA 2024 full | 완료 |
| Held-out test set | OSD-515 thymus | GAS held-out 미결 | |
| Tier 1 (Classical ML) | 완료 | — | |
| Tier 2 (Geneformer) | HPC 대기 | scGPT 등 v2.0 | |
| Tier 3 (Text LLM) | GPT-4o zero-shot | Claude, Llama v2.0 | |

### v2.0 범위 (v1.0 출판 후)
- 추가 조직 (Soleus/EDL, Skin, Brain)
- Category E: 다종 (C. elegans, 인간)
- Category F: 단세포 RNA-seq
- scGPT, BioGPT 등 추가 Foundation Model

### Geneformer 개발 전략 (MacBook Air M-series + Cornell Cayuga HPC)

**환경 제약**: Apple Silicon MPS, 24GB unified memory, 4-bit quantization 불가 (float32만).
**HPC**: Cornell Cayuga — Slurm v25.05.0, `scu-gpu` partition, A40 (48GB) / A100 (80GB).

| 단계 | 환경 | 내용 | 상태 |
|------|------|------|------|
| 1. 토크나이제이션 | MacBook (CPU) | ENSMUSG → ENSG ortholog 매핑, gene rank 생성 | ✅ 완료 |
| 2. 추론 테스트 | MacBook (MPS, float32) | 10M 파라미터 모델 로딩 + dry-run | ✅ 완료 |
| 3. Dev fine-tuning | MacBook (MPS, float32) | RR-9 fold 5-epoch 테스트 (진행 중) | 🔄 진행 중 |
| 4. HPC fine-tuning | Cayuga A40 GPU | LOMO 5-fold, 10 epoch, batch=16 | ⏳ 대기 |
| 5. 평가 | MacBook | 저장된 checkpoint → AUROC 계산 | ⏳ 대기 |

**구현 결과** (2026-03-01):
- 모델: Geneformer-V1-10M (10,263,298 파라미터, BERT 6-layer 256-dim)
- 토크나이제이션: `scripts/geneformer_tokenize.py` — 5 folds × 2 splits 완료
- Ortholog coverage: 57-64% (mouse ENSMUSG → human ENSG via Ensembl BioMart)
- 평균 token 길이: 618-2048 (fold/mission별 상이 — RR-9 test가 가장 짧음)
- HPC 스크립트: `scripts/hpc_submit_geneformer.sh` (Cayuga scu-gpu, sbatch array)
- 필요 패키지: `requirements_geneformer.txt`

**Cayuga HPC 사용법**:
```bash
# Login node — tokenize (CPU, fast)
python scripts/geneformer_tokenize.py --task A4 --model-version v1

# Submit array job (5 folds × A40 GPU)
sbatch scripts/hpc_submit_geneformer.sh
```

**Geneformer 선택 이유**: 마우스 단세포 데이터 사전학습 (30M 세포), 오픈소스, HuggingFace 제공, 가장 작은 Foundation Model (~10M 파라미터) → MPS 개발 가능, HPC fine-tuning 효율적.

---

## 13. 저널 타겟

| 버전 | 1순위 | 2순위 | 포지셔닝 |
|---|---|---|---|
| **v1.0** | **Genome Biology** | Nature Communications | 마우스 다조직 spaceflight 전사체 일반화 benchmark |
| **v2.0** | **Cell Systems** | npj Systems Biology | 다종·다해상도 우주생물학 벤치마크 |
| **v3.0** | **Nature Communications** | Bioinformatics | 방사선-오믹스 통합 + 파이프라인 표준화 |

> Nature Methods는 v2.0·v3.0 통합 후 big paper로 도전 검토.

---

## 13. "Unassailable" 전략

### 다른 그룹의 태클 예상 포인트 및 대응

| 예상 비판 | 대응 전략 |
|---|---|
| "샘플 수가 너무 적다" | Power analysis 명시. N 한계를 design으로 투명하게 공개. HU 아날로그 추가로 샘플 확대 |
| "배치효과와 생물학 분리 불가" | LFC feature 사용으로 mission-내 정규화. control_type 두 트랙. J 카테고리에서 배치보정 효과 직접 측정 |
| "SOMA와 중복" | SOMA는 인간 atlas, 이 벤치마크는 마우스 다조직 ML task. 인간 데이터는 v2.0 보조로만 |
| "파이프라인 선택에 결과가 종속" | J 카테고리에서 직접 측정·보고 — 이것이 오히려 강점 |
| "gold standard 없음" | 알려진 spaceflight signature (Cell 2020, SOMA 2024)를 외부 검증 세트로 명시 |

### 핵심 강점

1. **LOMO full cross-mission matrix (B1)**: 6×6 미션 전이 성능표 → 어디서 일반화되고 어디서 안 되는지 한눈에
2. **Control type 두 트랙**: 기존 연구들이 무시한 GC vs VC 차이를 명시적으로 정량
3. **Confounder 정량화 (D4, D5)**: 계통·하드웨어 효과를 숨기지 않고 측정
4. **AWG 내부 협력**: NASA 내부 전문가 검토 → 외부 비판 방어막

---

## 14. 참고 자료

### 주요 논문
- **Cell 2020**: Comprehensive Multi-omics / Mitochondrial stress hub (13조직, 59 astronauts)
- **SOMA (Nature 2024)**: 인간 멀티미션 멀티오믹스 atlas
- **npj Microgravity 2023**: Arabidopsis 15연구 메타분석
- **ISSOP (Patterns 2020)**: 국제 space omics 처리 표준
- **NASA OSDR NAR 2025**: OSDR 전체 데이터 유형 표 (Table 1)
- **SpaceOmicsBench**: 인간 멀티오믹스 benchmark 설계 참조

### 데이터 접근
- OSDR Search: https://osdr.nasa.gov/bio/repo/search
- GeneLab API: https://visualization.osdr.nasa.gov/biodata/api/v2/
- RadLab: https://visualization.osdr.nasa.gov/radlab/
- EDA: https://visualization.osdr.nasa.gov/eda/
- OSDR Data Curation Templates: https://github.com/nasa/OSDR_Data_Curation
- GeneLab Processing Pipeline: https://github.com/nasa/GeneLab_Data_Processing

### 관련 도구
- SpaceOmicsBench: /Users/jak4013/Dropbox/Bioinformatics/Claude/SpaceOmicsBench/v2_public/scripts/

---

## 15. Changelog

| 버전 | 날짜 | 주요 변경 |
|---|---|---|
| **v1.1** | 2026-03-01 | 벤치마크 강화: NC1 permutation summary (28 entries), NC2 housekeeping control (5 tissues), Cell 2020 external validation (71.7% pathway concordance), J3 batch correction comparison (H2 STRONGLY_SUPPORTED, mean delta 0.01), D4 strain + D5 hardware tasks, B 5-tissue summary 통합, utils.py 코드 중복 제거, §13 v1.0 scope 업데이트 (C/D/J5 모두 v1.0 포함). 새 스크립트: `aggregate_negative_controls.py`, `housekeeping_control.py`, `cell2020_validation.py`, `batch_correction_eval.py`, `utils.py`. 새 문서: `docs/BIOLOGICAL_GROUND_TRUTH.md`. |
| **v1.0** | 2026-03-01 | Category D 구현 (D3 liver 6-class, D6 MHU-2 3-class) + J5 gene vs pathway 체계적 비교 (12 comparisons). 핵심 발견: D3 gene F1=1.0 (배치 효과 증거), pathway F1=0.06 (배치 불변 확인); Kidney Rescue (gene 0.43→pathway 0.74). Phase 2 체크리스트 B/C/D/J5 완료 표시. 스크립트: `condition_prediction.py`, `gene_vs_pathway_comparison.py`. |
| **v0.9** | 2026-03-01 | Geneformer 파이프라인 구현 완료: `geneformer_tokenize.py` (Ensembl ortholog 매핑, 5 folds 토크나이제이션), `geneformer_finetune.py` (MPS/CUDA, BertForSequenceClassification), `hpc_submit_geneformer.sh` (Cornell Cayuga scu-gpu A40). §13 Geneformer 전략 업데이트 (Cayuga HPC 세부사항 포함). |
| **v0.8** | 2026-03-01 | §13 추가: v1.0 MVP 범위 재정의 (Category C/D/J → v2.0). OSD-515 held-out 결정 확정. Geneformer 개발 전략 (MacBook MPS + HPC). A4 within-LOO 4-fold 완료 기록. |
| **v0.7** | 2026-03-01 | §12 추가: 타겟 모델 3-tier (Classical ML / Gene FM / Text LLM) + v1.0 Dataset Freeze 정책. MHU-1 Track 2b 재분류 반영. 버그 수정 기록 (run_baselines.py B1+B2). |
| **v0.6** | 2026-02-28 | Pathway analysis 통합: §4.4 추가 (fGSEA+GSVA 로컬 실행 정책), Category C 방법 C 구체화 (fGSEA top-20 + GSVA score), J5 task 추가 (gene vs pathway feature 비교), 디렉토리 구조에 fgsea/pathway_scores 추가, DD-15 참조 연결 |
| **v0.5** | 2026-02-28 | 3개 리뷰 (PLAN_REVIEW.md, PLAN_REVIEW_2026-02-28.md, REVIEW.md) 통합 개선. P0 전항목 반영: Category A feature leakage 수정 (LFC→normalized counts), LOMO feature selection leakage 경고 추가, mission-stratified split 명시, GLDS-168 제외, composite score 정규화 수식 수정, apeglm shrinkage 확정, H1/H2/H3 예상 발견 가설 추가. D1 ordinal 재설계, Category C 3방법 병렬 설계, 배치보정 Harmony→limma/RUVseq 교체, 마우스 계통 Track 정책 추가 |
| **v0.4** | 2026-02-27 | 모듈식 출판 전략 도입 (v1/v2/v3), 저널 타겟 재조정, H 카테고리 v3.0으로 이동 |
| **v0.3** | 2026-02-27 | 방사선(RadLab)·환경(EDA)·미생물(Category G) 데이터 포함, 미션 메타데이터 JSON 스키마 추가 |
| **v0.2** | 2026-02-27 | 크로스미션(Category B)·크로스티슈(Category C) task 설계, control_type 두 트랙 분리 |
| **v0.1** | 2026-02-27 | 초기 계획 수립 (NASA OSDR 데이터 범위, 차별화 포인트) |

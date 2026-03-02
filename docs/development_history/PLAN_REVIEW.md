# GeneLab_benchmark PLAN.md — 상세 리뷰

> 리뷰 날짜: 2026-02-27
> 리뷰어: Claude Sonnet 4.6 (AI assistant)
> 리뷰 대상: PLAN.md v0.4

---

## 총평

**전체 평가**: 매우 잘 구성된 과학 벤치마크 계획서. 아래 사항을 수정하면 Genome Biology 제출 수준의 벤치마크가 될 수 있음.

**특히 우수한 점**:
- 대조군 유형(GC/VC/BC) 명시적 구분 — 우주생물학 문헌에서 흔히 무시되는 핵심 포인트
- LFC feature를 cross-mission 표현으로 선택한 근거가 명확
- Confounder task(D4, D5)를 명시적으로 설계 — 정직한 과학적 접근
- LOMO full 6×6 cross-mission matrix(B1) — 벤치마크의 핵심 기여

**리뷰 구조**:
1. Critical Issues (반드시 수정)
2. Important Improvements (권장 수정)
3. Minor Suggestions
4. Missing Sections
5. Task Count 검증

---

## 1. Critical Issues (반드시 수정)

### 1.1 LFC 계산 방법 미명시 — 결과 재현 불가 위험

**문제**: "LFC (log fold-change vs mission GC/VC)" 라고만 쓰여 있음.
n=3/그룹의 소규모 데이터에서 DESeq2 LFC 추정치는 apeglm/ashr 수축(shrinkage) 없이는 **분산이 매우 큼**.

```
RR-1 kidney (n=3/group) 예시:
- 수축 없는 LFC: 일부 유전자가 ±10 이상 (통계적 노이즈)
- apeglm 수축 LFC: 작은 n에서 LFC를 신호 크기에 비례해 수축
```

**필수 명시 사항**:
```python
# preprocess_mouse_rnaseq.py에 명시할 것
LFC_METHOD = "apeglm"          # "apeglm" | "ashr" | "none" (plain log2 ratio)
LFC_REFERENCE = "GC"           # "GC" | "VC" (Track별로 다름)
LFC_PSEUDOCOUNT = 0.5          # counts + 0.5 before log2 (for Track 2 only)
```

**대응**: `DESIGN_DECISIONS.md`에 apeglm vs plain LFC 비교 테이블 추가.
J 카테고리에 J0: "LFC shrinkage 방법 비교" task 추가 검토.

---

### 1.2 Category D1 (비행 기간 회귀) — 설계 결함

**문제**: 비행 기간이 실질적으로 4~5개 이산 클러스터임.

```
실제 OSDR 마우스 미션 기간 분포:
- ~13일: STS-135 (A2, B2에 포함)
- ~30-39일: RR-1 (37d), RR-3 (33d), RR-6 (39d), RR-8 (30d 내외)
- ~75일: RR-7
- ~92일: ISS 장기 (아직 데이터 제한적)
```

이 5개 포인트로 "연속 회귀"를 정당화하기 어려움.

**대응**:
```
Option A (권장): D1을 "5-class ordinal classification"으로 재설계
  → 13d / 22-30d / 31-40d / 41-75d / >75d

Option B: D1 유지하되 논문에서 명시
  → "5-point dose-response curve (not true continuous)"
  → Pearson r 대신 Spearman r을 primary metric으로
  → 비선형 모델(polynomial, spline) baseline 추가

Option C: D1 제거하고 D2 (3-class) 만 유지
```

---

### 1.3 Feature Selection 데이터 누출 위험 (Category B/C)

**문제**: 3.5절의 feature selection 정책이 **전체 데이터셋 기준**으로 되어 있음.
Cross-mission (B) 또는 Cross-tissue (C)에서 holdout mission의 정보가 feature selection에 새어 들어갈 수 있음.

```python
# 잘못된 방법 (현재 암묵적 설계):
variance_filtered_genes = compute_variance_filter(ALL_missions_combined)
# → holdout mission 포함 → 데이터 누출

# 올바른 방법:
for train_missions, test_mission in lomo_splits:
    variance_filtered_genes = compute_variance_filter(train_missions_only)
    # holdout mission에서 같은 유전자 세트 사용
    model.fit(train_missions[filtered_genes])
    model.predict(test_mission[filtered_genes])
```

**대응**: `generate_tasks.py`에 feature selection을 LOMO 루프 내부로 명시적으로 이동.
`DESIGN_DECISIONS.md`에 누출 방지 정책 섹션 추가.

---

### 1.4 마우스 계통(Strain) 혼합이 A/B/C를 오염시킬 수 있음

**문제**: D4에서 계통을 confounder로 측정하지만, A/B/C task에서는 계통을 통제하지 않음.

```
GLDS별 마우스 계통:
- RR-1 (GLDS-48 liver): C57BL/6J, Female, 10주령
- MHU-2 (GLDS-617 liver): C57BL/6J 또는 다른 계통? → 확인 필요
- STS-135: BALB/c Female (B2 task에 포함 예정)
```

STS-135 (BALB/c) → ISS missions (C57BL/6J) 전이(B2)는 계통 차이로 인해
spaceflight 신호와 계통 신호가 완전히 뒤엉킴.

**대응 옵션**:
```
Option A: A/B/C를 C57BL/6J만으로 한정 (STS-135는 B2에서 제외)
Option B: 계통 정보를 covariate로 모델에 포함 (e.g., linear mixed model)
Option C: B2 유지하되 "이 task에는 strain confound 있음"을 명시,
          D4와 B2 AUROC를 함께 해석하도록 논문 작성
```

Option C가 실제로 **더 흥미로운 발견**을 줄 수 있음 (strain effect vs spaceflight effect 상대적 크기).

---

### 1.5 NASA AWG 데이터 접근성 — 공개 벤치마크 조건 충족 여부

**문제**: "내부 데이터 접근"이 언급되어 있으나, 공개 벤치마크는 **모든 데이터가 공개 접근 가능**해야 함.

**확인 필요 사항**:
1. 리스트된 모든 GLDS ID가 OSDR에서 누구나 다운로드 가능한가?
2. 특히 최신 GLDS (GLDS-617, 638, 664, 674, 689)가 이미 공개 릴리즈됐는가?
3. AWG와 공동연구한 "내부" 분석 결과가 논문에 들어가는 경우, 해당 데이터도 공개돼야 함.

**대응**: `DATA_INVENTORY.md`에 각 GLDS의 공개 상태 (Public / Embargo / Restricted) 컬럼 추가.

---

## 2. Important Improvements (권장 수정)

### 2.1 Batch Correction 방법론 업데이트 필요

**현재 계획**: ComBat-seq vs Harmony 비교 (J3)

**문제**: Harmony는 주로 단일세포 데이터용. Bulk RNA-seq 배치보정에는:
- `ComBat-seq` (R, sva 패키지) ← 적절
- `limma::removeBatchEffect` ← 적절
- `RUVseq` ← 특히 부정적 대조군 있을 때
- `Harmony` ← 단일세포에 적합 (bulk에도 쓰이지만 standard는 아님)

**대응**: J3를 아래로 개정:
```
J3: 배치보정 없음 vs ComBat-seq vs limma::removeBatchEffect vs RUVseq
    + Harmony는 "단일세포용이지만 비교용으로 포함" 명시
```

---

### 2.2 Small n 환경의 ML 오버피팅 위험

**문제**: n=30 training samples, p=5,000~8,000 features → n << p.
현재 baselines (LogReg, RF, XGBoost, LightGBM)는 모두 과적합 위험.

**필수 추가**:
```python
# 강한 regularization 명시
logistic_regression_params = {
    "penalty": "l1",   # lasso (feature selection 효과)
    "C": [0.001, 0.01, 0.1],  # cross-validated
    "solver": "liblinear",
}

# 또는 Elastic Net
# XGBoost/LightGBM: max_depth=2~3, subsample 강화

# PCA 전처리 option 추가 (50~100 PCs)
# → raw features vs PCA(50) vs PCA(100) 비교
```

**추가 baseline 제안**:
- PLIER (Pathway Level Information Extractor) — bulk RNA-seq에 자주 쓰임
- Ridge + PCA pipeline

---

### 2.3 Category C Cross-Tissue Transfer — Feature 표준화 문제

**문제**: 조직별 LFC 값의 절대 스케일이 다름.

```
예시:
- 간 LFC: 유전자 A = +3.2 (간에서 강하게 상향)
- 신장 LFC: 유전자 A = +0.8 (신장에서 약하게 상향)
```

LFC를 직접 사용하면 조직 간 스케일 차이가 transfer 성능을 왜곡.

**대응**: C 카테고리에서 두 가지 표현 모두 비교:
```
C_variant_1: LFC 직접 사용 (조직별 스케일 차이 유지)
C_variant_2: gene rank (조직 내 상대 순위) 사용
             → 조직 간 스케일 차이 제거
```
Gene rank가 더 robust할 것으로 예측되지만, 이것이 benchmark 자체의 발견이 될 수 있음.

---

### 2.4 통계적 유의성 검증 체계 없음

**현재 계획**: AUROC, F1 등 점수만 보고. 방법 간 차이의 유의성 판단 기준 없음.

**추가 필요**:
```python
# classification task
from scipy.stats import wilcoxon
# LOMO AUROC 분포 (N_missions 개의 AUROC 값들)
# → paired Wilcoxon test (non-parametric, small n에 적합)
# → Bonferroni 또는 FDR correction (25개 task × 5개 모델)

# Transfer task
# → permutation test on transfer AUROC
# → null: mission 레이블 무작위 섞기 → permutation distribution
```

**최소 요건**: `metrics.py`에 다음 추가:
- Bootstrap CI for AUROC (n_bootstrap=1000)
- Permutation test p-value for each task

---

### 2.5 딥러닝 baseline 누락

**현재 baselines**: LogReg, RF, XGBoost, LightGBM, Random

**제안**: 적어도 다음 중 하나 추가:
```
1. MLP (2-layer, heavy dropout=0.5, L2 regularization)
   → n이 작아 성능은 낮겠지만 "deep learning on small n" 베이스라인 제공
2. AutoML (TabPFN) — small tabular data에 강함, 오믹스에도 최근 적용 사례
3. 1D-CNN on gene rank vector
```

없어도 논문은 되지만, 리뷰어가 "왜 DL을 비교 안 했나?" 물을 것.

---

### 2.6 데이터 공개 및 버전 관리 전략 없음

**현재 계획**: 데이터 저장 경로만 있음. 공개 배포 계획 없음.

**추가 필요**:
```
데이터 배포:
- 방법 A: Zenodo (DOI 획득, 재현성 보장) ← 권장
- 방법 B: HuggingFace datasets (datasets 라이브러리 호환)
- 방법 C: FigShare (용량 제한 있음)

버전 관리:
- benchmark_version: "1.0.0" (semver)
- git tag + Zenodo archive 연동
- task JSON schema version 명시

비교: SpaceOmicsBench가 데이터를 어떻게 배포하는지 참고
(/Users/jak4013/Dropbox/Bioinformatics/Claude/SpaceOmicsBench/v2_public/)
```

---

### 2.7 재현 환경 미명시

**문제**: Python과 R이 모두 필요한 프로젝트이나, 환경 설정이 전혀 명시되지 않음.

**추가 필요**:
```
environment.yml (conda):
  - R >= 4.3.0
    - DESeq2 (Bioconductor)
    - edgeR (Bioconductor)
    - limma (Bioconductor)
    - sva (ComBat-seq)
  - Python >= 3.12
    - scikit-learn
    - lightgbm, xgboost
    - harmonypy (Harmony)
    - pydeseq2 (optional, Python port)

또는: Docker image (Bioconductor/tidyverse + Python)
```

---

## 3. Minor Suggestions (소규모 개선)

### 3.1 GLDS ID 검증 우선순위 조정

계획서에 "첫 번째 작업: catalog_datasets.py"라고 쓰여 있는데, 이보다 먼저 다음을 해야 함:

```python
# Step 0 (catalog_datasets.py 실행 전):
# 이미 검증된 GLDS vs 추정/예상 GLDS 구분
# → 현재 표에서 (?)로 표시된 항목들:
#    - GLDS-638 (MHU-8 soleus)
#    - GLDS-664 (MHU-8 eye)
#    - GLDS-689 (RR-8 skin)
#    - GLDS-674 (kidney)
#    - GLDS-617 (MHU-2 liver)
# → 이 중 일부는 OSDR에 없거나 bulk RNA-seq이 아닐 수 있음
```

**대응**: 표에 "Status" 컬럼 추가 (Verified / Unverified / Assumed).

### 3.2 Task 수 정확성

계획서에 "~25개"라고 했지만 실제 count:

| Category | 확정 Tasks | 조건부 Tasks |
|---------|-----------|------------|
| A | A1~A6 = 6 | A7 (brain, TBD) |
| B | B1~B5 = 5 | B6 (HU, TBD) |
| C | C1~C5 = 5 | - |
| D | D1~D6 = 6 | - |
| J | J1~J4 = 4 | - |
| **합계** | **26** | **2 추가** |

→ 확정 26개, 조건부 포함 28개. "~25개"는 undercount.

### 3.3 A2 (근육) 통합 vs 분리 정책 명확화

현재 Soleus (A2a 보류)와 Gastrocnemius (A2b ✓)가 "A2"라는 같은 번호를 공유.
Soleus와 Gastrocnemius의 spaceflight 반응은 반대 방향일 수 있으므로 절대 합치면 안 됨.

**제안**: A2a와 A2b를 A2, A3으로 re-numbering (기존 A3→A4, 순서 밀기).
또는 A2a를 Supplementary로 완전히 분리.

### 3.4 Category B 배치보정 비교 위치

현재 "B 카테고리 각각에 배치보정 비교 실험" + "→ Category J와 연결"이라고 되어 있음.
이 실험이 B에 있는지 J3에 있는지 혼재. 명확화:

```
J3: 배치보정 없음 vs ComBat-seq vs limma::removeBatchEffect
    평가: B1 transfer AUROC에 대한 효과 (B1이 evaluation 기준)
    → B1 결과 자체는 보정 없음으로 보고
    → J3에서 보정 효과를 secondary 분석으로 보고
```

### 3.5 HU 아날로그 데이터 위치

현재 "3.7 Hindlimb Unloading" 섹션이 있고 B6 task가 있는데,
HU 연구의 GLDS ID가 한 개도 명시되지 않음.

**추가 필요**: HU bulk RNA-seq GLDS 목록 (OSDR에 다수 존재 주장이므로):
```
예상 HU 데이터:
- GLDS-XXX (HU liver, 기간?)
- GLDS-XXX (HU gastrocnemius?)
→ catalog_datasets.py 시 HU 키워드로 검색 필요
```

### 3.6 Changelog 없음

버전 0.4까지 왔지만 무엇이 바뀌었는지 없음. 추가 권장:

```markdown
## Changelog
- v0.4 (2026-02-27): Control type 두 트랙 명시, H 카테고리 v3.0으로 이동, J 카테고리 축소
- v0.3: ...
- v0.2: ...
- v0.1: ...
```

---

## 4. Missing Sections (추가해야 할 섹션)

### 4.1 Power Analysis — "샘플 수가 너무 적다" 비판 대응

섹션 13에서 비판 예상 포인트로 언급했지만, 실제 power analysis가 없음.

**추가 내용**:
```
최소 샘플 수 근거:
- SpaceOmicsBench의 PBMC task: n=40~60으로 AUROC 0.75 달성
- Reference: Buxbaum 2023, power analysis for AUROC in small n
- 우리 데이터: n_flight / n_GC = 5~10 per mission
- Power: α=0.05, power=0.8, expected AUROC=0.70 가정 시 최소 n=?
  (계산: AUROC for small n — Mason & Graham 2002)
- 결론: n≥50 total이 LOMO에서 의미 있는 최소 요건 (근거 인용 필요)
```

### 4.2 Expected Timeline — 언제 무엇이 완성되는가

섹션 10에 Phase 1/2/3가 있지만 달력 기준 타임라인이 없음.

```
Phase 1 (데이터 인프라): 2026-03 ~ 2026-05
  - catalog_datasets.py: 2주
  - 간 7개 GLDS 다운로드 + QC: 3주
  - Task A1 동작 확인: 1주
  - 완료 체크포인트: A1 AUROC > 0.7

Phase 2 (다조직 + 전체 task): 2026-05 ~ 2026-09
  - 나머지 5개 조직: 2개월
  - Category C, D, J 구현: 2개월

Phase 3 (논문): 2026-09 ~ 2027-03
  - 논문 초안: 3개월
  - 공동저자 검토 + 수정: 3개월
  - 제출: 2027-Q1
```

### 4.3 Authorship & Contribution — NASA AWG 공동연구 명시

NASA AWG와 공동연구한다면:
- 누가 제1저자인가?
- AWG 담당자의 역할은?
- 데이터 접근 협약(DUA) 또는 협력 협약이 있는가?
- 논문 제출 전 NASA 내부 승인이 필요한가?

이것이 타임라인에 영향을 줄 수 있음.

### 4.4 Ethical & Legal Considerations

- OSDR 데이터의 재배포 라이선스 조건 확인
- 마우스 실험 데이터는 일반적으로 공개 가능하지만, OSDR 이용약관 확인 필요
- 향후 인간 데이터(v2.0) 포함 시 별도의 data use agreement 필요 가능성

---

## 5. Category별 생물학적 타당성 추가 검토

### A 카테고리 — 생물학적 ground truth 매핑

계획서에 A1 (간), A2 (근육), A3 (신장) 검증 유전자가 있음. 나머지 보완:

| Task | 검증 유전자 세트 | 출처 |
|------|---------------|------|
| A1 Liver | ANGPTL4, PCK1, CYP7A1 (담즙산) | Cell 2020 |
| A2 Muscle | TRIM63/MuRF1, FBXO32/MAFbx | 근위축 문헌 표준 |
| A3 Kidney | Nrf2 경로, Slc34a1 (phosphate reabsorption) | OSD-253 |
| A4 Thymus | T세포 발달 유전자 (CD3e, Rag1, IL-7R) | 면역 표준 |
| A5 Skin | 콜라겐/MMP 경로 (UV/스트레스 반응) | 확인 필요 |
| A6 Eye | VEGFA, PEDF (망막 혈관), 산화스트레스 | 우주 안과 문헌 |

A5 (Skin)의 spaceflight 검증 유전자가 없음 → `BIOLOGICAL_GROUND_TRUTH.md`에 문헌 조사 후 추가 필요.

### C 카테고리 — 예상 AUROC 근거 없음

현재 계획:
- C1 (Liver→Kidney): 0.60~0.70
- C2 (Liver→Gastrocnemius): 0.55~0.65

이 예상치가 어디서 왔는지 불명확. Phase 1 완료 후 파일럿 데이터로 보정해야 하지만,
논문에서 이 예상치를 hypothesize로 명시해야 함.

**제안**: C 카테고리 도입부에 다음 추가:
```
예상 AUROC 근거:
- Pathway overlap (GSEA NES correlation): liver-kidney DEG에서 산화스트레스 pathway 공유 (Zhang 2023)
- 유전자 수준 Jaccard: top-200 DEG에서 15~25% 예상 (비공식 추정)
→ Phase 1 파일럿에서 실제 Jaccard 측정 후 보정
```

---

## 6. 최종 체크리스트 (우선순위 순)

### 즉시 수정 (착수 전)
- [ ] **LFC shrinkage 방법 명시** (apeglm vs none, DESIGN_DECISIONS.md)
- [ ] **GLDS ID 전수 검증** (특히 MHU-8 데이터 존재 여부)
- [ ] **Feature selection 누출 방지** (generate_tasks.py 설계)
- [ ] **마우스 계통 처리 방침** (STS-135 BALB/c 처리)
- [ ] **모든 GLDS 공개 접근 가능 여부 확인**

### Phase 1 완료 전 수정
- [ ] D1 회귀 task 재설계 (이산 vs 연속 논의)
- [ ] ComBat-seq + limma::removeBatchEffect로 J3 개정
- [ ] Power analysis 추가
- [ ] Bootstrap CI + permutation test 설계
- [ ] GLDS status 컬럼 추가 (Verified/Unverified)

### Phase 2 착수 전 수정
- [ ] 데이터 공개 플랫폼 결정 (Zenodo 강력 권장)
- [ ] conda environment.yml 또는 Docker 작성
- [ ] HU 아날로그 GLDS 목록 확보
- [ ] Authorship & NASA 내부 승인 일정 확인

### 논문 작성 전 확인
- [ ] Changelog 추가
- [ ] Task re-numbering (A2a/A2b 정리)
- [ ] Expected AUROC 예상 근거 문헌화
- [ ] BIOLOGICAL_GROUND_TRUTH.md (A5 Skin 포함 전체 완성)

---

## 7. 전체 평가 요약

| 측면 | 평가 | 비고 |
|------|------|------|
| 과학적 독창성 | ★★★★★ | GC/VC 구분, LOMO 크로스미션, confounder task |
| 방법론 엄밀성 | ★★★☆☆ | LFC 방법 미명시, feature selection 누출, strain 혼합 |
| 실행 가능성 | ★★★★☆ | 대부분 공개 데이터, 기존 도구 재활용 |
| 출판 가능성 | ★★★★☆ | Genome Biology 적합, 위 수정 후 |
| 재현성 | ★★★☆☆ | 환경 명시 필요, 데이터 배포 계획 없음 |
| 커뮤니티 기여 | ★★★★☆ | SpaceOmicsBench 연동 시 강화됨 |

**최우선 액션**: `catalog_datasets.py` 실행 전에 LFC 방법과 feature selection 누출 방지 설계를 확정할 것.
데이터 품질(실제 GLDS 내용)이 확인되기 전까지는 task 수와 예상 AUROC를 잠정치로만 처리.

# GeneLab_benchmark — PLAN 개선 계획
**근거**: PLAN_REVIEW_2026-02-28.md + PLAN_REVIEW.md + REVIEW.md 종합
**작성일**: 2026-02-28
**적용 대상**: PLAN.md v0.4 → v1.0

---

## 개선 계획 개요

세 리뷰 파일에서 공통으로 지적된 문제와 각각이 지적한 문제를 통합하여
**6개 테마 × 우선순위** 구조로 정리한다.

```
테마 1: 데이터 누출(Leakage) — 논문 신뢰성의 핵심
테마 2: 통계적 엄밀성 — 저널 수락의 최소 기준
테마 3: 방법론 선택 — 방어 가능한 설계
테마 4: 과학적 내러티브 — 출판 임팩트
테마 5: 인프라 — 분야 표준 벤치마크의 조건
테마 6: 구현 세부사항 — 재현성과 일관성
```

우선순위:
- **P0**: 코드 착수 전 설계 확정 필수 (구현 중 수정 불가)
- **P1**: Phase 1 완료 전 반영 필수 (논문 결과에 영향)
- **P2**: Phase 2 착수 전 반영 (논문 방법론 섹션 필요)
- **P3**: 논문 작성 전 반영 (출판 품질)
- **P4**: v1.0 이후 (living benchmark 관련)

---

## 테마 1: 데이터 누출(Leakage) 제거

*세 리뷰 모두 가장 치명적 문제로 지적. P0 전부.*

### 1-A. Category A feature leakage 수정 [P0]

**현재**: "LFC (log fold-change vs mission GC/VC)"를 Category A 입력으로 사용
**문제**: LFC는 flight/ground 그룹 레이블을 알아야 계산 가능 → 분류 목표 사전 사용
**3개 리뷰 일치**: PLAN_REVIEW_2026-02-28 #1, PLAN_REVIEW #1.2, REVIEW Part1-1-2

**확정 결정**:
```
Category A 입력:  log2(DESeq2 size-factor normalized counts) per sample
                 (sample별 정규화만. 그룹 평균 사용 금지)

Category B, C 입력: mission-level LFC
                   (학습 미션에서만 계산. 테스트 미션 leakage 없음)

적용 대상 PLAN.md 섹션: 3.4 Feature Representation 표 전면 교체
```

**PLAN.md 수정 내용**:
```
[현재] Category A: LFC (log fold-change vs mission GC/VC, log2)
[수정] Category A: log2(DESeq2 normalized counts), per-sample
       Mission label은 train/test 분리 후 각각 독립 처리
```

---

### 1-B. LOMO 내 Feature Selection 누출 수정 [P0]

**현재**: feature selection(variance filter)이 전체 데이터 기준으로 적용됨
**문제**: holdout 미션 정보가 feature selection에 노출됨
**근거**: PLAN_REVIEW #1.3, PLAN_REVIEW_2026-02-28 #2

**확정 결정**:
```python
# 잘못된 방법 (현재):
filtered_genes = variance_filter(ALL_missions)  # holdout 포함

# 올바른 방법 (수정 후):
for train_missions, test_mission in lomo_splits:
    filtered_genes = variance_filter(train_missions_only)  # holdout 제외
    train_X = train_data[filtered_genes]
    test_X  = test_data[filtered_genes]   # 같은 유전자 세트 적용
```

**PLAN.md 수정 내용**: 섹션 3.5 Feature Selection에 아래 주의사항 추가
```
[추가] ⚠️ Feature selection은 반드시 LOMO 루프 내부에서 수행.
       전체 데이터 기준 필터링은 holdout 미션 정보 노출로 결과 과대평가를 유발.
       generate_tasks.py에서 split_aware_feature_filter() 함수로 구현.
```

---

### 1-C. Animal/Cage 수준 독립성 보장 [P0]

**현재**: 샘플 독립성 단위 미명시
**문제**: 동일 케이지/그룹 마우스의 샘플이 train/test에 분리될 경우 성능 과대 추정
**근거**: PLAN_REVIEW_2026-02-28 #3

**확정 결정**:
```
독립성 단위: mission (미션 전체가 train 또는 test)
             → LOMO는 이미 mission 단위 분리이므로 구조적으로 cage-level leak 없음

단, Category D 80/20 split에서는:
  - 동일 mission의 샘플이 train/test에 섞이지 않도록
  - mission-stratified split 사용 (랜덤 샘플 분리 금지)

Category J에서는:
  - 동일 GLDS 샘플을 train과 test에 분리 불가
  - GLDS 단위 split 사용
```

**PLAN.md 수정 내용**: 섹션 6.1 Split 전략 표에 "독립성 단위" 컬럼 추가
```
| Split 유형 | 독립성 단위 | 적용 Task |
|---|---|---|
| LOMO | Mission (구조적으로 cage-level 보장) | A, B |
| LOTO | Tissue | C |
| Mission-stratified 80/20 | Mission (같은 미션은 같은 fold) | D |
| GLDS-stratified | GLDS study | J |
```

---

### 1-D. GLDS-168 중복 샘플 처리 확정 [P0]

**현재**: A1 Liver에 GLDS-48(RR-1), GLDS-137(RR-3), GLDS-168(RR-1+RR-3) 동시 포함
**문제**: GLDS-168이 GLDS-48, 137의 샘플을 재처리한 것 → 동일 샘플 중복
**근거**: REVIEW Part1-1-3

**확정 결정**:
```
GLDS-168: Category A1에서 완전 제외
           Category J에만 포함 (동일 샘플의 단독 vs 통합 처리 파이프라인 비교)

A1 Liver 미션 구성 (수정 후):
  GLDS-48  (RR-1)
  GLDS-137 (RR-3)
  GLDS-245 (RR-6)
  GLDS-379 (RR-8)
  GLDS-242 (RR-9)
  GLDS-617 (MHU-2)
  GLDS-25  (STS-135)  ← 별도 확인 필요 (BALB/c 계통 주의)
```

**PLAN.md 수정 내용**: 섹션 4.1 표에서 GLDS-168 제거, 주석 추가
```
[추가 주석] GLDS-168은 GLDS-48+137 통합 재처리 데이터.
            Category A에서 제외하여 샘플 중복 방지.
            Category J (J1)에서 파이프라인 비교용으로만 활용.
```

---

### 1-E. D4/D5 Confounder Task의 해석적 취약성 [P1]

**문제**: 계통(D4)·하드웨어(D5)가 미션/기간과 공선성(collinearity)이 있을 경우,
높은 예측 성능이 isolated confounder 효과인지 미션 전체 차이인지 구분 불가
**근거**: PLAN_REVIEW_2026-02-28 #6

**수정 방향**:
```
D4, D5를 "confounder 정량화"로만 쓰되, 해석 한계를 명시:

논문 Methods에 추가:
  "D4(strain)와 D5(hardware)의 예측 성능은 해당 변수의
   독립적 효과가 아닌, 관련 미션 배치 차이의 상한선으로 해석해야 한다.
   예: BALB/c는 RR-3에서만 사용되었으므로 D4 AUROC = D3 AUROC의 subset."

추가 분석: partial correlation (미션 효과 residualize 후 계통 효과 재계산)
```

**PLAN.md 수정 내용**: D4, D5 task 설명에 "해석 주의사항" 추가

---

## 테마 2: 통계적 엄밀성

*저널 리뷰어의 최소 요건. P0~P1.*

### 2-A. 불확실성 정량화 구현 계획 [P0]

**현재**: Point estimate만 보고
**문제**: AUROC 0.73 vs 0.71의 유의미한 차이 주장 불가
**3개 리뷰 일치**: PLAN_REVIEW_2026-02-28 #7, PLAN_REVIEW #2.4, REVIEW Part1-1-5

**확정 결정**:
```python
# metrics.py에 추가할 함수들:

def auroc_with_ci(y_true, y_score, n_bootstrap=1000, alpha=0.05):
    """Bootstrap 95% CI for AUROC"""

def delong_test(y_true, y_score_a, y_score_b):
    """DeLong's test for comparing two AUROCs"""

def permutation_test(y_true, y_score, n_permutations=1000):
    """Permutation p-value: 랜덤 레이블로 얻을 AUROC와 비교"""

def lomo_summary(fold_aurocs):
    """LOMO fold별 AUROC → mean ± SD"""
```

**보고 형식 표준화**:
```
AUROC = 0.81 [95% CI: 0.74–0.87], permutation p < 0.001
(mean ± SD across LOMO folds: 0.81 ± 0.06)
```

**PLAN.md 수정 내용**: 섹션 6.2 평가 지표 표에 "CI 방법" 컬럼 추가
```
| Task 유형 | Primary | CI 방법 | 모델 간 비교 |
|---|---|---|---|
| 이진 분류 | AUROC | Bootstrap 95% CI | DeLong's test |
| 다중 클래스 | macro-F1 | Bootstrap 95% CI | Wilcoxon signed-rank |
| 회귀 | Spearman r | Bootstrap 95% CI | - |
```

---

### 2-B. Composite Score 정규화 수식 오류 수정 [P0]

**현재**: 다중 클래스 primary = macro-F1이지만 random baseline = 1/K accuracy로 기술
**문제**: metric과 baseline이 다른 스케일 → 정규화 수식이 수학적으로 잘못됨
**근거**: PLAN_REVIEW_2026-02-28 #4

**확정 결정**:
```python
# 올바른 정규화:
def random_baseline(task_type, n_classes=None, class_weights=None):
    if task_type == "binary":
        return 0.5                    # AUROC random = 0.5
    elif task_type == "multiclass":
        # macro-F1 random baseline (클래스 균형 가정):
        return 1.0 / n_classes        # ≠ 1/K accuracy
    elif task_type == "regression":
        return 0.0                    # Spearman r random ≈ 0

def normalize_score(score, task_type, **kwargs):
    baseline = random_baseline(task_type, **kwargs)
    return (score - baseline) / (1.0 - baseline)
```

**PLAN.md 수정 내용**: 섹션 6.3 수식을 task_type별로 분리 명시

---

### 2-C. 다중 검정 보정 정책 확정 [P1]

**현재**: 없음
**근거**: PLAN_REVIEW #2.4, REVIEW Part3-3-3

**확정 결정**:
```
Level 1 (task 내 모델 비교):
  DeLong's test for AUROC, n_comparisons = n_models × (n_models-1) / 2
  Bonferroni correction within task

Level 2 (task 간 비교):
  Benjamini-Hochberg FDR correction across all tasks
  q < 0.05 threshold

Level 3 (composite score):
  Bootstrap CI로 95% CI 제시
  "유의미하게 높다"는 주장 = CI 비중첩
```

**PLAN.md 수정 내용**: 섹션 6.2에 "다중 검정 보정" 서브섹션 추가

---

### 2-D. Phase 1 체크포인트 기준 개선 [P1]

**현재**: "A1 AUROC > 0.7" (단일 임계값, CI 없음)
**문제**: CI 없는 단일 임계값은 샘플 수 변동에 취약
**근거**: PLAN_REVIEW_2026-02-28 #7

**수정 방향**:
```
Phase 1 체크포인트 (수정):
  [통과] A1 AUROC > 0.7, 95% CI 하한 > 0.6, permutation p < 0.05
  [통과] Permutation null AUROC = 0.50 ± 0.03 (신호 존재 확인)
  [통과] 알려진 spaceflight 유전자 (ANGPTL4, PCK1) SHAP top-50 이내
  [실패 시] 데이터 QC 재검토 → 미션별 샘플 수·mapping rate 확인
```

**PLAN.md 수정 내용**: 섹션 10 Phase 1 체크포인트를 3-조건 AND로 교체

---

## 테마 3: 방법론 선택

*방어 가능한 구현 결정들. P0~P2.*

### 3-A. LFC Shrinkage 방법 확정 [P0]

**현재**: "LFC (log fold-change vs mission GC/VC)"만 기술, 수축 방법 미명시
**문제**: n=3/그룹 소규모 데이터에서 shrinkage 없는 LFC는 분산이 매우 큼
**근거**: PLAN_REVIEW #1.1

**확정 결정**:
```python
# DESIGN_DECISIONS.md에 명시할 것:

LFC_PRIMARY_METHOD = "apeglm"     # DESeq2 apeglm shrinkage (기본)
LFC_COMPARISON_METHOD = "ashr"    # 비교용 (J 카테고리)
LFC_PSEUDOCOUNT = None            # apeglm 사용 시 불필요

# Category A에서는 LFC 사용 안 함 (1-A 수정 결과)
# Category B, C에서만 shrinkage LFC 사용
```

**PLAN.md 수정 내용**: 섹션 3.4 Feature Representation에 shrinkage 방법 추가

---

### 3-B. 배치 보정 방법 교체 [P1]

**현재**: "ComBat-seq vs Harmony" 비교
**문제**: Harmony는 single-cell 임베딩용. Bulk RNA-seq 표준이 아님
**근거**: PLAN_REVIEW_2026-02-28 #5, PLAN_REVIEW #2.1

**확정 결정**:
```
J3 배치 보정 비교 (수정):
  비교 대상 1: 보정 없음 (baseline)
  비교 대상 2: ComBat-seq (sva 패키지, negative binomial 모델)
  비교 대상 3: limma::removeBatchEffect (log-CPM 입력)
  비교 대상 4: RUVseq (negative control genes 있을 경우)
  참고용 비교: Harmony (단일세포용, "bulk에 적용 시 비교" 명시)

평가: Category B1 transfer AUROC에 대한 각 방법의 효과
```

**PLAN.md 수정 내용**: 섹션 5 Category J의 J3 행 교체

---

### 3-C. D1 비행 기간 회귀 재설계 [P1]

**현재**: D1 = 연속 회귀 (Pearson r)
**문제**: 실제 기간 포인트가 4~5개 클러스터 → 연속 회귀 의미 없음
**근거**: PLAN_REVIEW #1.2, PLAN_REVIEW_2026-02-28 (간접), REVIEW Part3-3-2

**확정 결정**:
```
D1 (수정): 5-class ordinal classification으로 재설계
  Class 1: ≤14일  (STS-135: 13일)
  Class 2: 15~25일 (해당 미션 없으면 제거)
  Class 3: 26~40일 (RR-1~3,6,8,9: 30~39일)
  Class 4: 41~74일 (해당 미션 없으면 3-class로 축소)
  Class 5: ≥75일  (RR-7: 75일)

실제 가능한 클래스: 3-class (≤14d / 26~40d / ≥75d)
Primary metric: macro-F1 (ordinal 특성 반영: Spearman r between predicted/true class)
Random baseline: 1/3 macro-F1

D1을 이렇게 재정의하면 D2(기존 3-class)와 차별화:
  D2: 비행 기간에 따른 분자 반응 존재 여부 (단기/중기/장기 분류)
  D1(수정): 비행 기간의 정밀 예측 (5-class ordinal)
```

**PLAN.md 수정 내용**: 섹션 5 Category D의 D1 행 및 설명 전면 교체

---

### 3-D. Category C Cross-Tissue 방법 명확화 [P1]

**현재**: "한 조직으로 학습 → 다른 조직으로 테스트"가 두 가지 해석 가능
**문제**: Covariate shift 방법(A)은 조직별 발현 차이로 실패 보장
**근거**: PLAN_REVIEW #2.3, REVIEW Part3-3-1

**확정 결정**:
```
Category C의 3가지 방법을 병렬 실행하여 비교:

방법 A (Gene-level covariate shift):
  Liver 학습 분류기 가중치 → Kidney에 직접 적용
  예상: AUROC 낮음 (조직 특이 유전자 문제)
  역할: 하한선(lower bound) 확인

방법 B (DEG feature transfer):
  Liver top-200 DEG를 feature로 선택 → Kidney에서 동일 유전자로 새 모델 학습
  예상: 방법 A보다 높음
  역할: 중간 baseline

방법 C (Pathway-level transfer) [PRIMARY]:
  Liver에서 GSEA 실행 → top-20 enriched pathway 선택
  Kidney에서 같은 pathway 점수 계산 → 분류
  예상: 가장 높음 (조직 불문 pathway 보존)
  역할: 생물학적으로 가장 의미 있음

C의 발견 스토리:
  "방법 A → B → C 순으로 성능 증가"는
  "유전자 수준보다 경로 수준에서 spaceflight 반응이 더 보존된다"는 발견이 됨
```

**PLAN.md 수정 내용**: 섹션 5 Category C에 3가지 방법 병렬 설계 명시

---

### 3-E. 마우스 계통 처리 방침 확정 [P1]

**현재**: STS-135 (BALB/c)를 Category B2에 포함, 계통 통제 없음
**문제**: BALB/c vs C57BL/6J 차이가 spaceflight 신호보다 클 수 있음
**근거**: PLAN_REVIEW #1.4, PLAN_REVIEW_2026-02-28 #6

**확정 결정**:
```
정책: "계통 혼합은 허용하되, 명시적으로 레이블링하고 분석"

Track 2a (C57BL/6J only): Category A, B, C의 primary 분석
Track 2b (모든 계통 포함): Category A, B, C의 secondary 분석

B2 (ISS → STS-135) 특별 처리:
  - "cross-strain" 실험으로 명시적으로 reframe
  - B2 결과 = spaceflight 신호 + strain 신호 합산의 상한선
  - D4 (strain prediction AUROC)로 strain 효과 크기 보정 후 해석

PLAN.md에 각 GLDS의 mouse_strain을 missions.json 필수 필드로 명시
```

**PLAN.md 수정 내용**: 섹션 3 핵심 설계 원칙에 "3.8 마우스 계통 처리 정책" 서브섹션 추가

---

### 3-F. 소규모 N 환경 ML 정규화 정책 [P1]

**현재**: regularization 파라미터 미명시
**문제**: n ≈ 30~60, p ≈ 5,000~8,000 → n << p → 과적합 위험
**근거**: PLAN_REVIEW #2.2

**확정 결정**:
```python
# baselines/run_baselines.py에 명시:

BASELINE_HYPERPARAMS = {
    "logistic_regression": {
        "penalty": "elasticnet",
        "l1_ratio": [0.1, 0.5, 0.9],
        "C": [0.001, 0.01, 0.1],
        "solver": "saga",
        # cross-validated within training fold
    },
    "random_forest": {
        "n_estimators": 100,
        "max_features": "sqrt",  # 유전자 수의 sqrt
        "max_depth": [3, 5, None],
        "class_weight": "balanced",
    },
    "xgboost": {
        "max_depth": [2, 3],       # 의도적으로 얕게 (과적합 방지)
        "subsample": 0.7,
        "colsample_bytree": 0.3,
        "n_estimators": 100,
    },
    "lightgbm": {
        "num_leaves": [15, 31],    # 기본값 31 유지 또는 15로 축소
        "subsample": 0.7,
        "colsample_bytree": 0.3,
    }
}

# 추가 baseline:
# PCA (50 components) + LogReg → "low-dim representation" baseline
# Ridge regression (multi-output) for D1
```

---

## 테마 4: 과학적 내러티브 수립

*출판 임팩트를 결정하는 핵심. 코드 전에 확정.*

### 4-A. 논문 핵심 발견 스토리 확정 [P0]

**현재**: 없음
**문제**: "벤치마크 만들었습니다"만으로는 상위 저널 수락 불가
**근거**: REVIEW Part1-1-1

**확정 가설 (데이터 탐색 후 검증)**:

```
Primary Hypothesis (H1):
  "마우스 우주비행 전사체 서명의 미션 간 일반화 능력은
   조직별로 유의미하게 다르다.
   간(liver)은 가장 일관된 서명을 보이고,
   근육(muscle)은 미션 특이적 반응이 지배적이다."

Secondary Hypothesis (H2):
  "ML 모델의 cross-mission 전이 실패는
   배치효과보다 생물학적 다양성에 기인한다.
   (ComBat-seq 보정 후 B1 transfer AUROC 개선 < 5%)"

Tertiary Hypothesis (H3):
  "우주비행 전사체 반응은 유전자 수준보다
   경로(pathway) 수준에서 더 잘 보존된다.
   (Category C: 방법 A < 방법 B << 방법 C)"

Expected Figure Structure:
  Figure 1: LOMO 6×6 cross-mission AUROC matrix (B1)
            + 미션별 방사선량·기간·계통 metadata 연결
  Figure 2: Tissue-specific classification performance (A1~A6) 비교
            + GC track vs VC track 차이
  Figure 3: Category C pathway-level transfer vs gene-level
  Figure 4: Confounder 정량화 (D4: strain, D5: hardware vs A1 spaceflight AUROC)
```

**PLAN.md 수정 내용**: 섹션 1 비전에 "Expected Findings" 서브섹션 추가

---

### 4-B. 생물학적 Ground Truth 검증 체계 강화 [P1] — 🔄 부분 완료

**현재**: feature importance 기반 sanity check (불안정)
**문제**: 모델 의존적, 재현 어려움
**근거**: REVIEW Part3-3-4, PLAN_REVIEW #5

**확정 결정**:
```
생물학적 검증 3단계:

Stage 1: GSEA on model features  ← ✅ 부분 완료 (fGSEA 파이프라인으로 대체)
  top-100 유전자 (SHAP values 기준) → GSEA
  기대: FOXO, Nrf2, mitochondria pathway 농축
  실패 시: 모델이 생물학적으로 의미 없는 feature 학습 의심

Stage 2: Direction concordance  ← ✅ 완료 (contrast direction 검증)
  Cell 2020 마우스 13조직 DEG 방향 (up/down) vs 이 벤치마크 LFC 방향
  기대: 공통 spaceflight 유전자의 방향이 80%+ 일치
  지표: Direction concordance score (%) per tissue

Stage 3: External DEG overlap  ← ⏳ 미완 (Cell 2020 데이터 필요)
  Cell 2020 top-500 DEG ∩ 이 벤치마크 top-500 feature
  기대: Jaccard index > 0.20 (>10% 기대 시 유의미)
  계층: liver > kidney ≈ thymus > muscle (예측)

BIOLOGICAL_GROUND_TRUTH.md:
  각 조직별 검증 유전자 세트 + 출처 논문 테이블
  A5 Skin: 콜라겐 리모델링(Col1a1, Col3a1), MMP (Mmp2, Mmp9) → 확인 필요
```

**✅ 구현 상태 (2026-03-01)**:
- **Stage 1 대체**: fGSEA 결과로 5 tissues 생물학적 검증 완료 (DD-15 참조)
  - Liver: OXIDATIVE_PHOSPHORYLATION, FATTY_ACID_METABOLISM enriched ✓
  - Thymus: E2F_TARGETS, G2M_CHECKPOINT, IFN-gamma ✓
  - Gastrocnemius: OXIDATIVE_PHOSPHORYLATION, MYOGENESIS ✓
  - Kidney: MTORC1_SIGNALING ✓, Eye: OXIDATIVE_PHOSPHORYLATION dominant ✓
- **Stage 2**: Contrast direction (GC-first) 검증 완료 — Pck1, Angptl4, Cyp2e1 방향 확인
- **Stage 3**: Cell 2020 외부 DEG 비교 미완 → Phase 3 (문서화) 예정
- **SHAP→fGSEA**: SHAP top-100 genes → pathway enrichment는 미구현 (SHAP 분석 자체 미실행)

---

### 4-C. Negative Control Tasks 추가 [P1]

**현재**: 없음
**근거**: REVIEW Part1-1-7

**확정 결정**:
```
NC1: Permutation Test (모든 task에 적용)
  y_true를 1000회 무작위 섞기 → AUROC 분포 생성
  실제 AUROC가 permutation 분포 95th percentile보다 높아야 유의미
  → metrics.py의 permutation_test() 함수로 자동화

NC2: Housekeeping Gene Control (A 카테고리)
  우주비행 무관 하우스키핑 유전자 50개 (GAPDH, ACTB 등)만 사용
  기대 AUROC ≈ 0.50 → 실제 signal의 배경 노이즈 확인

NC3: Cross-species Alien Control (C 카테고리)
  Liver 학습 → Arabidopsis 유전자로 테스트 (ortholog 없음)
  기대 Transfer AUROC ≈ 0.50 → 실제 C1~C5의 상대적 평가 기준
```

**PLAN.md 수정 내용**: 섹션 5 각 Category에 NC task 1개씩 추가 (supplementary로 분류)

---

## 테마 5: 인프라 구축

*분야 표준 벤치마크의 조건. P2~P4.*

### 5-A. Held-out Test Set 설계 [P2]

**현재**: 없음
**근거**: REVIEW Part2-2-2

**확정 결정**:
```
Held-out 방법: 최신 미션 1개를 완전 holdout
  후보: RRRM-2 간 데이터 (또는 catalog_datasets.py 확인 후 가장 최신 미션)
  역할: 어떤 분석에도 사용 안 함 → 최종 리더보드 점수만 이 미션으로 계산

개발 셋 (public):
  나머지 미션 전체 → LOMO split 공개

테스트 셋 (hidden):
  Holdout 미션의 유전자 발현 X만 공개 (레이블 y 비공개)
  제출 시스템에 ŷ만 제출 → 서버에서 AUROC 계산

구현 방법 (최소 비용):
  HuggingFace Spaces + Gradio interface
  또는 간단한 GitHub Actions + evaluation script
```

**PLAN.md 수정 내용**: 섹션 8 (신규) "Held-out Test Set 및 제출 시스템"에 추가

---

### 5-B. Benchmark 버전 동결 및 데이터 스냅샷 [P2]

**현재**: 없음
**근거**: REVIEW Part2-2-4, PLAN_REVIEW #2.6

**확정 결정**:
```
버전 정책:
  v1.0 (출판 시): GLDS 목록 + processed CSV + task JSON 영구 동결
  → Zenodo에 업로드 (DOI 취득)
  → GitHub에 git tag v1.0

업데이트 정책:
  v1.0.x: 버그 수정만 (결과에 영향 없음)
  v1.1: 새 미션 추가 (새 task 번호 부여, 이전 task 결과 유지)
  v2.0: 새 카테고리 (E+F+G)

파일 형식:
  processed CSV: SHA-256 checksum 기록 (DATA_INVENTORY.md)
  task JSON: schema version 필드 포함
  metadata: benchmark_version: "1.0.0" (semver)
```

**PLAN.md 수정 내용**: 섹션 9 디렉토리 구조에 `CHECKSUMS.md` 추가,
섹션 14 참고자료에 버전 정책 링크 추가

---

### 5-C. Foundation Model Baseline 추가 [P2]

**현재**: LogReg, RF, XGBoost, LightGBM, Random
**문제**: 2026년 기준 SOTA 비교 없으면 리뷰어 지적 불가피
**근거**: REVIEW Part2-2-3, PLAN_REVIEW #2.5

**확정 결정**:
```
baselines/에 추가:

1. Geneformer (zero-shot):
   bulk RNA-seq → rank tokenization → Geneformer embedding → linear probe
   구현: theislab/geneformer (HuggingFace)

2. Fine-tuned Geneformer:
   각 Category A task에 대해 fine-tuning
   → "general pre-training이 spaceflight에 전이되는가?" 검증

3. PCA + LogReg (baseline ensemble):
   50 PCs + LogReg → "low-dimensional representation" 효과 확인

4. MLP (small, heavy regularization):
   256→64→1 레이어, dropout=0.5, 배치 정규화
   → "deep learning on small n" reference
```

---

### 5-D. 재현 환경 명세 [P2]

**현재**: 없음
**근거**: PLAN_REVIEW #2.7, PLAN_REVIEW_2026-02-28 #12

**확정 결정**:
```
프로젝트 루트에 추가:

environment.yml (conda):
  name: genelabench
  channels: [conda-forge, bioconda, defaults]
  dependencies:
    - python=3.12
    - r-base=4.3
    - bioconductor-deseq2
    - bioconductor-edger
    - bioconductor-limma
    - bioconductor-sva        # ComBat-seq
    - bioconductor-ruvseq     # RUVseq
    - bioconductor-apeglm     # LFC shrinkage
    - bioconductor-fgsea      # GSEA
    - scikit-learn
    - lightgbm
    - xgboost
    - shap                    # feature importance
    - scipy                   # DeLong's test
    - pydeseq2                # optional Python port

Docker image:
  FROM bioconductor/bioconductor_docker:RELEASE_3_18
  + Python 3.12 + 위 pip 패키지
  → Docker Hub: genelabench/base:1.0
```

---

### 5-E. 공개 리더보드 (최소 버전) [P3]

**현재**: 없음
**근거**: REVIEW Part2-2-1

**확정 결정**:
```
Phase 2 완료 시점에 구축 목표:

최소 구현 (HuggingFace Spaces):
  1. 처리된 데이터: HuggingFace Datasets 업로드
  2. Spaces: 예측 파일(CSV) 업로드 → 자동 AUROC 계산 → 점수 표시
  3. Leaderboard 탭: 제출 기록 + 모델명 + 점수 + 날짜

제출 형식:
  CSV 파일: sample_id, y_pred_proba (0~1)
  JSON 헤더: model_name, description, paper_url

장기 목표: EvalAI 또는 CodaLab 이전
```

---

## 테마 6: 구현 세부사항 개선

*재현성과 일관성. P1~P3.*

### 6-A. PLAN.md 미션 수 불일치 수정 [P1]

**현재**: A1 Liver "7+ 미션", B1 "6×6 matrix" → 불일치
**근거**: PLAN_REVIEW_2026-02-28 #8

**확정 결정**:
```
A1 Liver 미션 목록 동결 (catalog_datasets.py 검증 후):
  [검증 필요] GLDS-48 (RR-1), GLDS-137 (RR-3), GLDS-245 (RR-6),
              GLDS-379 (RR-8), GLDS-242 (RR-9), GLDS-617 (MHU-2), GLDS-25 (STS-135)
  → GLDS-25 (STS-135): BALB/c → Track 2b에만 포함 (Track 2a 제외)
  → B1 6×6 = C57BL/6J 미션만: GLDS-48, 137, 245, 379, 242, 617 (= 6개)
  → B2 = ISS 6개 → STS-135 (BALB/c, 추가)

PLAN.md 해당 섹션에서 "7+" → "6 (Track 2a) + 1 BALB/c (Track 2b)"로 수정
```

---

### 6-B. Track별 독립 포함 기준 명시 [P1]

**현재**: 전체 기준(N_mission ≥ 3, N_total ≥ 50)만 명시
**문제**: Track 1 (GC only)은 전체 기준을 통과해도 샘플 부족 가능
**근거**: PLAN_REVIEW_2026-02-28 #9

**수정 방향**:
```python
def is_eligible(missions, track):
    if track == "track1_gc":
        # GC 보유 미션만 사용
        gc_missions = [m for m in missions if "GC" in m["control_types"]]
        return (len(gc_missions) >= 3 and
                sum(m["n_gc"] for m in gc_missions) >= 30)
    elif track == "track2_vc":
        return (len(missions) >= 3 and
                sum(m["n_vc"] + m["n_flight"] for m in missions) >= 50)
```

**PLAN.md 수정 내용**: 섹션 3.3 Task 포함 최소 기준 표에 Track별 기준 추가

---

### 6-C. 중복 샘플 탐지 방법 개선 [P1]

**현재**: "Pearson r > 0.99" 단일 기준
**문제**: 생물학적으로 유사한 샘플(같은 조직, 같은 미션)도 r > 0.99 가능
**근거**: PLAN_REVIEW_2026-02-28 #10

**수정 방향**:
```python
def detect_duplicates(samples):
    """
    단계적 중복 탐지:
    1. SHA-256 checksum: 완전 동일 파일 탐지
    2. GLDS × sample_id 메타데이터 매칭: 동일 study 내 같은 샘플 탐지
    3. Pearson r > 0.999 (더 엄격) + 같은 GLDS 또는 cross-GLDS
    4. GLDS technical_replicate 필드 확인
    """
```

---

### 6-D. GLDS ID 검증 상태 테이블 추가 [P1]

**현재**: GLDS ID가 "확인된 것"과 "추정한 것"이 섞여 있음
**근거**: PLAN_REVIEW #3.1

**수정 방향**: DATA_CATALOG.md에 Status 컬럼 추가

```markdown
| GLDS ID | Tissue | Mission | Status | Bulk RNA-seq | Notes |
|---|---|---|---|---|---|
| GLDS-48 | Liver | RR-1 | ✅ Verified | Yes | GeneLab processed |
| GLDS-617 | Liver | MHU-2 | ⚠️ Unverified | TBD | catalog_datasets.py로 확인 필요 |
| GLDS-638 | Soleus | MHU-8 | ⚠️ Unverified | TBD | 존재 여부 확인 필요 |
```

---

### 6-E. Task 번호 체계 정리 [P2]

**현재**: A2a/A2b가 같은 번호 공유
**근거**: PLAN_REVIEW #3.3

**확정 결정**:
```
근육 유형 처리:
  Gastrocnemius: A2 (확정)
  Soleus: Phase 1 이후 추가 여부 결정 → A7 또는 Supplementary

현재 A3~A6 번호 유지 (A3 Kidney, A4 Thymus, A5 Skin, A6 Eye)
A7 Brain: 데이터 확인 후 추가

B6 (HU analog): HU GLDS 목록 확보 후 포함 여부 결정
```

---

### 6-F. Changelog 추가 [P2]

**현재**: 없음
**근거**: PLAN_REVIEW #3.6

**PLAN.md에 추가할 섹션**:
```markdown
## Changelog

- **v1.0** (예정, 2027-Q1): 최종 공개 버전
- **v0.5** (2026-02-28): 3개 리뷰 통합 개선. Feature leakage 수정, 통계 체계 추가
- **v0.4** (2026-02-27): 모듈식 출판 전략 도입, 범위 v1.0으로 조정
- **v0.3** (2026-02-27): 방사선·환경·미생물 데이터 포함
- **v0.2** (2026-02-27): 크로스미션·크로스티슈 task 설계
- **v0.1** (2026-02-27): 초기 계획 수립
```

---

### 6-G. HU Analog GLDS 목록 수집 [P1]

**현재**: B6 task 있으나 구체적 GLDS ID 없음
**근거**: PLAN_REVIEW #3.5

**행동**: `catalog_datasets.py`에 HU 키워드 검색 추가:
```python
HU_QUERY_TERMS = [
    "hindlimb unloading",
    "tail suspension",
    "antiorthostatic",
    "HU liver",
    "HU gastrocnemius",
]
# OSDR API: factor=hindlimb+unloading, organism=mus+musculus, assay=rna-seq
```

---

## 통합 수정 우선순위 요약

| 우선순위 | 항목 | 담당 파일 |
|---|---|---|
| **P0-1** | Category A feature: LFC → normalized counts | PLAN.md 3.4 |
| **P0-2** | LOMO 내 feature selection 누출 방지 | PLAN.md 3.5 |
| **P0-3** | Animal/cage 독립성 → mission-stratified split | PLAN.md 6.1 |
| **P0-4** | GLDS-168 제외 → Category J로 이동 | PLAN.md 4.1 |
| **P0-5** | Composite score 정규화 수식 오류 수정 | PLAN.md 6.3 |
| **P0-6** | LFC shrinkage 방법 확정 (apeglm) | PLAN.md 3.4, DESIGN_DECISIONS.md |
| **P0-7** | 논문 핵심 발견 스토리 (H1, H2, H3) | PLAN.md 1절 |
| **P1-1** | Bootstrap CI + DeLong's test 구현 계획 | PLAN.md 6.2 |
| **P1-2** | 다중 검정 보정 정책 | PLAN.md 6.2 |
| **P1-3** | Phase 1 체크포인트 3-조건으로 교체 | PLAN.md 10 |
| **P1-4** | D1 ordinal classification 재설계 | PLAN.md 5 Category D |
| **P1-5** | Category C 3방법 병렬 설계 | PLAN.md 5 Category C |
| **P1-6** | 배치 보정: Harmony → limma/RUVseq | PLAN.md 5 Category J |
| **P1-7** | 마우스 계통 처리 정책 (Track 2a/2b) | PLAN.md 3.8 신규 |
| **P1-8** | ML 정규화 파라미터 명시 | PLAN.md, baselines/ |
| **P1-9** | Negative control tasks (NC1, NC2, NC3) | PLAN.md 5 각 Category |
| **P1-10** | D4/D5 해석 한계 명시 | PLAN.md 5 Category D |
| **P1-11** | GLDS ID 검증 상태 테이블 | DATA_CATALOG.md |
| **P1-12** | HU analog GLDS 목록 수집 | catalog_datasets.py |
| **P2-1** | 생물학적 Ground Truth 3단계 검증 | BIOLOGICAL_GROUND_TRUTH.md |
| **P2-2** | Foundation model baseline (Geneformer) | baselines/ |
| **P2-3** | Held-out test set 설계 | PLAN.md 8 신규 |
| **P2-4** | Benchmark 버전 동결 + Zenodo DOI | PLAN.md, CHECKSUMS.md |
| **P2-5** | conda environment.yml + Docker | environment.yml |
| **P2-6** | Track별 독립 포함 기준 | PLAN.md 3.3 |
| **P2-7** | Task 번호 정리 (A2a/A2b) | PLAN.md 5 |
| **P3-1** | HuggingFace 리더보드 프로토타입 | 별도 |
| **P3-2** | pip 패키지 기초 구조 | genelabench/ |
| **P3-3** | Changelog 추가 | PLAN.md |
| **P4-1** | Living benchmark 절차 문서화 | docs/ |
| **P4-2** | ESA 데이터 포함 여부 결정 | PLAN.md |

---

## 다음 즉시 행동 (P0 → 코드 착수)

```
✅ STEP 1: PLAN.md v0.5 업데이트 (이 개선 계획 적용)        [2026-02-28]
✅ STEP 2: DESIGN_DECISIONS.md 초안 작성                     [2026-02-28]
✅ STEP 3: catalog_datasets.py 착수 + DATA_CATALOG.md 생성   [2026-02-28]
✅ STEP 4: quality_filter.py 대체 → LOMO pipeline 직접 구현  [2026-02-28]
✅ STEP 5: fGSEA + GSVA 파이프라인 (5 tissues × 3 DBs)      [2026-02-28]
✅ STEP 6: Category B (cross-mission transfer, 5 tissues)     [2026-03-01]
✅ STEP 7: Category C (cross-tissue transfer, 4 pairs)        [2026-03-01]
✅ STEP 8: Category D (D3 mission ID + D6 gravity)            [2026-03-01]
✅ STEP 9: J5 gene vs pathway comparison (12 comparisons)     [2026-03-01]
✅ STEP 10: Geneformer tokenization + HPC finetune script     [2026-03-01]

현재 남은 작업:
  [ ] Geneformer HPC fine-tuning 실행 (Cayuga A40)
  [ ] Cell 2020 DEG overlap 비교 (Stage 3, 4-B)
  [ ] BIOLOGICAL_GROUND_TRUTH.md 작성
  [ ] Category J1-J4 (pipeline comparison)
  [ ] 논문 초안 (Genome Biology)
```

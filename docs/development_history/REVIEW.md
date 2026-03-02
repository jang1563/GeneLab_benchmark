# GeneLab_benchmark — 비판적 종합 리뷰
**기준**: 고품질 학술 논문 출판 + 분야 주요 벤치마크 확립
**검토 버전**: PLAN.md v0.4 (2026-02-27)
**리뷰 날짜**: 2026-02-27

---

## 요약 판정

| 평가 축 | 현재 상태 | 목표 대비 |
|---|---|---|
| 과학적 질문의 명확성 | 보통 | 핵심 내러티브 부재 |
| Task 설계 엄밀성 | 보통~낮음 | 여러 방법론 오류 잔존 |
| 논문 출판 가능성 (v1.0) | 보통 | 구조 재정비 필요 |
| 분야 표준 벤치마크 가능성 | 낮음 | 인프라 근본 요소 누락 |
| 재현성·투명성 | 보통 | CI, 유의성 검정 없음 |
| 경쟁 내성 | 보통 | Foundation model 대응 없음 |

**결론**: 현재 계획은 Genome Biology 수준 논문의 골격은 갖추었으나,
분야 표준 벤치마크로 확립되기 위한 인프라·방법론·내러티브 층이 부족하다.
아래 문제를 해결하지 않으면 출판 후에도 영향력이 제한될 것이다.

---

## Part 1. 논문 출판 기준 리뷰

### 1-1. 과학적 내러티브의 부재 [심각]

**문제**: v1.0의 핵심 질문은 "우주비행 전사체 서명이 미션 간·조직 간 일반화되는가?"이다.
이 질문은 **벤치마크 구성 이유**는 되지만, **논문의 과학적 기여**가 되지 않는다.

Nature Methods, Genome Biology 등 상위 저널에서 벤치마크 논문이 수락되려면
인프라 기술 이상의 **분야를 이해시키는 발견**이 있어야 한다.

성공한 벤치마크 논문의 구조:
```
GLUE (Wang et al., 2018):
  인프라: NLP 다중 태스크 평가 프레임워크
  발견: "현재 모델은 언어 이해의 특정 측면을 체계적으로 실패한다"
  → 이 발견이 GPT-2, BERT 등의 개선 방향을 제시

Quartet (Nat Biotechnology, 2023):
  인프라: 다중오믹스 참조 자료
  발견: "플랫폼별 배치효과가 생물학적 신호를 30~70% 가린다"
  → 구체적 수치로 문제의 크기를 보여줌
```

현재 계획에서 **예상 발견 스토리가 없다**.

**필요한 것**: 데이터 탐색 전에 가설 기반 예상 결과를 먼저 정립하라.

예시 내러티브 (제안):
> "우리는 간이 미션 간 가장 일관된 spaceflight 서명을 보이지만,
> 근육과 신장은 미션 특이적 반응이 지배적임을 보인다.
> ML 모델의 cross-mission 전이 실패는 배치효과보다 생물학적 다양성에서
> 기인함을 보인다 (ComBat 보정 후에도 성능 개선 < 5%).
> 이는 우주비행 생물학 연구에서 단일 미션 결과의 일반화 한계를 정량화한다."

이 수준의 예상 발견이 있어야 논문이 완성된다.

---

### 1-2. LFC Feature의 방법론적 순환성 문제 [심각]

**문제**: Category A (Spaceflight Detection)에서 LFC (log fold-change vs control)를 feature로 사용하면 **label leakage**가 발생한다.

```
LFC 계산:
  LFC_gene_i = log2(mean_flight) - log2(mean_ground)

Category A task:
  입력: LFC_vector (유전자별 LFC)
  출력: "이 샘플이 flight인가?"

문제:
  LFC는 그룹 평균으로 계산 → 어떤 샘플이 flight인지 이미 알아야 계산 가능
  → 분류 목표(flight/ground label)를 미리 사용해서 feature를 만드는 것
```

**이것은 Category A의 근본적 설계 오류다.**

LFC는 **mission-level aggregate feature**이지, sample-level feature가 아니다.
Category A (개별 샘플 분류)에서는 반드시 **sample-level raw counts** 또는
**normalized counts (sample별 size-factor 보정만 적용)**를 사용해야 한다.

**수정 방향**:
```
Category A: 입력 = log2(DESeq2 size-factor normalized counts) per sample
Category B, C: 입력 = mission-level LFC (per-mission 집계, 학습 미션에서만 계산)
```

---

### 1-3. GLDS-168 데이터 중복 오염 위험 [심각]

**문제**: Category A1 (Liver)에 다음을 동시에 포함할 경우:
- GLDS-48: RR-1 간
- GLDS-137: RR-3 간
- GLDS-168: RR-1 + RR-3 간 **통합 분석**

GLDS-168은 GLDS-48과 GLDS-137의 샘플을 **재처리한 것**이다.
즉 동일 샘플이 데이터셋에 두 번 들어갈 수 있다.

LOMO에서 "RR-1 holdout"으로 설정해도,
GLDS-168 내의 RR-1 샘플이 학습 세트에 포함될 수 있다.

**수정 방향**: GLDS-168은 Category A1에서 제외.
대신 GLDS-168을 Category J (파이프라인 비교: 동일 샘플의 단독 처리 vs 통합 처리 결과 비교)에 사용.

---

### 1-4. LOMO Split의 구현 모호성 [중요]

**문제**: LOMO에서 훈련 세트는 여러 미션의 샘플을 합친다.
이때 다음이 명시되지 않았다:

(a) **미션 간 샘플 수 불균형 처리**: RR-1 (n=20) vs STS-135 (n=6). 합칠 때 가중치?

(b) **훈련 세트 내 정규화**: 다른 미션에서 온 샘플들의 normalized counts를 그냥 합쳐도 되는가?
    → 안 된다. 미션별 size factor는 다르다.

(c) **Cross-validation 반복**: LOMO는 fold 수 = 미션 수. 6미션이면 6-fold.
    각 fold에서 독립적으로 hyperparameter 튜닝하는가?

(d) **훈련 시 미션 레이블 사용 금지**: 훈련 세트에서 어떤 샘플이 어느 미션인지 알면 안 됨 (Category A에서는).
    실제 구현에서 이를 보장하는 방법이 없음.

**수정 방향**: LOMO 구현 의사코드를 PLAN에 명시:
```python
for test_mission in all_missions:
    train_missions = [m for m in all_missions if m != test_mission]

    # 각 미션 내에서 독립적으로 size-factor 정규화 (이미 완료된 경우)
    # 미션 간 추가 정규화: quantile normalization 또는 ComBat
    train_X = pool_and_renormalize(train_missions)
    train_y = flight_labels(train_missions)

    # 훈련 (미션 ID feature 사용 금지)
    model.fit(train_X, train_y)

    # 테스트 미션은 독립적으로 정규화
    test_X = normalize_independently(test_mission)
    test_y = flight_labels(test_mission)
    auroc = evaluate(model, test_X, test_y)
```

---

### 1-5. 불확실성 정량화(Uncertainty Quantification) 완전 누락 [중요]

**문제**: 현재 계획은 모든 metric을 point estimate로만 보고한다.
이 경우:
- "Model A AUROC = 0.73, Model B AUROC = 0.71" → 통계적으로 유의미한 차이인가?
- 샘플 수가 적어 AUROC 추정치 자체의 신뢰구간이 넓을 수 있음

Genome Biology, Nature Methods 리뷰어는 95% 신뢰구간 없는 결과를 받아들이지 않는다.

**필요한 것**:
```
1. AUROC 신뢰구간: Bootstrap (n=1000), 95% CI
2. 모델 간 비교: DeLong's test for AUROC, McNemar's test for binary predictions
3. LOMO fold 간 분산: mean ± SD across folds
4. 다중 task 비교 시: FDR correction (Benjamini-Hochberg)
```

---

### 1-6. 수학적 Task 명세의 부재 [중요]

**문제**: 상위 저널 벤치마크 논문은 task를 수학적으로 정의한다.
현재 계획은 자연어 설명만 있다.

**필요한 형식** (Category A1 예시):
```
Task A1: Liver Spaceflight Detection

Given:
  X ∈ ℝ^(n × p) : normalized gene expression matrix
    n = number of samples, p = number of genes after filtering
  y ∈ {0,1}^n : binary flight label (1 = spaceflight, 0 = ground control)
  M ∈ {m₁,...,mₖ} : mission assignment per sample

Predict:
  ŷ = f(X) ∈ [0,1]^n : predicted probability of spaceflight

Evaluation:
  AUROC(y, ŷ) averaged over LOMO folds:
  Score(A1) = (1/K) Σₖ AUROC(y_test_k, f_k(X_test_k))
  where f_k is trained on {m : m ≠ mₖ}

Constraints:
  - f_k must not use mission labels as input features
  - X must be normalized independently per test mission
  - Ground control type (GC/VC) must be recorded and reported separately
```

---

### 1-7. 음성 대조군(Negative Control) 부재 [중요]

**문제**: 벤치마크 검증에 필수적인 negative control이 없다.

Negative control의 역할:
- "이 AUROC가 실제로 의미 있는 수치인가?" 확인
- 노이즈 수준 측정

**필요한 negative controls**:
```
NC1: 무작위 레이블 섞기 (permutation test)
  - y를 무작위로 섞은 후 동일 pipeline 실행
  - 기대값: AUROC ≈ 0.50 ± small CI
  - 실제값이 이보다 유의미하게 높아야 task에 신호가 있는 것

NC2: 무관한 유전자 집합 사용
  - 우주비행과 무관한 하우스키핑 유전자만으로 분류 시도
  - 기대값: AUROC ≈ 0.50
  - 실제 feature selection의 효과 검증

NC3: 잘못된 조직 쌍 (Category C용)
  - 간으로 학습, 완전히 다른 종(예: Arabidopsis)으로 테스트
  - 기대값: Transfer AUROC ≈ 0.50
```

---

## Part 2. 분야 표준 벤치마크 기준 리뷰

### 2-1. 공개 Leaderboard 및 제출 시스템 부재 [심각]

**문제**: 현재 계획에는 벤치마크 제출 시스템이 없다.

GLUE, SuperGLUE, BEIR, CASP 등 분야를 대표하는 벤치마크의 공통 요소:
```
1. 공개 리더보드 (EvalAI, CodaLab, HuggingFace Spaces)
2. 표준화된 제출 형식 (JSON / parquet)
3. 자동 평가 시스템
4. 제출별 성능 기록 및 공개
```

Leaderboard 없이는:
- 다른 그룹들이 벤치마크에 참여할 구조가 없다
- "benchmark"라고 부를 수 없다 ("dataset" 또는 "evaluation framework"에 불과)
- 논문 인용 후 실제 사용자가 생기지 않는다

**최소 요건**: HuggingFace Datasets로 처리된 데이터 + HuggingFace Spaces 리더보드.
이것이 현재 가장 낮은 비용으로 구축 가능한 방법이다.

---

### 2-2. Held-out Test Set 없음 [심각]

**문제**: 현재 모든 split 정보가 공개된다.
연구자들이 테스트 셋 레이블을 보고 모델을 튜닝할 수 있다.

분야 표준 벤치마크는 두 단계로 나눈다:
```
개발 셋 (public): 모델 개발·튜닝용. 정답 레이블 공개.
테스트 셋 (hidden): 최종 평가용. 정답 레이블 비공개.
           → 제출 시스템에 예측값만 제출 → 서버에서 점수 계산
```

우주비행 데이터 특성상 구현 방법:
```
옵션 1: 가장 최근 미션 1개를 완전 hold-out
  예: RRRM-1 간 데이터 → 어떤 분석에도 사용 안 함
  → 최종 리더보드 점수만 이 미션으로 계산

옵션 2: 각 task별 20% 샘플을 서버에만 보관
  → 개발 셋 AUROC와 서버 테스트 AUROC 둘 다 리더보드에 표시
```

---

### 2-3. Foundation Model 지원 없음 [중요]

**문제**: 현재 baseline은 LogReg, RF, XGBoost, LightGBM이다.
2025~2026년 기준으로 유전자 발현 분야의 state-of-the-art는 이미 foundation model이다:

- **Geneformer** (Theodoris et al., Nature 2023): 3천만 세포로 pre-trained transformer, rank-based tokenization
- **scGPT** (Cui et al., Nature Methods 2024): single-cell foundation model
- **Universal Cell Embeddings (UCE)**: 33M cell pre-trained
- **AIDO.Cell**: multi-organism cell representation

분야 표준 벤치마크가 되려면 이들을 평가 대상에 포함해야 한다.
두 가지 이유:

1. **현실적 이유**: 논문 출판 시 reviewer가 "현재 SOTA인 foundation model과 비교했는가?" 묻는다.
2. **임팩트 이유**: foundation model이 이 벤치마크에서 실패하면 → 우주비행 데이터의 특수성 증명. 성공하면 → 일반 생물학적 표현이 우주비행에 전이된다는 발견.

**필요한 것**:
- Geneformer 제로샷(zero-shot) 평가를 baseline에 추가
- Fine-tuned Geneformer와 비교
- "spaceflight pre-training이 필요한가?"를 category B로 검증

---

### 2-4. 버전 관리 및 Living Benchmark 메커니즘 없음 [중요]

**문제**: NASA OSDR에는 새로운 미션 데이터가 지속적으로 추가된다.
현재 계획에는 새 데이터 추가 시 처리 방법이 없다.

성공한 living benchmark 사례:
```
CASP (단백질 구조 예측):
  - 2년마다 새 단백질 구조로 blind challenge
  - 과거 버전과 현재 버전 결과 분리

HELM (언어 모델 벤치마크):
  - 버전 번호 엄격 관리
  - 이전 버전 결과는 이전 데이터에서만 유효
```

**필요한 메커니즘**:
```
v1.0 (2026): GLDS 목록 고정, 영구적으로 동결
  → DOI 부여, Zenodo에 데이터 스냅샷 저장
v1.1 (2027): 새 미션 추가 가능하지만 v1.0 task 결과와 비교 불가
  → 새 task만 추가, 이전 task 점수 유지
v2.0: 새 카테고리 추가 (E+F+G)
  → 별도 논문, 별도 리더보드
```

이 없이는 1년 후 데이터가 바뀌면 리더보드가 무의미해진다.

---

### 2-5. pip 설치 가능한 패키지 없음 [중요]

**문제**: 분야 표준이 되려면 사용하기 쉬워야 한다.
현재 계획은 스크립트 모음이다.

BEIR, MTEB 등 성공한 벤치마크는 모두 pip 패키지:
```bash
pip install beir
pip install mteb
```

**목표 사용 경험 (이상적)**:
```python
from genelabench import load_task, evaluate

# Task 로드
task = load_task("A1_liver_detection")
train_data, test_data = task.get_splits()

# 모델 실행
predictions = my_model.predict(test_data.X)

# 평가
results = evaluate(task, predictions)
# → {"auroc": 0.81, "auroc_ci": [0.74, 0.87], "auprc": 0.79}

# 리더보드 제출
task.submit(predictions, model_name="my_model")
```

이 수준의 API가 없으면 다른 그룹이 사용하지 않는다.

---

## Part 3. 방법론 심층 리뷰

### 3-1. Cross-Tissue Transfer (Category C) 설계 개선 필요

**현재 설계의 문제**:
"Liver feature로 학습 → Kidney에 적용"이 두 가지 다른 방법으로 해석된다.

**방법 A (Covariate shift)**: Liver와 Kidney가 같은 유전자 공간을 공유.
Liver에서 학습한 분류기의 가중치(weight)를 Kidney에 그대로 적용.
```
문제: 간세포 특이적 마커(ALB, CYP3A4)가 신장 샘플에서 발현이 0이면
     모델이 이 유전자에 의존하면 → 완전히 실패
```

**방법 B (Feature transfer)**: Liver에서 중요한 유전자 집합(DEG list)을 추출.
해당 유전자만 사용해서 Kidney에서 새 모델 학습 (또는 평가).
```
더 합리적이지만:
  - "새 모델 학습"이면 transfer가 아니라 feature selection
  - task 정의가 애매해짐
```

**방법 C (Pathway-level)**: 유전자 수준이 아닌 pathway 점수(GSEA NES) 수준에서 전이.
Liver pathway 점수 → 같은 pathway를 Kidney에서 계산 → 분류.
```
가장 생물학적으로 의미 있음
공통 pathway (미토콘드리아, Nrf2 등)는 조직 불문 존재
```

**권장**: 방법 C를 Category C의 주요 방법으로 채택. 방법 A는 supplementary.

---

### 3-2. Category D 회귀 task (D1)의 구조적 한계

**문제**: D1 (비행 기간 회귀)의 실제 데이터 포인트:

```
STS-135:  13일
BION-M1:  30일
RR-6:     30~55일 (가변)
RR-8:     39일
RR-9:     33일
RR-1:     37일
RR-3:     39일
RR-7:     75일 (일부)
```

실질적으로 세 클러스터: 단기(13d), 중기(30~40d), 장기(75d).
이것은 연속 회귀가 아닌 **3-class 분류**와 사실상 동일하다.

Pearson r로 측정하면 "장기 미션과 단기 미션을 구분"하는 것만으로도 높은 r을 얻을 수 있어, 실제 비행 기간 예측 능력을 과대평가한다.

**수정 방향**:
- D1을 "비행 기간 예측"이 아닌 "단기/중기/장기 비행 분류"로 재정의 (D2와 통합)
- 또는 D1 유지하되: 기간 포인트 클러스터링 한계를 Methods에 명시하고, within-cluster 예측 실패를 분석

---

### 3-3. 다중 검정(Multiple Testing) 문제

**문제**: 25개 task × 5개 baseline × 2개 primary metric = 250개 비교.
이 중 우연에 의해 유의한 것처럼 보이는 결과가 상당수 나올 수 있다.

현재 계획에 다중 검정 보정이 없다.

**필요한 것**:
```
1. Task 내 모델 비교: DeLong's test (AUROC) or paired bootstrap
2. Task 간 종합 순위: Composite score의 bootstrap CI
3. 다중 task에 걸친 주장: Bonferroni or FDR correction
4. 보고 기준: p < 0.05 (bonferroni-corrected) 또는 95% CI 비중첩
```

---

### 3-4. 생물학적 Ground Truth 검증 설계 미흡

**현재 계획**: 특정 유전자(ANGPTL4, MuRF1 등)의 feature importance 상위 여부를 "sanity check"로 사용.

**문제**: Feature importance는 모델 의존적이고 불안정하다.
- Random Forest feature importance: split 기반
- Linear coefficient: L2 정규화에 민감
- SHAP values: 더 안정적이지만 computationally expensive

**개선된 설계**:
```
생물학적 검증 방법 1: Gene Set Enrichment
  top-100 feature genes로 GSEA 실행 → spaceflight pathway 농축 확인
  (FOXO, Nrf2, mitochondria, oxidative stress gene sets 사용)

생물학적 검증 방법 2: 외부 검증 셋
  Cell 2020 (마우스 13조직 spaceflight DEG) → 이 벤치마크의 top feature와 overlap 계산
  SOMA 2024 (인간 spaceflight DEG) → cross-species 보존도 확인

생물학적 검증 방법 3: Direction concordance
  각 DEG의 up/down 방향이 알려진 spaceflight 반응과 일치하는가?
  (예: 미세중력에서 근위축 유전자는 반드시 상향)
```

---

## Part 4. 경쟁 환경 분석 및 포지셔닝

### 4-1. SOMA 팀 선점 위험 [중요]

**현황**: Mason Lab (SOMA, Nature 2024)은 인간 멀티미션 데이터의 소유자이며
지속적으로 후속 분석을 출판 중.

**위험**: 이 벤치마크 v1.0 출판 전에 Mason Lab이 유사한 ML 평가 논문을 출판할 수 있다.

**차별화 핵심**: 이 벤치마크의 차별점은 인간이 아닌 **마우스 다조직 크로스미션 일반화**다.
SOMA는 이 방향으로 움직이지 않았다.

**권장 행동**: Phase 1 완료 시점에 preprint (bioRxiv) 즉시 업로드 → 선점권 확보.
저널 리뷰 대기 중에도 커뮤니티에 공개하여 사용자 확보.

---

### 4-2. GeneLab 공식 분석 그룹과의 관계

**문제**: AWG AI/ML 그룹이 내부적으로 유사한 작업을 진행 중일 수 있다.
공동저자로 포함시키거나 참조 논문으로 인용하는 것이 "태클" 방어에 유리하다.

**AWG 협력 활용 방법**:
- 내부 미션 메타데이터 직접 접근 (RadLab, EDA 정확한 미션 매핑)
- GeneLab processing pipeline 버전 이력 정확한 정보
- 미공개 QC 리포트 접근
- 논문 endorsement ("NASA GeneLab AI/ML AWG의 공식 협력 연구"로 표기)

---

### 4-3. ESA / JAXA 데이터 미포함 [보통]

**현황**: MHU (JAXA) 데이터는 포함되어 있으나, ESA 우주생물학 데이터는 없다.

ESA Open Science: https://esamultimedia.esa.int/docs/SME/ESA_Open_Data_Policy.pdf
ESA도 OSDR와 유사한 데이터 공유 정책을 가지고 있으며, 일부 데이터가 OSDR에 공동 기탁되어 있다.

국제 표준 벤치마크를 목표로 한다면 ESA 데이터 포함 검토 필요.
v2.0 이상에서 포함하더라도 PLAN에 언급하는 것이 전략적으로 유리하다.

---

## Part 5. 현재 계획의 강점 (유지해야 할 것)

다음 설계 결정은 이 벤치마크를 기존 작업보다 우수하게 만든다.
변경하지 말 것.

### 강점 1: Control Type 두 트랙 (GC vs VC) [핵심 기여]
이것은 이 분야에서 아무도 체계적으로 다루지 않았다.
Track 1 (GC) vs Track 2 (VC)의 성능 차이가 곧 "하드웨어 효과의 크기"가 된다.
논문의 Figure 1급 결과가 될 수 있다.

### 강점 2: LOMO Full Cross-Table (B1)
6미션 × 6미션 transfer AUROC 행렬은 시각적으로 강력하다.
"어느 미션 쌍이 서로 일반화되는가?"가 한눈에 보인다.
마우스 계통, 미션 기간, 방사선량 등과 cross-table을 연결하면 생물학적 발견이 된다.

### 강점 3: Confounder 정량화 (D4, D5)
D4 (계통 분류), D5 (하드웨어 분류)는 "우리는 confounder를 숨기지 않는다"는 메시지다.
이것이 없으면 reviewer가 "계통 차이가 우주비행 차이보다 크지 않냐"고 공격할 것이다.
D4, D5로 confounder 크기를 정량화하면 이 공격을 선제적으로 막는다.

### 강점 4: HU Analog 포함 (B6)
ISS flight vs Hindlimb Unloading은 "미세중력 특이 반응 vs 우주비행 전체 반응"을 분리한다.
이것이 성공하면 "진짜 spaceflight unique signature"를 정의하는 데 기여한다.

### 강점 5: AWG 내부 협력
공식 NASA 협력은 데이터 접근, 리뷰어 방어, 출판 홍보 모두에서 유리하다.

---

## Part 6. 누락된 구성요소 목록

### 논문 구조상 누락
- [ ] Motivation figure: "왜 cross-mission generalization이 중요한가" 실증 데이터
- [ ] Expected findings narrative: 출판 전 가설 기반 예상 결과
- [ ] Limitation section: 소샘플, 마우스 중심, 인간 관련성의 한계 명시
- [ ] Data statement: OSDR 데이터 접근 정책, 재배포 제한

### 방법론 누락
- [ ] 수학적 task 명세 (formal notation)
- [ ] Bootstrap CI for all metrics
- [ ] Permutation-based p-values (NC1)
- [ ] Multiple testing correction 정책
- [ ] LOMO 구현 의사코드
- [ ] Foundation model (Geneformer) zero-shot baseline

### 인프라 누락
- [ ] Held-out test set 설계
- [ ] 제출 시스템 (최소: HuggingFace Spaces)
- [ ] pip 패키지 (`genelabench`)
- [ ] Benchmark versioning 정책 (v1.0 영구 동결 + Zenodo DOI)
- [ ] Living benchmark 업데이트 절차

### 데이터 누락
- [ ] GLDS-168 중복 처리 정책 확정
- [ ] LFC vs normalized counts 용도별 명확 분리
- [ ] Negative control tasks (permutation, housekeeping genes)
- [ ] HU analog GLDS 목록 (catalog_datasets.py 결과 후)

---

## Part 7. 수정 우선순위

### 즉시 수정 (코드 착수 전 필수)

| 우선순위 | 항목 | 이유 |
|---|---|---|
| P0 | Category A feature: LFC → normalized counts | 방법론 오류 |
| P0 | GLDS-168 중복 처리 정책 | 데이터 오염 |
| P0 | LOMO 구현 pseudocode 확정 | 미션 간 정규화 방법 결정 |
| P0 | 수학적 task 명세 작성 | 구현 기준 |
| P1 | Bootstrap CI 구현 계획 | 출판 필수 |
| P1 | Permutation test (NC1) 추가 | 신호 검증 |
| P1 | 예상 발견 스토리 확정 | 논문 내러티브 |

### 단기 수정 (Phase 1 중)

| 우선순위 | 항목 |
|---|---|
| P2 | DeLong's test 구현 |
| P2 | Geneformer zero-shot baseline 추가 |
| P2 | Category C: Pathway-level 방법 설계 |
| P2 | D1 (비행 기간 회귀) 재정의 |
| P3 | HuggingFace 리더보드 프로토타입 |
| P3 | bioRxiv 업로드 일정 확정 |

### 중기 수정 (Phase 2 중)

| 우선순위 | 항목 |
|---|---|
| P3 | Benchmark 버전 동결 및 Zenodo DOI |
| P3 | pip 패키지 기초 구조 |
| P4 | ESA 데이터 포함 여부 결정 |
| P4 | Living benchmark 절차 문서화 |

---

## 결론

현재 GeneLab_benchmark 계획은 **아이디어와 방향은 옳지만**,
고품질 논문과 분야 표준 벤치마크가 되기 위한 핵심 구성요소가 여러 층에서 누락되어 있다.

**가장 치명적인 3가지**:
1. Category A의 LFC feature leakage → 즉시 수정 필수
2. Held-out test set + 제출 시스템 없음 → "벤치마크"가 아닌 "평가 프레임워크"에 그침
3. 예상 발견 내러티브 없음 → 상위 저널 수락 불가

**가장 강력한 3가지** (바꾸지 말 것):
1. Control type 두 트랙 (GC vs VC) — 새로운 생물학적 기여
2. LOMO full cross-table — 시각적으로 강력한 핵심 결과
3. Confounder 정량화 (D4, D5) — 선제적 비판 방어

이 리뷰에서 지적한 P0 항목들을 수정한 후
실제 구현을 시작하는 것을 권장한다.

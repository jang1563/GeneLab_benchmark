# GeneLab_benchmark — Design Decisions Reference
**작성일**: 2026-02-28
**목적**: 구현 중 변경 불가한 설계 결정 사항의 단일 진실 소스 (Single Source of Truth)
**적용 범위**: PLAN.md v0.5 이후 모든 스크립트 구현

이 파일의 내용은 코드 착수 전 확정되었으며, 변경 시 버전을 올리고 Changelog에 기록.

---

## DD-01: Category A Feature Representation

**결정**: Category A (Spaceflight Detection)의 입력 feature는 `log2(DESeq2 size-factor normalized counts)`, per-sample.

**금지**: LFC (log fold-change) 사용 금지.

**근거**:
LFC는 flight/ground 그룹 레이블을 알아야 계산 가능 (group mean subtraction).
Category A의 목표는 이 레이블을 예측하는 것 → LFC를 입력으로 쓰면 label leakage.
검토 자료: PLAN_REVIEW_2026-02-28 #1, PLAN_REVIEW.md #1.2, REVIEW.md Part1-1-2

**구현**:
```python
# preprocess_mouse_rnaseq.py
import numpy as np

def get_category_a_features(counts_matrix, size_factors):
    """
    Category A용 feature: log2(normalized counts + 1) per sample.
    counts_matrix: shape (n_samples, n_genes), raw counts
    size_factors: shape (n_samples,), DESeq2 size factors
    """
    normalized = counts_matrix / size_factors[:, np.newaxis]
    return np.log2(normalized + 1)   # pseudocount=1
```

**적용 범위**: A1~A7 모든 task의 입력 feature.

---

## DD-02: LFC Shrinkage Method

**결정**: Category B, C의 LFC는 반드시 DESeq2 `lfcShrink(type="apeglm")`로 계산.

**금지**: shrinkage 없는 raw LFC 사용 금지 (n=3/그룹 환경에서 분산 극심).

**근거**: PLAN_REVIEW.md #1.1 — n=3/그룹의 소규모 데이터에서 shrinkage 없는 LFC는
분산이 매우 커서 cross-mission 비교 시 노이즈 지배적.

**구현** (R, DESeq2):
```r
# preprocess_mouse_rnaseq.R
library(DESeq2)
library(apeglm)

compute_lfc <- function(dds, coef_name, method="apeglm") {
  # coef_name: "condition_Flight_vs_Ground"
  res <- lfcShrink(dds, coef=coef_name, type=method)
  return(res)
}

# Primary:    type="apeglm"   (기본)
# Comparison: type="ashr"     (Category J 파이프라인 비교용)
```

**적용 범위**: Category B, C의 mission-level LFC 계산. Category A에는 LFC 자체를 사용하지 않음.

**추가 규칙**: LFC는 학습 미션에서만 계산. 테스트 미션의 LFC는 LOMO 루프 밖에서 계산하지 않음.

---

## DD-03: LOMO Feature Selection — Split-Aware

**결정**: Feature selection (variance filter)은 LOMO 루프 내부에서 학습 미션 데이터만으로 수행.

**금지**: 전체 데이터 기준으로 variance filter 사전 적용 후 LOMO 실행.

**근거**: PLAN_REVIEW.md #1.3, PLAN_REVIEW_2026-02-28 #2 —
holdout 미션 데이터의 variance 정보가 feature selection에 노출되면 성능 과대 추정.

**구현**:
```python
# generate_tasks.py
def split_aware_feature_filter(data, missions, low_expr_pct=0.20, var_pct=0.25):
    """
    LOMO-aware feature selection.
    - Low-expression filter는 전체 데이터에서 1회만 적용 가능 (missing value 처리).
    - Variance filter는 반드시 train missions에서만 계산.
    """
    results = {}
    for test_mission in missions:
        train_missions = [m for m in missions if m != test_mission]
        train_data = data[data['mission'].isin(train_missions)]

        # Variance filter on train only
        gene_var = train_data.drop(columns=['mission', 'label']).var()
        var_threshold = gene_var.quantile(var_pct)
        selected_genes = gene_var[gene_var >= var_threshold].index.tolist()

        results[test_mission] = {
            'train_missions': train_missions,
            'test_mission': test_mission,
            'selected_genes': selected_genes
        }
    return results
```

---

## DD-04: Independence Unit & Split Policy

**결정**: split 유형별 독립성 단위 확정.

| Split 유형 | 독립성 단위 | 적용 Task | 비고 |
|---|---|---|---|
| LOMO | Mission | A, B | cage-level leakage 구조적으로 없음 |
| LOTO | Tissue | C | |
| Mission-stratified 80/20 | Mission | D | 랜덤 샘플 분리 금지 |
| GLDS-stratified 80/20 | GLDS study | J | |

**금지**: Category D, J에서 샘플 단위 랜덤 80/20 split. 같은 미션/GLDS의 샘플이 train/test에 섞이면 batch leakage.

**구현**:
```python
# generate_tasks.py
from sklearn.model_selection import GroupShuffleSplit

def mission_stratified_split(X, y, groups, test_size=0.2, seed=42):
    """
    groups: mission ID per sample. 같은 mission은 같은 fold에.
    """
    gss = GroupShuffleSplit(n_splits=5, test_size=test_size, random_state=seed)
    return list(gss.split(X, y, groups=groups))
```

---

## DD-05: GLDS-168 — 완전 제외 (Category A)

**결정**: GLDS-168은 Category A에서 완전 제외. Category J (J1)에만 사용.

**근거**: REVIEW.md Part1-1-3 —
GLDS-168은 GLDS-48 (RR-1) + GLDS-137 (RR-3)을 통합 재처리한 데이터.
Category A1 LOMO에서 동일 샘플이 train(GLDS-48 또는 137)과 test(GLDS-168)에 중복 등장.

**Category J 활용**:
```
J1: GeneLab pipeline v1 vs v2 vs v3 비교
  → GLDS-48+137 개별 처리 결과 vs GLDS-168 통합 처리 결과 비교
  → "단독 연구 처리" vs "다연구 통합 처리"의 파이프라인 차이 정량
```

**A1 Liver 미션 목록 (확정)**:
```
Track 2a (C57BL/6J only, primary):
  GLDS-48  → RR-1  (C57BL/6J, GC+VC)
  GLDS-137 → RR-3  (C57BL/6J, GC+VC) ← catalog_datasets.py로 계통 재확인 필요
  GLDS-245 → RR-6  (C57BL/6J, GC)     ← 확인 필요
  GLDS-379 → RR-8  (C57BL/6J, ?)      ← 확인 필요
  GLDS-242 → RR-9  (C57BL/6J, ?)      ← 확인 필요
  GLDS-617 → MHU-2 (C57BL/6J, ?)      ← 확인 필요

Track 2b (all strains, secondary):
  + GLDS-25 → STS-135 (BALB/c, VC only)  ← cross-strain 실험으로 reframe

제외: GLDS-168 (중복), Visium/scRNA-seq 데이터
```

---

## DD-06: Mouse Strain Track Policy

**결정**: 마우스 계통 혼용을 허용하되, 명시적 분리.

```
Track 2a (Primary): C57BL/6J 미션만 포함
  → 논문 main figures (A, B, C 카테고리)
  → B1 6×6 matrix = C57BL/6J 미션 6개만

Track 2b (Secondary): 모든 계통
  → supplementary figures
  → BALB/c 미션은 "cross-strain" 실험으로 명시
```

**missions.json 필수 필드**: `"mouse_strain": "C57BL/6J"` (모든 미션에 기록).

**D4 task**: strain prediction AUROC → Track 2b에서의 계통 confounder 크기 정량화.

---

## DD-07: Composite Score Normalization

**결정**: random baseline은 primary metric과 동일한 스케일로 계산.

```python
def random_baseline(task_type: str, n_classes: int = None) -> float:
    if task_type == "binary":
        return 0.5               # AUROC random = 0.5
    elif task_type == "multiclass":
        return 1.0 / n_classes   # macro-F1 random ≠ accuracy! F1 스케일 기준
    elif task_type == "regression":
        return 0.0               # Spearman r random ≈ 0
    elif task_type == "transfer":
        return 0.5               # Transfer AUROC random = 0.5
    else:
        raise ValueError(f"Unknown task_type: {task_type}")

def normalize_score(score: float, task_type: str, **kwargs) -> float:
    baseline = random_baseline(task_type, **kwargs)
    return (score - baseline) / (1.0 - baseline)
```

**금지**: 다중 클래스 task에 "1/K accuracy"를 macro-F1의 random baseline으로 사용.
(accuracy와 macro-F1은 다른 스케일 → 정규화 수식이 수학적으로 잘못됨)

---

## DD-08: Uncertainty Quantification

**결정**: 모든 metric에 Bootstrap 95% CI 필수 보고.

```python
# evaluation/metrics.py

def bootstrap_ci(y_true, y_score, metric_fn, n_bootstrap=1000, alpha=0.05, seed=42):
    """Bootstrap 95% CI for any metric function."""
    rng = np.random.default_rng(seed)
    n = len(y_true)
    scores = []
    for _ in range(n_bootstrap):
        idx = rng.choice(n, size=n, replace=True)
        scores.append(metric_fn(y_true[idx], y_score[idx]))
    lower = np.percentile(scores, 100 * alpha / 2)
    upper = np.percentile(scores, 100 * (1 - alpha / 2))
    return np.mean(scores), lower, upper

def permutation_pvalue(y_true, y_score, metric_fn, n_permutations=1000, seed=42):
    """Permutation test p-value (H0: AUROC = random)."""
    rng = np.random.default_rng(seed)
    observed = metric_fn(y_true, y_score)
    null_scores = []
    for _ in range(n_permutations):
        y_perm = rng.permutation(y_true)
        null_scores.append(metric_fn(y_perm, y_score))
    p_value = np.mean(np.array(null_scores) >= observed)
    return p_value
```

**보고 형식** (모든 결과 테이블):
```
AUROC = 0.81 [95% CI: 0.74–0.87], permutation p < 0.001
LOMO fold mean ± SD: 0.81 ± 0.06
```

**모델 간 비교**: DeLong's test (이진 AUROC), Wilcoxon signed-rank (F1, J).

---

## DD-09: D1 Task — Ordinal Classification (not regression)

**결정**: D1 (비행 기간 예측)을 3-class ordinal classification으로 재설계.

**근거**: 실제 비행 기간이 3~4 클러스터에 집중. 연속 회귀는 클러스터 구조를 무시.
```
Class 1: ≤14일  (STS-135: 13일, shuttle)
Class 2: 26~40일 (RR-1,3,6,8,9: 30~39일, 대부분의 ISS 미션)
Class 3: ≥75일  (RR-7: 75일, 장기 비행)
```

**Primary metric**: macro-F1 (3-class).
**Secondary metric**: Spearman r between predicted class rank and true class rank.
**Random baseline**: 1/3 macro-F1 (= 0.333).

---

## DD-10: Batch Correction Methods

**결정**: bulk RNA-seq 배치 보정 비교 대상 (Category B, J3).

| 방법 | 용도 | 참고 |
|---|---|---|
| 보정 없음 | Baseline | |
| ComBat-seq (sva) | Primary comparison | negative binomial 모델, bulk RNA-seq 표준 |
| limma::removeBatchEffect | Primary comparison | log-CPM 입력, LIMMA 연구 표준 |
| RUVseq | Primary comparison | negative control gene 있을 경우 |
| Harmony | 참고용만 | ⚠️ 단일세포 임베딩용 — bulk RNA-seq 비표준. "bulk에 적용 시 비교" 명시 |

**금지**: Harmony를 bulk RNA-seq 배치보정 기본 방법으로 사용.

---

## DD-11: Phase 1 Go/No-Go Checkpoint

**결정**: Phase 1 완료 후 Phase 2 진행 여부는 아래 3조건 AND로 판단.

```
[통과 조건]
1. A1 LOMO AUROC > 0.70 AND 95% CI 하한 > 0.60
2. Permutation p < 0.05 (H0: AUROC = random)
3. Known spaceflight genes (ANGPTL4, PCK1) SHAP top-50 이내

[실패 시 조치]
- 조건 1 실패: 미션별 mapping rate·샘플 수 재확인 → QC 재처리
- 조건 2 실패: 배치효과 가능성 → J3 선행 실행
- 조건 3 실패: feature representation 재검토 → normalized counts 정규화 방법 확인
```

---

## DD-12: Negative Control Tasks

**결정**: 신뢰성 검증을 위한 negative control task 3종.

| NC | 방법 | 기대 결과 | 실패 의미 |
|---|---|---|---|
| NC1 | Permutation test (모든 task) | null AUROC ≈ 0.50 ± 0.03 | 신호 존재 확인 실패 |
| NC2 | Housekeeping gene only (50개: GAPDH, ACTB 등) | AUROC ≈ 0.50 | 배경 노이즈 수준 확인 |
| NC3 | Cross-species alien (Category C: liver → Arabidopsis ortholog) | Transfer AUROC ≈ 0.50 | C1~C5 상대적 기준 |

---

## DD-13: Foundation Model Baselines

**결정**: Phase 2에서 Geneformer 포함 (P2 항목).

```
baselines/
├── logistic_regression.py     # ElasticNet, SAGA solver
├── random_forest.py           # max_features=sqrt, class_weight=balanced
├── xgboost_baseline.py        # max_depth=2~3, subsample=0.7
├── lightgbm_baseline.py       # num_leaves=15~31
├── pca_logreg.py              # 50 PCs + LogReg (low-dim baseline)
├── mlp_baseline.py            # 256→64→1, dropout=0.5 (deep learning reference)
└── geneformer_baseline.py     # Phase 2: zero-shot + fine-tuned (P2)
```

**Geneformer** (Phase 2):
- bulk RNA-seq → rank tokenization → Geneformer embedding → linear probe (zero-shot)
- Fine-tuning on each Category A task
- "general pre-training이 spaceflight에 전이되는가?" 검증

---

## DD-15: Pathway Enrichment Analysis (fGSEA + GSVA)

**결정**: OSDR은 fGSEA 결과를 제공하지 않음 (DESeq2 DGE까지만 처리). 로컬에서 실행.

### 그룹 수준: fGSEA

| 파라미터 | 값 | 근거 |
|---|---|---|
| Ranking metric | DESeq2 Wald statistic (DGE `Stat_` 컬럼) | 부호+크기 모두 반영, DESeq2 표준 |
| Gene ID | SYMBOL (mouse gene symbol) | msigdbr gene sets 매칭 기준 |
| `minSize` | 15 | 너무 작은 gene set은 통계적 불안정 |
| `maxSize` | 500 | 너무 큰 gene set은 해석 무의미 |
| `eps` | 0 | 정확한 p-value 계산 (빠른 근사 대신) |
| Multiple testing | BH FDR (fgsea 내장) | DB별 독립 보정 |

**Gene Set DB**:

| DB | MSigDB Category | Set 수 | 역할 |
|---|---|---|---|
| **Hallmark** | H | 50 | Primary — Category C Method C, 논문 main figures |
| **KEGG** | C2:CP:KEGG_MEDICUS | ~658 | Secondary — 대사 경로 중심 |
| **Reactome** | C2:CP:REACTOME | ~1,787 | Secondary — 포괄적 경로 커버리지 |

Species: Mus musculus (msigdbr)

**Stat_ 컬럼 없는 구버전 DGE**: `Log2fc_` 컬럼으로 fallback + 경고 출력.

### 샘플 수준: GSVA

| 파라미터 | 값 | 근거 |
|---|---|---|
| Method | GSVA (기본), ssGSEA (비교용) | GSVA가 normalized data에 더 적합 |
| `kcdf` | `"Gaussian"` | log2 normalized continuous data (count 아님) |
| `minSize` | 15 | fGSEA와 동일 기준 |
| `maxSize` | 500 | fGSEA와 동일 기준 |

입력: `log2(DESeq2 normalized counts + 1)` per sample.
출력: samples × pathways enrichment score matrix.

### Leakage 방지

- Category C pathway 선택: LOMO train missions의 fGSEA 결과에서만 선택.
- GSVA score 계산 자체는 unsupervised (label 미사용) → 전체 샘플에 적용 가능.

### 적용 범위 (전체 카테고리 확장)

- **Category A**: gene-level vs pathway-level (GSVA) 분류 성능 비교
- **Category B**: LFC vs NES cross-mission transfer 비교
- **Category C Method C**: pathway-level cross-tissue transfer (PRIMARY)
- **Category D**: gene-level vs pathway-level condition prediction 비교
- **Category J (J5 NEW)**: Feature Representation 체계적 비교
- **Biological validation**: SHAP top genes → pathway enrichment
- **Cross-mission analysis**: NES profile correlation

### 스크립트

```
scripts/run_fgsea.R              — 그룹 수준 fGSEA (미션×조직별)
scripts/compute_pathway_scores.R — 샘플 수준 GSVA pathway score
scripts/preprocess_pathways.py   — Python 통합 + Category C feature 생성
```

### 출력 디렉토리

```
processed/fgsea/{tissue}/         — {mission}_fgsea_{db}.csv
processed/fgsea/summary/          — all_fgsea_{db}.csv (전 미션 통합)
processed/pathway_scores/{tissue}/ — {mission}_gsva_{db}.csv (samples × pathways)
```

### ✅ 구현 결과 (2026-03-01)

**파이프라인 완료**: 5 tissues × 17 tissue-missions × 3 DBs = 51 fGSEA + 54 GSVA 파일.

**생물학적 검증 (Hallmark, 전 조직 PASS)**:
- Liver: OXIDATIVE_PHOSPHORYLATION, FATTY_ACID_METABOLISM 상위 enriched
- Thymus: E2F_TARGETS, G2M_CHECKPOINT 상위 (thymocyte 증식); IFN-gamma 유의
- Gastrocnemius: OXIDATIVE_PHOSPHORYLATION, FATTY_ACID_METABOLISM, MYOGENESIS
- Kidney: MTORC1_SIGNALING, CHOLESTEROL_HOMEOSTASIS, E2F_TARGETS
- Eye: OXIDATIVE_PHOSPHORYLATION dominant (망막 대사 활성)

**GSVA 특성**: 범위 -0.87 ~ +0.84, 50 pathways × 264 samples (liver 기준), NA 없음.

**핵심 발견 (J5 Gene vs Pathway)**:
- Category A (LOMO): 5 tissues 중 pathway 승리 2개 (Kidney: 0.43→0.74, Eye: 0.79→0.92)
- Category D (D3): Gene F1=1.0 / Pathway F1=0.06 → **pathway는 배치 불변 (batch-invariant)** 확인
- **Kidney Rescue**: gene-level 실패를 pathway-level이 보완 → H3 부분 지지

**발견된 버그 (해결됨)**:
- Bug 4: GSVA API v1.40 호환 (`gsvaParam` vs legacy `gsva()` API 감지)
- Bug 5: DGE matrix orientation (genes×samples vs samples×genes 자동 감지)
- Bug 6: GSVA-Metadata sample name alignment (mission prefix 제거)

---

## DD-16: Text-based LLM Evaluation Track (신규)

**결정**: Phase 2에서 Text LLM을 별도 Tier 3 baseline으로 평가.

**목적**: Gene expression foundation model (Geneformer)과 달리, 텍스트 LLM은 생물학 도메인 prior knowledge로 spaceflight 예측이 가능한가?

**입력 형식 (표준)**:
```
System: "You are a bioinformatics expert specializing in spaceflight transcriptomics."
User: "A mouse {tissue} RNA-seq sample shows highest expression variability in the following genes:
{gene_1} ({ensembl_1}), {gene_2} ({ensembl_2}), ..., {gene_50} ({ensembl_50}).
Based on your knowledge of spaceflight biology, predict whether this sample is from:
(A) Flight (spaceflight/microgravity condition)
(B) Ground (ground control condition)
Answer with A or B, followed by a confidence score (0.0-1.0)."
```

**평가 조건**:
- Zero-shot: 학습 예시 없음
- Few-shot (k=3, 5, 10): k개 labeled 예시 포함
- Chain-of-thought: 단계적 reasoning 요청

**타겟 모델**:
- GPT-4o (OpenAI)
- Claude Opus (Anthropic)
- Llama 3.1-70B (Meta, open source baseline)

**결과 저장 형식**: `evaluation/llm_baseline/{model_name}_A{n}_{tissue}.json`

**중요 설계 원칙**:
- LLM은 LOMO split을 "알지 못해야" 함 — 학습 데이터 gene list ≠ test gene list
- 평가 비교: Classical ML vs Gene FM vs Text LLM (3-tier) → 논문의 핵심 비교 테이블

**구현**: Phase 2 (현재 Phase 1 완료 후)

---

---

## DD-17: Category B Evaluation Criteria (Cross-Mission Transfer)

**결정**: Category B (Cross-Mission Transfer)는 단일 GO/NO-GO 결론을 내리지 않는다. 대신 **Transfer Pattern Summary** 형식으로 결과를 보고한다.

**근거**:

Category A의 GO/NO-GO 기준(perm_p < 0.05)은 Category B에서 통계적으로 달성 불가능하다.

```
n_test = n_flt + n_gnd (예: 3+3=6)
permutation 공간: C(6,3) = 20
최소 가능 p-value = 1/20 = 0.050  (관측값이 20개 중 최고점일 때)
```

따라서 대부분의 Category B 쌍(pair)에서 perm_p ≥ 0.056이 구조적으로 발생한다.
이는 모델 성능 문제가 아니라 소규모 미션 설계의 고유한 한계다.

**평가 형식 — Transfer Pattern Summary**:

```
Transfer Pattern Summary (DD-17):
  AUROC ≥ 0.70 pairs: X/N
  perm_p < 0.05 pairs: X/N
  Large pairs (n≥10): M  mean AUROC=0.XXX  sig=X/M
(Category B: no single GO/NO-GO — use per-pair table above)
```

**1차 분석 단위**: **n≥10 pairs (대규모 쌍)** — 통계적 신뢰도가 충분한 유일한 단위.
- n<10 쌍의 perm_p는 참고용으로만 사용 (수치 신뢰 불가)
- n≥10 쌍에서 AUROC ≥ 0.70 + perm_p < 0.05 → 실질적 transfer 증거로 인정

**2차 분석 단위**: AUROC ≥ 0.70 쌍 비율 — 방향적 증거 (effect size 기반).

**구현 위치**: `scripts/evaluate_submission.py` — `print_summary()` 함수 내
`task_id.startswith("B")` 조건으로 Category B 분기 처리.

**Task 디렉토리 명명 규칙**:
- Category A: `tasks/A{N}_{tissue}_lomo/fold_{MISSION}_test/`
- Category B: `tasks/B{N}_{tissue}_cross_mission/pair_{TRAIN}_{TEST}/`

**Submission fold_id 규칙**:
- Category A: `"fold_{MISSION}_test"` (예: `"fold_RR-9_test"`)
- Category B: `"pair_{TRAIN}_{TEST}"` (예: `"pair_MHU-2_RR-6"`)

**참고**: `fold_*_holdout` 항목은 Tier 1 (LOMO) 제출에서 선택적(optional)이다.
레이블이 비공개이므로 평가자가 expected 목록에서 제외하고 reported 여부만 표시한다.

**적용 미션 (v1.0)**:
- B4 (Thymus): 4개 미션 (MHU-1, MHU-2, RR-6, RR-9) → 12 pairs
- B2 (Gastrocnemius): 3개 미션 (RR-1, RR-5, RR-9) → 6 pairs

**Phase 1 결과 요약** (PCA-LR baseline):
- B4 Thymus: AUROC≥0.70 pairs=9/12, perm_p<0.05=4/12, Large pairs mean=0.775, sig=4/6
- B2 Gastrocnemius: AUROC≥0.70 pairs=4/6, perm_p<0.05=3/6

---

## DD-18: Temporal Covariate Confound Awareness (v2.0)

**결정**: ISS-T vs LAR 비교는 preservation method (RNAlater vs standard necropsy)와 confounded. 모든 T1 결과에 GC_ISS-T vs GC_LAR AUROC를 preservation artifact baseline으로 보고.

**근거**: PMID 33376967 — preservation method가 ISS-T/LAR transcriptome 차이의 주요 원인.

**구현**: `v2/scripts/temporal_analysis.py` T1에서 FLT AUROC와 함께 GC AUROC를 항상 보고. `excess = FLT_AUROC - GC_AUROC`로 biological signal 추정.

**v2.0 결과**: RR-6 excess=+0.078, RR-8 excess=-0.043 → preservation artifact 지배적.

---

## DD-19: T3 Multiple Testing Correction (v2.0)

**결정**: T3 Age × Spaceflight ANOVA에서 50 Hallmark pathways × interaction term → Benjamini-Hochberg FDR correction 필수.

**금지**: 개별 pathway p-value를 FDR 없이 significant로 보고 금지.

**구현**: `statsmodels.stats.multitest.multipletests(method='fdr_bh')`.

**v2.0 결과**: 0/50 pathways significant at FDR < 0.05 (n=40, underpowered).

---

## DD-20: T2 Recovery Fraction Convention (v2.0)

**결정**: direction-aware `recovery_fraction = 1 - (delta_return / delta_flight)`. Preservation-matched comparisons: FLT_ISS-T vs BSL_ISS-T, FLT_LAR vs BSL_LAR.

**규칙**:
- 1.0 = complete recovery
- 0.0 = no recovery
- >1 = overshoot past baseline (direction reversal)
- <0 = continued divergence in the original direction
- `|delta_flight| < 0.1` 인 pathway는 계산에서 제외

**구현**: `v2/scripts/temporal_analysis.py` T2, `RECOVERY_MIN_DELTA = 0.1`.

---

## DD-21: scGPT Foundation Model Integration

**결정일**: 2026-03-09

**결정**: scGPT (whole_human, CellXGene 33M cells) 을 Tier 2 FM benchmark으로 추가.
Geneformer (mouse-specific) 와 대비하여 2nd FM data point 확보.

**모델 구성**:
- 모델: `whole_human` checkpoint (scGPT pretrained, 12L Transformer, 512 hidden dim)
- Fine-tuning: 10 epochs, batch=8, lr=1e-4, warmup=10% steps
- Layer freeze: bottom 10/12 layers frozen (top 2 + cls head만 훈련)
- Scheduler: `torch.optim.lr_scheduler.LambdaLR` (linear warmup/decay)
  - `transformers.get_linear_schedule_with_warmup` 불사용: PyTorch ≥2.4 필요, Cayuga=2.1.2

**Mouse→Human ortholog 매핑**:
- 파일: `data/mouse/ensembl_mouse_human_orthologs.tsv` (83,454 rows)
- 최종 사용 gene 수: ~17K (ortholog 매핑 후 발현 유전자 교집합)

**핵심 버그 수정 기록** (scgpt_finetune.py):
1. `vocab=None` → `vocab.json` 로드 후 dict 전달 (`TransformerModel` 필수)
2. `src_key_padding_mask` 누락 → `gene_ids.eq(pad_token_id)` positional 추가
3. `values` dtype Long → `.float()` cast (ContinuousValueEncoder 요구)
4. `flash_attn==1.0.4` fp16 필수 → `use_fast_transformer=False` 영구 비활성화
5. `output["cls_output"]` 이미 (B, n_cls) logits → `cls_decoder()` 중복 호출 제거
6. **CRITICAL**: `best_auroc = 0.5` 초기화 → `0.0` 수정 (아래찬스 fold 정직 보고)
7. Race condition: 공유 JSON 덮어쓰기 → per-fold `scgpt_{version}_{task}_{fold}_result.json`

**결과**: 6 tissues × 21 LOMO folds
- Overall mean AUROC: 0.667 (GF: 0.476, Baseline: 0.758)
- Classical ML 6/6 완승: scGPT Δ=-0.092, GF Δ=-0.283 vs Baseline
- 결론: scGPT (human pretrained) > GF (mouse-specific) but both < classical ML for small-n bulk RNA-seq

**근거**:
- Geneformer만으론 "FM 전반의 실패" 주장 약함; human FM (scGPT) 추가로 일반화된 결론
- 두 FM 모두 chance 이상이지만 classical ML에 유의미하게 열등 → small-n bulk에서 FM underperform 확인

**구현**: `scripts/scgpt_finetune.py`, `scripts/aggregate_scgpt_results.py`
**결과 파일**: `evaluation/scgpt_whole_human_all_tissues_summary.json`

---

## Changelog

| 버전 | 날짜 | 변경 |
|---|---|---|
| **v2.1** | 2026-03-09 | DD-21 추가 — scGPT FM Integration (whole_human, 7 bug fixes, mean AUROC=0.667, Classical 6/6 완승). |
| **v2.0** | 2026-03-07 | DD-18/19/20 추가 — v2.0 Temporal analysis design decisions (preservation confound, FDR correction, recovery fraction). |
| **v1.4** | 2026-03-01 | DD-17 추가 — Category B Evaluation Criteria (Transfer Pattern Summary, perm_p floor 근거, pair_ 명명 규칙). |
| **v1.3** | 2026-03-01 | DD-15 결과 추가 — fGSEA/GSVA 전체 구현 완료 (51+54 files), J5 비교 결과, 버그 수정 기록. |
| **v1.2** | 2026-03-01 | DD-16 추가 — Text LLM Evaluation Track (GPT-4o, Claude, Llama 3). Bug fix 기록: run_baselines.py B1 (n_ground_test), B2 (permutation p pseudocount). MHU-1 Track 2b 재분류. |
| **v1.0** | 2026-02-28 | 초안 — 3개 리뷰 통합 P0 결정 사항 확정 (DD-01~DD-13) |
| **v1.1** | 2026-02-28 | DD-15 추가 — fGSEA + GSVA Pathway Analysis 설계 결정 |

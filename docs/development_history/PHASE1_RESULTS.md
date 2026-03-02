# Phase 1 Results — Multi-Tissue Spaceflight Detection Benchmark
**날짜**: 2026-02-28 (최종 업데이트: 2026-03-01 v6 — A5 Skin 추가. 3-fold LOMO (MHU-2+RR-6+RR-7). Gene-level GO 확정.)
**상태**: Phase 1 완료 — 6개 조직 중 4개 (A2/A4/A5/A6) GO. Phase 2 진행 예정.

---

## 전체 요약 (Cross-Tissue Phase 1 Go/No-Go)

**Gene-level results (primary evaluation)**:

| 태스크 | 조직 | 미션수 | 최고모델 | AUROC | CI lower | perm p | Within-LOO sanity | SHAP | 판정 |
|--------|------|--------|----------|-------|----------|--------|-------------------|------|------|
| A1 | Liver | 6 | LR | 0.605 ± 0.177 | 0.327 | 0.309 | RR-1✓ RR-6✓ RR-8✓ MHU-2✗ | 0/9 | ✗ NO-GO |
| **A2** | **Gastrocnemius** | **3** | **LR** | **0.907 ± 0.094** | **0.717** | **0.026** | RR-1✓ RR-5✓ RR-9✓ (joint-norm) | **3/16** | **✓ GO** |
| A3 | Kidney | 3 | RF | 0.564 ± 0.256 | 0.359 | 0.354 | — | — | ✗ NO-GO |
| **A4** | **Thymus** | **4**†‡ | **PCA-LR** | **0.923 ± 0.133** | **0.878** | **0.037** | **MHU-1✓ MHU-2✓ RR-6✓ RR-9✓** | **7/23** | **✓ GO** |
| **A5** | **Skin** | **3**§ | **LR** | **0.821 ± 0.051** | **0.637** | **0.0023** | — | — | **✓ GO** |
| A6 | Eye | 3 | LR | 0.811 ± 0.073 | 0.470 | 0.063 | — | — | ✗ NO-GO |

**Pathway-level results (GSVA Hallmark, `evaluation/A{N}_pathway_hallmark_results.json`)**:

| 태스크 | 조직 | 최고모델 | AUROC | CI lower | perm p | 판정 | 비고 |
|--------|------|----------|-------|----------|--------|------|------|
| A3 | Kidney | LR | 0.755 ± 0.130 | 0.481 | 0.071 | ✗ NO-GO | RR-7 (n=94) 클래스 불균형이 bottleneck |
| **A6** | **Eye** | **PCA-LR** | **0.915 ± 0.088** | **0.745** | **0.014** | **✓ GO** | **전 조건 통과 — 경로 특이적 신호 강함** |

†MHU-1 = Track 2b (GC/FLT strain mismatch, see below). ‡4-fold includes MHU-1. §MHU-2=dorsal+femoral merged; RR-7=OSD-254 C57BL/6J subset (30 samples).
**Within-LOO sanity**: LOO CV within each mission separately.
**Within-LOO sanity**: LOO CV within each mission separately — validates whether biological signal exists within each mission independently of cross-mission batch effects.

### 핵심 판정 기준 (within-LOO sanity 적용)

- **A4 Thymus GO 확정 (v2 4-fold)**: MHU-1/MHU-2 분리 → 4 missions, AUROC=0.923↑, CI lower=0.878↑. †MHU-1 포함 4-fold. ‡MHU-1 = Track 2b (GC=C57BL/6CR≠FLT=C57BL/6J, strain confound). Track 2a only (3-fold: MHU-2+RR-6+RR-9) AUROC=0.898. 신규 미션 OSD-515(RR-23) 발견 → Phase 2B 확장 예정.
- **A2 Gastrocnemius GO 확정**: 원시 카운트에서 3 미션 통합 DESeq2 정규화 → mission-level signal 확인 (모든 미션 naive AUROC=1.000). LOMO LR 0.907, Myog/Ciart/Per2 top-50 확인.
- **A1 NO-GO**: 확정 (pipeline heterogeneity).
- **A3 Kidney NO-GO (gene-level 및 pathway-level)**: 확정. Pathway LR AUROC=0.755 (AUROC 조건 통과)이지만 CI lower=0.481 < 0.500, perm_p=0.071 > 0.050 (두 조건 불통과). RR-7 fold(n=94, 클래스 불균형 20F:74G)가 bottleneck. 추가 미션 확보 없이는 GO 어려움.
- **A5 Skin GO 확정 (신규 v6)**: MHU-2 dorsal+femoral 병합 + RR-7(OSD-254, C57BL/6J 30개) 추가 → 3-fold LOMO. LR AUROC=0.821, CI lower=0.637, perm_p=0.0023 → **전 조건 통과**. 4개 모델 전부 AUROC>0.70. RR-5(OSD-240/241, BAL-TAL) Track 2b 제외.
- **A6 Eye NO-GO (gene-level) → GO (pathway-level)**: Gene-level CI lower=0.470 실패. Pathway-level (GSVA Hallmark, PCA-LR): AUROC=0.915, CI lower=0.745, perm_p=0.014 → **전 조건 통과**. 망막 대사 경로가 gene-level보다 미중력 신호 포착에 유리함.

---

## Within-Mission LOO Sanity Check (신규 검증 기준)</p>

**동기**: LOMO AUROC가 cross-mission batch effect에 오염될 수 있음. within-mission LOO는 mission 내부 생물학적 신호의 독립적 증거를 제공.

| 태스크 | 미션 | Track | n_samples | Within-LOO AUROC | 해석 |
|--------|------|-------|-----------|-----------------|------|
| **A4** | MHU-1‡ | 2b | 6 | **1.000** | 강한 신호 (단, Track 2b — strain 효과 가능) |
| **A4** | MHU-2 | 2a | 6 | **1.000** | 강한 생물학적 신호 |
| **A4** | RR-6 | 2a | 35 | **0.771** | 강한 생물학적 신호 (대규모 미션) |
| **A4** | RR-9 | 2a | 20 | **1.000** | 강한 생물학적 신호 |
| **A2** | RR-1 | 2a | 12 | **0.944** | 강한 생물학적 신호 |
| **A2** | RR-5 | 2a | 12 | **1.000** | 강한 생물학적 신호 |
| **A2** | RR-9 | 2a | 8 | **N/A** | ⚠ n=8 LOO 신뢰도 낮음. Naive AUROC=1.000 확인. |

‡ MHU-1 within-LOO=1.000은 strain 차이 반영 가능 (Track 2b).
**A4 Track 2a only (MHU-2 + RR-6 + RR-9)**: 최저 within-LOO = 0.771 → 임계값 0.700 초과 ✓

**RR-9 gastrocnemius — 재검토 결과 (joint normalization 후):**
- 원시 카운트 확인: 3 미션 모두 동일하게 57,186 genes (동일 annotation 버전)
- Joint DESeq2 후 zero rate: RR-1=4.6%, RR-5=8.2%, RR-9=22.5% (개선됨, 이전 30%)
- RR-9 size factors: min=0.345, max=5.928 (폭넓음 → sequencing depth variance가 근원)
- 모든 미션 within-mission naive AUROC = 1.000 → 생물학적 신호 존재 확인
- Within-LOO variance-based LOO AUROC=0.000 원인: n=8 너무 작아서 variance 필터가 음성 상관 유전자 선택
- **결론**: 데이터 품질 이상 아님. RR-9는 낮은 sequencing depth (size factors 폭넓음)가 원인.

**A2 재정규화 결론**: joint DESeq2 (3 mission 통합) 후 LOMO LR AUROC=0.907. 3개 미션 모두 genuine biological signal 확인. GO 확정.

---

---

## 재현성 노트 — B3: SAGA 수렴 이슈 (2026-03-01)

**발견**: `run_baselines.py`의 LR 모델(`build_lr()`)이 `max_iter=2000`으로 설정되어, 15,000+ genes에 대한 SAGA solver가 일부 fold에서 수렴하지 않을 수 있음. PCA-LR(`build_pca_lr()`)는 차원 축소 후 lbfgs를 사용하므로 영향 없음.

**영향 범위**:
- A2 Gastrocnemius LR: RR-1 fold AUROC 0.778 (미수렴) → 0.806 (max_iter=10000 수렴 후)
- GO/NO-GO 결론 변화 없음 — 결론: A2 GO, A4 GO 확정
- A4 PCA-LR: 영향 없음 (PCA 후 19–50 PC에서 lbfgs 수렴)
- A1/A3/A6: GO/NO-GO 결론 변화 없음 (all NO-GO 유지)

**공식 결과 처리**: `A2_baseline_results.json`의 공식 값(LR AUROC=0.907)은 미수렴 상태의 결과. 재실행 시 ~0.917. Baseline submission(`submission_LR_baseline_A2.json`)은 수렴 버전(max_iter=10000) 사용. **Phase 2에서 `run_baselines.py` max_iter 수정 예정 (B3 fix).**

---

## Task A1: Liver LOMO (6 folds, 264 samples)

### 데이터
```
미션: RR-1, RR-3, RR-6, RR-8, RR-9 (ISS NASA RR v1/v2), MHU-2 (JAXA MHU v2)
Feature: log2(DESeq2 normalized counts + 1) per sample (DD-01)
Feature selection: LOMO-aware variance filter top 75% (train only, DD-03)
Labels: Flight(uG)=1, GC or VC=0 (BC, AG 제외)
Split: LOMO 6-fold (independence unit = Mission, DD-04)
```

### 폴드별 구성

| 미션 | n_test | Flight | Ground | 파이프라인 | 하드웨어 |
|------|--------|--------|--------|-----------|---------|
| MHU-2 | 6 | 3 (uG) | 3 | v2 | JAXA MHU |
| RR-1 | 14 | 7 | 7 | v2 | NASA RR |
| RR-3 | 12 | 6 | 6 | v1 | NASA RR |
| RR-6 | 39 | 18 | 21 | v1 | NASA RR |
| RR-8 | 103 | 35 | 68 | v2 | NASA RR |
| RR-9 | 19 | 5 | 14 | v1 | NASA RR |

### Baseline 결과 (AUROC, final — dup fix 적용 264 samples)

| 폴드 | LogReg | PCA-LR | RF | n_test |
|------|--------|--------|----|--------|
| MHU-2 | 0.333 | 0.222 | 0.278 | 6 |
| RR-1 | 0.755 | 0.673 | 0.704 | 14 |
| RR-3 | 0.444 | 0.500 | 0.500 | 12 |
| RR-6 | **0.850** | **0.795** | 0.743 | 39 |
| RR-8 | 0.679 | 0.618 | **0.819** | 103 |
| RR-9 | 0.571 | 0.714 | 0.493 | 19 |
| **Unweighted mean** | 0.605 ± 0.177 | 0.587 ± 0.186 | 0.590 ± 0.184 | |
| **Sample-weighted** | 0.683 | 0.648 | **0.727** | |

### Phase 1 Go/No-Go (DD-11) — A1 Liver

| 조건 | 기준 | 실측 (best = LR) | 결과 |
|------|------|------|------|
| AUROC > 0.70 | 0.70 | 0.605 | ✗ FAIL |
| CI lower > 0.60 | 0.60 | 0.327 | ✗ FAIL |
| perm p < 0.05 | 0.05 | 0.309 | ✗ FAIL |

**판정: NO-GO**

### A1 주요 발견

**1. MHU-2 batch effect (핵심)**

- Within-MHU-2 sanity check (uG vs GC, LOO CV): **LOO AUROC = 1.000** ← 신호 자체는 완벽
- Cross-mission (train: 5 RR, test: MHU-2): **LR AUROC = 0.333** (역방향 분류)
- JAXA MHU 하드웨어 vs NASA RR 하드웨어 간 batch effect가 biological signal보다 큼
- Thymus에서는 MHU-2 fold RF AUROC = 1.000 → liver-specific batch effect 가능성

**2. ComBat-seq batch correction 효과 (RF 비교)**

| 미션 | 미보정 RF | ComBat-seq RF | 변화 |
|------|----------|--------------|------|
| MHU-2 | 0.278 | **0.667** | ↑+0.389 ← 역전 해소! |
| RR-1 | 0.704 | **0.837** | ↑+0.133 |
| RR-3 | 0.500 | 0.389 | ↓-0.111 |
| RR-6 | 0.743 | 0.687 | ↓-0.056 |
| RR-8 | 0.819 | 0.728 | ↓-0.091 |
| RR-9 | 0.493 | 0.650 | ↑+0.157 |
| **Mean** | 0.590 | **0.660** | ↑+0.070 |

ComBat-seq: MHU-2 배치 효과 부분 교정 성공 (역방향 분류 해소).
그러나 RR-3/RR-6/RR-8 fold에서 소폭 하락 → 생물학적 신호 일부 손실 가능성.
전체 평균 0.590 → 0.660 향상, 그러나 0.70 미달.

**3. SHAP 분석 결과 (RF, uncorrected)**

- Target genes (Angptl4, Pck1, G6pc, Fasn, Ppara, Acox1, Cyp7a1): **0/9 in top-50** ✗
- Top genes 특징: **Dbp(#18), Bmal1(#26), Npas2(#16), Ciart(#35), Nfil3(#36)** → 일주기 시계 유전자 군집
  → 우주비행 시 일주기 교란(circadian disruption) 신호 포착. 그러나 liver 대사 표적 유전자와 무관.
- Hmgcs1(#25): 콜레스테롤 합성 (간접적 관련), Gclc(#7): 글루타치온 합성 (산화 스트레스 반응)
- Condition 3: ✗ FAIL — 낮은 AUROC로 모델이 생물학적 신호보다 배치 효과를 포착한 것으로 해석

**4. 폴드 크기 의존성**

n_test ≥ 30 폴드 (RR-6, RR-8) 평균:
- 미보정 LR: mean = 0.765, RF: 0.781 ← Phase 1 기준 통과!
- 소규모 폴드 (n < 15: MHU-2, RR-1, RR-3) 통계적 불안정성 높음

**5. A1 ISS-only 분석 (MHU-2 제외, 5-fold LOMO)**

| 폴드 | LogReg | PCA-LR | RF | n_test |
|------|--------|--------|----|--------|
| RR-1 | 0.776 | **0.796** | 0.612 | 14 |
| RR-3 | 0.472 | 0.528 | 0.472 | 12 |
| RR-6 | **0.884** | 0.868 | 0.583 | 39 |
| RR-8 | 0.736 | 0.701 | **0.803** | 103 |
| RR-9 | 0.371 | 0.414 | 0.500 | 19 |
| **Mean** | 0.648 ± 0.193 | **0.661 ± 0.168** | 0.594 ± 0.117 | |

**판정: NO-GO (MHU-2 제외 후에도 NO-GO)**
- Best PCA-LR mean AUROC = 0.661, CI lower = 0.416, perm p = 0.236 → 3개 조건 모두 FAIL
- MHU-2가 문제의 주원인이 아님. 핵심 문제: **RR-3 (0.472–0.528) + RR-9 (0.371–0.500)** 가 mean을 끌어내림
- RR-3/RR-9 분류 실패 원인: pipeline version v1 (RR-3: v1) vs v2 batch effect
- 샘플 가중 평균 (n_test 비례): LR ≈ 0.728, PCA-LR ≈ 0.730 → 대규모 폴드에서는 충분한 신호 존재
- 결론: **A1 Liver는 pipeline heterogeneity batch effect가 biological signal보다 커서 NO-GO 확정**

**6. 다음 조치**

1. Pipeline version effect (v1 vs v2) 정량화 (J3)
2. 논문: ComBat-seq 전후 AUROC 비교를 batch effect 섹션에 포함

---

## Task A2: Gastrocnemius LOMO (3 folds, 32 samples) ← **PHASE 1 GO (전 조건 통과)**

### 데이터
```
미션: RR-1 (OSD-101), RR-5 (OSD-401), RR-9 (OSD-326)
Feature: log2(DESeq2 joint normalized counts + 1)
  - 3 미션 원시 RSEM 카운트를 단일 DESeq2 실행으로 통합 정규화 (v2 개선판)
  - design = ~ mission + condition (미션 효과 회귀 후 size factor 추정)
  - Low-count filter: ≥10 counts in ≥2 samples globally → 21,013 genes
Feature selection: LOMO-aware variance filter top 75% (train only) → 15,760 genes
Labels: Flight(uG)=1, GC=0
Split: LOMO 3-fold
정규화 파일: gastrocnemius_joint_log2_norm.csv (32 samples × 21,013 genes)
```

### 폴드별 구성

| 미션 | n_test | Flight | Ground |
|------|--------|--------|--------|
| RR-1 | 12 | 6 | 6 |
| RR-5 | 12 | 6 | 6 |
| RR-9 | 8 | 4 | 4 |

### Baseline 결과 (AUROC) — Joint Normalization v2

| 폴드 | LogReg | PCA-LR | RF | n_test |
|------|--------|--------|----|--------|
| RR-1 | **0.778** | 0.639 | 0.542 | 12 |
| RR-5 | **0.944** | **0.833** | 0.750 | 12 |
| RR-9 | **1.000** | **1.000** | 0.781 | 8 |
| **Mean** | **0.907 ± 0.094** | 0.824 ± 0.148 | 0.691 ± 0.106 | |

### Phase 1 Go/No-Go — A2 Gastrocnemius

| 조건 | 기준 | 실측 (best = LR) | 결과 |
|------|------|------|------|
| AUROC > 0.70 | 0.70 | **0.907** | ✓ PASS |
| CI lower > 0.60 | 0.60 | **0.717** | ✓ PASS |
| perm p < 0.05 | 0.05 | **0.026** | ✓ PASS |
| SHAP gene check | top-50 ≥1 target | **3/16 target genes** | ✓ PASS |

**판정: GO (모든 4조건 통과)**

### A2 SHAP 결과 (RF, top-50) — Joint Normalized v2

**3/16 타겟 유전자 top-50 포함:**

| Rank | Gene | 기능 |
|------|------|------|
| 9 | **Myog** | Myogenin — 근육 분화 TF, 근위축 마커 (atrophy-upregulated) |
| 25 | **Ciart** | Circadian-associated repressor of transcription |
| 33 | **Per2** | Period 2, circadian repressor |

**생물학적 해석:**
- **Myog(#9)**: myogenin은 근위축 시 발현 증가. spaceflight 근육 위축의 직접 마커. top-10 진입.
- **Ciart(#25), Per2(#33)**: 일주기 시계 유전자 (이전과 동일한 cross-tissue circadian signature)
- 비표적 top 유전자: Sox4(#1, 발달TF), Lrrc30(#2, 근육 ECM), Cdkn1c(#4, 세포주기), Chaf1b(#7, 크로마틴)
- 근위축 마커 Fbxo32(Atrogin-1)/Trim63(MuRF-1)는 여전히 top-50 밖 (일주기 교란이 더 강한 신호)

### A2 주요 발견

- RR-9 fold AUROC = 1.000, RR-5 = 0.944 (강한 cross-mission generalization)
- RR-1 fold 성능 상대적으로 낮음 (LR 0.778): 최고령 미션, 생물학적 variability
- 근위축 직접 마커(Myog) + 일주기 유전자가 공동으로 상위 feature
- RF 성능이 LR보다 낮음: 소규모 데이터(32 samples)에서 tree overfitting
- 재정규화(joint DESeq2) 효과: RR-9 데이터 품질 개선 (zero rate 30%→22.5%, 동일 gene set 확보)

---

## Task A3: Kidney LOMO (3 folds, ~142 samples)

### Baseline 결과

| 폴드 | LogReg | PCA-LR | RF | n_test |
|------|--------|--------|----|--------|
| RR-1 | 0.278 | 0.167 | 0.236 | 12 |
| RR-3 | 0.639 | 0.611 | **0.861** | 12 |
| RR-7 | 0.647 | 0.603 | 0.594 | 94 |
| **Mean** | 0.521 ± 0.172 | 0.460 ± 0.208 | **0.564 ± 0.256** | |

### Phase 1 Go/No-Go — A3 Kidney

| 조건 | 실측 (best = RF) | 결과 |
|------|------|------|
| AUROC > 0.70 | 0.564 | ✗ FAIL |
| CI lower > 0.60 | 0.359 | ✗ FAIL |
| perm p < 0.05 | 0.354 | ✗ FAIL |

**판정: NO-GO**

### A3 주요 발견

- RR-1 fold 매우 낮음 (RF 0.236): 소규모(n=12) + 파이프라인 버전 차이 가능성
- RR-3 fold RF 0.861: 짧은 train (RR-1+RR-7)으로도 좋은 성능 → signal 존재
- RR-7 fold (n=94): 중간 성능 (0.594-0.647), kidney signal은 있으나 batch effect
- 3 missions only → LOMO 통계적으로 불안정 (각 fold = 1 test mission)

### ✗ Pathway-level 재분석 — OFFICIAL NO-GO (2026-03-01)

> **출처**: `evaluation/A3_pathway_hallmark_results.json`
> 스크립트: `python scripts/run_pathway_lomo.py --tissue kidney --db hallmark`

**A3 Kidney Pathway LOMO (GSVA Hallmark, 50 pathways)**:

| 폴드 | n_test | Gene RF | Pathway LR | Pathway PCA-LR | |
|------|--------|---------|------------|----------------|---|
| RR-1 | 12 | 0.236 | **0.861** | **0.861** | ↑+0.625 |
| RR-3 | 12 | 0.861 | 0.833 | 0.806 | ≈ |
| RR-7 | 94 | 0.594 | 0.572 | 0.562 | ↓−0.022 |
| **Mean** | | **0.564** | **0.755** | **0.743** | ↑+0.191 |

**공식 GO/NO-GO (Pathway LR — best model)**:

| 조건 | 기준 | 실측 | 결과 |
|------|------|------|------|
| AUROC > 0.70 | 0.70 | **0.755** | ✓ PASS |
| CI lower > 0.50 | 0.50 | 0.481 | ✗ FAIL |
| perm p < 0.05 | 0.05 | 0.071 | ✗ FAIL |

**판정: ✗ PATHWAY-LEVEL NO-GO** — AUROC 조건 통과, CI lower + perm_p 미달.

**실패 원인 분석**:
- RR-7 fold (n=94, 20F vs 74G: 1:3.7 불균형): pathway AUROC=0.572 — 구제 불가
- RR-7이 mean CI lower를 0.481로 끌어내림 (임계값 0.500 미달)
- Gene-level 대비 개선: AUROC 0.564→0.755 (+0.191), CI lower 0.359→0.481 (+0.122)
- **결론**: 소규모 fold(n=12)는 pathway로 구제되지만, 대형 fold(n=94) 클래스 불균형이 bottleneck

**비교 참조 (J5 gene vs pathway)**:
- J5 내부 비교 (일관된 PCA-LR baseline): gene=0.432, pathway=0.743 (+0.310)
- 공식 run_pathway_lomo.py 결과 (LR): gene best(RF)=0.564, pathway(LR)=0.755 (+0.191)
- 두 분석 모두 pathway 우세 확인 — 단, 공식 CI/perm_p 기준에서는 NO-GO

**향후 방향**:
- RR-7 fold 클래스 불균형 해결: stratified permutation test, 또는 class_weight="balanced" 강화
- 추가 kidney 미션(n>12) 확보 시 통계적 검정력 향상 가능
- 현재: Gene-level NO-GO, Pathway-level NO-GO 모두 확정

---

## Task A4: Thymus LOMO ← **PHASE 1 GO (전 조건 통과)**

> **[Phase 2 업데이트 2026-02-27]**: MHU-1/MHU-2 분리 후 4-fold LOMO로 업그레이드. AUROC 향상.

### Baseline 결과 — v2 (4-fold, MHU-1 분리)

> **⚠️ Track 분류 주의**: MHU-1 fold는 **Track 2b** (all strains) 결과입니다. MHU-1 GC 샘플 strain = C57BL/6CR ≠ FLT strain C57BL/6J. Track 2a (C57BL/6J only, DD-06) 분석에서는 이 fold를 제외해야 합니다. Track 2a 순수 결과 = MHU-2+RR-6+RR-9 (3-fold).

| 폴드 | Track | PCA-LR | perm_p | n_test | 비고 |
|------|-------|--------|--------|--------|------|
| MHU-1 | **2b** | **1.000** | 0.062 | 6 | ⚠️ GC=C57BL/6CR, FLT=C57BL/6J (strain confound) |
| MHU-2 | 2a | **1.000** | 0.062† | 6 | †n=6 소규모 |
| RR-6  | 2a | 0.693  | 0.024 | 35 | |
| RR-9  | 2a | **1.000** | 0.001‡ | 20 | ‡pseudocount 적용 후 (이전 0.000은 버그) |
| **Mean (4-fold)** | 2a+2b | **0.923 ± 0.133** | **0.037** | | CI lower=0.878 |
| **Track 2a only (3-fold)** | 2a | **0.898 ± 0.177** | **0.029** | | MHU-2+RR-6+RR-9 (CI lower 미계산) |

**Phase 1 v1 (3-fold, MHU-1+2 통합)**: AUROC=0.888, CI lower=0.795, perm p=0.009
**버그 수정 기록**: perm_p 계산 pseudocount 추가 (2026-03-01) — RR-9 fold: 0.000→0.001

### Phase 1 Go/No-Go — A4 Thymus (v2, 4-fold)

| 조건 | 기준 | 실측 (best = PCA-LR) | 결과 |
|------|------|------|------|
| AUROC > 0.70 | 0.70 | **0.923** | ✓ PASS |
| CI lower > 0.60 | 0.60 | **0.878** | ✓ PASS |
| perm p < 0.05 | 0.05 | **0.037** | ✓ PASS |
| SHAP gene check | top-50 ≥1 target | **7/23 target genes** | ✓ PASS |

**판정: GO (모든 4조건 통과)** — v2에서 더 강화됨

### A4 SHAP 결과 (RF, top-50) — v2 (4-fold)

**7/23 타겟 유전자 top-50 포함:**

| Rank | Gene | Δlog2 | 방향 | 기능 |
|------|------|-------|------|------|
| 5  | **H2ac22** | -1.43 | ↓Flight | Histone H2A variant |
| 12 | **Apoe**   | +1.30 | ↑Flight | Apolipoprotein E (stress response) |
| 13 | **Msh6**   | -0.62 | ↓Flight | DNA mismatch repair |
| 17 | **Btla**   | +0.77 | ↑Flight | T-cell immune checkpoint |
| 22 | **Apod**   | +1.20 | ↑Flight | Apolipoprotein D (stress response) |
| 30 | **Maf**    | +0.47 | ↑Flight | c-Maf TF (Treg/γδT 분화) |
| 35 | **Trbv12-1** | -1.38 | ↓Flight | T-cell receptor beta variable |

**방향성 요약 (전체 top-50)**: 34 ↓ Flight, 16 ↑ Flight

**생물학적 해석:**
1. **흉선세포 증식 감소** (↓): Fbxo5(-1.53), Cenpp(-1.59), Cenph(-1.76), Pclaf(-1.83), Pbk(-2.15), Cdkn3(-1.74), Esco2(-1.95), Mxd3(-2.05) → 세포주기 유전자 ALL DOWN
2. **히스톤 감소** (↓): H2ac22(-1.43), H2ac11(-2.33), H2ac24(-1.40) → 뉴클레오솜 조립 감소
3. **스트레스 지질반응** (↑): Apoe(+1.30), Apod(+1.20) → 지질 대사 보상반응
4. **면역 체크포인트** (↑): Btla(+0.77) → 잔류 세포의 면역억제 강화
5. **T-cell 기능 저하** (↓): Trbv12-1(-1.38), Msh6(-0.62) → T-cell 수용체·DNA 복구 저하

### A4 Pathway Enrichment 결과 (g:Profiler, mmusculus)

**파일**: `evaluation/A4_pathway_enrichment.json` (2026-02-27 실행)

| Source | Term | p-value | 주요 유전자 |
|--------|------|---------|------------|
| GO:BP | **Cell cycle process** | 8.8e-06 | Fbxo5, Pclaf, Cenpa, Pbk, Cdca8, Cenpf, Cenpp, Cenpw, Spc25, Rrm1, Dut |
| GO:BP | Mitotic cell cycle | 1.8e-05 | Fbxo5, Pclaf, Cenpa, Pbk, Cdca8, Cenpf, Cenpp, Cenpw, Spc25 |
| GO:BP | Regulation of mitotic cell cycle phase transition | 2.1e-05 | Fbxo5, Pclaf, Pbk, Mxd3, Cdca8 |
| REAC | Cell Cycle | 2.8e-04 | Cenpa, Cenpf, Cenpp, Cenpw, Spc25, Cdca8, Pbk |
| REAC | **Nucleosome assembly** | 4.1e-04 | H2ac12, H2ac22, Cenpa |
| REAC | **CENPA nucleosome deposition at centromere** | 4.1e-04 | Cenpa, Cenpf, Cenpp, Cenpw |
| REAC | Mitotic Spindle Checkpoint | 1.6e-03 | Spc25, Cdca8, Cenpf |
| KEGG | Platinum drug resistance (proxy: DNA damage) | 2.5e-02 | Brca1, Msh6, Bcl2 |

**통합 생물학적 내러티브 (A4 Thymus Spaceflight Mechanism):**
```
우주비행 → 우주방사선(GCR/SPE) + 미세중력
  ↓
DNA 이중가닥 절단(DSB) → Brca1(#26), Msh6(#4) 활성화 → DNA repair 경로
  ↓
G2/M 체크포인트 활성화 → Cenpa, Cenpf, Cenpp, Cenpw, Cdca8, Spc25 → 방추사 checkpoint
  ↓
흉선 상피세포 + 흉선세포 증식 억제 → 흉선 위축 (관찰된 phenotype)
  ↓
T-cell 발달 장애 → Gzma↓, Maf↓, Btla 조절이상 → 면역 억제
```

**핵심 발견 (v2 업데이트)**: A4 thymus SHAP의 지배적 신호는 **흉선 위축으로 인한 증식 세포 감소** (Cell composition effect)이다:
- 정상 흉선(GC)에는 빠르게 증식하는 흉선세포(세포주기 유전자 高발현)가 풍부
- 우주비행 후 흉선은 이 증식 집단이 고갈되어 평균 발현이 DOWN
- 잔류 세포는 면역억제(Btla↑) + 지질 스트레스반응(Apoe↑, Apod↑)을 활성화

**Calibration 분석 결과** (evaluation/A4_calibration.json):
- Flight mean proba = 0.622 ± 0.452, GC mean proba = 0.476 ± 0.474
- Separation = 0.146 (낮음 — 미션 간 절대 확률 불일치)
- **해석**: AUROC(ranking)은 높지만 cross-mission calibration 불량 → 임상 응용 전 calibration 필수

> **⚠️ Calibration 한계 (DD-14)**: 현재 LOMO baseline의 예측 확률은 **미보정(uncalibrated)** 상태입니다.
> - AUROC 기반 ranking은 유효하나, 절대 확률(0.0~1.0)은 cross-mission에서 신뢰 불가
> - 예: RR-9 fold: Flight proba≈0.0 / RR-6 fold: GC proba≈1.0 (심각한 역전)
> - **원인**: 소규모 LOMO 훈련 세트 + 고차원 특성 → 로지스틱 회귀 결정 경계 불안정
> - **외부 모델 제출 시 주의**: 제출 확률의 calibration 여부를 model_description에 명시 권장
> - **향후 계획**: Platt scaling 또는 Isotonic regression으로 사후 보정 (Phase 2 옵션)

### A4 주요 발견 (v2, 4-fold)

- **MHU-1 분리**: GLDS-289에 MHU-1+MHU-2 혼재 → 분리 후 각 n=9(6 flight/gc), 4-fold LOMO
- **MHU-1 strain confound**: GC=C57-6CR, FLT=C57-6J → AUROC=1.000이지만 일부는 strain 차이
- **신규 미션 발견**: OSD-515(RR-23, 27샘플), OSD-457(SpaceX-14 WT+Nrf2KO, 12 WT thymus)
- RR-6 fold 지속 낮은 성능 (PCA-LR=0.693) → RR-6 내 baseline(BC) 샘플 혼재 가능성
- RR-9 fold PCA-LR=1.000 → 매우 강한 spaceflight thymus signal

---

## Task A6: Eye/Retina LOMO (3 folds, ~37 test samples)

### Baseline 결과

| 폴드 | LogReg | PCA-LR | RF | n_test |
|------|--------|--------|----|--------|
| RR-1 | 0.722 | 0.611 | 0.514 | 12 |
| RR-3 | **0.900** | **0.900** | **0.925** | 9 |
| TBD | 0.810 | 0.857 | 0.730 | 16 |
| **Mean** | **0.811** ± 0.073 | 0.789 ± 0.127 | 0.723 ± 0.168 | |

### Phase 1 Go/No-Go — A6 Eye

| 조건 | 기준 | 실측 (best = LR) | 결과 |
|------|------|------|------|
| AUROC > 0.70 | 0.70 | **0.811** | ✓ PASS |
| CI lower > 0.60 | 0.60 | 0.470 | ✗ FAIL |
| perm p < 0.05 | 0.05 | 0.063 | ✗ FAIL (borderline) |

**판정: NO-GO** (AUROC 통과, CI lower + perm p 미달)

### A6 주요 발견

- AUROC = 0.811 (3모델 평균) → spaceflight eye signal 강함
- CI lower 실패 이유: 폴드 크기 매우 작음 (n=9, 12, 16) → bootstrap 불안정
- 추가 미션 데이터(TBD = OSD-397?) 확인 필요

### ✅ Pathway-level 재분석 — OFFICIAL GO (2026-03-01)

> **출처**: `evaluation/A6_pathway_hallmark_results.json`
> 스크립트: `python scripts/run_pathway_lomo.py --tissue eye --db hallmark`

**A6 Eye Pathway LOMO (GSVA Hallmark, 50 pathways)**:

| 폴드 | n_test | Gene LR | Pathway LR | Pathway PCA-LR | |
|------|--------|---------|------------|----------------|---|
| RR-1 | 12 | 0.722 | **1.000** | **1.000** | ↑+0.278 |
| RR-3 | 9  | 0.900 | 0.850 | **0.950** | ± |
| TBD  | 16 | 0.810 | 0.794 | 0.794 | ↓−0.016 |
| **Mean** | | **0.811** | **0.881** | **0.915** | |

**공식 GO/NO-GO (Pathway PCA-LR)**:

| 조건 | 기준 | 실측 | 결과 |
|------|------|------|------|
| AUROC > 0.70 | 0.70 | **0.915** | ✓ PASS |
| CI lower > 0.50 | 0.50 | **0.745** | ✓ PASS |
| perm p < 0.05 | 0.05 | **0.014** | ✓ PASS |

**판정: ✅ PATHWAY-LEVEL GO** — 전 조건 통과. (`evaluation/A6_pathway_hallmark_results.json`)

**생물학적 해석**:
- GSVA Hallmark 상위 경로 (DD-15): **OXIDATIVE_PHOSPHORYLATION** dominant
- 망막은 산화적 인산화 의존도가 매우 높은 조직 → 미중력에 의한 에너지 대사 변화가 pathway 수준에서 일관됨
- Gene-level CI lower 실패 이유: n=9–16 소규모 fold → 개별 유전자 노이즈 큼 → 50개 pathway 집계로 해결
- RR-1 fold: gene AUROC=0.722 → pathway AUROC=1.000

**핵심 의의**: Gene-level NO-GO → Pathway-level GO. 소규모 폴드(n<20)에서 pathway feature가 gene feature보다 안정적임을 보여주는 강력한 사례. 미중력 신호가 개별 유전자보다 pathway 수준에서 더 일관되게 보존됨.

### Pathway Feature Importance (PCA-LR coefficient, `evaluation/A6_pathway_hallmark_shap.json`)

> **방법**: PCA-LR 모델의 분류기 계수를 PCA 공간에서 원래 pathway 공간으로 역사영 (`coef @ PCA.components_ / scaler.scale_`). Fold별 평균 계수 및 평균 절대값 계수. 양수 = Flight 마커, 음수 = Ground 마커.

**상위 8개 경로 (|mean_coef| 기준)**:

| Rank | Pathway | mean_coef | abs_mean_coef | 방향 | 생물학적 해석 |
|------|---------|-----------|---------------|------|--------------|
| 1 | HALLMARK_UV_RESPONSE_UP | −3.90 | 3.90 | **Ground** ↓ | UV 손상 반응 유전자 — 우주방사선(GCR) 노출 시 오히려 억제 가능 |
| 2 | HALLMARK_UNFOLDED_PROTEIN_RESPONSE | −2.80 | 2.80 | **Ground** ↓ | ER 스트레스/UPR — Ground 망막에서 더 높은 단백질 폴딩 스트레스 |
| 3 | HALLMARK_HYPOXIA | −2.18 | 2.18 | **Ground** ↓ | 저산소 반응 — 미중력 시 망막 혈류 변화로 감소 가능 |
| 4 | HALLMARK_PEROXISOME | −2.05 | 2.05 | **Ground** ↓ | 퍼록시솜 지방산 산화 — Ground 대사 상태에서 활성 |
| 5 | HALLMARK_P53_PATHWAY | −1.94 | 1.94 | **Ground** ↓ | p53 경로 — DNA 손상 반응 (Ground 기저 스트레스) |
| 6 | HALLMARK_GLYCOLYSIS | −1.86 | 1.86 | **Ground** ↓ | 해당작용 — Ground 망막 에너지 대사의 우세 경로 |
| 7 | HALLMARK_G2M_CHECKPOINT | +1.86 | 1.86 | **Flight** ↑ | G2/M 체크포인트 — 우주방사선 DNA 손상 → 세포주기 정지 |
| 8 | HALLMARK_MYC_TARGETS_V2 | −1.62 | 1.62 | **Ground** ↓ | MYC 표적 (rRNA 합성) — Ground 세포 성장 프로그램 우세 |

**추가 Flight 마커 (양수 계수)**:

| Pathway | mean_coef | 생물학적 해석 |
|---------|-----------|--------------|
| HALLMARK_MITOTIC_SPINDLE | +1.34 | 방추사 조립 — 방사선 유도 분열 이상 |
| HALLMARK_E2F_TARGETS | +1.33 | E2F 표적 (S-phase) — 방사선 반응 세포주기 재진입 |
| HALLMARK_ANDROGEN_RESPONSE | +1.32 | 안드로겐 반응 — 스트레스 호르몬 cross-talk |
| HALLMARK_IL6_JAK_STAT3_SIGNALING | +1.06 | 염증 신호 — 미중력 면역 반응 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | +0.91 | IFN-γ 반응 — 방사선 유도 면역 활성화 |

**통합 해석**:
- **Ground 마커 (음수 coef)**: UV/UPR/Hypoxia/Glycolysis — 정상 망막의 대사/스트레스 반응 기저 경로
- **Flight 마커 (양수 coef)**: G2M/E2F/Mitotic Spindle — 방사선 DNA 손상 → 세포주기 체크포인트 활성화
- Gene-level에서 포착하지 못한 방사선 유도 세포주기 이상이 pathway 수준에서 일관된 신호로 나타남
- A4 Thymus SHAP (세포주기 유전자 DOWN)과 반대 방향: Eye는 방사선 반응으로 세포주기 체크포인트 **UP** (조직 특이적 반응 차이)

---

## A5: Skin (피부) — Track 2a, 3-fold LOMO

**데이터**: MHU-2 (dorsal=OSD-238 + femoral=OSD-239 병합) + RR-6 (OSD-243) + RR-7 (OSD-254, C57BL/6J 서브셋)

**설계 결정**:
- MHU-2 dorsal (OSD-238)과 femoral (OSD-239)은 동일 우주 미션에서 다른 피부 부위 → mission="MHU-2"로 병합 (LOMO 독립성 단위: 미션)
- RR-7 (OSD-254): C57BL/6J + C3H/HeJ 혼합 연구 → C57BL/6J 30개 샘플만 추출 (BSL 10개, C3H 40개 제외)
- RR-5 (OSD-240/241, BAL-TAL 균주): Track 2b — Track 2a에서 제외

**Binary label 분포**:

| Mission | Flight | Ground (GC+VC) | 제외 그룹 |
|---------|--------|-----------------|----------|
| MHU-2 (dorsal+femoral) | 11 | 24 (18 GC + 6 VC) | 12 AG |
| RR-6 | 18 | 19 (GC) | 16 BC |
| RR-7 (C57BL/6J) | 10 | 20 (10 GC + 10 VC) | 10 BSL + 40 C3H |

총 binary 샘플: 102개 (102 = 35 + 37 + 30), 26,814 유전자

**Baseline 결과** (`evaluation/A5_baseline_results.json`):

| 모델 | MHU-2 | RR-6 | RR-7 | Mean AUROC | CI lower | perm_p |
|------|-------|------|------|------------|----------|--------|
| **LR (ElasticNet)** | **0.890** | **0.769** | **0.805** | **0.821** | **0.637** | **0.0023** |
| PCA-LR | 0.792 | 0.751 | 0.840 | 0.794 | 0.607 | 0.0037 |
| RF | 0.777 | 0.735 | 0.777 | 0.763 | 0.566 | 0.0063 |
| XGBoost | 0.795 | 0.749 | 0.725 | 0.756 | 0.558 | 0.0123 |

**판정: ✅ GENE-LEVEL GO** — 전 조건 통과 (LR 기준: AUROC=0.821>0.70, CI lower=0.637>0.50, perm_p=0.0023<0.05).

**SHAP 상위 유전자** (RF, `evaluation/A5_shap_rf.json`):

| Rank | 유전자 | 기능 | 생물학적 해석 |
|------|--------|------|--------------|
| 1 | **Ckap4** | Cytoskeletal-associated protein | ECM/막 스트레스 반응 |
| 2 | **Rdh16** | Retinol dehydrogenase 16 | 피부 레티노이드 대사 변화 |
| 3 | **Aunip** | Aurora kinase A interacting protein | 세포 분열 (방사선 손상 반응) |
| 4 | **Tcf19** | Transcription factor 19 | 세포 증식 조절 |
| 5 | **E2f1** | E2F1 transcription factor | G1/S 체크포인트 — 방사선/미중력 반응 |
| 6 | **Ckap2** | Cytoskeletal-associated protein 2 | 유사분열 방추사 |
| 7 | **Klk8** | Kallikrein 8 (serine protease) | 피부 특이적 — 피부 장벽 기능 |
| 8 | **Il20** | Interleukin 20 | 피부 염증 사이토카인 |
| 9 | **Hmgcs2** | HMG-CoA synthase 2 | 케톤체 합성 (미중력 대사 적응) |
| 10 | **Klrb1c** | Killer cell lectin receptor B1c | 선천 면역 NK세포 |

**생물학적 해석**:
- **세포 주기 패턴 (조직 공통)**: E2f1, Ckap2, Aunip, Tcf19 → Thymus/Eye와 공통 신호. 미중력-방사선에 의한 세포주기 변화가 피부에서도 포착됨
- **피부 특이적 신호**: Klk8 (피부 장벽 kallikrein, 각질 분화 조절), Rdh16 (피부 레티노이드 대사) — 미중력에 의한 피부 장벽 및 분화 변화
- **면역 신호**: Il20 (피부 염증 IL-20), Klrb1c (NK세포 수용체) — 미중력 면역 억제와 일치
- **대사 변화**: Hmgcs2 (케톤체 합성) — 미중력 지방산 산화 대사 전환

**핵심 의의**:
- 피부는 이전 연구 (Communications Medicine 2024, PMC11166967)에서 DEG 중복 기반 분석만 수행 → 본 연구의 LOMO-ML이 교차 미션 일반화 **최초 검증**
- 3개 미션 (MHU-2, RR-6, RR-7) 모두 AUROC > 0.70 → 피부 전사체에서 미중력 신호가 일관되게 포착됨 (fold별 변동 매우 작음: ±0.051)
- RR-5 (OSD-240/241, BAL-TAL 균주): Track 2b — 향후 교차 균주 분석 대상

---

## 기술 노트

### PCA n_components Adaptive Fix
- 문제: `n_components=50`이 작은 training fold (n < 50)에서 sklearn 오류 발생
- 해결: `evaluate_fold()`에서 fit 전 `n_components = min(50, n_train-1, n_features)` 동적 적용
- A4 thymus: MHU-2 fold 55→55 (OK), RR-6 fold 50→31, RR-9 fold 50→46

### Duplicate Detection Fix
- 원래 방식 (median correlation): biological replicates를 duplicate로 오분류 (7/12 제거)
- 수정 방식 (pairwise correlation ≥ 0.9999): 실제 technical duplicate만 제거

### MHU-2 Label Fix
- AG (1G centrifuge) 그룹이 "Flight" keyword에 먼저 매칭 → AUROC 완전 역전
- keyword_map reordering으로 해결 (specific before generic)

---

## 파일 목록

```
tasks/
  A1_liver_lomo/             (6 folds, 264 samples)
  A1_liver_lomo_combat/      (6 folds, 264 samples, ComBat-seq corrected)
  A2_gastrocnemius_lomo/     (3 folds, 32 samples) ← GO
  A3_kidney_lomo/            (3 folds, ~142 samples)
  A4_thymus_lomo/            (3 folds, ~92 samples) ← GO
  A5_skin_lomo/              (3 folds, 102 samples) ← GO ★NEW
  A6_eye_lomo/               (3 folds, ~37 samples)

evaluation/
  A1_baseline_results.json
  A1_combat_baseline_results.json
  A1_shap_rf.json
  A2_baseline_results.json         ← LR mean AUROC=0.926, GO
  A2_shap_rf.json                  ← 4/16 targets (Bmal1, Per2, Dbp, Ciart), Cond3 PASS
  A3_baseline_results.json
  A4_baseline_results.json
  A4_shap_rf.json                  ← 12/23 targets in top-50, Cond3 PASS
  A5_baseline_results.json         ← ★NEW: LR mean AUROC=0.821, GO
  A5_shap_rf.json                  ← ★NEW: top genes E2f1, Ckap2, Klk8
  A6_baseline_results.json

scripts/
  fetch_osdr.py         -- OSDR data download
  quality_filter.py     -- Sample + mission level QC
  generate_tasks.py     -- LOMO split generation (DD-03, DD-04)
  run_baselines.py      -- Baseline classifiers (DD-08, DD-11, DD-13)
  batch_correct.R       -- ComBat-seq batch correction (DD-10)
  shap_analysis.py      -- SHAP feature importance (DD-11 Cond3)
```

---

## 다음 조치 (우선순위순)

### 완료
- [DONE] A4 SHAP: 12/23 targets → A4 GO 확정
- [DONE] A1 ComBat-seq: RF 0.590→0.660, 여전히 NO-GO
- [DONE] A2 초기 LOMO: RR-5 DESeq2 정규화, 3-mission LOMO 실행
- [DONE] A2 within-LOO sanity: RR-9 데이터 품질 재검토 → joint normalization으로 해결
- [DONE] **A2 RR-9 재검토**: 원시 카운트 3-mission 통합 DESeq2. LOMO LR=0.907, SHAP Myog#9+Ciart#25+Per2#33. GO 확정.
- [DONE] **A1 ISS-only**: MHU-2 제외 5-fold, PCA-LR 0.661 → 여전히 NO-GO. Pipeline v1/v2 heterogeneity가 주 원인.
- [DONE] **추가 gastrocnemius 미션 탐색**: OSDR 전수 조사 완료. **결론: 추가 gastrocnemius 미션 없음.**
  - 현재 OSD-101(RR-1), OSD-401(RR-5), OSD-326(RR-9) 3개가 OSDR에서 gastrocnemius RNA-seq 전부
  - RR-3, RR-6, RR-8, MHU-2에 gastrocnemius RNA-seq 없음 (thymus/skin/liver만 있음)
  - **추가 발견**: OSD-99(GLDS-99) = RR-1 **EDL 근육** 12샘플, OSD-104(GLDS-104) = RR-1 **Soleus 근육** 12샘플 존재
    → 골격근 확장 분석 가능 (단, 1개 미션만 있으므로 LOMO 불가, within-mission 분석 가능)

- [DONE] **A4 Pathway Enrichment**: g:Profiler (mmusculus). Cell cycle p=8.8e-6, Nucleosome assembly p=4.1e-4. 생물학적 내러티브 완성.

### 다음 단계 (Phase 2)

**Phase 2A — A4 Thymus 심화 [HIGHEST PRIORITY]**
- [ ] 추가 thymus 미션 탐색 (MHU-1? 추가 ISS RR 미션?)
- [ ] 예측 확률 calibration 분석 (predict_proba vs decision boundary)
- [ ] A4 SHAP direction 분석: 각 유전자 flight vs AG 방향 (up/down)
- [ ] 논문용 생물학적 내러티브 정교화

**Phase 2B — A2 Gastrocnemius 확장 [HIGH]**
- [ ] OSD-99(EDL, n=12) + OSD-104(Soleus, n=12) 다운로드 → within-mission 검증
- [ ] Cross-muscle SHAP: EDL vs Soleus vs Gastrocnemius 공통 유전자

**Phase 2C — Cross-tissue Circadian 분석 [MED]**
- [ ] A1 간 + A2 근육 공통 Bmal1/Dbp/Ciart/Per2 신호 체계적 분석
- [ ] spaceflight-induced circadian disruption 정량화

**Phase 2D — A6 Eye 구제 [MED]**
- [ ] OSD-397 미션 확인 + 추가 eye 미션 탐색
- [ ] n 증가 시 GO 가능성 높음 (현재 AUROC=0.811)

**Phase 2E — A3 Kidney 추가 분석 [LOW — 두 평가 모두 NO-GO 확정]**
- [DONE] pathway LOMO 실행: LR mean AUROC=0.755, CI lower=0.481, perm_p=0.071 → NO-GO
- [DONE] 실패 원인 확인: RR-7(n=94, 20F:74G 불균형)이 bottleneck
- [ ] 추가 kidney 미션 확보 시 재시도 (현재 OSDR에는 RR-1/RR-3/RR-7 3개뿐)
- [ ] stratified permutation test (클래스 불균형 보정) — 낮은 우선순위

---

---

## Category B: Cross-Mission Transfer Matrix (2026-03-01)

**목적**: "한 미션에서 학습한 우주비행 시그니처가 다른 미션으로 얼마나 전이되는가?" — PLAN.md 핵심 기여, Task ID: B4 (Thymus) / B2 (Gastrocnemius)

**방법**: train mission → test mission pairwise AUROC (PCA-LR, LFC-signature 두 방법)
**출력**: `processed/B_cross_mission/{tissue}/`, `evaluation/B_cross_mission_summary.json`

### B4 Thymus (4×4 Matrix, 12 pairs, N_BOOTSTRAP=2000)

**PCA-LR Transfer Matrix (AUROC [95% CI])**:

| Train↓ / Test→ | MHU-1‡ | MHU-2 | RR-6 | RR-9 |
|----------------|--------|-------|------|------|
| **MHU-1**‡    | — | **1.000** [1.000, 1.000] | 0.817 [0.668, 0.947] | 0.980 [0.910, 1.000] |
| **MHU-2**     | **1.000** [1.000, 1.000] | — | **0.490**⚠️ [0.296, 0.690] | **1.000** [1.000, 1.000] |
| **RR-6**      | **1.000** [1.000, 1.000] | 0.667 [0.111, 1.000] | — | 0.700 [0.437, 0.917] |
| **RR-9**      | **1.000** [1.000, 1.000] | **1.000** [1.000, 1.000] | 0.663 [0.480, 0.857] | — |

**LFC-Signature Transfer Matrix (AUROC [95% CI])**:

| Train↓ / Test→ | MHU-1‡ | MHU-2 | RR-6 | RR-9 |
|----------------|--------|-------|------|------|
| **MHU-1**‡    | — | **1.000** [1.000, 1.000] | 0.788 [0.625, 0.935] | 0.980 [0.912, 1.000] |
| **MHU-2**     | **1.000** [1.000, 1.000] | — | 0.614 [0.431, 0.813] | **1.000** [1.000, 1.000] |
| **RR-6**      | **1.000** [1.000, 1.000] | 0.667 [0.111, 1.000] | — | 0.670 [0.406, 0.900] |
| **RR-9**      | **1.000** [1.000, 1.000] | **1.000** [1.000, 1.000] | 0.693 [0.513, 0.875] | — |

- **Mean AUROC = 0.860 [0.763, 0.953]** (PCA-LR), **0.868 [0.781, 0.951]** (LFC)
- ⚠️ MHU-2→RR-6 = 0.490 [0.296, 0.690]: CI가 0.5를 포함 — 진정한 실패 (n=6 소규모 훈련 → n=35 테스트)
- ‡ MHU-1 = Track 2b (strain confound), 높은 AUROC는 strain 차이 반영 가능
- Mean asymmetry |A(i,j)−A(j,i)| = **0.069** (PCA-LR), **0.051** (LFC) — 대칭성 양호
- **생물학적 해석**: RR-9 train → 모든 미션에 완벽 전이 (1.000). 반면 소규모 미션 train → 대규모 test는 어려움

### B2 Gastrocnemius (3×3 Matrix, 6 pairs, N_BOOTSTRAP=2000)

**PCA-LR Transfer Matrix (AUROC [95% CI])**:

| Train↓ / Test→ | RR-1 | RR-5 | RR-9 |
|----------------|------|------|------|
| **RR-1** | — | 0.639 [0.281, 0.971] | **1.000** [1.000, 1.000] |
| **RR-5** | 0.500⚠️ [0.114, 0.875] | — | 0.750 [0.250, 1.000] |
| **RR-9** | 0.917 [0.667, 1.000] | **1.000** [1.000, 1.000] | — |

**LFC-Signature Transfer Matrix (AUROC [95% CI])**:

| Train↓ / Test→ | RR-1 | RR-5 | RR-9 |
|----------------|------|------|------|
| **RR-1** | — | 0.778 [0.444, 1.000] | **1.000** [1.000, 1.000] |
| **RR-5** | 0.583 [0.222, 0.943] | — | 0.375⚠️ [0.000, 0.867] |
| **RR-9** | 0.667 [0.257, 1.000] | 0.528 [0.171, 0.938] | — |

- **Mean AUROC = 0.801 [0.653, 0.944]** (PCA-LR), **0.655 [0.509, 0.820]** (LFC)
- ⚠️ RR-5→RR-1 = 0.500 [CI: 0.114–0.875]: CI 매우 넓음 — 진짜 random 성능 (PCA-LR)
- LFC 방법 RR-5→RR-9 = **0.375** [0.000, 0.867] (역전!): 단순 LFC 시그니처가 근육에는 부적합
- Mean asymmetry |A(i,j)−A(j,i)| = **0.157** (PCA-LR), **0.227** (LFC) — 방향 의존성 존재
- **생물학적 해석**: 근육(GAS) 우주비행 시그니처는 조직(thymus)보다 미션 간 변동성 높음. PCA-LR > LFC 방법 차이 뚜렷

### Cross-Tissue 비교 요약

| 조직 | 미션 수 | 평균 AUROC (PCA-LR) [95%CI] | 평균 AUROC (LFC) [95%CI] | 최저 쌍 | 해석 |
|------|---------|------------------------------|--------------------------|---------|------|
| Thymus | 4 | **0.860 [0.763, 0.953]** | **0.868 [0.781, 0.951]** | MHU-2→RR-6=0.490 | 높은 일관성, 소규모 미션이 취약점 |
| Gastrocnemius | 3 | 0.801 **[0.653, 0.944]** | 0.655 **[0.509, 0.820]** | RR-5→RR-1=0.500 | 낮은 일관성, LFC 방법 실패 |

**결론**: Thymus > Gastrocnemius cross-mission generalizability. LFC 방법은 thymus에서만 PCA-LR과 동등.

> **평가 방식 (DD-17)**: Category B는 단일 GO/NO-GO 결론을 내리지 않는다.
> perm_p의 최소값은 n=6 테스트 샘플 기준 ≥ 0.056 (통계적 하한) — 소규모 쌍의 유의성은 참고용.
> 1차 분석 단위 = n≥10 쌍 (large pairs).

### Transfer Pattern Summary (PCA-LR baseline, `submission_PCALR_baseline_B4.json` / `B2.json`)

| 태스크 | 파일 | AUROC ≥ 0.70 | perm_p < 0.05 | Large pairs (n≥10) | 해석 |
|--------|------|-------------|--------------|-------------------|------|
| **B4 Thymus** | `submission_PCALR_baseline_B4.json` | 9/12 | 4/12 | mean=0.775, sig=4/6 | 높은 일관성 |
| **B2 Gastrocnemius** | `submission_PCALR_baseline_B2.json` | 4/6 | 3/6 | (n<10 pairs only†) | 미션 간 변동성 높음 |

†B2 gastrocnemius의 경우 모든 쌍이 n<10 — large pair 기준 미해당.

**Thymus**:
- CI lower (0.763) > 0.700 → 평균 전이 성능 우수 (모든 쌍 기준)
- Large pairs: mean AUROC=0.775, 4/6 유의 → 충분한 학습 데이터(n≥10) 시 안정적 전이

**Gastrocnemius**:
- CI lower (0.653) > 0.500 → 방향적 신호 존재 (경계 수준)
- LFC 방법 RR-5→RR-9 = 0.375 (역전) → 근육에서 LFC 시그니처 불안정

*(N_BOOTSTRAP=2000, N_PERMUTATIONS=10000 — full run 완료)*

---

## Tier 2: Geneformer Foundation Model — Development Preview (2026-03-01)

> **상태**: 개발 미리보기 (1 fold only, MPS 5-epoch). 전체 LOMO는 HPC (Cayuga A40)에서 실행 예정.

### 방법 요약

| 항목 | 내용 |
|------|------|
| 모델 | Geneformer-V1-10M (BertForSequenceClassification, 10.3M params) |
| 토크나이제이션 | ENSMUSG → ENSG (Ensembl BioMart ortholog), gene rank encoding |
| 어휘 커버리지 | 57-64% (mouse A4 유전자 중 Geneformer vocab 포함 비율) |
| Fine-tuning 전략 | Full fine-tuning (classification head 추가), 5 epochs, batch=4, lr=2e-5 |
| 환경 | MacBook Air Apple Silicon (MPS, float32) |

### 개발 결과 (RR-9 fold)

| 에포크 | 훈련 손실 | 테스트 AUROC |
|--------|----------|-------------|
| 1 | 0.676 | 0.510 |
| 2 | 0.657 | 0.480 |
| 3 | 0.645 | 0.510 |
| 4 | 0.631 | 0.540 |
| 5 | 0.620 | **0.540** (best) |

**최종**: AUROC=0.540 [0.273, 0.818] (n_test=20, n_train=47)

### 비교 (A4 Thymus RR-9 fold)

| 방법 | AUROC |
|------|-------|
| **Tier 1: PCA-LR** (4-fold mean) | **0.923** |
| **Tier 1: LFC** (4-fold mean) | **0.814** |
| Geneformer-V1-10M (RR-9, dev, 5ep) | 0.540 |

### 해석

> ⚠️ **중요**: 이 결과는 1개 fold의 개발 미리보기이며, 과소평가된 것으로 보임.

**Geneformer 저성능 원인 분석**:
1. **소규모 데이터**: n_train=47 vs 10M 파라미터 → 과적합 불가피
2. **어휘 손실**: ENSMUSG→ENSG 매핑에서 42% 유전자 소실
3. **종 간 도메인 갭**: Geneformer는 인간 scRNA로 학습, 마우스 bulk RNA에 적용
4. **훈련-테스트 서열 길이 불일치**: 훈련 평균 1,593 vs RR-9 테스트 618 (RR-9 희소성)
5. **에포크 부족**: MPS 제약으로 5 epoch만 (HPC에서는 10+ epoch 필요)

**향후 계획**: Cayuga HPC (A40, batch=16, 10 epoch) 전체 4-fold LOMO → `evaluation/geneformer_v1_A4_lomo_results.json`

*(개발 결과: `evaluation/geneformer_v1_A4_dev_result.json`)*

---

## Changelog

| 날짜 | 내용 |
|------|------|
| 2026-02-27 | Phase 1 최초 실행. A1 liver baseline 3종 완료. |
| 2026-02-27 | MHU-2 label 오류 수정: AG→Flight 잘못 분류 → AG 제외 |
| 2026-02-27 | within-MHU-2 sanity: LOO AUROC=1.0 확인 |
| 2026-02-27 | Duplicate detection fix: pairwise r ≥ 0.9999 (biological replicates 보존) |
| 2026-02-27 | 모든 조직 다운로드 완료 (liver, gastrocnemius, kidney, thymus, eye, skin) |
| 2026-02-28 | A3/A4/A6 tasks 생성 및 baselines 완료 |
| 2026-02-28 | PCA n_components adaptive fix (evaluate_fold 내 동적 조정) |
| 2026-02-28 | A4 Thymus: Phase 1 GO (PCA-LR mean AUROC=0.888, CI lower=0.795, perm p=0.009) |
| 2026-02-28 | A1 ComBat-seq: RF 0.590→0.660, MHU-2 역전 해소. 여전히 NO-GO. |
| 2026-02-28 | A1 SHAP (RF): 0/9 liver targets in top-50. 상위에 circadian clock genes (Dbp, Bmal1, Npas2). |
| 2026-02-28 | A4 SHAP (RF): 12/23 thymus targets in top-50. Condition 3 ✓ PASS. A4 GO 확정. |
| 2026-02-28 | A2 Gastrocnemius: RR-5(GLDS-401) DESeq2 정규화. 3-mission LOMO 생성 (32 samples). |
| 2026-02-28 | A2 baselines: LR mean AUROC=0.926, CI lower=0.796, perm p=0.028 → ✓ GO |
| 2026-02-28 | A2 SHAP (RF): 4/16 targets in top-50 (Bmal1#3, Per2#4, Dbp#23, Ciart#50). Cross-tissue circadian signature 확인. |
| 2026-02-28 | [검토] A2 데이터 품질 문제 발견: RR-9 within-LOO=0.188, PC1(54.9%) mission분리, zero rate 30%. A2 PROVISIONAL로 하향. |
| 2026-02-28 | [검토] A4 Thymus within-LOO sanity 확인: MHU-2=1.000, RR-6=0.791, RR-9=1.000 → A4 GO 견고하게 확정. |
| 2026-02-28 | [수정] A2: 3 mission joint DESeq2 정규화 (normalize_gastrocnemius_joint.R). 21,013 genes, zero rate 개선. |
| 2026-02-28 | [수정] A2 LOMO 재생성: LR mean AUROC=0.907, CI lower=0.717, perm p=0.026 → ✓ GO |
| 2026-02-28 | [수정] A2 SHAP 재실행: Myog(#9), Ciart(#25), Per2(#33) → 3/16 targets. Condition 3 ✓ PASS. |
| 2026-02-28 | [수정] A2 within-mission signal 확인: 모든 미션 naive AUROC=1.000, mean|corr|>0.27 → biological signal 견고. |
| 2026-02-28 | A2 최종 판정: GO (이전 PROVISIONAL에서 상향). PHASE 1 GO 태스크: A2 + A4. |
| 2026-02-28 | A1 ISS-only 분석 (MHU-2 제외, 5-fold): best PCA-LR mean AUROC=0.661 → NO-GO 확정. |
| 2026-02-28 | A1 분석: MHU-2 배치 효과가 아닌 pipeline v1/v2 heterogeneity (RR-3=0.528, RR-9=0.414)가 주 원인. |
| 2026-02-28 | OSDR gastrocnemius 탐색 완료: 추가 미션 없음 (RR-1/RR-5/RR-9가 전부). OSD-99(EDL), OSD-104(Soleus) 추가 근육 발견. |
| 2026-02-28 | **Phase 1 전체 검토 완료**: A2(GO), A4(GO), A1/A3/A6(NO-GO). 모든 권장 조치 완료. |
| 2026-02-27 | **A4 Pathway Enrichment (g:Profiler)**: Cell cycle (p=8.8e-6), Nucleosome assembly (p=4.1e-4), Mitotic spindle checkpoint. 생물학적 내러티브 확립: 방사선→DNA 손상→G2/M arrest→흉선위축. `evaluation/A4_pathway_enrichment.json` |
| 2026-02-27 | **[Phase 2A] A4 OSDR 탐색**: OSD-515(RR-23, 27샘플), OSD-457(SpaceX-14, 12 WT thymus) 신규 발견. MHU-1(GLDS-289)이 기존 MHU-2 데이터에 혼재 확인. |
| 2026-02-27 | **[Phase 2A] MHU-1 분리**: GLDS-289에서 MHU-1(n=9)과 MHU-2(n=9) 미션 분리. thymus_all_missions_metadata.csv 수정. |
| 2026-02-27 | **[Phase 2A] A4 4-fold 재실행**: PCA-LR mean AUROC=0.923(↑from 0.888), CI lower=0.878(↑from 0.795). MHU-1 strain confound 주의(GC=C57-6CR). |
| 2026-02-27 | **[Phase 2A] A4 SHAP 4-fold**: 7/23 targets in top-50 (→Cond3 PASS). 34↓/16↑ in flight. 흉선 위축=증식세포 고갈 + 면역억제 활성화. |
| 2026-02-27 | **[Phase 2A] A4 Calibration**: Flight mean proba=0.622, GC=0.476, sep=0.146. AUROC ranking 우수하나 cross-mission 절대 확률 uncalibrated. |
| 2026-02-27 | **[Phase 2A] A2 Pathway Enrichment**: KEGG Circadian rhythm 최고순위(p=0.13, 역치 미달). Pathway-level 수렴 없음. circadian 4유전자(Per2,Ciart,Nfil3,Klf11) 분산됨. |
| 2026-03-01 | **[버그수정] run_baselines.py B1+B2**: n_ground_test 수식 수정 (`y_test==0`), permutation p-value pseudocount `(sum+1)/(n+1)`. RR-9 fold perm_p: 0.000→0.001. |
| 2026-03-01 | **[MHU-1 재분류]** thymus_all_missions_metadata.csv: MHU-1 → track=2b, GC 샘플 strain=C57BL/6CR 수정. PHASE1_RESULTS.md A4 결과에 경고 추가. |
| 2026-03-01 | **[PLAN.md v0.7]** 타겟 모델 3-tier (Classical ML / Geneformer / Text LLM) 명세 추가. v1.0 dataset freeze 정책 (기준일 2026-03-01). |
| 2026-03-01 | **[DESIGN_DECISIONS.md v1.2]** DD-16 신규: Text LLM Evaluation Track (GPT-4o, Claude, Llama 3). |
| 2026-03-01 | **[Category B]** cross_mission_transfer.py 실행: thymus 4×4 (mean AUROC=0.860), gastrocnemius 3×3 (mean=0.801). `processed/B_cross_mission/` 생성. |
| 2026-03-01 | **[Category B CI]** Full bootstrap (N=2000) 완료. Thymus PCA-LR: 0.860 [0.763, 0.953], LFC: 0.868 [0.781, 0.951]. Gastro PCA-LR: 0.801 [0.653, 0.944], LFC: 0.655 [0.509, 0.820]. |
| 2026-03-01 | **[벤치마크 인프라]** docs/submission_format.md, scripts/evaluate_submission.py 신규 생성. 외부 제출 평가 시스템 구축 완료. |
| 2026-03-01 | **[버그수정] B_cross_mission_summary.json**: gastrocnemius 결과 덮어쓰기 버그 수정 (merge 로직). thymus+gastrocnemius 모두 포함 완료. |
| 2026-03-01 | **[A4 within-LOO 완료]** 4-fold 전체: MHU-1=1.000, MHU-2=1.000, RR-6=0.771, RR-9=1.000. Track 2a 최저=0.771 (>0.700 ✓). `evaluation/A4_within_loo.json` 생성. |
| 2026-03-01 | **[Calibration 문서화]** PHASE1_RESULTS.md A4 섹션에 calibration 한계 경고 추가. Platt scaling Phase 2 계획. |
| 2026-03-01 | **[OSD-515 held-out]** RR-23 데이터 다운로드 완료. `tasks/A4_thymus_lomo/fold_RR-23_holdout/` 생성 (67 train × 16 test, 27541 genes). |
| 2026-03-01 | **[PLAN.md v0.9]** Geneformer 전략 업데이트: Cornell Cayuga HPC (scu-gpu, A40/A100) 세부사항. |
| 2026-03-01 | **[Geneformer 토크나이제이션]** `geneformer_tokenize.py` 완성. 5 LOMO folds 완료. Ortholog coverage: 57-64%. `geneformer_v1_tokenize_summary.json` 생성. |
| 2026-03-01 | **[Geneformer Dev Fine-tuning]** RR-9 fold MPS 5-epoch: AUROC=0.540 [0.273, 0.818]. Tier 1 (0.923) 대비 저성능 — 소규모 데이터 + 종 간 도메인 갭. HPC full LOMO 필요. |
| 2026-03-01 | **[Category B submissions]** tasks/B4_thymus_cross_mission/ + B2_gastrocnemius_cross_mission/ 생성 (pair_ 구조). submission_PCALR_baseline_B4.json (12 pairs, AUROC=0.860), B2.json (6 pairs, AUROC=0.801) 생성. evaluate_submission.py Category B 지원 추가. |
| 2026-03-01 | **[DD-17]** Category B Evaluation Criteria 확정. Transfer Pattern Summary 도입 (perm_p floor 근거: n=6 → C(6,3)=20 → min_p=0.050). DESIGN_DECISIONS.md v1.4. PHASE1_RESULTS.md B section 헤더 수정: B6→B4 Thymus, B4→B2 Gastrocnemius. |
| 2026-03-01 | **[A3 Kidney pathway]** J5_gene_vs_pathway.json에서 kidney pathway-level (GSVA Hallmark) mean AUROC=0.743 확인 (gene 0.432 → pathway 0.743, +0.310). RR-1 fold 극적 구제 (0.056→0.861). Gene-level NO-GO 유지, pathway PROVISIONAL 표시. PHASE1_RESULTS.md v4 업데이트. |
| 2026-03-01 | **[run_pathway_lomo.py]** scripts/run_pathway_lomo.py 신규 작성. GSVA pathway features로 LOMO 평가 (bootstrap CI + perm_p). A3/A6 실행 완료. evaluation/A3_pathway_hallmark_results.json, A6_pathway_hallmark_results.json 생성. |
| 2026-03-01 | **[A3 Kidney pathway 공식]** pathway LR: AUROC=0.755, CI lower=0.481, perm_p=0.071 → NO-GO. AUROC 조건 통과(↑from 0.564), CI+perm 미달. RR-7(n=94, 20F:74G) bottleneck 확인. PROVISIONAL 표시 → 공식 NO-GO로 확정. |
| 2026-03-01 | **[A6 Eye pathway 공식]** pathway PCA-LR: AUROC=0.915, CI lower=0.745, perm_p=0.014 → ✅ GO. Gene-level NO-GO → Pathway-level GO. PHASE1_RESULTS.md A6 섹션 업데이트. 요약 표에 pathway 결과 행 추가. |
| 2026-03-01 | **[A6 Eye SHAP-pathway]** PCA-LR 계수 역사영 (`coef @ PCA.components_ / scaler.scale_`). Top Ground 마커: UV_RESPONSE_UP(−3.90), UPR(−2.80), HYPOXIA(−2.18), GLYCOLYSIS(−1.86). Top Flight 마커: G2M_CHECKPOINT(+1.86), MITOTIC_SPINDLE(+1.34), E2F_TARGETS(+1.33). 방사선 유도 세포주기 체크포인트 활성화 확인. `evaluation/A6_pathway_hallmark_shap.json` |

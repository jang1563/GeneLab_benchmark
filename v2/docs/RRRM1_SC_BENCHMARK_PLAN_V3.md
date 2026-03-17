# RRRM-1 scRNA-seq Benchmark Plan (v3)

Date: 2026-03-12 (updated 2026-03-13)
Author: GeneLab Benchmark project

---

## 1. 왜 scRNA-seq 분석을 시작했는가

### 1.1 PLAN.md §11 Category F 근거

`PLAN.md` v1.4, §11 "Future Extensions":

> **Category F (Single-Cell & Spatial)**
> - F2: RRRM-1 다조직 scRNA-seq (GLDS-746~762) — 세포 조성 변화 (검증 후)

v1.0 벤치마크 (6 tissues, ~450 samples, bulk RNA-seq)는 **조직 수준** spaceflight 탐지를 수행한다.
v2.0 Category F의 목적은 이 신호를 **세포 유형 수준으로 분해**하는 것이다:

- "어느 세포 유형이 bulk spaceflight 신호를 주도하는가?"
- "세포 조성 자체가 우주비행으로 변하는가?"
- "v1.0에서 발견된 pathway 신호가 특정 세포 유형에서만 나타나는가?"

### 1.2 v1.0과의 연결

```
v1.0 (bulk):  조직 전체 → FLT vs. GC  → AUROC (Category A)
                                            pathway (Category J/GSVA)
                ↓ Cell-type 분해
v2.0 (scRNA): 조직 내 세포 유형별 → FLT vs. GC per cell type
                                     → 세포 조성 변화
                                     → pathway 신호 귀속
```

### 1.3 선택한 4개 조직이 v1.0과 동일한 이유

| scRNA OSD | Tissue  | v1.0 대응 Task |
|-----------|---------|----------------|
| OSD-918   | Blood   | (A6, eye 인접) — immune context |
| OSD-920   | Eye     | A6 Eye (2 missions, ⚠️ 미션 수 부족) 보완 |
| OSD-924   | Muscle  | A2 Gastrocnemius (bulk muscle) 비교 |
| OSD-934   | Skin    | A5 Skin + held-out RR-7 비교 |

---

## 2. 현재 파이프라인 상태 (2026-03-13 기준)

### 2.1 완료된 단계

```
STARsolo (GLDS-746~762, CellRanger2.2 barcode filter)
  ↓ Solo.out/GeneFull/filtered ✅
rrrm1_h5ad_convert.py
  ↓ per-OSD .h5ad ✅
rrrm1_initial_scanpy.py (QC: min_genes=200, pct_mt<25%)
  ↓ *_processed.h5ad ✅ (Mar 12)
rrrm1_broad_annotate.py (marker score → broad label)
  ↓ *_annotated.h5ad ✅ (Mar 12)
rrrm1_singlecell_hardening.py (doublet removal via scrublet)
  ↓ *_hardened.h5ad ✅ (Mar 13)
```

### 2.2 세포 수 및 구성 (hardened 기준)

| Tissue  | 총 세포   | 주요 세포 유형                                    | Doublet 제거 |
|---------|---------|------------------------------------------------|------------|
| Blood   | 4,377   | erythroid 83%, B cell 7%, T cell 5%             | 18 (0.4%)  |
| Eye     | 2,197   | retinal neuronal 52%, Muller glia 18%            | 9 (0.4%)   |
| Muscle  | 15,669  | macrophage/myeloid 38%, T/NK 25%, FAP 13%       | 349 (2.2%) |
| Skin    | 15,838  | immune myeloid 37%, basal keratinocyte 30%       | 38 (0.2%)  |

### 2.3 현재 실행 중인 Jobs (Cayuga)

- `bioreview_sft` (Job 2702484): RUNNING 23h+ on g0002 — Codex v3 작업
- `figures` (Job 2702039): DependencyNeverSatisfied ← **취소 검토 필요**
- `bioreview_inf/src` (Job 2702592, 2702810, 2703139): Dependency pending

---

## 3. 핵심 전제: FLT/GC 레이블

### 3.1 현재 갭

모든 downstream 분석의 전제는 **각 세포에 FLT/GC 조건 레이블**이 있어야 한다.
현재 h5ad의 메타데이터에 `condition` (FLT vs. GC) 컬럼이 있는지 확인 필요.

RRRM-1 (OSDR GLDS-746~762):
- 10x Chromium 3' v3, 멀티플렉싱 여부 확인 필요
- 각 SRX run = 1마리? 또는 pooled?
- Typical RRRM-1 설계: 6 FLT + 6 GC mice (n=6 각 조건)

### 3.2 확인 방법

```python
import scanpy as sc
adata = sc.read_h5ad("RRRM1_blood_hardened.h5ad")
print(adata.obs.columns.tolist())
print(adata.obs["condition"].value_counts())   # or "group", "treatment" 등
```

**Phase 0 (최우선)**: FLT/GC 레이블이 obs에 있는지, animal_id가 있는지 확인.
없으면 OSDR metadata에서 SRX → FLT/GC 매핑을 추가해야 함.

---

## 4. Benchmark Task 정의

RRRM-1은 **단일 미션**이므로 v1.0의 LOMO(Leave-One-Mission-Out)는 불가능.
대신 **LOAO(Leave-One-Animal-Out)** cross-validation 사용.

### Task F2-A: Cell-type Composition Shift

**질문**: 우주비행이 조직 내 세포 유형 비율을 바꾸는가?

```
입력: per-animal 세포 유형 비율 (FLT k동물, GC k동물)
출력: 각 세포 유형의 FLT/GC proportion 차이 + 통계
평가 지표: Wilcoxon rank-sum p-value, effect size (Cohen's d)
다중 검정 보정: BH FDR
```

| Tissue  | FLT 동물 | GC 동물 | 분석 가능 여부       |
|---------|--------|--------|-----------------|
| Blood   | 확인 필요 | 확인 필요 | OSDR metadata 확인 |
| Eye     | 확인 필요 | 확인 필요 | OSDR metadata 확인 |
| Muscle  | 확인 필요 | 확인 필요 | OSDR metadata 확인 |
| Skin    | 확인 필요 | 확인 필요 | OSDR metadata 확인 |

**연결점**: v1.0 D3 (mission prediction) — 세포 조성 자체가 배치 효과 원천인지 확인 가능.

---

### Task F2-B: Cell-type-specific Pathway Activity (pseudo-bulk fGSEA)

**질문**: 어느 세포 유형이 bulk spaceflight pathway 신호를 주도하는가?

```
파이프라인:
1. 각 tissue별, 세포 유형별, 동물별 pseudo-bulk 합산 (sum of raw counts)
2. DESeq2 FLT vs. GC (per cell type, per tissue)
3. fGSEA on Hallmark DB (NES, padj)
4. v1.0 bulk NES와 비교 (Spearman r per pathway)
```

**기준 비교**:
- Muscle: v1.0에서 MuRF1/MAFbx (근위축 E3 리가제) 탐지 여부 → 어느 세포 유형에서?
- Skin: v1.0 held-out RR-7 AUROC=0.885 → 어느 세포 유형이 신호 원천?
- Blood: I4 PBMC (F1, snRNA-seq)와의 비교 — 같은 immune 세포 유형에서 일치하는가?

**연결점**: v2.0 E1/E2/E3 (cross-species) — 마우스 scRNA muscle vs. human I4 PBMC.

---

### Task F2-C: Cell-level Spaceflight Classifier (핵심 benchmark task)

**질문**: 세포 수준에서 FLT vs. GC를 얼마나 잘 구분할 수 있는가? 세포 유형마다 다른가?

```
방법: PCA-LR (v1.0 baseline과 동일)
입력: single-cell normalized expression (log1p)
검증: Leave-One-Animal-Out (LOAO) cross-validation
평가: AUROC, Bootstrap 95% CI (n=1000), permutation p-value
단위: tissue × cell type
```

**설계 세부사항**:
```python
# per-tissue, per-cell-type LOAO loop
for tissue in [blood, eye, muscle, skin]:
    for cell_type in adata.obs["broad_celltype"].unique():
        subset = adata[adata.obs["broad_celltype"] == cell_type]
        for animal_out in subset.obs["animal_id"].unique():
            train = subset[subset.obs["animal_id"] != animal_out]
            test  = subset[subset.obs["animal_id"] == animal_out]
            # PCA(50) → LR(ElasticNet) → predict_proba
        # aggregate LOAO predictions → AUROC
```

**최소 세포 수 기준**: 각 조건(FLT/GC) × 각 세포 유형당 ≥ 20 cells (power threshold).
혈액 erythroid는 83%로 충분하지만 NK (14세포)는 제외 가능성 있음.

**예상 산출물**:
```
F2-C_results.json:
  blood:
    erythroid: {auroc: ?, ci_lo: ?, ci_hi: ?, p_perm: ?}
    b_cell: {auroc: ?, ...}
    ...
  muscle:
    macrophage_myeloid: {auroc: ?, ...}
    ...
```

**v1.0 비교**: Category A bulk AUROC vs. F2-C per-cell-type AUROC
→ "어느 세포 유형의 scRNA 신호가 bulk 신호보다 강한가?"

---

### Task F2-D: Cross-species Pathway Concordance (mouse RRRM-1 ↔ human I4)

**질문**: RRRM-1 마우스 scRNA pathway 신호가 I4 인간 PBMC fGSEA 신호와 일치하는가?

```
방법:
1. RRRM-1 blood monocyte pseudo-bulk fGSEA NES (FLT vs. GC)
2. I4 CD14+ monocyte fGSEA NES (flight vs. pre)
3. Spearman r on 50 Hallmark pathways
4. 비교: r(RRRM-1 monocyte, I4 monocyte) vs. r(RRRM-1 bulk, I4 cfRNA) = 0.352
```

**연결점**:
- E1 (mouse liver bulk vs. JAXA cfRNA r=0.352) → 세포 유형 매칭하면 더 높아지는가?
- F1 (I4 PBMC 10 cell types fGSEA) → 대응하는 마우스 세포 유형과 비교

---

## 5. 전체 실행 계획

### Phase 0: FLT/GC 레이블 확인 (즉시)

```bash
# Cayuga에서
python3 - << 'EOF'
import scanpy as sc
for tissue in ["blood", "eye", "muscle", "skin"]:
    h5ad = f"downstream_initial/hardening/objects/RRRM1_{tissue}_hardened.h5ad"
    adata = sc.read_h5ad(h5ad)
    print(f"\n=== {tissue} ===")
    print("obs columns:", adata.obs.columns.tolist())
    for col in ["condition", "group", "treatment", "sample", "animal_id"]:
        if col in adata.obs.columns:
            print(f"  {col}:", adata.obs[col].value_counts().to_dict())
EOF
```

**If labels missing**: OSDR GLDS-746~762 metadata에서 SRX → FLT/GC 매핑 추가 필요.

---

### Phase 1: F2-A Composition (1일)

**Script**: `v2/scripts/rrrm1_composition_analysis.py`

```python
# 개요
# 1. per-tissue: obs["animal_id"] × obs["broad_celltype"] → proportion matrix
# 2. FLT vs. GC Wilcoxon test per cell type
# 3. stacked bar (FLT vs. GC) per tissue figure
# 출력: F2A_composition_results.json, F2A_composition_figure.html
```

**HPC**: login node에서 실행 가능 (연산 경량)

---

### Phase 2: F2-B Pseudo-bulk fGSEA (2일)

**Script**: `v2/scripts/rrrm1_pseudobulk_fgsea.R`

```r
# 개요 (DESeq2 + fGSEA)
# 1. AnnData → pseudo-bulk (per animal, per cell type)
# 2. DESeq2: ~ condition (FLT vs. GC)
# 3. fGSEA: Hallmark DB (mouse gene sets)
# 4. NES matrix: pathways × cell types
# 출력: {tissue}_celltype_fgsea_hallmark.csv
```

**HPC**: sbatch, ~1-2h per tissue on scu-gpu

---

### Phase 3: F2-C Cell Classifier (1-2일)

**Script**: `v2/scripts/rrrm1_celltype_classifier.py`

```python
# 개요
# 1. per-tissue, per-cell-type LOAO PCA-LR
# 2. AUROC + bootstrap CI + permutation p
# 3. comparison table vs. v1.0 bulk AUROC
# 출력: F2C_classifier_results.json, F2C_auroc_comparison.html
```

**HPC**: GPU 불필요, sbatch --cpus-per-task=8 --mem=32G

---

### Phase 4: F2-D Cross-species (0.5일)

**Script**: `v2/scripts/rrrm1_crossspecies_concordance.py`

```python
# 개요
# 1. RRRM-1 blood monocyte/macrophage pseudo-bulk NES
# 2. I4 PBMC monocyte NES (기존 v2/processed/F1_scrna/i4_snrnaseq_celltype_fgsea.csv)
# 3. Spearman r: RRRM-1 mouse cell type ↔ I4 human cell type
# 4. 비교: E1 r=0.352 (bulk liver ↔ cfRNA)
# 출력: F2D_crossspecies_results.json
```

---

### Phase 5: Figure 통합 (1일)

새 Fig4 (scRNA-seq summary):
- Panel A: Cell-type composition FLT vs. GC (4 tissues stacked bar)
- Panel B: F2-B top cell types with significant pathway changes (heatmap)
- Panel C: F2-C AUROC per cell type vs. v1.0 bulk AUROC
- Panel D: F2-D cross-species concordance (RRRM-1 monocyte ↔ I4 monocyte)

---

## 6. 의존성 및 제약

### 6.1 단일 미션 제약

RRRM-1 = 1개 미션 → LOMO 불가 → LOAO 사용
- LOAO 유효성: 동물 간 독립성 가정 (cage effect 있을 수 있음)
- 결론에 "single-mission, LOAO validation" 명시 필요

### 6.2 세포 수 제약

| Tissue | 문제 세포 유형              | n cells | 대응                   |
|--------|------------------------|---------|----------------------|
| Blood  | NK cell                | 14      | F2-C에서 제외 (n<20)    |
| Blood  | monocyte_macrophage    | 25      | 경계선 — F2-B만          |
| Eye    | lens_crystallin        | 26      | F2-B/C 제외 고려         |

### 6.3 FLT/GC 레이블 (Phase 0 블로커)

FLT/GC 레이블 없이 Phase 1-4 전체 불가.
→ Phase 0을 즉시 실행하여 확인.

### 6.4 v1.0 baseline 비교를 위한 데이터

F2-C 비교용 v1.0 AUROC:
- A2 Gastrocnemius: bulk gene AUROC (기존 결과)
- A5 Skin: bulk gene AUROC = 0.885 (A5 held-out)
→ `evaluation/` 디렉토리에서 로드

---

## 7. 성공 기준

| Task | 성공 기준                                              |
|------|-----------------------------------------------------|
| F2-A | ≥1 tissue에서 유의한 FLT/GC 세포 유형 비율 차이 (BH p<0.05) |
| F2-B | ≥1 cell type에서 v1.0 bulk pathway 신호 재현 (r>0.3)   |
| F2-C | ≥1 cell type에서 AUROC > v1.0 bulk baseline           |
| F2-D | RRRM-1 monocyte ↔ I4 monocyte r > E1 bulk r=0.352   |

---

## 8. 즉시 해야 할 것

### Step 1: bioreview_sft job (2702484) 확인
23h+ running — 무엇을 하고 있는지 확인 필요.
```bash
sacct -j 2702484 --format=JobID,JobName,State,Elapsed
cat /athena/.../logs/bioreview_sft_2702484.out | tail -20
```

### Step 2: Phase 0 실행 (FLT/GC 레이블 확인)
```bash
# Cayuga에서 hardened h5ad obs 컬럼 확인
```

### Step 3: figures job (2702039) 취소 고려
DependencyNeverSatisfied → 원인 파악 또는 취소.

### Step 4: raw FASTQ cleanup
308 GB → 삭제 가능 조건: manifest frozen + h5ad readable
현재 조건 충족됨 (manifest 2026-03-12 confirmed).

---

## 9. 버전 관리

이 문서는 `v2/docs/RRRM1_SC_BENCHMARK_PLAN_V3.md`로 저장.
이전 계획 문서: `RRRM1_DOWNSTREAM_PLAN_2026-03-11.md`, `RRRM1_NEXT_STEPS_FOR_BENCHMARK_2026-03-12.md`

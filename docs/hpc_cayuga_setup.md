# Cayuga HPC — Geneformer LOMO Setup Guide

*Cornell CAC Cayuga cluster — Slurm v25.05.0*
*Contact: scu@med.cornell.edu (subject: "Cayuga")*

---

## 1. 클러스터 구성

| 노드 | GPU | VRAM | CPU | RAM | 파티션 |
|------|-----|------|-----|-----|--------|
| g0001 | 4× A100 PCIe | 80GB each | 128 | 1024GB | scu-gpu |
| g0002–g0003 | 4× A40 PCIe each | 48GB each | 128 | 1024GB | scu-gpu |
| g0004 | 4× H100 | — | 128 | 1032GB | **restricted** (별도 신청) |
| c0001–c0011 | — | — | 112 | 768GB | scu-cpu |
| c0012–c0023 | — | — | 128 | 512GB | scu-cpu |

**파티션 최대 시간**: 7일 (168h)
**GRES 문법**:
```
--gres=gpu:a40:1    # A40 1개
--gres=gpu:a100:1   # A100 1개
```

---

## 2. 스토리지

| 위치 | 경로 | 특성 |
|------|------|------|
| 홈 디렉토리 | `/home/fs01/<cwid>` | NFS, 백업 없음, 용량 소 |
| Athena scratch | `/athena/cayuga_XXXX/scratch/$USER/` | 병렬 파일시스템 7.6PB, **모든 작업 데이터 여기** |
| Lab symlink | `/athena/<labname>` → `/athena/cayuga_XXXX` | |

> **규칙**: 모든 conda env, 프로젝트 파일, 체크포인트는 **Athena scratch**에 저장.
> 홈 디렉토리는 설정파일(.bashrc 등)만 사용.

---

## 3. 원타임 환경 설정 (로그인 노드에서 실행)

```bash
# 1. Athena에 Miniconda 설치 ($HOME이 아님!)
mkdir -p /athena/cayuga_XXXX/scratch/$USER/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
     -O /tmp/miniconda.sh
bash /tmp/miniconda.sh -b -u -p /athena/cayuga_XXXX/scratch/$USER/miniconda3
source /athena/cayuga_XXXX/scratch/$USER/miniconda3/etc/profile.d/conda.sh

# 2. geneformer 환경 생성 (CUDA 12.4용 PyTorch)
conda create -n geneformer python=3.11 -y
conda activate geneformer
pip install torch --index-url https://download.pytorch.org/whl/cu124
pip install transformers datasets accelerate scikit-learn safetensors tqdm huggingface_hub

# 3. 프로젝트 복사
cp -r /path/to/GeneLab_benchmark /athena/cayuga_XXXX/scratch/$USER/

# 4. 사전 토큰화 (CPU, 로그인 노드 OK, ~5분)
cd /athena/cayuga_XXXX/scratch/$USER/GeneLab_benchmark
conda activate geneformer
python scripts/geneformer_tokenize.py \
    --task A4 \
    --task-dir A4_thymus_lomo \
    --model-version v1
# → tasks/A4_thymus_lomo/fold_*/geneformer_tokens/v1/ 에 저장됨

# Mouse-Geneformer 사용 시 (local dictionary 경로 명시)
python scripts/geneformer_tokenize.py \
    --task A4 \
    --task-dir A4_thymus_lomo \
    --model-version mouse_gf \
    --mouse-gf-base /athena/cayuga_XXXX/scratch/$USER/GeneLab_benchmark/models
```

---

## 4. 작업 제출

```bash
# hpc_submit_geneformer.sh 상단 변수 편집 먼저:
# ATHENA_SCRATCH="/athena/cayuga_XXXX/scratch/${USER}"  ← XXXX 수정

# LOMO 배열 작업 제출 (4 fold 동시 실행)
cd /athena/cayuga_XXXX/scratch/$USER/GeneLab_benchmark
sbatch scripts/hpc_submit_geneformer.sh

# 상태 확인
squeue -u $USER
sacct -j <JOBID> --format=JobID,State,ExitCode,Elapsed,MaxRSS

# 인터랙티브 디버그 (A40 1개, 1시간)
srun -p scu-gpu --gres=gpu:a40:1 --mem=32G --time=01:00:00 --pty bash
```

---

## 5. 설정 요약 (hpc_submit_geneformer.sh)

| 파라미터 | 값 | 이유 |
|----------|-----|------|
| `--partition` | `scu-gpu` | GPU 파티션 |
| `--gres` | `gpu:a40:1` | A40 48GB (충분) |
| `--mem` | `32G` | 안전 마진 (노드: 1024GB) |
| `--cpus-per-task` | `4` | DataLoader num_workers=4 |
| `--time` | `04:00:00` | fold당 ~2h 예상 |
| `--array` | `0-3` | 4 LOMO fold (RR-23 제외) |
| `EPOCHS` | `10` | 개발 5 → 프로덕션 10 |
| `BATCH_SIZE` | `16` | A40 48GB에서 안정적 |
| `FREEZE_LAYERS` | `4` | n~50: 하위 4층 고정, 상위 2층 + head 학습 |

---

## 6. 결과 수집

각 fold 완료 후:
```
tasks/A4_thymus_lomo/fold_*/geneformer_tokens/v1/
├── checkpoints/best_model/   ← 최적 체크포인트
└── finetune_result.json      ← fold별 AUROC, CI, 학습 곡선
```

전체 LOMO 요약:
```
evaluation/geneformer_v1_A4_lomo_results.json
```

로컬에서 결과 동기화 (HPC → Mac):
```bash
rsync -avz <cwid>@cayuga.med.cornell.edu:/athena/cayuga_XXXX/scratch/$USER/GeneLab_benchmark/evaluation/ \
    /path/to/local/GeneLab_benchmark/evaluation/
```

---

## 7. 다음 단계 (HPC 완료 후)

1. `evaluation/geneformer_v1_A4_lomo_results.json` 로컬 동기화
2. `PHASE2_RESULTS.md` 작성 — 4-fold LOMO mean AUROC vs Tier 1 비교
3. Layer freezing ablation 결과 정리 (freeze=0 vs 4 vs 6)

# cell-eval Integration Note — GeneLab Benchmark

**Created:** 2026-05-04
**Trigger:** Mehran Karimzadeh (Xaira) flagged cell-eval (Arc Institute, May 2025) as the de facto evaluation standard for perturbation prediction models. While GeneLab Benchmark targets bulk transcriptomics + cross-mission classification (different scope), selected cell-eval design principles are worth porting to v8+.

---

## Context: Why This Note Exists

GeneLab Benchmark v1–v7 evaluates models on **mission-level Flight vs. Ground classification** using LOMO CV. cell-eval evaluates models on **per-perturbation prediction quality** using PDS + DES + MAE.

The two frameworks have **non-overlapping primary use cases**:

| | GeneLab Benchmark | cell-eval (Arc) |
|---|-------------------|-----------------|
| Data modality | Bulk RNA-seq (mouse OSDR) | Single-cell Perturb-seq |
| Task | Cross-mission classification (LOMO) | Perturbation effect prediction |
| Primary metric | Macro-F1, balanced accuracy | PDS, DES (overlap_at_N), MAE |
| Generalization unit | Mission (RR-1, RR-3, ...) | Cell type / perturbation |
| Model targets | Geneformer, scGPT, UCE, scFoundation, scPRINT2, text-LLMs | GEARS, CPA, scVI, X-Cell, STATE |

**However:** Some cell-eval design principles are transferable to GeneLab Benchmark, especially for the v8+ scope that includes scRNA-seq RRRM-1/RRRM-2 datasets.

---

## What cell-eval Is

- Arc Institute Python package: https://github.com/ArcInstitute/cell-eval
- Created: Yusuf Roohani, May 2025
- Companion to the Virtual Cell Challenge (VCC, *Cell* 2025)
- Used by Xaira's X-Cell (March 2026), 300+ VCC teams

**Citation:**
Roohani YH et al. "Virtual Cell Challenge: Toward a Turing test for the virtual cell." *Cell* 2025; 188(13):3370–3374.

---

## Three Transferable Design Principles

### 1. Discrimination over accuracy (PDS principle)

**cell-eval's PDS** asks: "Can the model tell perturbation A apart from perturbation B?" — not just "Is prediction A accurate in isolation?"

**Translation to GeneLab Benchmark:** Add a **mission discrimination** metric — for each held-out mission, rank model embeddings/predictions by similarity to all training missions. Score = where the held-out mission's nearest neighbor falls (1.0 if always correctly identified, 0.5 if random).

Currently GeneLab measures **classification accuracy** (Flight vs. Ground per tissue). A discrimination metric would test whether the model captures **mission-specific signatures**, not just the binary label. This is a stronger generalization claim.

**Implementation pattern:**
```python
# For each model embedding of mission M:
#   1. Compute distance to centroids of all other missions
#   2. Rank: where does mission M's own centroid fall in the sorted list?
#   3. Normalize to [0, 1]; aggregate across missions
```

### 2. Multi-metric profiles (cell-eval's `vcc`/`minimal`/`full`)

cell-eval registers 15+ metrics and exposes pre-configured profiles. GeneLab Benchmark currently reports macro-F1 + balanced accuracy as primary; supplementary metrics scattered across categories.

**Translation:** Define explicit metric profiles in v8+:
- `genelab_minimal` = macro-F1, balanced accuracy, mission discrimination
- `genelab_full` = + per-tissue F1, cross-mission transfer matrix entropy, confounder-conditional accuracy
- `genelab_fm` = + embedding-level metrics (silhouette, clustering agreement vs ground truth labels)

### 3. DE-level evaluation for scRNA-seq sub-track

GeneLab v6+ includes scRNA-seq tasks (RRRM-1, RRRM-2). For these, cell-eval's **`overlap_at_N`** (Wilcoxon-based DE gene set overlap) is directly applicable.

**Translation:** When evaluating FMs on scRNA-seq spaceflight tasks (cell-type-specific Flight vs. Ground DE), use cell-eval's DE metrics directly:
- `overlap_at_N` for top DE gene recovery
- `de_direction_match` for sign agreement
- `pr_auc` for FDR-calibrated DE recovery

This brings GeneLab Benchmark v8+ scRNA-seq evaluations into alignment with the broader virtual cell field.

---

## What GeneLab Benchmark Has That cell-eval Lacks

GeneLab's **LOMO design** (Leave-One-Mission-Out CV) and **tissue × mission matrix** are independent contributions that have no analog in cell-eval. The framework solves a different problem: can FMs trained on Earth-based data generalize to the spaceflight domain shift?

The GeneLab v7 finding — **PCA-LR (0.753) > scGPT (0.667) > Mouse-Geneformer (0.476)** — is a domain-shift failure pattern that cell-eval's perturbation-prediction framing wouldn't surface.

---

## Suggested v8+ Integration Path

### Optional drop-in (low priority)
- Install `cell-eval` for the scRNA-seq sub-track only (RRRM-1, RRRM-2)
- Apply `overlap_at_N` + `de_direction_match` to per-cell-type DE gene predictions
- Document as "harmonized scRNA-seq evaluation alongside cell-eval"

### Conceptual port (higher priority)
- Add **mission discrimination metric** (PDS-inspired) to v8+ as a third primary metric alongside macro-F1 and balanced accuracy
- Define explicit metric profiles (`genelab_minimal`, `genelab_full`, `genelab_fm`)
- Reframe v8+ paper: positioning section can note alignment with cell-eval design principles while preserving GeneLab's distinct LOMO+mission focus

### What NOT to do
- **Don't replace** macro-F1 / balanced accuracy with PDS — they answer different questions (classification accuracy vs. discrimination ranking)
- **Don't apply cell-eval to bulk RNA-seq tasks** — DES requires single-cell distributions; bulk-only models can't run Wilcoxon properly
- **Don't merge frameworks** — GeneLab and cell-eval have distinct, complementary scopes

---

## Why This Matters Strategically

The virtual cell field's evaluation standards are converging on cell-eval (Arc) + complementary axis-design frameworks (e.g., ORBIT-Cell). GeneLab Benchmark sits in a different niche (cross-mission domain shift on bulk + sc data) but should:

1. **Acknowledge the alignment** in the v8+ paper's related-work section
2. **Adopt PDS-style discrimination metrics** where biologically meaningful
3. **Use cell-eval directly for scRNA-seq sub-tasks** to maximize community comparability
4. **Defend the LOMO + mission matrix design** as the orthogonal contribution

---

## Reference

Internal source note: the original private prep memo is intentionally omitted
from the public repository. Public implementation decisions should rely on the
sources below and on `docs/V8_IMPLEMENTATION_RESEARCH.md`.

**Key papers:**
- cell-eval / VCC: Roohani et al. *Cell* 2025; 188(13):3370–3374
- STATE model: biorxiv.org/content/10.1101/2025.06.26.661135v1
- VCBench (12-metric benchmark): Mao et al. arXiv:2604.27646 (April 2026)
- Evaluation difficulty: Csendes et al. biorxiv:2026.02.14.705879 (Feb 2026)

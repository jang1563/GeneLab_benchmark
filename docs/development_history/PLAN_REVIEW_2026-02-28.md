# PLAN.md Detailed Review (2026-02-28)

## Findings (Severity Ordered)

1. **[Critical] Feature/label leakage risk (LFC definition can make classification trivially easier)**  
   Ref: `PLAN.md` L236, L545  
   If `Flight vs Ground` uses `LFC vs mission control` as input, control samples become structurally near zero, which can leak label signal. Also, inference requires mission-specific control in the test mission, reducing deployment realism.

2. **[Critical] Split design allows confounder/batch leakage**  
   Ref: `PLAN.md` L344, L306  
   `80/20 stratified` sample split for D/J can place same mission/study samples in both train/test, inflating performance. This is especially problematic for D3 (mission ID) and D4/D5 (confounder characterization).

3. **[Critical] Independence unit undefined (animal/cage-level leakage possible)**  
   Ref: `PLAN.md` L285, L413  
   Without group constraints, multi-tissue samples from the same animal or same cage/batch may cross split boundaries and overestimate transfer performance.

4. **[High] Composite-score normalization mismatches metric definitions**  
   Ref: `PLAN.md` L351, L370  
   Multi-class primary metric is `macro-F1`, but random baseline is written as `1/K accuracy`; normalization is not mathematically aligned.

5. **[High] Harmony as a bulk RNA-seq batch-correction baseline is weakly justified**  
   Ref: `PLAN.md` L276, L468  
   Harmony is mainly embedding-oriented (often single-cell). For bulk count workflows, ComBat/ComBat-seq, limma, RUV-family methods are more defensible.

6. **[High] D4/D5 confounder interpretation is causally fragile**  
   Ref: `PLAN.md` L311, L315  
   If strain/hardware is collinear with mission/duration/year, high predictive performance may reflect mission identification rather than isolated confounder effect.

7. **[High] Phase checkpoint ignores uncertainty**  
   Ref: `PLAN.md` L522  
   A single threshold (`A1 AUROC > 0.7`) is fragile without confidence intervals or permutation tests.

8. **[Medium] Mission-count/table inconsistency**  
   Ref: `PLAN.md` L165, L266  
   Liver is listed as `7+` missions while B1 is fixed `6x6`. Freeze criteria for mission inclusion are unclear.

9. **[Medium] Track-specific eligibility rules are missing**  
   Ref: `PLAN.md` L82, L100  
   Track 1 (GC-only) can fail sample requirements even when overall inclusion passes. Eligibility should be defined per track.

10. **[Medium] Duplicate-sample rule may create false positives**  
    Ref: `PLAN.md` L397  
    `Pearson r > 0.99` alone can flag biologically similar samples. Should combine checksum, metadata, technical replicate tags, and batch info.

11. **[Medium] Assay exclusion policy is not yet machine-enforced**  
    Ref: `PLAN.md` L88, L180  
    Category A exclusions are clear textually, but automatic whitelist/blacklist rules for `catalog_datasets.py` are not specified.

12. **[Medium] Reproducibility contract is missing**  
    Ref: `PLAN.md` L499  
    Benchmark credibility needs explicit env lock, data snapshot/versioning, and deterministic seed policy.

## Open Questions

1. For A/B/C, is the target primarily sample-level prediction or mission/group-level signature comparison?  
2. Should D4/D5 remain confounder quantification tasks, or shift to residualized prediction after adjustment?  
3. At v1.0 freeze, will liver remain fixed at 6 missions or expand to 7+ (with explicit inclusion priority)?

## Short Overall Assessment

The concept and scientific scope are strong. The largest practical risk is inflated performance due to leakage/split design. Fixing items 1–4 first will substantially improve paper defensibility and reproducibility.

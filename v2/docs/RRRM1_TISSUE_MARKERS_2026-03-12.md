# RRRM-1 Tissue Marker Reference

Date: 2026-03-12
Scope: first-pass broad annotation for RRRM-1 scRNA-seq downstream analysis

## Blood

Broad cell classes and starter markers:

- Erythroid: `Hba-a1`, `Hba-a2`, `Hbb-bs`, `Hbb-bt`, `Alas2`, `Klf1`
- B cell: `Cd79a`, `Cd79b`, `Ms4a1`, `Cd74`, `H2-Aa`
- T cell: `Cd3d`, `Cd3e`, `Trac`, `Ltb`, `Il7r`
- NK cell: `Nkg7`, `Klrd1`, `Klrb1c`, `Prf1`, `Ccl5`
- Monocyte / macrophage: `Lyz2`, `Lgals3`, `Ctss`, `Tyrobp`, `Fcgr3`
- Neutrophil: `S100a8`, `S100a9`, `Retnlg`, `Lcn2`, `Mmp8`
- Platelet / megakaryocyte: `Ppbp`, `Pf4`, `Gng11`, `Itga2b`, `Gp9`
- Dendritic cell: `Flt3`, `Ccr7`, `Itgax`, `H2-Ab1`, `Xcr1`

## Eye

Broad cell classes and starter markers:

- Epithelial: `Krt8`, `Krt18`, `Krt19`, `Epcam`, `Krt14`
- Corneal / conjunctival epithelial: `Krt12`, `Krt13`, `Krt14`, `Krt15`, `Tacstd2`
- Lens / crystallin-rich: `Cryaa`, `Cryab`, `Cryba1`, `Crybb2`, `Mip`
- Retinal neuronal: `Rbfox3`, `Snap25`, `Tubb3`, `Elavl4`, `Syt1`
- Photoreceptor-like: `Rho`, `Gnat1`, `Pde6a`, `Rcvrn`, `Crx`
- Muller glia / retinal glia: `Rlbp1`, `Slc1a3`, `Glul`, `Aqp4`, `Vim`
- Stromal / fibroblast: `Col1a1`, `Col1a2`, `Dcn`, `Lum`, `Pdgfra`
- Endothelial / perivascular: `Kdr`, `Klf2`, `Pecam1`, `Cdh5`, `Rgs5`
- Immune: `Ptprc`, `Lyz2`, `Tyrobp`, `Ctss`, `H2-Ab1`

## Muscle

Broad cell classes and starter markers:

- Myonuclear / muscle structural: `Acta1`, `Tnnt3`, `Tpm1`, `Myh1`, `Ckm`
- Satellite cell / myogenic progenitor: `Pax7`, `Myf5`, `Myod1`, `Myog`, `Vcam1`
- Fibro-adipogenic progenitor: `Pdgfra`, `Col1a1`, `Dcn`, `Lum`, `Pi16`
- Fibroblast / matrix: `Col1a1`, `Col3a1`, `Fn1`, `Dcn`, `Postn`
- Endothelial: `Pecam1`, `Cdh5`, `Kdr`, `Emcn`, `Klf2`
- Pericyte / smooth muscle: `Rgs5`, `Myl9`, `Acta2`, `Cspg4`, `Pdgfrb`
- Macrophage / myeloid: `Lyz2`, `Adgre1`, `Tyrobp`, `Ctss`, `Lgals3`
- T / NK lymphocyte: `Cd3d`, `Cd3e`, `Trac`, `Nkg7`, `Ccl5`

## Skin

Broad cell classes and starter markers:

- Basal keratinocyte: `Krt5`, `Krt14`, `Krt15`, `Dst`, `Trp63`
- Suprabasal / differentiating keratinocyte: `Krt1`, `Krt10`, `Lor`, `Flg`, `Sprr1b`
- Hair follicle / appendage epithelial: `Krt17`, `Krt79`, `Lhx2`, `Sox9`, `Shh`
- Fibroblast / dermal stromal: `Col1a1`, `Col1a2`, `Dcn`, `Lum`, `Pdgfra`
- Endothelial: `Pecam1`, `Cdh5`, `Kdr`, `Emcn`, `Klf2`
- Pericyte / vascular smooth muscle: `Rgs5`, `Acta2`, `Myl9`, `Tagln`, `Cspg4`
- Immune myeloid: `Lyz2`, `Tyrobp`, `Ctss`, `Fcgr3`, `Adgre1`
- T / NK lymphocyte: `Cd3d`, `Cd3e`, `Trac`, `Nkg7`, `Ccl5`
- Melanocyte-like: `Pmel`, `Dct`, `Mlana`, `Tyr`, `Tyrp1`

## Annotation Notes

- Use these as first-pass markers only. Final annotation should be confirmed from cluster-level differential expression, not single genes.
- For mixed clusters, prioritize coherent marker programs over one-off high-expression genes.
- If a tissue contains strong erythroid contamination, label it explicitly rather than forcing assignment to the resident tissue atlas.

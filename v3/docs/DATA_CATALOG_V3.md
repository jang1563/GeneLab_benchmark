# GeneLabBench v3 — Data Catalog

**Generated**: 2026-03-18
**Method**: OSDR API automated verification (`v3/scripts/catalog_v3_datasets.py`)

## Summary

| OSD ID | Organism | Assay | n | Title | Phase | Go? |
|--------|----------|-------|---|-------|-------|-----|
| OSD-112 | Caenorhabditis elega | Microarray | 9 | Microgravity effect on C. elegans N2/VC (CERI... | 1 | ⚠️  |
| OSD-113 | Caenorhabditis elega | Microarray | 6 | Expression Data from International C. elegans... | 1 | ✅ M |
| OSD-120 | Arabidopsis thaliana | RNA-seq, Imaging | 90 | Genetic dissection of the Arabidopsis spacefl... | 1 | ✅ R |
| OSD-207 | Drosophila melanogas | RNA-seq, Proteomics | 32 | Correlated Gene and Protein Expression in hea... | 1 | ✅ R |
| OSD-208 | Arabidopsis thaliana | RNA-seq, Microarray | 6 | Comparing RNA-Seq and microarray gene express... | 1 | ✅ R |
| OSD-37 | Arabidopsis thaliana | RNA-seq | 56 | Comparison of the spaceflight transcriptome o... | 1 | ✅ R |
| OSD-38 | Arabidopsis thaliana | Proteomics, RNA-seq | 36 | Proteomics and Transcriptomics analysis of Ar... | 1 | ✅ R |
| OSD-41 | ? | Microarray | 18 | Microgravity effect on C. elegans N2/VC (CERI... | 1 | ⚠️  |
| OSD-7 | Arabidopsis thaliana | Microarray | 36 | The Arabidopsis spaceflight transcriptome: a ... | 1 | ⚠️  |
| OSD-270 | Mus musculus | Spatial (Visium), RNA-seq | 46 | Bulk and spatially resolved transcriptional a... | 2 | ✅ R |
| OSD-402 | Mus musculus | RNA-seq | 16 | Single cell transcriptional profiling of femu... | 3 | ✅ R |
| OSD-403 | Mus musculus | RNA-seq | 16 | Single cell transcriptional profiling of hume... | 3 | ✅ R |
| OSD-404 | Mus musculus | RNA-seq | 16 | Single cell transcriptional profiling of peri... | 3 | ✅ R |
| OSD-405 | Mus musculus | RNA-seq | 16 | Single cell transcriptional profiling of sple... | 3 | ✅ R |
| OSD-240 | Mus musculus | RNA-seq | 20 | Transcriptional analysis of dorsal skin from ... | 4 | ✅ R |
| OSD-241 | Mus musculus | RNA-seq | 19 | Transcriptional analysis of femoral skin from... | 4 | ✅ R |
| OSD-243 | Mus musculus | RNA-seq | 53 | Transcriptional analysis of dorsal skin from ... | 4 | ✅ R |
| OSD-247 | Mus musculus | RNA-seq | 54 | Transcriptional analysis of colon from mice f... | 4 | ✅ R |
| OSD-248 | Mus musculus | RNA-seq | 58 | Transcriptional analysis of lung from mice fl... | 4 | ✅ R |
| OSD-253 | Mus musculus | RNA-seq | 112 | Transcriptional analysis of kidneys from mice... | 4 | ✅ R |
| OSD-109 | Mus musculus | Microarray | 15 | Delayed Cardiomyocyte Response to Total Body ... | 5 | ✅ M |
| OSD-202 | Mus musculus | RNA-seq, Microarray, Methylation | 84 | Low dose (0.04 Gy) irradiation (LDR) and hind... | 5 | ✅ R |
| OSD-211 | Mus musculus | RNA-seq | 21 | Transcriptomic analysis of spleens from mice ... | 5 | ✅ R |
| OSD-237 | Mus musculus | RNA-seq | 21 | Transcriptomic analysis of skin from mice sub... | 5 | ✅ R |
| OSD-73 | Homo sapiens | Microarray | 134 | Transcriptome Profiles in Normal Human Bronch... | 5 | ✅ M |
| OSD-804 | Mus musculus | MicroCT | 170 | 37-Day microgravity exposure in 16-Week femal... | 6 | ⚠️  |

## Phase 1 Details

### OSD-207 — ✅ RNA-seq
- **Title**: Correlated Gene and Protein Expression in heads from Drosophila reared in microgravity
- **Organism**: Drosophila melanogaster
- **Assay types**: RNA-seq, Proteomics
- **Sample count**: 32
- **Factors**: Spaceflight (Space Flight), Genotype (genotype), Strain (strain)
  - Strain: Canton-S, w1118
  - Spaceflight: Ground Control, Space Flight
  - Genotype: KCNQ370, Sei ts1, Wild Type
- **Expected**: Drosophila head RNA-seq, spaceflight
- **Description**: Omics analyses of RNA and protein isolated from heads of microgravity reared adult Drosophila.
- **First data file**: `GLDS-207_rna-seq_01Mix_1-KCNQ_370-Female.tar.gz`

### OSD-7 — ⚠️ Microarray (not RNA-seq)
- **Title**: The Arabidopsis spaceflight transcriptome: a comparison of whole plants to discrete root, hypocotyl and shoot responses to the orbital environment
- **Organism**: Arabidopsis thaliana
- **Assay types**: Microarray
- **Sample count**: 36
- **Factors**: Spaceflight (Space Flight), Organism Part (organism part)
  - Organism Part: Plant Roots, Plant Shoots, Whole Organism, hypocotyl
  - Spaceflight: Ground Control, Space Flight
- **Expected**: Arabidopsis RNA-seq, spaceflight (ABRS)
- **Description**: Arabidopsis thaliana was evaluated for its response to the spaceflight environment in three replicated experiments on the International Space Station. Two approaches were used; GFP reporter genes were used to collect gene expression data in real time within unique GFP imaging hardware, and plants we
- **First data file**: `GLDS-7_microarray_8311_FLTAGES2_25_Ferl.Paul_(ATH1-121501).CEL`

### OSD-37 — ✅ RNA-seq
- **Title**: Comparison of the spaceflight transcriptome of four commonly used Arabidopsis thaliana ecotypes
- **Organism**: Arabidopsis thaliana
- **Assay types**: RNA-seq
- **Sample count**: 56
- **Factors**: Spaceflight (Space Flight), Ecotype (ecotype)
  - Spaceflight: Ground Control, Space Flight
  - Ecotype: Col-0, Cvi-0, Ler-0, Ws-2
- **Expected**: Arabidopsis RNA-seq, spaceflight (BRIC)
- **Description**: This experiment compared the spaceflight transcriptomes of four commonly used natural variants (ecotypes) of Arabidopsis thaliana using RNAseq. In nature, Arabidopsis is a native of Europe/Asia/Northwestern Africa and is found across the globe growing in a wide range of environments. The geographica
- **First data file**: `GLDS-37_rna-seq_Rep1-FL-A1.tar.gz`

### OSD-38 — ✅ RNA-seq
- **Title**: Proteomics and Transcriptomics analysis of Arabidopsis Seedlings in Microgravity
- **Organism**: Arabidopsis thaliana
- **Assay types**: Proteomics, RNA-seq
- **Sample count**: 36
- **Factors**: Spaceflight (Space Flight), Sample Preservation Method (Preservation, Biological)
  - Sample Preservation Method: Liquid Nitrogen, RNALater
  - Spaceflight: Ground Control, Space Flight
- **Expected**: Arabidopsis RNA-seq + proteomics
- **Description**: On Earth plants are constantly exposed to a gravitational field of 1G. Gravity affects a plant in every step of its development. Germinating seedlings orient their radicle and hypocotyl and growing plants position organs at a specific Gravitropic Set-point Angle dictated by the asymmetric distributi
- **First data file**: `GLDS-38_membrane_proteome_WyattLab_membrane.csv,GLDS-38_membrane_proteome_WyattLab_membrane.xlsx`

### OSD-120 — ✅ RNA-seq
- **Title**: Genetic dissection of the Arabidopsis spaceflight transcriptome: Are some responses dispensable for the physiological adaptation of plants to spaceflight?
- **Organism**: Arabidopsis thaliana
- **Assay types**: RNA-seq, Imaging
- **Sample count**: 90
- **Factors**: Ecotype (Ecotype), Spaceflight (Space Flight), Treatment (treatment)
  - Treatment: Dark Treatment, Light Treatment
  - Ecotype: Col-0, Col-0 PhyD, Wassilewskija ecotype
  - Spaceflight: Ground Control, Space Flight
- **Expected**: Arabidopsis RNA-seq, spaceflight (ISS)
- **Description**: Experimentation on the International Space Station has reached the stage where repeated and nuanced transcriptome studies are beginning to illuminate the structural and metabolic differences between plants grown in space compared to plants on the Earth. Genes that are important in establishing the s
- **First data file**: `GLDS-120_rna-seq_GSM2493759.tar`

### OSD-208 — ✅ RNA-seq
- **Title**: Comparing RNA-Seq and microarray gene expression data in two zones of the Arabidopsis root apex relevant to spaceflight.
- **Organism**: Arabidopsis thaliana
- **Assay types**: RNA-seq, Microarray
- **Sample count**: 6
- **Factors**: Tissue (tissue)
  - Tissue: 
- **Expected**: Arabidopsis root RNA-seq (NOT spaceflight)
- **Description**: Premise of the study: The root apex is an important region involved in environmental sensing, but comprises a very small part of the root. Obtaining root apex transcriptomes is therefore challenging when the samples are limited. The feasibility of using tiny root sections for transcriptome analysis 
- **First data file**: `GLDS-208_rna-seq_SRR7287361_1.fastq.gz, GLDS-208_rna-seq_SRR7287361_2.fastq.gz, GLDS-208_rna-seq_SRR7287362_1.fastq.gz, GLDS-208_rna-seq_SRR7287362_2.fastq.gz, GLDS-208_rna-seq_SRR7287363_1.fastq.gz, GLDS-208_rna-seq_SRR7287363_2.fastq.gz, GLDS-208_rna-seq_SRR7287364_1.fastq.gz, GLDS-208_rna-seq_SRR7287364_2.fastq.gz`

### OSD-113 — ✅ Microarray (as expected)
- **Title**: Expression Data from International C. elegans Experiment 1st (ICE-FIRST)
- **Organism**: Caenorhabditis elegans
- **Assay types**: Microarray
- **Sample count**: 6
- **Factors**: Spaceflight (Space Flight)
  - Spaceflight: Ground control, Space Flight
- **Expected**: C. elegans spaceflight (microarray?)
- **Description**: The effect of microgravity on gene expression in C.elegans was comprehensively analysed by DNA microarray (Agilent). Overall design: Mix staged C. elegans N2 was exposed microgravity for 10 days.
- **First data file**: `GLDS-113_microarray_GSE71771_processed.zip`

### OSD-112 — ⚠️ Microarray (not RNA-seq)
- **Title**: Microgravity effect on C. elegans N2/VC (CERISE, 4 days)
- **Organism**: Caenorhabditis elegans
- **Assay types**: Microarray
- **Sample count**: 9
- **Factors**: Altered Gravity (Gravity, Altered), Spaceflight (Space Flight)
  - Altered Gravity: 1G by centrifugation, 1G on Earth, uG
  - Spaceflight: Ground Control, Space Flight
- **Expected**: C. elegans spaceflight
- **Description**: Microgravity effect on C. elegans gene expression was analysed by whole genome microarray. The worms were cultivated under microgravity for 4 days in the Japanese Module of the International Space Station. C. elegans N2 was exposed microgravity for 4 days. The worms synchronously were cultivated fro
- **First data file**: `GLDS-112_microarray_GSM1845018_US84100219_252018610315_S01_GE1_107_Sep09_1_1.txt.gz`

### OSD-41 — ⚠️ Microarray (not RNA-seq)
- **Title**: Microgravity effect on C. elegans N2/VC (CERISE, 8days)
- **Organism**: 
- **Assay types**: Microarray
- **Sample count**: 18
- **Factors**: gravity (Gravitation)
- **Expected**: C. elegans spaceflight
- **Description**: Microgravity effect on C. elegans gene expression was analysed by whole genome microarray. The worms were cultivated under microgravity for 8 days in the Japanese Module of the International Space Station. The samples of this study were divided three experimental groups: 1. microgravity for 8 days 2
- **First data file**: `GSM675866.txt`


## Phase 2 Details

### OSD-270 — ✅ RNA-seq
- **Title**: Bulk and spatially resolved transcriptional analysis of hearts from mice flown on the RR-3
- **Organism**: Mus musculus
- **Assay types**: Spatial (Visium), RNA-seq
- **Sample count**: 46
- **Factors**: Spaceflight (Space Flight)
  - Spaceflight: Ground Control, Space Flight
- **Expected**: RR-3 heart Visium spatial
- **Description**: The Rodent Research-3 (RR-3) mission was sponsored by the pharmaceutical company Eli Lilly and Co. and the Center for the Advancement of Science in Space to study the effectiveness of a potential countermeasure for the loss of muscle and bone mass that occurs during spaceflight. Twenty BALB/c, 12-we
- **First data file**: `GLDS-270_SpatialTranscriptomics_Lev1_FLT_Rep1_TecRep1_V19L29_043_C1_F1_R1_raw.fastq.gz,GLDS-270_SpatialTranscriptomics_Lev1_FLT_Rep1_TecRep1_V19L29_043_C1_F1_R2_raw.fastq.gz`


## Phase 3 Details

### OSD-402 — ✅ RNA-seq
- **Title**: Single cell transcriptional profiling of femur bone marrow from mice flown on Rodent Research Reference Mission-2 (RRRM-2)
- **Organism**: Mus musculus
- **Assay types**: RNA-seq
- **Sample count**: 16
- **Factors**: Spaceflight (Space Flight), Age (age)
  - Age: 
  - Spaceflight: Ground Control, Space Flight
- **Expected**: RRRM-2 femur bone marrow scRNA-seq
- **Description**: In the Rodent Research Reference Mission (RRRM-2), forty female C57BL/6NTac mice were flown on the International Space Station. To assess differences in outcomes due to age, twenty 12 week-old and twenty 29 week-old mice were flown, respectively. To directly assess spaceflight effects, half of the y
- **First data file**: `GLDS-402_scRNA-Seq_RRRM2_Femur_BM_FLT_LAR_OLD_FO16_R1_raw.fastq.gz,GLDS-402_scRNA-Seq_RRRM2_Femur_BM_FLT_LAR_OLD_FO16_R2_raw.fastq.gz`

### OSD-403 — ✅ RNA-seq
- **Title**: Single cell transcriptional profiling of humerus bone marrow from mice flown on Rodent Research Reference Mission-2 (RRRM-2)
- **Organism**: Mus musculus
- **Assay types**: RNA-seq
- **Sample count**: 16
- **Factors**: Spaceflight (Space Flight), Age (age)
  - Age: 
  - Spaceflight: Ground Control, Space Flight
- **Expected**: RRRM-2 humerus bone marrow scRNA-seq
- **Description**: In the Rodent Research Reference Mission (RRRM-2), forty female C57BL/6NTac mice were flown on the International Space Station. To assess differences in outcomes due to age, twenty 12 week-old and twenty 29 week-old mice were flown, respectively. To directly assess spaceflight effects, half of the y
- **First data file**: `GLDS-403_scRNA-Seq_RRRM2_Humerus_BM_FLT_LAR_OLD_FO16_R1_raw.fastq.gz,GLDS-403_scRNA-Seq_RRRM2_Humerus_BM_FLT_LAR_OLD_FO16_R2_raw.fastq.gz`

### OSD-404 — ✅ RNA-seq
- **Title**: Single cell transcriptional profiling of peripheral blood mononuclear cells (PBMCs) from mice flown on Rodent Research Reference Mission-2 (RRRM-2)
- **Organism**: Mus musculus
- **Assay types**: RNA-seq
- **Sample count**: 16
- **Factors**: Spaceflight (Space Flight), Age (age)
  - Age: 
  - Spaceflight: Ground Control, Space Flight
- **Expected**: RRRM-2 PBMCs scRNA-seq
- **Description**: In the Rodent Research Reference Mission (RRRM-2), forty female C57BL/6NTac mice were flown on the International Space Station. To assess differences in outcomes due to age, twenty 12 week-old and twenty 29 week-old mice were flown, respectively. To directly assess spaceflight effects, half of the y
- **First data file**: `GLDS-404_scRNA-Seq_RRRM2_PBMC_FLT_LAR_OLD_FO16_R1_raw.fastq.gz,GLDS-404_scRNA-Seq_RRRM2_PBMC_FLT_LAR_OLD_FO16_R2_raw.fastq.gz`

### OSD-405 — ✅ RNA-seq
- **Title**: Single cell transcriptional profiling of spleens from mice flown on Rodent Research Reference Mission-2 (RRRM-2)
- **Organism**: Mus musculus
- **Assay types**: RNA-seq
- **Sample count**: 16
- **Factors**: Spaceflight (Space Flight), Age (age)
  - Age: 
  - Spaceflight: Ground Control, Space Flight
- **Expected**: RRRM-2 spleen scRNA-seq
- **Description**: In the Rodent Research Reference Mission (RRRM-2), forty female C57BL/6NTac mice were flown on the International Space Station. To assess differences in outcomes due to age, twenty 12 week-old and twenty 29 week-old mice were flown, respectively. To directly assess spaceflight effects, half of the y
- **First data file**: `GLDS-405_scRNA-Seq_RRRM2_SPL_FLT_LAR_OLD_FO16_R1_raw.fastq.gz,GLDS-405_scRNA-Seq_RRRM2_SPL_FLT_LAR_OLD_FO16_R2_raw.fastq.gz`


## Phase 4 Details

### OSD-253 — ✅ RNA-seq
- **Title**: Transcriptional analysis of kidneys from mice flown on the RR-7 mission
- **Organism**: Mus musculus
- **Assay types**: RNA-seq
- **Sample count**: 112
- **Factors**: Spaceflight (Space Flight), Strain (strain), Duration (duration), Treatment (treatment)
  - Strain: C3H/HeJ, C57BL/6J
  - Spaceflight: Basal Control, Ground Control, Ground Control Rerun, Space Flight, Vivarium Control
  - Duration: 
  - Treatment: Blue Light, Not Applicable, White light
- **Expected**: Kidney RR-7 RNA-seq (extends kidney)
- **Description**: The objective of the Rodent Research-7 mission (RR-7) was to study the impact of the space environment on the gut microbiota of two strains of mice and how any changes in-turn affect the immune system, metabolic system, and circadian or daily rhythms. To this end, ten 11-week-old female C57BL/6J and
- **First data file**: `GLDS-253_rna-seq_Mmus_C3H-HeJ_KDN_BSL_0days_Rep1_B2_R1_raw.fastq.gz, GLDS-253_rna-seq_Mmus_C3H-HeJ_KDN_BSL_0days_Rep1_B2_R2_raw.fastq.gz`

### OSD-240 — ✅ RNA-seq
- **Title**: Transcriptional analysis of dorsal skin from mice flown on the RR-5 mission
- **Organism**: Mus musculus
- **Assay types**: RNA-seq
- **Sample count**: 20
- **Factors**: Spaceflight (Space Flight)
  - Spaceflight: Ground Control, Space Flight
- **Expected**: Dorsal skin RR-5 RNA-seq
- **Description**: The objective of the Rodent Research-5 (RR-5) study was to evaluate bone loss in mice during spaceflight and to determine if treatment with a modified version of NEL-like molecule-1 (NELL-1) can reduce or prevent bone loss that would otherwise occur during spaceflight. To this end, a cohort of forty
- **First data file**: `GLDS-240_rna-seq_Mmus_BAL-TAL_DSKN_FLT_Rep10_F10_R1_raw.fastq.gz, GLDS-240_rna-seq_Mmus_BAL-TAL_DSKN_FLT_Rep10_F10_R2_raw.fastq.gz`

### OSD-241 — ✅ RNA-seq
- **Title**: Transcriptional analysis of femoral skin from mice flown on the RR-5 mission
- **Organism**: Mus musculus
- **Assay types**: RNA-seq
- **Sample count**: 19
- **Factors**: Spaceflight (Space Flight)
  - Spaceflight: Ground Control, Space Flight
- **Expected**: Femoral skin RR-5 RNA-seq
- **Description**: The objective of the Rodent Research-5 (RR-5) study was to evaluate bone loss in mice during spaceflight and to determine if treatment with a modified version of NEL-like molecule-1 (NELL-1) can reduce or prevent bone loss that would otherwise occur during spaceflight. To this end, a cohort of forty
- **First data file**: `GLDS-241_rna-seq_Mmus_BAL-TAL_FSKN_FLT_Rep1_F1_R1_raw.fastq.gz, GLDS-241_rna-seq_Mmus_BAL-TAL_FSKN_FLT_Rep1_F1_R2_raw.fastq.gz`

### OSD-243 — ✅ RNA-seq
- **Title**: Transcriptional analysis of dorsal skin from mice flown on the RR-6 mission
- **Organism**: Mus musculus
- **Assay types**: RNA-seq
- **Sample count**: 53
- **Factors**: Spaceflight (Space Flight), Duration (duration), Dissection Condition (dissection), Euthanasia Location (Euthanasia Location)
  - Euthanasia Location: On Earth, On ISS
  - Dissection Condition: Carcass, Upon euthanasia
  - Duration: 
  - Spaceflight: Basal Control, Ground Control, Space Flight
- **Expected**: Dorsal skin RR-6 RNA-seq
- **Description**: The objective of the Rodent Research-6 (RR-6) study was to evaluate muscle atrophy in mice during spaceflight and to test the efficacy of a novel therapeutic to mitigate muscle wasting. The experiment involved an implantable subcutaneous nanochannel delivery system (nDS; between scapula), which deli
- **First data file**: `GLDS-243_rna-seq_Mmus_C57-6T_DSKN_BSL_ISS-T_Rep1_B2_R1_raw.fastq.gz, GLDS-243_rna-seq_Mmus_C57-6T_DSKN_BSL_ISS-T_Rep1_B2_R2_raw.fastq.gz`

### OSD-248 — ✅ RNA-seq
- **Title**: Transcriptional analysis of lung from mice flown on the RR-6 mission
- **Organism**: Mus musculus
- **Assay types**: RNA-seq
- **Sample count**: 58
- **Factors**: Spaceflight (Space Flight), Duration (duration), Dissection Condition (dissection), Euthanasia Location (Euthanasia Location)
  - Euthanasia Location: On Earth, On ISS
  - Dissection Condition: Carcass, Upon euthanasia
  - Duration: 
  - Spaceflight: Basal Control, Ground Control, Space Flight
- **Expected**: Lung RR-6 RNA-seq
- **Description**: The objective of the Rodent Research-6 (RR-6) study was to evaluate muscle atrophy in mice during spaceflight and to test the efficacy of a novel therapeutic to mitigate muscle wasting. The experiment involved an implantable subcutaneous nanochannel delivery system (nDS; between scapula), which deli
- **First data file**: `GLDS-248_rna-seq_Mmus_C57-6T_LNG_BSL_ISS-T_Rep1_B2_R1_raw.fastq.gz, GLDS-248_rna-seq_Mmus_C57-6T_LNG_BSL_ISS-T_Rep1_B2_R2_raw.fastq.gz`

### OSD-247 — ✅ RNA-seq
- **Title**: Transcriptional analysis of colon from mice flown on the RR-6 mission
- **Organism**: Mus musculus
- **Assay types**: RNA-seq
- **Sample count**: 54
- **Factors**: Spaceflight (Space Flight), Duration (duration), Dissection Condition (dissection), Euthanasia Location (Euthanasia Location)
  - Euthanasia Location: On Earth, On ISS
  - Dissection Condition: Carcass, Upon euthanasia
  - Duration: 
  - Spaceflight: Basal Control, Ground Control, Space Flight
- **Expected**: Colon RR-6 RNA-seq (new tissue)
- **Description**: The objective of the Rodent Research-6 (RR-6) study was to evaluate muscle atrophy in mice during spaceflight and to test the efficacy of a novel therapeutic to mitigate muscle wasting. The experiment involved an implantable subcutaneous nanochannel delivery system (nDS; between scapula), which deli
- **First data file**: `GLDS-247_rna-seq_Mmus_C57-6T_CLN_BSL_ISS-T_Rep1_B1_R1_raw.fastq.gz, GLDS-247_rna-seq_Mmus_C57-6T_CLN_BSL_ISS-T_Rep1_B1_R2_raw.fastq.gz`


## Phase 5 Details

### OSD-202 — ✅ RNA-seq
- **Title**: Low dose (0.04 Gy) irradiation (LDR) and hindlimb unloading (HLU) microgravity in mice: brain transcriptomic and epigenomic data
- **Organism**: Mus musculus
- **Assay types**: RNA-seq, Microarray, Methylation
- **Sample count**: 84
- **Factors**: Ionizing Radiation (ionizing radiation), Hindlimb Unloading (Hindlimb Suspension), Time of Sample Collection After Treatment (Time of Sample Collection After Treatment)
  - Hindlimb Unloading: Hindlimb Unloaded, Normally Loaded Control
  - Ionizing Radiation: cobalt-57 gamma radiation, non-irradiated
  - Time of Sample Collection After Treatment: 
- **Expected**: Retina+brain, LDR+HLU, RNA-seq
- **Description**: The purpose of the present study was to evaluate damage in brain and eye in a ground-based model for spaceflight which includes prolonged unloading and low-dose radiation. Low-dose/Low-dose-rate (LDR) gamma-radiation using 57Co plates (0.04 Gy) was delivered whole-body to mature, 6-month old female 
- **First data file**: `GLDS-202_rna-seq_CFG2006_S1_L001_R1_001.fastq.gz`

### OSD-211 — ✅ RNA-seq
- **Title**: Transcriptomic analysis of spleens from mice subjected to chronic low-dose radiation, hindlimb unloading or a combination of both
- **Organism**: Mus musculus
- **Assay types**: RNA-seq
- **Sample count**: 21
- **Factors**: Ionizing Radiation (Ionizing Radiation), Hindlimb Unloading (Hindlimb Suspension)
  - Ionizing Radiation: Irradiated with 0.04 Gy 57Co, non-irradiated
  - Hindlimb Unloading: Hindlimb Unloaded, Normally Loaded Control
- **Expected**: Spleen, LDR+HLU, RNA-seq
- **Description**: The purpose of this study was to evaluate transcriptional changes in mouse spleens using a ground-based model for spaceflight. This model includes prolonged unloading and low-dose irradiation. Low-dose-rate gamma-radiation was delivered to 6-month old female C57BL/6J mice using 57Co plates (0.04 Gy)
- **First data file**: `GLDS-211_rna-seq_Mmus_C57-6J_SPL_HLLC_IRC_Rep1_M18_R1_raw.fastq.gz, GLDS-211_rna-seq_Mmus_C57-6J_SPL_HLLC_IRC_Rep1_M18_R2_raw.fastq.gz`

### OSD-237 — ✅ RNA-seq
- **Title**: Transcriptomic analysis of skin from mice subjected to chronic low-dose radiation, hindlimb unloading or a combination of both
- **Organism**: Mus musculus
- **Assay types**: RNA-seq
- **Sample count**: 21
- **Factors**: Ionizing Radiation (Ionizing Radiation), Hindlimb Unloading (Hindlimb Suspension)
  - Ionizing Radiation: Irradiated with 0.04 Gy 57Co, non-irradiated
  - Hindlimb Unloading: Hindlimb Unloaded, Normally Loaded Control
- **Expected**: Skin, LDR+HLU, RNA-seq
- **Description**: The purpose of this study was to evaluate transcriptional changes in mouse skin using a ground-based model for spaceflight. This model includes prolonged unloading and low-dose irradiation. Low-dose-rate gamma-radiation was delivered to 6-month old female C57BL/6J mice using 57Co plates (0.04 Gy) to
- **First data file**: `GLDS-237_rna-seq_Mmus_C57-6J_DSKN_HLLC_IRC_Rep1_M22_R1_raw.fastq.gz, GLDS-237_rna-seq_Mmus_C57-6J_DSKN_HLLC_IRC_Rep1_M22_R2_raw.fastq.gz`

### OSD-73 — ✅ Microarray (as expected)
- **Title**: Transcriptome Profiles in Normal Human Bronchial Epithelial Cells after Exposure to gamma-rays and different HZE particles
- **Organism**: Homo sapiens
- **Assay types**: Microarray
- **Sample count**: 134
- **Factors**: Ionizing Radiation (ionizing radiation), Absorbed Radiation Dose (absorbed radiation dose), Time of Sample Collection After Treatment (Time of Sample Collection After Treatment)
  - Time of Sample Collection After Treatment: 
  - Absorbed Radiation Dose: 
  - Ionizing Radiation: Cesium-137 gamma radiation, Fe-56 ion radiation, Si-28 ion radiation, non-irradiated, sham-irradiated
- **Expected**: Radiation, microarray (verify)
- **Description**: Distinct transcriptome profiles in response to low-LET and high-LET, and different radiation qualities of HZE particles. Total RNA obtained from HBEC3KT cells after 1, 4, 12 and 24 hours of radiation. Mock-irradiated samples at each time point and control samples before radiation (0 hour) were also 

### OSD-109 — ✅ Microarray (as expected)
- **Title**: Delayed Cardiomyocyte Response to Total Body Particle Radiation Exposure - Identification of Regulatory Gene Network [iron]
- **Organism**: Mus musculus
- **Assay types**: Microarray
- **Sample count**: 15
- **Factors**: Time of Sample Collection After Treatment (Time of Sample Collection After Treatment), Ionizing Radiation (ionizing radiation)
  - Ionizing Radiation: Fe-56 ion radiation, sham-irradiated
  - Time of Sample Collection After Treatment: 
- **Expected**: Radiation, microarray (verify)
- **Description**: We examined molecular responses using transcriptome profiling in isolated left ventricular murine cardiomyocytes to 90 cGy, 1 GeV proton (1H) and 15 cGy, 1 GeV/nucleon (n) iron (56Fe) particles 1, 3, 7, 14 and 28 days after exposure. Unsupervised clustering analysis of gene expression segregated sam
- **First data file**: `GLDS-109_microarray_GSM1684875_DG_A4_IRR_24h_7.CEL`


## Phase 6 Details

### OSD-804 — ⚠️ Non-omics (MicroCT)
- **Title**: 37-Day microgravity exposure in 16-Week female C57BL/6J mice during the NASA Rodent Research 1 mission is associated with bone loss specific to weight-bearing skeletal sites (femur and vertebrae, micro computed tomography)
- **Organism**: Mus musculus
- **Assay types**: MicroCT
- **Sample count**: 170
- **Factors**: Spaceflight (Space Flight)
  - Spaceflight: Basal Control, Ground Control, Space Flight, Vivarium Control
- **Expected**: RR-1 MicroCT bone imaging
- **Description**: Exposure to weightlessness in microgravity and elevated space radiation are associated with rapid bone loss in mammals, but questions remain about their mechanisms of action and relative importance. In this study, we tested the hypothesis that bone loss during spaceflight in Low Earth Orbit is prima
- **First data file**: `Not Available`

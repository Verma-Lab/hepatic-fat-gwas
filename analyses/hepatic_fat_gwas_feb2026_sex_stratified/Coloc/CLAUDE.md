# Coloc: GTEx Colocalization Analysis (PMBB Hepatic Fat, Sex-Stratified)

## Overview
Colocalization of sex-stratified/combined hepatic fat GWAS signals against GTEx v8 eQTLs. See `../LeadSNPs/CLAUDE.md` for lead SNP discovery and full sumstats provenance. Same input sumstats as the LeadSNPs analysis: `PMBB_ALL_F`/`PMBB_ALL_M` (from this same `hepatic_fat_gwas_feb2026_sex_stratified/` directory) and `PMBB_ALL_ALL` (from `../../hepatic_fat_gwas_jan2026/`).

Originally developed locally at `/Users/agaro/Documents/PhD/meno_effects_liver/`; this copy on remote lives at `analyses/hepatic_fat_gwas_feb2026_sex_stratified/Coloc/`.

## Status (as of 2026-08-18)
- `coloc_hepatic_fat_gtex.py` / `.R` — GTEx v8 eQTL colocalization, run 2026-05-11–17
- `coloc_results_v2.tsv` is the current/final version (an earlier `coloc_results.tsv` was superseded and not carried over to remote)
- `build_eqtl_locus_table.py` — builds the locus table used as input to the coloc scripts, from lead SNPs/loci identified in `../LeadSNPs/`
- Last activity on this analysis before a 3-month gap; picked back up 2026-09-04 to consolidate local work onto remote

(Scripts live flat at the top of this directory, not in a nested `scripts/` folder — matches how the rest of `hepatic_fat_gwas_feb2026_sex_stratified/` puts single scripts at the top level rather than nesting them.)

## GTEx v8 tissues used (31, same eQTL set as the FUMA runs — see `../LeadSNPs/CLAUDE.md`)
Adipose Subcutaneous, Adipose Visceral Omentum, Adrenal Gland, Cells EBV-transformed lymphocytes, Whole Blood, Artery Aorta, Artery Coronary, Artery Tibial, Brain Amygdala, Brain Anterior cingulate cortex BA24, Brain Caudate basal ganglia, Brain Cerebellar Hemisphere, Brain Cerebellum, Brain Cortex, Brain Frontal Cortex BA9, Brain Hippocampus, Brain Hypothalamus, Brain Nucleus accumbens basal ganglia, Brain Putamen basal ganglia, Brain Spinal cord cervical c-1, Brain Substantia nigra, Breast Mammary Tissue, Colon Sigmoid, Colon Transverse, Esophagus Gastroesophageal Junction, Esophagus Mucosa, Esophagus Muscularis, Heart Atrial Appendage, Heart Left Ventricle, Kidney Cortex, Liver

## File index

| File | Description | Status |
|------|-------------|--------|
| `coloc_results_v2.tsv` | Colocalization results, current version | Keep |
| `coloc_hepatic_fat_gtex.py` / `.R` | Coloc analysis scripts | Keep |
| `build_eqtl_locus_table.py` | Builds locus table input for coloc | Keep |

## Pipeline timeline
See `../LeadSNPs/CLAUDE.md` for the full reconstructed pipeline timeline (EDA → phenotyping → covariates → GWAS → FUMA prep → this colocalization step).

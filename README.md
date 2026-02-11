# Sex-Specific and Menopause-Specific Genetic Effects on Liver Health

## Overview
This project investigates sex-specific and menopause-specific genetic effects on hepatic fat accumulation and liver health using CT imaging-derived phenotypes from the Penn Medicine Biobank. The study examines how genetic variants influencing hepatic fat may exhibit different effects between males and females, and how these effects change across menopausal status in women.

## Background
- Sexual dimorphism in liver disease is well-documented, with distinct patterns of metabolic dysfunction-associated steatotic liver disease (MASLD) between men and women
- Women have a 2.4-fold increased risk of MASLD development after menopause compared to pre-menopausal women
- Estrogen's protective effects on liver health decline after menopause, creating a critical biological transition period
- Menopause represents a unique window to study sex-specific genetic architecture and hormone-modulated genetic effects on liver health

## Research Hypothesis
Genetic variants known to influence hepatic fat accumulation (e.g., PNPLA3, TM6SF2, GPAM, APOE) will exhibit:
  - **Baseline GWAS**: All participants to identify hepatic fat-associated variants
  - **Sex-stratified GWAS**: Separate analyses in males vs. females to identify sex-specific effects
  - **SNP × menopause interaction testing**: Women-only analysis testing for differential genetic effects between pre-menopausal and post-menopausal women

## Methods

### GWAS Pipeline
- **Tool**: SAIGE (Scalable and Accurate Implementation of GEneralized mixed model)
- **Genotype data**: Penn Medicine Biobank imputed genotypes
- **Sample size**: ~18,000 participants with both genotype and imaging data
- **Analysis approach**: 
  - **Baseline GWAS**: All participants to identify hepatic fat-associated variants
  - **Sex-stratified GWAS**: Separate analyses in males vs. females to identify sex-specific effects
  - **SNP × menopause interaction testing**: Women-only analysis testing for differential genetic effects between pre-menopausal and post-menopausal women
  - Adjustment for age, genetic ancestry PCs, and technical covariates

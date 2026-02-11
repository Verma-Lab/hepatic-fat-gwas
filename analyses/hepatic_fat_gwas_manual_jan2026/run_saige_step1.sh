#!/bin/bash

set -euo pipefail

SIF=saige_1.5.0.sif
SINGULARITY_CMD="singularity exec \
  -B /project/verma_shared:/project \
  -B /static:/static \
  ${SIF}"

# Create output directory
OUTDIR=/project/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_manual_jan2026/PMBB_ALL_ALL
mkdir -p ${OUTDIR}

# Step 1: Fit null GLMM model
${SINGULARITY_CMD} /app/.pixi/envs/default/bin/Rscript /usr/local/bin/step1_fitNULLGLMM.R \
  --plinkFile=/static/PMBB/PMBB-Release-2024-3.0/Imputed/common_snps_LD_pruned/PMBB-Release-2024-3.0_genetic_imputed.commonsnps \
  --phenoFile=/project/projects/Liver_IDPs/covariates/hepatic_fat_all_covariates.csv \
  --phenoCol=hepatic_fat \
  --covarColList=Sex,Sample_age,PC1,PC2,PC3,PC4,PC5,PC6 \
  --qCovarColList=Sex \
  --sampleIDColinphenoFile=PMBB_ID \
  --traitType=quantitative \
  --outputPrefix=/project/projects/Liver_IDPs/analyses/hepatic_fat_gwas_manual_jan2026/PMBB_ALL_ALL/step1_null_model \
  --nThreads=4


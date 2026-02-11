#!/bin/bash
set -euo pipefail

SIF=saige_1.5.0.sif
BGEN_PREFIX=/static/PMBB/PMBB-Release-2024-3.0/Imputed/chunked_bgen_files/PMBB-Release-2024-3.0_genetic_imputed
CHUNK_LIST=/project/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_jan2026/bgen_chunk_list.txt
OUTDIR=/project/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_manual_jan2026/PMBB_ALL_ALL

# Get the chunk for this array job based on LSB_JOBINDEX
CHUNK=$(sed -n "${LSB_JOBINDEX}p" ${CHUNK_LIST})

if [ -z "${CHUNK}" ]; then
    echo "ERROR: Could not read chunk for index ${LSB_JOBINDEX}"
    exit 1
fi

echo "[INFO] Processing chunk ${LSB_JOBINDEX}/980: ${CHUNK}"

# Extract chromosome from chunk name
CHR=$(echo ${CHUNK} | cut -d'_' -f1 | sed 's/chr//')
CHUNK_OUT=${OUTDIR}/hepatic_fat_${CHUNK}_results.txt

# Check if output already exists
if [ -s "${CHUNK_OUT}" ]; then
    echo "[SKIP] ${CHUNK_OUT} already exists and is non-empty"
    exit 0
fi

# Run saige step 2 
echo "[START] Running SAIGE step 2 for ${CHUNK}"
singularity exec \
    -B /project/verma_shared:/project \
    -B /static:/static \
    "${SIF}" \
    /app/.pixi/envs/default/bin/Rscript /usr/local/bin/step2_SPAtests.R \
    --bgenFile=${BGEN_PREFIX}.${CHUNK}.bgen \
    --bgenFileIndex=${BGEN_PREFIX}.${CHUNK}.bgen.bgi \
    --sampleFile=/static/PMBB/PMBB-Release-2024-3.0/Imputed/chunked_bgen_files/PMBB-Release-2024-3.0_genetic_imputed.commonsnps.samplelist.txt \
    --AlleleOrder=ref-first \
    --SAIGEOutputFile=${CHUNK_OUT} \
    --chrom=chr${CHR} \
    --minMAF=0 \
    --minMAC=40 \
    --GMMATmodelFile=${OUTDIR}/step1_null_model.rda \
    --varianceRatioFile=${OUTDIR}/step1_null_model.varianceRatio.txt \
    --LOCO=TRUE \
    --is_Firth_beta=TRUE \
    --pCutoffforFirth=0.1 \
    --is_output_moreDetails=TRUE
    
echo "[DONE] Completed chunk ${CHUNK}"

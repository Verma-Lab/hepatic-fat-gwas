#!/bin/bash
set -euo pipefail

SIF=saige_1.5.0.sif
BGEN_PREFIX=/static/PMBB/PMBB-Release-2024-3.0/Imputed/chunked_bgen_files/PMBB-Release-2024-3.0_genetic_imputed
CHUNK_LIST=/project/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_jan2026/test_one_chunk.txt
OUTDIR=/project/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_manual_jan2026/PMBB_ALL_ALL

# Create output directory
mkdir -p ${OUTDIR}


# Step 2: Perform association testing
while read CHUNK; do
    CHR=$(echo ${CHUNK} | cut -d'_' -f1 | sed 's/chr//')
    CHUNK_OUT=${OUTDIR}/hepatic_fat_${CHUNK}_results.txt

    if [ -s "${CHUNK_OUT}" ]; then
	echo "[SKIP] ${CHUNK_OUT} already exists and is non-empty"
        continue
    fi	

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
done < ${CHUNK_LIST}

# Merge chunks per chromosome
# for chr in {1..22}; do
#    files=(${OUTDIR}/hepatic_fat_chr${chr}_*_results.txt)
#    if [ -e "${files[0]}" ]; then
#        awk 'FNR==1 && NR!=1 { next } { print }' "${files[@]}" \
#            > ${OUTDIR}/hepatic_fat_chr${chr}_results.txt
#    fi
# done

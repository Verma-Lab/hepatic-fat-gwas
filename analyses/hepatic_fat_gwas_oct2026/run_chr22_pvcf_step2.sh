#!/bin/bash
# chr22 SAIGE 1.5.1 step 2 on the release pVCF (DS dosages), reusing the SL_transformed step 1 model.
# Flags match the PGEN step 2 run. No variant list (SAIGE's VCF marker IDs didn't match the
# release IDs); the same MAF/MAC filters apply, and runs are matched on variants afterwards.
#
# Usage: bash run_chr22_pvcf_step2.sh <EUR|AFR> <hepatic_fat.rda> <hepatic_fat.varianceRatio.txt>
set -euo pipefail

ANC=$1
RDA=$(readlink -f "$2")
VR=$(readlink -f "$3")

BASE=/project/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_oct2026
SIF=/project/verma_shared/projects/Liver_IDPs/tools/pmbb-nf-toolkit-saige-family_oct2026/saige.sif
VCF_DIR=/static/PMBB/PMBB-Release-2026-4.0/Imputed/chunked_pVCF_files
OUT=$BASE/test_chr22_pvcf_saige151/$ANC
mkdir -p "$OUT/jobs" "$OUT/vcf"

# Variants tested in the PGEN run on chr22 (kept for matching the runs afterwards)
zcat "$BASE/SL_transformed/PMBB_$ANC/Sumstats/PMBB_$ANC.hepatic_fat.gwas.saige.gz" \
  | awk -F'\t' 'NR>1 && ($2=="22" || $2=="chr22") {print $4}' > "$OUT/pgen_chr22_ids.txt"
echo "$ANC: $(wc -l < "$OUT/pgen_chr22_ids.txt") chr22 variants from the PGEN run"
echo "model: $RDA"

# One bash job file per chunk (a bare bsub command string runs in a POSIX shell that
# chokes on the module system's exported functions)
for vcf in "$VCF_DIR"/*_chr22_*.vcf.gz; do
  chunk=$(basename "$vcf" .vcf.gz | sed 's/.*_imputed\.//')
  # SAIGE needs a .csi index; the release only has .tbi, so index a symlink in our folder
  link=$OUT/vcf/$(basename "$vcf")
  ln -sf "$vcf" "$link"
  job=$OUT/jobs/$chunk.sh
  cat > "$job" <<EOF
#!/bin/bash
#BSUB -q epistasis_normal
#BSUB -M 16000
#BSUB -R "rusage[mem=16000]"
#BSUB -J pvcf_${ANC}_$chunk
#BSUB -oo $OUT/$chunk.lsf.log

module load apptainer/1.5.3 bcftools/1.21
set -e
if [ ! -s $link.csi ]; then
  bcftools index -c -f $link
fi
apptainer exec -B /home/,/project/,/static/ $SIF \\
  /usr/local/bin/step2_SPAtests.R \\
    --vcfFile=$link --vcfFileIndex=$link.csi --vcfField=DS \\
    --chrom=chr22 \\
    --GMMATmodelFile=$RDA --varianceRatioFile=$VR \\
    --SAIGEOutputFile=$OUT/PMBB_$ANC.$chunk.txt \\
    --minMAF=0.01 --minMAC=20 --LOCO=TRUE \\
    --is_imputed_data=FALSE --minInfo=0 --is_output_moreDetails=TRUE \\
  > $OUT/$chunk.saige.log 2>&1
EOF
  bsub < "$job"
done



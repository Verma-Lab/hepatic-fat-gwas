#!/bin/bash
# chr22 SAIGE 1.5.1 step 2 on the PGEN chunks with --is_noadjCov=FALSE, reusing the SL_transformed
# step 1 model. Flags otherwise match the PGEN step 2 run; only is_noadjCov changes.
#
# Usage: bash run_chr22_pgen_noadjcov_step2.sh <EUR|AFR> <hepatic_fat.rda> <hepatic_fat.varianceRatio.txt>
set -euo pipefail

ANC=$1
RDA=$(readlink -f "$2")
VR=$(readlink -f "$3")

BASE=/project/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_oct2026
SIF=/project/verma_shared/projects/Liver_IDPs/tools/pmbb-nf-toolkit-saige-family_oct2026/saige.sif
PGEN_DIR=/static/PMBB/PMBB-Release-2026-4.0/Imputed/chunked_pgen_files
OUT=$BASE/test_chr22_pgen_noadjcovF_saige151/$ANC
mkdir -p "$OUT/jobs"
echo "model: $RDA"

# One bash job file per chunk (a bare bsub command string runs in a POSIX shell that
# chokes on the module system's exported functions)
for pgen in "$PGEN_DIR"/*_chr22_*.pgen; do
  prefix=${pgen%.pgen}
  chunk=$(basename "$prefix" | sed 's/.*_imputed\.//')
  # chromosome name exactly as written in the .pvar (22 or chr22), as the pipeline does
  chrom=$(grep -vm1 '^#' "$prefix.pvar" | cut -f1)
  job=$OUT/jobs/$chunk.sh
  cat > "$job" <<EOF
#!/bin/bash
#BSUB -q epistasis_normal
#BSUB -M 16000
#BSUB -R "rusage[mem=16000]"
#BSUB -J noadj_${ANC}_$chunk
#BSUB -oo $OUT/$chunk.lsf.log

module load apptainer/1.5.3
apptainer exec -B /home/,/project/,/static/ $SIF \\
  /usr/local/bin/step2_SPAtests.R \\
    --pgenPrefix=$prefix --AlleleOrder=ref-first \\
    --chrom=$chrom \\
    --GMMATmodelFile=$RDA --varianceRatioFile=$VR \\
    --SAIGEOutputFile=$OUT/PMBB_$ANC.$chunk.txt \\
    --minMAF=0.01 --minMAC=20 --LOCO=TRUE \\
    --is_imputed_data=FALSE --minInfo=0 --is_output_moreDetails=TRUE \\
    --is_noadjCov=FALSE \\
  > $OUT/$chunk.saige.log 2>&1
EOF
  bsub < "$job"
done
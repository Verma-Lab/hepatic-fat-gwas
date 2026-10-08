#!/usr/bin/env bash
# Ancestry- and sex-stratified LDSC heritability for hepatic fat (ΔHU).
# Runs the full pipeline (see README.md for the 6-stage breakdown) for each
# of EUR_F, EUR_M, AFR_F, AFR_M.
#
# Writes everything to WORK_DIR 
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RAW_ROOT="/home/agaro/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_feb2026_sex_stratified"
WORK_DIR="/home/agaro/ldsc_heritability_outputs"

LDSC_TOOLS="/home/agaro/ldsc39"                                          
LDSC39_PYTHON="/home/agaro/mambaforge/envs/ldsc39/bin/python"           # direct path
HM3_SNPLIST="/project/verma_shared/for_Trust/personal/jwoerner/meta/data/alleles/w_hm3.snplist"
LD_PANEL_ROOT="/project/verma_shared/projects/LDSC/AOU_Genetic_Corr_Files/1KG_LD"

mkdir -p "${WORK_DIR}"/{pre_munge,munged,clean,output}

# trait_key -> ancestry (which 1000G LD panel to use)
declare -A ANCESTRY=(
  [EUR_F]=EUR
  [EUR_M]=EUR
  [AFR_F]=AFR
  [AFR_M]=AFR
)

for key in EUR_F EUR_M AFR_F AFR_M; do
  ancestry="${ANCESTRY[$key]}"
  raw_file="${RAW_ROOT}/PMBB_${key}/Sumstats/PMBB_${key}.hepatic_fat.gwas.saige.gz"

  echo "=== ${key} (ancestry=${ancestry}) ==="

  if [[ ! -f "${raw_file}" ]]; then
    echo "[WARN] Missing raw file: ${raw_file} -- skipping ${key}"
    continue
  fi

  # Reformat: attach rsIDs, restrict to HapMap3 (plain python3 + pandas)
  pre_munge="${WORK_DIR}/pre_munge/${key}.pre_munge.tsv.gz"
  python3 "${SCRIPT_DIR}/attach_rsid.py" "${raw_file}" "${pre_munge}"

  # Munge: standardize + QC into ldsc.py's expected shape (needs ldsc39 env)
  munged_prefix="${WORK_DIR}/munged/${key}"
  "${LDSC39_PYTHON}" "${LDSC_TOOLS}/munge_sumstats.py" \
    --sumstats "${pre_munge}" \
    --out "${munged_prefix}" \
    --merge-alleles "${HM3_SNPLIST}" \
    --snp SNP --a1 A1 --a2 A2 --p P \
    --signed-sumstats BETA,0 \
    --N-col N

  # Relabel + clean: SNP -> chr:pos0:A1:A2, drop missing Z (plain python3)
  clean_file="${WORK_DIR}/clean/${key}.clean.sumstats.gz"
  python3 "${SCRIPT_DIR}/relabel_and_clean.py" "${munged_prefix}.sumstats.gz" "${clean_file}"

  # Estimate h2 against the ancestry-matched 1000G LD panel (needs ldsc39 env)
  "${LDSC39_PYTHON}" "${LDSC_TOOLS}/ldsc.py" \
    --h2 "${clean_file}" \
    --ref-ld-chr "${LD_PANEL_ROOT}/${ancestry}/ldscore_gz_L2/subset." \
    --w-ld-chr "${LD_PANEL_ROOT}/${ancestry}/ldscore_gz_L2/subset." \
    --not-M-5-50 \
    --out "${WORK_DIR}/output/h2_${key}"

  echo "=== ${key} done -- see ${WORK_DIR}/output/h2_${key}.log ==="
done

echo "All done. h2 results in ${WORK_DIR}/output/h2_*.log"

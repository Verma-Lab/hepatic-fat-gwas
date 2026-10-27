#!/usr/bin/env bash
# LDSC SNP-heritability, EUR only, for the Oct 2026 GWAS runs (EUR 1000G LD panel).
# Same steps as run_ldsc_heritability.sh, but uses attach_rsid_oct2026.py (checks
# which HapMap3 BED column matches the GWAS positions) and relabel_and_clean_oct2026.py
# (matches SNP IDs to the panel in either allele order).
#
# Usage (from scripts/ldsc/, on a compute node or via bsub):
#   ./run_ldsc_h2_eur_oct2026.sh                                                    # SL_transformed_saige145_bed
#   ./run_ldsc_h2_eur_oct2026.sh SL_transformed_saige145_bed SL_untransformed_saige145_bed
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNS=("${@:-SL_transformed_saige145_bed}")

GWAS_ROOT="/project/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_oct2026"
WORK_DIR="/home/agaro/ldsc_heritability_outputs/oct2026_EUR"     # /home, as in the Feb run (shared volume was full)

LDSC_TOOLS="/home/agaro/ldsc39"
LDSC39_PYTHON="/home/agaro/mambaforge/envs/ldsc39/bin/python"
HM3_SNPLIST="/home/agaro/ldsc_reference/w_hm3.snplist"   # old for_Trust copy was deleted; set to wherever your copy is
LD_PANEL="/project/verma_shared/projects/LDSC/AOU_Genetic_Corr_Files/1KG_LD/EUR/ldscore_gz_L2/subset."

[[ -f "${HM3_SNPLIST}" ]] || { echo "[ERROR] HapMap3 SNP list not found: ${HM3_SNPLIST}"; exit 1; }
mkdir -p "${WORK_DIR}"/{pre_munge,munged,clean,output}

for run in "${RUNS[@]}"; do
  key="${run}_EUR"
  raw_file="${GWAS_ROOT}/${run}/PMBB_EUR/Sumstats/PMBB_EUR.hepatic_fat.gwas.saige.gz"
  echo "=== ${key} ==="
  if [[ ! -f "${raw_file}" ]]; then
    echo "[ERROR] Missing: ${raw_file}"; exit 1
  fi

  # 1. rsIDs + HapMap3 filter (checks the position match on the data; stops if unclear)
  pre_munge="${WORK_DIR}/pre_munge/${key}.pre_munge.tsv.gz"
  python3 "${SCRIPT_DIR}/attach_rsid_oct2026.py" "${raw_file}" "${pre_munge}" "${HM3_SNPLIST}"

  # 2. munge
  munged_prefix="${WORK_DIR}/munged/${key}"
  "${LDSC39_PYTHON}" "${LDSC_TOOLS}/munge_sumstats.py" \
    --sumstats "${pre_munge}" --out "${munged_prefix}" \
    --merge-alleles "${HM3_SNPLIST}" \
    --snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats BETA,0 --N-col N

  # 3. relabel SNP IDs to the LD panel's chr:pos:a1:a2 IDs (either allele order), drop bad Z
  clean_file="${WORK_DIR}/clean/${key}.clean.sumstats.gz"
  python3 "${SCRIPT_DIR}/relabel_and_clean_oct2026.py" "${munged_prefix}.sumstats.gz" "${clean_file}" "${LD_PANEL}"

  # 4. h2 against the EUR 1000G LD scores
  "${LDSC39_PYTHON}" "${LDSC_TOOLS}/ldsc.py" \
    --h2 "${clean_file}" \
    --ref-ld-chr "${LD_PANEL}" --w-ld-chr "${LD_PANEL}" \
    --not-M-5-50 \
    --out "${WORK_DIR}/output/h2_${key}"

  echo "=== ${key} done: ${WORK_DIR}/output/h2_${key}.log ==="
done
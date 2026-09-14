#!/bin/bash
# Step 1 of FUMA prep, adapted for GWAMA meta-analysis output (vs. the original raw-SAIGE
# version documented in LeadSNPs/CLAUDE.md).
#
# GWAMA .out columns:
#   rs_number(chr_pos_ref_alt) reference_allele(=effect allele, confirmed from GWAMA source
#   marker.cpp) other_allele eaf beta se beta_95L beta_95U z p-value -log10_p-value
#   q_statistic q_p-value i2 n_studies n_samples effects
#
# No effect_allele_count column available (GWAMA doesn't report it) -- omitted; FUMA does
# not strictly require it.

set -euo pipefail

for SEX in M F; do
  IN="PMBB_META_EUR_AFR_${SEX}_hepatic_fat.out"
  OUT="PMBB_META_EUR_AFR_${SEX}.hepatic_fat.gwama.fuma.v2.gz"

  awk -F'\t' 'BEGIN{OFS="\t"; print "chromosome","base_pair_location","variant_id","p_value","effect_allele","other_allele","effect_allele_frequency","beta","standard_error"}
  NR>1 {
    n = split($1, parts, "_")
    chr = parts[1]
    sub(/^chr/, "", chr)
    pos = parts[2]
    print chr, pos, $1, $10, $2, $3, $4, $5, $6
  }' "$IN" | gzip > "$OUT"

  echo "Wrote $OUT"
done

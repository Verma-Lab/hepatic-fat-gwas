#!/usr/bin/env python3
"""
liftover_gwama_hg38_to_hg19.py

Liftover hg38 -> hg19 for the FUMA-formatted GWAMA meta-analysis output
(Step 1 output), ahead of the MAF >= 1% filter (Step 3, done separately via
awk/gzcat, same as the baseline pipeline -- see LeadSNPs/CLAUDE.md).

Adapted from the user's original liftover_hg38_to_hg19.py (run locally for the
baseline PMBB_ALL_F/M/ALL_ALL analysis) -- only BASE_DIR, CHAIN, and the file
list differ; loop structure, single shared LiftOver object, and leaving the
chromosome field untouched (only base_pair_location is rewritten) all match
the original verbatim.

Input:  *.gwama.fuma.v2.gz   (9 cols, hg38, tab-separated, header row --
        no effect_allele_count column, since GWAMA doesn't report one)
Output: *.gwama.fuma.hg19.gz (same 9 cols, base_pair_location replaced with hg19 pos)

Run with:
    /project/verma_shared/projects/Liver_IDPs/tools/pyliftover_env/bin/python liftover_gwama_hg38_to_hg19.py
"""

import gzip
from pyliftover import LiftOver

BASE_DIR = "/project/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_feb2026_sex_stratified/Meta"
CHAIN    = "/project/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_feb2026_sex_stratified/LeadSNPs/hg38ToHg19.over.chain.gz"

lo = LiftOver(CHAIN)

for infile, outfile in [
    (f"{BASE_DIR}/PMBB_META_EUR_AFR_M.hepatic_fat.gwama.fuma.v2.gz",
     f"{BASE_DIR}/PMBB_META_EUR_AFR_M.hepatic_fat.gwama.fuma.hg19.gz"),
    (f"{BASE_DIR}/PMBB_META_EUR_AFR_F.hepatic_fat.gwama.fuma.v2.gz",
     f"{BASE_DIR}/PMBB_META_EUR_AFR_F.hepatic_fat.gwama.fuma.hg19.gz"),
]:
    with gzip.open(infile, 'rt') as fin, gzip.open(outfile, 'wt') as fout:
        fout.write(fin.readline())  # header
        failed = converted = 0
        for line in fin:
            fields = line.strip().split('\t')
            result = lo.convert_coordinate('chr' + fields[0], int(fields[1]) - 1)
            if result:
                fields[1] = str(result[0][1] + 1)
                fout.write('\t'.join(fields) + '\n')
                converted += 1
            else:
                failed += 1
        print(f"{outfile}: Converted {converted:,}, Failed {failed:,}")

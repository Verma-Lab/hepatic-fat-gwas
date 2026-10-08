"""
LDSC (Reformat): attach HapMap3 rsIDs to a raw SAIGE GWAS file,
keeping only SNPs present in the HapMap3 reference set.

Match logic: matches directly on base_pair_location == hm3 start0, since this
raw SAIGE output's position is already 0-based. Matches on chr+position only,
allele QC happens later in munge_sumstats.py.

Input:  raw SAIGE gwas.gz (15 cols, hg38, chr-prefixed chromosome column)
Output: {out}  with columns SNP(rsid) A1 A2 BETA P N, the input munge_sumstats.py 
        expects via --snp SNP --a1 A1 --a2 A2 --p P
        --signed-sumstats BETA,0 --N-col N

Run with:
    python attach_rsid.py <raw_gwas.gz> <out.pre_munge.tsv.gz>
"""
import sys
import pandas as pd

HM3_BED = "/project/pmbb_codeworks/datasets/Hapmap_SNPs/BED_Files/hapmap_phase3_snps_b38.bed"


def main(infile, outfile):
    hm3 = pd.read_csv(HM3_BED, sep='\t', header=None,
                       names=['chr', 'start0', 'end1', 'rsid', 'hm3_a1', 'hm3_a2'],
                       dtype={'chr': str})
    hm3['chr'] = hm3['chr'].str.replace('^chr', '', regex=True)
    # Force int64 explicitly 
    hm3['start0'] = pd.to_numeric(hm3['start0'], errors='coerce').astype('Int64')
    hm3 = hm3.dropna(subset=['start0'])
    hm3['start0'] = hm3['start0'].astype('int64')
    hm3 = hm3.drop_duplicates(subset=['chr', 'start0'])

    gwas = pd.read_csv(infile, sep='\t', compression='gzip',
                        usecols=['chromosome', 'base_pair_location', 'effect_allele',
                                 'other_allele', 'beta', 'p_value', 'N'],
                        dtype={'chromosome': str})
    gwas['chr'] = gwas['chromosome'].str.replace('^chr', '', regex=True)
    gwas['start0'] = gwas['base_pair_location'].astype('int64')

    merged = gwas.merge(hm3[['chr', 'start0', 'rsid']], on=['chr', 'start0'], how='inner')

    out = merged.rename(columns={
        'rsid': 'SNP', 'effect_allele': 'A1', 'other_allele': 'A2',
        'beta': 'BETA', 'p_value': 'P',
    })[['SNP', 'A1', 'A2', 'BETA', 'P', 'N']]

    out.to_csv(outfile, sep='\t', index=False, compression='gzip')
    print(f"{infile}: {len(gwas):,} input SNPs -> {len(out):,} matched to HapMap3 -> {outfile}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <raw_gwas.gz> <out.pre_munge.tsv.gz>")
    main(sys.argv[1], sys.argv[2])

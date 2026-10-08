"""
Relabel + clean: rewrite munge_sumstats.py's rsID-keyed SNP column into 
chr:pos0:A1:A2 convention, to match the SNP IDs used by the 1KG_LD
reference LD-score panel, then drop rows with missing/invalid Z, which
can't contribute to the h2 regression. Run this after munge_sumstats.py and
before ldsc.py --h2 (see README.md for the full pipeline order).

Same rsid->position lookup source as attach_rsid.py (the HapMap3 BED file),
just used in the reverse direction: rsID -> (chr, position) instead of
(chr, position) -> rsID.

Input:  {trait}.sumstats.gz  (munge_sumstats.py output: SNP[rsid] A1 A2 N Z)
Output: {trait}.clean.sumstats.gz (same columns, SNP relabeled, bad Z dropped)

Run with:
    python relabel_and_clean.py <in.sumstats.gz> <out.clean.sumstats.gz>
"""
import sys
import pandas as pd

HM3_BED = "/project/pmbb_codeworks/datasets/Hapmap_SNPs/BED_Files/hapmap_phase3_snps_b38.bed"


def main(infile, outfile):
    hm3 = pd.read_csv(HM3_BED, sep='\t', header=None,
                       names=['chr', 'start0', 'end1', 'rsid', 'hm3_a1', 'hm3_a2'],
                       dtype={'chr': str})
    hm3['chr'] = hm3['chr'].str.replace('^chr', '', regex=True)
    hm3 = hm3.drop_duplicates(subset=['rsid'])

    ss = pd.read_csv(infile, sep='\t', compression='gzip')
    n_before = len(ss)

    merged = ss.merge(hm3[['rsid', 'chr', 'start0']], left_on='SNP', right_on='rsid', how='left')
    merged['SNP'] = ('chr' + merged['chr'] + ':' + merged['start0'].astype('Int64').astype(str)
                      + ':' + merged['A1'] + ':' + merged['A2'])

    merged = merged.dropna(subset=['SNP', 'Z'])
    out = merged[['SNP', 'A1', 'A2', 'N', 'Z']]

    out.to_csv(outfile, sep='\t', index=False, compression='gzip')
    print(f"{infile}: {n_before:,} -> {len(out):,} after relabel+clean -> {outfile}")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(f"Usage: {sys.argv[0]} <in.sumstats.gz> <out.clean.sumstats.gz>")
    main(sys.argv[1], sys.argv[2])

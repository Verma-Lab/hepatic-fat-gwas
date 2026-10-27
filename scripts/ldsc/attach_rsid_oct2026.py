#!/usr/bin/env python3
"""
Attach HapMap3 rsIDs to an Oct 2026 SAIGE GWAS file (SAIGE 1.5.1, PGEN,
GRCh38) and keep only SNPs in w_hm3.snplist; output feeds munge_sumstats.py.

The GWAS positions are 1-based (FUMA's dbSNP lookup named the SNPs at our
positions rs738408/rs1883350). The HapMap3 BED's 2nd column turned out to hold
that same 1-based position (rs1883350: BED 43932163, GWAS 43932163), so it is
not standard 0-based BED. Rather than assume, this script matches GWAS
positions to both BED columns, compares allele agreement with HapMap3, uses the
column that clearly wins, and stops if neither does.

Output columns: SNP(rsid) A1 A2 BETA P N, for munge_sumstats.py with
--snp SNP --a1 A1 --a2 A2 --p P --signed-sumstats BETA,0 --N-col N

Run with:
    python3 attach_rsid_oct2026.py <raw_gwas.gz> <out.pre_munge.tsv.gz> [w_hm3.snplist]
"""
import sys
import pandas as pd

HM3_BED = "/project/pmbb_codeworks/datasets/Hapmap_SNPs/BED_Files/hapmap_phase3_snps_b38.bed"
HM3_SNPLIST = "/project/verma_shared/for_Trust/personal/jwoerner/meta/data/alleles/w_hm3.snplist"
MIN_ALLELE_AGREEMENT = 0.90      # among position matches, share whose alleles agree with HapMap3
COMPLEMENT = str.maketrans("ACGT", "TGCA")
KNOWN_SNP = ("22", 43932163, "rs1883350")   # GWAS position; FUMA mapped this position to rs1883350


def allele_agreement(m):
    """Share of matched rows whose {A1, A2} equal the HapMap3 alleles (either order or strand)."""
    gwas = m["effect_allele"].str.upper() + "/" + m["other_allele"].str.upper()
    hm3 = m["hm3_a1"].str.upper() + "/" + m["hm3_a2"].str.upper()
    hm3_rev = m["hm3_a2"].str.upper() + "/" + m["hm3_a1"].str.upper()
    ok = ((gwas == hm3) | (gwas == hm3_rev)
          | (gwas == hm3.str.translate(COMPLEMENT)) | (gwas == hm3_rev.str.translate(COMPLEMENT)))
    return ok.mean() if len(m) else 0.0


def main(infile, outfile, snplist_path=HM3_SNPLIST):
    hm3 = pd.read_csv(HM3_BED, sep="\t", header=None,
                      names=["chr", "col2", "col3", "rsid", "hm3_a1", "hm3_a2"], dtype=str)
    hm3["chr"] = hm3["chr"].str.replace("^chr", "", regex=True)
    for col in ["col2", "col3"]:
        hm3[col] = pd.to_numeric(hm3[col], errors="coerce")
    hm3 = hm3.dropna(subset=["col2", "col3"])
    hm3[["col2", "col3"]] = hm3[["col2", "col3"]].astype("int64")
    print(f"HapMap3 BED: {len(hm3):,} rows")

    gwas = pd.read_csv(infile, sep="\t", compression="gzip", dtype={"chromosome": str},
                       usecols=["chromosome", "base_pair_location", "effect_allele",
                                "other_allele", "beta", "p_value", "N"])
    gwas["chr"] = gwas["chromosome"].str.replace("^chr", "", regex=True)
    gwas["pos"] = gwas["base_pair_location"].astype("int64")
    print(f"GWAS: {len(gwas):,} variants")

    # match on both BED position columns; let the alleles decide which is right
    results = {}
    for col in ["col2", "col3"]:
        ref = hm3.drop_duplicates(subset=["chr", col])[["chr", col, "rsid", "hm3_a1", "hm3_a2"]]
        m = gwas.merge(ref.rename(columns={col: "pos"}), on=["chr", "pos"], how="inner")
        results[col] = (m, allele_agreement(m))
        print(f"  GWAS position == BED {col}: {len(m):,} matches, alleles agree with HapMap3 for {results[col][1]:.1%}")

    best, other = sorted(results, key=lambda c: results[c][1], reverse=True)
    if results[best][1] < MIN_ALLELE_AGREEMENT or results[best][1] - results[other][1] < 0.2:
        sys.exit("Coordinate check FAILED: neither BED column clearly matches the GWAS positions.\n"
                 "Check the positions before running LDSC (see the numbers above).")
    print(f"Coordinate check: matching GWAS positions to BED {best}.")
    m = results[best][0]

    hit = m[(m["chr"] == KNOWN_SNP[0]) & (m["pos"] == KNOWN_SNP[1])]
    print(f"Known-SNP check: {KNOWN_SNP[0]}:{KNOWN_SNP[1]} -> "
          f"{hit['rsid'].tolist() if len(hit) else 'not matched'} (expect {KNOWN_SNP[2]})")

    out = m.rename(columns={"rsid": "SNP", "effect_allele": "A1", "other_allele": "A2",
                            "beta": "BETA", "p_value": "P"})[["SNP", "A1", "A2", "BETA", "P", "N"]]
    n_matched = len(out)
    snplist = set(pd.read_csv(snplist_path, sep=r"\s+", usecols=["SNP"])["SNP"])
    out = out[out["SNP"].isin(snplist)]
    out.to_csv(outfile, sep="\t", index=False, compression="gzip")
    print(f"{infile}: {len(gwas):,} variants -> {n_matched:,} matched to HapMap3 "
          f"-> {len(out):,} also in w_hm3.snplist -> {outfile}")


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        sys.exit(f"Usage: {sys.argv[0]} <raw_gwas.gz> <out.pre_munge.tsv.gz> [w_hm3.snplist]")
    main(*sys.argv[1:])
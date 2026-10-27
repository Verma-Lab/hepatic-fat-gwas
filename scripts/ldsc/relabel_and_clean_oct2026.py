#!/usr/bin/env python3
"""
Relabel munged LDSC sumstats from rsID to the LD panel's chr:pos:a1:a2 IDs and
drop rows with missing Z. Unlike relabel_and_clean.py, the allele order is
taken from the panel itself: each SNP gets whichever order (A1:A2 or A2:A1)
exists in the panel, so SNPs aren't lost to an order mismatch.

Run with:
    python3 relabel_and_clean_oct2026.py <in.sumstats.gz> <out.clean.sumstats.gz> <ld_panel_prefix>
    (ld_panel_prefix e.g. .../1KG_LD/EUR/ldscore_gz_L2/subset. ; files <prefix><chr>.l2.ldscore.gz)
"""
import sys

import pandas as pd

HM3_BED = "/project/pmbb_codeworks/datasets/Hapmap_SNPs/BED_Files/hapmap_phase3_snps_b38.bed"


def main(infile, outfile, panel_prefix):
    # BED column 2 holds the 1-based GRCh38 position (checked 2026-10-14), as in the panel IDs
    hm3 = pd.read_csv(HM3_BED, sep="\t", header=None, usecols=[0, 1, 3],
                      names=["chr", "pos", "rsid"], dtype=str)
    hm3["chr"] = hm3["chr"].str.replace("^chr", "", regex=True)
    hm3 = hm3.drop_duplicates(subset=["rsid"])

    panel_ids = set()
    for c in range(1, 23):
        panel_ids.update(pd.read_csv(f"{panel_prefix}{c}.l2.ldscore.gz", sep="\t", usecols=["SNP"])["SNP"])
    print(f"LD panel: {len(panel_ids):,} SNP IDs")

    ss = pd.read_csv(infile, sep="\t", compression="gzip")
    n_in = len(ss)
    ss = ss.dropna(subset=["Z"])
    m = ss.merge(hm3, left_on="SNP", right_on="rsid", how="inner")

    base = "chr" + m["chr"] + ":" + m["pos"] + ":"
    id_fwd = base + m["A1"] + ":" + m["A2"]
    id_rev = base + m["A2"] + ":" + m["A1"]
    fwd, rev = id_fwd.isin(panel_ids), id_rev.isin(panel_ids)
    m["SNP"] = id_fwd.where(fwd, id_rev.where(rev, id_fwd))
    print(f"Matched to panel: as A1:A2 {fwd.sum():,} | as A2:A1 {(rev & ~fwd).sum():,} | "
          f"not in panel {(~fwd & ~rev).sum():,}")

    out = m[["SNP", "A1", "A2", "N", "Z"]]
    out.to_csv(outfile, sep="\t", index=False, compression="gzip")
    print(f"{infile}: {n_in:,} rows -> {len(ss):,} with Z -> {len(out):,} relabeled -> {outfile}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(f"Usage: {sys.argv[0]} <in.sumstats.gz> <out.clean.sumstats.gz> <ld_panel_prefix>")
    main(*sys.argv[1:])

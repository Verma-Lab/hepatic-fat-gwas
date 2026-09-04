#!/usr/bin/env python3
"""
make_male_top10_leadSNPs.py

Generate PMBB_male_top10_leadSNPS_FUMA.tsv from FUMA output files.

Required inputs (adjust paths below):
  - leadSNPs.txt  : FUMA lead SNPs table
  - snps.txt      : FUMA all SNPs table (used for annotation of lead SNPs)
  - eqtl.txt      : FUMA eQTL mapping table

Output:
  - PMBB_male_top10_leadSNPS_FUMA.tsv
"""

import pandas as pd
import os

# ── paths ──────────────────────────────────────────────────────────────────────
FUMA_DIR   = "/Users/agaro/Documents/meno_effects_liver/fuma_results/fuma_male"
OUT_FILE   = "/Users/agaro/Documents/meno_effects_liver/fuma_results/PMBB_male_top10_leadSNPS_FUMA.tsv"

LEAD_SNPS  = os.path.join(FUMA_DIR, "leadSNPs.txt")
SNPS       = os.path.join(FUMA_DIR, "snps.txt")
EQTL       = os.path.join(FUMA_DIR, "eqtl.txt")
# ───────────────────────────────────────────────────────────────────────────────

# 1. Load lead SNPs and select top 10 by p-value
lead = pd.read_csv(LEAD_SNPS, sep="\t")
lead = lead.sort_values("p").head(10).reset_index(drop=True)
lead["rank"] = lead.index + 1

# 2. Load snps.txt — keep only the lead SNP rows (r2 == 1 with themselves)
#    matched on uniqID so we pull beta, SE, MAF, annotation columns
snps = pd.read_csv(SNPS, sep="\t")
lead_annot = snps[snps["uniqID"].isin(lead["uniqID"])][
    ["uniqID", "MAF", "beta", "se", "nearestGene", "dist", "func"]
].copy()

# 3. Load eqtl.txt — collect unique mapped gene symbols per GenomicLocus
#    Join eqtl SNPs to snps.txt to get their GenomicLocus, then group
eqtl = pd.read_csv(EQTL, sep="\t")
eqtl_locus = (
    eqtl[eqtl["eqtlMapFilt"] == 1]
    .merge(snps[["uniqID", "GenomicLocus"]], on="uniqID", how="left")
    .dropna(subset=["GenomicLocus"])
)
eqtl_genes = (
    eqtl_locus.groupby("GenomicLocus")["symbol"]
    .apply(lambda x: ";".join(sorted(set(x))))
    .reset_index()
    .rename(columns={"symbol": "eQTL_gene"})
)
eqtl_genes["GenomicLocus"] = eqtl_genes["GenomicLocus"].astype(int)

# 4. Merge everything together
out = (
    lead[["rank", "GenomicLocus", "rsID", "chr", "pos", "p"]]
    .merge(lead_annot, on="uniqID", how="left")
    .merge(eqtl_genes, on="GenomicLocus", how="left")
)
out["eQTL_gene"] = out["eQTL_gene"].fillna("—")

# 5. Round numeric columns and rename for clarity
out = out.rename(columns={
    "GenomicLocus": "locus",
    "p":            "p_value",
    "se":           "SE",
    "nearestGene":  "nearest_gene",
    "dist":         "dist_bp",
})
out["p_value"] = out["p_value"].map(lambda x: f"{x:.2e}")
out["beta"]    = out["beta"].round(4)
out["SE"]      = out["SE"].round(5)
out["MAF"]     = out["MAF"].round(4)

# Final column order
out = out[[
    "rank", "locus", "rsID", "chr", "pos",
    "p_value", "beta", "SE", "MAF",
    "nearest_gene", "dist_bp", "func", "eQTL_gene"
]]

out.to_csv(OUT_FILE, sep="\t", index=False)
print(f"Written: {OUT_FILE}")
print(out.to_string(index=False))

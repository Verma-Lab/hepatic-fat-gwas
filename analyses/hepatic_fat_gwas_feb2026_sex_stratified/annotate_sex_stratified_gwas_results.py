"""
annotate_sex_stratified_gwas_results.py

Filters and annotates sex-stratified SAIGE GWAS summary statistics with gene names
and rsIDs from a biofilter annotation file, then exports results for each sex separately.

Usage:
    python annotate_sex_stratified_gwas_results.py \
        --results_dir_m <path/to/male/sumstats> \
        --results_dir_f <path/to/female/sumstats> \
        --annot_dir <path/to/annotations> \
        [--prefix_m PMBB_ALL_M] [--prefix_f PMBB_ALL_F] \
        [--annot_file GWAS_biofilter_genes_rsids.csv] [--pval_threshold 5e-8]

Output:
    <results_dir_m>/<prefix_m>.hepatic_fat.gwas.annotated.txt
    <results_dir_f>/<prefix_f>.hepatic_fat.gwas.annotated.txt
"""

import pandas as pd
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser(description="Annotate sex-stratified GWAS significant hits")
parser.add_argument("--results_dir_m", required=True, help="Directory containing male GWAS summary stats")
parser.add_argument("--results_dir_f", required=True, help="Directory containing female GWAS summary stats")
parser.add_argument("--annot_dir",     required=True, help="Directory containing biofilter annotation file")
parser.add_argument("--prefix_m", default="PMBB_ALL_M", help="File prefix for male summary stats (default: PMBB_ALL_M)")
parser.add_argument("--prefix_f", default="PMBB_ALL_F", help="File prefix for female summary stats (default: PMBB_ALL_F)")
parser.add_argument("--annot_file",    default="GWAS_biofilter_genes_rsids.csv", help="Annotation filename (default: GWAS_biofilter_genes_rsids.csv)")
parser.add_argument("--pval_threshold",type=float, default=5e-8, help="P-value threshold (default: 5e-8)")
args = parser.parse_args()

results_file_m = os.path.join(args.results_dir_m, f"{args.prefix_m}.hepatic_fat.gwas.saige.gz")
results_file_f = os.path.join(args.results_dir_f, f"{args.prefix_f}.hepatic_fat.gwas.saige.gz")
annot_file     = os.path.join(args.annot_dir, args.annot_file)

# Load annotations once
print("Loading gene annotations...")
annot = pd.read_csv(annot_file)

def load_and_annotate(results_file, label, results_dir, prefix):
    print(f"\nLoading GWAS results ({label})...")
    df = pd.read_csv(results_file, sep='\t', compression='gzip')
    print(f"Total variants: {len(df):,}")
    print(f"Columns: {list(df.columns)}")

    df = df.merge(annot[['Var_ID', 'Gene', 'RSID']],
                  left_on='variant_id',
                  right_on='Var_ID',
                  how='left')
    print(f"Variants with gene annotations: {df['Gene'].notna().sum():,}")

    df_sorted = df.sort_values('p_value')
    sig_hits = df_sorted[df_sorted['p_value'] < args.pval_threshold]
    print(f"Hits passing threshold (p < {args.pval_threshold}): {len(sig_hits)}")

    if len(sig_hits) > 0:
        print(sig_hits[['chromosome', 'base_pair_location', 'variant_id',
                         'Gene', 'RSID', 'effect_allele', 'beta', 'p_value']].to_string())

    output_file = os.path.join(results_dir, f"{prefix}.hepatic_fat.gwas.annotated.txt")
    sig_hits.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved {len(sig_hits)} hits to: {output_file}")

load_and_annotate(results_file_m, "Male",   args.results_dir_m, args.prefix_m)
load_and_annotate(results_file_f, "Female", args.results_dir_f, args.prefix_f)
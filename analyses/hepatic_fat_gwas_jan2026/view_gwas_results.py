import pandas as pd
import numpy as np
import os

results_dir = "/home/agaro/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_manual_jan2026/PMBB_ALL_ALL/Sumstats"
annot_dir = "/home/agaro/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_jan2026/Annotations"

results_file = os.path.join(results_dir, "PMBB_ALL_ALL.hepatic_fat.gwas.saige.gz")

# Read the gzipped file
print("Loading GWAS results...")
df = pd.read_csv(results_file, sep='\t', compression='gzip')
annot_file = os.path.join(annot_dir, "GWAS_biofilter_genes_rsids.csv")

print(f"Total variants: {len(df):,}")
print(f"\nColumns: {list(df.columns)}")

# Load and merge the annotations
print("\nLoading gene annotations...")
annot = pd.read_csv(annot_file)
df = df.merge(annot[['Var_ID', 'Gene', 'RSID']], 
              left_on='variant_id', 
              right_on='Var_ID', 
              how='left')

print(f"Variants with gene annotations: {df['Gene'].notna().sum():,}")


# Calculate Odds Ratio (OR = e^beta)
df['OR'] = np.exp(df['beta'])

# Sort by p_value
df_sorted = df.sort_values('p_value')

# View genome-wide significant hits (p < 5e-8)
sig_hits = df_sorted[df_sorted['p_value'] < 5e-8]
print(f"\n\nGenome-wide significant hits (p < 5e-8): {len(sig_hits)}")
if len(sig_hits) > 0:
    print(sig_hits[['chromosome', 'base_pair_location', 'variant_id',
                    'Gene', 'RSID','effect_allele', 'OR', 'beta', 'p_value']].to_string())
output_file = os.path.join(results_dir, "PMBB_ALL_ALL.hepatic_fat.gwas.significant_annotated.txt")

sig_hits.to_csv(output_file, sep = '\t', index = False)

print(f"\n\nSaved {len(sig_hits)} significant hits with annotations to: {output_file}")


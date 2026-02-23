import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from manhattan_plot import ManhattanPlot

# ── Argument parsing ────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description="Generate Manhattan and QQ plots from GWAS summary statistics.")

parser.add_argument("--cohort",      required=True,  help="Cohort name (e.g. PMBB_ALL_ALL)")
parser.add_argument("--pheno",       required=True,  help="Phenotype name (e.g. hepatic_fat)")
parser.add_argument("--results_dir", required=True,  help="Directory containing GWAS summary stats .gz file")
parser.add_argument("--annot_dir",   required=False, default=None, help="Directory containing annotation CSV file")
parser.add_argument("--output_dir",  required=True,  help="Directory to save plots")
parser.add_argument("--pheno_table", required=True,  help="Full path to pheno_summaries.csv")
parser.add_argument("--colnames",    required=True,  help="Full path to colnames.txt")

args = parser.parse_args()

# ── Assign variables ─────────────────────────────────────────────────────────
cohort      = args.cohort
pheno       = args.pheno
results_dir = args.results_dir
annot_dir   = args.annot_dir
output_dir  = args.output_dir

sumstats_file   = os.path.join(results_dir, f"{cohort}.{pheno}.gwas.saige.gz")
annot_file      = os.path.join(annot_dir, "GWAS_biofilter_genes_rsids.csv") if annot_dir else None
pheno_table_file = args.pheno_table
colnames_file   = args.colnames

# Load column mapping
print("Loading column mappings...")
columns_map = {}
with open(colnames_file, 'r') as file:
    for line in file:
        key, value = line.split('=')
        key = key.strip()
        value = value.strip().strip("'")
        columns_map[key] = value

# Reverse the columns map
columns_map_inv = {v: k for k, v in columns_map.items()}
print(f"Column mappings loaded: {columns_map}")

# Load phenotype data 
print("\nLoading phenotype information...")
pheno_df = pd.read_csv(pheno_table_file)
pheno_df = pheno_df[pheno_df['PHENO'] == pheno]
pheno_df = pheno_df[pheno_df['COHORT'] == cohort]


# Determine trait type 
trait_type = 'bin' if pheno_df['Cases'].count() != 0 else 'quant'
print(f"Trait type: {trait_type}")

# Set up output path
output_manhattan = os.path.join(output_dir, f'{cohort}.{pheno}.manhattan_vertical.png')
output_qq = os.path.join(output_dir, f'{cohort}.{pheno}.qq.png')

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Create plot title
plot_title = f'GWAS {cohort}: {pheno.replace("_", " ")}'

if trait_type == 'bin':
    plot_title += f'\nCases = {pheno_df.iloc[0].Cases:,.0f}, Controls = {pheno_df.iloc[0].Controls:,.0f}'
else:
    plot_title += f'\nN = {pheno_df.iloc[0].N:,.0f}'

print(f"\nPlot title: {plot_title}")

# Load and process data
print("\nInitializing ManhattanPlot...")
mp = ManhattanPlot(sumstats_file, test_rows=None, title=plot_title)
mp.load_data()

print("Processing chromosome and variant IDs...")
mp.df['chromosome_noCHR'] = mp.df['chromosome'].astype(str).str.replace('chr', '').astype(int)
mp.df['variant_id'] = mp.df['variant_id'].str.replace('_', ':')

print("Cleaning data with column mappings...")
mp.clean_data(col_map={
    'chromosome_noCHR': '#CHROM', 
    columns_map['POS']: 'POS', 
    columns_map['MarkerID']: 'ID', 
    columns_map['p.value']: 'P'
})

print(f"Data shape: {mp.df.shape}")
print(mp.df.head())

# Annotations
if annot_file is not None and os.path.exists(annot_file):
    print("\nLoading and adding annotations...")
    annot_df = pd.read_csv(annot_file)
    annot_df['ID'] = annot_df['Gene']
    annot_df.iloc[:, 0] = annot_df.iloc[:, 0].str.replace('_', ':')
    mp.add_annotations(annot_df, extra_cols=['RSID'])
    print(f"Annotations added: {len(annot_df)} entries")
else:
    print("\nNo annotation file found, skipping annotations...")


# Thin data 
print("\nThinning data for plotting...")
mp.get_thinned_data()
print(f"Thinned data points: {len(mp.thinned)}")
print(mp.thinned.head())

# Set annotation threshold
annot_thresh = 5E-8 if np.any(mp.thinned['P'].min() < 5E-8) else np.nanquantile(mp.thinned['P'], 10 / len(mp.thinned))
print(f"\nAnnotation threshold: {annot_thresh:.2e}")

# Update plotting parameters
print("Setting plotting parameters...")
mp.update_plotting_parameters(
    vertical=True, 
    sig=annot_thresh if not np.any(mp.thinned['P'] < 5E-8) else 5E-8, 
    sug=annot_thresh, 
    annot_thresh=annot_thresh, 
    merge_genes=True
)

# Create Manhattan plot
print("\nGenerating Manhattan plot...")
mp.full_plot(
    save=output_manhattan, 
    rep_boost=False, 
    extra_cols={'beta': 'BETA', 'RSID': 'RSID'}, 
    number_cols=['BETA'], 
    keep_chr_pos=False, 
    with_table_bg=False, 
    with_table_grid=False
)

plt.clf()
print(f"✓ Manhattan plot saved to: {output_manhattan}")

# Create QQ Plot 
print("\nGenerating QQ plot...")
mp.qq_plot(save=output_qq)
print(f"✓ QQ plot saved to: {output_qq}")

# Create plots manifest
print("\nCreating plots manifest...")

print("\n" + "="*50)
print("DONE! All plots generated successfully.")
print("="*50)
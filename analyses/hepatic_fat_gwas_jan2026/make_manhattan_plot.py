import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Paths
results_dir = "/home/agaro/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_manual_jan2026/PMBB_ALL_ALL/Sumstats"
annot_dir = "/home/agaro/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_jan2026/Annotations"
plots_dir = "/home/agaro/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_manual_jan2026/Plots"

results_file = os.path.join(results_dir, "PMBB_ALL_ALL.hepatic_fat.gwas.saige.gz")
annot_file = os.path.join(annot_dir, "GWAS_biofilter_genes_rsids.csv")

# Read data
print("Loading GWAS results...")
df = pd.read_csv(results_file, sep='\t', compression='gzip')

# Load annotations
print("Loading gene annotations...")
annot = pd.read_csv(annot_file)
df = df.merge(annot[['Var_ID', 'Gene', 'RSID']], 
              left_on='variant_id', 
              right_on='Var_ID', 
              how='left')

# Clean data
df = df[df['p_value'].notna()].copy()
df['neg_log_p'] = -np.log10(df['p_value'])

# Convert chromosome to numeric
df['chr_num'] = df['chromosome'].astype(str).str.replace('chr', '')
df = df[df['chr_num'].str.isnumeric()]
df['chr_num'] = df['chr_num'].astype(int)
df = df.sort_values(['chr_num', 'base_pair_location'])

print(f"Total variants before thinning: {len(df):,}")

# Thin data to reduce overplotting - keep all significant + sample of non-significant
def thin_data(df, sig_threshold=5e-8, target_points=50000):
    """Keep all significant variants and thin non-significant ones"""
    sig = df[df['p_value'] < sig_threshold].copy()
    nonsig = df[df['p_value'] >= sig_threshold].copy()
    
    # Randomly sample non-significant variants
    if len(nonsig) > target_points:
        nonsig = nonsig.sample(n=target_points, random_state=42)
    
    thinned = pd.concat([sig, nonsig]).sort_values(['chr_num', 'base_pair_location'])
    return thinned

df_thinned = thin_data(df, sig_threshold=5e-8, target_points=50000)
print(f"Variants after thinning: {len(df_thinned):,}")

# Calculate cumulative position
df_thinned['cumulative_pos'] = 0
chr_centers = []

for chrom in range(1, 23):
    chr_mask = df_thinned['chr_num'] == chrom
    if chr_mask.sum() > 0:
        if chrom == 1:
            chr_start = 0
        else:
            chr_start = df_thinned[df_thinned['chr_num'] < chrom]['cumulative_pos'].max()
        
        df_thinned.loc[chr_mask, 'cumulative_pos'] = chr_start + df_thinned.loc[chr_mask, 'base_pair_location']
        chr_centers.append(df_thinned[chr_mask]['cumulative_pos'].median())

# Create plot
fig, ax = plt.subplots(figsize=(18, 6))

# Plot by chromosome with alternating colors
colors = ['#4169E1', '#FFA500']  # Blue and Orange
for i, chrom in enumerate(range(1, 23)):
    chr_df = df_thinned[df_thinned['chr_num'] == chrom]
    if len(chr_df) > 0:
        ax.scatter(chr_df['cumulative_pos'], chr_df['neg_log_p'], 
                  c=colors[i % 2], s=10, alpha=0.7, edgecolors='none')

# Significance lines
ax.axhline(y=-np.log10(5e-8), color='red', linestyle='--', linewidth=1.5, 
          label='Genome-wide sig. (p=5e-8)', zorder=3)
ax.axhline(y=-np.log10(1e-5), color='gray', linestyle='--', linewidth=1, 
          alpha=0.5, label='Suggestive (p=1e-5)', zorder=3)

# Label top hits
sig_hits = df_thinned[df_thinned['p_value'] < 5e-8].sort_values('p_value')
labeled_genes = set()

for idx, row in sig_hits.head(20).iterrows():
    if pd.notna(row['Gene']) and row['Gene'] not in labeled_genes:
        ax.annotate(row['Gene'], 
                   xy=(row['cumulative_pos'], row['neg_log_p']),
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=9, fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                   arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
        labeled_genes.add(row['Gene'])

# Formatting
ax.set_xlabel('Chromosome', fontsize=14, fontweight='bold')
ax.set_ylabel('-log₁₀(p-value)', fontsize=14, fontweight='bold')
ax.set_title('GWAS Manhattan Plot: Hepatic Fat\nN = 18,079', 
            fontsize=16, fontweight='bold', pad=20)

# X-axis
ax.set_xticks(chr_centers)
ax.set_xticklabels(range(1, 23), fontsize=11)
ax.set_xlim(df_thinned['cumulative_pos'].min(), df_thinned['cumulative_pos'].max())

# Y-axis
ax.set_ylim(0, df_thinned['neg_log_p'].max() * 1.05)
ax.grid(True, alpha=0.2, axis='y', linestyle=':')

# Legend
ax.legend(loc='upper right', fontsize=10, framealpha=0.9)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

# Save
output_file = os.path.join(plots_dir, 'hepatic_fat_manhattan_plot.png')
plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
print(f"\nManhattan plot saved to: {output_file}")

# Print summary
print(f"\nGenome-wide significant hits (p < 5e-8): {len(sig_hits)}")
print("\nTop 10 significant loci:")
print(sig_hits[['chr_num', 'base_pair_location', 'Gene', 'RSID', 'p_value', 'beta']].head(10).to_string())

plt.close()
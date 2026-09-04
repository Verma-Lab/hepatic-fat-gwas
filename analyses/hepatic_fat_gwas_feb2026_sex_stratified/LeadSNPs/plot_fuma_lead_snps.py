#!/usr/bin/env python3
"""
Mirrored forest plot — FUMA lead SNPs (sex-stratified hepatic fat GWAS)

Method:
  1. Load lead SNPs from male and female FUMA runs
  2. Union positions; deduplicate cross-sex near-duplicates using 250 kb window
     (FUMA already defined independent loci via LD clumping — this only merges
     cases where male and female runs picked slightly different positions for the
     same locus, e.g. PNPLA3 at 44,324,727 vs 44,324,855)
  3. Pick representative SNP per locus (min P across either sex)
  4. Look up cross-sex stats from raw GWAS files
  5. Compute het Z-test: Z = (β_M - β_F) / sqrt(SE_M² + SE_F²)
  6. Classify: shared (P < 1e-5 in both, same direction) vs sex-specific
  7. Asterisk (*) on label = p_het < 0.05
  8. Display top 10 loci by min P

Gene labels: taken directly from FUMA nearest_gene column (positional mapping).
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Paths ──────────────────────────────────────────────────────────────────────
MALE_GWAS    = '/Users/agaro/Documents/meno_effects_liver/gwas_processed/PMBB_ALL_M.hepatic_fat.gwas.saige.fuma.maf01.hg19.gz'
FEMALE_GWAS  = '/Users/agaro/Documents/meno_effects_liver/gwas_processed/PMBB_ALL_F.hepatic_fat.gwas.saige.fuma.maf01.hg19.gz'
MALE_LEADS   = '/Users/agaro/Documents/meno_effects_liver/fuma_results/PMBB_male_top10_leadSNPS_FUMA.tsv'
FEMALE_LEADS = '/Users/agaro/Documents/meno_effects_liver/fuma_results/PMBB_female_top10_loci_FUMA.tsv'
OUT_PNG      = '/Users/agaro/Documents/meno_effects_liver/sex_specific_loci_fuma_leads.png'

# ── Parameters ─────────────────────────────────────────────────────────────────
P_THRESH   = 1e-5
GWS_THRESH = 5e-8
MERGE_DIST = 250_000
CHUNKSIZE  = 300_000

MALE_COL   = '#2166AC'
FEMALE_COL = '#E75480'

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})


# ══ STEP 1 — load FUMA lead SNPs ══════════════════════════════════════════════
print('\n=== Step 1: load FUMA lead SNPs ===')

# SNPs removed after MAF/population-specificity QC:
# Male:  rs114584016 (PDE4D, MAF<1%), rs77770701 (ZC3H8, inconsistent databases),
#         rs2688717 (TRIM48, absent in Europeans)
# Female: rs6892286 (HCN1, imputation artifact), rs114018234 (C2orf73, inflated beta at MAF~1%),
#          rs144305284 (RP11-189B4.6, MAF<1%), rs115045876 (ADRA2C, population-specific)
EXCLUDE_RSIDS = {
    'rs114584016', 'rs77770701', 'rs2688717',
    'rs6892286', 'rs114018234', 'rs144305284', 'rs115045876',
}

male_leads = pd.read_csv(MALE_LEADS, sep='\t')
male_leads = male_leads[~male_leads['rsID'].isin(EXCLUDE_RSIDS)]
male_leads = male_leads.rename(columns={
    'chr': 'chromosome',
    'pos': 'base_pair_location',
    'nearest_gene': 'gene',
})
male_leads['sex'] = 'male'
male_leads = male_leads[['chromosome', 'base_pair_location', 'gene', 'sex']]

female_leads = pd.read_csv(FEMALE_LEADS, sep='\t')
female_leads = female_leads[~female_leads['rsID'].isin(EXCLUDE_RSIDS)]
female_leads = female_leads.rename(columns={
    'Chr': 'chromosome',
    'Pos_hg19': 'base_pair_location',
    'Nearest_Gene': 'gene',
})
female_leads['sex'] = 'female'
female_leads = female_leads[['chromosome', 'base_pair_location', 'gene', 'sex']]

all_leads = pd.concat([male_leads, female_leads], ignore_index=True)
print(f'  Male leads: {len(male_leads)}, Female leads: {len(female_leads)}')


# ══ STEP 2 — deduplicate cross-sex near-duplicates ════════════════════════════
# FUMA already defined independent loci via LD clumping. The only merging needed
# is when male and female runs independently identified lead SNPs at the same
# locus but at slightly different positions (e.g. PNPLA3: 44,324,727 vs
# 44,324,855). Without this step those would appear as two rows in the plot.
# We use 250 kb to match FUMA's mergeDist, but this is deduplication, not
# locus definition.
print('\n=== Step 2: deduplicate cross-sex near-duplicates (250 kb) ===')

all_leads = all_leads.sort_values(['chromosome', 'base_pair_location']).reset_index(drop=True)

locus_ids = []
locus_id  = 0
prev_chr  = None
prev_pos  = -np.inf
for _, row in all_leads.iterrows():
    c, p = row['chromosome'], row['base_pair_location']
    if c != prev_chr or (p - prev_pos) > MERGE_DIST:
        locus_id += 1
    locus_ids.append(locus_id)
    prev_chr, prev_pos = c, p
all_leads['locus_id'] = locus_ids

n_before = len(all_leads)
n_loci   = locus_id
print(f'  {n_before} lead SNPs → {n_loci} unique loci '
      f'({n_before - n_loci} cross-sex near-duplicate(s) merged)')


# ══ STEP 3 — look up stats for all lead positions from both GWAS files ════════
print('\n=== Step 3: look up lead SNP stats from GWAS files ===')

pos_lookup = all_leads[['chromosome', 'base_pair_location']].drop_duplicates()


def lookup_positions(filepath, pos_df, label):
    print(f'  [{label}] scanning GWAS file ...')
    keep = ['chromosome', 'base_pair_location', 'p_value', 'beta', 'standard_error']
    chunks = []
    for chunk in pd.read_csv(filepath, sep='\t', compression='gzip',
                             usecols=keep, chunksize=CHUNKSIZE):
        matched = chunk.merge(pos_df, on=['chromosome', 'base_pair_location'], how='inner')
        if len(matched):
            chunks.append(matched)
    result = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=keep)
    result = result.drop_duplicates(subset=['chromosome', 'base_pair_location'])
    print(f'         → {len(result):,} positions found')
    return result


male_stats   = lookup_positions(MALE_GWAS,   pos_lookup, 'male')
female_stats = lookup_positions(FEMALE_GWAS, pos_lookup, 'female')


# ══ STEP 4 — merge cross-sex stats + het test ══════════════════════════════════
print('\n=== Step 4: merge and het test ===')

combined = (
    male_stats
    .rename(columns={'p_value': 'p_M', 'beta': 'beta_M', 'standard_error': 'se_M'})
    .merge(
        female_stats
        .rename(columns={'p_value': 'p_F', 'beta': 'beta_F', 'standard_error': 'se_F'}),
        on=['chromosome', 'base_pair_location'],
        how='outer'
    )
)

Z_het = (combined['beta_M'] - combined['beta_F']) / \
        np.sqrt(combined['se_M']**2 + combined['se_F']**2)
combined['Z_het'] = Z_het
combined['p_het'] = 2 * stats.norm.sf(np.abs(Z_het))
combined['min_p'] = combined[['p_M', 'p_F']].min(axis=1)

# Attach locus_id and gene label
combined = combined.merge(
    all_leads[['chromosome', 'base_pair_location', 'locus_id', 'gene']],
    on=['chromosome', 'base_pair_location'],
    how='left'
)

# Per locus: keep the row with min P as the representative
lead_idx = combined.groupby('locus_id')['min_p'].idxmin()
leads = combined.loc[lead_idx].copy().reset_index(drop=True)

print(f'  {len(leads)} representative loci')


# ══ STEP 5 — classify ══════════════════════════════════════════════════════════
def classify_locus(row):
    m_sig    = pd.notna(row['p_M']) and row['p_M'] < P_THRESH
    f_sig    = pd.notna(row['p_F']) and row['p_F'] < P_THRESH
    same_dir = (pd.notna(row['beta_M']) and pd.notna(row['beta_F']) and
                np.sign(row['beta_M']) == np.sign(row['beta_F']))
    if m_sig and f_sig and same_dir:
        return 'shared'
    elif m_sig and f_sig and not same_dir:
        return 'opposite'
    elif m_sig:
        return 'male'
    else:
        return 'female'

leads['class']   = leads.apply(classify_locus, axis=1)
leads['het_sig'] = leads['p_het'] < 0.05

print('\n=== Step 5: locus classification ===')
print(leads[['chromosome', 'base_pair_location', 'gene',
             'p_M', 'p_F', 'p_het', 'class', 'het_sig']].to_string())
print('\nClass counts:', leads['class'].value_counts().to_dict())
print('Het sig (p_het < 0.05):', leads['het_sig'].sum())


# ══ STEP 6 — top 10, sort, plot ═══════════════════════════════════════════════
leads = leads.nsmallest(15, 'min_p').reset_index(drop=True)
print(f'\nRestricted to top 15 loci for display.')

sort_order = {'shared': 0, 'opposite': 1, 'male': 2, 'female': 2}
leads['_sort'] = leads['class'].map(sort_order)
leads = leads.sort_values(['_sort', 'min_p']).reset_index(drop=True)

# Only display bars for P < 0.05; absence of bar means P >= 0.05 in that sex
DISPLAY_THRESH = 0.05
leads['logp_M'] = leads['p_M'].apply(
    lambda p: -np.log10(p) if pd.notna(p) and p < DISPLAY_THRESH else 0.0)
leads['logp_F'] = leads['p_F'].apply(
    lambda p: -np.log10(p) if pd.notna(p) and p < DISPLAY_THRESH else 0.0)

n     = len(leads)
y_pos = np.arange(n - 1, -1, -1)

fig_h = max(7, n * 0.55 + 3.0)
fig, ax = plt.subplots(figsize=(11, fig_h))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

thresh_sug = -np.log10(P_THRESH)
thresh_gws = -np.log10(GWS_THRESH)
x_max = max(leads['logp_M'].max(), leads['logp_F'].max()) * 1.10 + 0.5

for i, row in leads.iterrows():
    y = y_pos[i]

    if row['logp_M'] > 0:
        ax.plot([0, -row['logp_M']], [y, y], color=MALE_COL, lw=2.2,
                solid_capstyle='round', zorder=2, alpha=0.88)
        mk = '^' if pd.notna(row['beta_M']) and row['beta_M'] > 0 else 'v'
        ax.scatter(-row['logp_M'], y, color=MALE_COL, marker=mk, s=130,
                   zorder=3, edgecolors='white', lw=0.8)
        if pd.notna(row['beta_M']):
            ax.text(-row['logp_M'] - 0.15, y, f"{row['beta_M']:+.2f}",
                    color=MALE_COL, fontsize=10.5, ha='right', va='center', zorder=4)

    if row['logp_F'] > 0:
        ax.plot([0, row['logp_F']], [y, y], color=FEMALE_COL, lw=2.2,
                solid_capstyle='round', zorder=2, alpha=0.88)
        mk = '^' if pd.notna(row['beta_F']) and row['beta_F'] > 0 else 'v'
        ax.scatter(row['logp_F'], y, color=FEMALE_COL, marker=mk, s=130,
                   zorder=3, edgecolors='white', lw=0.8)
        if pd.notna(row['beta_F']):
            ax.text(row['logp_F'] + 0.15, y, f"{row['beta_F']:+.2f}",
                    color=FEMALE_COL, fontsize=10.5, ha='left', va='center', zorder=4)

    label = f'$\\it{{{row["gene"]}}}$'
    ax.text(0, y, label, ha='center', va='center', fontsize=9.5, zorder=4,
            bbox=dict(facecolor='white', edgecolor='none', pad=2.0))

ax.axvline(0, color='#333333', lw=0.9, zorder=1)
for sign in (-1, 1):
    ax.axvline(sign * thresh_sug, color='#636363', ls='--', lw=1.3, alpha=0.85, zorder=0)
    ax.axvline(sign * thresh_gws, color='#c0392b', ls='--', lw=1.3, alpha=0.85, zorder=0)

# Threshold labels — horizontal, placed just above the top data row
y_label = n + 0.18
for sign in (-1, 1):
    ax.text(sign * thresh_sug, y_label, '$P{=}1{\\times}10^{-5}$',
            color='#636363', fontsize=10.5, ha='center', va='bottom')
    ax.text(sign * thresh_gws, y_label, '$P{=}5{\\times}10^{-8}$',
            color='#c0392b', fontsize=10.5, ha='center', va='bottom')

ax.set_xlim(-x_max, x_max)
ax.set_ylim(-0.8, n + 0.3)
ax.set_yticks([])

max_tick = int(np.floor(x_max))
ticks = list(range(-max_tick, max_tick + 1))
ax.set_xticks(ticks)
ax.set_xticklabels([str(abs(t)) for t in ticks], fontsize=10)
ax.spines[['left', 'top', 'right']].set_visible(False)
ax.spines['bottom'].set_linewidth(1.2)

ax.set_title(
    'Mirrored forest plot \u2014 FUMA lead SNPs (hepatic fat GWAS)\n'
    'Same locus shown side-by-side: male (blue) left, female (pink) right',
    fontsize=12, fontweight='bold', pad=14, loc='center')

# Matplotlib fills ncol=3 top-to-bottom per column, so order is:
# col0: Male, Female  |  col1: β>0, β<0  |  col2: P=1e-5, P=5e-8
legend_elements = [
    mpatches.Patch(facecolor=MALE_COL,   edgecolor='none', label='Male signal'),
    mpatches.Patch(facecolor=FEMALE_COL, edgecolor='none', label='Female signal'),
    plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='#555',
               markersize=10, label=r'$\beta > 0$ (increases fat)'),
    plt.Line2D([0], [0], marker='v', color='w', markerfacecolor='#555',
               markersize=10, label=r'$\beta < 0$ (decreases fat)'),
    plt.Line2D([0], [0], color='#636363', ls='--', lw=1.2,
               label='$P = 1\\times10^{-5}$'),
    plt.Line2D([0], [0], color='#c0392b', ls='--', lw=1.2,
               label='$P = 5\\times10^{-8}$'),
]

plt.tight_layout(rect=[0, 0.11, 1, 1])
ax_bottom = ax.get_position().y0
tick_gap  = 0.028
fig.text(0.27, ax_bottom - tick_gap,
         r'$\leftarrow$ Male $-\log_{10}(P\text{-value})$',
         ha='center', va='top', fontsize=11)
fig.text(0.73, ax_bottom - tick_gap,
         r'Female $-\log_{10}(P\text{-value})$ $\rightarrow$',
         ha='center', va='top', fontsize=11)

ax.legend(handles=legend_elements, loc='lower center',
          bbox_to_anchor=(0.5, 0.005), bbox_transform=fig.transFigure,
          ncol=3, fontsize=10.5, framealpha=0.95, edgecolor='#cccccc', borderpad=0.8)

plt.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'\nSaved → {OUT_PNG}')

#!/usr/bin/env python3
"""
Mirrored forest plot — sex-stratified hepatic fat GWAS

Method:
  1. Filter each GWAS (maf01.hg19) to P < 1e-5
  2. Union positions across both sexes
  3. Re-scan both files to get cross-sex stats for every union position
  4. Compute heterogeneity Z-test: Z = (β_M - β_F) / sqrt(SE_M² + SE_F²)
  5. Cluster into loci by 250 kb merge; pick lead SNP (min P across sexes)
  6. Classify: shared (P < 1e-5 in both, same direction)
               male / female (P < 1e-5 in one sex only)
               opposite (P < 1e-5 in both, opposite direction)
  7. Asterisk (*) on label = p_het < 0.05 (significant effect-size heterogeneity)

Note: alleles are assumed concordant between sexes (same PMBB pipeline/VCF).
"""

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Paths ─────────────────────────────────────────────────────────────────────
MALE_GWAS   = '/Users/agaro/Documents/meno_effects_liver/gwas_processed/PMBB_ALL_M.hepatic_fat.gwas.saige.fuma.maf01.hg19.gz'
FEMALE_GWAS = '/Users/agaro/Documents/meno_effects_liver/gwas_processed/PMBB_ALL_F.hepatic_fat.gwas.saige.fuma.maf01.hg19.gz'
MALE_GENES   = '/Users/agaro/Documents/meno_effects_liver/fuma_results/FUMA_male_genes.txt'
FEMALE_GENES = '/Users/agaro/Documents/meno_effects_liver/fuma_results/FUMA_female_genes.txt'
OUT_PNG     = '/Users/agaro/Documents/meno_effects_liver/sex_specific_loci_lollipop.png'

# ── Parameters ────────────────────────────────────────────────────────────────
P_THRESH   = 1e-5
GWS_THRESH = 5e-8
MERGE_DIST = 250_000   # bp — matches FUMA mergeDist
CHUNKSIZE  = 300_000

# Loci excluded after MAF/population-specificity QC (chr, center_pos)
# Entire ±250 kb window around each flagged position is dropped
# Male:   PDE4D (MAF<1%), ZC3H8 (inconsistent databases), TRIM48 (absent in Europeans)
# Female: HCN1 (imputation artifact), C2orf73 (inflated beta at MAF~1%),
#         RP11-189B4.6 (MAF<1%), ADRA2C (population-specific)
EXCLUDE_REGIONS = [
    (5,  58_294_396),   # PDE4D      male
    (2,  112_977_265),  # ZC3H8      male
    (11, 54_907_648),   # TRIM48     male
    (5,  45_928_575),   # HCN1       female
    (2,  54_605_940),   # C2orf73    female
    (13, 47_034_337),   # RP11-189B4.6 female
    (4,  3_776_272),    # ADRA2C     female
    (6,  137_859_147),  # SNORD112   female - not in FUMA top 10
    (22, 40_092_645),   # CACNA1I    female - not in FUMA top 10
]

def in_excluded_region(chrom, pos):
    for exc_chr, exc_pos in EXCLUDE_REGIONS:
        if chrom == exc_chr and abs(pos - exc_pos) <= MERGE_DIST:
            return True
    return False

MALE_COL   = '#2166AC'
FEMALE_COL = '#E75480'
SHARED_COL = '#6A4C93'

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})


# ── Helper: stream GWAS file, keep rows below p threshold ─────────────────────
def filter_gwas(filepath, p_thresh, label):
    print(f'  [{label}] scanning for P < {p_thresh:.0e} ...')
    keep = ['chromosome', 'base_pair_location', 'p_value', 'beta', 'standard_error']
    chunks = []
    for chunk in pd.read_csv(filepath, sep='\t', compression='gzip',
                             usecols=keep, chunksize=CHUNKSIZE):
        sub = chunk[chunk['p_value'] < p_thresh]
        if len(sub):
            mask = sub.apply(
                lambda r: not in_excluded_region(r['chromosome'], r['base_pair_location']),
                axis=1)
            sub = sub[mask]
        if len(sub):
            chunks.append(sub)
    result = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=keep)
    print(f'         → {len(result):,} variants')
    return result


# ── Helper: stream GWAS file, keep rows matching a position set ───────────────
def lookup_positions(filepath, pos_df, label):
    print(f'  [{label}] cross-sex stat lookup ...')
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


# ══ STEP 1 — significant variants per sex ════════════════════════════════════
print('\n=== Step 1: filter to P < 1e-5 ===')
male_sig   = filter_gwas(MALE_GWAS,   P_THRESH, 'male')
female_sig = filter_gwas(FEMALE_GWAS, P_THRESH, 'female')

# ══ STEP 2 — union position set ══════════════════════════════════════════════
pos_union = pd.concat([
    male_sig[['chromosome', 'base_pair_location']],
    female_sig[['chromosome', 'base_pair_location']]
]).drop_duplicates().reset_index(drop=True)
print(f'\n=== Step 2: {len(pos_union):,} unique positions in union ===')

# ══ STEP 3 — cross-sex stats for all union positions ═════════════════════════
print('\n=== Step 3: cross-sex lookup ===')
male_stats   = lookup_positions(MALE_GWAS,   pos_union, 'male')
female_stats = lookup_positions(FEMALE_GWAS, pos_union, 'female')

# ══ STEP 4 — merge + het test ═════════════════════════════════════════════════
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

# ══ STEP 5 — cluster into loci by 250 kb merge ═══════════════════════════════
combined = combined.sort_values(['chromosome', 'base_pair_location']).reset_index(drop=True)

locus_ids = []
locus_id  = 0
prev_chr  = None
prev_pos  = -np.inf
for _, row in combined.iterrows():
    c, p = row['chromosome'], row['base_pair_location']
    if c != prev_chr or (p - prev_pos) > MERGE_DIST:
        locus_id += 1
    locus_ids.append(locus_id)
    prev_chr, prev_pos = c, p
combined['locus_id'] = locus_ids

print(f'\n=== Step 5: {locus_id} loci after 250 kb merge ===')

# ══ STEP 6 — lead SNP per locus (min P across either sex) ════════════════════
combined['min_p'] = combined[['p_M', 'p_F']].min(axis=1)
lead_idx = combined.groupby('locus_id')['min_p'].idxmin()
leads    = combined.loc[lead_idx].copy().reset_index(drop=True)

# ══ STEP 7 — classify loci ════════════════════════════════════════════════════
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

print('\n=== Step 7: locus classification ===')
print(leads[['chromosome', 'base_pair_location', 'p_M', 'p_F',
             'beta_M', 'beta_F', 'p_het', 'class', 'het_sig']].to_string())
print('\nClass counts:', leads['class'].value_counts().to_dict())
print('Het sig (p_het < 0.05):', leads['het_sig'].sum())

# ── Restrict to top 10 loci by min P across either sex ───────────────────────
leads = leads.nsmallest(10, 'min_p').reset_index(drop=True)
print(f'\nRestricted to top 10 loci for display.')
print(leads[['chromosome','base_pair_location','p_M','p_F','beta_M','beta_F','p_het','class','het_sig']].to_string())
print(f"Top 10 het_sig count: {leads['het_sig'].sum()}")
# print gene-annotated table after labels are assigned — done below

# ══ STEP 8 — gene annotation from FUMA genes.txt files ════════════════════════
# Each genes.txt row covers one gene mapped to a locus; we pick the best label
# per locus: protein-coding genes preferred, then closest gene to lead SNP.
def load_genes(filepath):
    df = pd.read_csv(filepath, sep='\t', usecols=['symbol', 'chr', 'start', 'end'])
    df['mid'] = (df['start'] + df['end']) // 2
    return df

genes_all = pd.concat([load_genes(MALE_GENES), load_genes(FEMALE_GENES)]) \
              .drop_duplicates(subset=['symbol', 'chr']).reset_index(drop=True)


def nearest_gene(chr_, pos):
    """Return the gene whose midpoint is closest to pos on chr_."""
    same = genes_all[genes_all['chr'] == chr_].copy()
    if same.empty:
        return f'chr{chr_}:{pos}'
    same['dist'] = (same['mid'] - pos).abs()
    return same.sort_values('dist').iloc[0]['symbol']


leads['gene'] = leads.apply(
    lambda r: nearest_gene(r['chromosome'], r['base_pair_location']), axis=1)

LABEL_OVERRIDES = {
    (8,  9_183_596):  'TNKS',
    (14, 21_587_865): 'ARHGEF40',
    (22, 48_615_615): 'MIR3201',
}
leads['gene'] = leads.apply(
    lambda r: LABEL_OVERRIDES.get((r['chromosome'], r['base_pair_location']), r['gene']),
    axis=1)

print('\n=== Top 10 loci with gene labels ===')
print(leads[['gene','chromosome','p_M','p_F','beta_M','beta_F','p_het','class','het_sig']].to_string())

# ══ STEP 9 — sort: shared first, then by min_p ═══════════════════════════════
sort_order = {'shared': 0, 'opposite': 1, 'male': 2, 'female': 2}
leads['_sort'] = leads['class'].map(sort_order)
leads = leads.sort_values(['_sort', 'min_p']).reset_index(drop=True)

# ══ STEP 10 — plot ════════════════════════════════════════════════════════════
leads['logp_M'] = leads['p_M'].apply(
    lambda p: -np.log10(p) if pd.notna(p) and p > 0 else 0.0)
leads['logp_F'] = leads['p_F'].apply(
    lambda p: -np.log10(p) if pd.notna(p) and p > 0 else 0.0)

n     = len(leads)
y_pos = np.arange(n - 1, -1, -1)   # most significant at top

fig_h = max(7, n * 0.55 + 3.0)
fig, ax = plt.subplots(figsize=(11, fig_h))
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

thresh_sug = -np.log10(P_THRESH)    # 5.0
thresh_gws = -np.log10(GWS_THRESH)  # ≈7.3

x_max = max(leads['logp_M'].max(), leads['logp_F'].max()) * 1.10 + 0.5

for i, row in leads.iterrows():
    y = y_pos[i]

    # Male bar — extends LEFT, always blue
    if row['logp_M'] > 0:
        ax.plot([0, -row['logp_M']], [y, y], color=MALE_COL, lw=2.2,
                solid_capstyle='round', zorder=2, alpha=0.88)
        mk = '^' if pd.notna(row['beta_M']) and row['beta_M'] > 0 else 'v'
        ax.scatter(-row['logp_M'], y, color=MALE_COL, marker=mk, s=130,
                   zorder=3, edgecolors='white', lw=0.8)

    # Female bar — extends RIGHT, always pink
    if row['logp_F'] > 0:
        ax.plot([0, row['logp_F']], [y, y], color=FEMALE_COL, lw=2.2,
                solid_capstyle='round', zorder=2, alpha=0.88)
        mk = '^' if pd.notna(row['beta_F']) and row['beta_F'] > 0 else 'v'
        ax.scatter(row['logp_F'], y, color=FEMALE_COL, marker=mk, s=130,
                   zorder=3, edgecolors='white', lw=0.8)

    # Gene label centered; * if het significant
    label = f'$\\it{{{row["gene"]}}}$'
    if row['het_sig']:
        label += ' *'
    ax.text(0, y, label, ha='center', va='center', fontsize=9.5, zorder=4,
            bbox=dict(facecolor='white', edgecolor='none', pad=2.0))

# Center axis
ax.axvline(0, color='#333333', lw=0.9, zorder=1)

# Threshold reference lines (both sides)
for sign in (-1, 1):
    ax.axvline(sign * thresh_sug, color='#636363', ls='--', lw=1.3,
               alpha=0.85, zorder=0)
    ax.axvline(sign * thresh_gws, color='#c0392b', ls='--', lw=1.3,
               alpha=0.85, zorder=0)

# Threshold labels at top of plot
y_top = n - 0.15
ax.text(-thresh_sug, y_top, '$P{=}1{\\times}10^{-5}$',
        color='#636363', fontsize=8, ha='center', va='bottom')
ax.text( thresh_sug, y_top, '$P{=}1{\\times}10^{-5}$',
        color='#636363', fontsize=8, ha='center', va='bottom')
ax.text(-thresh_gws, y_top, '$P{=}5{\\times}10^{-8}$',
        color='#c0392b', fontsize=8, ha='center', va='bottom')
ax.text( thresh_gws, y_top, '$P{=}5{\\times}10^{-8}$',
        color='#c0392b', fontsize=8, ha='center', va='bottom')

# Axes
ax.set_xlim(-x_max, x_max)
ax.set_ylim(-0.8, n + 0.3)
ax.set_yticks([])

max_tick = int(np.floor(x_max))
ticks = list(range(-max_tick, max_tick + 1))
ax.set_xticks(ticks)
ax.set_xticklabels([str(abs(t)) for t in ticks], fontsize=10)
ax.spines[['left', 'top', 'right']].set_visible(False)
ax.spines['bottom'].set_linewidth(1.2)

# Title
ax.set_title(
    'Mirrored forest plot \u2014 sex-specific lead SNPs (hepatic fat GWAS)\n'
    'Same locus shown side-by-side: male (blue) left, female (pink) right',
    fontsize=12, fontweight='bold', pad=14, loc='center')

# Legend
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
    plt.Line2D([0], [0], color='none', marker='$*$', markerfacecolor='black',
               markersize=11, label='$p_{\\mathrm{het}} < 0.05$'),
]

# Layout: tight_layout first, then place labels relative to actual axes position
plt.tight_layout(rect=[0, 0.11, 1, 1])

# Axis direction labels — placed just below tick labels using actual axes position
ax_bottom = ax.get_position().y0   # axes bottom in figure coords
tick_gap  = 0.028                  # approximate height of tick labels in figure coords
fig.text(0.27, ax_bottom - tick_gap,
         r'$\leftarrow$ Male $-\log_{10}(P\text{-value})$',
         ha='center', va='top', fontsize=11)
fig.text(0.73, ax_bottom - tick_gap,
         r'Female $-\log_{10}(P\text{-value})$ $\rightarrow$',
         ha='center', va='top', fontsize=11)

# Legend — anchored to figure bottom, below direction labels
ax.legend(handles=legend_elements, loc='lower center',
          bbox_to_anchor=(0.5, 0.005), bbox_transform=fig.transFigure,
          ncol=4, fontsize=9, framealpha=0.95, edgecolor='#cccccc', borderpad=0.8)
plt.savefig(OUT_PNG, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f'\nSaved → {OUT_PNG}')

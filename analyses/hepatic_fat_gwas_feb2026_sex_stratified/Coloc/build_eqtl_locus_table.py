#!/usr/bin/env python3
"""
Build per-locus eQTL summary table — sex-stratified hepatic fat GWAS (PMBB).

Inputs:
  - FUMA snps.txt from male and female runs
  - FUMA eqtl.txt from male and female runs
  - Raw GWAS .gz files (hg19, MAF-filtered) for cross-sex stat lookup
  - Curated lead SNP TSV files (already in the working directory)

Outputs:
  - locus_eqtl_summary.tsv  — one row per merged locus, with:
      Locus_ID, Chr, Lead_SNP, Nearest_Gene, Sex_detected,
      p_M, p_F, beta_M, beta_F, se_M, se_F, Z_het, p_het,
      Classification, eQTL_genes_liver, eQTL_genes_adipose, eQTL_genes_other

Classification (three-tier):
  Shared          — P < 1e-5 in both sexes, same direction, p_het > 0.05
  Sex-enriched    — P < 1e-5 in one sex only
  (no het test)   — het Z cannot be computed (one sex missing)

Z_het = (beta_M - beta_F) / sqrt(SE_M^2 + SE_F^2)
p_het = 2 * Phi(-|Z_het|)

Tissue priority for display: Liver > Adipose Visceral > Adipose Subcutaneous > Blood > other
"""

import numpy as np
import pandas as pd
from scipy import stats

# ── Paths — edit to match where you placed FUMA output files ──────────────────
MALE_GWAS    = '/Users/agaro/Documents/meno_effects_liver/gwas_processed/PMBB_ALL_M.hepatic_fat.gwas.saige.fuma.maf01.hg19.gz'
FEMALE_GWAS  = '/Users/agaro/Documents/meno_effects_liver/gwas_processed/PMBB_ALL_F.hepatic_fat.gwas.saige.fuma.maf01.hg19.gz'
MALE_LEADS   = '/Users/agaro/Documents/meno_effects_liver/fuma_results/PMBB_male_top10_leadSNPS_FUMA.tsv'
FEMALE_LEADS = '/Users/agaro/Documents/meno_effects_liver/fuma_results/PMBB_female_top10_loci_FUMA.tsv'

# FUMA output files — place downloaded files here (or update paths)
MALE_SNPS    = '/Users/agaro/Documents/meno_effects_liver/fuma_results/fuma_male/snps.txt'
MALE_EQTL    = '/Users/agaro/Documents/meno_effects_liver/fuma_results/fuma_male/eqtl.txt'
FEMALE_SNPS  = '/Users/agaro/Documents/meno_effects_liver/fuma_results/fuma_female/snps.txt'
FEMALE_EQTL  = '/Users/agaro/Documents/meno_effects_liver/fuma_results/fuma_female/eqtl.txt'

OUT_TSV = '/Users/agaro/Documents/meno_effects_liver/locus_eqtl_summary.tsv'

# ── Parameters ─────────────────────────────────────────────────────────────────
P_MAIN     = 1e-5     # primary threshold for three-tier classification
P_HET_SIG  = 0.05     # heterogeneity significance cutoff (shared requires p_het > 0.05)
MERGE_DIST = 250_000  # cross-sex deduplication window (matches FUMA mergeDist)
EQTL_DIST  = 500_000  # window to assign eQTL SNPs to a lead locus
CHUNKSIZE  = 300_000  # rows per chunk when scanning GWAS .gz files

# Tissue display priority (lower = more relevant to hepatic fat)
TISSUE_PRIORITY = [
    'Liver',
    'Adipose Visceral Omentum',
    'Adipose Subcutaneous',
    'Whole Blood',
]


# ── Helpers ────────────────────────────────────────────────────────────────────
def tissue_sort_key(tissue):
    tl = tissue.lower()
    for rank, name in enumerate(TISSUE_PRIORITY):
        if name.lower() in tl:
            return rank
    return len(TISSUE_PRIORITY)


def categorize_tissue(tissue):
    tl = tissue.lower()
    if 'liver' in tl:
        return 'liver'
    if 'adipose' in tl:
        return 'adipose'
    return 'other'


def format_gene_tissues(gene_tissue_dict, category):
    """
    For the given tissue category return a pipe-separated string:
        GENE1 (tissue_a; tissue_b) | GENE2 (tissue_c)
    Genes are sorted alphabetically; tissues by priority within category.
    Returns '—' when no genes map to that category.
    """
    entries = []
    for gene in sorted(gene_tissue_dict):
        tissues_in_cat = [t for t in gene_tissue_dict[gene]
                          if categorize_tissue(t) == category]
        if not tissues_in_cat:
            continue
        tissues_in_cat.sort(key=tissue_sort_key)
        # Shorten verbose GTEx tissue names for readability
        short = [t.replace('Adipose Visceral Omentum', 'Adipose Visceral')
                  .replace('Adipose Subcutaneous', 'Adipose Subcut.')
                  for t in tissues_in_cat]
        entries.append(f'{gene} ({"; ".join(short)})')
    return ' | '.join(entries) if entries else '—'


# ══ STEP 1 — Load curated lead SNPs from both sexes ═══════════════════════════
print('\n=== Step 1: Load lead SNPs ===')

male_leads = pd.read_csv(MALE_LEADS, sep='\t').rename(columns={
    'chr': 'chromosome',
    'pos': 'base_pair_location',
    'nearest_gene': 'gene',
})
male_leads['sex'] = 'male'

female_leads = pd.read_csv(FEMALE_LEADS, sep='\t').rename(columns={
    'Chr': 'chromosome',
    'Pos_hg19': 'base_pair_location',
    'Nearest_Gene': 'gene',
})
female_leads['sex'] = 'female'

# Minimal union for locus building
all_leads = pd.concat([
    male_leads[['chromosome', 'base_pair_location', 'gene', 'sex', 'rsID']],
    female_leads[['chromosome', 'base_pair_location', 'gene', 'sex', 'rsID']],
], ignore_index=True)

print(f'  Male: {len(male_leads)} leads, Female: {len(female_leads)} leads')


# ══ STEP 2 — Deduplicate cross-sex near-duplicates (250 kb) ══════════════════
# Same logic as plot_fuma_lead_snps.py — merges cases where male/female FUMA
# runs picked slightly different positions for the same locus (e.g. PNPLA3).
print('\n=== Step 2: Deduplicate cross-sex near-duplicates ===')

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

n_loci = locus_id
print(f'  {len(all_leads)} lead SNPs → {n_loci} unique loci')

# Locus reference table: one row per unique (chr, pos, locus_id)
locus_ref = (all_leads[['chromosome', 'base_pair_location', 'locus_id']]
             .drop_duplicates()
             .reset_index(drop=True))


# ══ STEP 3 — Load FUMA snps.txt (informational; confirms loci) ════════════════
print('\n=== Step 3: Load FUMA snps.txt ===')

male_snps   = pd.read_csv(MALE_SNPS,   sep='\t')
female_snps = pd.read_csv(FEMALE_SNPS, sep='\t')
print(f'  Male:   {len(male_snps):,} SNPs across {male_snps["GenomicLocus"].nunique()} loci')
print(f'  Female: {len(female_snps):,} SNPs across {female_snps["GenomicLocus"].nunique()} loci')


# ══ STEP 4 — Load FUMA eqtl.txt, filter to eqtlMapFilt == 1 ══════════════════
print('\n=== Step 4: Load FUMA eqtl.txt ===')

def load_eqtl(path, sex_label):
    df = pd.read_csv(path, sep='\t')
    print(f'  [{sex_label}] raw rows: {len(df):,}  columns: {list(df.columns)}')
    if 'eqtlMapFilt' in df.columns:
        df = df[df['eqtlMapFilt'] == 1].copy()
        print(f'           after eqtlMapFilt==1: {len(df):,}')
    # Ensure numeric chr/pos
    if 'chr' in df.columns:
        df['chr'] = pd.to_numeric(df['chr'], errors='coerce')
    if 'pos' in df.columns:
        df['pos'] = pd.to_numeric(df['pos'], errors='coerce')
    # Fall back to parsing uniqID if chr/pos missing
    elif 'uniqID' in df.columns:
        parts = df['uniqID'].str.split(':', expand=True)
        df['chr'] = pd.to_numeric(parts[0], errors='coerce')
        df['pos'] = pd.to_numeric(parts[1], errors='coerce')
    return df

male_eqtl   = load_eqtl(MALE_EQTL,   'male')
female_eqtl = load_eqtl(FEMALE_EQTL, 'female')


# ══ STEP 5 — Map eQTL genes to merged loci ════════════════════════════════════
print('\n=== Step 5: Map eQTL genes to loci ===')

def map_eqtl_to_loci(eqtl_df, locus_ref_df, window=EQTL_DIST):
    """
    For each eQTL SNP, find the nearest lead locus within `window` bp on the
    same chromosome.  Returns dict: locus_id → {symbol → set(tissues)}.

    Strategy: merge on chromosome, compute distance, filter to window, keep
    nearest locus per eQTL SNP.  Fast because locus_ref is tiny (~20-30 rows).
    """
    cols_needed = ['chr', 'pos', 'symbol', 'tissue']
    missing = [c for c in cols_needed if c not in eqtl_df.columns]
    if missing:
        raise ValueError(f'eqtl.txt missing expected columns: {missing}')

    sub = eqtl_df[cols_needed].dropna(subset=['chr', 'pos']).copy()
    sub['chr'] = sub['chr'].astype(int)

    # Cross-join on chromosome (locus_ref is tiny — no memory concern)
    ref = locus_ref_df.rename(columns={'chromosome': 'chr',
                                        'base_pair_location': 'locus_pos'})
    ref['chr'] = ref['chr'].astype(int)

    merged = sub.merge(ref, on='chr', how='left')
    merged['dist'] = (merged['pos'] - merged['locus_pos']).abs()
    merged = merged[merged['dist'] <= window].copy()

    # Keep only nearest locus per (chr, pos, symbol, tissue) tuple
    idx = merged.groupby(['chr', 'pos', 'symbol', 'tissue'])['dist'].idxmin()
    merged = merged.loc[idx]

    result = {}
    for _, row in merged.iterrows():
        lid = int(row['locus_id'])
        sym = str(row['symbol'])
        tis = str(row['tissue'])
        result.setdefault(lid, {}).setdefault(sym, set()).add(tis)
    return result

male_locus_eqtl   = map_eqtl_to_loci(male_eqtl,   locus_ref)
female_locus_eqtl = map_eqtl_to_loci(female_eqtl, locus_ref)

def _count(d):
    return sum(len(v) for v in d.values())

print(f'  eQTL gene-locus pairs: {_count(male_locus_eqtl)} (male), '
      f'{_count(female_locus_eqtl)} (female)')


# ══ STEP 6 — Look up cross-sex GWAS stats from raw files ══════════════════════
print('\n=== Step 6: Look up GWAS stats from .gz files ===')

pos_lookup = locus_ref[['chromosome', 'base_pair_location']].drop_duplicates()

def lookup_positions(filepath, pos_df, label):
    print(f'  [{label}] scanning {filepath.split("/")[-1]} ...')
    keep = ['chromosome', 'base_pair_location', 'p_value', 'beta', 'standard_error']
    chunks = []
    for chunk in pd.read_csv(filepath, sep='\t', compression='gzip',
                             usecols=keep, chunksize=CHUNKSIZE):
        matched = chunk.merge(pos_df, on=['chromosome', 'base_pair_location'], how='inner')
        if len(matched):
            chunks.append(matched)
    result = (pd.concat(chunks, ignore_index=True) if chunks
              else pd.DataFrame(columns=keep))
    result = result.drop_duplicates(subset=['chromosome', 'base_pair_location'])
    print(f'    → {len(result):,} positions found')
    return result

male_stats   = lookup_positions(MALE_GWAS,   pos_lookup, 'male')
female_stats = lookup_positions(FEMALE_GWAS, pos_lookup, 'female')


# ══ STEP 7 — Merge cross-sex stats, compute het Z, classify ══════════════════
print('\n=== Step 7: Merge stats and classify ===')

combined = (
    male_stats
    .rename(columns={'p_value': 'p_M', 'beta': 'beta_M', 'standard_error': 'se_M'})
    .merge(
        female_stats
        .rename(columns={'p_value': 'p_F', 'beta': 'beta_F', 'standard_error': 'se_F'}),
        on=['chromosome', 'base_pair_location'],
        how='outer',
    )
)

Z = ((combined['beta_M'] - combined['beta_F']) /
     np.sqrt(combined['se_M']**2 + combined['se_F']**2))
combined['Z_het'] = Z
combined['p_het'] = 2 * stats.norm.sf(np.abs(Z))
combined['min_p'] = combined[['p_M', 'p_F']].min(axis=1)

# Attach locus_id
combined = combined.merge(locus_ref, on=['chromosome', 'base_pair_location'], how='left')

# Per locus: representative position = min P
lead_idx = combined.groupby('locus_id')['min_p'].idxmin()
leads_df = combined.loc[lead_idx].copy().reset_index(drop=True)

# Attach gene label and sex from original lead list
# When both sexes reported a locus, prefer the one with lower p
gene_labels = (all_leads
               .sort_values('locus_id')
               .drop_duplicates(subset=['locus_id'], keep='first')
               [['locus_id', 'gene', 'sex', 'rsID']])
leads_df = leads_df.merge(gene_labels, on='locus_id', how='left')


def classify_locus(row):
    m_main = pd.notna(row['p_M']) and row['p_M'] < P_MAIN
    f_main = pd.notna(row['p_F']) and row['p_F'] < P_MAIN
    het    = pd.notna(row['p_het']) and row['p_het'] < P_HET_SIG
    same_d = (pd.notna(row['beta_M']) and pd.notna(row['beta_F']) and
              np.sign(row['beta_M']) == np.sign(row['beta_F']))

    if m_main and f_main and same_d and not het:
        return 'Shared'
    if m_main and f_main and not same_d:
        return 'Opposite direction'
    if m_main and not f_main:
        return 'Male sex-enriched'
    if f_main and not m_main:
        return 'Female sex-enriched'
    return 'Unclear'

leads_df['Classification'] = leads_df.apply(classify_locus, axis=1)

print('\nClassification counts:')
print(leads_df['Classification'].value_counts().to_string())
print(f'\nHet p < 0.05: {(leads_df["p_het"] < 0.05).sum()} / {len(leads_df)} loci')


# ══ STEP 8 — Build output table with eQTL annotations ════════════════════════
print('\n=== Step 8: Annotate and write output ===')

def get_merged_eqtl(locus_id):
    """Merge male + female eQTL gene-tissue dicts for a locus."""
    merged = {}
    for src in [male_locus_eqtl.get(locus_id, {}),
                female_locus_eqtl.get(locus_id, {})]:
        for gene, tissues in src.items():
            merged.setdefault(gene, set()).update(tissues)
    return merged

rows = []
for _, row in leads_df.iterrows():
    lid = int(row['locus_id'])
    gt  = get_merged_eqtl(lid)

    rows.append({
        'Locus_ID':            lid,
        'Chr':                 int(row['chromosome']),
        'Lead_SNP':            row.get('rsID', ''),
        'Nearest_Gene':        row.get('gene', ''),
        'Sex_detected':        row.get('sex', ''),
        'p_M':                 f"{row['p_M']:.3e}" if pd.notna(row['p_M']) else 'NA',
        'p_F':                 f"{row['p_F']:.3e}" if pd.notna(row['p_F']) else 'NA',
        'beta_M':              f"{row['beta_M']:.4f}" if pd.notna(row['beta_M']) else 'NA',
        'beta_F':              f"{row['beta_F']:.4f}" if pd.notna(row['beta_F']) else 'NA',
        'se_M':                f"{row['se_M']:.4f}" if pd.notna(row['se_M']) else 'NA',
        'se_F':                f"{row['se_F']:.4f}" if pd.notna(row['se_F']) else 'NA',
        'Z_het':               f"{row['Z_het']:.3f}" if pd.notna(row['Z_het']) else 'NA',
        'p_het':               f"{row['p_het']:.3f}" if pd.notna(row['p_het']) else 'NA',
        'Classification':      row['Classification'],
        'eQTL_genes_liver':    format_gene_tissues(gt, 'liver'),
        'eQTL_genes_adipose':  format_gene_tissues(gt, 'adipose'),
        'eQTL_genes_other':    format_gene_tissues(gt, 'other'),
    })

out = pd.DataFrame(rows)

# Sort: shared first, then male-specific, then female-specific; within group by p
sort_class = {
    'Shared': 0,
    'Male sex-enriched': 1,
    'Female sex-enriched': 2,
    'Opposite direction': 3,
    'Unclear': 4,
}
out['_sort'] = out['Classification'].map(sort_class).fillna(9)
# Use numeric min_p for sorting
out['_min_p'] = leads_df['min_p'].values
out = out.sort_values(['_sort', '_min_p']).drop(columns=['_sort', '_min_p']).reset_index(drop=True)

out.to_csv(OUT_TSV, sep='\t', index=False)
print(f'\nSaved → {OUT_TSV}')
print(f'Total loci: {len(out)}')
print('\n── Preview ──────────────────────────────────────────────────────────────────')
print(out[['Chr', 'Lead_SNP', 'Nearest_Gene', 'Classification',
           'p_M', 'p_F', 'p_het',
           'eQTL_genes_liver', 'eQTL_genes_adipose']].to_string(index=False))

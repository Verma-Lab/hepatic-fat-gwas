#!/usr/bin/env python3
"""
coloc_hepatic_fat_gtex.py
Formal colocalization: PMBB hepatic fat GWAS × GTEx v8 eQTL
Analyses: male-only, female-only, ALL combined
Method: coloc.abf (approximate Bayes factors, single causal variant)
eQTL source: eQTL Catalogue (GTEx v8 full nominal stats) via pysam tabix

Run with conda Python:
    /Users/agaro/opt/anaconda3/bin/python3 scripts/coloc_hepatic_fat_gtex.py

VARIANT MATCHING:
    GWAS variant_id encodes hg38 positions (e.g. chr22_43556374_G_A).
    eQTL Catalogue variant column uses the same format.
    Matching is done on normalized hg38 key — no liftover needed.
    hg19 positions in the GWAS file are used only for the GWAS window filter.

eQTL CATALOGUE NOTE:
    Files contain all nominal cis-eQTL associations (all tested variants,
    not just significant ones). This is required for valid coloc.abf.
"""

import os
import sys
import numpy as np
import pandas as pd
import pysam
import requests
import time
from scipy.special import logsumexp

BASE_DIR = "/Users/agaro/Documents/meno_effects_liver"

# ══ CONFIG ════════════════════════════════════════════════════════════════════

GWAS_META = {
    'male':   {'file': f'{BASE_DIR}/gwas_processed/PMBB_ALL_M.hepatic_fat.gwas.saige.fuma.maf01.hg19.gz',   'N': 9048},
    'female': {'file': f'{BASE_DIR}/gwas_processed/PMBB_ALL_F.hepatic_fat.gwas.saige.fuma.maf01.hg19.gz',   'N': 9031},
    'all':    {'file': f'{BASE_DIR}/gwas_processed/PMBB_ALL_ALL.hepatic_fat.gwas.saige.fuma.maf01.hg19.gz', 'N': 18079},
}

LEAD_CFG = {
    'male':   {'file': f'{BASE_DIR}/fuma_results/PMBB_male_top10_leadSNPS_FUMA.tsv',
               'chr': 'chr',  'pos': 'pos',      'rsid': 'rsID', 'ngene': 'nearest_gene', 'eqtl': 'eQTL_gene'},
    'female': {'file': f'{BASE_DIR}/fuma_results/PMBB_female_top10_loci_FUMA.tsv',
               'chr': 'Chr',  'pos': 'Pos_hg19', 'rsid': 'rsID', 'ngene': 'Nearest_Gene', 'eqtl': 'eQTL_Genes'},
    'all':    {'file': f'{BASE_DIR}/fuma_results/PMBB_ALL_top10_leadSNPs_FUMA.tsv',
               'chr': 'chr',  'pos': 'pos',      'rsid': 'rsID', 'ngene': 'nearest_gene', 'eqtl': 'eQTL_gene'},
}

# eQTL Catalogue GTEx v8 gene expression datasets
TISSUES = [
    {'name': 'Liver',                    'qtd': 'QTD000266', 'N': 208},
    {'name': 'Whole_Blood',              'qtd': 'QTD000356', 'N': 670},
    {'name': 'Adipose_Subcutaneous',     'qtd': 'QTD000116', 'N': 581},
    {'name': 'Adipose_Visceral_Omentum', 'qtd': 'QTD000121', 'N': 469},
    {'name': 'Muscle_Skeletal',          'qtd': 'QTD000281', 'N': 702},
    {'name': 'Pancreas',                 'qtd': 'QTD000296', 'N': 305},
]

EQTL_BASE  = "https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015"
WINDOW_BP  = 500_000
PP4_THRESH = 0.8
MIN_SNPS   = 50
W_PRIOR    = 0.15 ** 2   # ABF prior variance

# eQTL Catalogue column indices
# molecular_trait_id | chr | pos | ref | alt | variant | ma_samples | maf |
# pvalue | beta | se | type | ac | an | r2 | ... | rsid
COL_GENE, COL_POS, COL_VARIANT = 0, 2, 5
COL_MAF, COL_BETA, COL_SE      = 7, 9, 10

# ══ HELPERS ══════════════════════════════════════════════════════════════════

def norm_vid(v):
    """Normalize variant ID to bare chr_pos_ref_alt key (no chr prefix, no _b38)."""
    if v.startswith('chr'):
        v = v[3:]
    if v.endswith('_b38') or v.endswith('_b37'):
        v = v.rsplit('_', 1)[0]
    return v

def flip_vid(v):
    """Swap ref/alt alleles: 22_43423038_T_C → 22_43423038_C_T."""
    p = v.split('_')
    return f"{p[0]}_{p[1]}_{p[3]}_{p[2]}" if len(p) == 4 else v

def parse_hg38_pos(vid):
    """Extract hg38 position from normalized variant_id string."""
    try:
        return int(vid.split('_')[1])
    except (IndexError, ValueError):
        return None

# ══ GTEx GENE LOOKUP ════════════════════════════════════════════════════════

_gene_cache = {}

def get_gencode_id(symbol):
    """Convert gene symbol to unversioned Ensembl ID via GTEx API."""
    if symbol in _gene_cache:
        return _gene_cache[symbol]
    try:
        r = requests.get(
            'https://gtexportal.org/api/v2/reference/gene',
            params={'geneId': symbol},
            timeout=10
        )
        data = r.json().get('data', [])
        gid = data[0]['gencodeId'].split('.')[0] if data else None
    except Exception:
        gid = None
    _gene_cache[symbol] = gid
    time.sleep(0.1)
    return gid

# ══ GWAS EXTRACTION ═════════════════════════════════════════════════════════

_gwas_cache = {}

def extract_gwas(sex, chrom, center_hg19):
    key = (sex, chrom, center_hg19)
    if key in _gwas_cache:
        return _gwas_cache[key]

    lo, hi = center_hg19 - WINDOW_BP, center_hg19 + WINDOW_BP
    f = GWAS_META[sex]['file']
    cmd = f"gzcat '{f}' | awk 'NR==1 || ($1=={chrom} && $2>={lo} && $2<={hi})'"

    try:
        df = pd.read_csv(os.popen(cmd), sep='\t')
    except Exception as e:
        print(f"  GWAS read error: {e}")
        _gwas_cache[key] = None
        return None

    if df.empty:
        _gwas_cache[key] = None
        return None

    df['vid'] = df['variant_id'].apply(norm_vid)
    df['maf'] = df['effect_allele_frequency'].apply(lambda x: min(x, 1 - x))
    df = df[(df['standard_error'] > 0) & (df['maf'] > 0) & df['vid'].notna()]
    df = df.drop_duplicates('vid').set_index('vid')

    _gwas_cache[key] = df
    return df

# ══ eQTL CATALOGUE TABIX QUERY ══════════════════════════════════════════════

_tabix_handles = {}

def get_tabix(qtd):
    if qtd not in _tabix_handles:
        url = f"{EQTL_BASE}/{qtd}/{qtd}.all.tsv.gz"
        _tabix_handles[qtd] = pysam.TabixFile(url)
    return _tabix_handles[qtd]

def query_eqtl(qtd, chrom, hg38_start, hg38_end, gene_ensg_base):
    """Fetch full nominal eQTL stats for one gene in a genomic window."""
    try:
        tbx = get_tabix(qtd)
        rows = []
        for line in tbx.fetch(str(chrom), hg38_start, hg38_end):
            fields = line.split('\t')
            if len(fields) <= COL_SE:
                continue
            # Filter to the requested gene (strip version from both sides)
            if fields[COL_GENE].split('.')[0] != gene_ensg_base:
                continue
            try:
                beta = float(fields[COL_BETA])
                se   = float(fields[COL_SE])
                maf  = float(fields[COL_MAF])
                vid  = norm_vid(fields[COL_VARIANT])
            except ValueError:
                continue
            if se > 0 and 0 < maf < 1:
                rows.append({'vid': vid, 'beta_e': beta, 'se_e': se, 'maf_e': maf})
    except Exception as e:
        print(f"tabix error: {e}")
        return None

    if not rows:
        return None
    return pd.DataFrame(rows).drop_duplicates('vid').set_index('vid')

# ══ COLOC.ABF ════════════════════════════════════════════════════════════════

def log_abf(beta, se):
    z = beta / se
    V = se ** 2
    return 0.5 * np.log(V / (V + W_PRIOR)) + 0.5 * z**2 * W_PRIOR / (V + W_PRIOR)

def coloc_abf(beta1, se1, beta2, se2, p1=1e-4, p2=1e-4, p12=1e-5):
    """
    Approximate Bayes Factor colocalization (Giambartolomei et al. 2014).
    Returns dict of PP.H0 – PP.H4.
    PP.H4 = posterior probability of shared causal variant.
    """
    lbf1  = log_abf(np.asarray(beta1, float), np.asarray(se1, float))
    lbf2  = log_abf(np.asarray(beta2, float), np.asarray(se2, float))
    lbf12 = lbf1 + lbf2

    lsum1  = logsumexp(lbf1)
    lsum2  = logsumexp(lbf2)
    lsum12 = logsumexp(lbf12)

    # H3: log of (sum1 * sum2 - sum12) — different causal variants
    inner = np.exp(lsum1 + lsum2) - np.exp(lsum12)
    lH3   = np.log(inner) if inner > 0 else -np.inf

    log_pp = np.array([
        0,                                    # H0: no association
        np.log(p1)  + lsum1,                 # H1: GWAS only
        np.log(p2)  + lsum2,                 # H2: eQTL only
        np.log(p1)  + np.log(p2) + lH3,     # H3: both, different causal
        np.log(p12) + lsum12                 # H4: both, same causal
    ])

    log_norm = logsumexp(log_pp)
    pp = np.exp(log_pp - log_norm)
    return {f'PP.H{i}': float(pp[i]) for i in range(5)}

# ══ LOAD LEAD SNPS ══════════════════════════════════════════════════════════

def load_leads(sex):
    cfg = LEAD_CFG[sex]
    df  = pd.read_csv(cfg['file'], sep='\t')
    leads = []
    for _, row in df.iterrows():
        raw = f"{row[cfg['ngene']]};{row[cfg['eqtl']]}"
        genes = list({g.strip() for g in raw.replace(',', ';').split(';')
                      if g.strip() and g.strip() not in ('—', 'NA', 'nan', '')})
        leads.append({
            'sex': sex, 'rsid': str(row[cfg['rsid']]),
            'chr': int(row[cfg['chr']]), 'pos_hg19': int(row[cfg['pos']]),
            'genes': genes
        })
    return leads

# ══ MAIN LOOP ════════════════════════════════════════════════════════════════

results = []

for sex in ['male', 'female', 'all']:
    leads = load_leads(sex)
    for locus in leads:
        chrom     = locus['chr']
        pos_hg19  = locus['pos_hg19']
        gene_syms = locus['genes']

        print(f"\n[{sex}] {locus['rsid']} | chr{chrom}:{pos_hg19} | genes: {', '.join(gene_syms)}")

        gwas_df = extract_gwas(sex, chrom, pos_hg19)
        if gwas_df is None or gwas_df.empty:
            print("  [skip] no GWAS data in window")
            continue
        print(f"  GWAS variants: {len(gwas_df)}")

        # Derive hg38 window from variant_ids embedded in GWAS data
        hg38_pos = gwas_df.index.map(parse_hg38_pos).dropna().astype(int)
        if hg38_pos.empty:
            print("  [skip] cannot determine hg38 window")
            continue
        hg38_lo = int(hg38_pos.min()) - 10_000
        hg38_hi = int(hg38_pos.max()) + 10_000

        for gene_sym in gene_syms:
            ensg_versioned = get_gencode_id(gene_sym)
            if not ensg_versioned:
                print(f"  [skip] {gene_sym}: not found in GTEx")
                continue
            ensg_base = ensg_versioned.split('.')[0]

            for tissue in TISSUES:
                print(f"  {gene_sym} × {tissue['name']} ...", end=' ', flush=True)

                eqtl_df = query_eqtl(tissue['qtd'], chrom, hg38_lo, hg38_hi, ensg_base)
                if eqtl_df is None or eqtl_df.empty:
                    print("no data")
                    continue
                print(f"{len(eqtl_df)} variants", end=' ', flush=True)

                if len(eqtl_df) < 100:
                    print("[few variants — eQTL may be absent in this tissue]", end=' ')

                # Match GWAS and eQTL on normalized hg38 variant key
                shared = gwas_df.index.intersection(eqtl_df.index)

                # Allele-flip fallback
                flipped = eqtl_df.copy()
                flipped.index = flipped.index.map(flip_vid)
                flipped['beta_e'] = -flipped['beta_e']
                extra = gwas_df.index.intersection(flipped.index).difference(shared)
                if len(extra):
                    eqtl_df = pd.concat([eqtl_df, flipped.loc[extra]])
                    shared = shared.append(extra)

                if len(shared) < MIN_SNPS:
                    print(f"only {len(shared)} shared variants — skip")
                    continue

                merged = gwas_df.loc[shared].join(eqtl_df.loc[shared], how='inner')
                merged = merged.dropna(subset=['beta', 'standard_error', 'beta_e', 'se_e'])
                if len(merged) < MIN_SNPS:
                    print(f"only {len(merged)} after merge — skip")
                    continue

                pp = coloc_abf(merged['beta'], merged['standard_error'],
                               merged['beta_e'], merged['se_e'])

                print(f"PP.H4={pp['PP.H4']:.3f}")
                results.append({
                    'gwas': sex, 'lead_rsid': locus['rsid'],
                    'chr': chrom, 'pos_hg19': pos_hg19,
                    'gene': gene_sym, 'gencode_id': ensg_base,
                    'tissue': tissue['name'], 'n_snps': len(merged),
                    **{k: round(v, 4) for k, v in pp.items()},
                    'coloc': pp['PP.H4'] >= PP4_THRESH
                })

# ══ OUTPUT ═══════════════════════════════════════════════════════════════════

if not results:
    print("\nNo colocalization results produced.")
    sys.exit(0)

out = pd.DataFrame(results).sort_values('PP.H4', ascending=False)
out_file = f"{BASE_DIR}/coloc_results_v2.tsv"
out.to_csv(out_file, sep='\t', index=False)
print(f"\nResults → {out_file} ({len(out)} rows)")

hits = out[out['coloc']]
print(f"\n── PP.H4 ≥ {PP4_THRESH}: {len(hits)} colocalized signals ──────────────────")
if len(hits):
    print(hits[['gwas', 'lead_rsid', 'gene', 'tissue', 'n_snps',
                'PP.H3', 'PP.H4']].to_string(index=False))

print(f"\nTop 20 by PP.H4:")
print(out[['gwas', 'lead_rsid', 'gene', 'tissue', 'n_snps', 'PP.H4']].head(20).to_string(index=False))

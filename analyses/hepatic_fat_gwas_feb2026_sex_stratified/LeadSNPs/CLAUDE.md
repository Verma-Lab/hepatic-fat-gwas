# LeadSNPs: FUMA Lead SNP / Locus Analysis (PMBB Hepatic Fat, Sex-Stratified)

## Overview
FUMA-based lead SNP / locus discovery for sex-stratified (M/F) and combined (ALL) SAIGE GWAS of CT-derived hepatic fat (ΔHU) in PMBB. Sample sizes: male-only N=9048, female-only N=9031, combined ALL N=18079.

Originally developed locally at `/Users/agaro/Documents/PhD/meno_effects_liver/`; this copy on remote lives at `analyses/hepatic_fat_gwas_feb2026_sex_stratified/LeadSNPs/`. See `../Coloc/CLAUDE.md` for the companion colocalization analysis, which uses the same input sumstats.

## Input provenance (resolved 2026-09-04)
- `PMBB_ALL_F` / `PMBB_ALL_M` raw sumstats: from this same analysis directory's `PMBB_ALL_F` and `PMBB_ALL_M` sex-stratified cohorts, generated ~2026-02-20, pulled to local Mac 2026-03-16.
- `PMBB_ALL_ALL` raw sumstats: from `../../hepatic_fat_gwas_jan2026/PMBB_ALL_ALL/Sumstats/`, merged 2026-02-12, pulled to local Mac 2026-02-19.
- This resolves a previously open question in this doc about whether the multi-ancestry Nextflow run under `hepatic_fat_gwas_feb2026_sex_stratified/` fed this analysis — confirmed yes, for the F/M cohorts. `hepatic_fat_gwas_manual_jan2026` was NOT the source — it only holds unmerged per-chunk SAIGE results, no merged sumstats file.

## Status of fuma_results/
**These 5 TSV/gene files are from superseded FUMA runs** (dated 2026-03-30/31): no MAF filter applied, inconsistent settings (see "Previous FUMA runs" below). Kept here for continuity/history, not as final results. Replace once the corrected `_hg19_maf0.01` FUMA jobs' outputs are confirmed and pulled — as of 2026-08-18 it was unclear whether those corrected outputs (`leadSNPs.txt` etc.) were ever pulled locally at all.

## plot_sex_specific_lollipop.py — known issue (as of 2026-08-18)
Computes p_het correctly (`Z_het = (β_male − β_female) / √(SE_male² + SE_female²)`, `p_het = 2Φ(−|Z_het|)`) but does **not** gate the shared/male/female/opposite classification on p_het — p_het is only overlaid as a `*` annotation, not used to define the class. Confirmed misclassifications:
- **TM6SF2** labeled "male" despite p_het=0.099 (not significant) — a well-established both-sexes liver-fat locus, so this is a false sex-specific call.
- **TNKS** labeled generically "male" but by the strict tier rules below (female nominal p<0.05 + p_het<0.05) is actually **sex-enriched**, not sex-specific.

`sex_specific_loci_lollipop.png` (generated 2026-05-17) is **provisional** — pre-dates this fix, do not treat as final.

(Scripts and the liftover chain file live flat at the top of this directory, not in a nested `scripts/`/`reference/` folder — matches how the rest of `hepatic_fat_gwas_feb2026_sex_stratified/` puts single scripts at the top level rather than nesting them.)

### Pending
- Fix `classify_locus()` to actually gate on p_het per the three-tier system below (currently a 4-way significance/direction split with p_het as cosmetic only)
- Re-check the full locus list (not just top 10) once tier logic is fixed — 42/46 loci in the full `combined` table have het_sig=True, worth confirming that's expected
- Regenerate lollipop plot and poster TSV once tier-consistent classification is in place
- Re-run baseline ALL FUMA with maf01 file and corrected settings (never confirmed done)

### Next immediate step
Fix the p_het-gated tier logic in `plot_sex_specific_lollipop.py`, rerun on the full locus set, regenerate the poster figure/TSV from the corrected classification.

## Sex-Specificity Classification — method
Heterogeneity Z-test + three-tier system (replaces a previous 1 Mb positional criterion):

| Tier | Criterion |
|------|-----------|
| Shared | P < 1×10⁻⁵ in both sexes, same direction, p_het > 0.05 |
| Sex-enriched | P < 1×10⁻⁵ in one sex, p < 0.05 in other, p_het < 0.05 |
| Sex-specific | P < 1×10⁻⁵ in one sex, p > 0.05 in other, p_het < 0.05 |

p_het < 0.05 required to exclude power-difference artifacts (N_M ≈ N_F ≈ 9,000).

Threshold corrected 2026-08-18: previously stated as P < 1×10⁻⁶ here, which was inconsistent with the `leadP = 1e-5` used for the actual sex-stratified FUMA submissions and with `P_THRESH` in `plot_sex_specific_lollipop.py`. 1×10⁻⁵ is correct.

## FUMA parameters — corrected runs (use for all new submissions)

| Parameter | Sex-stratified (M/F) | Baseline ALL |
|-----------|---------------------|--------------|
| Build | hg19 | hg19 |
| leadP | 1e-5 | 5e-8 |
| gwasP | 0.05 | 0.05 |
| r2 | 0.6 | 0.6 |
| r2_2 | 0.1 | 0.1 |
| mergeDist | 250 kb | 250 kb |
| ref panel | 1KG/Phase3 ALL | 1KG/Phase3 ALL |
| MHC | excluded | excluded |
| genetype | all | all |
| eQTL | ON (GTEx v8, 31 tissues) | ON (same) |
| Input file | `*.maf01.hg19.gz` | `*.maf01.hg19.gz` |

GTEx v8 tissues (31): Adipose Subcutaneous, Adipose Visceral Omentum, Adrenal Gland, Cells EBV-transformed lymphocytes, Whole Blood, Artery Aorta, Artery Coronary, Artery Tibial, Brain Amygdala, Brain Anterior cingulate cortex BA24, Brain Caudate basal ganglia, Brain Cerebellar Hemisphere, Brain Cerebellum, Brain Cortex, Brain Frontal Cortex BA9, Brain Hippocampus, Brain Hypothalamus, Brain Nucleus accumbens basal ganglia, Brain Putamen basal ganglia, Brain Spinal cord cervical c-1, Brain Substantia nigra, Breast Mammary Tissue, Colon Sigmoid, Colon Transverse, Esophagus Gastroesophageal Junction, Esophagus Mucosa, Esophagus Muscularis, Heart Atrial Appendage, Heart Left Ventricle, Kidney Cortex, Liver

## Previous FUMA runs (superseded — no MAF filter + inconsistent settings, do not use)
- Female (job 718915): eQTL OFF, protein_coding only
- Male (job 718921): eQTL ON, protein_coding only
- Baseline (job ~): leadP=5e-8, genetype=all, eQTL ON (no MAF filter)
- The files in `fuma_results/` (dated 2026-03-30/31) are from these superseded runs.

## gwas_processed/ — file provenance

### Raw GWAS format (15 cols, hg38, chr prefix)
`phenotype chromosome base_pair_location variant_id other_allele effect_allele effect_allele_count effect_allele_frequency missing_rate beta standard_error t_statistic variance p_value N`

### Step 1 — Reformat to FUMA v2 (10 cols, hg38, chr stripped)
```bash
gzcat PMBB_ALL_F.hepatic_fat.gwas.saige.gz \
| awk 'BEGIN{OFS="\t"; print "chromosome\tbase_pair_location\tvariant_id\tp_value\teffect_allele\tother_allele\teffect_allele_frequency\teffect_allele_count\tbeta\tstandard_error"}
NR>1{sub(/^chr/, "", $2); print $2, $3, $4, $14, $6, $5, $8, $7, $10, $11}' \
| gzip > PMBB_ALL_F.hepatic_fat.gwas.saige.fuma.v2.gz
# Same for M; ALL used .fuma.gz as source with awk to strip chr from $1
```
Column mapping: $2→chr, $3→pos, $4→variant_id, $14→p_value, $6→EA, $5→OA, $8→EAF, $7→EAC, $10→beta, $11→SE

### Step 2 — Liftover hg38 → hg19 (pyliftover; CrossMap v0.7.3 failed — no gwas subcommand)
```bash
pip install pyliftover
# Python script using LiftOver('hg38ToHg19.over.chain.gz')
# converts chr+pos-1 (0-based) → hg19, writes pos+1 back
```
Results: F=19.5M converted/43.8K failed, M=18.0M/41.0K, ALL=1.1M/3.4K

### Step 3 — MAF filter ≥ 1% (applied to hg19 files)
```bash
gzcat PMBB_ALL_F.hepatic_fat.gwas.saige.fuma.hg19.gz \
| awk 'BEGIN{OFS="\t"}
NR==1 {print $0, "maf"; next}
{
  maf = ($7 < 0.5 ? $7 : 1 - $7)
  if (maf >= 0.01) print $0, maf
}' \
| gzip > PMBB_ALL_F.hepatic_fat.gwas.saige.fuma.maf01.hg19.gz
# Same for M and ALL
```
`$7` = effect_allele_frequency. Appends maf column (FUMA ignores extra columns). Removes variants with MAF < 1% (e.g. GJD2 MAF=0.0016 drops out).

### SAIGE config note
`min_maf = 0` in the original SAIGE run allowed variants down to MAF ~0.22% (MAC≥40, N~9000). GWAS does NOT need rerun — the reformat step (Step 1) used maf=0.01 correctly; the MAF filter on summary stats (Step 3) is sufficient.

## File index

| File | Description | Status |
|------|-------------|--------|
| `gwas_processed/*.fuma.v2.gz` | FUMA-formatted hg38, MAF unfiltered | Keep (reproducibility) |
| `gwas_processed/*.fuma.hg19.gz` | hg19 lifted, MAF unfiltered | Keep |
| `gwas_processed/*.fuma.maf01.hg19.gz` | hg19, MAF ≥ 1% — **submitted to FUMA** | Keep |
| `hg38ToHg19.over.chain.gz` | Liftover chain used for Step 2 | Keep |
| `sex_specific_loci_lollipop.png` | Lollipop plot, generated 2026-05-17 | Provisional — pre-dates p_het-gated tier fix, do not treat as final |
| `plot_sex_specific_lollipop.py` | Plot script — p_het computed correctly, class labels not yet gated by it | Needs fix (see above) |
| `fuma_results/*` | Top10 loci/gene TSVs, dated 2026-03-30/31 | Superseded — see "Previous FUMA runs" above |

## Broader pipeline timeline (reconstructed 2026-08-18)
1. **EDA** (`../../../EDA/`, Dec 2025–Jan 6 2026) — raw liver/spleen attenuation distributions
2. **Phenotyping** (`../../../phenotyping/`, Jan 12–30 2026) — menopause cohort defined, merged with genotype data
3. **Covariates** (`../../../covariates/`, Jan 29 2026) — `hepatic_fat_all_covariates.csv`
4. **Pipeline setup** (`../../../tools/pmbb-nf-toolkit-saige-family`, Feb 9 2026) — Nextflow SAIGE toolkit; git repo formalized Feb 11–12
5. **First GWAS runs** (`../../hepatic_fat_gwas_manual_jan2026` Feb 13, `../../hepatic_fat_gwas_jan2026` Feb 18–22 2026) — ALL-ancestry, sex-stratified SAIGE run; source of the `PMBB_ALL_ALL` sumstats used here
6. **Multi-ancestry expansion** (this directory, `hepatic_fat_gwas_feb2026_sex_stratified/`, Mar 3–Apr 29 2026) — AFR/EUR/ALL × M/F, source of the `PMBB_ALL_F`/`PMBB_ALL_M` sumstats used here
7. **FUMA prep, local** (Mar 16–30 2026) — reformat → hg19 liftover → MAF≥1% filter (see gwas_processed/ provenance above)
8. **Colocalization**, local (May 11–17 2026) — see `../Coloc/CLAUDE.md`

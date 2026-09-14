## Compare cross-ancestry GWAMA meta-analysis lead SNPs against the baseline
## combined-ancestry GWAS (PMBB_ALL_M/F), and pull each candidate locus's
## EUR-only and AFR-only effect size to directly test whether any "new"
## locus is genuinely AFR-driven (i.e. actually recovered, diluted signal)
## rather than an artifact of PMBB_ALL_M/F containing ~350-450 participants
## beyond the EUR+AFR strata (see caveats below).
##
## CAVEATS (read before trusting any "candidate_new" result):
## 1. [RESOLVED 2026-09-11] Baseline lead SNP lists now come from the
##    corrected, full (non-top-10), MAF-filtered FUMA runs
##    (LeadSNPs/fuma_all_sex_stratified/PMBB_sex_stratified_{male,female}_hg19_maf0.01/leadSNPs.txt),
##    generated 2026-03-30 but only pulled/wired in on this date. This
##    replaces the old top-10-only, no-MAF-filter superseded list
##    (fuma_results/PMBB_*_top10_*.tsv), which could not reliably distinguish
##    "genuinely new locus" from "just ranked below #10."
## 2. PMBB_ALL_M/F is NOT simply EUR+AFR pooled -- it has ~350-450 more
##    participants than EUR_M+AFR_M / EUR_F+AFR_F combined (other/unclassified
##    ancestry). Some meta-vs-baseline differences could reflect this
##    population difference, not just the pooled-vs-meta statistical method.
##    The EUR-alone/AFR-alone lookup below is what actually distinguishes a
##    genuine ancestry-specific effect from this confound.

library(data.table)

MERGE_WINDOW_KB <- 250  # matches FUMA's own mergeDist setting used for both runs

setwd("/project/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_feb2026_sex_stratified")

## ---- 1. Load meta lead SNPs ----

meta_lead <- rbindlist(list(
  M = fread("LeadSNPs/fuma_results_meta/eur_afr_m_FUMA_job769143/leadSNPs.txt")[, sex := "M"],
  F = fread("LeadSNPs/fuma_results_meta/eur_afr_f_FUMA_job769145/leadSNPs.txt")[, sex := "F"]
))
setnames(meta_lead, c("chr", "pos", "p", "rsID"), c("chr", "pos_hg19", "meta_p", "meta_rsID"), skip_absent = TRUE)

## ---- 2. Load baseline lead SNPs (corrected, full, MAF-filtered run) ----

base_m <- fread("LeadSNPs/fuma_all_sex_stratified/PMBB_sex_stratified_male_hg19_maf0.01/leadSNPs.txt")
base_m[, sex := "M"]

base_f <- fread("LeadSNPs/fuma_all_sex_stratified/PMBB_sex_stratified_female_hg19_maf0.01/leadSNPs.txt")
base_f[, sex := "F"]

baseline_lead <- rbindlist(list(base_m, base_f), fill = TRUE)
baseline_lead[, chr := as.character(chr)]
meta_lead[, chr := as.character(chr)]

## ---- 3. Classify each meta locus: shared with baseline (within window) vs candidate-new ----

classify_vs_baseline <- function(chr_i, pos_i, sex_i) {
  cand <- baseline_lead[chr == chr_i & sex == sex_i]
  if (nrow(cand) == 0) return(list(status = "candidate_new", baseline_p = NA_real_, baseline_rsid = NA_character_))
  dists <- abs(cand$pos - pos_i)
  hit <- which(dists <= MERGE_WINDOW_KB * 1000)
  if (length(hit) == 0) return(list(status = "candidate_new", baseline_p = NA_real_, baseline_rsid = NA_character_))
  best <- cand[hit][which.min(dists[hit])]
  list(status = "shared", baseline_p = best$p, baseline_rsid = best$rsID)
}

meta_lead[, c("vs_baseline_status", "baseline_p", "baseline_rsID") :=
  {
    res <- Map(classify_vs_baseline, chr, pos_hg19, sex)
    list(sapply(res, `[[`, "status"),
         sapply(res, `[[`, "baseline_p"),
         sapply(res, `[[`, "baseline_rsid"))
  }]

cat(sprintf("Meta loci: %d total (%d shared with baseline top10, %d candidate_new)\n",
            nrow(meta_lead), sum(meta_lead$vs_baseline_status == "shared"),
            sum(meta_lead$vs_baseline_status == "candidate_new")))

## ---- 4. Get the hg38 variant_id for each meta lead SNP (exact hg19 position match
##         against the meta's own hg19.gz file, which still carries the original
##         hg38-based variant_id string) ----

get_variant_ids <- function(sex_i) {
  f <- sprintf("Meta/PMBB_META_EUR_AFR_%s.hepatic_fat.gwama.fuma.hg19.gz", sex_i)
  d <- fread(cmd = sprintf("zcat %s", f),
             select = c("chromosome", "base_pair_location", "variant_id", "effect_allele_frequency"))
  d[, chromosome := as.character(chromosome)]
  d[, meta_maf := pmin(effect_allele_frequency, 1 - effect_allele_frequency)]
  d
}
lookup_m <- get_variant_ids("M")
lookup_f <- get_variant_ids("F")

meta_lookup_result <- rbindlist(lapply(seq_len(nrow(meta_lead)), function(i) {
  lu <- if (meta_lead$sex[i] == "M") lookup_m else lookup_f
  hit <- lu[chromosome == meta_lead$chr[i] & base_pair_location == meta_lead$pos_hg19[i]]
  if (nrow(hit) > 0) hit[1, .(variant_id, meta_maf)] else data.table(variant_id = NA_character_, meta_maf = NA_real_)
}))
meta_lead[, variant_id := meta_lookup_result$variant_id]
meta_lead[, meta_maf := meta_lookup_result$meta_maf]

## ---- 5. Pull EUR-alone and AFR-alone beta/SE/p for each locus's variant_id ----
## NOTE: GWAMA aligns effect alleles across studies, which can flip a
## variant's allele order (ref/alt) relative to how it's written in the raw
## per-ancestry sumstats -- an exact string match on variant_id alone misses
## these. We look up BOTH the as-is id and its allele-swapped counterpart,
## and negate beta when the swapped id is what matched (since that row's
## beta is relative to the opposite allele).

parse_variant <- function(vid) {
  parts <- tstrsplit(vid, "_")
  list(chr = parts[[1]], pos = parts[[2]], a1 = parts[[3]], a2 = parts[[4]])
}

variant_ids_needed <- unique(na.omit(meta_lead$variant_id))
vp <- parse_variant(variant_ids_needed)
flipped_ids <- sprintf("%s_%s_%s_%s", vp$chr, vp$pos, vp$a2, vp$a1)
all_ids_to_search <- unique(c(variant_ids_needed, flipped_ids))
writeLines(all_ids_to_search, "LeadSNPs/tmp_variant_ids_needed.txt")

id_map <- data.table(variant_id = variant_ids_needed, flipped_id = flipped_ids)

extract_ancestry_stats <- function(ancestry_sex) {
  infile <- sprintf("PMBB_%s/Sumstats/PMBB_%s.hepatic_fat.gwas.saige.gz", ancestry_sex, ancestry_sex)
  cmd <- sprintf(
    "zcat %s | awk -F'\\t' 'NR==FNR{ids[$1]=1; next} FNR==1 || ($4 in ids)' LeadSNPs/tmp_variant_ids_needed.txt -",
    infile
  )
  raw <- fread(cmd = cmd)[, .(variant_id, beta, standard_error, p_value, effect_allele_frequency)]
  raw[, maf := pmin(effect_allele_frequency, 1 - effect_allele_frequency)]

  ## resolve back to the ORIGINAL (as-in-meta) variant_id, negating beta
  ## whenever only the flipped-allele-order id was found (SE and MAF are
  ## unaffected by which allele is "effect")
  direct <- merge(id_map[, .(variant_id)], raw, by = "variant_id")
  via_flip <- merge(id_map, raw, by.x = "flipped_id", by.y = "variant_id")
  if (nrow(via_flip) > 0) {
    via_flip[, beta := -beta]
    via_flip <- via_flip[, .(variant_id, beta, standard_error, p_value, maf)]
    ## don't double-count a locus already resolved directly
    via_flip <- via_flip[!variant_id %in% direct$variant_id]
  }
  rbindlist(list(direct[, .(variant_id, beta, standard_error, p_value, maf)], via_flip), fill = TRUE)
}

eur_m <- extract_ancestry_stats("EUR_M")[, .(variant_id, EUR_beta = beta, EUR_se = standard_error, EUR_p = p_value, EUR_maf = maf)]
afr_m <- extract_ancestry_stats("AFR_M")[, .(variant_id, AFR_beta = beta, AFR_se = standard_error, AFR_p = p_value, AFR_maf = maf)]
eur_f <- extract_ancestry_stats("EUR_F")[, .(variant_id, EUR_beta = beta, EUR_se = standard_error, EUR_p = p_value, EUR_maf = maf)]
afr_f <- extract_ancestry_stats("AFR_F")[, .(variant_id, AFR_beta = beta, AFR_se = standard_error, AFR_p = p_value, AFR_maf = maf)]

meta_m <- merge(meta_lead[sex == "M"], eur_m, by = "variant_id", all.x = TRUE)
meta_m <- merge(meta_m, afr_m, by = "variant_id", all.x = TRUE)
meta_f <- merge(meta_lead[sex == "F"], eur_f, by = "variant_id", all.x = TRUE)
meta_f <- merge(meta_f, afr_f, by = "variant_id", all.x = TRUE)

result <- rbindlist(list(meta_m, meta_f), fill = TRUE)

## ---- 5b. Pull n_studies/effects from the raw GWAMA .out to distinguish
##          "AFR-only tested" (n_studies=1, invisible to any EUR-inclusive
##          analysis by definition) from "tested in both, AFR-driven"
##          (n_studies=2, real dilution-in-pooled-regression case) ----

get_n_studies <- function(sex_i) {
  f <- sprintf("Meta/PMBB_META_EUR_AFR_%s_hepatic_fat.out", sex_i)
  fread(f, select = c("rs_number", "n_studies", "effects"))
}
nstud_m <- get_n_studies("M")
nstud_f <- get_n_studies("F")
nstud <- rbindlist(list(nstud_m, nstud_f))
setnames(nstud, "rs_number", "variant_id")

result <- merge(result, unique(nstud, by = "variant_id"), by = "variant_id", all.x = TRUE)
result[, ancestry_category := fifelse(
  n_studies == 1 & !is.na(AFR_p) & is.na(EUR_p), "AFR-only tested (invisible to EUR-inclusive analysis)",
  fifelse(n_studies == 1 & !is.na(EUR_p) & is.na(AFR_p), "EUR-only tested",
  fifelse(n_studies == 2, "tested in both ancestries", "unknown")))]

setorder(result, sex, meta_p)

## ---- 6. Save + summarize ----

fwrite(result, "LeadSNPs/meta_vs_baseline_ancestry_comparison.csv")
file.remove("LeadSNPs/tmp_variant_ids_needed.txt")

cat("\n=== Candidate-new loci (not in the full corrected baseline lead SNP list) ===\n")
print(result[vs_baseline_status == "candidate_new",
             .(sex, chr, pos_hg19, meta_rsID, meta_p, ancestry_category, EUR_beta, EUR_p, AFR_beta, AFR_p)])

cat("\n=== Breakdown of candidate-new loci by ancestry category ===\n")
print(result[vs_baseline_status == "candidate_new", .N, by = .(sex, ancestry_category)])

cat("\nSaved full comparison table to LeadSNPs/meta_vs_baseline_ancestry_comparison.csv\n")

## ---- 7. Presentation-ready trimmed table: AFR-only + tested-in-both loci,
##         with mapped gene symbol and per-ancestry MAF ----

genes_m <- fread("LeadSNPs/fuma_results_meta/eur_afr_m_FUMA_job769143/genes.txt")[, sex := "M"]
genes_f <- fread("LeadSNPs/fuma_results_meta/eur_afr_f_FUMA_job769145/genes.txt")[, sex := "F"]
genes <- rbindlist(list(genes_m, genes_f), fill = TRUE)
gene_by_locus <- genes[, .(gene = paste(unique(symbol), collapse = ";")), by = .(sex, GenomicLocus)]

result <- merge(result, gene_by_locus, by.x = c("sex", "GenomicLocus"), by.y = c("sex", "GenomicLocus"), all.x = TRUE)

## IMPORTANT: must also require vs_baseline_status == "candidate_new" here --
## PMBB_ALL_M/F is a COMBINED-ancestry cohort that includes AFR participants,
## so a variant absent from EUR-only data can still appear in the baseline's
## own lead SNP list via its AFR subset. Without this filter, "AFR-only
## tested" rows that are already present in baseline (not novel) get mixed
## in with genuinely new ones -- caught 2026-09-11 when a manual row count
## (22) didn't match the script's own reported breakdown count (18).
presentation_table <- result[
  vs_baseline_status == "candidate_new" &
  ancestry_category %in% c("AFR-only tested (invisible to EUR-inclusive analysis)", "tested in both ancestries"),
  .(sex, chr, pos_hg19, gene, meta_rsID, meta_p, meta_maf, ancestry_category,
    EUR_beta, EUR_se, EUR_p, EUR_maf, AFR_beta, AFR_se, AFR_p, AFR_maf)
]
setorder(presentation_table, sex, meta_p)
fwrite(presentation_table, "LeadSNPs/meta_vs_baseline_presentation_table.csv")

cat("\n=== Presentation table (AFR-only + tested-in-both, with gene + MAF) ===\n")
print(presentation_table)
cat("\nSaved to LeadSNPs/meta_vs_baseline_presentation_table.csv\n")

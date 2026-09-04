#!/usr/bin/env Rscript
# coloc_hepatic_fat_gtex.R
# Formal colocalization: PMBB hepatic fat GWAS × GTEx v8 eQTL
# Analyses: male-only, female-only, ALL combined
# Method: coloc.abf (Wakefield approximate Bayes factors; single causal variant)
# Tissues: Liver, Adipose Visceral Omentum, Adipose Subcutaneous, Whole Blood
#
# Install once:
#   install.packages(c("coloc", "httr2", "data.table"))
#
# VARIANT MATCHING NOTE:
#   GWAS base_pair_location is hg19 (post-liftover), but variant_id retains the
#   original SAIGE hg38 encoding (e.g., chr22_43556374_G_A). GTEx API returns
#   variantId in the same hg38 format (chr22_43556374_G_A_b38). After stripping
#   the chr prefix and _b38 suffix both resolve to 22_43556374_G_A — no liftover
#   needed for variant matching. Positional window filtering still uses hg19.
#
# API COMPLETENESS NOTE:
#   The GTEx v2 singleTissueEqtl endpoint should return all nominal cis-eQTL
#   associations (all tested variants, not just significant ones). If <100
#   variants are returned per gene×tissue, the API may be filtering to
#   significant hits, which biases coloc. In that case switch to the GTEx bulk
#   "all_associations" files available at gtexportal.org/home/datasets.

suppressPackageStartupMessages({
  library(coloc)
  library(httr2)
  library(data.table)
})

setwd("/Users/agaro/Documents/meno_effects_liver")

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ══ CONFIG ════════════════════════════════════════════════════════════════════

WINDOW_BP     <- 500000L   # ±500 kb around each lead SNP (hg19)
GTEX_BASE     <- "https://gtexportal.org/api/v2"
DATASET_ID    <- "gtex_v8"
PP4_THRESHOLD <- 0.8       # PP.H4 threshold for reporting colocalization
MIN_SNPS      <- 50L       # minimum shared variants required to run coloc
MAX_PER_PAGE  <- 2000L     # GTEx API items per page (increase if API allows)

# GWAS files and sample sizes
GWAS_META <- list(
  male   = list(file = "gwas_processed/PMBB_ALL_M.hepatic_fat.gwas.saige.fuma.maf01.hg19.gz",   N = 9048L),
  female = list(file = "gwas_processed/PMBB_ALL_F.hepatic_fat.gwas.saige.fuma.maf01.hg19.gz",   N = 9031L),
  all    = list(file = "gwas_processed/PMBB_ALL_ALL.hepatic_fat.gwas.saige.fuma.maf01.hg19.gz", N = 18079L)
)

# Lead SNP file configs — UPDATE these paths when new FUMA runs complete
LEAD_CFG <- list(
  male = list(
    file     = "fuma_results/PMBB_male_top10_leadSNPS_FUMA.tsv",
    chr      = "chr",          pos      = "pos",
    rsid     = "rsID",         ngene    = "nearest_gene",
    eqtl_col = "eQTL_gene"
  ),
  female = list(
    file     = "fuma_results/PMBB_female_top10_loci_FUMA.tsv",
    chr      = "Chr",          pos      = "Pos_hg19",
    rsid     = "rsID",         ngene    = "Nearest_Gene",
    eqtl_col = "eQTL_Genes"
  ),
  all = list(
    file     = "fuma_results/PMBB_ALL_top10_leadSNPs_FUMA.tsv",
    chr      = "chr",          pos      = "pos",
    rsid     = "rsID",         ngene    = "nearest_gene",
    eqtl_col = "eQTL_gene"
  )
)

# GTEx v8 tissue IDs and approximate sample sizes (verify at gtexportal.org)
TISSUES <- data.frame(
  id = c("Liver", "Adipose_Visceral_Omentum", "Adipose_Subcutaneous",
         "Artery_Aorta", "Artery_Tibial"),
  N  = c(208L,    581L,                        581L,
         387L,          584L),
  stringsAsFactors = FALSE
)

# ══ VARIANT KEY NORMALIZATION ═════════════════════════════════════════════════

# Normalize to bare chr_pos_ref_alt key for matching across GWAS and GTEx
# "chr22_43556374_G_A"      → "22_43556374_G_A"
# "chr22_43556374_G_A_b38"  → "22_43556374_G_A"
norm_vid <- function(v) {
  v <- sub("^chr", "", v)
  v <- sub("_(b38|b37)$", "", v)
  v
}

# Flip ref/alt alleles in a normalized key (handles strand ambiguity)
# "22_43556374_G_A" → "22_43556374_A_G"
flip_vid <- function(v) {
  parts <- strsplit(v, "_", fixed = TRUE)[[1]]
  if (length(parts) == 4L) paste(parts[1], parts[2], parts[4], parts[3], sep = "_") else v
}

# ══ LOAD LEAD SNPS ════════════════════════════════════════════════════════════

load_leads <- function(sex) {
  cfg <- LEAD_CFG[[sex]]
  dt  <- fread(cfg$file, sep = "\t")
  raw <- paste(dt[[cfg$ngene]], dt[[cfg$eqtl_col]], sep = ";")
  gene_lists <- lapply(raw, function(g) {
    syms <- unique(trimws(unlist(strsplit(g, "[;,]"))))
    syms[nchar(syms) > 0L & syms != "—" & syms != "NA" & syms != "—"]
  })
  data.table(
    sex       = sex,
    rsid      = dt[[cfg$rsid]],
    chr       = as.integer(dt[[cfg$chr]]),
    pos       = as.integer(dt[[cfg$pos]]),
    gene_list = gene_lists
  )
}

leads <- rbindlist(lapply(c("male", "female", "all"), load_leads), fill = TRUE)
cat(sprintf("Loaded %d lead SNP entries across 3 GWAS\n", nrow(leads)))

# ══ GTEx API FUNCTIONS ════════════════════════════════════════════════════════

gencode_cache <- list()

get_gencode_id <- function(symbol) {
  if (!is.null(gencode_cache[[symbol]])) return(gencode_cache[[symbol]])
  resp <- tryCatch(
    request(paste0(GTEX_BASE, "/reference/gene")) |>
      req_url_query(geneSymbol = symbol, referenceGenomeId = "GRCh38/hg38") |>
      req_retry(max_tries = 3) |>
      req_perform(),
    error = function(e) NULL
  )
  gid <- NA_character_
  if (!is.null(resp)) {
    genes <- resp_body_json(resp)$gene
    if (length(genes) > 0L) gid <- genes[[1L]]$gencodeId %||% NA_character_
  }
  gencode_cache[[symbol]] <<- gid
  Sys.sleep(0.1)
  gid
}

query_gtex_eqtl <- function(gencode_id, tissue_id) {
  all_rows <- list()
  page     <- 0L
  repeat {
    resp <- tryCatch(
      request(paste0(GTEX_BASE, "/association/singleTissueEqtl")) |>
        req_url_query(
          gencodeId          = gencode_id,
          tissueSiteDetailId = tissue_id,
          datasetId          = DATASET_ID,
          itemsPerPage       = MAX_PER_PAGE,
          page               = page
        ) |>
        req_retry(max_tries = 3) |>
        req_perform(),
      error = function(e) NULL
    )
    if (is.null(resp)) break
    body    <- resp_body_json(resp)
    records <- body$data
    if (!length(records)) break
    all_rows <- c(all_rows, records)
    if (page >= (body$numPages %||% 1L) - 1L) break
    page <- page + 1L
    Sys.sleep(0.15)
  }
  if (!length(all_rows)) return(NULL)

  dt <- rbindlist(lapply(all_rows, function(r) {
    data.table(
      vid  = norm_vid(r$variantId %||% ""),
      beta = r$nes    %||% NA_real_,
      se   = r$se     %||% NA_real_,
      maf  = r$maf    %||% NA_real_
    )
  }), fill = TRUE)

  dt <- dt[nchar(vid) > 0L & !is.na(beta) & !is.na(se) & se > 0 & !is.na(maf) & maf > 0]
  dt[dt[, .I[1L], by = vid]$V1]  # deduplicate
}

# ══ GWAS EXTRACTION ══════════════════════════════════════════════════════════

gwas_cache <- list()

extract_gwas <- function(sex, chrom, center_pos) {
  cache_key <- paste(sex, chrom, center_pos, sep = "_")
  if (!is.null(gwas_cache[[cache_key]])) return(gwas_cache[[cache_key]])

  lo  <- center_pos - WINDOW_BP
  hi  <- center_pos + WINDOW_BP
  cmd <- sprintf(
    "gzcat '%s' | awk 'NR==1 || ($1==%d && $2>=%d && $2<=%d)'",
    GWAS_META[[sex]]$file, chrom, lo, hi
  )
  dt <- tryCatch(
    fread(cmd = cmd, sep = "\t", showProgress = FALSE),
    error = function(e) { message(sprintf("GWAS read error: %s", conditionMessage(e))); NULL }
  )
  if (is.null(dt) || nrow(dt) == 0L) { gwas_cache[[cache_key]] <<- NULL; return(NULL) }

  dt[, vid     := norm_vid(variant_id)]
  dt[, maf_gwas := pmin(effect_allele_frequency, 1 - effect_allele_frequency)]
  dt <- dt[nchar(vid) > 0L & !is.na(beta) & !is.na(standard_error) & standard_error > 0 &
             maf_gwas > 0 & maf_gwas < 1]
  dt <- dt[dt[, .I[which.min(p_value)], by = vid]$V1]

  gwas_cache[[cache_key]] <<- dt
  dt
}

# ══ COLOC ════════════════════════════════════════════════════════════════════

run_coloc <- function(gwas_dt, eqtl_dt, gwas_n, eqtl_n) {
  # Primary match on exact normalized hg38 key
  shared <- intersect(gwas_dt$vid, eqtl_dt$vid)

  # Extend with allele-flipped matches (handles strand differences)
  eqtl_flipped        <- copy(eqtl_dt)
  eqtl_flipped[, vid  := sapply(vid, flip_vid)]
  eqtl_flipped[, beta := -beta]
  extra_shared <- intersect(gwas_dt$vid, eqtl_flipped$vid)
  extra_shared <- setdiff(extra_shared, shared)
  if (length(extra_shared)) {
    eqtl_dt <- rbindlist(list(eqtl_dt, eqtl_flipped[vid %in% extra_shared]))
    shared   <- c(shared, extra_shared)
  }

  if (length(shared) < MIN_SNPS) return(NULL)

  gw <- gwas_dt[vid %in% shared, .(vid, beta,   se = standard_error, maf = maf_gwas)]
  eq <- eqtl_dt[vid %in% shared, .(vid, beta_e = beta, se_e = se, maf_e = maf)]
  m  <- merge(gw, eq, by = "vid")
  m  <- m[!is.na(beta) & !is.na(beta_e)]
  if (nrow(m) < MIN_SNPS) return(NULL)

  d1 <- list(snp = m$vid, beta = m$beta,   varbeta = m$se^2,
             N = gwas_n, MAF = m$maf,   type = "quant")
  d2 <- list(snp = m$vid, beta = m$beta_e, varbeta = m$se_e^2,
             N = eqtl_n, MAF = m$maf_e, type = "quant")

  res <- tryCatch(suppressMessages(coloc.abf(d1, d2)), error = function(e) NULL)
  if (is.null(res)) return(NULL)

  pp <- as.list(res$summary)
  list(nsnps = nrow(m),
       PP.H0 = pp$PP.H0.abf, PP.H1 = pp$PP.H1.abf,
       PP.H2 = pp$PP.H2.abf, PP.H3 = pp$PP.H3.abf,
       PP.H4 = pp$PP.H4.abf)
}

# ══ MAIN LOOP ════════════════════════════════════════════════════════════════

results <- list()

for (i in seq_len(nrow(leads))) {
  row       <- leads[i]
  sex       <- row$sex
  rsid      <- row$rsid
  chrom     <- row$chr
  pos       <- row$pos
  gene_syms <- row$gene_list[[1L]]

  cat(sprintf("\n[%d/%d] %s | %s | chr%d:%d | genes: %s\n",
              i, nrow(leads), sex, rsid, chrom, pos,
              paste(gene_syms, collapse = ", ")))

  gwas_dt <- extract_gwas(sex, chrom, pos)
  if (is.null(gwas_dt)) { cat("  [skip] no GWAS variants in window\n"); next }
  cat(sprintf("  GWAS variants in window: %d\n", nrow(gwas_dt)))

  for (gene_sym in gene_syms) {
    gencode_id <- get_gencode_id(gene_sym)
    if (is.na(gencode_id)) {
      cat(sprintf("  [skip] %s: gencodeId not found in GTEx\n", gene_sym))
      next
    }

    for (ti in seq_len(nrow(TISSUES))) {
      tissue_id <- TISSUES$id[ti]
      tissue_n  <- TISSUES$N[ti]

      cat(sprintf("  querying %-15s × %-30s ... ", gene_sym, tissue_id))
      eqtl_dt <- query_gtex_eqtl(gencode_id, tissue_id)

      if (is.null(eqtl_dt)) { cat("no data\n"); next }
      cat(sprintf("%d variants", nrow(eqtl_dt)))

      if (nrow(eqtl_dt) < 100L) {
        cat(" [WARNING: very few variants — API may return significant-only; consider GTEx bulk files]")
      }
      cat("\n")

      res <- run_coloc(gwas_dt, eqtl_dt, GWAS_META[[sex]]$N, tissue_n)
      if (!is.null(res)) {
        results[[length(results) + 1L]] <- data.table(
          gwas       = sex,
          lead_rsid  = rsid,
          chr        = chrom,
          pos_hg19   = pos,
          gene       = gene_sym,
          gencode_id = gencode_id,
          tissue     = tissue_id,
          n_snps     = res$nsnps,
          PP.H0      = round(res$PP.H0, 4),
          PP.H1      = round(res$PP.H1, 4),
          PP.H2      = round(res$PP.H2, 4),
          PP.H3      = round(res$PP.H3, 4),
          PP.H4      = round(res$PP.H4, 4),
          coloc      = res$PP.H4 >= PP4_THRESHOLD
        )
      }
    }
  }
}

# ══ OUTPUT ════════════════════════════════════════════════════════════════════

if (!length(results)) {
  cat("\nNo colocalization results produced (check warnings above).\n")
} else {
  out <- rbindlist(results)
  out <- out[order(-PP.H4)]
  fwrite(out, "coloc_results.tsv", sep = "\t")
  cat(sprintf("\nResults written → coloc_results.tsv (%d rows)\n", nrow(out)))

  hits <- out[coloc == TRUE]
  cat(sprintf(
    "\n── Colocalized signals (PP.H4 ≥ %.1f): %d ───────────────────────\n",
    PP4_THRESHOLD, nrow(hits)
  ))
  if (nrow(hits) > 0L) {
    print(hits[, .(gwas, lead_rsid, gene, tissue, n_snps,
                   PP.H3 = round(PP.H3, 3), PP.H4 = round(PP.H4, 3))])
  }
  cat(sprintf(
    "\nAll PP.H4 values (top 20):\n"
  ))
  print(out[1:min(20L, nrow(out)),
            .(gwas, lead_rsid, gene, tissue, n_snps, PP.H4)])
}

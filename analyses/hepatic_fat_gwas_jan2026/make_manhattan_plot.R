library(data.table)
library(qqman)
library(dplyr)

# Paths
results_file <- "/home/agaro/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_manual_jan2026/PMBB_ALL_ALL/Sumstats/PMBB_ALL_ALL.hepatic_fat.gwas.saige.gz"
annot_file <- "/home/agaro/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_jan2026/Annotations/GWAS_biofilter_genes_rsids.csv"
output_plot <- "/home/agaro/verma_shared/projects/Liver_IDPs/analyses/hepatic_fat_gwas_manual_jan2026/Plots/hepatic_fat_manhattan_qqman.png"

# Load data
cat("Loading GWAS results...\n")
gwas <- fread(results_file)

# Load annotations
cat("Loading gene annotations...\n")
annot <- fread(annot_file)

# Merge
gwas <- merge(gwas, annot[, .(Var_ID, Gene, RSID)], 
              by.x = "variant_id", by.y = "Var_ID", 
              all.x = TRUE)

# Clean chromosome names
gwas$CHR <- as.integer(gsub("chr", "", gwas$chromosome))
gwas <- gwas[!is.na(CHR) & CHR %in% 1:22]

# Rename columns for qqman
gwas$BP <- gwas$base_pair_location
gwas$P <- gwas$p_value
gwas$SNP <- gwas$variant_id

# Remove missing p-values
gwas <- gwas[!is.na(P) & P > 0]

cat(sprintf("Total variants: %s\n", format(nrow(gwas), big.mark=",")))

# Separate significant and non-significant
sig_threshold <- 5e-8
sig <- gwas[P < sig_threshold]
nonsig <- gwas[P >= sig_threshold]

cat(sprintf("Significant variants: %s\n", nrow(sig)))
cat(sprintf("Non-significant variants: %s\n", nrow(nonsig)))

# Clump only significant variants (to avoid PNPLA3 repetition)
sig <- sig[order(CHR, BP, P)]
sig$region <- paste0(sig$CHR, "_", floor(sig$BP / 1e6))
sig_clumped <- sig %>%
  group_by(region) %>%
  slice_min(P, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  as.data.table()

cat(sprintf("Clumped significant to: %s independent signals\n", nrow(sig_clumped)))

# Thin non-significant variants for plotting speed (keep random 50k)
set.seed(42)
if (nrow(nonsig) > 50000) {
  nonsig_thinned <- nonsig[sample(.N, 50000)]
} else {
  nonsig_thinned <- nonsig
}

# Combine clumped significant + thinned non-significant
gwas_plot <- rbind(sig_clumped, nonsig_thinned, fill = TRUE)
gwas_plot <- gwas_plot[order(CHR, BP)]

# Add gene names for labeling
gwas_plot$LABEL <- ifelse(!is.na(gwas_plot$Gene) & gwas_plot$P < sig_threshold, 
                          gwas_plot$Gene, 
                          NA)

cat(sprintf("Total variants for plotting: %s\n", nrow(gwas_plot)))

# Create Manhattan plot
png(output_plot, width = 1800, height = 600, res = 100)

manhattan(gwas_plot,
         chr = "CHR",
         bp = "BP", 
         p = "P",
         snp = "LABEL",
         col = c("#4169E1", "#FFA500"),
         suggestiveline = -log10(1e-5),
         genomewideline = -log10(5e-8),
         annotatePval = 5e-8,
         annotateTop = TRUE,
         main = "GWAS Manhattan Plot: Hepatic Fat (N = 18,079)",
         ylim = c(0, max(-log10(gwas_plot$P)) * 1.05))

dev.off()

cat(sprintf("\nManhattan plot saved to: %s\n", output_plot))
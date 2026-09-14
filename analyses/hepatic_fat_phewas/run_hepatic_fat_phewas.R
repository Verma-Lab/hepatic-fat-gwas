## Hepatic Fat PheWAS (Aim 2D)
## Exposure: CT-derived hepatic fat deltaHU (continuous, from hepatic_fat_all_covariates.csv)
## Outcome: Phecode 1.2 case/control status, built from PMBB-Release-2026-4.0 condition occurrences
## Covariates: Age, Sex, PC1-6 
##
## Data provenance:
## - cohort_phecode12_occurrences.tsv: pre-filtered from
##   /static/PMBB/PMBB-Release-2026-4.0/Phenotype/4.0/PMBB-Release-2026-4.0_phenotype_conditions_phecode_12.txt
##   (99.9% coverage of this cohort)
## - condition_source_value in that release file is already phecode-mapped (not raw ICD)

library(data.table)
library(PheWAS)

MIN_CASES <- 200  # per-lab precedent (Sudha's imaging-guided PheWAS pipeline), stricter than the commonly-used >=50

## ---- 1. Load inputs ----

occ <- fread("cohort_phecode12_occurrences.tsv", header = FALSE,
             col.names = c("PMBB_ID", "phecode"),
             colClasses = c("character", "character"))
# NOTE: colClasses is required here -- fread auto-infers phecode as numeric since every
# value (e.g. "773") happens to be numeric-parseable. Left as numeric, `pc` in the loop
# below becomes a bare number, so `wide[[pc]] <- ...` gets read as *positional* indexing
# (e.g. wide[[773]], "the 773rd column") instead of assigning a new column by name --
# this caused an "[[ ]] subscript out of bounds" failure when this was first run.

covars <- fread("../../covariates/hepatic_fat_all_covariates.csv")
# Expected columns: PMBB_ID, Sample_age, hepatic_fat, PC1-PC6, Sex

## ---- 2. Build wide phecode case/control matrix (with proper exclusion-range logic) ----
## Case = >=1 occurrence of the tested phecode.
## Excluded (NA, neither case nor control) = no occurrence of the tested phecode itself, but
## >=1 occurrence of a phecode listed in its exclusion range (PheWAS::phecode_exclude) --
## e.g. someone with Type 1 diabetes is excluded from the control pool when testing Type 2
## diabetes, rather than being counted as a "clean" control.
## Control = everyone else.
##
## NOTE on implementation: PheWAS::mapPhecodesToExclusions() is the package's built-in
## function for this, but it's broken in this environment -- it calls the removed
## dplyr::tbl_df() (dplyr 1.2.1 here; tbl_df() was deprecated in dplyr 1.0.0 and has since
## been fully removed). Reimplemented the same logic manually below directly against the
## PheWAS::phecode_exclude reference table instead. Validated against a small synthetic
## example (case/excluded/control all resolve as expected) before use here.

occ <- occ[phecode != ""]  # drop condition records that didn't map to a valid Phecode 1.2 entry
                            # (646,672 of the pre-filtered rows, across 16,353 of 18,079
                            # participants, had a blank condition_source_value in the release
                            # file -- left in, these would get tested as a bogus "" phecode)
occ <- unique(occ)  # guard against duplicate (id, phecode) rows from multiple diagnosis events
setkey(occ, phecode)

case_counts <- occ[, .(n_cases = uniqueN(PMBB_ID)), by = phecode]
test_phecodes <- case_counts[n_cases >= MIN_CASES, phecode]
cat(sprintf("Phecodes with >=%d cases: %d (of %d total observed)\n",
            MIN_CASES, length(test_phecodes), nrow(case_counts)))

excl <- as.data.table(phecode_exclude)  # columns: code, exclusion_criteria

all_ids <- unique(covars$PMBB_ID)
wide <- data.table(PMBB_ID = all_ids)

for (pc in test_phecodes) {
  case_ids <- unique(occ[.(pc), PMBB_ID])
  excl_codes <- excl[code == pc, exclusion_criteria]
  excl_ids <- if (length(excl_codes) > 0) {
    setdiff(unique(occ[.(excl_codes), on = "phecode", nomatch = NULL, PMBB_ID]), case_ids)
  } else {
    character(0)
  }
  wide[[pc]] <- ifelse(all_ids %in% case_ids, TRUE,
                 ifelse(all_ids %in% excl_ids, NA, FALSE))
  # NOTE: must be logical (TRUE/FALSE/NA), not integer (1L/0L/NA_integer_) -- phewas()
  # only runs logistic regression (and only populates n_cases/n_controls) for logical
  # outcome columns; an integer column gets silently treated as continuous, running
  # linear regression instead. Confirmed via the `type` column in a first real run:
  # all 944 completed tests came back "linear" before this fix.
}

## ---- 3. Merge phenotype + covariates ----

analysis_df <- merge(covars, wide, by = "PMBB_ID")
cat(sprintf("Final analysis N: %d (of %d cohort participants)\n",
            nrow(analysis_df), length(all_ids)))

phecode_cols <- setdiff(names(wide), "PMBB_ID")

## ---- 4. Run PheWAS ----
## `predictors` here is the continuous exposure (hepatic fat deltaHU), not a genotype --
## PheWAS::phewas() treats it generically as a predictor of interest regardless.

results <- phewas(
  phenotypes = phecode_cols,
  predictors = "hepatic_fat",
  data = analysis_df,
  covariates = c("Sample_age", "Sex", paste0("PC", 1:6)),
  cores = 1,
  additive.genotypes = FALSE  # predictor is a continuous IDP, not a genotype -- avoids
                               # the "Genotype is not coded 0,1" warning on every row and
                               # skips the (inapplicable) HWE/allele-frequency attempt
)

results <- addPhecodeInfo(results)
setDT(results)  # phewas() and addPhecodeInfo() both return plain data.frames -- the
                # `:=` syntax below requires an actual data.table, and failed with
                # "Check that is.data.table(DT) == TRUE" without this conversion.

## ---- 5. Multiple testing + output ----

n_tested <- length(phecode_cols)
bonferroni_p <- 0.05 / n_tested
results[, bonferroni_threshold := bonferroni_p]
results[, significant := p < bonferroni_threshold]

setorder(results, p)

fwrite(results, "hepatic_fat_phewas_results.csv")

cat(sprintf("\nTested %d phecodes. Bonferroni threshold: %.2e\n", n_tested, bonferroni_p))
cat(sprintf("Significant associations: %d\n", sum(results$significant, na.rm = TRUE)))

# ============================================================
# GSE213001 W6.6 — Biological interpretation preflight
#
# Purpose:
#   Validate the frozen W5/W6 inputs required for downstream
#   gene prioritisation and enrichment before inspecting
#   individual gene identities.
#
# No biological interpretation is performed.
# No individual gene identities are printed.
# ============================================================

source("renv/activate.R")

# ------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------

de <- read.csv(
  "results/differential_expression/primary_IPF_vs_NDC_annotated.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("")
)

age <- read.csv(
  "results/differential_expression/sensitivity_age_gene_concordance.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("\n============================================================\n")
cat("GSE213001 W6.6 — BIOLOGICAL INTERPRETATION PREFLIGHT\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 2. Required-column audit
# ------------------------------------------------------------

required_de <- c(
  "gene_id",
  "logFC",
  "AveExpr",
  "t",
  "P.Value",
  "adj.P.Val",
  "B",
  "SYMBOL",
  "ENTREZID",
  "GENENAME",
  "n_symbol",
  "n_entrezid",
  "n_genename",
  "n_mapping_records",
  "mapping_status"
)

required_age <- c(
  "gene_id",
  "logFC_reduced",
  "logFC_age_adjusted",
  "delta_logFC_age_minus_reduced",
  "abs_delta_logFC",
  "direction_concordant",
  "FDR_reduced",
  "FDR_age_adjusted",
  "significant_reduced",
  "significant_age_adjusted"
)

stopifnot(
  all(required_de %in% names(de)),
  all(required_age %in% names(age))
)

cat("\n=== REQUIRED COLUMNS ===\n")
cat("Primary DE columns present:     TRUE\n")
cat("Age-sensitivity columns present: TRUE\n")

# ------------------------------------------------------------
# 3. Frozen universe invariants
# ------------------------------------------------------------

stopifnot(
  nrow(de) == 15012,
  nrow(age) == 15012,
  length(unique(de$gene_id)) == 15012,
  length(unique(age$gene_id)) == 15012,
  setequal(de$gene_id, age$gene_id)
)

cat("\n=== FROZEN GENE UNIVERSE ===\n")
cat("Primary DE genes:        ", nrow(de), "\n", sep = "")
cat("Age-sensitivity genes:   ", nrow(age), "\n", sep = "")
cat("Same gene universe:      TRUE\n")

# ------------------------------------------------------------
# 4. Primary statistical completeness
# ------------------------------------------------------------

stat_cols <- c(
  "logFC",
  "AveExpr",
  "t",
  "P.Value",
  "adj.P.Val",
  "B"
)

stat_missing <- vapply(
  de[stat_cols],
  function(x) sum(is.na(x)),
  integer(1)
)

stopifnot(
  all(stat_missing == 0L),
  all(de$P.Value >= 0 & de$P.Value <= 1),
  all(de$adj.P.Val >= 0 & de$adj.P.Val <= 1)
)

direction_agreement <- all(
  sign(de$logFC[de$logFC != 0]) ==
    sign(de$t[de$logFC != 0])
)

stopifnot(direction_agreement)

cat("\n=== PRIMARY STATISTICAL INTEGRITY ===\n")
cat("Missing statistical values:  ", sum(stat_missing), "\n", sep = "")
cat("P-values within [0,1]:       TRUE\n")
cat("FDR values within [0,1]:     TRUE\n")
cat("logFC / t sign agreement:    TRUE\n")

# ------------------------------------------------------------
# 5. Frozen descriptive priority counts
# ------------------------------------------------------------

primary_sig <- de$adj.P.Val < 0.05
large_effect <- abs(de$logFC) >= 1
priority <- primary_sig & large_effect

cat("\n=== PRIMARY DESCRIPTIVE SETS ===\n")
cat("FDR < 0.05:                     ", sum(primary_sig), "\n", sep = "")
cat("|logFC| >= 1:                   ", sum(large_effect), "\n", sep = "")
cat("FDR < 0.05 and |logFC| >= 1:    ", sum(priority), "\n", sep = "")
cat(
  "  higher in IPF:                ",
  sum(priority & de$logFC > 0),
  "\n",
  sep = ""
)
cat(
  "  lower in IPF:                 ",
  sum(priority & de$logFC < 0),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 6. Age-sensitivity integrity
# ------------------------------------------------------------

idx_age <- match(de$gene_id, age$gene_id)

stopifnot(!anyNA(idx_age))

age_aligned <- age[idx_age, , drop = FALSE]

stopifnot(
  identical(de$gene_id, age_aligned$gene_id),
  all(!is.na(age_aligned$logFC_reduced)),
  all(!is.na(age_aligned$logFC_age_adjusted)),
  all(!is.na(age_aligned$FDR_reduced)),
  all(!is.na(age_aligned$FDR_age_adjusted))
)

age_robust <- (
  age_aligned$FDR_reduced < 0.05 &
  age_aligned$FDR_age_adjusted < 0.05 &
  age_aligned$direction_concordant
)

cat("\n=== AGE-SENSITIVITY STRUCTURE ===\n")
cat(
  "Direction concordant:           ",
  sum(age_aligned$direction_concordant),
  "\n",
  sep = ""
)
cat(
  "FDR < 0.05 reduced:             ",
  sum(age_aligned$FDR_reduced < 0.05),
  "\n",
  sep = ""
)
cat(
  "FDR < 0.05 age-adjusted:        ",
  sum(age_aligned$FDR_age_adjusted < 0.05),
  "\n",
  sep = ""
)
cat(
  "FDR significant in both:        ",
  sum(
    age_aligned$FDR_reduced < 0.05 &
    age_aligned$FDR_age_adjusted < 0.05
  ),
  "\n",
  sep = ""
)
cat(
  "Strong age-robust set:          ",
  sum(age_robust),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Enrichment-ID eligibility audit
#
# Frozen rule:
#   ENTREZID must be present, uniquely mapped for the Ensembl
#   gene, and must occur only once in the eligible universe.
# ------------------------------------------------------------

entrez_present <- (
  !is.na(de$ENTREZID) &
  nzchar(trimws(as.character(de$ENTREZID)))
)

single_entrez <- (
  !is.na(de$n_entrezid) &
  de$n_entrezid == 1
)

candidate_entrez <- entrez_present & single_entrez

candidate_ids <- as.character(de$ENTREZID[candidate_entrez])

entrez_frequency <- table(candidate_ids)

unique_candidate_id <- candidate_ids %in%
  names(entrez_frequency)[entrez_frequency == 1L]

enrichment_eligible <- rep(FALSE, nrow(de))
enrichment_eligible[candidate_entrez] <- unique_candidate_id

n_missing_entrez <- sum(!entrez_present)
n_ambiguous_entrez <- sum(entrez_present & !single_entrez)
n_shared_entrez_rows <- sum(candidate_entrez) -
  sum(enrichment_eligible)

stopifnot(
  !anyDuplicated(
    as.character(de$ENTREZID[enrichment_eligible])
  )
)

cat("\n=== ENRICHMENT IDENTIFIER AUDIT ===\n")
cat("Frozen DE universe:              ", nrow(de), "\n", sep = "")
cat("Missing ENTREZID:                ", n_missing_entrez, "\n", sep = "")
cat("Ambiguous ENTREZID mapping:      ", n_ambiguous_entrez, "\n", sep = "")
cat("Rows with shared ENTREZID:       ", n_shared_entrez_rows, "\n", sep = "")
cat(
  "One-to-one enrichment universe: ",
  sum(enrichment_eligible),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 8. ORA eligibility counts
# ------------------------------------------------------------

ora_up <- enrichment_eligible &
  de$adj.P.Val < 0.05 &
  de$logFC >= 1

ora_down <- enrichment_eligible &
  de$adj.P.Val < 0.05 &
  de$logFC <= -1

cat("\n=== FROZEN ORA INPUT COUNTS ===\n")
cat(
  "Eligible enrichment background: ",
  sum(enrichment_eligible),
  "\n",
  sep = ""
)
cat(
  "Higher-in-IPF foreground:        ",
  sum(ora_up),
  "\n",
  sep = ""
)
cat(
  "Lower-in-IPF foreground:         ",
  sum(ora_down),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 9. Ranked-analysis readiness
# ------------------------------------------------------------

rank_ready <- enrichment_eligible & !is.na(de$t)

stopifnot(
  !anyDuplicated(
    as.character(de$ENTREZID[rank_ready])
  )
)

cat("\n=== RANKED ENRICHMENT READINESS ===\n")
cat("Ranking statistic:                moderated t\n")
cat(
  "Genes eligible for ranking:       ",
  sum(rank_ready),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 10. Machine-readable audit
# ------------------------------------------------------------

audit <- data.frame(
  metric = c(
    "frozen_DE_universe",
    "age_sensitivity_universe",
    "primary_FDR_lt_0.05",
    "primary_abs_logFC_ge_1",
    "primary_FDR_lt_0.05_and_abs_logFC_ge_1",
    "age_direction_concordant",
    "age_FDR_significant_both",
    "age_strong_robust_set",
    "missing_ENTREZID",
    "ambiguous_ENTREZID_mapping",
    "rows_with_shared_ENTREZID",
    "enrichment_background",
    "ORA_higher_IPF",
    "ORA_lower_IPF",
    "ranked_enrichment_genes"
  ),
  value = c(
    nrow(de),
    nrow(age),
    sum(primary_sig),
    sum(large_effect),
    sum(priority),
    sum(age_aligned$direction_concordant),
    sum(
      age_aligned$FDR_reduced < 0.05 &
      age_aligned$FDR_age_adjusted < 0.05
    ),
    sum(age_robust),
    n_missing_entrez,
    n_ambiguous_entrez,
    n_shared_entrez_rows,
    sum(enrichment_eligible),
    sum(ora_up),
    sum(ora_down),
    sum(rank_ready)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  audit,
  "metadata/qc_interpretation_preflight.csv",
  row.names = FALSE
)

cat("\n=== OUTPUT ===\n")
cat("metadata/qc_interpretation_preflight.csv\n")

cat("\nIMPORTANT:\n")
cat("No individual gene identities were printed.\n")
cat("No genes were removed from the frozen DE table.\n")
cat("No biological interpretation was performed.\n")
cat("No enrichment analysis was performed.\n")

cat("\n============================================================\n")
cat("W6.6 BIOLOGICAL INTERPRETATION PREFLIGHT: PASSED\n")
cat("============================================================\n")

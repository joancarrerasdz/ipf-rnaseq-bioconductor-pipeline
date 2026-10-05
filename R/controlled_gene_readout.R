# ============================================================
# GSE213001 W6.7 — First controlled biological gene readout
#
# Ranking rules were frozen before gene identities were
# inspected for biological prioritisation.
#
# Primary descriptive set:
#   BH FDR < 0.05 AND |logFC| >= 1
#
# Deterministic order within direction:
#   1. adj.P.Val ascending
#   2. |logFC| descending
#   3. |moderated t| descending
#   4. gene_id ascending
#
# Age robustness and annotation status are descriptive only.
# They do not modify the primary ranking.
# ============================================================

source("renv/activate.R")

# ------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------

de <- read.csv(
  "results/differential_expression/primary_IPF_vs_NDC_annotated.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

age <- read.csv(
  "results/differential_expression/sensitivity_age_gene_concordance.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)

cat("\n============================================================\n")
cat("GSE213001 W6.7 — FIRST CONTROLLED BIOLOGICAL READOUT\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 2. Frozen-universe checks
# ------------------------------------------------------------

stopifnot(
  nrow(de) == 15012,
  nrow(age) == 15012,
  !anyDuplicated(de$gene_id),
  !anyDuplicated(age$gene_id)
)

idx_age <- match(de$gene_id, age$gene_id)

stopifnot(!anyNA(idx_age))

age <- age[idx_age, , drop = FALSE]

stopifnot(
  identical(
    as.character(de$gene_id),
    as.character(age$gene_id)
  )
)

# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

nonempty <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}

logical_flag <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  tolower(trimws(as.character(x))) %in%
    c("true", "t", "1")
}

# ------------------------------------------------------------
# 4. Age robustness
# ------------------------------------------------------------

age_direction_concordant <-
  logical_flag(age$direction_concordant)

strong_age_robust <-
  age$FDR_reduced < 0.05 &
  age$FDR_age_adjusted < 0.05 &
  age_direction_concordant

# ------------------------------------------------------------
# 5. Annotation status
# ------------------------------------------------------------

symbol_present <- nonempty(de$SYMBOL)

annotation_clean <-
  symbol_present &
  !is.na(de$n_symbol) &
  de$n_symbol == 1 &
  !is.na(de$n_mapping_records) &
  de$n_mapping_records == 1

# Reproduce frozen one-to-one ENTREZ enrichment eligibility.

entrez_present <- nonempty(de$ENTREZID)

single_entrez <-
  entrez_present &
  !is.na(de$n_entrezid) &
  de$n_entrezid == 1

candidate_ids <-
  as.character(de$ENTREZID[single_entrez])

entrez_frequency <- table(candidate_ids)

unique_entrez <- candidate_ids %in%
  names(entrez_frequency)[entrez_frequency == 1L]

enrichment_eligible <- rep(FALSE, nrow(de))
enrichment_eligible[single_entrez] <- unique_entrez

stopifnot(
  sum(enrichment_eligible) == 14827,
  !anyDuplicated(
    as.character(de$ENTREZID[enrichment_eligible])
  )
)

# ------------------------------------------------------------
# 6. Construct complete reporting table
# ------------------------------------------------------------

report <- data.frame(
  gene_id = de$gene_id,
  SYMBOL = de$SYMBOL,
  GENENAME = de$GENENAME,
  ENTREZID = de$ENTREZID,

  logFC = de$logFC,
  AveExpr = de$AveExpr,
  t = de$t,
  P.Value = de$P.Value,
  adj.P.Val = de$adj.P.Val,
  B = de$B,

  mapping_status = de$mapping_status,
  n_mapping_records = de$n_mapping_records,
  n_symbol = de$n_symbol,
  n_entrezid = de$n_entrezid,

  annotation_clean = annotation_clean,
  enrichment_eligible = enrichment_eligible,

  age_logFC_reduced = age$logFC_reduced,
  age_logFC_adjusted = age$logFC_age_adjusted,
  age_delta_logFC = age$delta_logFC_age_minus_reduced,
  age_abs_delta_logFC = age$abs_delta_logFC,
  age_direction_concordant = age_direction_concordant,
  age_FDR_reduced = age$FDR_reduced,
  age_FDR_adjusted = age$FDR_age_adjusted,
  strong_age_robust = strong_age_robust,

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 7. Frozen primary descriptive set
# ------------------------------------------------------------

priority <-
  report$adj.P.Val < 0.05 &
  abs(report$logFC) >= 1

priority_report <- report[priority, , drop = FALSE]

stopifnot(nrow(priority_report) == 1861)

# ------------------------------------------------------------
# 8. Deterministic ranking function
# ------------------------------------------------------------

rank_direction <- function(x, direction_name) {

  if (direction_name == "higher") {
    x <- x[x$logFC >= 1, , drop = FALSE]
    direction_label <- "Higher in IPF"
  } else {
    x <- x[x$logFC <= -1, , drop = FALSE]
    direction_label <- "Lower in IPF"
  }

  ord <- order(
    x$adj.P.Val,
    -abs(x$logFC),
    -abs(x$t),
    x$gene_id,
    na.last = TRUE
  )

  x <- x[ord, , drop = FALSE]

  x$rank_direction <- seq_len(nrow(x))
  x$direction <- direction_label

  x <- x[
    ,
    c(
      "rank_direction",
      "direction",
      "gene_id",
      "SYMBOL",
      "GENENAME",
      "ENTREZID",
      "logFC",
      "AveExpr",
      "t",
      "P.Value",
      "adj.P.Val",
      "B",
      "mapping_status",
      "n_mapping_records",
      "n_symbol",
      "n_entrezid",
      "annotation_clean",
      "enrichment_eligible",
      "age_logFC_reduced",
      "age_logFC_adjusted",
      "age_delta_logFC",
      "age_abs_delta_logFC",
      "age_direction_concordant",
      "age_FDR_reduced",
      "age_FDR_adjusted",
      "strong_age_robust"
    )
  ]

  x
}

higher <- rank_direction(priority_report, "higher")
lower  <- rank_direction(priority_report, "lower")

stopifnot(
  nrow(higher) == 1387,
  nrow(lower) == 474
)

# ------------------------------------------------------------
# 9. Complete frozen reporting outputs
# ------------------------------------------------------------

all_priority <- rbind(higher, lower)

write.csv(
  all_priority,
  "results/interpretation/W6_priority_genes_all.csv",
  row.names = FALSE,
  na = ""
)

write.csv(
  higher,
  "results/interpretation/W6_priority_higher_in_IPF.csv",
  row.names = FALSE,
  na = ""
)

write.csv(
  lower,
  "results/interpretation/W6_priority_lower_in_IPF.csv",
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 10. First controlled readout
# ------------------------------------------------------------

top_n <- 20L

first_readout <- rbind(
  head(higher, top_n),
  head(lower, top_n)
)

write.csv(
  first_readout,
  "results/interpretation/W6_first_readout_top20_each_direction.csv",
  row.names = FALSE,
  na = ""
)

# ------------------------------------------------------------
# 11. Quantitative summary
# ------------------------------------------------------------

summary_df <- data.frame(
  metric = c(
    "priority_genes_total",
    "priority_higher_IPF",
    "priority_lower_IPF",
    "priority_annotation_clean",
    "priority_enrichment_eligible",
    "priority_strong_age_robust",
    "higher_strong_age_robust",
    "lower_strong_age_robust",
    "higher_annotation_clean",
    "lower_annotation_clean"
  ),
  value = c(
    nrow(all_priority),
    nrow(higher),
    nrow(lower),
    sum(all_priority$annotation_clean),
    sum(all_priority$enrichment_eligible),
    sum(all_priority$strong_age_robust),
    sum(higher$strong_age_robust),
    sum(lower$strong_age_robust),
    sum(higher$annotation_clean),
    sum(lower$annotation_clean)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  summary_df,
  "metadata/qc_W6_controlled_gene_readout.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Console certificate
# ------------------------------------------------------------

cat("\n=== FROZEN PRIMARY DESCRIPTIVE SET ===\n")
cat("Total:                         ", nrow(all_priority), "\n", sep = "")
cat("Higher in IPF:                 ", nrow(higher), "\n", sep = "")
cat("Lower in IPF:                  ", nrow(lower), "\n", sep = "")
cat(
  "Annotation-clean:              ",
  sum(all_priority$annotation_clean),
  "\n",
  sep = ""
)
cat(
  "Enrichment-eligible:           ",
  sum(all_priority$enrichment_eligible),
  "\n",
  sep = ""
)
cat(
  "Strongly age-robust:           ",
  sum(all_priority$strong_age_robust),
  "\n",
  sep = ""
)

console_cols <- c(
  "rank_direction",
  "SYMBOL",
  "gene_id",
  "logFC",
  "adj.P.Val",
  "strong_age_robust",
  "mapping_status"
)

cat("\n============================================================\n")
cat("FIRST GENE-IDENTITY READOUT — HIGHER IN IPF\n")
cat("============================================================\n")

print(
  head(higher[, console_cols, drop = FALSE], top_n),
  row.names = FALSE
)

cat("\n============================================================\n")
cat("FIRST GENE-IDENTITY READOUT — LOWER IN IPF\n")
cat("============================================================\n")

print(
  head(lower[, console_cols, drop = FALSE], top_n),
  row.names = FALSE
)

cat("\n=== OUTPUTS ===\n")
cat("results/interpretation/W6_priority_genes_all.csv\n")
cat("results/interpretation/W6_priority_higher_in_IPF.csv\n")
cat("results/interpretation/W6_priority_lower_in_IPF.csv\n")
cat("results/interpretation/W6_first_readout_top20_each_direction.csv\n")
cat("metadata/qc_W6_controlled_gene_readout.csv\n")

cat("\nIMPORTANT:\n")
cat("The W5 statistical model was not recalculated.\n")
cat("The frozen 15,012-gene universe was not modified.\n")
cat("Age robustness did not affect ranking.\n")
cat("Annotation content did not affect ranking.\n")
cat("No pathway enrichment was performed.\n")

cat("\n============================================================\n")
cat("W6.7 FIRST CONTROLLED BIOLOGICAL READOUT: COMPLETED\n")
cat("============================================================\n")

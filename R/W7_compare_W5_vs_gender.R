# ============================================================
# GSE213001
# W7.3e - Frozen W5 vs gender-adjusted W7 comparison
#
# PURPOSE
# Compare the frozen W5 primary IPF-vs-NDC differential-
# expression result against the prespecified W7 gender-adjusted
# sensitivity result.
#
# IMPORTANT
# - No differential-expression model is fitted here.
# - No voom transformation is recalculated.
# - No duplicateCorrelation estimate is recalculated.
# - No W5 result is modified.
# - No W6 result is modified.
# - No W7 gender model is modified.
# - This script performs deterministic post-hoc comparison only.
# ============================================================

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# 1. Frozen inputs
# ------------------------------------------------------------

w5_file <- paste0(
  "results/differential_expression/",
  "primary_IPF_vs_NDC_all_genes.csv"
)

w7_file <- paste0(
  "results/sensitivity/",
  "W7_gender_IPF_vs_NDC_all_genes.csv"
)

output_dir <- "results/sensitivity"

comparison_file <- file.path(
  output_dir,
  "W7_gender_vs_W5_gene_comparison.csv"
)

summary_file <- file.path(
  output_dir,
  "W7_gender_vs_W5_summary.csv"
)

qc_file <- "metadata/qc_W7_gender_vs_W5_comparison.csv"

stopifnot(
  file.exists(w5_file),
  file.exists(w7_file)
)

# ------------------------------------------------------------
# 2. Read frozen results
# ------------------------------------------------------------

w5 <- read.csv(
  w5_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

w7 <- read.csv(
  w7_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "gene_id",
  "logFC",
  "P.Value",
  "adj.P.Val"
)

# ------------------------------------------------------------
# 3. Frozen-input QC
# ------------------------------------------------------------

stopifnot(
  all(required_columns %in% names(w5)),
  all(required_columns %in% names(w7)),

  nrow(w5) == 15012L,
  nrow(w7) == 15012L,

  length(unique(w5$gene_id)) == 15012L,
  length(unique(w7$gene_id)) == 15012L,

  !anyDuplicated(w5$gene_id),
  !anyDuplicated(w7$gene_id),

  setequal(w5$gene_id, w7$gene_id),
  identical(w5$gene_id, w7$gene_id),

  all(is.finite(w5$logFC)),
  all(is.finite(w7$logFC)),

  all(is.finite(w5$P.Value)),
  all(is.finite(w7$P.Value)),

  all(is.finite(w5$adj.P.Val)),
  all(is.finite(w7$adj.P.Val)),

  all(w5$P.Value >= 0 & w5$P.Value <= 1),
  all(w7$P.Value >= 0 & w7$P.Value <= 1),

  all(w5$adj.P.Val >= 0 & w5$adj.P.Val <= 1),
  all(w7$adj.P.Val >= 0 & w7$adj.P.Val <= 1),

  sum(w5$adj.P.Val < 0.05) == 7898L,
  sum(w7$adj.P.Val < 0.05) == 7745L
)

# ------------------------------------------------------------
# 4. Gene-level comparison
# ------------------------------------------------------------

w5_fdr <- w5$adj.P.Val < 0.05
w7_fdr <- w7$adj.P.Val < 0.05

delta_logFC <- w7$logFC - w5$logFC
abs_delta_logFC <- abs(delta_logFC)

direction_concordant <- sign(w5$logFC) == sign(w7$logFC)

FDR_status <- ifelse(
  w5_fdr & w7_fdr,
  "shared_FDR",
  ifelse(
    w5_fdr & !w7_fdr,
    "W5_only_FDR",
    ifelse(
      !w5_fdr & w7_fdr,
      "W7_only_FDR",
      "neither_FDR"
    )
  )
)

comparison <- data.frame(
  gene_id = w5$gene_id,

  W5_logFC = w5$logFC,
  W7_logFC = w7$logFC,
  delta_logFC = delta_logFC,
  abs_delta_logFC = abs_delta_logFC,

  W5_P_Value = w5$P.Value,
  W7_P_Value = w7$P.Value,

  W5_adj_P_Val = w5$adj.P.Val,
  W7_adj_P_Val = w7$adj.P.Val,

  W5_FDR = w5_fdr,
  W7_FDR = w7_fdr,

  FDR_status = FDR_status,
  direction_concordant = direction_concordant,

  check.names = FALSE
)

# ------------------------------------------------------------
# 5. FDR-set metrics
# ------------------------------------------------------------

shared_FDR <- sum(FDR_status == "shared_FDR")
W5_only_FDR <- sum(FDR_status == "W5_only_FDR")
W7_only_FDR <- sum(FDR_status == "W7_only_FDR")
neither_FDR <- sum(FDR_status == "neither_FDR")

FDR_union <- shared_FDR + W5_only_FDR + W7_only_FDR

FDR_Jaccard <- shared_FDR / FDR_union

W5_FDR_retained_in_W7 <- shared_FDR /
  (shared_FDR + W5_only_FDR)

W7_FDR_already_present_in_W5 <- shared_FDR /
  (shared_FDR + W7_only_FDR)

# ------------------------------------------------------------
# 6. Genome-wide effect stability
# ------------------------------------------------------------

Pearson_logFC <- cor(
  comparison$W5_logFC,
  comparison$W7_logFC,
  method = "pearson"
)

Spearman_logFC <- cor(
  comparison$W5_logFC,
  comparison$W7_logFC,
  method = "spearman"
)

direction_concordance_all <- mean(
  comparison$direction_concordant
)

median_abs_delta_logFC <- median(
  comparison$abs_delta_logFC
)

mean_abs_delta_logFC <- mean(
  comparison$abs_delta_logFC
)

q95_abs_delta_logFC <- as.numeric(
  quantile(
    comparison$abs_delta_logFC,
    probs = 0.95,
    names = FALSE,
    type = 7
  )
)

q99_abs_delta_logFC <- as.numeric(
  quantile(
    comparison$abs_delta_logFC,
    probs = 0.99,
    names = FALSE,
    type = 7
  )
)

max_abs_delta_logFC <- max(
  comparison$abs_delta_logFC
)

# ------------------------------------------------------------
# 7. Shared-FDR effect stability
# ------------------------------------------------------------

shared_idx <- comparison$FDR_status == "shared_FDR"

stopifnot(sum(shared_idx) == shared_FDR)

shared_FDR_direction_concordance <- mean(
  comparison$direction_concordant[shared_idx]
)

shared_FDR_direction_discordant <- sum(
  !comparison$direction_concordant[shared_idx]
)

shared_FDR_Pearson_logFC <- cor(
  comparison$W5_logFC[shared_idx],
  comparison$W7_logFC[shared_idx],
  method = "pearson"
)

shared_FDR_Spearman_logFC <- cor(
  comparison$W5_logFC[shared_idx],
  comparison$W7_logFC[shared_idx],
  method = "spearman"
)

shared_FDR_median_abs_delta_logFC <- median(
  comparison$abs_delta_logFC[shared_idx]
)

shared_FDR_q95_abs_delta_logFC <- as.numeric(
  quantile(
    comparison$abs_delta_logFC[shared_idx],
    probs = 0.95,
    names = FALSE,
    type = 7
  )
)

# ------------------------------------------------------------
# 8. Deterministic summary
# ------------------------------------------------------------

summary <- data.frame(
  metric = c(
    "genes",
    "W5_FDR_0_05",
    "W7_FDR_0_05",
    "shared_FDR",
    "W5_only_FDR",
    "W7_only_FDR",
    "neither_FDR",
    "FDR_union",
    "FDR_Jaccard",
    "W5_FDR_retained_in_W7",
    "W7_FDR_already_present_in_W5",
    "Pearson_logFC",
    "Spearman_logFC",
    "direction_concordance_all",
    "shared_FDR_direction_concordance",
    "shared_FDR_direction_discordant",
    "median_abs_delta_logFC",
    "mean_abs_delta_logFC",
    "q95_abs_delta_logFC",
    "q99_abs_delta_logFC",
    "max_abs_delta_logFC",
    "shared_FDR_Pearson_logFC",
    "shared_FDR_Spearman_logFC",
    "shared_FDR_median_abs_delta_logFC",
    "shared_FDR_q95_abs_delta_logFC"
  ),

  value = c(
    nrow(comparison),
    sum(w5_fdr),
    sum(w7_fdr),
    shared_FDR,
    W5_only_FDR,
    W7_only_FDR,
    neither_FDR,
    FDR_union,
    FDR_Jaccard,
    W5_FDR_retained_in_W7,
    W7_FDR_already_present_in_W5,
    Pearson_logFC,
    Spearman_logFC,
    direction_concordance_all,
    shared_FDR_direction_concordance,
    shared_FDR_direction_discordant,
    median_abs_delta_logFC,
    mean_abs_delta_logFC,
    q95_abs_delta_logFC,
    q99_abs_delta_logFC,
    max_abs_delta_logFC,
    shared_FDR_Pearson_logFC,
    shared_FDR_Spearman_logFC,
    shared_FDR_median_abs_delta_logFC,
    shared_FDR_q95_abs_delta_logFC
  ),

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 9. QC metadata
# ------------------------------------------------------------

qc <- data.frame(
  metric = c(
    "W5_rows",
    "W7_rows",
    "W5_unique_gene_ids",
    "W7_unique_gene_ids",
    "same_gene_set",
    "same_gene_order",
    "W5_FDR_0_05",
    "W7_FDR_0_05",
    "comparison_rows",
    "summary_rows",
    "shared_FDR_direction_discordant"
  ),

  value = c(
    as.character(nrow(w5)),
    as.character(nrow(w7)),
    as.character(length(unique(w5$gene_id))),
    as.character(length(unique(w7$gene_id))),
    as.character(setequal(w5$gene_id, w7$gene_id)),
    as.character(identical(w5$gene_id, w7$gene_id)),
    as.character(sum(w5_fdr)),
    as.character(sum(w7_fdr)),
    as.character(nrow(comparison)),
    as.character(nrow(summary)),
    as.character(shared_FDR_direction_discordant)
  ),

  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 10. Output-level invariants
# ------------------------------------------------------------

stopifnot(
  nrow(comparison) == 15012L,
  length(unique(comparison$gene_id)) == 15012L,
  !anyDuplicated(comparison$gene_id),

  identical(comparison$gene_id, w5$gene_id),
  identical(comparison$gene_id, w7$gene_id),

  sum(comparison$FDR_status == "shared_FDR") +
    sum(comparison$FDR_status == "W5_only_FDR") ==
      7898L,

  sum(comparison$FDR_status == "shared_FDR") +
    sum(comparison$FDR_status == "W7_only_FDR") ==
      7745L,

  sum(
    comparison$FDR_status %in% c(
      "shared_FDR",
      "W5_only_FDR",
      "W7_only_FDR",
      "neither_FDR"
    )
  ) == 15012L,

  all(is.finite(comparison$delta_logFC)),
  all(is.finite(comparison$abs_delta_logFC)),

  Pearson_logFC >= -1,
  Pearson_logFC <= 1,

  Spearman_logFC >= -1,
  Spearman_logFC <= 1,

  direction_concordance_all >= 0,
  direction_concordance_all <= 1,

  shared_FDR_direction_concordance >= 0,
  shared_FDR_direction_concordance <= 1,

  shared_FDR_direction_discordant >= 0,

  FDR_Jaccard >= 0,
  FDR_Jaccard <= 1,

  W5_FDR_retained_in_W7 >= 0,
  W5_FDR_retained_in_W7 <= 1,

  W7_FDR_already_present_in_W5 >= 0,
  W7_FDR_already_present_in_W5 <= 1,

  nrow(summary) == 25L,
  nrow(qc) == 11L
)

# ------------------------------------------------------------
# 11. Write outputs
# ------------------------------------------------------------

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  comparison,
  comparison_file,
  row.names = FALSE
)

write.csv(
  summary,
  summary_file,
  row.names = FALSE
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Console certificate
# ------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "GSE213001 W7.3e - W5 vs GENDER-ADJUSTED W7 COMPARISON\n"
)

cat(
  "============================================================\n"
)

cat("\n===== INPUT =====\n")
cat("Genes:", nrow(comparison), "\n")
cat("W5 FDR < 0.05:", sum(w5_fdr), "\n")
cat("W7 FDR < 0.05:", sum(w7_fdr), "\n")

cat("\n===== FDR OVERLAP =====\n")
cat("Shared FDR:", shared_FDR, "\n")
cat("W5 only:", W5_only_FDR, "\n")
cat("W7 only:", W7_only_FDR, "\n")
cat("Neither:", neither_FDR, "\n")
cat("FDR union:", FDR_union, "\n")
cat("FDR Jaccard:", signif(FDR_Jaccard, 7), "\n")
cat(
  "W5 FDR retained in W7:",
  signif(W5_FDR_retained_in_W7, 7),
  "\n"
)
cat(
  "W7 FDR already present in W5:",
  signif(W7_FDR_already_present_in_W5, 7),
  "\n"
)

cat("\n===== GENOME-WIDE EFFECT STABILITY =====\n")
cat("Pearson logFC:", signif(Pearson_logFC, 7), "\n")
cat("Spearman logFC:", signif(Spearman_logFC, 7), "\n")
cat(
  "Direction concordance:",
  signif(direction_concordance_all, 7),
  "\n"
)
cat(
  "Median |delta logFC|:",
  signif(median_abs_delta_logFC, 7),
  "\n"
)
cat(
  "Mean |delta logFC|:",
  signif(mean_abs_delta_logFC, 7),
  "\n"
)
cat(
  "95th percentile |delta logFC|:",
  signif(q95_abs_delta_logFC, 7),
  "\n"
)
cat(
  "99th percentile |delta logFC|:",
  signif(q99_abs_delta_logFC, 7),
  "\n"
)
cat(
  "Max |delta logFC|:",
  signif(max_abs_delta_logFC, 7),
  "\n"
)

cat("\n===== SHARED-FDR EFFECT STABILITY =====\n")
cat("Shared FDR genes:", shared_FDR, "\n")
cat(
  "Pearson logFC:",
  signif(shared_FDR_Pearson_logFC, 7),
  "\n"
)
cat(
  "Spearman logFC:",
  signif(shared_FDR_Spearman_logFC, 7),
  "\n"
)
cat(
  "Direction concordance:",
  signif(shared_FDR_direction_concordance, 7),
  "\n"
)
cat(
  "Direction discordant:",
  shared_FDR_direction_discordant,
  "\n"
)
cat(
  "Median |delta logFC|:",
  signif(shared_FDR_median_abs_delta_logFC, 7),
  "\n"
)
cat(
  "95th percentile |delta logFC|:",
  signif(shared_FDR_q95_abs_delta_logFC, 7),
  "\n"
)

cat("\n===== OUTPUTS =====\n")
cat(comparison_file, "\n")
cat(summary_file, "\n")
cat(qc_file, "\n")

cat("\nIMPORTANT:\n")
cat("No differential-expression model was fitted.\n")
cat("No W5 result was recalculated or modified.\n")
cat("No W6 result was recalculated or modified.\n")
cat("No W7 gender model was recalculated or modified.\n")
cat("This is a deterministic comparison of frozen results.\n")

cat(
  "\n============================================================\n"
)
cat(
  "W7.3e W5 vs GENDER COMPARISON: COMPLETED\n"
)
cat(
  "============================================================\n"
)

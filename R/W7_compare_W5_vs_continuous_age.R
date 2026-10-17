# ============================================================
# GSE213001
# W7.2d — W5 vs continuous-age W7 sensitivity comparison
#
# Purpose:
#   Compare the frozen W5 primary IPF-vs-NDC differential-
#   expression result with the frozen W7 continuous-age
#   sensitivity model.
#
# Important:
#   - W5 is the frozen primary reference.
#   - W7 is a sensitivity analysis only.
#   - No differential-expression model is fitted here.
#   - No voom transformation is recalculated.
#   - No duplicateCorrelation estimate is recalculated.
#   - No P value or adjusted P value is recalculated.
#   - W5 and W7 inputs are never modified.
#   - All comparison metrics are descriptive robustness metrics.
# ============================================================


# ------------------------------------------------------------
# 1. Frozen inputs
# ------------------------------------------------------------

w5_file <- file.path(
  "results",
  "differential_expression",
  "primary_IPF_vs_NDC_all_genes.csv"
)

w7_file <- file.path(
  "results",
  "sensitivity",
  "W7_continuous_age_IPF_vs_NDC_all_genes.csv"
)

output_dir <- file.path(
  "results",
  "sensitivity"
)

comparison_file <- file.path(
  output_dir,
  "W7_continuous_age_vs_W5_gene_comparison.csv"
)

summary_file <- file.path(
  output_dir,
  "W7_continuous_age_vs_W5_summary.csv"
)

qc_file <- file.path(
  "metadata",
  "qc_W7_continuous_age_vs_W5_comparison.csv"
)


# ------------------------------------------------------------
# 2. Input existence
# ------------------------------------------------------------

stopifnot(
  file.exists(w5_file),
  file.exists(w7_file)
)


# ------------------------------------------------------------
# 3. Read frozen results
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

stopifnot(
  all(required_columns %in% names(w5)),
  all(required_columns %in% names(w7)),
  nrow(w5) == 15012L,
  nrow(w7) == 15012L,
  !anyDuplicated(w5$gene_id),
  !anyDuplicated(w7$gene_id),
  identical(w5$gene_id, w7$gene_id),
  all(is.finite(w5$logFC)),
  all(is.finite(w7$logFC)),
  all(w5$P.Value >= 0 & w5$P.Value <= 1),
  all(w7$P.Value >= 0 & w7$P.Value <= 1),
  all(w5$adj.P.Val >= 0 & w5$adj.P.Val <= 1),
  all(w7$adj.P.Val >= 0 & w7$adj.P.Val <= 1)
)


# ------------------------------------------------------------
# 4. Frozen preflight expectations
# ------------------------------------------------------------

w5_fdr <- w5$adj.P.Val < 0.05
w7_fdr <- w7$adj.P.Val < 0.05

stopifnot(
  sum(w5_fdr) == 7898L,
  sum(w7_fdr) == 5277L,

  sum(w5_fdr & w5$logFC > 0) == 3964L,
  sum(w5_fdr & w5$logFC < 0) == 3934L,

  sum(w7_fdr & w7$logFC > 0) == 2610L,
  sum(w7_fdr & w7$logFC < 0) == 2667L
)


# ------------------------------------------------------------
# 5. Gene-level comparison table
# ------------------------------------------------------------

comparison <- data.frame(
  gene_id = w5$gene_id,

  W5_logFC = w5$logFC,
  W7_logFC = w7$logFC,

  delta_logFC = w7$logFC - w5$logFC,
  abs_delta_logFC = abs(w7$logFC - w5$logFC),

  W5_P_Value = w5$P.Value,
  W7_P_Value = w7$P.Value,

  W5_adj_P_Val = w5$adj.P.Val,
  W7_adj_P_Val = w7$adj.P.Val,

  W5_FDR_0_05 = w5_fdr,
  W7_FDR_0_05 = w7_fdr,

  stringsAsFactors = FALSE
)

comparison$W5_direction <- ifelse(
  comparison$W5_logFC > 0,
  "higher_in_IPF",
  ifelse(
    comparison$W5_logFC < 0,
    "lower_in_IPF",
    "zero"
  )
)

comparison$W7_direction <- ifelse(
  comparison$W7_logFC > 0,
  "higher_in_IPF",
  ifelse(
    comparison$W7_logFC < 0,
    "lower_in_IPF",
    "zero"
  )
)

comparison$direction_concordant <- (
  comparison$W5_direction ==
    comparison$W7_direction
)

comparison$FDR_status <- ifelse(
  comparison$W5_FDR_0_05 &
    comparison$W7_FDR_0_05,
  "shared_FDR",
  ifelse(
    comparison$W5_FDR_0_05 &
      !comparison$W7_FDR_0_05,
    "W5_only_FDR",
    ifelse(
      !comparison$W5_FDR_0_05 &
        comparison$W7_FDR_0_05,
      "W7_only_FDR",
      "neither_FDR"
    )
  )
)

comparison$shared_FDR_direction_concordant <- ifelse(
  comparison$FDR_status == "shared_FDR",
  comparison$direction_concordant,
  NA
)


# ------------------------------------------------------------
# 6. Global robustness metrics
# ------------------------------------------------------------

pearson_logFC <- cor(
  comparison$W5_logFC,
  comparison$W7_logFC,
  method = "pearson"
)

spearman_logFC <- cor(
  comparison$W5_logFC,
  comparison$W7_logFC,
  method = "spearman"
)

direction_concordance_all <- mean(
  comparison$direction_concordant
)

shared_fdr <- sum(
  comparison$FDR_status == "shared_FDR"
)

w5_only_fdr <- sum(
  comparison$FDR_status == "W5_only_FDR"
)

w7_only_fdr <- sum(
  comparison$FDR_status == "W7_only_FDR"
)

neither_fdr <- sum(
  comparison$FDR_status == "neither_FDR"
)

union_fdr <- shared_fdr +
  w5_only_fdr +
  w7_only_fdr

fdr_jaccard <- shared_fdr / union_fdr

shared_fdr_direction_concordance <- if (
  shared_fdr > 0L
) {
  mean(
    comparison$direction_concordant[
      comparison$FDR_status ==
        "shared_FDR"
    ]
  )
} else {
  NA_real_
}

median_abs_delta <- median(
  comparison$abs_delta_logFC
)

mean_abs_delta <- mean(
  comparison$abs_delta_logFC
)

max_abs_delta <- max(
  comparison$abs_delta_logFC
)

q95_abs_delta <- unname(
  quantile(
    comparison$abs_delta_logFC,
    probs = 0.95
  )
)

q99_abs_delta <- unname(
  quantile(
    comparison$abs_delta_logFC,
    probs = 0.99
  )
)


# ------------------------------------------------------------
# 7. Summary table
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
    "Pearson_logFC",
    "Spearman_logFC",
    "direction_concordance_all",
    "shared_FDR_direction_concordance",
    "median_abs_delta_logFC",
    "mean_abs_delta_logFC",
    "q95_abs_delta_logFC",
    "q99_abs_delta_logFC",
    "max_abs_delta_logFC"
  ),

  value = c(
    nrow(comparison),
    sum(w5_fdr),
    sum(w7_fdr),
    shared_fdr,
    w5_only_fdr,
    w7_only_fdr,
    neither_fdr,
    union_fdr,
    fdr_jaccard,
    pearson_logFC,
    spearman_logFC,
    direction_concordance_all,
    shared_fdr_direction_concordance,
    median_abs_delta,
    mean_abs_delta,
    q95_abs_delta,
    q99_abs_delta,
    max_abs_delta
  ),

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 8. QC metadata
# ------------------------------------------------------------

qc <- data.frame(
  metric = c(
    "W5_rows",
    "W7_rows",
    "gene_order_identical",
    "W5_FDR_0_05",
    "W7_FDR_0_05",
    "shared_FDR",
    "W5_only_FDR",
    "W7_only_FDR",
    "neither_FDR",
    "direction_concordant_genes",
    "direction_discordant_genes"
  ),

  value = c(
    nrow(w5),
    nrow(w7),
    identical(w5$gene_id, w7$gene_id),
    sum(w5_fdr),
    sum(w7_fdr),
    shared_fdr,
    w5_only_fdr,
    w7_only_fdr,
    neither_fdr,
    sum(comparison$direction_concordant),
    sum(!comparison$direction_concordant)
  ),

  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 9. Structural assertions
# ------------------------------------------------------------

stopifnot(
  nrow(comparison) == 15012L,
  !anyDuplicated(comparison$gene_id),

  shared_fdr +
    w5_only_fdr ==
    7898L,

  shared_fdr +
    w7_only_fdr ==
    5277L,

  shared_fdr +
    w5_only_fdr +
    w7_only_fdr +
    neither_fdr ==
    15012L,

  is.finite(pearson_logFC),
  is.finite(spearman_logFC),

  pearson_logFC >= -1 &
    pearson_logFC <= 1,

  spearman_logFC >= -1 &
    spearman_logFC <= 1,

  direction_concordance_all >= 0 &
    direction_concordance_all <= 1,

  fdr_jaccard >= 0 &
    fdr_jaccard <= 1,

  nrow(summary) == 18L,
  nrow(qc) == 11L
)


# ------------------------------------------------------------
# 10. Write outputs
# ------------------------------------------------------------

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  comparison,
  comparison_file,
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  summary,
  summary_file,
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE,
  quote = TRUE
)


# ------------------------------------------------------------
# 11. Console certificate
# ------------------------------------------------------------

cat("============================================================\n")
cat("GSE213001 W7.2d — W5 vs CONTINUOUS-AGE W7 COMPARISON\n")
cat("============================================================\n")

cat("\n===== INPUT =====\n")
cat("Genes: ", nrow(comparison), "\n", sep = "")
cat("W5 FDR < 0.05: ", sum(w5_fdr), "\n", sep = "")
cat("W7 FDR < 0.05: ", sum(w7_fdr), "\n", sep = "")

cat("\n===== EFFECT-SIZE STABILITY =====\n")
cat(
  "Pearson logFC: ",
  signif(pearson_logFC, 7),
  "\n",
  sep = ""
)
cat(
  "Spearman logFC: ",
  signif(spearman_logFC, 7),
  "\n",
  sep = ""
)
cat(
  "Direction concordance: ",
  signif(direction_concordance_all, 7),
  "\n",
  sep = ""
)
cat(
  "Median |delta logFC|: ",
  signif(median_abs_delta, 7),
  "\n",
  sep = ""
)
cat(
  "95th percentile |delta logFC|: ",
  signif(q95_abs_delta, 7),
  "\n",
  sep = ""
)

cat("\n===== FDR OVERLAP =====\n")
cat("Shared FDR: ", shared_fdr, "\n", sep = "")
cat("W5-only FDR: ", w5_only_fdr, "\n", sep = "")
cat("W7-only FDR: ", w7_only_fdr, "\n", sep = "")
cat("Neither FDR: ", neither_fdr, "\n", sep = "")
cat(
  "FDR Jaccard: ",
  signif(fdr_jaccard, 7),
  "\n",
  sep = ""
)
cat(
  "Shared-FDR direction concordance: ",
  signif(
    shared_fdr_direction_concordance,
    7
  ),
  "\n",
  sep = ""
)

cat("\n===== OUTPUTS =====\n")
cat(comparison_file, "\n")
cat(summary_file, "\n")
cat(qc_file, "\n")

cat("\nIMPORTANT:\n")
cat("W5 remains the frozen primary reference.\n")
cat("W7 remains a sensitivity analysis only.\n")
cat("No differential-expression model was fitted.\n")
cat("No statistical result was recalculated.\n")
cat("No W5 or W7 input was modified.\n")
cat("No biological conclusion is generated by this script.\n")

cat("\n============================================================\n")
cat("W7.2d W5 vs CONTINUOUS-AGE COMPARISON: COMPLETED\n")
cat("============================================================\n")

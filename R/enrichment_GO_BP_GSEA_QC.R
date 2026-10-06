# ============================================================
# GSE213001 - W6.12d GO:BP ranked enrichment result QC
#
# Independent structural and numerical audit of the frozen
# W6.12c ranked-enrichment output.
#
# IMPORTANT:
#   - no pathway descriptions are printed
#   - no biological interpretation is performed
#   - W5 statistics are not recalculated
# ============================================================

rank_file <-
  "results/enrichment/input/W6_ranked_enrichment_moderated_t.csv"

all_file <-
  "results/enrichment/GO_BP/W6_GO_BP_GSEA_all_terms.csv"

sig_file <-
  "results/enrichment/GO_BP/W6_GO_BP_GSEA_significant.csv"

summary_file <-
  "metadata/qc_W6_GO_BP_GSEA.csv"

qc_file <-
  "metadata/qc_W6_GO_BP_GSEA_result_qc.csv"

stopifnot(
  file.exists(rank_file),
  file.exists(all_file),
  file.exists(sig_file),
  file.exists(summary_file)
)

ranked <- read.csv(
  rank_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

all <- read.csv(
  all_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sig <- read.csv(
  sig_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

summary_df <- read.csv(
  summary_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

summary_values <- setNames(
  summary_df$value,
  summary_df$metric
)

num_summary <- function(x) {
  as.numeric(summary_values[[x]])
}

cat("\n============================================================\n")
cat("GSE213001 W6.12d - GO:BP RANKED ENRICHMENT RESULT QC\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. Frozen ranked-input audit
# ------------------------------------------------------------

rank_required <- c(
  "gene_id",
  "ENTREZID",
  "t"
)

rank_required_present <-
  all(rank_required %in% names(ranked))

stopifnot(rank_required_present)

rank_rows <- nrow(ranked)
rank_unique_gene_id <- length(unique(ranked$gene_id))
rank_unique_entrez <- length(unique(ranked$ENTREZID))
rank_missing_entrez <- sum(
  is.na(ranked$ENTREZID) |
    ranked$ENTREZID == ""
)
rank_missing_t <- sum(is.na(ranked$t))
rank_nonfinite_t <- sum(!is.finite(ranked$t))
rank_positive_t <- sum(ranked$t > 0)
rank_negative_t <- sum(ranked$t < 0)
rank_zero_t <- sum(ranked$t == 0)
rank_decreasing <-
  identical(order(ranked$t, decreasing = TRUE), seq_len(nrow(ranked)))

cat("\n=== FROZEN RANKED INPUT ===\n")
cat("Rows:                  ", rank_rows, "\n", sep = "")
cat("Unique gene_id:        ", rank_unique_gene_id, "\n", sep = "")
cat("Unique ENTREZID:       ", rank_unique_entrez, "\n", sep = "")
cat("Positive t:            ", rank_positive_t, "\n", sep = "")
cat("Negative t:            ", rank_negative_t, "\n", sep = "")
cat("Zero t:                ", rank_zero_t, "\n", sep = "")
cat("Missing ENTREZID:      ", rank_missing_entrez, "\n", sep = "")
cat("Missing t:             ", rank_missing_t, "\n", sep = "")
cat("Non-finite t:          ", rank_nonfinite_t, "\n", sep = "")
cat("Already decreasing t: ", rank_decreasing, "\n", sep = "")

stopifnot(
  rank_rows == 14827L,
  rank_unique_gene_id == 14827L,
  rank_unique_entrez == 14827L,
  rank_missing_entrez == 0L,
  rank_missing_t == 0L,
  rank_nonfinite_t == 0L,
  rank_positive_t == 7334L,
  rank_negative_t == 7493L,
  rank_zero_t == 0L,
  rank_decreasing
)

# ------------------------------------------------------------
# 2. Result structure
# ------------------------------------------------------------

required_cols <- c(
  "ID",
  "setSize",
  "NES",
  "pvalue",
  "p.adjust",
  "rank",
  "core_enrichment"
)

required_all <-
  all(required_cols %in% names(all))

required_sig <-
  all(required_cols %in% names(sig))

cat("\n=== RESULT STRUCTURE ===\n")
cat("All terms:             ", nrow(all), "\n", sep = "")
cat("Significant terms:     ", nrow(sig), "\n", sep = "")
cat("Required columns/all:  ", required_all, "\n", sep = "")
cat("Required columns/sig:  ", required_sig, "\n", sep = "")

stopifnot(
  required_all,
  required_sig,
  nrow(all) == 5823L,
  nrow(sig) == 1054L
)

# ------------------------------------------------------------
# 3. ID integrity
# ------------------------------------------------------------

all_unique_ids <- length(unique(all$ID))
sig_unique_ids <- length(unique(sig$ID))
duplicated_all_ids <- sum(duplicated(all$ID))
duplicated_sig_ids <- sum(duplicated(sig$ID))

cat("\n=== ID INTEGRITY ===\n")
cat("Unique all IDs:        ", all_unique_ids, "\n", sep = "")
cat("Unique significant IDs:", sig_unique_ids, "\n", sep = "")
cat("Duplicated all IDs:    ", duplicated_all_ids, "\n", sep = "")
cat("Duplicated sig IDs:    ", duplicated_sig_ids, "\n", sep = "")

stopifnot(
  all_unique_ids == nrow(all),
  sig_unique_ids == nrow(sig),
  duplicated_all_ids == 0L,
  duplicated_sig_ids == 0L
)

# ------------------------------------------------------------
# 4. Numerical integrity
# ------------------------------------------------------------

missing_nes <- sum(is.na(all$NES))
missing_raw_p <- sum(is.na(all$pvalue))
missing_adj_p <- sum(is.na(all$p.adjust))

nonfinite_nes <- sum(!is.finite(all$NES))
nonfinite_raw_p <- sum(!is.finite(all$pvalue))
nonfinite_adj_p <- sum(!is.finite(all$p.adjust))

raw_p_valid <-
  all(all$pvalue >= 0 & all$pvalue <= 1)

adj_p_valid <-
  all(all$p.adjust >= 0 & all$p.adjust <= 1)

positive_nes <- sum(all$NES > 0)
negative_nes <- sum(all$NES < 0)
zero_nes <- sum(all$NES == 0)

cat("\n=== NUMERICAL INTEGRITY ===\n")
cat("Missing NES:           ", missing_nes, "\n", sep = "")
cat("Missing raw P:         ", missing_raw_p, "\n", sep = "")
cat("Missing adjusted P:    ", missing_adj_p, "\n", sep = "")
cat("Non-finite NES:        ", nonfinite_nes, "\n", sep = "")
cat("Raw P valid [0,1]:     ", raw_p_valid, "\n", sep = "")
cat("Adjusted P valid [0,1]:", adj_p_valid, "\n", sep = "")
cat("Positive NES:          ", positive_nes, "\n", sep = "")
cat("Negative NES:          ", negative_nes, "\n", sep = "")
cat("Zero NES:              ", zero_nes, "\n", sep = "")

stopifnot(
  missing_nes == 0L,
  missing_raw_p == 0L,
  missing_adj_p == 0L,
  nonfinite_nes == 0L,
  nonfinite_raw_p == 0L,
  nonfinite_adj_p == 0L,
  raw_p_valid,
  adj_p_valid,
  positive_nes == 3099L,
  negative_nes == 2724L,
  zero_nes == 0L
)

# ------------------------------------------------------------
# 5. Gene-set-size integrity
# ------------------------------------------------------------

min_set_size <- min(all$setSize)
max_set_size <- max(all$setSize)
set_size_valid <-
  all(all$setSize >= 10 & all$setSize <= 500)

cat("\n=== GENE-SET SIZE ===\n")
cat("Minimum setSize:       ", min_set_size, "\n", sep = "")
cat("Maximum setSize:       ", max_set_size, "\n", sep = "")
cat("All within [10,500]:   ", set_size_valid, "\n", sep = "")

stopifnot(set_size_valid)

# ------------------------------------------------------------
# 6. Significant-table identity
# ------------------------------------------------------------

expected_sig <- all[all$p.adjust < 0.05, , drop = FALSE]

significant_count_identical <-
  nrow(expected_sig) == nrow(sig)

significant_ids_identical <-
  setequal(expected_sig$ID, sig$ID)

sig_all_below_threshold <-
  all(sig$p.adjust < 0.05)

cat("\n=== SIGNIFICANCE FILTER ===\n")
cat("Computed significant:  ", nrow(expected_sig), "\n", sep = "")
cat("Written significant:   ", nrow(sig), "\n", sep = "")
cat("Count identical:       ", significant_count_identical, "\n", sep = "")
cat("IDs identical:         ", significant_ids_identical, "\n", sep = "")
cat("All sig adj.P < 0.05:  ", sig_all_below_threshold, "\n", sep = "")

stopifnot(
  significant_count_identical,
  significant_ids_identical,
  sig_all_below_threshold,
  nrow(sig) == 1054L
)

# ------------------------------------------------------------
# 7. Summary consistency
# ------------------------------------------------------------

summary_consistent <- all(c(
  num_summary("ranked_input_genes") == rank_rows,
  num_summary("positive_moderated_t") == rank_positive_t,
  num_summary("negative_moderated_t") == rank_negative_t,
  num_summary("zero_moderated_t") == rank_zero_t,
  num_summary("GO_BP_terms_returned") == nrow(all),
  num_summary("BH_adjusted_P_lt_0.05") == nrow(sig),
  num_summary("positive_NES_terms") == positive_nes,
  num_summary("negative_NES_terms") == negative_nes,
  num_summary("missing_NES") == missing_nes,
  num_summary("missing_raw_P") == missing_raw_p,
  num_summary("missing_adjusted_P") == missing_adj_p,
  num_summary("minGSSize") == 10,
  num_summary("maxGSSize") == 500,
  num_summary("seed") == 20261006
))

min_raw_matches <-
  isTRUE(all.equal(
    min(all$pvalue),
    num_summary("minimum_raw_P"),
    tolerance = 1e-12
  ))

min_adj_matches <-
  isTRUE(all.equal(
    min(all$p.adjust),
    num_summary("minimum_adjusted_P"),
    tolerance = 1e-12
  ))

cat("\n=== SUMMARY CONSISTENCY ===\n")
cat("Count metrics consistent: ", summary_consistent, "\n", sep = "")
cat("Minimum raw P matches:    ", min_raw_matches, "\n", sep = "")
cat("Minimum adjusted P matches:", min_adj_matches, "\n", sep = "")

stopifnot(
  summary_consistent,
  min_raw_matches,
  min_adj_matches
)

# ------------------------------------------------------------
# 8. Audit output
# ------------------------------------------------------------

qc <- data.frame(
  metric = c(
    "ranked_input_rows",
    "ranked_unique_gene_id",
    "ranked_unique_ENTREZID",
    "positive_moderated_t",
    "negative_moderated_t",
    "zero_moderated_t",
    "GO_BP_terms_returned",
    "BH_adjusted_P_lt_0.05",
    "positive_NES_terms",
    "negative_NES_terms",
    "missing_NES",
    "missing_raw_P",
    "missing_adjusted_P",
    "minimum_setSize",
    "maximum_setSize",
    "significant_IDs_identical",
    "summary_consistent",
    "minimum_raw_P_matches",
    "minimum_adjusted_P_matches"
  ),
  value = c(
    rank_rows,
    rank_unique_gene_id,
    rank_unique_entrez,
    rank_positive_t,
    rank_negative_t,
    rank_zero_t,
    nrow(all),
    nrow(sig),
    positive_nes,
    negative_nes,
    missing_nes,
    missing_raw_p,
    missing_adj_p,
    min_set_size,
    max_set_size,
    significant_ids_identical,
    summary_consistent,
    min_raw_matches,
    min_adj_matches
  ),
  stringsAsFactors = FALSE
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE
)

cat("\n=== OUTPUT ===\n")
cat(qc_file, "\n")

cat("\nIMPORTANT:\n")
cat("No pathway descriptions were printed.\n")
cat("No biological interpretation was performed.\n")
cat("No W5 differential-expression statistics were recalculated.\n")
cat("The frozen ranked input remained unchanged.\n")

cat("\n============================================================\n")
cat("W6.12d GO:BP RANKED ENRICHMENT RESULT QC: PASSED\n")
cat("============================================================\n")

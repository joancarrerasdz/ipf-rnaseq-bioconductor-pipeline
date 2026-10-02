# ============================================================
# GSE213001 — W5.5 age sensitivity analysis
#
# Frozen framework:
#   Primary biological contrast: IPF vs NDC
#   Frozen gene universe: 15,012 genes
#   Sensitivity cohort: age-complete primary cohort
#   Same-cohort models:
#
#     Reduced:
#       ~ diseasegroup + lunglocation
#
#     Age-adjusted:
#       ~ diseasegroup + lunglocation + age
#
#   Repeated measurements:
#       block = donorid
#
#   Normalization:
#       TMM recalculated in the 99-sample cohort
#
#   Inference:
#       two-pass voom + duplicateCorrelation
#       eBayes(robust = TRUE, trend = FALSE)
#       Benjamini-Hochberg FDR
#
# IMPORTANT:
#   This implementation was written before inspection of
#   age-adjusted gene-level results.
# ============================================================

source("renv/activate.R")

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

cat("\n============================================================\n")
cat("GSE213001 W5.5 — AGE SENSITIVITY ANALYSIS\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------

metadata <- read.csv(
  "metadata/sample_metadata.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_filter <- read.csv(
  "metadata/qc_gene_filtering.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

counts_df <- read.csv(
  gzfile(
    "data/source/GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv.gz"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

gene_id <- counts_df[[1]]

counts <- as.matrix(
  counts_df[, -1, drop = FALSE]
)

rownames(counts) <- gene_id
storage.mode(counts) <- "numeric"

# ------------------------------------------------------------
# 2. Frozen primary cohort
# ------------------------------------------------------------

primary <- metadata[
  metadata$diseasegroup %in% c("IPF", "NDC") &
    !is.na(metadata$lunglocation),
  ,
  drop = FALSE
]

stopifnot(
  nrow(primary) == 101,
  length(unique(primary$donorid)) == 34,
  sum(primary$diseasegroup == "IPF") == 61,
  sum(primary$diseasegroup == "NDC") == 40
)

# ------------------------------------------------------------
# 3. Age-complete sensitivity cohort
# ------------------------------------------------------------

sensitivity <- primary[
  !is.na(primary$age),
  ,
  drop = FALSE
]

sensitivity$age <- as.numeric(sensitivity$age)

stopifnot(
  nrow(sensitivity) == 99,
  length(unique(sensitivity$donorid)) == 33,
  !anyNA(sensitivity$age),
  !anyDuplicated(sensitivity$sample_title),
  all(sensitivity$sample_title %in% colnames(counts))
)

missing_age_donors <- unique(
  primary$donorid[is.na(primary$age)]
)

stopifnot(
  length(missing_age_donors) == 1,
  identical(missing_age_donors, "ALF017")
)

# Explicit reference levels.

sensitivity$diseasegroup <- factor(
  sensitivity$diseasegroup,
  levels = c("NDC", "IPF")
)

sensitivity$lunglocation <- factor(
  sensitivity$lunglocation,
  levels = c("Base", "Apex")
)

# ------------------------------------------------------------
# 4. Frozen W4 gene universe
# ------------------------------------------------------------

stopifnot(
  nrow(gene_filter) == 15065,
  "gene_id" %in% names(gene_filter),
  "keep_primary" %in% names(gene_filter)
)

keep_ids <- gene_filter$gene_id[
  gene_filter$keep_primary
]

stopifnot(
  length(keep_ids) == 15012,
  !anyDuplicated(keep_ids),
  all(keep_ids %in% rownames(counts))
)

counts_sensitivity <- counts[
  keep_ids,
  sensitivity$sample_title,
  drop = FALSE
]

stopifnot(
  nrow(counts_sensitivity) == 15012,
  ncol(counts_sensitivity) == 99,
  identical(
    colnames(counts_sensitivity),
    sensitivity$sample_title
  )
)

cat("\n=== SENSITIVITY COHORT ===\n")
cat("Samples:                 ", nrow(sensitivity), "\n", sep = "")
cat(
  "Donors:                  ",
  length(unique(sensitivity$donorid)),
  "\n",
  sep = ""
)
cat("Removed missing-age donor: ALF017\n")
cat("Frozen retained genes:    ", nrow(counts_sensitivity), "\n", sep = "")

# ------------------------------------------------------------
# 5. TMM normalization on identical 99-sample cohort
# ------------------------------------------------------------

dge <- DGEList(
  counts = counts_sensitivity
)

dge <- calcNormFactors(
  dge,
  method = "TMM"
)

stopifnot(
  ncol(dge) == 99,
  nrow(dge) == 15012,
  all(is.finite(dge$samples$norm.factors))
)

cat("\n=== TMM ===\n")
print(summary(dge$samples$norm.factors))
cat(
  "Product of normalization factors: ",
  signif(prod(dge$samples$norm.factors), 7),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 6. Designs
# ------------------------------------------------------------

design_reduced <- model.matrix(
  ~ diseasegroup + lunglocation,
  data = sensitivity
)

design_age <- model.matrix(
  ~ diseasegroup + lunglocation + age,
  data = sensitivity
)

stopifnot(
  qr(design_reduced)$rank == ncol(design_reduced),
  qr(design_age)$rank == ncol(design_age),
  "diseasegroupIPF" %in% colnames(design_reduced),
  "diseasegroupIPF" %in% colnames(design_age)
)

block <- sensitivity$donorid

cat("\n=== DESIGN MATRICES ===\n")
cat(
  "Reduced:      ",
  paste(colnames(design_reduced), collapse = " | "),
  "\n",
  sep = ""
)
cat(
  "Rank:         ",
  qr(design_reduced)$rank,
  "/",
  ncol(design_reduced),
  "\n",
  sep = ""
)

cat(
  "Age-adjusted: ",
  paste(colnames(design_age), collapse = " | "),
  "\n",
  sep = ""
)
cat(
  "Rank:         ",
  qr(design_age)$rank,
  "/",
  ncol(design_age),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 7. Donor-aware model helper
# ------------------------------------------------------------

fit_donor_model <- function(dge, design, block) {

  v1 <- voom(
    dge,
    design,
    plot = FALSE
  )

  dc1 <- duplicateCorrelation(
    v1,
    design,
    block = block
  )

  v2 <- voom(
    dge,
    design,
    plot = FALSE,
    block = block,
    correlation = dc1$consensus.correlation
  )

  dc2 <- duplicateCorrelation(
    v2,
    design,
    block = block
  )

  fit <- lmFit(
    v2,
    design,
    block = block,
    correlation = dc2$consensus.correlation
  )

  fit <- eBayes(
    fit,
    robust = TRUE,
    trend = FALSE
  )

  list(
    v1 = v1,
    v2 = v2,
    rho1 = dc1$consensus.correlation,
    rho2 = dc2$consensus.correlation,
    fit = fit
  )
}

# ------------------------------------------------------------
# 8. Same-cohort reduced model
# ------------------------------------------------------------

cat("\n=== REDUCED 99-SAMPLE MODEL ===\n")

reduced_fit <- fit_donor_model(
  dge = dge,
  design = design_reduced,
  block = block
)

cat(
  "duplicateCorrelation 1: ",
  signif(reduced_fit$rho1, 7),
  "\n",
  sep = ""
)

cat(
  "duplicateCorrelation 2: ",
  signif(reduced_fit$rho2, 7),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 9. Age-adjusted model
# ------------------------------------------------------------

cat("\n=== AGE-ADJUSTED 99-SAMPLE MODEL ===\n")

age_fit <- fit_donor_model(
  dge = dge,
  design = design_age,
  block = block
)

cat(
  "duplicateCorrelation 1: ",
  signif(age_fit$rho1, 7),
  "\n",
  sep = ""
)

cat(
  "duplicateCorrelation 2: ",
  signif(age_fit$rho2, 7),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# 10. Complete gene-level tables
# ------------------------------------------------------------

extract_table <- function(fit_object) {

  x <- topTable(
    fit_object,
    coef = "diseasegroupIPF",
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  )

  x$gene_id <- rownames(x)

  x <- x[
    ,
    c(
      "gene_id",
      setdiff(names(x), "gene_id")
    ),
    drop = FALSE
  ]

  rownames(x) <- NULL

  x
}

de_reduced <- extract_table(
  reduced_fit$fit
)

de_age <- extract_table(
  age_fit$fit
)

stopifnot(
  nrow(de_reduced) == 15012,
  nrow(de_age) == 15012,
  identical(de_reduced$gene_id, de_age$gene_id),
  all(is.finite(de_reduced$logFC)),
  all(is.finite(de_age$logFC)),
  all(is.finite(de_reduced$adj.P.Val)),
  all(is.finite(de_age$adj.P.Val))
)

# ------------------------------------------------------------
# 11. Same-cohort concordance
# ------------------------------------------------------------

logfc_reduced <- de_reduced$logFC
logfc_age <- de_age$logFC

sig_reduced <- de_reduced$adj.P.Val < 0.05
sig_age <- de_age$adj.P.Val < 0.05

direction_concordant <-
  sign(logfc_reduced) == sign(logfc_age)

pearson_logfc <- cor(
  logfc_reduced,
  logfc_age,
  method = "pearson"
)

spearman_logfc <- cor(
  logfc_reduced,
  logfc_age,
  method = "spearman"
)

direction_pct <-
  100 * mean(direction_concordant)

median_abs_delta <-
  median(
    abs(logfc_age - logfc_reduced)
  )

n_sig_reduced <- sum(sig_reduced)
n_sig_age <- sum(sig_age)

n_sig_overlap <- sum(
  sig_reduced & sig_age
)

n_sig_reduced_only <- sum(
  sig_reduced & !sig_age
)

n_sig_age_only <- sum(
  !sig_reduced & sig_age
)

fdr_union <- sum(
  sig_reduced | sig_age
)

fdr_jaccard <-
  if (fdr_union == 0) {
    NA_real_
  } else {
    n_sig_overlap / fdr_union
  }

# Pre-specified |logFC| >= 1 robustness summaries.
#
# Report both:
#   1) union: >=1 in either model
#   2) intersection: >=1 in both models
#
# This avoids selecting one model post hoc as the denominator.

strong_union <-
  abs(logfc_reduced) >= 1 |
  abs(logfc_age) >= 1

strong_both <-
  abs(logfc_reduced) >= 1 &
  abs(logfc_age) >= 1

strong_union_direction_pct <-
  if (sum(strong_union) == 0) {
    NA_real_
  } else {
    100 * mean(
      direction_concordant[strong_union]
    )
  }

strong_both_direction_pct <-
  if (sum(strong_both) == 0) {
    NA_real_
  } else {
    100 * mean(
      direction_concordant[strong_both]
    )
  }

# ------------------------------------------------------------
# 12. Gene-level concordance file
# ------------------------------------------------------------

concordance <- data.frame(
  gene_id = de_reduced$gene_id,
  logFC_reduced = logfc_reduced,
  logFC_age_adjusted = logfc_age,
  delta_logFC_age_minus_reduced =
    logfc_age - logfc_reduced,
  abs_delta_logFC =
    abs(logfc_age - logfc_reduced),
  direction_concordant =
    direction_concordant,
  FDR_reduced =
    de_reduced$adj.P.Val,
  FDR_age_adjusted =
    de_age$adj.P.Val,
  significant_reduced =
    sig_reduced,
  significant_age_adjusted =
    sig_age,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 13. Summary table
# ------------------------------------------------------------

summary_table <- data.frame(
  metric = c(
    "samples",
    "donors",
    "tested_genes",
    "rho2_reduced",
    "rho2_age_adjusted",
    "pearson_logFC",
    "spearman_logFC",
    "direction_concordance_percent",
    "median_absolute_logFC_change",
    "FDR_lt_0.05_reduced",
    "FDR_lt_0.05_age_adjusted",
    "FDR_overlap",
    "FDR_reduced_only",
    "FDR_age_adjusted_only",
    "FDR_jaccard",
    "strong_union_genes",
    "strong_union_direction_concordance_percent",
    "strong_both_genes",
    "strong_both_direction_concordance_percent"
  ),
  value = c(
    nrow(sensitivity),
    length(unique(sensitivity$donorid)),
    nrow(de_reduced),
    reduced_fit$rho2,
    age_fit$rho2,
    pearson_logfc,
    spearman_logfc,
    direction_pct,
    median_abs_delta,
    n_sig_reduced,
    n_sig_age,
    n_sig_overlap,
    n_sig_reduced_only,
    n_sig_age_only,
    fdr_jaccard,
    sum(strong_union),
    strong_union_direction_pct,
    sum(strong_both),
    strong_both_direction_pct
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# 14. Outputs
# ------------------------------------------------------------

dir.create(
  "results/differential_expression",
  recursive = TRUE,
  showWarnings = FALSE
)

reduced_file <-
  "results/differential_expression/sensitivity_99_reduced_all_genes.csv"

age_file <-
  "results/differential_expression/sensitivity_99_age_adjusted_all_genes.csv"

concordance_file <-
  "results/differential_expression/sensitivity_age_gene_concordance.csv"

summary_file <-
  "results/differential_expression/sensitivity_age_summary.csv"

write.csv(
  de_reduced,
  reduced_file,
  row.names = FALSE
)

write.csv(
  de_age,
  age_file,
  row.names = FALSE
)

write.csv(
  concordance,
  concordance_file,
  row.names = FALSE
)

write.csv(
  summary_table,
  summary_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. Global-only reporting
# ------------------------------------------------------------

cat("\n=== SAME-COHORT AGE SENSITIVITY SUMMARY ===\n")

cat(
  "Pearson logFC correlation:           ",
  sprintf("%.6f", pearson_logfc),
  "\n",
  sep = ""
)

cat(
  "Spearman logFC correlation:          ",
  sprintf("%.6f", spearman_logfc),
  "\n",
  sep = ""
)

cat(
  "Direction concordance:               ",
  sprintf("%.2f%%", direction_pct),
  "\n",
  sep = ""
)

cat(
  "Median absolute logFC change:        ",
  sprintf("%.6f", median_abs_delta),
  "\n",
  sep = ""
)

cat(
  "FDR < 0.05 — reduced model:          ",
  n_sig_reduced,
  "\n",
  sep = ""
)

cat(
  "FDR < 0.05 — age-adjusted model:     ",
  n_sig_age,
  "\n",
  sep = ""
)

cat(
  "FDR-significant overlap:             ",
  n_sig_overlap,
  "\n",
  sep = ""
)

cat(
  "Significant only before age:         ",
  n_sig_reduced_only,
  "\n",
  sep = ""
)

cat(
  "Significant only after age:          ",
  n_sig_age_only,
  "\n",
  sep = ""
)

cat(
  "FDR-set Jaccard index:               ",
  sprintf("%.6f", fdr_jaccard),
  "\n",
  sep = ""
)

cat(
  "|logFC| >= 1 in either model:        ",
  sum(strong_union),
  "\n",
  sep = ""
)

cat(
  "Direction concordance in union:      ",
  sprintf("%.2f%%", strong_union_direction_pct),
  "\n",
  sep = ""
)

cat(
  "|logFC| >= 1 in both models:         ",
  sum(strong_both),
  "\n",
  sep = ""
)

cat(
  "Direction concordance in both:       ",
  sprintf("%.2f%%", strong_both_direction_pct),
  "\n",
  sep = ""
)

cat("\n=== OUTPUTS ===\n")
cat(reduced_file, "\n")
cat(age_file, "\n")
cat(concordance_file, "\n")
cat(summary_file, "\n")

cat("\nIMPORTANT:\n")
cat("No individual gene identities were printed by this analysis.\n")
cat("Biological interpretation remains deferred.\n")

cat("\n============================================================\n")
cat("W5.5 AGE SENSITIVITY ANALYSIS: COMPLETED\n")
cat("============================================================\n")

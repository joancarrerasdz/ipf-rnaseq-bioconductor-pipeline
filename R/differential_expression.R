# ============================================================
# GSE213001 — Differential expression
#
# Primary analysis:
#   IPF vs NDC
#
# Fixed effects:
#   ~ diseasegroup + lunglocation
#
# Repeated measurements:
#   donorid handled with limma::duplicateCorrelation()
#
# Normalization:
#   TMM
#
# Filtering:
#   Frozen W4 keep_primary rule
#
# Sensitivity:
#   ~ diseasegroup + lunglocation + age
#
# IMPORTANT:
#   Statistical implementation was frozen before inspection
#   of formal gene-level differential-expression results.
# ============================================================

source("renv/activate.R")

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

# ------------------------------------------------------------
# Inputs
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

storage.mode(counts) <- "integer"
rownames(counts) <- gene_id

# ------------------------------------------------------------
# Frozen W4 filter
# ------------------------------------------------------------

stopifnot(
  nrow(gene_filter) == 15065,
  "gene_id" %in% names(gene_filter),
  "keep_primary" %in% names(gene_filter),
  is.logical(gene_filter$keep_primary),
  !anyNA(gene_filter$keep_primary),
  sum(gene_filter$keep_primary) == 15012,
  sum(!gene_filter$keep_primary) == 53
)

# ------------------------------------------------------------
# Primary cohort
# ------------------------------------------------------------

primary <- metadata[
  metadata$diseasegroup %in% c("IPF", "NDC") &
    !is.na(metadata$lunglocation) &
    metadata$lunglocation %in% c("Apex", "Base"),
  ,
  drop = FALSE
]

primary$diseasegroup <- relevel(
  factor(primary$diseasegroup),
  ref = "NDC"
)

primary$lunglocation <- relevel(
  factor(primary$lunglocation),
  ref = "Base"
)

stopifnot(
  nrow(primary) == 101,
  length(unique(primary$donorid)) == 34,
  sum(primary$diseasegroup == "IPF") == 61,
  sum(primary$diseasegroup == "NDC") == 40,
  sum(primary$lunglocation == "Apex") == 52,
  sum(primary$lunglocation == "Base") == 49
)

# ------------------------------------------------------------
# Align counts to primary metadata
# ------------------------------------------------------------

stopifnot(
  all(primary$sample_title %in% colnames(counts))
)

counts_primary <- counts[
  ,
  primary$sample_title,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(counts_primary),
    primary$sample_title
  )
)

idx <- match(
  gene_filter$gene_id,
  rownames(counts_primary)
)

stopifnot(!anyNA(idx))

counts_primary <- counts_primary[
  idx,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(counts_primary),
    gene_filter$gene_id
  )
)

counts_filtered <- counts_primary[
  gene_filter$keep_primary,
  ,
  drop = FALSE
]

stopifnot(
  nrow(counts_filtered) == 15012,
  ncol(counts_filtered) == 101
)

# ------------------------------------------------------------
# Primary design
# ------------------------------------------------------------

design <- model.matrix(
  ~ diseasegroup + lunglocation,
  data = primary
)

block <- factor(primary$donorid)

stopifnot(
  qr(design)$rank == ncol(design),
  identical(
    colnames(design),
    c(
      "(Intercept)",
      "diseasegroupIPF",
      "lunglocationApex"
    )
  )
)

# ------------------------------------------------------------
# TMM
# ------------------------------------------------------------

y <- DGEList(
  counts = counts_filtered
)

y <- calcNormFactors(
  y,
  method = "TMM"
)

stopifnot(
  all(is.finite(y$samples$norm.factors)),
  all(y$samples$norm.factors > 0)
)

# ------------------------------------------------------------
# Donor-aware voom / duplicateCorrelation
# ------------------------------------------------------------

v1 <- voom(
  y,
  design = design,
  plot = FALSE
)

corfit1 <- duplicateCorrelation(
  v1,
  design = design,
  block = block
)

rho1 <- corfit1$consensus.correlation

v2 <- voom(
  y,
  design = design,
  plot = FALSE,
  block = block,
  correlation = rho1
)

corfit2 <- duplicateCorrelation(
  v2,
  design = design,
  block = block
)

rho2 <- corfit2$consensus.correlation

stopifnot(
  is.finite(rho1),
  is.finite(rho2),
  rho1 > -1,
  rho1 < 1,
  rho2 > -1,
  rho2 < 1
)

# Dry-run value observed before formal DE:
# rho1 ≈ 0.315821
# rho2 ≈ 0.315828
#
# The final inferential model will use:
#
#   expression ~ diseasegroup + lunglocation
#
# with:
#
#   block       = donorid
#   correlation = rho2
#
# Primary coefficient:
#
#   diseasegroupIPF
#
# Interpretation:
#
#   IPF - NDC, adjusted for lung region.
#
# Formal lmFit/eBayes code is intentionally added only
# after this implementation has been committed/frozen.

cat("\n============================================================\n")
cat("W5 DE FRAMEWORK PRE-INFERENCE CHECK: PASSED\n")
cat("============================================================\n")
cat("Genes:                 ", nrow(counts_filtered), "\n", sep = "")
cat("Samples:               ", ncol(counts_filtered), "\n", sep = "")
cat("Donors:                ", length(unique(primary$donorid)), "\n", sep = "")
cat("Design:                 ~ diseasegroup + lunglocation\n")
cat("Primary coefficient:    diseasegroupIPF\n")
cat("Blocking unit:          donorid\n")
cat("duplicateCorrelation 1: ", signif(rho1, 7), "\n", sep = "")
cat("duplicateCorrelation 2: ", signif(rho2, 7), "\n", sep = "")
cat("lmFit/eBayes:           NOT PERFORMED\n")
cat("============================================================\n")

# ============================================================
# W5.4 — FORMAL PRIMARY DIFFERENTIAL-EXPRESSION ANALYSIS
# ============================================================

cat("\n============================================================\n")
cat("GSE213001 W5.4 — PRIMARY DONOR-AWARE DE\n")
cat("============================================================\n")

# ------------------------------------------------------------
# 1. Formal donor-aware model fit
# ------------------------------------------------------------

fit <- lmFit(
  v2,
  design,
  block = block,
  correlation = rho2
)

fit <- eBayes(
  fit,
  robust = TRUE,
  trend = FALSE
)

stopifnot(
  "diseasegroupIPF" %in% colnames(fit$coefficients),
  nrow(fit$coefficients) == 15012
)

# ------------------------------------------------------------
# 2. Complete primary result table
#
# IMPORTANT:
# sort.by = "none" preserves the tested-gene ordering.
# No gene ranking is printed at this stage.
# ------------------------------------------------------------

de_primary <- topTable(
  fit,
  coef = "diseasegroupIPF",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

de_primary$gene_id <- rownames(de_primary)

de_primary <- de_primary[
  ,
  c(
    "gene_id",
    setdiff(names(de_primary), "gene_id")
  ),
  drop = FALSE
]

rownames(de_primary) <- NULL

stopifnot(
  nrow(de_primary) == 15012,
  !anyDuplicated(de_primary$gene_id),
  all(is.finite(de_primary$logFC)),
  all(is.finite(de_primary$P.Value)),
  all(is.finite(de_primary$adj.P.Val)),
  all(de_primary$P.Value >= 0 & de_primary$P.Value <= 1),
  all(de_primary$adj.P.Val >= 0 & de_primary$adj.P.Val <= 1)
)

# ------------------------------------------------------------
# 3. Global inferential summaries
# ------------------------------------------------------------

n_tested <- nrow(de_primary)

n_raw_005 <- sum(
  de_primary$P.Value < 0.05
)

n_fdr_005 <- sum(
  de_primary$adj.P.Val < 0.05
)

n_fdr_up <- sum(
  de_primary$adj.P.Val < 0.05 &
    de_primary$logFC > 0
)

n_fdr_down <- sum(
  de_primary$adj.P.Val < 0.05 &
    de_primary$logFC < 0
)

n_abs_logfc_1 <- sum(
  abs(de_primary$logFC) >= 1
)

n_fdr_abs_logfc_1 <- sum(
  de_primary$adj.P.Val < 0.05 &
    abs(de_primary$logFC) >= 1
)

# ------------------------------------------------------------
# 4. Result directories
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  "results/differential_expression",
  showWarnings = FALSE,
  recursive = TRUE
)

# ------------------------------------------------------------
# 5. Write complete result table
# ------------------------------------------------------------

primary_result_file <-
  "results/differential_expression/primary_IPF_vs_NDC_all_genes.csv"

write.csv(
  de_primary,
  primary_result_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Write model-level summary
# ------------------------------------------------------------

de_summary <- data.frame(
  metric = c(
    "tested_genes",
    "primary_samples",
    "primary_donors",
    "consensus_correlation",
    "raw_p_lt_0.05",
    "fdr_lt_0.05",
    "fdr_lt_0.05_logFC_positive",
    "fdr_lt_0.05_logFC_negative",
    "abs_logFC_ge_1",
    "fdr_lt_0.05_and_abs_logFC_ge_1",
    "min_raw_p",
    "min_fdr"
  ),
  value = c(
    n_tested,
    ncol(v2$E),
    length(unique(primary$donorid)),
    rho2,
    n_raw_005,
    n_fdr_005,
    n_fdr_up,
    n_fdr_down,
    n_abs_logfc_1,
    n_fdr_abs_logfc_1,
    min(de_primary$P.Value),
    min(de_primary$adj.P.Val)
  ),
  stringsAsFactors = FALSE
)

summary_file <-
  "results/differential_expression/primary_IPF_vs_NDC_summary.csv"

write.csv(
  de_summary,
  summary_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Final model-level gate
# ------------------------------------------------------------

stopifnot(
  n_fdr_up + n_fdr_down <= n_fdr_005
)

cat("\n=== MODEL ===\n")
cat("Contrast:               IPF vs NDC\n")
cat("Coefficient:            diseasegroupIPF\n")
cat("Adjusted for:           lunglocation\n")
cat("Blocking unit:          donorid\n")
cat(
  "Consensus correlation: ",
  signif(rho2, 7),
  "\n",
  sep = ""
)
cat("Empirical Bayes:        robust=TRUE, trend=FALSE\n")
cat("Multiple testing:       Benjamini-Hochberg\n")

cat("\n=== GLOBAL INFERENCE SUMMARY ===\n")
cat("Genes tested:                    ", n_tested, "\n", sep = "")
cat("Raw P < 0.05:                    ", n_raw_005, "\n", sep = "")
cat("FDR < 0.05:                      ", n_fdr_005, "\n", sep = "")
cat("  positive logFC:                ", n_fdr_up, "\n", sep = "")
cat("  negative logFC:                ", n_fdr_down, "\n", sep = "")
cat("|logFC| >= 1:                    ", n_abs_logfc_1, "\n", sep = "")
cat(
  "FDR < 0.05 and |logFC| >= 1:    ",
  n_fdr_abs_logfc_1,
  "\n",
  sep = ""
)
cat(
  "Minimum raw P:                  ",
  format(min(de_primary$P.Value), scientific = TRUE, digits = 4),
  "\n",
  sep = ""
)
cat(
  "Minimum FDR:                    ",
  format(min(de_primary$adj.P.Val), scientific = TRUE, digits = 4),
  "\n",
  sep = ""
)

cat("\n=== OUTPUTS ===\n")
cat(primary_result_file, "\n")
cat(summary_file, "\n")

cat("\nIMPORTANT:\n")
cat("No individual gene identities were printed by this stage.\n")
cat("Biological gene-level interpretation has not yet been performed.\n")

cat("\n============================================================\n")
cat("W5.4 PRIMARY DONOR-AWARE DE: COMPLETED\n")
cat("============================================================\n")

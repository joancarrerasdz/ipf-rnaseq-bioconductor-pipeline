# ============================================================
# GSE213001 — W7 gender sensitivity analysis
#
# Frozen W5 reference:
#   IPF vs NDC
#
# W7 sensitivity design:
#   ~ 0 + group + gender
#
# Cohort:
#   n = 101 samples
#   n = 34 donors
#
# Gender:
#   F = 39 samples
#   M = 62 samples
#   no missing values
#
# Repeated measurements:
#   donorid handled with limma::duplicateCorrelation()
#
# Normalization:
#   TMM
#
# Gene universe:
#   frozen W4/W5 keep_primary = 15,012 genes
#
# IMPORTANT:
#   This is a W7 sensitivity model only.
#   No W5, W6, or previous W7 output is overwritten.
# ============================================================

source("renv/activate.R")

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

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

storage.mode(counts) <- "integer"
rownames(counts) <- gene_id

# ------------------------------------------------------------
# 2. Frozen gene filter
# ------------------------------------------------------------

stopifnot(
  nrow(gene_filter) == 15065L,
  "gene_id" %in% names(gene_filter),
  "keep_primary" %in% names(gene_filter),
  is.logical(gene_filter$keep_primary),
  !anyNA(gene_filter$keep_primary),
  sum(gene_filter$keep_primary) == 15012L,
  sum(!gene_filter$keep_primary) == 53L
)

# ------------------------------------------------------------
# 3. Reconstruct frozen W5 primary cohort
# ------------------------------------------------------------

primary <- metadata[
  metadata$diseasegroup %in% c("IPF", "NDC") &
    !is.na(metadata$lunglocation) &
    metadata$lunglocation %in% c("Apex", "Base"),
  ,
  drop = FALSE
]

stopifnot(
  nrow(primary) == 101L,
  length(unique(primary$donorid)) == 34L
)

# ------------------------------------------------------------
# 4. Gender audit and frozen encoding
# ------------------------------------------------------------

stopifnot(
  "gender" %in% names(primary),
  sum(is.na(primary$gender) | primary$gender == "") == 0L
)

primary$gender <- factor(
  primary$gender,
  levels = c("F", "M")
)

stopifnot(
  !anyNA(primary$gender),
  identical(
    as.integer(table(primary$gender)),
    c(39L, 62L)
  )
)

gender_by_donor <- tapply(
  as.character(primary$gender),
  primary$donorid,
  function(z) length(unique(z[!is.na(z) & z != ""]))
)

stopifnot(
  length(gender_by_donor) == 34L,
  all(gender_by_donor == 1L)
)

# ------------------------------------------------------------
# 5. Frozen group encoding
# ------------------------------------------------------------

primary$group <- factor(
  primary$group,
  levels = c(
    "IPF.Apex",
    "IPF.Base",
    "NDC.Apex",
    "NDC.Base"
  )
)

stopifnot(
  !anyNA(primary$group),
  identical(
    as.integer(table(primary$group)),
    c(30L, 31L, 22L, 18L)
  )
)

group_gender <- table(
  primary$group,
  primary$gender
)

stopifnot(
  identical(
    unname(unclass(group_gender)),
    matrix(
      c(
        11L, 11L, 9L, 8L,
        19L, 20L, 13L, 10L
      ),
      nrow = 4L,
      ncol = 2L
    )
  )
)

# ------------------------------------------------------------
# 6. Align counts to W7 gender cohort
# ------------------------------------------------------------

stopifnot(
  all(primary$sample_title %in% colnames(counts))
)

counts_gender <- counts[
  ,
  primary$sample_title,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(counts_gender),
    primary$sample_title
  )
)

idx <- match(
  gene_filter$gene_id,
  rownames(counts_gender)
)

stopifnot(!anyNA(idx))

counts_gender <- counts_gender[
  idx,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(counts_gender),
    gene_filter$gene_id
  )
)

counts_filtered <- counts_gender[
  gene_filter$keep_primary,
  ,
  drop = FALSE
]

stopifnot(
  nrow(counts_filtered) == 15012L,
  ncol(counts_filtered) == 101L
)

# ------------------------------------------------------------
# 7. Prespecified W7 gender design
# ------------------------------------------------------------

design <- model.matrix(
  ~ 0 + group + gender,
  data = primary
)

block <- factor(primary$donorid)

stopifnot(
  qr(design)$rank == ncol(design),
  identical(
    colnames(design),
    c(
      "groupIPF.Apex",
      "groupIPF.Base",
      "groupNDC.Apex",
      "groupNDC.Base",
      "genderM"
    )
  )
)

# ------------------------------------------------------------
# 8. Prespecified IPF-vs-NDC contrast
# ------------------------------------------------------------

contrast_matrix <- makeContrasts(
  IPF_vs_NDC =
    (
      groupIPF.Apex +
      groupIPF.Base
    ) / 2 -
    (
      groupNDC.Apex +
      groupNDC.Base
    ) / 2,
  levels = design
)

stopifnot(
  ncol(contrast_matrix) == 1L,
  colnames(contrast_matrix) == "IPF_vs_NDC"
)

# ------------------------------------------------------------
# 9. TMM normalization
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
# 10. Donor-aware voom / duplicateCorrelation
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

# ------------------------------------------------------------
# 11. Sensitivity-model fit
# ------------------------------------------------------------

fit <- lmFit(
  v2,
  design,
  block = block,
  correlation = rho2
)

fit <- contrasts.fit(
  fit,
  contrast_matrix
)

fit <- eBayes(
  fit,
  robust = TRUE,
  trend = FALSE
)

stopifnot(
  "IPF_vs_NDC" %in% colnames(fit$coefficients),
  nrow(fit$coefficients) == 15012L
)

# ------------------------------------------------------------
# 12. Complete result table
#
# Preserve frozen gene order. No ranking-based selection.
# ------------------------------------------------------------

de_gender <- topTable(
  fit,
  coef = "IPF_vs_NDC",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

de_gender$gene_id <- rownames(de_gender)

de_gender <- de_gender[
  ,
  c(
    "gene_id",
    setdiff(names(de_gender), "gene_id")
  ),
  drop = FALSE
]

rownames(de_gender) <- NULL

stopifnot(
  nrow(de_gender) == 15012L,
  !anyDuplicated(de_gender$gene_id),
  identical(
    de_gender$gene_id,
    gene_filter$gene_id[gene_filter$keep_primary]
  ),
  all(is.finite(de_gender$logFC)),
  all(is.finite(de_gender$P.Value)),
  all(is.finite(de_gender$adj.P.Val)),
  all(de_gender$P.Value >= 0 & de_gender$P.Value <= 1),
  all(de_gender$adj.P.Val >= 0 & de_gender$adj.P.Val <= 1)
)

# ------------------------------------------------------------
# 13. Global summaries
# ------------------------------------------------------------

n_tested <- nrow(de_gender)

n_raw_005 <- sum(
  de_gender$P.Value < 0.05
)

n_fdr_005 <- sum(
  de_gender$adj.P.Val < 0.05
)

n_fdr_up <- sum(
  de_gender$adj.P.Val < 0.05 &
    de_gender$logFC > 0
)

n_fdr_down <- sum(
  de_gender$adj.P.Val < 0.05 &
    de_gender$logFC < 0
)

n_abs_logfc_1 <- sum(
  abs(de_gender$logFC) >= 1
)

n_fdr_abs_logfc_1 <- sum(
  de_gender$adj.P.Val < 0.05 &
    abs(de_gender$logFC) >= 1
)

# ------------------------------------------------------------
# 14. Output directories
# ------------------------------------------------------------

dir.create(
  "results/sensitivity",
  showWarnings = FALSE,
  recursive = TRUE
)

result_file <-
  "results/sensitivity/W7_gender_IPF_vs_NDC_all_genes.csv"

summary_file <-
  "results/sensitivity/W7_gender_IPF_vs_NDC_summary.csv"

qc_file <-
  "metadata/qc_W7_gender_sensitivity.csv"

# ------------------------------------------------------------
# 15. Write outputs
# ------------------------------------------------------------

write.csv(
  de_gender,
  result_file,
  row.names = FALSE
)

de_summary <- data.frame(
  metric = c(
    "tested_genes",
    "analysis_samples",
    "analysis_donors",
    "female_samples",
    "male_samples",
    "duplicateCorrelation_pass1",
    "duplicateCorrelation_pass2",
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
    nrow(primary),
    length(unique(primary$donorid)),
    sum(primary$gender == "F"),
    sum(primary$gender == "M"),
    rho1,
    rho2,
    n_raw_005,
    n_fdr_005,
    n_fdr_up,
    n_fdr_down,
    n_abs_logfc_1,
    n_fdr_abs_logfc_1,
    min(de_gender$P.Value),
    min(de_gender$adj.P.Val)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  de_summary,
  summary_file,
  row.names = FALSE
)

qc <- data.frame(
  metric = c(
    "primary_reference_samples",
    "analysis_samples",
    "primary_reference_donors",
    "analysis_donors",
    "missing_gender_samples",
    "female_samples",
    "male_samples",
    "tested_genes",
    "design_columns",
    "design_rank",
    "gender_reference",
    "contrast"
  ),
  value = c(
    "101",
    "101",
    "34",
    "34",
    "0",
    "39",
    "62",
    "15012",
    as.character(ncol(design)),
    as.character(qr(design)$rank),
    "F",
    "0.5*(IPF.Apex+IPF.Base)-0.5*(NDC.Apex+NDC.Base)"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  qc,
  qc_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 16. Final gates
# ------------------------------------------------------------

stopifnot(
  n_fdr_up + n_fdr_down <= n_fdr_005,
  nrow(primary) == 101L,
  length(unique(primary$donorid)) == 34L,
  sum(primary$gender == "F") == 39L,
  sum(primary$gender == "M") == 62L,
  nrow(de_gender) == 15012L,
  qr(design)$rank == ncol(design)
)

# ------------------------------------------------------------
# 17. Console certificate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GSE213001 W7.3 — GENDER SENSITIVITY MODEL\n")
cat("============================================================\n")

cat("\n=== COHORT ===\n")
cat("Samples:                     ", nrow(primary), "\n", sep = "")
cat(
  "Donors:                      ",
  length(unique(primary$donorid)),
  "\n",
  sep = ""
)
cat("Female samples:               ", sum(primary$gender == "F"), "\n", sep = "")
cat("Male samples:                 ", sum(primary$gender == "M"), "\n", sep = "")
cat("Missing gender:               0\n")

cat("\n=== DESIGN ===\n")
cat("Formula:   ~ 0 + group + gender\n")
cat("Gender reference: F\n")
cat("Full rank: ", qr(design)$rank == ncol(design), "\n", sep = "")
cat("Block:     donorid\n")
cat(
  "Contrast:  0.5*(IPF.Apex + IPF.Base) - ",
  "0.5*(NDC.Apex + NDC.Base)\n",
  sep = ""
)

cat("\n=== CORRELATION ===\n")
cat("duplicateCorrelation pass 1: ", signif(rho1, 7), "\n", sep = "")
cat("duplicateCorrelation pass 2: ", signif(rho2, 7), "\n", sep = "")

cat("\n=== INFERENCE ===\n")
cat("Genes tested:                 ", n_tested, "\n", sep = "")
cat("Raw P < 0.05:                 ", n_raw_005, "\n", sep = "")
cat("FDR < 0.05:                   ", n_fdr_005, "\n", sep = "")
cat("FDR positive logFC:           ", n_fdr_up, "\n", sep = "")
cat("FDR negative logFC:           ", n_fdr_down, "\n", sep = "")

cat("\n=== OUTPUTS ===\n")
cat(result_file, "\n")
cat(summary_file, "\n")
cat(qc_file, "\n")

cat("\nIMPORTANT:\n")
cat("This is a W7 sensitivity analysis only.\n")
cat("The frozen W5 primary analysis remains the reference.\n")
cat("No W5 or W6 output is overwritten.\n")

cat("\n============================================================\n")
cat("W7.3 GENDER SENSITIVITY MODEL: COMPLETED\n")
cat("============================================================\n")

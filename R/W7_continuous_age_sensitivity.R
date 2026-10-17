# ============================================================
# GSE213001 — W7 continuous-age sensitivity analysis
#
# Frozen W5 reference:
#   IPF vs NDC
#
# W7 sensitivity design:
#   ~ 0 + group + age_numeric
#
# Continuous-age complete-case cohort:
#   n = 99 samples
#   n = 33 donors
#
# Missing continuous age:
#   GSM6568354 — ALF017 — NDC.Apex
#   GSM6568343 — ALF017 — NDC.Apex
#
# No age imputation.
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
#   No W5 or W6 output is overwritten.
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
  nrow(gene_filter) == 15065,
  "gene_id" %in% names(gene_filter),
  "keep_primary" %in% names(gene_filter),
  is.logical(gene_filter$keep_primary),
  !anyNA(gene_filter$keep_primary),
  sum(gene_filter$keep_primary) == 15012,
  sum(!gene_filter$keep_primary) == 53
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
  nrow(primary) == 101,
  length(unique(primary$donorid)) == 34
)

# ------------------------------------------------------------
# 4. Continuous-age complete-case cohort
# ------------------------------------------------------------

primary$age_numeric <- suppressWarnings(
  as.numeric(primary$age)
)

stopifnot(
  sum(is.na(primary$age_numeric)) == 2L
)

missing_age <- primary[
  is.na(primary$age_numeric),
  c(
    "geo_accession",
    "sample_title",
    "donorid",
    "diseasegroup",
    "lunglocation",
    "group",
    "age"
  ),
  drop = FALSE
]

stopifnot(
  nrow(missing_age) == 2L,
  all(missing_age$donorid == "ALF017"),
  all(missing_age$group == "NDC.Apex")
)

primary_age <- primary[
  !is.na(primary$age_numeric),
  ,
  drop = FALSE
]

stopifnot(
  nrow(primary_age) == 99L,
  length(unique(primary_age$donorid)) == 33L,
  !"ALF017" %in% primary_age$donorid,
  all(is.finite(primary_age$age_numeric))
)

# ------------------------------------------------------------
# 5. Frozen W7 group encoding
# ------------------------------------------------------------

primary_age$group <- factor(
  primary_age$group,
  levels = c(
    "IPF.Apex",
    "IPF.Base",
    "NDC.Apex",
    "NDC.Base"
  )
)

stopifnot(
  !anyNA(primary_age$group),
  identical(
    as.integer(table(primary_age$group)),
    c(30L, 31L, 20L, 18L)
  )
)

# ------------------------------------------------------------
# 6. Align counts to W7 age cohort
# ------------------------------------------------------------

stopifnot(
  all(primary_age$sample_title %in% colnames(counts))
)

counts_age <- counts[
  ,
  primary_age$sample_title,
  drop = FALSE
]

stopifnot(
  identical(
    colnames(counts_age),
    primary_age$sample_title
  )
)

idx <- match(
  gene_filter$gene_id,
  rownames(counts_age)
)

stopifnot(!anyNA(idx))

counts_age <- counts_age[
  idx,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(counts_age),
    gene_filter$gene_id
  )
)

counts_filtered <- counts_age[
  gene_filter$keep_primary,
  ,
  drop = FALSE
]

stopifnot(
  nrow(counts_filtered) == 15012L,
  ncol(counts_filtered) == 99L
)

# ------------------------------------------------------------
# 7. Prespecified W7 design
# ------------------------------------------------------------

design <- model.matrix(
  ~ 0 + group + age_numeric,
  data = primary_age
)

block <- factor(primary_age$donorid)

stopifnot(
  qr(design)$rank == ncol(design),
  identical(
    colnames(design),
    c(
      "groupIPF.Apex",
      "groupIPF.Base",
      "groupNDC.Apex",
      "groupNDC.Base",
      "age_numeric"
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
# Preserve frozen gene order. No ranking-based selection here.
# ------------------------------------------------------------

de_age <- topTable(
  fit,
  coef = "IPF_vs_NDC",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

de_age$gene_id <- rownames(de_age)

de_age <- de_age[
  ,
  c(
    "gene_id",
    setdiff(names(de_age), "gene_id")
  ),
  drop = FALSE
]

rownames(de_age) <- NULL

stopifnot(
  nrow(de_age) == 15012L,
  !anyDuplicated(de_age$gene_id),
  identical(
    de_age$gene_id,
    gene_filter$gene_id[gene_filter$keep_primary]
  ),
  all(is.finite(de_age$logFC)),
  all(is.finite(de_age$P.Value)),
  all(is.finite(de_age$adj.P.Val)),
  all(de_age$P.Value >= 0 & de_age$P.Value <= 1),
  all(de_age$adj.P.Val >= 0 & de_age$adj.P.Val <= 1)
)

# ------------------------------------------------------------
# 13. Global summaries
# ------------------------------------------------------------

n_tested <- nrow(de_age)

n_raw_005 <- sum(
  de_age$P.Value < 0.05
)

n_fdr_005 <- sum(
  de_age$adj.P.Val < 0.05
)

n_fdr_up <- sum(
  de_age$adj.P.Val < 0.05 &
    de_age$logFC > 0
)

n_fdr_down <- sum(
  de_age$adj.P.Val < 0.05 &
    de_age$logFC < 0
)

n_abs_logfc_1 <- sum(
  abs(de_age$logFC) >= 1
)

n_fdr_abs_logfc_1 <- sum(
  de_age$adj.P.Val < 0.05 &
    abs(de_age$logFC) >= 1
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
  "results/sensitivity/W7_continuous_age_IPF_vs_NDC_all_genes.csv"

summary_file <-
  "results/sensitivity/W7_continuous_age_IPF_vs_NDC_summary.csv"

qc_file <-
  "metadata/qc_W7_continuous_age_sensitivity.csv"

# ------------------------------------------------------------
# 15. Write outputs
# ------------------------------------------------------------

write.csv(
  de_age,
  result_file,
  row.names = FALSE
)

de_summary <- data.frame(
  metric = c(
    "tested_genes",
    "analysis_samples",
    "analysis_donors",
    "excluded_missing_age_samples",
    "excluded_missing_age_donors",
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
    nrow(primary_age),
    length(unique(primary_age$donorid)),
    nrow(missing_age),
    length(unique(missing_age$donorid)),
    rho1,
    rho2,
    n_raw_005,
    n_fdr_005,
    n_fdr_up,
    n_fdr_down,
    n_abs_logfc_1,
    n_fdr_abs_logfc_1,
    min(de_age$P.Value),
    min(de_age$adj.P.Val)
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
    "age_complete_samples",
    "primary_reference_donors",
    "age_complete_donors",
    "missing_age_samples",
    "missing_age_donor",
    "tested_genes",
    "design_columns",
    "design_rank",
    "contrast"
  ),
  value = c(
    "101",
    "99",
    "34",
    "33",
    "2",
    "ALF017",
    "15012",
    as.character(ncol(design)),
    as.character(qr(design)$rank),
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
  nrow(primary_age) == 99L,
  length(unique(primary_age$donorid)) == 33L,
  nrow(de_age) == 15012L,
  qr(design)$rank == ncol(design)
)

# ------------------------------------------------------------
# 17. Console certificate
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GSE213001 W7.2 — CONTINUOUS-AGE SENSITIVITY MODEL\n")
cat("============================================================\n")

cat("\n=== COHORT ===\n")
cat("Samples:                     ", nrow(primary_age), "\n", sep = "")
cat(
  "Donors:                      ",
  length(unique(primary_age$donorid)),
  "\n",
  sep = ""
)
cat("Missing-age samples removed: ", nrow(missing_age), "\n", sep = "")
cat("Missing-age donor removed:   ALF017\n")
cat("Age imputation:              NONE\n")

cat("\n=== DESIGN ===\n")
cat("Formula:   ~ 0 + group + age_numeric\n")
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
cat("The frozen W5 primary model was not modified.\n")
cat("No W6 enrichment or biological interpretation was modified.\n")
cat("The frozen W4/W5 15,012-gene universe was preserved.\n")
cat("No age value was imputed.\n")
cat("This is a W7 sensitivity analysis only.\n")

cat("\n============================================================\n")
cat("W7.2 CONTINUOUS-AGE SENSITIVITY MODEL: COMPLETED\n")
cat("============================================================\n")

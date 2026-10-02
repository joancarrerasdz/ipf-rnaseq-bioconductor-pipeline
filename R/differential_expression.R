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
